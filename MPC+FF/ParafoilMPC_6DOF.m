classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_6DOF_FF
    % フィードフォワード(FF)則を組み込んだ数値線形化MPC
    
    properties
        Ts                % 制御周期
        N                 % 予測ホライゾン
        ParamsList        % 物理パラメータ
        
        % 線形モデル (離散時間)
        Ad
        Bd
        
        % 重み行列
        Q, R, P
        
        % 状態保持
        LastU
        RefFullTable      % 参照軌道テーブル
        
        % 線形化オブジェクト
        Linearizer 
        DynamicsModel 
        LinearizationMode = 'Fixed';
        SoftMinK = 30;
        
        % ★追加: 誘導則(FF)用パラメータ
        FF_Gain_Factor
        Base_Trim
    end
    
    methods
        %% コンストラクタ
        function obj = ParafoilMPC_6DOF(ts, horizon, params)
            if nargin < 1, ts = 0.1; end
            if nargin < 2, horizon = 20; end
            
            obj.Ts = ts;
            obj.N  = horizon;
            obj.LastU = [0; 0];
            obj.ParamsList = params; 
            
            % --- ★追加: FFパラメータの計算 ---
            % K = - (b * Cnr) / (2 * Cnda)
            % 構造体 params から必要な係数を取得 (フィールド名は環境に合わせて調整してください)
            b = params.b;
            Cnr = params.C_n_r;      
            Cnda = params.C_n_delta_a;
            
            if abs(Cnda) < 1e-9, error('C_n_delta_a is too small.'); end
            obj.FF_Gain_Factor = - (b * Cnr) / (2 * Cnda);
            obj.Base_Trim = 0.05; % 滑空用トリム (デフォルト)
            
            % --- 重み行列 ---
            % 誤差に対する重み (Q)
            err_vec = [2.0, 1.0, 1.0, ...
                       deg2rad(10), deg2rad(10), deg2rad(10), ...
                       deg2rad(15), deg2rad(15), deg2rad(5), ...
                       5.0, 5.0, 10.0];
            w_vec = 1 ./ (err_vec.^2);
            obj.Q = diag(w_vec);
            obj.P = obj.Q;
            
            % 入力に対する重み (R)
            % ★変更: ここでのRは「FF値からの偏差」に対するペナルティとなる
            % FF項のおかげで、ここを小さくしなくても追従性が確保できる
            max_u_dev = 1; % FF値から0.4くらいズレても許容
            r_val  = 1 / max_u_dev^2;
            obj.R = diag([r_val, r_val]);
            
            % --- 線形化オブジェクト ---
            obj.DynamicsModel = ParafoilDynamics(params); 
            obj.Linearizer = ParafoilLinearizerNumeric(obj.DynamicsModel);
            if isprop(obj.Linearizer, 'soft_min_k')
                obj.Linearizer.soft_min_k = obj.SoftMinK;
            end
            
            obj.Ad = eye(12); 
            obj.Bd = zeros(12, 2);
        end
        
        %% Missionデータロード (変更なし)
        function load_reference(obj, missionObj)
            stateData = missionObj.export_6dof_state(); 
            trajGround = missionObj.export_detailed_trajectory('Ground');
            timeVec = trajGround.Time;
            
            obj.RefFullTable = [table(timeVec, 'VariableNames', {'Time'}), stateData];
            fprintf('MPC: Reference trajectory loaded (%d steps).\n', height(obj.RefFullTable));
            
            obj.update_linear_model(obj.RefFullTable(1, :));
        end
    
        %% Stepメソッド (QPソルバー実行)
        function [u_cmd, info] = step(obj, x_curr, current_time)
            
            if isempty(obj.RefFullTable)
                error('Reference trajectory not loaded.');
            end
            
            % 1. 線形モデル更新
            if strcmpi(obj.LinearizationMode, 'Adaptive')
                obj.update_linear_model_from_state(x_curr, obj.LastU);
            end
            
            % 2. ホライゾン参照値の取得
            % xref_horizon: [12 x N]
            xref_horizon = obj.get_horizon_ref(current_time, x_curr(9));
            xref_vec = xref_horizon(:); % [12N x 1]
            
            % --- ★追加: FF入力列の計算 ---
            % ホライゾン内の各時刻における理想入力 u_FF を計算
            U_ff_seq = obj.calc_horizon_ff(xref_horizon); % [2N x 1]
            
            % 3. QP行列構築 (FF項を f に反映)
            [H, f] = obj.construct_qp_matrices(x_curr, xref_vec, U_ff_seq);
            
            % 4. 制約条件
            lb = zeros(2 * obj.N, 1);
            ub = ones(2 * obj.N, 1);
            
            % ウォームスタート: 初期推定値を「FF入力」にするのが収束に良い
            u0_qp = U_ff_seq; 
            
            % 5. QP実行
            options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
            [U_seq, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, u0_qp, options);
            
            if exitflag < 0
                % 解なし時はFF値をそのまま使う (安全策)
                u_cmd = U_ff_seq(1:2)';
                obj.LastU = u_cmd(:);
            else
                u_cmd = U_seq(1:2)';
                obj.LastU = u_cmd(:);
            end
            
            info.ExitFlag = exitflag;
            info.RefState = xref_horizon(:, 1);
            info.FF_Input = U_ff_seq(1:2)'; % デバッグ用にFF値も返す
        end
        
        %% ★追加: ホライゾンFF入力計算
        function U_ff_vec = calc_horizon_ff(obj, xref_horizon)
            % xref_horizon: [12 x N]
            % 各予測ステップの参照状態から u_FF を計算し、縦ベクトルに結合する
            
            U_ff_vec = zeros(2 * obj.N, 1);
            
            for k = 1:obj.N
                % 状態抽出
                u_vel = xref_horizon(1, k);
                v_vel = xref_horizon(2, k);
                w_vel = xref_horizon(3, k);
                phi_ref = xref_horizon(7, k);
                theta_ref = xref_horizon(8, k);
                
                V_air = sqrt(u_vel^2 + v_vel^2 + w_vel^2);
                if V_air < 1.0, V_air = 15.0; end
                
                % 単一ステップのFF計算呼び出し
                u_k = obj.calc_guidance_single(phi_ref, V_air, theta_ref);
                
                % ベクトルに格納 [u_L; u_R]
                idx = (k-1)*2 + 1;
                U_ff_vec(idx : idx+1) = u_k(:);
            end
        end
        
        %% ★追加: 単一時刻の誘導則計算
        function u_val = calc_guidance_single(obj, phi_ref, V, theta)
            % 理論式: delta_a = K * (g/V^2) * cos(theta) * sin(phi)
            g = 9.81;
            K_dyn = (g / V^2) * cos(theta);
            delta_a = obj.FF_Gain_Factor * K_dyn * sin(phi_ref);
            
            trim = obj.Base_Trim;
            
            if delta_a > 0 % 右旋回 -> 右を引く
                u_L = trim;
                u_R = trim + delta_a;
            else % 左旋回
                u_L = trim + abs(delta_a);
                u_R = trim;
            end
            
            % クランプ
            u_val = max(0, min(1, [u_L; u_R]));
        end

        %% 内部: QP行列構築 (FF対応版)
        function [H, f] = construct_qp_matrices(obj, x0, xref_vec, U_ff_seq)
            n = 12; m = 2; Np = obj.N;
            
            % 予測行列 F, Phi の構築 (変更なし)
            F = zeros(n*Np, n);
            Phi = zeros(n*Np, m*Np);
            
            A_pow = eye(n);
            for i = 1:Np
                for j = 1:i
                    if i == j, term = obj.Bd;
                    else, term = (obj.Ad^(i-j)) * obj.Bd;
                    end
                    Phi((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = term;
                end
                A_pow = obj.Ad * A_pow;
                F((i-1)*n+1 : i*n, :) = A_pow;
            end
            
            % --- コスト関数の構築 ---
            % J = (X - Xref)^T Q_bar (X - Xref) + (U - U_ff)^T R_bar (U - U_ff)
            
            Q_bar = kron(eye(Np), obj.Q);
            R_bar = kron(eye(Np), obj.R);
            
            % 二次の項 (Hessian): U^T (Phi^T Q Phi + R) U
            H = 2 * (Phi' * Q_bar * Phi + R_bar);
            H = (H + H') / 2; 
            
            % 一次の項 (Linear term): f^T U
            % f = 2 * (x0^T F^T - Xref^T) Q Phi - 2 * U_ff^T R
            
            Y_free = F * x0;       
            E = Y_free - xref_vec; % (F x0 - Xref)
            
            % 部分1: 追従誤差による勾配
            f_track = (2 * E' * Q_bar * Phi)';
            
            % 部分2: ★FF入力との偏差による勾配 (- 2 R U_ff)
            f_ff = -2 * R_bar * U_ff_seq;
            
            % 合計
            f = f_track + f_ff;
        end
        
        %% (その他のメソッドは変更なし)
        function update_linear_model(obj, ref_row)
            y_trim = zeros(12, 1);
            y_trim(1) = ref_row.u; y_trim(2) = ref_row.v; y_trim(3) = ref_row.w;
            y_trim(4) = ref_row.p; y_trim(5) = ref_row.q; y_trim(6) = ref_row.r;
            y_trim(7) = ref_row.phi; y_trim(8) = ref_row.theta; y_trim(9) = ref_row.psi;
            y_trim(10)= ref_row.N; y_trim(11)= ref_row.E; y_trim(12)= ref_row.D;
            u_trim = [0; 0]; 
            extra.GAMMA = 0; extra.wind_I = [0;0;0];
            [A_c, B_c] = obj.Linearizer.get_linear_model(0, y_trim, u_trim, extra);
            sys_d = c2d(ss(A_c, B_c, eye(12), zeros(12,2)), obj.Ts);
            obj.Ad = sys_d.A; obj.Bd = sys_d.B;
        end
        
        function update_linear_model_from_state(obj, x_vec, u_prev)
            y_trim = x_vec; u_trim = u_prev;
            extra.GAMMA = 0; extra.wind_I = [0;0;0];
            [A_c, B_c] = obj.Linearizer.get_linear_model(0, y_trim, u_trim, extra);
            sys_d = c2d(ss(A_c, B_c, eye(12), zeros(12,2)), obj.Ts);
            obj.Ad = sys_d.A; obj.Bd = sys_d.B;
        end

        function xref_mat = get_horizon_ref(obj, t_curr, psi_curr)
            t_query = t_curr + (1:obj.N)' * obj.Ts;
            ref_time = obj.RefFullTable.Time;
            ref_data = obj.RefFullTable{:, 2:13}; 
            xref_mat = interp1(ref_time, ref_data, t_query, 'linear', 'extrap');
            xref_mat = xref_mat'; 
            idx_psi = 9;
            ref_psi = xref_mat(idx_psi, :);
            d_psi = ref_psi(1) - psi_curr;
            d_psi_wrapped = atan2(sin(d_psi), cos(d_psi));
            offset = psi_curr + d_psi_wrapped - ref_psi(1);
            xref_mat(idx_psi, :) = ref_psi + offset;
        end
    end
end