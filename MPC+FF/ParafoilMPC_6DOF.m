classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_SCRATCH
    % MPC Toolboxを使用せず、quadprogで直接QPを解く実装
    % 
    % 特徴:
    % - カルマンフィルタや外乱推定器を持たないため、特異行列エラーが起きない
    % - 数式が透明 (x_k+1 = Ax + Bu)
    
    properties
        Ts                % 制御周期
        N                 % 予測ホライゾン (Prediction Horizon)
        ParamsList
        
        % 線形モデル (離散時間)
        Ad
        Bd
        
        % 重み行列
        Q  % 状態重み (12x12)
        R  % 入力重み (2x2)
        P  % 終端重み (12x12)
        
        % 直前の入力
        LastU
        
        % 参照軌道
        RefFullTable
    end
    
    methods
        %% コンストラクタ
        function obj = ParafoilMPC_6DOF(ts, horizon, params)
            if nargin < 1, ts = 0.1; end
            if nargin < 2, horizon = 20; end
            
            obj.Ts = ts;
            obj.N  = horizon;
            obj.LastU = [0; 0];
            
            % --- パラメータ展開 ---
            m = params.m_total; S = params.S_c; b = params.b;
            if isfield(params, 'Ixx'), Ixx = params.Ixx; else, Ixx = m*(b/4)^2; end
            if isfield(params, 'Iyy'), Iyy = params.Iyy; else, Iyy = m*(b/5)^2; end
            if isfield(params, 'Izz'), Izz = params.Izz; else, Izz = m*(b/3)^2; end
            if isfield(params, 'Ixz'), Ixz = params.Ixz; else, Ixz = 0; end
            
            rho = 1.225; g = 9.81;
            CL0=params.C_L_0; CLa=params.C_L_alpha;
            CD0=params.C_D_0; CDa2=params.C_D_alpha;
            Clp=params.C_l_p; Cnr=params.C_n_r; Cld=params.C_l_delta_a; Cnd=params.C_n_delta_a;
            
            obj.ParamsList = [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CLa, CD0, CDa2, Clp, Cnr, Cld, Cnd];
            
            % --- 線形モデル構築 (初期化時に固定) ---
            % LTI-MPCとするため、コンストラクタで一度だけ離散化する
            [plant_d, ~] = obj.get_linear_model_discrete(ts);
            obj.Ad = plant_d.A;
            obj.Bd = plant_d.B;
            
            % --- 重み行列の設定 ---
            % Q: 状態偏差へのペナルティ
            % [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            w_vec = zeros(1, 12);
            w_vec(7)  = 10.0;  % phi
            w_vec(9)  = 100.0; % psi
            w_vec(10) = 10.0;  % N
            w_vec(11) = 10.0;  % E
            
            obj.Q = diag(w_vec);
            obj.P = obj.Q;         % 終端コスト (簡易的にQと同じ)
            obj.R = diag([1.0, 1.0]); % 入力コスト
        end
        
        %% Stepメソッド (QPソルバー実行)
        function [u_cmd, info] = step(obj, x_curr, current_time)
            
            if isempty(obj.RefFullTable)
                error('Reference trajectory not set.');
            end
            
            % 1. ホライゾン分の参照軌道を取得
            % [12 x N] の行列を作成
            xref_horizon = obj.get_horizon_ref(current_time, x_curr(9));
            
            % 参照軌道を縦ベクトルに変換 [12N x 1]
            xref_vec = xref_horizon(:);
            
            % 2. QP問題の行列構築 (Condensing)
            % 評価関数: J = 0.5 * U'*H*U + f'*U
            [H, f] = obj.construct_qp_matrices(x_curr, xref_vec);
            
            % 3. 制約条件
            % 入力制約: 0 <= u <= 1
            n_inputs = 2;
            lb = zeros(n_inputs * obj.N, 1);
            ub = ones(n_inputs * obj.N, 1);
            
            % オプション設定 (非表示)
            options = optimoptions('quadprog', 'Display', 'off');
            
            % 4. ソルバー実行
            % U_seq = [u_0; u_1; ...; u_{N-1}]
            [U_seq, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, [], options);
            
            if exitflag < 0
                warning('QP Solver failed. Using previous input.');
                u_cmd = obj.LastU';
            else
                % 最初のステップのみ採用
                u_cmd = U_seq(1:2)';
                obj.LastU = u_cmd(:);
            end
            
            info.ExitFlag = exitflag;
            info.FullSequence = U_seq;
        end
        
        %% 参照軌道設定
        function set_reference_trajectory(obj, time_vec, pos_ned, alpha_trim)
            % Toolbox版と同じロジックでテーブルを作成
             if nargin < 4, alpha_trim = 0.1; end
            
            dt = mean(diff(time_vec));
            V_N = gradient(pos_ned(:,1), dt);
            V_E = gradient(pos_ned(:,2), dt);
            V_D = gradient(pos_ned(:,3), dt);
            V_horiz = sqrt(V_N.^2 + V_E.^2);
            psi_ref = unwrap(atan2(V_E, V_N));
            dot_psi = gradient(psi_ref, dt);
            phi_ref = atan((V_horiz .* dot_psi) / 9.81);
            theta_ref = atan2(-V_D, V_horiz) + alpha_trim;
            
            u_ref = V_horiz; v_ref = zeros(size(u_ref)); w_ref = V_D;
            p_ref = zeros(size(u_ref)); q_ref = zeros(size(u_ref)); r_ref = dot_psi;
            
            data_mat = [u_ref, v_ref, w_ref, p_ref, q_ref, r_ref, ...
                        phi_ref, theta_ref, psi_ref, ...
                        pos_ned(:,1), pos_ned(:,2), pos_ned(:,3)];
            
            var_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            obj.RefFullTable = array2table([time_vec, data_mat], 'VariableNames', ['Time', var_names]);
        end
        
        %% 内部: QP行列構築 (最も重要な計算部分)
        function [H, f] = construct_qp_matrices(obj, x0, xref_vec)
            % 標準的なCondensing手法で H, f を構築します。
            % J = (Y - Yref)' Q_bar (Y - Yref) + U' R_bar U
            % Y = F*x0 + Phi*U
            
            n = 12; % 状態数
            m = 2;  % 入力数
            N_steps = obj.N;
            
            % 巨大行列の事前割り当て
            F   = zeros(n*N_steps, n);
            Phi = zeros(n*N_steps, m*N_steps);
            
            Q_bar = kron(eye(N_steps), obj.Q);
            % 終端コストを入れる場合:
            % Q_bar(end-n+1:end, end-n+1:end) = obj.P;
            
            R_bar = kron(eye(N_steps), obj.R);
            
            % 行列の構築 (Loop)
            A_pow = eye(n);
            for i = 1:N_steps
                A_pow = obj.Ad * A_pow;
                F((i-1)*n+1 : i*n, :) = A_pow;
                
                for j = 1:i
                    if i == j
                        coeff = obj.Bd;
                    else
                        coeff = (obj.Ad^(i-j)) * obj.Bd;
                    end
                    Phi((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = coeff;
                end
            end
            
            % QP形式: 0.5 * U' * H * U + f' * U
            % H = 2 * (Phi' * Q_bar * Phi + R_bar)
            % f = 2 * (x0' * F' - xref') * Q_bar * Phi
            
            H = 2 * (Phi' * Q_bar * Phi + R_bar);
            H = (H + H') / 2; % 対称性保証
            
            Y_free = F * x0; % 入力なしの応答
            E = Y_free - xref_vec; % 誤差
            
            f = (2 * E' * Q_bar * Phi)';
        end
        
        %% 内部: ホライゾン参照値取得
        function xref_mat = get_horizon_ref(obj, t_curr, psi_curr)
            horizon_times = t_curr + (1 : obj.N)' * obj.Ts; % 次のステップから
            
            ref_all = obj.RefFullTable{:, 2:end};
            ref_t   = obj.RefFullTable.Time;
            
            % 補間
            xref_mat = interp1(ref_t, ref_all, horizon_times, 'linear', 'extrap');
            
            % 転置して [12 x N] にする (ベクトル化しやすいように)
            xref_mat = xref_mat'; 
            
            % Psi Unwrap (Toolbox版と同じロジック)
            idx_psi = 9;
            ref_psi_vec = xref_mat(idx_psi, :);
            
            diff_psi = ref_psi_vec(1) - psi_curr;
            diff_psi = atan2(sin(diff_psi), cos(diff_psi));
            offset = psi_curr + diff_psi - ref_psi_vec(1);
            
            xref_mat(idx_psi, :) = ref_psi_vec + offset;
        end
        
        %% 内部: 線形モデル取得
        function [plant, plant_c] = get_linear_model_discrete(obj, Ts)
             V_trim = 15.0; 
            x_trim = zeros(12, 1); x_trim(1) = V_trim; 
            u_trim = [0; 0];
            
            [A_full, B_full] = get_jacobians_sym(x_trim, u_trim, obj.ParamsList);
            C_full = eye(12); D_full = zeros(12, 2);
            
            plant_c = ss(A_full, B_full, C_full, D_full);
            plant = c2d(plant_c, Ts);
        end
    end
end