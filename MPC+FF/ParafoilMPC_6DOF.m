classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_INTEGRATED (Full Version)
    % 
    % 機能:
    % 1. ブライソンの法則を用いた重み付けによる12DOF MPC制御
    % 2. 位置座標のみの参照軌道から、速度・姿勢・角速度を自動計算して補完
    %
    % 依存ファイル: get_jacobians_sym.m (generate_linear_model.mで生成)
    
    properties
        MpcObj            % MPCオブジェクト (Model Predictive Control Toolbox)
        MpcState          % MPC内部状態推定器
        Ts                % 制御周期 [s]
        PredictionHorizon % 予測ホライゾン [steps]
        ParamsList        % 機体パラメータ配列
        CurrentU          % 直前の入力 [delta_L, delta_R]
        
        RefFullTable      % 補完された参照軌道データ (Table型)
    end
    
    methods
        %% 1. コンストラクタ (初期化と重み設定)
        function obj = ParafoilMPC_6DOF(ts, prediction_horizon, params)
            % 引数:
            %   ts: 制御周期 (例: 0.1)
            %   prediction_horizon: 予測ステップ数 (例: 20)
            %   params: 機体諸元構造体
            
            if nargin < 1, ts = 0.1; end
            if nargin < 2, prediction_horizon = 20; end
            
            obj.Ts = ts;
            obj.PredictionHorizon = prediction_horizon;
            obj.CurrentU = [0, 0];
            
            % --- 依存ファイルの確認 ---
            if exist('get_jacobians_sym', 'file') ~= 2
                error('エラー: get_jacobians_sym.m が見つかりません。generate_linear_model.m を実行してください。');
            end
            
            % --- パラメータ展開 ---
            % ※ユーザー環境に合わせてフィールド名を調整してください
            m = params.m_total; S = params.S_c; b = params.b;
            if isfield(params, 'Ixx'), Ixx = params.Ixx; else, Ixx = m*(b/4)^2; end
            if isfield(params, 'Iyy'), Iyy = params.Iyy; else, Iyy = m*(b/5)^2; end
            if isfield(params, 'Izz'), Izz = params.Izz; else, Izz = m*(b/3)^2; end
            if isfield(params, 'Ixz'), Ixz = params.Ixz; else, Ixz = 0; end
            rho = 1.225; g = 9.81;
            
            % 空力係数
            CL0=params.C_L_0; CLa=params.C_L_alpha;
            CD0=params.C_D_0; CDa2=params.C_D_alpha;
            Clp=params.C_l_p; Cnr=params.C_n_r; Cld=params.C_l_delta_a; Cnd=params.C_n_delta_a;
            
            obj.ParamsList = [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CLa, CD0, CDa2, Clp, Cnr, Cld, Cnd];
            
            % --- 線形モデル取得 (12状態) ---
            plant = obj.get_symbolic_linear_model_full(ts);
            
            % --- MPCオブジェクト作成 ---
            obj.MpcObj = mpc(plant, ts);
            obj.MpcObj.PredictionHorizon = prediction_horizon;
            obj.MpcObj.ControlHorizon = 3;
            
            % --- 制約条件 (Constraints) ---
            % 入力範囲: 0.0 (全解除) ~ 1.0 (フルブレーク)
            obj.MpcObj.MV(1).Min = 0; obj.MpcObj.MV(1).Max = 1.0;
            obj.MpcObj.MV(2).Min = 0; obj.MpcObj.MV(2).Max = 1.0;
            
            % --- 重みづけ (Bryson's Rule) ---
            % Q_ii = 1 / (max_error)^2
            
            % 1. 入力 (Manipulated Variables)
            u_max = 1.0;            % 入力の最大振幅
            du_max_sec = 1.0;       % 1秒間に許容する最大変化量 (1.0/s)
            du_max_step = du_max_sec * ts; 
            
            w_u  = 1 / u_max^2; 
            w_du = 1 / du_max_step^2;
            
            % 2. 出力 (Output Variables: 12 states)
            % 許容誤差の設定
            pos_err_max = 10.0;          % 位置許容誤差 [m]
            psi_err_max = deg2rad(5.0);  % 方位許容誤差 [rad]
            phi_err_max = deg2rad(10.0); % ロール許容誤差 [rad] (少し緩める)
            
            w_pos = 1 / pos_err_max^2;
            w_psi = 1 / psi_err_max^2;
            w_phi = 1 / phi_err_max^2;
            w_ignore = 0.0; % 制御しない変数は重みゼロ
            
            % 出力重みベクトルの構成
            % 順序: [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            weights_vec = [ ...
                w_ignore, w_ignore, w_ignore, ... % u, v, w (速度は成行)
                w_ignore, w_ignore, w_ignore, ... % p, q, r (角速度は成行)
                w_phi,    w_ignore, w_psi,    ... % phi(制御), theta, psi(最重要)
                w_pos,    w_pos,    w_ignore      % N(重要), E(重要), D(成行)
            ];
            
            % 重みの適用
            obj.MpcObj.Weights.ManipulatedVariables = [w_u, w_u];
            obj.MpcObj.Weights.ManipulatedVariablesRate = [w_du, w_du]; % 変化率抑制
            obj.MpcObj.Weights.OutputVariables = weights_vec;
            
            % --- 状態推定器初期化 ---
            obj.MpcState = mpcstate(obj.MpcObj);
        end
        
        %% 2. 参照軌道セット & 物理量自動計算
        function set_reference_trajectory(obj, time_vec, pos_ned, alpha_trim)
            % set_reference_trajectory
            % 入力: 
            %   time_vec: 時間 [T x 1]
            %   pos_ned : 位置 [T x 3] (North, East, Down)
            %   alpha_trim: トリム迎え角 [rad] (通常 0.1 程度)
            
            if nargin < 4, alpha_trim = 0.1; end
            
            g = 9.81;
            dt = mean(diff(time_vec));
            if dt <= 0, error('Time vector must be monotonically increasing.'); end
            
            % A. 速度ベクトル (数値微分)
            V_N = gradient(pos_ned(:,1), dt);
            V_E = gradient(pos_ned(:,2), dt);
            V_D = gradient(pos_ned(:,3), dt);
            V_horiz = sqrt(V_N.^2 + V_E.^2);
            
            % B. 方位角(psi) と 経路角(gamma)
            psi_ref = unwrap(atan2(V_E, V_N)); % 2piジャンプを除去
            gamma_ref = atan2(-V_D, V_horiz);  % 上昇角 (パラフォイルは負)
            
            % C. 定常旋回理論によるロール角(phi)計算
            % tan(phi) = (V * dot_psi) / g
            dot_psi = gradient(psi_ref, dt);
            phi_ref = atan((V_horiz .* dot_psi) / g);
            
            % ピッチ角 (gamma + alpha)
            theta_ref = gamma_ref + alpha_trim;
            
            % D. 機体座標系速度 (u,v,w) への変換
            u_ref = zeros(size(time_vec));
            v_ref = zeros(size(time_vec));
            w_ref = zeros(size(time_vec));
            
            % E. 角速度 (p,q,r) の計算
            dot_phi = gradient(phi_ref, dt);
            dot_theta = gradient(theta_ref, dt);
            p_ref = zeros(size(time_vec));
            q_ref = zeros(size(time_vec));
            r_ref = zeros(size(time_vec));
            
            for i = 1:length(time_vec)
                ph = phi_ref(i); th = theta_ref(i); ps = psi_ref(i);
                
                % DCM (NED -> Body)
                R_nb = [ cos(th)*cos(ps), cos(th)*sin(ps), -sin(th);
                         sin(ph)*sin(th)*cos(ps)-cos(ph)*sin(ps), sin(ph)*sin(th)*sin(ps)+cos(ph)*cos(ps), sin(ph)*cos(th);
                         cos(ph)*sin(th)*cos(ps)+sin(ph)*sin(ps), cos(ph)*sin(th)*sin(ps)-sin(ph)*cos(ps), cos(ph)*cos(th) ];
                
                V_ned_i = [V_N(i); V_E(i); V_D(i)];
                V_b = R_nb * V_ned_i;
                u_ref(i) = V_b(1); v_ref(i) = V_b(2); w_ref(i) = V_b(3);
                
                % Body Rates (Kinematics)
                p_ref(i) = dot_phi(i) - dot_psi(i)*sin(th);
                q_ref(i) = dot_theta(i)*cos(ph) + dot_psi(i)*cos(th)*sin(ph);
                r_ref(i) = -dot_theta(i)*sin(ph) + dot_psi(i)*cos(th)*cos(ph);
            end
            
            % F. データの格納
            % [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            data_mat = [u_ref, v_ref, w_ref, p_ref, q_ref, r_ref, ...
                        phi_ref, theta_ref, psi_ref, ...
                        pos_ned(:,1), pos_ned(:,2), pos_ned(:,3)];
            
            var_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            obj.RefFullTable = array2table([time_vec, data_mat], 'VariableNames', ['Time', var_names]);
            
            fprintf('Ref Trajectory Updated: %d steps. (Full State Calculated)\n', length(time_vec));
        end
        
        %% 3. Step実行 (MPC計算)
        function [u_cmd, info] = step(obj, current_x_full, current_time)
            % current_x_full: [u,v,w, p,q,r, phi,theta,psi, N,E,D] (12x1)
            
            if isempty(obj.RefFullTable)
                error('参照軌道がセットされていません。set_reference_trajectoryを実行してください。');
            end
            
            % 観測値 (Full State Feedback)
            y_measure = current_x_full;
            
            % 予測ホライゾン分の時刻作成
            horizon_times = current_time + (0 : obj.PredictionHorizon)' * obj.Ts;
            
            % 参照データの補間
            % Tableから行列を取り出し (Time列を除く2列目以降)
            ref_all = obj.RefFullTable{:, 2:end};
            ref_t   = obj.RefFullTable.Time;
            
            % 線形補間
            ref_signal = interp1(ref_t, ref_all, horizon_times, 'linear', 'extrap');
            
            % --- 方位角(psi)の連続性補正 (Unwrap) ---
            % 9列目が psi と仮定
            idx_psi = 9;
            psi_curr = current_x_full(idx_psi);
            ref_psi_vec = ref_signal(:, idx_psi);
            
            % 現在値と参照値の差分を -pi~pi に収める
            diff_psi = ref_psi_vec(1) - psi_curr;
            diff_psi = atan2(sin(diff_psi), cos(diff_psi));
            
            % オフセットを適用して滑らかに繋ぐ
            start_ref_psi = psi_curr + diff_psi;
            ref_psi_relative = ref_psi_vec - ref_psi_vec(1);
            ref_signal(:, idx_psi) = start_ref_psi + ref_psi_relative;
            
            % --- MPC計算 ---
            obj.MpcState.Plant = current_x_full; % 内部モデルの状態更新
            [u_opt, info] = mpcmove(obj.MpcObj, obj.MpcState, y_measure, ref_signal);
            
            obj.CurrentU = u_opt;
            u_cmd = u_opt;
        end
        
        %% 4. 線形モデル生成 (Internal)
        function plant = get_symbolic_linear_model_full(obj, Ts)
            % トリム条件 (動作点)
            V_trim = 15.0; % [m/s]
            x_trim = zeros(12, 1);
            x_trim(1) = V_trim; 
            u_trim = [0; 0];
            
            % 生成済みのシンボリック関数呼び出し
            [A_full, B_full] = get_jacobians_sym(x_trim, u_trim, obj.ParamsList);
            
            % 12状態フル観測モデル
            C_full = eye(12);
            D_full = zeros(12, 2);
            
            plant_c = ss(A_full, B_full, C_full, D_full);
            
            st_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            plant_c.StateName = st_names;
            plant_c.OutputName = st_names;
            plant_c.InputName = {'delta_L', 'delta_R'};
            
            % 離散化
            plant = c2d(plant_c, Ts);
        end
    end
end