classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_6DOF (Integrated Complete Version)
    % 
    % 機能要約:
    % 1. 12自由度 (u,v,w, p,q,r, phi,theta,psi, N,E,D) フル状態MPC
    % 2. 参照軌道(位置)から物理的に整合する姿勢・速度を自動生成 (Feedforward効果)
    % 3. ブライソンの法則による重み調整
    % 4. Kalmanフィルタエラー回避 (Output Disturbance Model除去)
    %
    % 依存ファイル: get_jacobians_sym.m (generate_linear_model.mで生成)
    
    properties
        MpcObj            % MPCオブジェクト (Model Predictive Control Toolbox)
        MpcState          % MPC内部状態推定器
        Ts                % 制御周期 [s]
        PredictionHorizon % 予測ホライゾン [steps]
        ParamsList        % 機体パラメータ配列
        CurrentU          % 直前の入力 [delta_L, delta_R]
        
        RefFullTable      % 補完・計算済みの参照軌道データ (Table型)
    end
    
    methods
        %% 1. コンストラクタ
        function obj = ParafoilMPC_6DOF(ts, prediction_horizon, params)
            % コンストラクタ
            % ts: 制御周期 (例: 0.1)
            % prediction_horizon: 予測ステップ数 (例: 20)
            % params: 機体諸元構造体
            
            if nargin < 1, ts = 0.1; end
            if nargin < 2, prediction_horizon = 20; end
            
            obj.Ts = ts;
            obj.PredictionHorizon = prediction_horizon;
            obj.CurrentU = [0, 0];
            
            % --- 依存ファイルの確認 ---
            if exist('get_jacobians_sym', 'file') ~= 2
                error('エラー: "get_jacobians_sym.m" が見つかりません。先に "generate_linear_model.m" を実行してください。');
            end
            
            % --- パラメータリストの構築 ---
            m = params.m_total; S = params.S_c; b = params.b;
            % 慣性モーメント
            if isfield(params, 'Ixx'), Ixx = params.Ixx; else, Ixx = m*(b/4)^2; end
            if isfield(params, 'Iyy'), Iyy = params.Iyy; else, Iyy = m*(b/5)^2; end
            if isfield(params, 'Izz'), Izz = params.Izz; else, Izz = m*(b/3)^2; end
            if isfield(params, 'Ixz'), Ixz = params.Ixz; else, Ixz = 0; end
            
            % 空力・環境定数
            rho = 1.225; g = 9.81;
            CL0=params.C_L_0; CLa=params.C_L_alpha;
            CD0=params.C_D_0; CDa2=params.C_D_alpha;
            Clp=params.C_l_p; Cnr=params.C_n_r; Cld=params.C_l_delta_a; Cnd=params.C_n_delta_a;
            
            obj.ParamsList = [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CLa, CD0, CDa2, Clp, Cnr, Cld, Cnd];
            
            % --- 線形モデルの取得 (12状態フルモデル) ---
            plant = obj.get_symbolic_linear_model_full(ts);
            
            % --- MPCオブジェクト作成 ---
            obj.MpcObj = mpc(plant, ts);
            
            % --- 基本設定 ---
            obj.MpcObj.PredictionHorizon = prediction_horizon;
            obj.MpcObj.ControlHorizon = 3; 
            
          % =================================================================
            % ★状態推定器の無効化 (Custom化)
            % =================================================================
            % ブライソンの法則は「制御則(Q,R)」には効いていますが、
            % 「推定器(Kalman Filter)」の内部スケーリングには影響しません。
            % 桁違いの状態量(位置vs角度)による推定エラーを防ぐため、
            % 推定器をバイパス(custom)設定にします。
            
            setEstimator(obj.MpcObj, 'custom');
            
            % (以前の setoutdist の記述は削除してください)
            % =================================================================
            
            % --- 制約条件 (Constraints) ---
            % 入力: 0.0 (全解除) ~ 1.0 (フルブレーク)
            obj.MpcObj.MV(1).Min = 0; obj.MpcObj.MV(1).Max = 1.0;
            obj.MpcObj.MV(2).Min = 0; obj.MpcObj.MV(2).Max = 1.0;
            
            % レートリミット (急操作の制限)
            obj.MpcObj.MV(1).RateMin = -2.0; obj.MpcObj.MV(1).RateMax = 2.0;
            obj.MpcObj.MV(2).RateMin = -2.0; obj.MpcObj.MV(2).RateMax = 2.0;
            
            % =================================================================
            % ★ブライソンの法則 (Bryson's Rule) による重みづけ
            % =================================================================
            
            % 1. 入力 (MV) のスケーリング
            % ---------------------------
            max_u_val  = 1.0;            % 入力最大値
            max_du_sec = 1.0;            % 1秒間の許容操作量 (1.0/s)
            max_du_step = max_du_sec * ts; 
            
            w_u   = 1 / (max_u_val^2);
            w_du  = 1 / (max_du_step^2);
            
            % 2. 出力 (OV) のスケーリング
            % ---------------------------
            % 制御したい変数の許容誤差を設定します
            
            % 位置 (N, E): 10m 以内の誤差を許容
            max_pos_err = 10.0; 
            w_pos = 1 / (max_pos_err^2);
            
            % 方位角 (psi): 5度 以内の誤差を許容
            max_psi_err = deg2rad(5.0);
            w_psi = 1 / (max_psi_err^2);
            
            % ロール角 (phi): 10度 以内のふらつきを許容 (安定性重視)
            max_phi_err = deg2rad(10.0);
            w_phi = 1 / (max_phi_err^2);
            
            % その他 (速度 u,v,w や 高度 D): 制御せず物理法則に任せる
            w_ignore = 0.0; 
            
            % --- 重みベクトルの構築 (12要素) ---
            % 順序: [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            weights_vec = [ ...
                w_ignore, w_ignore, w_ignore, ... % u, v, w
                w_ignore, w_ignore, w_ignore, ... % p, q, r
                w_phi,    w_ignore, w_psi,    ... % phi(制御), theta, psi(最重要)
                w_pos,    w_pos,    w_ignore      % N(重要), E(重要), D(成行)
            ];
            
            % --- 重みの適用 ---
            % 入力(MV): 小さくして動きやすくする
            obj.MpcObj.Weights.ManipulatedVariables = [w_u, w_u];
            % 入力変化率(MVRate): ガタつきを抑える
            obj.MpcObj.Weights.ManipulatedVariablesRate = [w_du, w_du];
            obj.MpcObj.Weights.OutputVariables = weights_vec;
            
            % --- 状態推定器の初期化 ---
            obj.MpcState = mpcstate(obj.MpcObj);
        end
        
        %% 2. 参照軌道セット & 物理量自動計算 (統合メソッド)
        function set_reference_trajectory(obj, time_vec, pos_ned, alpha_trim)
            % set_reference_trajectory
            % 位置だけの軌道データから、定常旋回理論に基づいて
            % 速度・姿勢・角速度を逆算し、内部テーブルに保存します。
            % これによりフィードフォワード的な先行制御が可能になります。
            %
            % Inputs:
            %   time_vec:   時間 [T x 1]
            %   pos_ned:    位置 [T x 3] (North, East, Down)
            %   alpha_trim: トリム迎え角 [rad] (default: 0.1)
            
            if nargin < 4, alpha_trim = 0.1; end
            
            g = 9.81;
            dt = mean(diff(time_vec));
            if dt <= 0, error('Time vector must be monotonically increasing.'); end
            
            % A. 速度ベクトルの計算
            V_N = gradient(pos_ned(:,1), dt);
            V_E = gradient(pos_ned(:,2), dt);
            V_D = gradient(pos_ned(:,3), dt);
            V_horiz = sqrt(V_N.^2 + V_E.^2);
            
            % B. 方位角(psi) と 経路角(gamma)
            psi_ref = unwrap(atan2(V_E, V_N));
            gamma_ref = atan2(-V_D, V_horiz);
            
            % C. 旋回率とロール角 (Coordinated Turn Logic)
            % tan(phi) = (V * dot_psi) / g
            dot_psi = gradient(psi_ref, dt);
            phi_ref = atan((V_horiz .* dot_psi) / g);
            
            % ピッチ角
            theta_ref = gamma_ref + alpha_trim;
            
            % D. 機体座標系速度 (u,v,w) への変換
            u_ref = zeros(size(time_vec));
            v_ref = zeros(size(time_vec));
            w_ref = zeros(size(time_vec));
            
            % E. 角速度 (p,q,r) の計算 (Kinematics)
            dot_phi = gradient(phi_ref, dt);
            dot_theta = gradient(theta_ref, dt);
            p_ref = zeros(size(time_vec));
            q_ref = zeros(size(time_vec));
            r_ref = zeros(size(time_vec));
            
            for i = 1:length(time_vec)
                ph = phi_ref(i); th = theta_ref(i); ps = psi_ref(i);
                
                % DCM (NED -> Body) 3-2-1 Rotation
                R_nb = [ cos(th)*cos(ps), cos(th)*sin(ps), -sin(th);
                         sin(ph)*sin(th)*cos(ps)-cos(ph)*sin(ps), sin(ph)*sin(th)*sin(ps)+cos(ph)*cos(ps), sin(ph)*cos(th);
                         cos(ph)*sin(th)*cos(ps)+sin(ph)*sin(ps), cos(ph)*sin(th)*sin(ps)-sin(ph)*cos(ps), cos(ph)*cos(th) ];
                     
                V_ned_i = [V_N(i); V_E(i); V_D(i)];
                V_b = R_nb * V_ned_i;
                u_ref(i)=V_b(1); v_ref(i)=V_b(2); w_ref(i)=V_b(3);
                
                % Angular Rates
                p_ref(i) = dot_phi(i) - dot_psi(i)*sin(th);
                q_ref(i) = dot_theta(i)*cos(ph) + dot_psi(i)*cos(th)*sin(ph);
                r_ref(i) = -dot_theta(i)*sin(ph) + dot_psi(i)*cos(th)*cos(ph);
            end
            
            % F. データの格納
            data_mat = [u_ref, v_ref, w_ref, p_ref, q_ref, r_ref, ...
                        phi_ref, theta_ref, psi_ref, ...
                        pos_ned(:,1), pos_ned(:,2), pos_ned(:,3)];
            
            var_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            obj.RefFullTable = array2table([time_vec, data_mat], 'VariableNames', ['Time', var_names]);
            
            fprintf('Ref Trajectory Updated: %d steps processed.\n', length(time_vec));
        end
        
        %% 3. Step実行メソッド (Run MPC)
        function [u_cmd, info] = step(obj, current_x_full, current_time)
            % step
            % 内部参照テーブルを用いてMPC計算を実行
            
            if isempty(obj.RefFullTable)
                error('Reference trajectory not set. Call set_reference_trajectory() first.');
            end
            
            y_measure = current_x_full;
            
            % 予測ホライゾン分の参照値取得
            horizon_times = current_time + (0 : obj.PredictionHorizon)' * obj.Ts;
            
            % データ補間 (Time列を除く全データ)
            ref_all = obj.RefFullTable{:, 2:end};
            ref_t   = obj.RefFullTable.Time;
            
            ref_signal = interp1(ref_t, ref_all, horizon_times, 'linear', 'extrap');
            
            % --- 方位角(psi)の Unwrap (2pi問題の解決) ---
            % 9列目が psi である前提
            idx_psi = 9;
            psi_curr = current_x_full(idx_psi);
            ref_psi_vec = ref_signal(:, idx_psi);
            
            % 現在値に近い角度にシフト (-pi ~ pi)
            diff_psi = ref_psi_vec(1) - psi_curr;
            diff_psi = atan2(sin(diff_psi), cos(diff_psi));
            start_ref_psi = psi_curr + diff_psi;
            
            % 全体にシフトを適用
            ref_psi_relative = ref_psi_vec - ref_psi_vec(1);
            ref_signal(:, idx_psi) = start_ref_psi + ref_psi_relative;
            
            % --- MPC計算実行 ---
            obj.MpcState.Plant = current_x_full;
            [u_opt, info] = mpcmove(obj.MpcObj, obj.MpcState, y_measure, ref_signal);
            
            obj.CurrentU = u_opt;
            u_cmd = u_opt;
        end
        
        %% 4. 線形モデル取得 (内部メソッド)
        function plant = get_symbolic_linear_model_full(obj, Ts)
            % トリム条件 (動作点)
            V_trim = 15.0; 
            x_trim = zeros(12, 1); x_trim(1) = V_trim; 
            u_trim = [0; 0];
            
            % ヤコビアン取得
            [A_full, B_full] = get_jacobians_sym(x_trim, u_trim, obj.ParamsList);
            
            % 12状態フル観測モデル
            % C = I (全状態出力), D = 0
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