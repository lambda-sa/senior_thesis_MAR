classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_6DOF (Integrated Version)
    % 
    % 機能:
    % 1. ブライソンの法則を用いた重み付けによる12DOF MPC制御
    % 2. set_reference_trajectory による物理量補完 (位置→全状態)
    %
    % ※ファイル名は ParafoilMPC_6DOF.m のままでOKです
    
    properties
        MpcObj            % MPCオブジェクト
        MpcState          % MPC内部状態推定器
        Ts                % 制御周期 [s]
        PredictionHorizon % 予測ホライゾン [steps]
        ParamsList        % 機体パラメータ配列
        CurrentU          % 直前の入力 [delta_L, delta_R]
        
        RefFullTable      % ★補完された参照軌道データ (Table型)
    end
    
    methods
        %% 1. コンストラクタ (クラス名と一致させる)
        function obj = ParafoilMPC_6DOF(ts, prediction_horizon, params)
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
            
            % --- 線形モデル取得 (12状態) ---
            plant = obj.get_symbolic_linear_model_full(ts);
            
            % --- MPCオブジェクト作成 ---
            obj.MpcObj = mpc(plant, ts);
            obj.MpcObj.PredictionHorizon = prediction_horizon;
            obj.MpcObj.ControlHorizon = 3;

            % =================================================================
            % ★【重要修正】Kalmanフィルタエラーの回避
            % =================================================================
            % 12状態すべてを観測値(Output)としている場合、MPCがデフォルトで追加する
            % 「出力外乱モデル(積分器)」によってシステムが可観測でなくなり、
            % ss/kalman でエラーが発生します。
            % これを防ぐため、デフォルトの外乱モデルを削除('remove')します。
            setoutdist(obj.MpcObj, 'model', 'remove');
            % =================================================================

            % --- 制約条件 ---
            obj.MpcObj.MV(1).Min = 0; obj.MpcObj.MV(1).Max = 1.0;
            
            % --- 制約条件 ---
            obj.MpcObj.MV(1).Min = 0; obj.MpcObj.MV(1).Max = 1.0;
            obj.MpcObj.MV(2).Min = 0; obj.MpcObj.MV(2).Max = 1.0;
            
            % --- 重みづけ (Bryson's Rule) ---
            u_max = 1.0;            
            du_max_sec = 1.0;       
            du_max_step = du_max_sec * ts; 
            
            w_u  = 1 / u_max^2; 
            w_du = 1 / du_max_step^2;
            
            % 許容誤差の設定 (チューニングポイント)
            pos_err_max = 10.0;          % 位置 [m]
            psi_err_max = deg2rad(5.0);  % 方位 [rad]
            phi_err_max = deg2rad(10.0); % ロール [rad]
            
            w_pos = 1 / pos_err_max^2;
            w_psi = 1 / psi_err_max^2;
            w_phi = 1 / phi_err_max^2;
            w_ignore = 0.0; 
            
            % 重みベクトル: [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            weights_vec = [ ...
                w_ignore, w_ignore, w_ignore, ... 
                w_ignore, w_ignore, w_ignore, ... 
                w_phi,    w_ignore, w_psi,    ... 
                w_pos,    w_pos,    w_ignore ];
            
            obj.MpcObj.Weights.ManipulatedVariables = [w_u, w_u];
            obj.MpcObj.Weights.ManipulatedVariablesRate = [w_du, w_du];
            obj.MpcObj.Weights.OutputVariables = weights_vec;
            
            % --- 状態推定器初期化 ---
            obj.MpcState = mpcstate(obj.MpcObj);
        end
        
        %% 2. 参照軌道セット & 物理量自動計算 (★ここです！)
        function set_reference_trajectory(obj, time_vec, pos_ned, alpha_trim)
            if nargin < 4, alpha_trim = 0.1; end
            
            g = 9.81;
            dt = mean(diff(time_vec));
            if dt <= 0, error('Time vector must be monotonically increasing.'); end
            
            V_N = gradient(pos_ned(:,1), dt);
            V_E = gradient(pos_ned(:,2), dt);
            V_D = gradient(pos_ned(:,3), dt);
            V_horiz = sqrt(V_N.^2 + V_E.^2);
            
            psi_ref = unwrap(atan2(V_E, V_N)); 
            gamma_ref = atan2(-V_D, V_horiz);  
            
            dot_psi = gradient(psi_ref, dt);
            phi_ref = atan((V_horiz .* dot_psi) / g);
            theta_ref = gamma_ref + alpha_trim;
            
            u_ref = zeros(size(time_vec)); v_ref = zeros(size(time_vec)); w_ref = zeros(size(time_vec));
            p_ref = zeros(size(time_vec)); q_ref = zeros(size(time_vec)); r_ref = zeros(size(time_vec));
            
            dot_phi = gradient(phi_ref, dt);
            dot_theta = gradient(theta_ref, dt);
            
            for i = 1:length(time_vec)
                ph = phi_ref(i); th = theta_ref(i); ps = psi_ref(i);
                R_nb = [ cos(th)*cos(ps), cos(th)*sin(ps), -sin(th);
                         sin(ph)*sin(th)*cos(ps)-cos(ph)*sin(ps), sin(ph)*sin(th)*sin(ps)+cos(ph)*cos(ps), sin(ph)*cos(th);
                         cos(ph)*sin(th)*cos(ps)+sin(ph)*sin(ps), cos(ph)*sin(th)*sin(ps)-sin(ph)*cos(ps), cos(ph)*cos(th) ];
                V_ned_i = [V_N(i); V_E(i); V_D(i)];
                V_b = R_nb * V_ned_i;
                u_ref(i) = V_b(1); v_ref(i) = V_b(2); w_ref(i) = V_b(3);
                
                p_ref(i) = dot_phi(i) - dot_psi(i)*sin(th);
                q_ref(i) = dot_theta(i)*cos(ph) + dot_psi(i)*cos(th)*sin(ph);
                r_ref(i) = -dot_theta(i)*sin(ph) + dot_psi(i)*cos(th)*cos(ph);
            end
            
            data_mat = [u_ref, v_ref, w_ref, p_ref, q_ref, r_ref, ...
                        phi_ref, theta_ref, psi_ref, ...
                        pos_ned(:,1), pos_ned(:,2), pos_ned(:,3)];
            
            var_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            obj.RefFullTable = array2table([time_vec, data_mat], 'VariableNames', ['Time', var_names]);
            
            fprintf('Ref Trajectory Updated: %d steps.\n', length(time_vec));
        end
        
        %% 3. Step実行 (MPC計算)
        function [u_cmd, info] = step(obj, current_x_full, current_time)
            if isempty(obj.RefFullTable)
                error('参照軌道がセットされていません。set_reference_trajectoryを実行してください。');
            end
            
            y_measure = current_x_full;
            horizon_times = current_time + (0 : obj.PredictionHorizon)' * obj.Ts;
            
            ref_all = obj.RefFullTable{:, 2:end};
            ref_t   = obj.RefFullTable.Time;
            ref_signal = interp1(ref_t, ref_all, horizon_times, 'linear', 'extrap');
            
            % 方位角Unwrap
            idx_psi = 9;
            psi_curr = current_x_full(idx_psi);
            ref_psi_vec = ref_signal(:, idx_psi);
            diff_psi = ref_psi_vec(1) - psi_curr;
            diff_psi = atan2(sin(diff_psi), cos(diff_psi));
            start_ref_psi = psi_curr + diff_psi;
            ref_psi_relative = ref_psi_vec - ref_psi_vec(1);
            ref_signal(:, idx_psi) = start_ref_psi + ref_psi_relative;
            
            obj.MpcState.Plant = current_x_full; 
            [u_opt, info] = mpcmove(obj.MpcObj, obj.MpcState, y_measure, ref_signal);
            
            obj.CurrentU = u_opt;
            u_cmd = u_opt;
        end
        
        %% 4. 線形モデル生成 (Internal)
        function plant = get_symbolic_linear_model_full(obj, Ts)
            V_trim = 15.0; 
            x_trim = zeros(12, 1); x_trim(1) = V_trim; 
            u_trim = [0; 0];
            [A_full, B_full] = get_jacobians_sym(x_trim, u_trim, obj.ParamsList);
            
            plant_c = ss(A_full, B_full, eye(12), zeros(12,2));
            st_names = {'u','v','w', 'p','q','r', 'phi','theta','psi', 'N','E','D'};
            plant_c.StateName = st_names; plant_c.OutputName = st_names;
            plant_c.InputName = {'delta_L', 'delta_R'};
            plant = c2d(plant_c, Ts);
        end
    end
end