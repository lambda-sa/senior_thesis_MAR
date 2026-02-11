classdef ParafoilMission < handle
    % PARAFOILMISSION (Updated for Pure Dubins & Run-up Integration)
    
    properties
        ExcelFileName
        Params        % パラメータ構造体
        SimSettings   % 初期設定
        
        AtmoModel     % 大気モデル
        Planner       % 経路計画クラス (ParafoilPathPlanner)
        PhysicsParams % 物理パラメータ
        
        % 軌道ログ保存用
        TrajectoryLogAir table
        TrajectoryLogGround table
        
        ManualSpeedCorrectionFlag = [];
    end
    
    methods
        function obj = ParafoilMission(excelFileName)
            obj.ExcelFileName = excelFileName;
            
            % 大気モデル初期化
            try
                obj.AtmoModel = AtmoTempPressRho();
            catch
                obj.AtmoModel = AtmoModel_Standard(); 
            end
            
            % Planner初期化
            obj.Planner = ParafoilPathPlanner(obj.AtmoModel);
        end
        
        % 実行スクリプト等からフラグを設定するための関数
        function set_speed_correction(obj, enable_flag)
            obj.ManualSpeedCorrectionFlag = logical(enable_flag);
            if ~isempty(obj.Planner)
                obj.Planner.set_speed_correction(obj.ManualSpeedCorrectionFlag);
            end
            fprintf('Mission Setting: Bank Speed Correction manually set to %d\n', obj.ManualSpeedCorrectionFlag);
        end
        
        % 物理パラメータ計算（Excel読み込み含む）
        function compute_physics_parameters(obj)
            [params, sim_settings] = load_params_from_excel(obj.ExcelFileName);
            obj.Params = params;
            obj.SimSettings = sim_settings;
            
            % 速度補正フラグ設定
            if ~isempty(obj.ManualSpeedCorrectionFlag)
                use_cor = obj.ManualSpeedCorrectionFlag;
            elseif isfield(sim_settings, 'enable_bank_speed_correction')
                use_cor = logical(sim_settings.enable_bank_speed_correction);
            else
                use_cor = false; 
            end
            obj.Planner.EnableBankSpeedCorrection = use_cor;
            
            % 初期高度・位置
            if isfield(sim_settings, 'h_init'), h_start = sim_settings.h_init;
            elseif isfield(sim_settings, 'z_initial'), h_start = -sim_settings.z_initial;
            else, h_start = 1000; end
            
            sx = -500; sy = 0; syaw = 0; land_dir = 180;
            if isfield(sim_settings, 'X_initial'), sx = sim_settings.X_initial; end
            if isfield(sim_settings, 'Y_initial'), sy = sim_settings.Y_initial; end
            if isfield(sim_settings, 'psi_initial_deg'), syaw = deg2rad(sim_settings.psi_initial_deg); end
            if isfield(sim_settings, 'landing_direction'), land_dir = sim_settings.landing_direction; end
            
            % トリム計算
            dyn = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, obj.AtmoModel);
            mu_rigging = 0;
            if isfield(sim_settings, 'ang_gamma_rigging'), mu_rigging = sim_settings.ang_gamma_rigging; end
            
            try
                [~, gamma_trim, V_trim_TAS, ~, ~] = dyn.calculate_initial_trim(h_start, 0, 0, mu_rigging);
            catch
                V_trim_TAS = 15.0; gamma_trim = deg2rad(-15);
            end
            
            % EAS換算
            rho_start = obj.AtmoModel.get_density(h_start/1000);
            rho_0 = 1.225;
            V_EAS = V_trim_TAS * sqrt(rho_start / rho_0);
            gr = -1.0 / tan(gamma_trim);
            
            % 旋回半径計算 (TAS基準)
            g = obj.AtmoModel.get_gravity(h_start);
            if isfield(params.guidance, 'nominal_bank_deg')
                phi_plan = deg2rad(params.guidance.nominal_bank_deg);
            else
                phi_plan = deg2rad(15); % デフォルトバンク角
            end
            R_calc = V_trim_TAS^2 / (g * tan(phi_plan));
            
            % 保存
            pp.start_pos = [sx, sy, h_start, syaw];
            pp.landing_direction = land_dir;
            pp.V_air = V_EAS;
            pp.glide_ratio = gr;
            pp.min_turn_radius = R_calc;
            
            obj.PhysicsParams = pp;
            
            fprintf('Physics Parameters Computed:\n');
            fprintf('  Start Alt: %.0fm, Run-up: %.0fm\n', h_start, obj.Planner.RunUpDistance);
            fprintf('  V_TAS: %.2f m/s, L/D: %.2f, Turn Radius: %.1f m\n', V_trim_TAS, gr, R_calc);
        end
        
        % =========================================================
        % ★重要: 助走付きPure Dubinsシミュレーション実行メソッド
        % =========================================================
        function run_wind_simulation(obj, target_pos, L_final, wind_vec_2d)
            % 1. 物理パラメータ更新
            obj.compute_physics_parameters();
            
            % 2. Plannerに風を設定
            obj.Planner.WindVector = [wind_vec_2d(1); wind_vec_2d(2); 0];
            
            % 3. パラメータ展開
            pp = obj.PhysicsParams;
            land_dir = pp.landing_direction;

            fprintf('Running Pure Dubins Planning (with Run-up)...\n');
            
            % 4. 経路生成 (calc_pure_dubins を呼び出し)
            % ※ 旋回方向は探索ロジック(plan_optimal_trajectory)経由が理想ですが、
            %    ここではシンプルに直接Dubins計算を呼ぶか、Plannerのメインメソッドを呼びます。
            %    (Planner側で plan_optimal_trajectory が calc_pure_dubins を呼ぶように修正済みと想定)
            
            % Plannerの plan_optimal_trajectory を呼ぶ (内部で calc_pure_dubins が呼ばれる)
            obj.Planner.plan_optimal_trajectory(...
                pp.start_pos, ...
                target_pos, ...
                land_dir, ...
                L_final, ...
                pp.V_air, ...
                pp.glide_ratio, ...
                pp.min_turn_radius);
                
            % 5. ログ生成 (新しいデータ構造に対応)
            obj.compute_dual_trajectories();
        end
        
        % =========================================================
        % ★重要: ログ生成 (runup, loiter, dubins, final 対応版)
        % =========================================================
        function compute_dual_trajectories(obj)
            d = obj.Planner.ResultData;
            if isempty(d), error('No Path Data'); end
            
            % データの結合 (runup が含まれていることに注意)
            fs = {'runup', 'loiter', 'dubins', 'final'};
            t=[]; x=[]; y=[]; z=[]; psi=[]; kappa=[];
            
            for i=1:length(fs)
                f=fs{i};
                if isfield(d, f) && ~isempty(d.(f).x)
                    t=[t; d.(f).t(:)]; 
                    x=[x; d.(f).x(:)]; y=[y; d.(f).y(:)]; z=[z; d.(f).z(:)];
                    psi=[psi; d.(f).psi(:)];
                    kappa=[kappa; d.(f).kappa(:)];
                end
            end
            
            if isempty(t), return; end
            
            % 重複削除
            [t_unique, idx] = unique(t, 'stable');
            x=x(idx); y=y(idx); z=z(idx); psi=psi(idx); kappa=kappa(idx);
            
            num_steps = length(t_unique);
            
            % 風ベクトル
            Wx = obj.Planner.WindVector(1); Wy = obj.Planner.WindVector(2);
            W_vec = repmat([Wx, Wy, 0], num_steps, 1);
            
            % 物理計算
            Va_vec = zeros(num_steps, 3);
            Eul = zeros(num_steps, 3);
            
            g = 9.80665;
            rho0 = 1.225;
            V0_EAS = obj.PhysicsParams.V_air;
            
            for i=1:num_steps
                h = max(0, z(i));
                rho = obj.AtmoModel.get_density(h/1000);
                V_TAS = V0_EAS * sqrt(rho0 / rho);
                
                % 対気速度ベクトル (方向はpsi, 大きさはTAS)
                % ※厳密には経路角gamma成分があるが、水平速度成分主体で近似
                u_air = V_TAS * cos(psi(i));
                v_air = V_TAS * sin(psi(i));
                
                % 沈下速度 w_air の推定
                if i < num_steps
                    dz = z(i+1) - z(i); dt = t_unique(i+1) - t_unique(i);
                    w_air = dz/dt; % 近似
                else
                    w_air = Va_vec(i-1, 3);
                end
                
                Va_vec(i,:) = [u_air, v_air, w_air];
                
                % バンク角 (kappaから復元)
                % tan(phi) = V^2 * kappa / g
                phi_est = atan( (V_TAS^2 * kappa(i)) / g );
                
                Eul(i, :) = [phi_est, 0, psi(i)]; % thetaは簡易的に0
            end
            
            % 対地速度 (Va + Wind)
            Vg_vec = Va_vec + W_vec;
            
            % 位置 (x,y,z は既に対地)
            Pos = [x, y, z];
            
            % テーブル作成
            obj.TrajectoryLogAir = table(t_unique, Pos, Vg_vec, Va_vec, zeros(num_steps,1), W_vec, Eul, ...
                'VariableNames',{'Time','Position','V_Ground','V_Air','Alpha','Wind','Euler_RPY'});
            
            obj.TrajectoryLogGround = obj.TrajectoryLogAir; % 内容は同じ(風考慮済み座標なので)
            
            fprintf('  -> Dual Trajectories Calculated (Steps: %d)\n', num_steps);
        end
        
        % 詳細軌道データのエクスポート
        function trajTable = export_detailed_trajectory(obj, mode)
            if nargin < 2, mode = 'Ground'; end
            if isempty(obj.TrajectoryLogGround), obj.compute_dual_trajectories(); end
            trajTable = obj.TrajectoryLogGround;
        end
        
        % 6DOF用状態量の算出
        function stateTable = export_6dof_state(obj)
            % EXPORT_6DOF_STATE
            % シミュレーション結果から、以下の12変数を計算してテーブルで返す
            % 1:u, 2:v, 3:w (Body速度)
            % 4:p, 5:q, 6:r (Body角速度)
            % 7:phi, 8:theta, 9:psi (Euler角)
            % 10:N, 11:E, 12:D (NED位置)
            
            % 1. 基本データの取得 (Ground Frame)
            T = obj.export_detailed_trajectory('Ground');
            if isempty(T), error('データがありません'); end
            
            time = T.Time;
            
            % --- 位置 (NED) ---
            % T.Position = [North, East, Altitude]
            % D = -Altitude (Depth)
            x = T.Position(:,1);
            y = T.Position(:,2);
            z = -T.Position(:,3);
            
            % --- 姿勢角 (Euler) ---
            % T.Euler_RPY = [phi, theta, psi]
            Phi   = T.Euler_RPY(:,1);
            Theta = T.Euler_RPY(:,2);
            Psi   = T.Euler_RPY(:,3);
            
            % --- 対気速度ベクトル (NED Frame) ---
            % ※u,v,wは「空気に対する」機体の挙動なので、V_Airを使う
            V_Air_NED = T.V_Air; % [Vn, Ve, Vd]
            
            % データ点数
            n_steps = length(time);
            
            % 結果格納用配列
            u = zeros(n_steps, 1);
            v = zeros(n_steps, 1);
            w = zeros(n_steps, 1);
            p = zeros(n_steps, 1);
            q = zeros(n_steps, 1);
            r = zeros(n_steps, 1);
            
            % --- 微分による角速度算出用 ---
            % Euler角の時間微分 (dPhi/dt, dTheta/dt, dPsi/dt)
            dPhi   = gradient(unwrap(Phi), time);
            dTheta = gradient(unwrap(Theta), time);
            dPsi   = gradient(unwrap(Psi), time);
            
            % ループ計算
            for i = 1:n_steps
                % 現在の姿勢
                ph = Phi(i); th = Theta(i); ps = Psi(i);
                
                % ----------------------------------------------------
                % 1. 機体速度 (u, v, w) の計算
                % ----------------------------------------------------
                % NED座標系の速度ベクトルを、機体座標系へ回転させる
                % V_body = R_nb' * V_ned  (R_nb: Body to NED DCM)
                
                cph = cos(ph); sph = sin(ph);
                cth = cos(th); sth = sin(th);
                cps = cos(ps); sps = sin(ps);
                
                % Body -> NED 回転行列 (Z-Y-X順)
                R_nb = [ ...
                    cth*cps,  sph*sth*cps - cph*sps,  cph*sth*cps + sph*sps;
                    cth*sps,  sph*sth*sps + cph*cps,  cph*sth*sps - sph*cps;
                   -sth,      sph*cth,                cph*cth ...
                ];
            
                % NED -> Body (転置)
                R_bn = R_nb';
                
                % 変換実行
                V_ned_vec = V_Air_NED(i, :)'; % 列ベクトル
                V_body_vec = R_bn * V_ned_vec;
                
                u(i) = V_body_vec(1);
                v(i) = V_body_vec(2); % 横滑りなし仮定ならほぼ0になるはず
                w(i) = V_body_vec(3);
                
                % ----------------------------------------------------
                % 2. 機体角速度 (p, q, r) の計算
                % ----------------------------------------------------
                % Euler角速度 -> 機体角速度 の変換式 (Kinematic Equations)
                % p = dphi - dpsi * sin(theta)
                % q = dtheta * cos(phi) + dpsi * cos(theta) * sin(phi)
                % r = -dtheta * sin(phi) + dpsi * cos(theta) * cos(phi)
                
                d_ph = dPhi(i); d_th = dTheta(i); d_ps = dPsi(i);
                
                p(i) = d_ph - d_ps * sth;
                q(i) = d_th * cph + d_ps * cth * sph;
                r(i) = -d_th * sph + d_ps * cth * cph;
            end
            
            % --- テーブル作成 ---
            stateTable = table(u, v, w, p, q, r, Phi, Theta, Psi, x, y, z, ...
                'VariableNames', {'u','v','w', 'p','q','r', 'phi','theta','psi', 'x','y','z'});
            
            % 単位情報の付与 (Properties.VariableUnits)
            stateTable.Properties.VariableUnits = { ...
                'm/s', 'm/s', 'm/s', ... % u, v, w
                'rad/s', 'rad/s', 'rad/s', ... % p, q, r
                'rad', 'rad', 'rad', ... % phi, theta, psi
                'm', 'm', 'm' ... % N, E, D
            };
            
            fprintf('  -> 6-DOF State Variables (12 states) computed.\n');
        end
    end
end