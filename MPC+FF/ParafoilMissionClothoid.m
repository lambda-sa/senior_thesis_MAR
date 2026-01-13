classdef ParafoilMissionClothoid < ParafoilMission
    % PARAFOILMISSIONCLOTHOID (12-State Export Version)
    %
    % 追加機能:
    % - export_6dof_state: 参照軌道として使える12変数を計算して出力
    %   [u, v, w, p, q, r, phi, theta, psi, x, y, z]
    
    properties
        TrajectoryLogAir table     % 対気座標系ログ
        TrajectoryLogGround table  % 対地座標系ログ
    end
    
    methods
        function obj = ParafoilMissionClothoid(excelFileName)
            obj@ParafoilMission(excelFileName);
            obj.Planner = ParafoilPathPlannerClothoid(obj.AtmoModel);
        end
        
        % オプション設定
        function set_clothoid_options(obj, max_roll_rate, look_ahead, run_up_dist)
            if isa(obj.Planner, 'ParafoilPathPlannerClothoid')
                if nargin > 1, obj.Planner.MaxRollRate = max_roll_rate; end
                if nargin > 2, obj.Planner.LookAheadFactor = look_ahead; end
                if nargin > 3, obj.Planner.RunUpDistance = run_up_dist; end
                fprintf('Clothoid Options Updated.\n');
            end
        end

        % 詳細軌道データのエクスポート
        function trajTable = export_detailed_trajectory(obj, mode)
            if nargin < 2, mode = 'Ground'; end
            if isempty(obj.TrajectoryLogGround), obj.compute_dual_trajectories(); end
            if strcmpi(mode, 'Air'), trajTable = obj.TrajectoryLogAir; else, trajTable = obj.TrajectoryLogGround; end
        end

        % =========================================================
        % ★追加: 12変数 (6-DOF) 状態量の計算と出力
        % =========================================================
        function stateTable = export_6dof_state(obj)
            % EXPORT_6DOF_STATE
            % シミュレーション結果(Ground)から、以下の12変数を計算する
            % 1-3: u, v, w (Body Frame Velocity)
            % 4-6: p, q, r (Body Frame Angular Rate)
            % 7-9: phi, theta, psi (Euler Angles)
            % 10-12: x, y, z (NED Position, z is -Altitude)
            
            % 1. ログ取得
            T = obj.export_detailed_trajectory('Ground');
            if isempty(T), error('ログデータがありません。'); end
            
            time = T.Time;
            n_steps = length(time);
            
            % 2. 状態量の抽出
            % 位置 (NED: zは下向き正なので、高度のマイナスを取る)
            x_ned = T.Position(:,1); % North
            y_ned = T.Position(:,2); % East
            z_ned = -T.Position(:,3); % Down (-Altitude)
            
            % 姿勢角
            Phi   = T.Euler_RPY(:,1);
            Theta = T.Euler_RPY(:,2);
            Psi   = T.Euler_RPY(:,3);
            
            % 対気速度ベクトル (NED Frame)
            % ※制御では対気速度ベースの挙動が重要
            V_Air_NED = T.V_Air; 
            
            % 3. 微分によるレート計算
            % gradientを使って数値微分 (dtが一定でない場合も考慮される)
            dPhi   = gradient(unwrap(Phi), time);
            dTheta = gradient(unwrap(Theta), time);
            dPsi   = gradient(unwrap(Psi), time);
            
            % 結果格納用
            u = zeros(n_steps, 1); v = zeros(n_steps, 1); w = zeros(n_steps, 1);
            p = zeros(n_steps, 1); q = zeros(n_steps, 1); r = zeros(n_steps, 1);
            
            for i = 1:n_steps
                ph = Phi(i); th = Theta(i); ps = Psi(i);
                
                % --- (A) 角速度 (p, q, r) の計算 ---
                % Euler角速度 -> 機体角速度 の変換 (Kinematics)
                % p = dphi - dpsi*sin(theta)
                % q = dtheta*cos(phi) + dpsi*cos(theta)*sin(phi)
                % r = -dtheta*sin(phi) + dpsi*cos(theta)*cos(phi)
                
                sth = sin(th); cth = cos(th);
                sph = sin(ph); cph = cos(ph);
                
                p(i) = dPhi(i) - dPsi(i) * sth;
                q(i) = dTheta(i) * cph + dPsi(i) * cth * sph;
                r(i) = -dTheta(i) * sph + dPsi(i) * cth * cph;
                
                % --- (B) 機体速度 (u, v, w) の計算 ---
                % NED速度 -> Body速度 (回転行列の転置)
                % R_nb (Body to NED)
                cps = cos(ps); sps = sin(ps);
                
                R_nb = [ ...
                    cth*cps,  sph*sth*cps - cph*sps,  cph*sth*cps + sph*sps;
                    cth*sps,  sph*sth*sps + cph*cps,  cph*sth*sps - sph*cps;
                   -sth,      sph*cth,                cph*cth ...
                ];
                
                % V_body = R_nb' * V_ned
                V_ned_vec = V_Air_NED(i, :)';
                V_body_vec = R_nb' * V_ned_vec;
                
                u(i) = V_body_vec(1);
                v(i) = V_body_vec(2);
                w(i) = V_body_vec(3);
            end
            
            % 4. テーブル作成
            stateTable = table(u, v, w, p, q, r, Phi, Theta, Psi, x_ned, y_ned, z_ned, ...
                'VariableNames', {'u','v','w', 'p','q','r', 'phi','theta','psi', 'x','y','z'});
                
            stateTable.Properties.VariableUnits = { ...
                'm/s','m/s','m/s', 'rad/s','rad/s','rad/s', 'rad','rad','rad', 'm','m','m'};
                
            fprintf('  -> 6-DOF Reference State (12 vars) Computed.\n');
        end

        % =========================================================
        % 風考慮シミュレーション (前回と同じ)
        % =========================================================
        function run_wind_simulation(obj, target_pos, L_final, wind_vec_2d)
            obj.Planner.apply_wind_effect(wind_vec_2d(1), wind_vec_2d(2));
            obj.compute_physics_parameters();
            wind_up_rad = atan2(-wind_vec_2d(2), -wind_vec_2d(1));
            landing_dir_deg_auto = rad2deg(wind_up_rad);
            
            fprintf('--- Wind Simulation (Clothoid) ---\n');
            target_pos_air = target_pos; 
            for iter = 1:5
                obj.Planner.plan_optimal_trajectory(...
                    obj.PhysicsParams.start_pos, target_pos_air, landing_dir_deg_auto, L_final, ...              
                    obj.PhysicsParams.V_air, obj.PhysicsParams.glide_ratio, obj.PhysicsParams.min_turn_radius);
                d = obj.Planner.ResultData;
                if isempty(d.final.t), break; end
                actual_time = d.final.t(end);
                Drift = wind_vec_2d * actual_time;
                ideal_target = target_pos - Drift;
                err = norm(ideal_target - target_pos_air);
                if err < 1.0, break; end
                target_pos_air = ideal_target;
            end
            obj.Planner.save_naive_path();
            obj.Planner.apply_wind_effect(wind_vec_2d(1), wind_vec_2d(2));
            obj.compute_dual_trajectories();
        end

        % =========================================================
        % ログ生成 (前回と同じ)
        % =========================================================
        function compute_dual_trajectories(obj)
            if isprop(obj.Planner, 'ResultDataNaive') && ~isempty(obj.Planner.ResultDataNaive)
                d_base = obj.Planner.ResultDataNaive;
            else, d_base = obj.Planner.ResultData; end
            if isempty(d_base), error('No Data'); end
            
            t=[]; x=[]; y=[]; z=[]; fields = {'runup', 'entry', 'loiter', 'dubins', 'final'};
            for i=1:length(fields), f=fields{i}; 
                if isfield(d_base, f) && ~isempty(d_base.(f).t)
                    t=[t; d_base.(f).t(:)]; x=[x; d_base.(f).x(:)]; y=[y; d_base.(f).y(:)]; z=[z; d_base.(f).z(:)];
                end
            end
            [t_unique, idx] = unique(t, 'stable'); x=x(idx); y=y(idx); z=z(idx);
            num_steps = length(t_unique); if num_steps==0, warning('No Data'); return; end
            
            [params, sim_settings] = load_params_from_excel(obj.ExcelFileName); atmo = obj.AtmoModel; 
            dyn = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, atmo);
            if isprop(obj.Planner, 'WindVector'), Wx=obj.Planner.WindVector(1); Wy=obj.Planner.WindVector(2); else, Wx=0; Wy=0; end
            if isprop(obj.Planner, 'V0'), V0=obj.Planner.V0; else, V0=15.0; end
            
            dx=gradient(x, t_unique); dy=gradient(y, t_unique); psi_v=atan2(dy,dx); dpsi_v=gradient(unwrap(psi_v),t_unique);
            Va=zeros(num_steps,3); alph=zeros(num_steps,1); phi=zeros(num_steps,1);
            if isfield(sim_settings, 'ang_gamma_rigging'), mu=sim_settings.ang_gamma_rigging; else, mu=0; end
            
            for i=1:num_steps
                h=max(0,z(i)); rho=atmo.get_density(h/1000); V_TAS=V0*sqrt(1.225/rho);
                dz=gradient(z,t_unique); dir=[dx(i),dy(i),dz(i)]; vn=norm(dir);
                if vn>1e-6, Va(i,:)=(dir/vn)*V_TAS; else, Va(i,:)=[V_TAS,0,0]; end
                phi(i)=atan((V_TAS*dpsi_v(i))/9.81);
                pl=params; pl.rho=rho; pl.g=atmo.get_gravity(h);
                av=dyn.solve_trim_alpha(pl,0,0,mu); if isnan(av), if i>1, av=alph(i-1); else, av=0; end; end
                alph(i)=av;
            end
            gam=atan2(Va(:,3), hypot(Va(:,1),Va(:,2))); the=gam+alph; Eul=[phi, the, psi_v];
            posG=[x+Wx*t_unique, y+Wy*t_unique, z];
            obj.TrajectoryLogAir=table(t_unique, [x,y,z], Va, Va, alph, zeros(num_steps,3), Eul, 'VariableNames',{'Time','Position','V_Ground','V_Air','Alpha','Wind','Euler_RPY'});
            obj.TrajectoryLogGround=table(t_unique, posG, Va+[Wx,Wy,0], Va, alph, repmat([Wx,Wy,0],num_steps,1), Eul, 'VariableNames',{'Time','Position','V_Ground','V_Air','Alpha','Wind','Euler_RPY'});
        end
        
        function plot_phases_static(obj), obj.Planner.plot_wind_comparison(); end % 簡易表示
        function plot_wind_comparison(obj), obj.Planner.plot_wind_comparison(); end
    end
end