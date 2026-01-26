classdef ParafoilMissionClothoid < ParafoilMission
    % PARAFOILMISSIONCLOTHOID (Robust Version)
    % ParafoilMissionWind の堅牢なロジックを移植した改良版
    
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

        % =========================================================
        % 風考慮シミュレーション (Wind版と同等のロジックへ更新)
        % =========================================================
        function run_wind_simulation(obj, target_pos, L_final, wind_vec_2d)
            % 1. 風速ベクトルをプランナーに登録
            obj.Planner.WindVector = [wind_vec_2d(1); wind_vec_2d(2); 0];
            
            % 2. 物理パラメータ計算
            obj.compute_physics_parameters();
            
            % 3. 風上方向の自動計算
            wind_up_rad = atan2(-wind_vec_2d(2), -wind_vec_2d(1));
            landing_dir_deg_auto = rad2deg(wind_up_rad);
            fprintf('  [Auto] Landing Direction set to Wind-Up: %.1f deg\n', landing_dir_deg_auto);
            
            fprintf('--- Iterative Wind Correction (Clothoid) ---\n');
            target_pos_air = target_pos; 
            max_iter = 5;
            tol_dist = 1.0;
            
            for iter = 1:max_iter
                % (A) 経路生成 (Air Frame)
                % ※前回の修正で、失敗しても安全な構造体が返ってくるようになっています
                obj.Planner.plan_optimal_trajectory(...
                    obj.PhysicsParams.start_pos, ...
                    target_pos_air, ...
                    landing_dir_deg_auto, ...
                    L_final, ...              
                    obj.PhysicsParams.V_air, ...
                    obj.PhysicsParams.glide_ratio, ...
                    obj.PhysicsParams.min_turn_radius);
                
                d = obj.Planner.ResultData;
                
                % (B) 時間取得 (失敗時は短時間または空)
                if isempty(d.final.t)
                    warning('Path generation failed (Empty Result).');
                    break; 
                end
                actual_time = d.final.t(end);
                
                % (C) ドリフト計算と目標更新
                Drift = [wind_vec_2d(:)', 0] * actual_time; % 3Dベクトル化
                ideal_target = target_pos - Drift(1:3);     % ターゲットを風上にずらす
                
                % (D) 収束判定
                % ※Z成分は無視して水平距離で判定
                err_dist = norm(ideal_target(1:2) - target_pos_air(1:2));
                
                fprintf('  Iter %d: Time=%.1fs, Err=%.1fm\n', iter, actual_time, err_dist);
                
                if err_dist < tol_dist
                    fprintf('  -> Converged!\n');
                    break;
                end
                
                % 更新 (Zは元のまま維持、XYのみ更新)
                target_pos_air(1) = ideal_target(1);
                target_pos_air(2) = ideal_target(2);
                
                % ★安全装置: ループの途中で高度不足で墜落している場合、
                % 無理にターゲットを遠ざけると永久に届かないので、Breakする判断もアリですが
                % ここでは前回のPlanner修正(マイナス高度許容)が入っている前提で回します。
            end
            
            % 4. 最終結果の保存と展開
            obj.Planner.save_naive_path();
            obj.Planner.apply_wind_effect(wind_vec_2d(1), wind_vec_2d(2));
            
            % 5. 詳細ログ生成
            obj.compute_dual_trajectories();
        end

        % =========================================================
        % ログ生成 (Wind版のロジック + Clothoid対応)
        % =========================================================
        function compute_dual_trajectories(obj)
            % ベースデータの取得 (Air Frame)
            if isprop(obj.Planner, 'ResultDataNaive') && ~isempty(obj.Planner.ResultDataNaive)
                d_air = obj.Planner.ResultDataNaive;
            else
                d_air = obj.Planner.ResultData; 
            end
            if isempty(d_air), error('No Data'); end
            
            % データの結合
            t=[]; x=[]; y=[]; z=[]; kappa=[];
            fields = {'runup', 'entry', 'loiter', 'dubins', 'final'};
            
            for i=1:length(fields)
                f=fields{i}; 
                if isfield(d_air, f) && ~isempty(d_air.(f).t)
                    t=[t; d_air.(f).t(:)]; 
                    x=[x; d_air.(f).x(:)]; 
                    y=[y; d_air.(f).y(:)]; 
                    z=[z; d_air.(f).z(:)];
                    
                    % kappa (曲率) の取得 (Planner修正版で追加されたフィールド)
                    if isfield(d_air.(f), 'kappa') && ~isempty(d_air.(f).kappa)
                        kappa = [kappa; d_air.(f).kappa(:)];
                    else
                        % kappaがない場合のフォールバック（ゼロ埋め）
                        kappa = [kappa; zeros(length(d_air.(f).t), 1)];
                    end
                end
            end
            
            [t_unique, idx] = unique(t, 'stable'); 
            if isempty(t_unique), warning('No valid trajectory steps.'); return; end
            
            x=x(idx); y=y(idx); z=z(idx); kappa=kappa(idx);
            num_steps = length(t_unique);
            
            % パラメータ準備
            [params, sim_settings] = load_params_from_excel(obj.ExcelFileName); 
            atmo = obj.AtmoModel; 
            dyn = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, atmo);
            
            Wx = obj.Planner.WindVector(1); 
            Wy = obj.Planner.WindVector(2);
            W_vec_ground = repmat([Wx, Wy, 0], num_steps, 1);
            
            if isprop(obj.Planner, 'V0'), V0 = obj.Planner.V0; else, V0=15.0; end
            
            % 配列初期化
            Va_vec = zeros(num_steps, 3);
            alph_vec = zeros(num_steps, 1);
            phi_vec = zeros(num_steps, 1);
            
            % 設定値
            mu = 0;
            if isfield(sim_settings, 'ang_gamma_rigging'), mu=sim_settings.ang_gamma_rigging; end
            
            g_std = 9.80665;
            rho_0 = 1.225;
            
            % --- 物理量計算ループ ---
            for i=1:num_steps
                % 高度と密度（マイナス高度対策）
                h = max(0, z(i)); 
                rho = atmo.get_density(h/1000);
                
                % 真対気速度 (TAS)
                V_TAS = V0 * sqrt(rho_0 / rho);
                
                % 1. 対気速度ベクトル (幾何形状から算出)
                % 数値微分を使用するが、方向ベクトルとしてのみ使用
                if i < num_steps
                    dx = x(i+1)-x(i); dy = y(i+1)-y(i); dz_step = z(i+1)-z(i);
                    dist = sqrt(dx^2 + dy^2 + dz_step^2);
                else
                    dx = x(i)-x(i-1); dy = y(i)-y(i-1); dz_step = z(i)-z(i-1);
                    dist = sqrt(dx^2 + dy^2 + dz_step^2);
                end
                
                if dist > 1e-6
                    dir = [dx, dy, dz_step] / dist;
                else
                    dir = [1, 0, 0];
                end
                Va_vec(i, :) = dir * V_TAS;
                
                % 2. バンク角 (phi) の計算
                % ★重要: 数値微分(dpsi)ではなく、Plannerが計画した曲率(kappa)から計算する
                % 公式: tan(phi) = (V^2 * kappa) / g
                k_curr = kappa(i);
                phi_vec(i) = atan( (V_TAS^2 * k_curr) / g_std );
                
                % 3. 迎角 (alpha) の計算 (トリム)
                pl = params; pl.rho = rho; pl.g = atmo.get_gravity(h);
                av = dyn.solve_trim_alpha(pl, 0, 0, mu); 
                if isnan(av)
                    if i>1, av=alph_vec(i-1); else, av=0; end
                end
                alph_vec(i) = av;
            end
            
            % 姿勢角の統合
            % Gamma (経路角)
            gam = atan2(Va_vec(:,3), hypot(Va_vec(:,1), Va_vec(:,2)));
            
            % Theta (ピッチ角) = Gamma + Alpha
            the = gam + alph_vec;
            
            % Psi (方位角)
            % データから計算（数値微分よりアークタンジェントの方が安定）
            % ただし、kappa積分で求めた psi が d_air にあるならそれを使いたい
            psi_vec = zeros(num_steps, 1);
            current_psi_idx = 1;
            for k=1:length(fields)
                f=fields{k};
                if isfield(d_air, f) && ~isempty(d_air.(f).psi)
                    len = length(d_air.(f).psi);
                    psi_vec(current_psi_idx : current_psi_idx+len-1) = d_air.(f).psi;
                    current_psi_idx = current_psi_idx + len;
                end
            end
            % 重複削除の影響でサイズが合わない場合の保険
            if length(psi_vec) > num_steps, psi_vec = psi_vec(idx); end
            if length(psi_vec) < num_steps, psi_vec = [psi_vec; repmat(psi_vec(end), num_steps-length(psi_vec), 1)]; end
            
            Eul = [phi_vec, the, psi_vec];
            
            % 対地位置 (Air位置 + 風オフセット)
            posG = [x + Wx*t_unique, y + Wy*t_unique, z];
            Vg_vec = Va_vec + W_vec_ground;
            
            % テーブル作成
            obj.TrajectoryLogAir = table(t_unique, [x,y,z], Va_vec, Va_vec, alph_vec, zeros(num_steps,3), Eul, ...
                'VariableNames',{'Time','Position','V_Ground','V_Air','Alpha','Wind','Euler_RPY'});
            
            obj.TrajectoryLogGround = table(t_unique, posG, Vg_vec, Va_vec, alph_vec, W_vec_ground, Eul, ...
                'VariableNames',{'Time','Position','V_Ground','V_Air','Alpha','Wind','Euler_RPY'});
            
            fprintf('  -> Dual Trajectories (Clothoid) Calculated.\n');
        end
        
        % 詳細軌道データのエクスポート
        function trajTable = export_detailed_trajectory(obj, mode)
            if nargin < 2, mode = 'Ground'; end
            if isempty(obj.TrajectoryLogGround), obj.compute_dual_trajectories(); end
            if strcmpi(mode, 'Air'), trajTable = obj.TrajectoryLogAir; else, trajTable = obj.TrajectoryLogGround; end
        end

        % =========================================================
        % 12変数 (6-DOF) 状態量の計算と出力 (Wind版と同じロジックでOK)
        % =========================================================
        function stateTable = export_6dof_state(obj)
            % 1. ログ取得
            T = obj.export_detailed_trajectory('Ground');
            if isempty(T), error('ログデータがありません。'); end
            
            time = T.Time; n_steps = length(time);
            x_ned = T.Position(:,1); y_ned = T.Position(:,2); z_ned = -T.Position(:,3);
            Phi = T.Euler_RPY(:,1); Theta = T.Euler_RPY(:,2); Psi = T.Euler_RPY(:,3);
            V_Air_NED = T.V_Air;
            
            % 微分
            dPhi = gradient(unwrap(Phi), time);
            dTheta = gradient(unwrap(Theta), time);
            dPsi = gradient(unwrap(Psi), time);
            
            u=zeros(n_steps,1); v=zeros(n_steps,1); w=zeros(n_steps,1);
            p=zeros(n_steps,1); q=zeros(n_steps,1); r=zeros(n_steps,1);
            
            for i = 1:n_steps
                ph=Phi(i); th=Theta(i); ps=Psi(i);
                sth=sin(th); cth=cos(th); sph=sin(ph); cph=cos(ph); cps=cos(ps); sps=sin(ps);
                
                % p, q, r
                p(i) = dPhi(i) - dPsi(i) * sth;
                q(i) = dTheta(i) * cph + dPsi(i) * cth * sph;
                r(i) = -dTheta(i) * sph + dPsi(i) * cth * cph;
                
                % u, v, w
                R_nb = [cth*cps, sph*sth*cps-cph*sps, cph*sth*cps+sph*sps;
                        cth*sps, sph*sth*sps+cph*cps, cph*sth*sps-sph*cps;
                       -sth,     sph*cth,             cph*cth];
                V_body = R_nb' * V_Air_NED(i, :)';
                u(i)=V_body(1); v(i)=V_body(2); w(i)=V_body(3);
            end
            
            stateTable = table(u, v, w, p, q, r, Phi, Theta, Psi, x_ned, y_ned, z_ned, ...
                'VariableNames', {'u','v','w', 'p','q','r', 'phi','theta','psi', 'x','y','z'});
        end
        
        % 可視化用ショートカット
        function plot_wind_comparison(obj), obj.Planner.plot_wind_comparison(); end
    end
end