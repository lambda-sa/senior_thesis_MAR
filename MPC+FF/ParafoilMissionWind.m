classdef ParafoilMissionWind < ParafoilMission
    % PARAFOILMISSIONWIND (決定版: 風上着陸 & 反復補正)
    
    properties
        TrajectoryLogAir table     % 対気座標系
        TrajectoryLogGround table  % 対地座標系
    end
    
    methods
        function obj = ParafoilMissionWind(excelFileName)
            obj@ParafoilMission(excelFileName);
            obj.Planner = ParafoilPathPlannerWind(obj.AtmoModel);
        end
        
        function set_turn_input_ratio(obj, ratio)
            obj.Planner.TurnInputRatio = ratio;
        end
        
        function set_minimum_orbit(obj, num_turns)
            obj.Planner.set_min_loiter_turn(num_turns);
        end

        function run_wind_simulation(obj, target_pos, L_final, wind_vec_2d)
            % 1. 風速ベクトル保存
            obj.Planner.WindVector = [wind_vec_2d(1); wind_vec_2d(2); 0];
            
            % 2. 物理パラメータ計算
            obj.compute_physics_parameters();
            
            % ========================================================
            % ★追加: 着陸方向を「風上」に自動設定
            % ========================================================
            % 風ベクトル (North, East) の逆方向を計算
            % atan2(Y, X) -> atan2(East, North) で方位角(北基準x軸からの偏角)を取得
            wind_up_rad = atan2(-wind_vec_2d(2), -wind_vec_2d(1));
            landing_dir_deg_auto = rad2deg(wind_up_rad);
            
            fprintf('  [Auto] Landing Direction set to Wind-Up: %.1f deg\n', landing_dir_deg_auto);
            
            % --------------------------------------------------------
            % 3. 目標点の風補正 (Iterative Wind Correction)
            % --------------------------------------------------------
            
            target_pos_air = target_pos; 
            max_iter = 5;
            tol_dist = 1.0;
            
            fprintf('--- Iterative Wind Correction Start ---\n');
            
            for iter = 1:max_iter
                % (A) 経路生成
                % ★ここで Excelの設定値ではなく、計算した「風上方位」を渡す
                obj.Planner.plan_optimal_trajectory(obj.PhysicsParams.start_pos, ...
                                                    target_pos_air, ...
                                                    landing_dir_deg_auto, ... % <--- 自動設定した方位
                                                    L_final, ...              % <--- 定常滑空区間の長さ
                                                    obj.PhysicsParams.V_air, ...
                                                    obj.PhysicsParams.glide_ratio, ...
                                                    obj.PhysicsParams.min_turn_radius);
                
                % (B) 時間取得
                d = obj.Planner.ResultData;
                if isempty(d.final.t), warning('経路生成失敗'); break; end
                actual_time = d.final.t(end);
                
                % (C) ドリフト計算
                Drift_X = wind_vec_2d(1) * actual_time;
                Drift_Y = wind_vec_2d(2) * actual_time;
                
                % (D) 理想の対気目標点
                ideal_target_air_X = target_pos(1) - Drift_X;
                ideal_target_air_Y = target_pos(2) - Drift_Y;
                
                % (E) 収束判定
                diff_X = ideal_target_air_X - target_pos_air(1);
                diff_Y = ideal_target_air_Y - target_pos_air(2);
                err_dist = sqrt(diff_X^2 + diff_Y^2);
                
                fprintf('  Iter %d: Time=%.1fs, Drift=[%.1f, %.1f], Err=%.1fm\n', ...
                    iter, actual_time, Drift_X, Drift_Y, err_dist);
                
                if err_dist < tol_dist
                    fprintf('  -> Converged!\n');
                    break;
                end
                
                % (F) 更新
                target_pos_air(1) = ideal_target_air_X;
                target_pos_air(2) = ideal_target_air_Y;
            end
            
            % --------------------------------------------------------
            
            % 4. 最終保存
            obj.Planner.save_naive_path();
            obj.Planner.apply_wind_effect(wind_vec_2d(1), wind_vec_2d(2));
            
            fprintf('  -> Path Planning Complete. (Final Err: %.2f m)\n', obj.Planner.FinalError);
            
            % 5. 軌道詳細計算
            obj.compute_dual_trajectories();
        end
        
        function compute_dual_trajectories(obj)
            d_air = obj.Planner.ResultData;
            if isempty(d_air), error('経路データがありません。'); end
            
            t = [d_air.loiter.t(:); d_air.dubins.t(:); d_air.final.t(:)];
            x_air = [d_air.loiter.x(:); d_air.dubins.x(:); d_air.final.x(:)];
            y_air = [d_air.loiter.y(:); d_air.dubins.y(:); d_air.final.y(:)];
            z     = [d_air.loiter.z(:); d_air.dubins.z(:); d_air.final.z(:)];
            
            [t_unique, idx] = unique(t, 'stable');
            x_air = x_air(idx); y_air = y_air(idx); z = z(idx);
            num_steps = length(t_unique);
            
            [params, sim_settings] = load_params_from_excel(obj.ExcelFileName); 
            atmo = obj.AtmoModel; 
            dyn_calculator = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, atmo);
            
            Wx = obj.Planner.WindVector(1); 
            Wy = obj.Planner.WindVector(2);
            W_vec_ground = repmat([Wx, Wy, 0], num_steps, 1);
            
            if isprop(obj.Planner, 'V0'), V_EAS = obj.Planner.V0; else, V_EAS = 15.0; end
            if isprop(obj.Planner, 'R_fixed'), R_turn = obj.Planner.R_fixed; else, R_turn = 50.0; end
            g = 9.81; rho_0 = 1.225;
            
            dx_air = gradient(x_air, t_unique); 
            dy_air = gradient(y_air, t_unique);
            psi_air = atan2(dy_air, dx_air);
            dpsi_air = gradient(unwrap(psi_air), t_unique);
            
            Va_vec = zeros(num_steps, 3);
            alpha_vec = zeros(num_steps, 1);
            phi_vec = zeros(num_steps, 1);
            
            delta_a_in = 0; delta_s_in = 0; 
            if isfield(sim_settings, 'ang_gamma_rigging'), mu_in = sim_settings.ang_gamma_rigging;
            elseif isfield(params, 'ang_psi_rigging'), mu_in = params.ang_psi_rigging;
            else, mu_in = 0; end

            for i = 1:num_steps
                h_curr = max(0, z(i));
                rho_curr = atmo.get_density(h_curr / 1000);
                
                V_TAS = V_EAS * sqrt(rho_0 / rho_curr);
                
                dz_geom = gradient(z, t_unique);
                dir_geom = [dx_air(i), dy_air(i), dz_geom(i)];
                dir_norm = norm(dir_geom);
                if dir_norm > 1e-6, dir_vec = dir_geom / dir_norm; else, dir_vec = [1,0,0]; end
                
                Va_vec(i, :) = dir_vec * V_TAS;
                
                phi_ideal = atan(V_TAS^2 / (g * R_turn));
                if abs(dpsi_air(i)) > 0.01
                    phi_vec(i) = sign(dpsi_air(i)) * phi_ideal;
                else
                    phi_vec(i) = 0;
                end
                
                p_loc = params; p_loc.rho = rho_curr; p_loc.g = atmo.get_gravity(h_curr);
                alpha_val = dyn_calculator.solve_trim_alpha(p_loc, delta_a_in, delta_s_in, mu_in);
                if isnan(alpha_val), if i>1, alpha_val=alpha_vec(i-1); else, alpha_val=0; end; end
                alpha_vec(i) = alpha_val;
            end
            
            gamma_vec = atan2(Va_vec(:,3), sqrt(Va_vec(:,1).^2 + Va_vec(:,2).^2));
            theta_vec = gamma_vec + alpha_vec;
            Euler_Air = [phi_vec, theta_vec, psi_air];
            Euler_Ground = [phi_vec, theta_vec, psi_air];
            
            x_g = x_air + Wx * t_unique;
            y_g = y_air + Wy * t_unique;
            Vg_vec = Va_vec + W_vec_ground;
            
            obj.TrajectoryLogAir = table(t_unique, [x_air, y_air, z], Va_vec, Va_vec, alpha_vec, zeros(num_steps, 3), Euler_Air, ...
                'VariableNames', {'Time', 'Position', 'V_Ground', 'V_Air', 'Alpha', 'Wind', 'Euler_RPY'});
            obj.TrajectoryLogAir.Properties.VariableUnits = {'s', 'm', 'm/s', 'm/s', 'rad', 'm/s', 'rad'};
            
            obj.TrajectoryLogGround = table(t_unique, [x_g, y_g, z], Vg_vec, Va_vec, alpha_vec, W_vec_ground, Euler_Ground, ...
                'VariableNames', {'Time', 'Position', 'V_Ground', 'V_Air', 'Alpha', 'Wind', 'Euler_RPY'});
            obj.TrajectoryLogGround.Properties.VariableUnits = {'s', 'm', 'm/s', 'm/s', 'rad', 'm/s', 'rad'};
            
            fprintf('  -> Dual Trajectories (Air/Ground) calculated and cached.\n');
        end
        
        function trajTable = export_detailed_trajectory(obj, mode)
            if nargin < 2, mode = 'Ground'; end
            if isempty(obj.TrajectoryLogGround)
                obj.compute_dual_trajectories();
            end
            if strcmpi(mode, 'Air')
                trajTable = obj.TrajectoryLogAir;
            else
                trajTable = obj.TrajectoryLogGround;
            end
        end
        
        function air_data = export_air_path_data(obj)
            d = obj.Planner.ResultData;
            if isempty(d), air_data = []; return; end
            t = [d.loiter.t(:); d.dubins.t(:); d.final.t(:)];
            x = [d.loiter.x(:); d.dubins.x(:); d.final.x(:)];
            y = [d.loiter.y(:); d.dubins.y(:); d.final.y(:)];
            z = [d.loiter.z(:); d.dubins.z(:); d.final.z(:)];
            [t_unique, idx] = unique(t, 'stable');
            air_data.time = t_unique;
            air_data.position = [x(idx), y(idx), z(idx)];
        end
        
        function plot_wind_comparison(obj)
            if isempty(obj.Planner.ResultDataWind), warning('No result data.'); return; end
            obj.Planner.plot_wind_comparison();
        end
        
        function animate_trajectory(obj, speed_factor, mode)
            if nargin < 2, speed_factor = 10.0; end
            if nargin < 3, mode = 'Ground'; end
            T = obj.export_detailed_trajectory(mode);
            t = T.Time; pos = T.Position; eul = T.Euler_RPY;
            figure('Name', sprintf('Trajectory Animation (%s Frame)', mode), 'Color', 'w');
            plot3(pos(:,2), pos(:,1), pos(:,3), 'Color', [0.7 0.7 0.7], 'LineStyle', '--');
            hold on; grid on; axis equal; view(3);
            if strcmpi(mode, 'Ground'), plot3(0,0,0, 'gx', 'MarkerSize',15); end
            xlabel('East [m]'); ylabel('North [m]'); zlabel('Altitude [m]');
            h_pt = plot3(0,0,0, 'bo', 'MarkerFaceColor','b', 'MarkerSize', 6);
            h_vec = plot3([0 0],[0 0],[0 0], 'r-', 'LineWidth', 2);
            h_title = title('');
            xlim([min(pos(:,2))-100, max(pos(:,2))+100]); ylim([min(pos(:,1))-100, max(pos(:,1))+100]); zlim([min(pos(:,3)), max(pos(:,3))+100]);
            tic;
            while true
                t_curr = toc * speed_factor;
                if t_curr > t(end), break; end
                idx = find(t >= t_curr, 1); if isempty(idx), idx = length(t); end
                P = pos(idx,:); E = eul(idx,:);
                set(h_pt, 'XData', P(2), 'YData', P(1), 'ZData', P(3));
                len = 50; 
                u = len * cos(E(3)) * cos(E(2)); v = len * sin(E(3)) * cos(E(2)); w = len * sin(E(2));
                set(h_vec, 'XData', [P(2), P(2)+v], 'YData', [P(1), P(1)+u], 'ZData', [P(3), P(3)+w]);
                set(h_title, 'String', sprintf('[%s] T: %.1fs / Alt: %.0fm', mode, t_curr, P(3)));
                drawnow limitrate;
            end
        end

        % ============================================================
        %  対地・対気 同時アニメーション
        % ============================================================
        function animate_dual_views(obj, speed_factor)
            if nargin < 2, speed_factor = 10.0; end
            TG = obj.export_detailed_trajectory('Ground');
            TA = obj.export_detailed_trajectory('Air');
            if isempty(TG) || isempty(TA), error('データなし'); end
            
            time = TG.Time; t_final = time(end);
            posG = TG.Position; posA = TA.Position; eul = TG.Euler_RPY;
            
            fig = figure('Name', 'Dual View Animation', 'Color', 'w', ...
                         'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.6]);
            
            pad = 50;
            all_pos = [posG; posA];
            xlim_common = [min(all_pos(:,2))-pad, max(all_pos(:,2))+pad];
            ylim_common = [min(all_pos(:,1))-pad, max(all_pos(:,1))+pad];
            zlim_common = [min(0), max(all_pos(:,3))+pad];
            vec_len = 100; 

            % Ground
            axG = subplot(1, 2, 1); hold on; grid on; axis equal; view(3);
            title('Ground Frame (Actual)');
            xlabel('East [m]'); ylabel('North [m]'); zlabel('Alt [m]');
            xlim(xlim_common); ylim(ylim_common); zlim(zlim_common);
            plot3(axG, posG(:,2), posG(:,1), posG(:,3), 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
            plot3(axG, 0, 0, 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
            hPointG = plot3(axG, nan, nan, nan, 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8);
            hHeadG  = plot3(axG, nan, nan, nan, 'r-', 'LineWidth', 2);

            % Air
            axA = subplot(1, 2, 2); hold on; grid on; axis equal; view(3);
            title('Air Frame (Aerodynamic)');
            xlabel('East [m]'); ylabel('North [m]'); zlabel('Alt [m]');
            xlim(xlim_common); ylim(ylim_common); zlim(zlim_common);
            plot3(axA, posA(:,2), posA(:,1), posA(:,3), 'Color', [0.7 0.7 0.7], 'LineStyle', ':');
            plot3(axA, posA(end,2), posA(end,1), posA(end,3), 'mo', 'MarkerSize', 10, 'LineWidth', 2); 
            hPointA = plot3(axA, nan, nan, nan, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
            hHeadA  = plot3(axA, nan, nan, nan, 'r-', 'LineWidth', 2);
            
            sgtitle(sprintf('Simulation Time: 0.0s / Speed: x%.1f', speed_factor));

            tic;
            while ishandle(fig)
                t_curr = toc * speed_factor;
                if t_curr > t_final, break; end
                idx = find(time >= t_curr, 1); if isempty(idx), idx = length(time); end
                
                PG = posG(idx, :); PA = posA(idx, :); E = eul(idx, :);
                u = vec_len * cos(E(2)) * cos(E(3)); 
                v = vec_len * cos(E(2)) * sin(E(3)); 
                w = vec_len * -sin(E(2));
                
                set(hPointG, 'XData', PG(2), 'YData', PG(1), 'ZData', PG(3));
                set(hHeadG,  'XData', [PG(2), PG(2)+v], 'YData', [PG(1), PG(1)+u], 'ZData', [PG(3), PG(3)-w]);
                
                set(hPointA, 'XData', PA(2), 'YData', PA(1), 'ZData', PA(3));
                set(hHeadA,  'XData', [PA(2), PA(2)+v], 'YData', [PA(1), PA(1)+u], 'ZData', [PA(3), PA(3)-w]);
                
                sgtitle(sprintf('Time: %.1fs / Alt: %.0fm / Bank: %.1f\circ', t_curr, PG(3), rad2deg(E(1))));
                drawnow limitrate;
            end
        end

        % ============================================================
        %  ★修正版: フェーズ別 静止画プロット (Ground Frame)
        % ============================================================
        function plot_phases_static(obj)
            d = obj.Planner.ResultDataWind; % 風適用後のデータ
            if isempty(d), error('データがありません。シミュレーションを実行してください。'); end
            
            % 全体範囲の取得（軸スケールを統一するため）
            all_x = [d.loiter.x, d.dubins.x, d.final.x];
            all_y = [d.loiter.y, d.dubins.y, d.final.y];
            lim_x = [min(all_y)-50, max(all_y)+50]; % East
            lim_y = [min(all_x)-50, max(all_x)+50]; % North

            figure('Name', 'Flight Phases (Ground Frame)', 'Color', 'w', 'Position', [100, 100, 1000, 800]);

            % --- 1. Loiter Phase (待機旋回) ---
            subplot(2, 2, 1);
            grid on; axis equal; hold on; % 設定を直接記述
            xlabel('East [m]'); ylabel('North [m]');
            title('Phase 1: Loiter (Orbit)', 'FontSize', 12, 'FontWeight', 'bold');
            plot(0, 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2); % Target
            
            % 背景に全体を薄く表示
            plot(all_y, all_x, 'Color', [0.9 0.9 0.9], 'LineWidth', 1);
            % メイン描画
            plot(d.loiter.y, d.loiter.x, 'b-', 'LineWidth', 2);
            % 始点・終点
            if ~isempty(d.loiter.x)
                plot(d.loiter.y(1), d.loiter.x(1), 'bo'); 
                plot(d.loiter.y(end), d.loiter.x(end), 'b^');
            end
            xlim(lim_x); ylim(lim_y);

            % --- 2. Dubins Phase (遷移経路) ---
            subplot(2, 2, 2);
            grid on; axis equal; hold on;
            xlabel('East [m]'); ylabel('North [m]');
            title('Phase 2: Dubins (Transition)', 'FontSize', 12, 'FontWeight', 'bold');
            plot(0, 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
            
            plot(all_y, all_x, 'Color', [0.9 0.9 0.9], 'LineWidth', 1);
            plot(d.dubins.y, d.dubins.x, 'r-', 'LineWidth', 2);
            if ~isempty(d.dubins.x)
                plot(d.dubins.y(1), d.dubins.x(1), 'ro');
                plot(d.dubins.y(end), d.dubins.x(end), 'r^');
            end
            xlim(lim_x); ylim(lim_y);

            % --- 3. Final Phase (着陸進入) ---
            subplot(2, 2, 3);
            grid on; axis equal; hold on;
            xlabel('East [m]'); ylabel('North [m]');
            title('Phase 3: Final Approach', 'FontSize', 12, 'FontWeight', 'bold');
            plot(0, 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
            
            plot(all_y, all_x, 'Color', [0.9 0.9 0.9], 'LineWidth', 1);
            plot(d.final.y, d.final.x, 'm-', 'LineWidth', 2);
            if ~isempty(d.final.x)
                plot(d.final.y(1), d.final.x(1), 'mo');
                plot(d.final.y(end), d.final.x(end), 'm^'); % 着地点
            end
            xlim(lim_x); ylim(lim_y);

            % --- 4. Total View (全体像) ---
            subplot(2, 2, 4);
            grid on; axis equal; hold on;
            xlabel('East [m]'); ylabel('North [m]');
            title('Total Trajectory', 'FontSize', 12, 'FontWeight', 'bold');
            plot(0, 0, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
            
            % 各色で結合
            h1 = plot(d.loiter.y, d.loiter.x, 'b-', 'LineWidth', 1.5);
            h2 = plot(d.dubins.y, d.dubins.x, 'r-', 'LineWidth', 1.5);
            h3 = plot(d.final.y, d.final.x, 'm-', 'LineWidth', 1.5);
            
            % 風向き矢印の表示 (右上に配置)
            wx = obj.Planner.WindVector(1); wy = obj.Planner.WindVector(2);
            if norm([wx, wy]) > 0
                arrow_center = [lim_x(2)-100, lim_y(2)-100];
                quiver(arrow_center(1), arrow_center(2), wy*20, wx*20, ...
                    'k', 'LineWidth', 2, 'MaxHeadSize', 0.5, 'AutoScale','off');
                text(arrow_center(1), arrow_center(2)-30, sprintf('Wind\n(%.1fm/s)', norm([wx,wy])), 'HorizontalAlignment','center');
            end
            
            legend([h1, h2, h3], {'Loiter', 'Dubins', 'Final'}, 'Location', 'best');
            xlim(lim_x); ylim(lim_y);
        end
    end
end