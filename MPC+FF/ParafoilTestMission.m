classdef ParafoilTestMission < ParafoilMissionWind
    % PARAFOILTESTMISSION: ParafoilMissionWindの12変数変換機能を継承
    
    methods
        function obj = ParafoilTestMission(excelFileName)
            obj@ParafoilMissionWind(excelFileName);
        end
        
        function traj12 = run_test_maneuver(obj, duration, bank_deg)
            % 1. 基本パラメータ計算 (親クラスのメソッド)
            obj.compute_physics_parameters();
            pp = obj.PhysicsParams;
            
            % 2. テスト用の幾何軌道を生成 (x, y, z, psi, t)
            dt = 0.05;
            t_vec = (0:dt:duration)';
            N = length(t_vec);
            
            % トリム情報
            alpha = pp.alpha_trim;
            gamma0 = -atan(1/pp.glide_ratio); 
            phi_ref = deg2rad(bank_deg);
            
            % 積分準備
            curr_pos = pp.start_pos(1:3); % [N, E, Alt]
            psi_ref = pp.start_pos(4);
            rho0 = 1.225;
            
            x = zeros(N,1); y = zeros(N,1); z = zeros(N,1); psi = zeros(N,1);
            
            for i = 1:N
                h = max(0, curr_pos(3));
                rho = obj.AtmoModel.get_density(h/1000);
                V_tas = pp.V_air * sqrt(rho0 / rho);
                
                gamma = atan(tan(gamma0) / cos(phi_ref));
                x(i) = curr_pos(1); y(i) = curr_pos(2); z(i) = -curr_pos(3);
                psi(i) = psi_ref;
                
                % 位置更新
                V_horiz = V_tas * cos(gamma);
                curr_pos(1) = curr_pos(1) + V_horiz * cos(psi_ref) * dt;
                curr_pos(2) = curr_pos(2) + V_horiz * sin(psi_ref) * dt;
                curr_pos(3) = curr_pos(3) + (V_tas * sin(gamma)) * dt;
                
                % 旋回率
                if abs(bank_deg) > 0.1
                    g = obj.AtmoModel.get_gravity(h);
                    psi_dot = (g / V_tas) * tan(phi_ref);
                    psi_ref = psi_ref + psi_dot * dt;
                end
            end
            
            % 3. 親クラスの計算機（Planner）にデータを無理やり詰め込む
            % これにより、親クラスの compute_dual_trajectories が動くようになる
            d = struct();
            d.loiter.x=[]; d.loiter.y=[]; d.loiter.z=[]; d.loiter.t=[];
            d.dubins.x=[]; d.dubins.y=[]; d.dubins.z=[]; d.dubins.t=[];
            d.final.x = x'; d.final.y = y'; d.final.z = z'; d.final.t = t_vec';
            d.final.psi = psi';
            obj.Planner.ResultData = d;
            obj.Planner.ResultDataWind = d; % 風なしテスト用
            
            % 4. ★親クラスの「12変数変換機能」を呼び出す
            obj.compute_dual_trajectories();
            
            % 5. ★親クラスの「12変数専用エクスポート」で表を出力
            % これが u, v, w, p, q, r を全て含んだ決定版の表です
            traj12 = obj.export_6dof_state();
        end
    end
end