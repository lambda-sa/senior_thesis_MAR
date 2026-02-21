classdef ParafoilPathPlanner < handle
    % PARAFOILPATHPLANNER (Original Logic Restored + Mid RunUp Added)
    
    properties
        AtmoModel
        ResultData      % .x, .y, .z, .t
        FinalError
        StartV
        
        V0 = 15.0; 
        R_fixed = 50.0;
        MinLoiterAngle = 0;
        
        % ★追加: 助走距離パラメータ
        RunUpDistance = 0.0;     % 初期助走
        MidRunUpDistance = 130.0; % 中間助走 (Loiter -> Dubinsの操舵反転バッファ)
        
        WindVector = [0; 0; 0];
        EnableBankSpeedCorrection = true;
    end
    
    methods
        function set_speed_correction(obj, enable_flag)
            obj.EnableBankSpeedCorrection = enable_flag;
        end
        function obj = ParafoilPathPlanner(atmo_model)
            obj.AtmoModel = atmo_model;
        end
        function set_min_loiter_turn(obj, num_turns)
            obj.MinLoiterAngle = num_turns * 2 * pi;
        end
        
        function chosen_dir = plan_optimal_trajectory(obj, start_3d, target_3d, land_deg, L_final, V0, glide_ratio, R_fixed)
            obj.V0 = V0;
            obj.R_fixed = R_fixed;
            fprintf('--- Path Planning Optimization (R=%.1fm) ---\n', R_fixed);
            
            [resLL, errLL, lenLL] = obj.try_combination('L', 'L', start_3d, target_3d, land_deg, L_final, V0, glide_ratio, R_fixed);
            [resRR, errRR, lenRR] = obj.try_combination('R', 'R', start_3d, target_3d, land_deg, L_final, V0, glide_ratio, R_fixed);
            
            validLL = abs(errLL) < 10; validRR = abs(errRR) < 10;
            if ~validLL && ~validRR
                [resLR, errLR, lenLR] = obj.try_combination('L', 'R', start_3d, target_3d, land_deg, L_final, V0, glide_ratio, R_fixed);
                [resRL, errRL, lenRL] = obj.try_combination('R', 'L', start_3d, target_3d, land_deg, L_final, V0, glide_ratio, R_fixed);
                candidates = [struct('res',resLL,'err',errLL,'len',lenLL,'valid',validLL), ...
                              struct('res',resRR,'err',errRR,'len',lenRR,'valid',validRR), ...
                              struct('res',resLR,'err',errLR,'len',lenLR,'valid',abs(errLR)<10), ...
                              struct('res',resRL,'err',errRL,'len',lenRL,'valid',abs(errRL)<10)];
            else
                candidates = [struct('res',resLL,'err',errLL,'len',lenLL,'valid',validLL), ...
                              struct('res',resRR,'err',errRR,'len',lenRR,'valid',validRR)];
            end
            
            valid_idx = find([candidates.valid]);
            if isempty(valid_idx)
                [~, best_idx] = min(abs([candidates.err]));
            else
                [~, min_len_idx] = min([candidates(valid_idx).len]);
                best_idx = valid_idx(min_len_idx);
            end
            
            best = candidates(best_idx);
            obj.ResultData = best.res.data;
            obj.FinalError = best.err;
            obj.StartV = best.res.data.V_start;
            chosen_dir = best.res.l_dir;
            fprintf('  >>> Selected: %s-Loiter (DubinsLen: %.1f, Err: %.1f)\n', chosen_dir, best.len, best.err);
        end
        
        function [res, final_err, dubins_len] = try_combination(obj, l_dir, d_dir, start_3d, target_3d, land_deg, Lf, V0, gr, R)
            h0 = start_3d(3); step = 0.5; target_z = target_3d(3); tan_g_str = -1.0 / gr;
            
            % ★修正: 初期助走＋中間助走の両方の高度ロスを考慮する
            z_run_loss = (obj.RunUpDistance + obj.MidRunUpDistance) * abs(tan_g_str);
            cost_func = @(alpha) obj.calc_core(alpha, start_3d, target_3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir) - target_z;
            
            best_alpha = 0;
            try
                min_a = obj.MinLoiterAngle;
                z_at_min = obj.calc_core(min_a, start_3d, target_3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir);
                if z_at_min >= target_z
                    h_avg = (h0 - z_run_loss + target_z)/2; % ここにもz_run_lossを適用
                    g = obj.AtmoModel.get_gravity(h_avg);
                    rho0 = obj.AtmoModel.get_density(h0/1000);
                    rho_avg = obj.AtmoModel.get_density(h_avg/1000);
                    V_avg = V0 * sqrt(rho0 / rho_avg);
                    sigma_avg = atan(V_avg^2 / (g * R));
                    tan_g_avg = tan_g_str / cos(sigma_avg);
                    guess_add = (z_at_min - target_z) / (R * abs(tan_g_avg));
                    guess = min_a + guess_add;
                    best_alpha = fzero(cost_func, [min_a, max(min_a + 4*pi, guess*5)], optimset('Display','off','TolX',1e-4));
                end
            catch
                best_alpha = 0;
            end
            [z_fin, data] = obj.calc_core(best_alpha, start_3d, target_3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir);
            res.data = data; res.l_dir = l_dir; res.d_dir = d_dir;
            final_err = z_fin - target_z;
            if isempty(data.dubins) || isempty(data.dubins.x)
                dubins_len = inf; 
            else
                dubins_len = sum(sqrt(diff(data.dubins.x).^2 + diff(data.dubins.y).^2));
            end
        end
        
        % =========================================================
        % ★★★ 修正版 calc_core ★★★
        % =========================================================
        function [final_z, data] = calc_core(obj, alpha, s3d, t3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir)
            
            rho0 = obj.AtmoModel.get_density(h0/1000);
            g_ref = obj.AtmoModel.get_gravity(h0);
            
            % 初期化
            curr = struct('x',s3d(1), 'y',s3d(2), 'z',s3d(3), 'yaw',s3d(4));
            current_time = 0;
            data.V_start = V0;
            
            % =====================================================
            % 0. Run-up (初期助走区間)
            % =====================================================
            rx=[]; ry=[]; rz=[]; rt=[]; rpsi=[]; rk=[];
            dist_run = 0;
            while dist_run < obj.RunUpDistance
                h_now = max(0, curr.z);
                rho = obj.AtmoModel.get_density(h_now/1000);
                V_TAS = V0 * sqrt(rho0 / rho);
                
                rx(end+1)=curr.x; ry(end+1)=curr.y; rz(end+1)=curr.z;
                rt(end+1)=current_time; rpsi(end+1)=curr.yaw; 
                rk(end+1)=0; % Kappa=0 (直線)
                
                curr.x = curr.x + step * cos(curr.yaw);
                curr.y = curr.y + step * sin(curr.yaw);
                curr.z = curr.z + step * tan_g_str;
                current_time = current_time + step / V_TAS;
                dist_run = dist_run + step;
            end
            data.runup = struct('x',rx, 'y',ry, 'z',rz, 't',rt, 'psi',rpsi, 'kappa',rk);
            
            % =====================================================
            % 1. Loiter (定常旋回)
            % =====================================================
            lx=[]; ly=[]; lz=[]; lt=[]; lpsi=[]; lk=[];
            
            if strcmp(l_dir, 'L')
                l_sign = -1;       
                k_loiter = -1.0/R; 
            else
                l_sign = 1; 
                k_loiter = 1.0/R;
            end
            
            tot = 0;
            while tot < alpha
                if curr.z <= 0, final_z = -9999; data.loiter=[]; data.runup_mid=[]; data.dubins=[]; data.final=[]; return; end
                
                h_km = max(0, curr.z/1000);
                rho = obj.AtmoModel.get_density(h_km);
                V = V0 * sqrt(rho0 / rho);
                
                tan_phi = abs(k_loiter) * V^2 / g_ref;
                cos_phi = 1.0 / sqrt(1 + tan_phi^2);
                tan_g = tan_g_str / cos_phi;
                
                d_ang = step / R;
                if tot + d_ang > alpha, d_ang = alpha - tot; move = d_ang * R; else, move = step; end
                
                V_horiz = V * cos_phi;
                dt = move / V_horiz;
                
                lx(end+1)=curr.x; ly(end+1)=curr.y; lz(end+1)=curr.z; 
                lt(end+1)=current_time; lpsi(end+1)=curr.yaw; lk(end+1)=k_loiter;
                
                if l_sign == -1, curr.yaw = curr.yaw - d_ang; else, curr.yaw = curr.yaw + d_ang; end
                
                curr.x = curr.x + move * cos(curr.yaw);
                curr.y = curr.y + move * sin(curr.yaw);
                curr.z = curr.z + move * tan_g;
                current_time = current_time + dt;
                
                tot = tot + d_ang;
            end
            data.loiter = struct('x',lx, 'y',ly, 'z',lz, 't',lt, 'psi',lpsi, 'kappa',lk);
            
            % =====================================================
            % 1.5. Mid Run-up (中間助走: 操舵反転バッファ) ★新規追加★
            % =====================================================
            mx=[]; my=[]; mz=[]; mt=[]; mpsi=[]; mk=[];
            dist_mid = 0;
            while dist_mid < obj.MidRunUpDistance
                h_now = max(0, curr.z);
                rho = obj.AtmoModel.get_density(h_now/1000);
                V_TAS = V0 * sqrt(rho0 / rho);
                
                mx(end+1)=curr.x; my(end+1)=curr.y; mz(end+1)=curr.z;
                mt(end+1)=current_time; mpsi(end+1)=curr.yaw; 
                mk(end+1)=0; % Kappa=0 (直線)
                
                curr.x = curr.x + step * cos(curr.yaw);
                curr.y = curr.y + step * sin(curr.yaw);
                curr.z = curr.z + step * tan_g_str;
                current_time = current_time + step / V_TAS;
                dist_mid = dist_mid + step;
            end
            data.runup_mid = struct('x',mx, 'y',my, 'z',mz, 't',mt, 'psi',mpsi, 'kappa',mk);

            % =====================================================
            % 2. Dubins
            % =====================================================
            psi_land = deg2rad(land_deg);
            fix_x = t3d(1) - Lf*cos(psi_land); 
            fix_y = t3d(2) - Lf*sin(psi_land);
            
            % Solver呼び出し (中間助走を終えたcurrの状態が自動的に始点になる)
            [px, py, pyaw, mode, lengths] = obj.dubins_solve(...
                [curr.x, curr.y, mod(curr.yaw, 2*pi)], ...
                [fix_x, fix_y, psi_land], ...
                1/R, step, d_dir);
            
            if isempty(px), final_z = -99999; data.dubins=[]; data.final=[]; return; end
            
            px = px(:); py = py(:); pyaw = pyaw(:);
            
            dz = zeros(size(px)); dz(1) = curr.z; 
            dt_dub = zeros(size(px)); dt_dub(1) = current_time;
            dk = zeros(size(px));
            
            d_s = sqrt([0; diff(px).^2 + diff(py).^2]);
            cum_dist = cumsum(d_s);
            L1 = lengths(1); L2 = lengths(2);
            
            for i=2:length(px)
                dist_now = cum_dist(i);
                if dist_now <= L1 + 1e-3
                    m = mode{1};
                elseif dist_now <= L1 + L2 + 1e-3
                    m = mode{2};
                else
                    m = mode{3};
                end
                
                if strcmp(m, 'L'), k_val = -1.0/R;
                elseif strcmp(m, 'R'), k_val = 1.0/R;
                else, k_val = 0; end
                
                dk(i) = k_val;
                is_turn = (k_val ~= 0);
                
                h_prev = max(0, dz(i-1));
                rho = obj.AtmoModel.get_density(h_prev/1000);
                V = V0 * sqrt(rho0 / rho);
                
                if ~is_turn
                    bg = tan_g_str;
                    cos_phi = 1.0;
                else
                    tan_phi = abs(k_val) * V^2 / g_ref;
                    cos_phi = 1.0 / sqrt(1 + tan_phi^2);
                    bg = tan_g_str / cos_phi;
                end
                
                dz(i) = dz(i-1) + d_s(i) * bg;
                dt_dub(i) = dt_dub(i-1) + d_s(i) / (V * cos_phi);
            end
            
            data.dubins = struct('x',px, 'y',py, 'z',dz, 't',dt_dub, 'psi',pyaw, 'kappa',dk);
            curr.z = dz(end);
            current_time = dt_dub(end);
            
            % =====================================================
            % 3. Final (直線区間)
            % =====================================================
            df = 0:step:Lf;
            fx = px(end) + df(:) * cos(psi_land); 
            fy = py(end) + df(:) * sin(psi_land);
            fpsi = repmat(psi_land, size(fx));
            fkappa = zeros(size(fx));
            
            tf = zeros(size(fx)); tf(1) = current_time; 
            zf = zeros(size(fx)); zf(1) = curr.z;
            
            for k=2:length(fx)
                dist = sqrt((fx(k)-fx(k-1))^2 + (fy(k)-fy(k-1))^2);
                h_prev = max(0, zf(k-1)); 
                rho = obj.AtmoModel.get_density(h_prev/1000);
                V = V0 * sqrt(rho0 / rho);
                
                tf(k) = tf(k-1) + dist/V; 
                zf(k) = zf(k-1) + dist*tan_g_str;
            end
            data.final = struct('x',fx, 'y',fy, 'z',zf, 't',tf, 'psi',fpsi, 'kappa',fkappa);
            final_z = zf(end);
        end
        
        function [px, py, pyaw, mode, lengths] = dubins_solve(~, s, e, c, step, start_cons)
             [px, py, pyaw, mode, lengths, ~] = dubins_path_planning(s, e, c, 1.0, step, start_cons);
        end
    end
end