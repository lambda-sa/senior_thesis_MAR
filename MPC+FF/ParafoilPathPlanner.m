classdef ParafoilPathPlanner < handle
    % PARAFOILPATHPLANNER (決定版: 幾何形状計算)
    
    properties
        AtmoModel
        ResultData      % .x, .y, .z, .t
        FinalError
        StartV
        
        % 経路計画パラメータ (初期値は仮置き。計算時に上書きされる)
        V0 = 15.0; 
        R_fixed = 50.0;
        MinLoiterAngle = 0;
        TurnInputRatio = 0.8;
        
        % 風速ベクトル (子クラスと共有)
        WindVector = [0; 0; 0]; 
    end
    
    methods
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
            
            % パターン探索 (L-L, R-R等)
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
            fprintf('  >>> 選択: %s-Loiter (DubinsLen: %.1f, Err: %.1f)\n', chosen_dir, best.len, best.err);
        end
        
        function [res, final_err, dubins_len] = try_combination(obj, l_dir, d_dir, start_3d, target_3d, land_deg, Lf, V0, gr, R)
            h0 = start_3d(3); step = 0.5; target_z = target_3d(3); tan_g_str = -1.0 / gr;
            cost_func = @(alpha) obj.calc_core(alpha, start_3d, target_3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir) - target_z;
            
            best_alpha = 0;
            try
                min_a = obj.MinLoiterAngle;
                z_at_min = obj.calc_core(min_a, start_3d, target_3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir);
                if z_at_min >= target_z
                    h_avg = (h0 + target_z)/2;
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
            if isempty(data.dubins.x), dubins_len = inf; else, dubins_len = sum(sqrt(sum(diff([data.dubins.x(:), data.dubins.y(:)]).^2, 2))); end
        end
        
        function [final_z, data] = calc_core(obj, alpha, s3d, t3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir)
            rho0 = obj.AtmoModel.get_density(h0/1000);
            g_ref = obj.AtmoModel.get_gravity(h0);
            sx = s3d(1); sy = s3d(2); sz = s3d(3); syaw = s3d(4);
            curr = struct('x',sx, 'y',sy, 'z',sz, 'yaw',syaw);
            lx=[sx]; ly=[sy]; lz=[sz]; lt=[0]; lpsi=[syaw];
            rho_st = obj.AtmoModel.get_density(max(0, sz/1000));
            data.V_start = V0 * sqrt(rho0 / rho_st);
            
            % 1. Loiter
            tot = 0; current_time = 0;
            if strcmp(l_dir, 'L'), l_sign = -1; else, l_sign = 1; end
            while tot < alpha
                if curr.z <= 0, final_z = -9999; data.loiter.x=[]; data.dubins.x=[]; data.final.x=[]; return; end
                h_km = curr.z/1000;
                rho = obj.AtmoModel.get_density(max(0, h_km));
                V = V0 * sqrt(rho0 / rho); % EAS -> TAS conversion for path integration
                sigma = atan(V^2 / (g_ref * R));
                tan_g = tan_g_str / cos(sigma);
                d_ang = step / R;
                if tot + d_ang > alpha, d_ang = alpha - tot; move = d_ang * R; else, move = step; end
                cos_gamma = 1 / sqrt(1 + tan_g^2); V_horiz = V * cos_gamma;
                dt = move / V_horiz; current_time = current_time + dt;
                if l_sign == -1, curr.yaw = curr.yaw + d_ang; else, curr.yaw = curr.yaw - d_ang; end
                curr.x = curr.x + move * cos(curr.yaw); curr.y = curr.y + move * sin(curr.yaw); curr.z = curr.z + move * tan_g;
                lx(end+1)=curr.x; ly(end+1)=curr.y; lz(end+1)=curr.z; lt(end+1)=current_time; lpsi(end+1)=curr.yaw;
                tot = tot + d_ang;
            end
            data.loiter.x=lx; data.loiter.y=ly; data.loiter.z=lz; data.loiter.t=lt; data.loiter.psi=lpsi;
            
            % 2. Dubins
            psi_land = deg2rad(land_deg);
            fix_x = t3d(1) - Lf*cos(psi_land); fix_y = t3d(2) - Lf*sin(psi_land);
            [px, py, pyaw, mode, lengths] = obj.dubins_solve([curr.x, curr.y, mod(curr.yaw, 2*pi)], [fix_x, fix_y, psi_land], 1/R, step, d_dir);
            if isempty(px), final_z = -99999; data.dubins.x=[]; data.final.x=[]; return; end
            dz = zeros(size(px)); dz(1) = curr.z; dt_dub = zeros(size(px)); dt_dub(1) = current_time;
            d_s = sqrt([0, diff(px).^2 + diff(py).^2]);
            cum = cumsum(d_s); L1 = lengths(1); L2 = lengths(2);
            for i=2:length(px)
                dist_step = d_s(i);
                if cum(i)<=L1, m=mode{1}; elseif cum(i)<=L1+L2, m=mode{2}; else, m=mode{3}; end
                is_turn = ~strcmp(m, 'S');
                h_prev = max(0, dz(i-1));
                rho = obj.AtmoModel.get_density(h_prev/1000);
                V = V0 * sqrt(rho0 / rho);
                if ~is_turn, bg = tan_g_str; else, sig = atan(V^2 / (g_ref * R)); bg = tan_g_str / cos(sig); end
                dz(i) = dz(i-1) + dist_step*bg; dt_dub(i) = current_time + dist_step/V; current_time = dt_dub(i);
            end
            data.dubins.x=px; data.dubins.y=py; data.dubins.z=dz; data.dubins.t=dt_dub; data.dubins.psi=pyaw;
            curr.z = dz(end);
            
            % 3. Final
            df = 0:step:Lf;
            data.final.x = px(end) + df*cos(psi_land); data.final.y = py(end) + df*sin(psi_land); data.final.psi = repmat(psi_land, size(df));
            tf = zeros(size(df)); tf(1) = current_time; zf = zeros(size(df)); zf(1) = curr.z;
            for k=2:length(df)
                dist = df(k)-df(k-1); h_prev = max(0, zf(k-1)); rho = obj.AtmoModel.get_density(h_prev/1000);
                V = V0 * sqrt(rho0 / rho);
                tf(k) = current_time + dist/V; current_time = tf(k);
                zf(k) = zf(k-1) + dist*tan_g_str;
            end
            data.final.z = zf; data.final.t = tf; final_z = zf(end);
        end
        
        function [px, py, pyaw, mode, lengths] = dubins_solve(~, s, e, c, step, start_cons)
             [px, py, pyaw, mode, lengths, ~] = dubins_path_planning(s, e, c, 1.0, step, start_cons);
        end
        
        function plot_result(obj)
            if isempty(obj.ResultData), return; end
            figure('Color','w'); hold on; grid on; axis equal; view(3);
            d = obj.ResultData;
            if ~isempty(d.loiter.x), plot3(d.loiter.x, d.loiter.y, d.loiter.z, 'c-', 'LineWidth', 2); end
            plot3(d.dubins.x, d.dubins.y, d.dubins.z, 'b-', 'LineWidth', 2);
            plot3(d.final.x, d.final.y, d.final.z, 'g-', 'LineWidth', 3);
            plot3(d.final.x(end), d.final.y(end), d.final.z(end), 'ro', 'MarkerFaceColor','r');
            title(sprintf('Optimal Path (Err: %.1fm)', obj.FinalError)); xlabel('X'); ylabel('Y'); zlabel('Z');
        end

        % --- ParafoilPathPlanner.m の methods ブロック内に追加 ---
        
        function plan_simple_physics_based(obj, type, duration, start_pos, V_EAS_trim, glide_ratio, bank_angle_deg)
            % PLAN_SIMPLE_PHYSICS_BASED
            % 入力バンク角から旋回半径を決定し、その半径を維持する軌道を生成する
            
            fprintf('--- Generating Constant Radius Trajectory (%s) ---\n', type);
            
            % 1. 初期化
            dt = 0.1;
            N = ceil(duration / dt) + 1;
            
            t = zeros(N, 1);
            x = zeros(N, 1); y = zeros(N, 1); z = zeros(N, 1);
            psi = zeros(N, 1);
            
            % 保存用: 実際に適用されたバンク角履歴
            phi_log = zeros(N, 1);
            
            % 初期値セット
            x(1) = start_pos(1);
            y(1) = start_pos(2);
            z(1) = start_pos(3);
            psi(1) = start_pos(4);
            
            % 物理定数
            rho_0 = 1.225;
            g_0 = 9.80665; % 標準重力（またはAtmoModelから取得しても良い）
            
            % --- 2. 旋回半径 R の固定 ---
            % 開始高度での状態を取得
            rho_start = obj.AtmoModel.get_density(max(0, z(1))/1000);
            V_TAS_start = V_EAS_trim * sqrt(rho_0 / rho_start);
            g_start = obj.AtmoModel.get_gravity(z(1));
            
            if strcmpi(type, 'Straight') || abs(bank_angle_deg) < 0.1
                R_origin = inf; % 直線
                target_phi_rad = 0;
                fprintf('  -> Mode: Straight (R = inf)\n');
            else
                % 入力バンク角から半径を決定
                phi_input = deg2rad(bank_angle_deg);
                % R = V^2 / (g * tan(phi))
                R_origin = V_TAS_start^2 / (g_start * tan(phi_input));
                
                % 符号（旋回方向）の保持
                turn_sign = sign(bank_angle_deg);
                R_origin = abs(R_origin) * turn_sign; % 左旋回ならマイナスの半径として扱うか、計算時に考慮
                
                fprintf('  -> Mode: Constant Radius Turn\n');
                fprintf('     Input Bank: %.1f deg, Start TAS: %.1f m/s\n', bank_angle_deg, V_TAS_start);
                fprintf('     => Fixed Radius: %.1f m\n', abs(R_origin));
            end
    
            phi_log(1) = deg2rad(bank_angle_deg);
            
            % 基本滑空角 (バンク0のとき)
            tan_gamma_base = -1.0 / glide_ratio;
            
            % --- 3. 積分ループ ---
            for i = 1:N-1
                t(i+1) = t(i) + dt;
                
                % 現在の環境
                h_curr = max(0, z(i));
                rho_curr = obj.AtmoModel.get_density(h_curr / 1000);
                g = obj.AtmoModel.get_gravity(h_curr);
                
                % 現在のTAS
                V_TAS = V_EAS_trim * sqrt(rho_0 / rho_curr);
                
                % --- 制御量の決定 ---
                if isinf(R_origin)
                    % 直線
                    dot_psi = 0;
                    phi_curr = 0;
                else
                    % 旋回: 固定半径 R を守るためのバンク角を逆算
                    % phi_curr = atan( V^2 / (g * R) )
                    % ※ R_fixedには符号が含まれていると想定
                    phi_curr = atan( V_TAS^2 / (g * R_origin) );
                    
                    % 旋回率
                    % dot_psi = V / R
                    dot_psi = V_TAS / R_origin;
                end
                
                phi_log(i+1) = phi_curr;
                
                % --- 物理挙動の更新 ---
                
                % バンク角による沈下率補正: tan(gamma) = tan(gamma_base) / cos(phi)
                % ※ phi_curr が大きくなると沈下も速くなる
                tan_gamma_curr = tan_gamma_base / cos(phi_curr);
                
                dz_dt = V_TAS * tan_gamma_curr; % 降下速度 (負)
                V_horiz = sqrt(max(0, V_TAS^2 - dz_dt^2));
                
                % 状態積分
                psi(i+1) = psi(i) + dot_psi * dt;
                x(i+1)   = x(i) + V_horiz * cos(psi(i)) * dt;
                y(i+1)   = y(i) + V_horiz * sin(psi(i)) * dt;
                z(i+1)   = z(i) + dz_dt * dt;
                
                if z(i+1) <= 0
                    z(i+1) = 0;
                    t = t(1:i+1); x = x(1:i+1); y = y(1:i+1); z = z(1:i+1); psi = psi(1:i+1); phi_log = phi_log(1:i+1);
                    break;
                end
            end
            
            % 結果格納
            d.loiter.x=[]; d.loiter.y=[]; d.loiter.z=[]; d.loiter.t=[];
            d.dubins.x=[]; d.dubins.y=[]; d.dubins.z=[]; d.dubins.t=[];
            d.final.x = x'; d.final.y = y'; d.final.z = z'; d.final.t = t';
            d.final.psi = psi';
            
            % ★ここ重要: 計算されたバンク角履歴を保存しておく場所がないため、
            % 一旦 ResultDataには入れず、軌道形状から compute_dual_trajectories で
            % 正しく再計算されることに期待します。
            % (直前の回答の「dpsi_airから逆算する修正」が入っていれば、これで完璧に整合します)
            
            obj.ResultData = d;
            obj.FinalError = 0;
            obj.StartV = V_EAS_trim;
            obj.V0 = V_EAS_trim; 
            
            % 旋回半径を保存（デバッグ用）
            if ~isinf(R_origin)
                obj.R_fixed = abs(R_origin);
            end
        end
    end
end