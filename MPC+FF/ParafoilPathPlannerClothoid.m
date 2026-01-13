classdef ParafoilPathPlannerClothoid < ParafoilPathPlanner
    % PARAFOILPATHPLANNERCLOTHOID (Bug Fixed Version)
    %
    % 修正点:
    %   calc_core内の異常終了時に data.dubins = [] としていた箇所を
    %   data.dubins.x = [] に変更。これにより親クラスでのドットアクセスエラーを回避。
    
    properties
        MaxRollRate = 0.5;     % [deg/s]
        LookAheadFactor = 0.5;  % 先行係数
        RunUpDistance = 0.0;   % 助走距離 [m]
        FineStepMultiplier = 5; % 積分精度
        
        ResultDataWind   % 対地データ
        ResultDataNaive  % 対気データ
    end
    
    methods
        function obj = ParafoilPathPlannerClothoid(atmo_model)
            obj@ParafoilPathPlanner(atmo_model);
            obj.WindVector = [0; 0; 0];
        end
        
        % =========================================================
        % メイン計算ロジック (Calc Core)
        % =========================================================
        function [final_z, data] = calc_core(obj, alpha, s3d, t3d, land_deg, Lf, step, V0, h0, R, tan_g_str, l_dir, d_dir)
            
            rho0 = obj.AtmoModel.get_density(h0/1000);
            g_ref = obj.AtmoModel.get_gravity(h0);
            
            curr = struct('x',s3d(1), 'y',s3d(2), 'z',s3d(3), 'yaw',s3d(4));
            current_time = 0;
            
            rho_st = obj.AtmoModel.get_density(max(0, curr.z/1000));
            data.V_start = V0 * sqrt(rho0 / rho_st);
            
            max_dkappa_ds = (9.80665 * deg2rad(obj.MaxRollRate)) / V0^2;%dκ/dsの最大値　1m進むごとに，どれくらい曲率を変えていいか
            k_loiter = (strcmp(l_dir, 'L') * 1 + strcmp(l_dir, 'R') * -1) / R;
            current_kappa = 0;
            
            % --- Phase -1: Run-up ---
            %空配列を回避したい→最初から正しい構造体フォーマット
            rx=[curr.x]; ry=[curr.y]; rz=[curr.z]; rt=[current_time]; rpsi=[curr.yaw];
            dist_run = 0;
            while dist_run < obj.RunUpDistance
                [curr, current_time] = obj.integrate_step(curr, 0, step, V0, rho0, g_ref, tan_g_str, current_time);
                rx(end+1)=curr.x; ry(end+1)=curr.y; rz(end+1)=curr.z; rt(end+1)=current_time; rpsi(end+1)=curr.yaw;
                dist_run = dist_run + step;
            end
            data.runup.x=rx; data.runup.y=ry; data.runup.z=rz; data.runup.t=rt; data.runup.psi=rpsi;
            
            % --- Phase 0: Entry ---　滑空→Loiterの遷移区間
            ex=[]; ey=[]; ez=[]; et=[]; epsi=[];
            while abs(current_kappa - k_loiter) > 1e-6
                % ★修正1: 異常終了時のデータ構造を維持する
                if curr.z <= 0
                    final_z = -9999; 
                    data.loiter.x=[]; data.dubins.x=[]; data.final.x=[]; % ここを修正
                    return; 
                end
                
                fine_step = step / obj.FineStepMultiplier;
                for k=1:obj.FineStepMultiplier
                    d_k = k_loiter - current_kappa;
                    current_kappa = current_kappa + sign(d_k) * min(abs(d_k), max_dkappa_ds * fine_step);
                    [curr, current_time] = obj.integrate_step(curr, current_kappa, fine_step, V0, rho0, g_ref, tan_g_str, current_time);
                end
                ex(end+1)=curr.x; ey(end+1)=curr.y; ez(end+1)=curr.z; et(end+1)=current_time; epsi(end+1)=curr.yaw;
            end
            data.entry.x=ex; data.entry.y=ey; data.entry.z=ez; data.entry.t=et; data.entry.psi=epsi;
            
            % --- Phase 2 (Dry Run) --- 一回リハーサルで走らせてみる
            psi_land = deg2rad(land_deg);
            fix_x = t3d(1) - Lf*cos(psi_land); fix_y = t3d(2) - Lf*sin(psi_land);%風のオフセットあり(一応)での目標点
            
            [~, ~, ~, mode, lengths] = obj.dubins_solve(...
                [curr.x, curr.y, mod(curr.yaw, 2*pi)], [fix_x, fix_y, psi_land], 1/R, step, d_dir);
            
            % ★修正2: Dubins解なしの場合
            if isempty(mode)
                final_z = -99999; 
                data.dubins.x=[]; data.final.x=[]; data.loiter.x=[]; % ここも念のため統一
                return; 
            end
            
            [~, ~, z_dry_end, ~] = obj.simulate_dubins_sequence(...
                curr, current_time, current_kappa, mode, lengths, step, V0, rho0, g_ref, tan_g_str, R, max_dkappa_ds);
            
            delta_z_dubins = curr.z - z_dry_end;
            delta_z_final = abs(Lf * tan_g_str);
            
            z_to_burn = (curr.z - t3d(3)) - (delta_z_dubins + delta_z_final);
            
            h_now = max(0, curr.z); V_TAS = V0 * sqrt(rho0 / obj.AtmoModel.get_density(h_now/1000));
            tan_phi = abs(current_kappa) * V_TAS^2 / g_ref;
            dz_ds_loiter = abs(tan_g_str / (1 / sqrt(1 + tan_phi^2)));
            
            if z_to_burn < -10, req_loiter_dist = 0; else, req_loiter_dist = z_to_burn / dz_ds_loiter; end %余剰高度の算出
            
            % --- Phase 1: Loiter ---
            lx=[]; ly=[]; lz=[]; lt=[]; lpsi=[];
            dist_loitered = 0;
            while dist_loitered < req_loiter_dist
                 % ★修正3: 異常終了時
                 if curr.z <= 0
                     final_z = -9999; 
                     data.loiter.x=[]; data.dubins.x=[]; data.final.x=[]; 
                     return; 
                 end
                 d_dist = step;
                 if dist_loitered + d_dist > req_loiter_dist, d_dist = req_loiter_dist - dist_loitered; end
                 [curr, current_time] = obj.integrate_step(curr, current_kappa, d_dist, V0, rho0, g_ref, tan_g_str, current_time);
                 lx(end+1)=curr.x; ly(end+1)=curr.y; lz(end+1)=curr.z; lt(end+1)=current_time; lpsi(end+1)=curr.yaw;
                 dist_loitered = dist_loitered + d_dist;
            end
            data.loiter.x=lx; data.loiter.y=ly; data.loiter.z=lz; data.loiter.t=lt; data.loiter.psi=lpsi;

            % --- Phase 2: Dubins Path (Replan) ---　高度消費が終わったタイミングで，本番の軌道を求める
            [~, ~, ~, mode, lengths] = obj.dubins_solve(...
                [curr.x, curr.y, mod(curr.yaw, 2*pi)], [fix_x, fix_y, psi_land], 1/R, step, d_dir);
             
             % ★修正4: Dubins解なしの場合
             if isempty(mode)
                 final_z = -99999; 
                 data.dubins.x=[]; data.final.x=[]; 
                 return; 
             end
            %Dubinsにクロゾイドを追加する
            [D_res, final_state] = obj.run_dubins_integration(...
                 curr, current_time, current_kappa, mode, lengths, step, V0, rho0, g_ref, tan_g_str, R, max_dkappa_ds);
             
            data.dubins = D_res;
            curr = final_state;
            current_time = D_res.t(end);

            % --- Phase 3: Final ---
            df = 0:step:Lf; final_yaw = curr.yaw; 
            data.final.x = curr.x + df*cos(final_yaw); 
            data.final.y = curr.y + df*sin(final_yaw); 
            data.final.psi = repmat(final_yaw, size(df));
            % ★★★ この行を追加してください！ ★★★
            data.final.kappa = zeros(size(df)); 
            % （解説: Final区間は直線なので、曲率はすべてゼロです）
            zf = zeros(size(df)); zf(1) = curr.z; tf = zeros(size(df)); tf(1) = current_time;
            for k=2:length(df)
                h_now = max(0, zf(k-1)); V_TAS = V0 * sqrt(rho0 / obj.AtmoModel.get_density(h_now/1000));
                tf(k) = tf(k-1) + step/V_TAS; 
                zf(k) = zf(k-1) + step * tan_g_str; 
            end
            data.final.z = zf; data.final.t = tf; final_z = zf(end);
        end
        
        % =========================================================
        % Dubins積分制御 (Safe Look-ahead Logic)　
        % =========================================================
        function [Data, EndState] = run_dubins_integration(obj, start_curr, start_time, start_kappa, mode, lengths, step, V0, rho0, g_ref, tan_g_str, R, max_dkappa_ds)
            
            curr = start_curr; current_time = start_time; current_kappa = start_kappa;
            
            % --- ★修正1: 曲率(kappa)も保存する配列を作る ---
            dx=[]; dy=[]; dz=[]; dt=[]; dpsi=[];
            dx(1)=curr.x; dy(1)=curr.y; dz(1)=curr.z; dt(1)=current_time; dpsi(1)=curr.yaw;
            dkappa = [current_kappa]; % 曲率履歴
            
            % --- ★修正2: 切り替え点を記録する変数を追加 ---
            seg_end_indices = zeros(1, 3);
            idx_clothoid_start = 0; % 第1旋回(Loiter)終了 -> クロソイド開始点
            
            est_bank = atan(V0^2 / (g_ref * R));
            theoretical_lead = (V0 * est_bank) / deg2rad(obj.MaxRollRate) * obj.LookAheadFactor;
            
            total_len = sum(lengths); dist_accum = 0; seg_idx = 1; dist_in_seg = 0;
            
            while dist_accum < total_len
                cur_len = lengths(seg_idx); cur_m = mode{seg_idx};
                if seg_idx < 3, nxt_m = mode{seg_idx+1}; else, nxt_m = 'F'; end
                base_k = obj.char2kappa(cur_m, R); next_k = obj.char2kappa(nxt_m, R);
                
                L_lead = min(theoretical_lead, cur_len * 0.45);
                
                % --- クロソイド判定 ---
                if (cur_len - dist_in_seg) < L_lead
                    target_k = next_k;
                    
                    % ★修正3: 「第1セグメント」かつ「クロソイドに入った瞬間」を記録
                    if seg_idx == 1 && idx_clothoid_start == 0
                        idx_clothoid_start = length(dx);
                    end
                else
                    target_k = base_k; 
                end
                
                % --- 積分 ---
                d_k = target_k - current_kappa;
                if abs(d_k) > 1e-6
                    fine_step = step / obj.FineStepMultiplier;
                    for k=1:obj.FineStepMultiplier
                        d_k_sub = target_k - current_kappa;
                        current_kappa = current_kappa + sign(d_k_sub) * min(abs(d_k_sub), max_dkappa_ds * fine_step);
                         [curr, current_time] = obj.integrate_step(curr, current_kappa, fine_step, V0, rho0, g_ref, tan_g_str, current_time);
                    end
                else
                     [curr, current_time] = obj.integrate_step(curr, current_kappa, step, V0, rho0, g_ref, tan_g_str, current_time);
                end
                
                % --- データ保存 ---
                dx(end+1)=curr.x; dy(end+1)=curr.y; dz(end+1)=curr.z; dt(end+1)=current_time; dpsi(end+1)=curr.yaw;
                
                % ★修正4: 曲率も保存
                dkappa(end+1) = current_kappa;
                
                dist_accum = dist_accum + step; dist_in_seg = dist_in_seg + step;
                
                % --- セグメント切り替え判定 ---
                if dist_in_seg >= cur_len
                    seg_end_indices(seg_idx) = length(dx);
                    
                    if seg_idx < 3
                        seg_idx = seg_idx+1; dist_in_seg=0; 
                    else
                        break; 
                    end
                end
            end
            
            % --- 結果出力 ---
            Data.x = dx; Data.y = dy; Data.z = dz; Data.t = dt; Data.psi = dpsi; 
            
            % ★修正5: 必要な情報を構造体に詰める
            Data.kappa = dkappa; % Autopilot用
            
            % Scheduler用: もし見つからなかったら(短すぎる場合など)セグメント終了点を使う
            if idx_clothoid_start == 0
                idx_clothoid_start = seg_end_indices(1);
            end
            Data.idx_clothoid_start = idx_clothoid_start; 
            Data.segment_end_idx = seg_end_indices;
            
            EndState = curr;
        end
        %simulate_dubins_sequence...結果を保存せずに、最終形態だけを知りたいときに用いる
        %最終到達点をget
        function [xf, yf, zf, tf] = simulate_dubins_sequence(obj, start_curr, start_time, start_kappa, mode, lengths, step, V0, rho0, g_ref, tan_g_str, R, max_dkappa_ds)
            [~, EndState] = obj.run_dubins_integration(start_curr, start_time, start_kappa, mode, lengths, step, V0, rho0, g_ref, tan_g_str, R, max_dkappa_ds);
            xf = EndState.x; yf = EndState.y; zf = EndState.z; tf = 0; 
        end

        function [next_curr, next_t] = integrate_step(obj, curr, kappa, step, V0, rho0, g_ref, tan_g_str, t_now)
            % 1. 現在高度での対気速度 (TAS) の計算
            h_now = max(0, curr.z); V_TAS = V0 * sqrt(rho0 / obj.AtmoModel.get_density(h_now/1000));

            % 2. 曲率(kappa)からバンク角(phi)を逆算
            % 物理式: tan(phi) = (V^2 * kappa) / g
            tan_phi = kappa * V_TAS^2 / g_ref; cos_phi = 1 / sqrt(1 + tan_phi^2);

            % 3. バンク角による沈下率の補正 ★重要
            % 旋回中は揚力が傾くため、高度維持に必要な揚力が増え、抗力も増える＝沈下が速くなる
            % 簡易モデルとして dz/ds = 基本沈下率 / cos(phi) で計算
            dz_ds = tan_g_str / cos_phi;
            next_curr.yaw = curr.yaw + kappa * step;
            next_curr.x = curr.x + step * cos(curr.yaw);
            next_curr.y = curr.y + step * sin(curr.yaw);
            next_curr.z = curr.z + dz_ds * step;
            next_t = t_now + step / (V_TAS * cos_phi);
        end
        
        function save_naive_path(obj), if ~isempty(obj.ResultData), obj.ResultDataNaive = obj.ResultData; end, end
        function apply_wind_effect(obj, wx, wy)
            obj.WindVector = [wx; wy; 0];
            if isempty(obj.ResultData), return; end
            d = obj.ResultData; w_data = d; 
            fields = {'runup', 'entry', 'loiter', 'dubins', 'final'};
            for i = 1:length(fields)
                f = fields{i}; if isfield(d, f) && ~isempty(d.(f).t)
                    w_data.(f).x = d.(f).x + wx * d.(f).t; w_data.(f).y = d.(f).y + wy * d.(f).t;
                end
            end
            obj.ResultDataWind = w_data;
        end
        function plot_wind_comparison(obj)
            if isempty(obj.ResultDataWind), warning('No Data'); return; end
            figure('Color', 'w', 'Name', 'Wind Effect'); hold on; grid on; axis equal; view(3);
            w = obj.ResultDataWind;
            if isfield(w, 'runup') && ~isempty(w.runup.x), plot3(w.runup.y, w.runup.x, w.runup.z, 'k-', 'LineWidth', 2, 'DisplayName', 'Run-up'); end
            if isfield(w, 'entry') && ~isempty(w.entry.x), plot3(w.entry.y, w.entry.x, w.entry.z, 'g-', 'LineWidth', 2, 'DisplayName', 'Entry'); end
            if ~isempty(w.loiter.x), plot3(w.loiter.y, w.loiter.x, w.loiter.z, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Loiter'); end
            if ~isempty(w.dubins.x), plot3(w.dubins.y, w.dubins.x, w.dubins.z, 'b-', 'LineWidth', 2, 'DisplayName', 'Dubins'); end
            if ~isempty(w.final.x), plot3(w.final.y, w.final.x, w.final.z, 'r-', 'LineWidth', 2, 'DisplayName', 'Final'); end
            legend;
        end
        function k = char2kappa(~, m, R), switch m, case 'L', k=1/R; case 'R', k=-1/R; otherwise, k=0; end; end
        function plan_simple_physics_based(obj, type, duration, start_pos, V_EAS_trim, glide_ratio, bank_angle_deg)
             % (RunUp=0対応版の単純シミュレーションコード。前回の回答と同じものを記述してください)
             fprintf('--- Generating Smooth Trajectory (Simple) ---\n');
            dt_base = 0.1; t = 0; x = start_pos(1); y = start_pos(2); z = start_pos(3); psi = start_pos(4);
            R_t=[t]; R_x=[x]; R_y=[y]; R_z=[z]; R_psi=[psi];
            rho_0 = 1.225; tan_gamma_base = -1.0 / glide_ratio;
            if strcmpi(type, 'Straight'), target_phi=0; else, target_phi=deg2rad(bank_angle_deg); end
            current_phi = 0; max_rr = deg2rad(obj.MaxRollRate);
            runup_time = obj.RunUpDistance / V_EAS_trim;
            while t < runup_time
                t = t + dt_base; h = max(0, z); V_TAS = V_EAS_trim * sqrt(rho_0 / obj.AtmoModel.get_density(h/1000));
                x = x + V_TAS * cos(psi) * dt_base; y = y + V_TAS * sin(psi) * dt_base; z = z + V_TAS * tan_gamma_base * dt_base;
                R_t(end+1)=t; R_x(end+1)=x; R_y(end+1)=y; R_z(end+1)=z; R_psi(end+1)=psi;
            end
            while t < duration
                err_phi = target_phi - current_phi;
                if abs(err_phi) > 1e-4, dt = dt_base / obj.FineStepMultiplier; else, dt = dt_base; end
                t = t + dt; h = max(0, z); V_TAS = V_EAS_trim * sqrt(rho_0 / obj.AtmoModel.get_density(h/1000));
                d_phi = sign(err_phi) * min(abs(err_phi), max_rr * dt);
                current_phi = current_phi + d_phi;
                if abs(current_phi) < 1e-6, k=0; else, k = (9.80665*tan(current_phi))/V_TAS^2; end
                psi = psi + k * V_TAS * dt; dz = V_TAS * (tan_gamma_base / cos(current_phi)) * dt;
                x = x + V_TAS * cos(psi) * dt; y = y + V_TAS * sin(psi) * dt; z = z + dz;
                R_t(end+1)=t; R_x(end+1)=x; R_y(end+1)=y; R_z(end+1)=z; R_psi(end+1)=psi;
                if z<=0, break; end
            end
            d.final.x=R_x; d.final.y=R_y; d.final.z=R_z; d.final.t=R_t; d.final.psi=R_psi;
            d.loiter.x=[]; d.dubins.x=[]; d.entry.x=[]; d.runup.x=[];
            obj.ResultData = d; obj.FinalError = 0; obj.StartV = V_EAS_trim; obj.V0 = V_EAS_trim;
        end
    end
end