classdef ParafoilAutopilot < handle
    % PARAFOILAUTOPILOT (Thesis Final Version)
    %
    % 概要:
    %   VF法(空間微分)、厳密な風補正(Kinematics)、および
    %   逆ダイナミクス(Inverse Dynamics)を統合したパラフォイル誘導制御クラス。
    %   クロソイド曲線を含む高度な経路追従に対応。
    
    properties
        % --- 機体・環境 ---
        Params
        Yaw_Factor
        AtmoModel
        Linearizer
        
        % --- ゲイン ---
        Gains
        
        % --- 状態管理 ---
        CurrentMode = 'Loiter'; 
        LoiterParams
        ReferencePathMatrix = []; % [x, y, z, t, psi, kappa] (6列構成)
        
        LastIndex = 1;
        WindowSize = 100;
        NeedGlobalSearch = true;
        
        PsiErrIntegral = 0; 
        LastUpdateTime = -1;
        UpdateInterval = 0.05; 
        
        LastDeltaR = 0;     
        LastDeltaL = 0;     
        IsReady = false;    
    end
    
    methods
        % =========================================================
        % コンストラクタ
        % =========================================================
        function obj = ParafoilAutopilot(params, linearizer)
            obj.Params = params;
            if nargin > 1, obj.Linearizer = linearizer; end
            
            try
                obj.AtmoModel = AtmoTempPressRho();
            catch
                obj.AtmoModel = [];
            end
            
            % --- ゲイン設定 ---
            obj.Gains.k_vf_loiter  = 0.08;
            obj.Gains.k_vf_mission = 0.05;
            obj.Gains.chi_inf      = pi/2;
            
            obj.Gains.K_guidance   = 0.0;    % レートFB (基本0で運用可能)
            obj.Gains.K_linear     = 1.0;    % 線形化FF有効
            
            % バンク角制御 (Inner Loop)
            obj.Gains.kp_phi = 0.0; % 論文の構成に合わせて調整してください
            
            % 補助PID
            obj.Gains.kp_psi = 1.0; 
            obj.Gains.ki_psi = 0.0; 
            
            % --- Yaw_Factor 計算 ---
            if isfield(params, 'prop'), b=params.prop.b; else, b=params.b; end
            if isfield(params, 'params')
                d=params.params.d; Cnr=params.params.C_n_r; Cnda=params.params.C_n_delta_a;
            else
                d=params.d; Cnr=params.C_n_r; Cnda=params.C_n_delta_a;
            end
            
            base_factor = -d * b * Cnr / (2 * Cnda);
            obj.Yaw_Factor = base_factor * 1.0;
            
            obj.IsReady = true;
        end
        
        % =========================================================
        % データセット
        % =========================================================
        function set_target_path(obj, path_matrix)
            if size(path_matrix, 2) < 2, return; end
            obj.ReferencePathMatrix = path_matrix;
            obj.LastIndex = 1;
            obj.NeedGlobalSearch = true;
        end
        
        function import_mission_data(obj, planner)
            % Missionデータの取り込み (kappa列対応)
            if isprop(planner, 'ResultData') && isfield(planner.ResultData, 'dubins')
                d = planner.ResultData;
                if isfield(d, 'dubins') && ~isempty(d.dubins.x)
                    dx = [d.dubins.x(:); d.final.x(:)];
                    dy = [d.dubins.y(:); d.final.y(:)];
                    dz = [d.dubins.z(:); d.final.z(:)];
                    dt = [d.dubins.t(:); d.final.t(:)];
                    dpsi = [d.dubins.psi(:); d.final.psi(:)];
                    
                    % ★修正: 曲率 kappa の読み込み
                    if isfield(d.dubins, 'kappa')
                        dk = [d.dubins.kappa(:); d.final.kappa(:)];
                    else
                        dk = zeros(size(dx)); % データがない場合は0 (直線扱い)
                    end
                    
                    % 6列構成でセット
                    obj.set_target_path([dx, dy, dz, dt, dpsi, dk]);
                end
            end
        end
        
        function set_mode(obj, mode_str)
            if any(strcmpi(mode_str, {'Loiter', 'Mission'}))
                if ~strcmpi(obj.CurrentMode, mode_str)
                    obj.CurrentMode = mode_str;
                    obj.PsiErrIntegral = 0;
                    if strcmpi(mode_str, 'Mission')
                        obj.NeedGlobalSearch = true;
                    end
                    fprintf('Autopilot: Mode switched to [%s]\n', upper(mode_str));
                end
            end
        end
        
        % =========================================================
        % ★★★ メイン更新ループ (Update) ★★★
        % =========================================================
        function [delta_R, delta_L, log_struct] = update(obj, t, current_state, wind_vec, ~)
            
            if ~obj.IsReady
                delta_R=0; delta_L=0; log_struct=[]; return;
            end
            if (obj.LastUpdateTime >= 0) && (t - obj.LastUpdateTime < obj.UpdateInterval)
                delta_R = obj.LastDeltaR; delta_L = obj.LastDeltaL; log_struct=[]; return;
            end
            dt = obj.UpdateInterval;
            
            % --- 状態量 ---
            u=current_state(1); v=current_state(2); w=current_state(3);
            phi=current_state(7); theta=current_state(8); psi=current_state(9);
            pos_ned = current_state(10:12);
            
            V_tas = sqrt(u^2+v^2+w^2);
            V_horiz = max(0.1, V_tas * cos(theta)); % 対気水平速度
            
            % 現在の対地速度・コース角の計算
            R_b2n = obj.rotation_matrix(phi, theta, psi);
            Vel_ned = R_b2n * [u;v;w];
            Vg_n = Vel_ned(1);
            Vg_e = Vel_ned(2) ;
            Vg_horiz = sqrt(Vg_n^2 + Vg_e^2);
            chi_curr = atan2(Vg_e, Vg_n);
            
            r_req = 0; cmd_chi = 0; cmd_psi = 0; % 初期化
            
            % =====================================================
            % Step 1: Geometry & Kinematics (誘導則)
            % =====================================================
            if strcmpi(obj.CurrentMode, 'Loiter') && ~isempty(obj.LoiterParams)
                % === Loiterモード ===
                p = obj.LoiterParams;
                dx = pos_ned(1) - p.xc; dy = pos_ned(2) - p.yc;
                d = sqrt(dx^2 + dy^2);
                phi_pos = atan2(dy, dx);
                
                ang_diff = chi_curr - phi_pos;
                dot_phi_pos = Vg_horiz * sin(ang_diff) / d;
                
                k = obj.Gains.k_vf_loiter;
                err = d - p.R;
                dot_d = Vg_horiz * cos(ang_diff);
                dot_chi_app = (p.lambda * k * dot_d) / (1 + (k * err)^2);
                
                dot_chi_cmd = dot_phi_pos + dot_chi_app;
                
                % 厳密風補正メソッドへ
                % Loiter中はターゲット方位が刻々と変わるため、簡便に chi_curr をターゲットとみなすか
                % または phi_pos + 90deg を使う
                current_chi_target = chi_curr; 
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_vec, psi, current_chi_target, phi, theta, V_tas);
                
            elseif strcmpi(obj.CurrentMode, 'Mission') && ~isempty(obj.ReferencePathMatrix)
                % === Missionモード (厳密版) ===
                
                idx = obj.find_closest_index(pos_ned(1), pos_ned(2));
                P = obj.ReferencePathMatrix;
                path_x = P(idx, 1); path_y = P(idx, 2); path_psi = P(idx, 5);
                
                % ★修正: 曲率 kappa の取得
                if size(P, 2) >= 6
                    kappa_ref = P(idx, 6);
                else
                    kappa_ref = 0;
                end
                
                dx = pos_ned(1) - path_x; dy = pos_ned(2) - path_y;
                e_cross = -sin(path_psi)*dx + cos(path_psi)*dy;
                chi_rel = obj.normalize_angle(chi_curr - path_psi);
                
                % (A) 経路追従項 (Path Term) - 厳密幾何補正
                % Scale Factor = 1 / (1 - kappa * e)
                % 特異点保護: 分母がゼロ付近にならないようにする
                denom = 1.0 - kappa_ref * e_cross;
                if abs(denom) < 0.1, denom = sign(denom) * 0.1; end
                scale_factor = 1.0 / denom;
                
                % dot_chi_path = Vg * kappa * scale * cos(chi_rel)
                dot_chi_path = Vg_horiz * kappa_ref * scale_factor * cos(chi_rel);
                
                % (B) 誤差収束項 (VF Term)
                k_vf = obj.Gains.k_vf_mission;
                chi_inf = obj.Gains.chi_inf;
                dot_e = Vg_horiz * sin(chi_rel);
                
                numerator = -(chi_inf * 2 / pi) * (k_vf * dot_e);
                denominator = 1 + (k_vf * e_cross)^2;
                dot_chi_app = numerator / denominator;
                
                dot_chi_cmd = dot_chi_path + dot_chi_app;
                
                % (C) 目標コース角の算出 (ログ・補正用)
                chi_app_val = -(chi_inf * 2 / pi) * atan(k_vf * e_cross);
                current_chi_target = path_psi + chi_app_val;
                
                % 厳密風補正メソッドへ
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_vec, psi, current_chi_target, phi, theta, V_tas);
            end
            
            % =====================================================
            % Step 2: Dynamics & Control (制御則)
            % =====================================================
            
            % 目標バンク角 (対気速度 V_horiz を使用)
            g = 9.81;
            phi_ref = atan( (V_horiz * r_req) / g );
            phi_ref = max(-0.8, min(0.8, phi_ref)); % リミッタ
            
            % 線形化フィードフォワード
            ref_state.V = V_tas;
            ref_state.phi = phi_ref;
            ref_state.theta = theta;
            ref_state.r_req = r_req;
            if abs(u) > 0.1, ref_state.alpha = atan2(w, u); else, ref_state.alpha = 0; end
            
            [delta_delta_a, da_ref] = obj.calc_linear_correction(current_state, ref_state);
            dda_term = obj.Gains.K_linear * delta_delta_a;
            
            % 残留誤差補正 (PID)
            err_psi = obj.normalize_angle(cmd_psi - psi);
            obj.PsiErrIntegral = max(-0.3, min(0.3, obj.PsiErrIntegral + err_psi * dt));
            
            phi_pid = obj.Gains.kp_psi * err_psi + obj.Gains.ki_psi * obj.PsiErrIntegral;
            phi_cmd_final = phi_ref + phi_pid;
            
            % バンク角FB
            da_fb = obj.Gains.kp_phi * (phi_cmd_final - phi);
            
            % ミキシング
            da_total = da_ref + dda_term + da_fb;
            [delta_R, delta_L] = obj.apply_mixing(da_total, 0);
            
            % 更新
            obj.LastUpdateTime = t;
            obj.LastDeltaR = delta_R; obj.LastDeltaL = delta_L;
            
            log_struct.chi_cmd = cmd_chi;
            log_struct.psi_cmd = cmd_psi;
            log_struct.r_req   = r_req;
            log_struct.phi_ref = phi_ref;
            log_struct.da_total= da_total;
            log_struct.mode    = obj.CurrentMode;
        end
    end
    
    methods (Access = private)
        
        % =========================================================
        % ★重要修正: 厳密物理モデルによる風補正 & 座標変換
        % =========================================================
        function [r_req, chi_cmd_out, psi_cmd_out] = apply_wind_and_feedback(obj, ...
                dot_chi_cmd, Vg, wind_vec, psi, chi_target_now, current_phi, current_theta, V_tas)
            
            % 1. 偏流角 (Crab Angle) eta の計算
            Wx = wind_vec(1); Wy = wind_vec(2);
            W_mag = norm([Wx, Wy]);
            chi_wind = atan2(Wy, Wx);
            
            chi_diff = chi_target_now - chi_wind;
            % eta = asin( - (W/Vg) * sin(chi - chi_wind) )
            % Vgが極小時のゼロ割防止
            sin_eta = -(W_mag / max(1.0, Vg)) * sin(chi_diff);
            sin_eta = max(-0.95, min(0.95, sin_eta));
            eta = asin(sin_eta);
            
            % 2. 偏流角変化率 dot_eta の予測 (Kinematics Prediction)
            % 理論式: dot_eta = (Vg / (Va_horiz * cos(eta)) - 1) * dot_chi
            Va_horiz = max(1.0, V_tas * cos(current_theta));
            
            term_drift = (Vg / (Va_horiz * cos(eta))) - 1.0;
            dot_eta_pred = term_drift * dot_chi_cmd;
            
            % 3. 座標変換 (Coordinate Transform)
            % 理論式: r_req = (cos(theta) / cos(phi)) * (dot_chi - dot_eta)
            % 特異点回避のため phi を制限
            safe_phi = max(-1.0, min(1.0, current_phi)); 
            coord_factor = cos(current_theta) / cos(safe_phi);
            
            % フィードフォワード要求レート
            r_ff = coord_factor * (dot_chi_cmd - dot_eta_pred);
            
            % 4. フィードバック項 (Rate FB - Optional)
            psi_cmd = chi_target_now + eta;
            err_psi = obj.normalize_angle(psi_cmd - psi);
            r_fb = obj.Gains.K_guidance * err_psi;
            
            r_req = r_ff + r_fb;
            
            chi_cmd_out = chi_target_now;
            psi_cmd_out = psi_cmd;
        end
        
        % =========================================================
        % ユーティリティ (変更なし)
        % =========================================================
        function idx = find_closest_index(obj, x, y)
            P = obj.ReferencePathMatrix;
            if isempty(P), idx=1; return; end
            if obj.NeedGlobalSearch
                d_sq = (P(:,1) - x).^2 + (P(:,2) - y).^2;
                [~, idx] = min(d_sq);
                obj.NeedGlobalSearch = false;
            else
                center = obj.LastIndex;
                start_i = max(1, center - 10);
                end_i   = min(size(P,1), center + obj.WindowSize);
                range_idx = start_i : end_i;
                sub_P = P(range_idx, 1:2);
                d_sq = (sub_P(:,1) - x).^2 + (sub_P(:,2) - y).^2;
                [~, loc] = min(d_sq);
                idx = range_idx(loc);
            end
            obj.LastIndex = idx;
        end
        
        function [delta_delta_a, da_ref] = calc_linear_correction(obj, current_state, ref_state)
            % (前回と同じ実装のため省略可だが、クラス完結のため記述)
            p = obj.Params;
            if isfield(p, 'prop'), b=p.prop.b; else, b=p.b; end
            if isfield(p, 'params')
                d=p.params.d; Cnr=p.params.C_n_r; Cnda=p.params.C_n_delta_a;
                Ixx=p.params.I_xx; Iyy=p.params.I_yy;
            else
                d=p.d; Cnr=p.C_n_r; Cnda=p.C_n_delta_a;
                Ixx=p.I_xx; Iyy=p.I_yy;
            end
            S = p.S_c;
            
            u_curr = current_state(1); w_curr = current_state(3);
            phi_curr = current_state(7); z_curr = current_state(12);
            
            V_ref = ref_state.V; phi_ref = ref_state.phi;
            theta_ref = ref_state.theta; alpha_ref = ref_state.alpha;
            u_ref = V_ref * cos(alpha_ref); w_ref = V_ref * sin(alpha_ref);
            g = 9.81;
            
            K_dyn = (g / (V_ref^2)) * cos(theta_ref);
            da_ref = obj.Yaw_Factor * K_dyn * tan(phi_ref); 
            
            psi_dot = (g / V_ref) * tan(phi_ref);
            r_ref_body = psi_dot * cos(phi_ref) * cos(theta_ref);
            p_ref_body = -psi_dot * sin(theta_ref);
            q_ref_body = psi_dot * sin(phi_ref) * cos(theta_ref);
            
            du = u_curr - u_ref; dw = w_curr - w_ref; dphi = phi_curr - phi_ref;
            
            if ~isempty(obj.AtmoModel)
                rho = obj.AtmoModel.get_density(-z_curr / 1000);
            else
                rho = 1.225 * exp(-(-z_curr)/8500);
            end
            
            K_damp = 0.25 * rho * S * b^2 * Cnr;
            K_ctl  = 0.5  * rho * S * b / d * Cnda;
            
            dN_dV = K_damp * r_ref_body + 2 * K_ctl * V_ref * da_ref;
            N_u = dN_dV * (u_ref / V_ref);
            N_w = dN_dV * (w_ref / V_ref);
            N_r = K_damp * V_ref;
            N_q = (Iyy - Ixx) * p_ref_body;
            S_phi = N_r * (-q_ref_body) + N_q * (r_ref_body);
            N_da = K_ctl * V_ref^2;
            
            if abs(N_da) < 1e-9
                delta_delta_a = 0;
            else
                term_vel = N_u * du + N_w * dw;
                term_phi = S_phi * dphi;
                delta_delta_a = - (1 / N_da) * (term_vel + term_phi);
            end
            delta_delta_a = max(-0.5, min(0.5, delta_delta_a));
        end
        
        function a = normalize_angle(~, a)
            a = mod(a + pi, 2*pi) - pi; 
        end
        
        function [dR, dL] = apply_mixing(~, da, ds)
            da = max(-1, min(1, da));
            if da > 0, dR = da + ds; dL = ds;
            else,      dR = ds;      dL = abs(da) + ds;
            end
            dR = max(0, min(1, dR)); dL = max(0, min(1, dL));
        end
        
        function R = rotation_matrix(~, phi, theta, psi)
            ct=cos(theta); st=sin(theta); cp=cos(phi); sp=sin(phi); cpsi=cos(psi); spsi=sin(psi);
            R = [ct*cpsi, sp*st*cpsi-cp*spsi, cp*st*cpsi+sp*spsi;
                 ct*spsi, sp*st*spsi+cp*cpsi, cp*st*spsi-sp*cpsi;
                 -st,     sp*ct,              cp*ct];
        end
    end
end