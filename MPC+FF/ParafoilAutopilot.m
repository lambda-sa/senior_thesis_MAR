classdef ParafoilAutopilot < handle
    % PARAFOILAUTOPILOT (Unified Trajectory Following - Fixed)
    %
    % 概要:
    %   モード分けを廃止し、Plannerが生成した連続軌道データを統合的に追従。
    %   set_nominal_glide_ratio メソッドを追加し、外部からのL/D設定に対応。
    
    properties
        % --- 機体・環境 ---
        Params
        Yaw_Factor
        AtmoModel
        Linearizer
        Dynamics
        
        % --- ゲイン ---
        Gains
        
        % --- 状態管理 ---
        ReferencePathMatrix = []; % [x, y, z, t, psi, kappa, v_ref]
        
        LastIndex = 1;
        WindowSize = 100;
        NeedGlobalSearch = true;
        
        PsiErrIntegral = 0; 
        LastUpdateTime = -1;
        UpdateInterval = 0.05; 
        
        LastDeltaR = 0;     
        LastDeltaL = 0;     
        IsReady = false;    
        
        % --- 制御フラグ ---
        EnableLongitudinalControl = true;
        
        % --- 環境・物理パラメータ ---
        BaseWindVector = [0; 0]; % 計画時の風
        NominalGlideRatio = 3.0; % ★追加: 基準滑空比

        % ▼▼▼▼▼▼ 追加・変更 ▼▼▼▼▼▼
        % 先読み時間 [秒] (推奨: 2.0 ~ 3.0)
        LookAheadTime = 0.0;       
        
        % 操舵反転フラグ (false:通常, true:反転)
        ReverseControlInput = false; 
        % ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
        
    end
    
    methods
        % =========================================================
        % コンストラクタ
        % =========================================================
        function obj = ParafoilAutopilot(params, dynamics_instance, linearizer)
            obj.Params = params;
            if nargin > 1, obj.Dynamics = dynamics_instance; end
            if nargin > 2, obj.Linearizer = linearizer; end
            
            try
                obj.AtmoModel = AtmoTempPressRho();
            catch
                obj.AtmoModel = [];
            end
            
            % --- ゲイン設定 ---
            obj.Gains.k_vf       = 0.06;   % 統合VFゲイン
            obj.Gains.chi_inf    = pi/2;
            obj.Gains.K_guidance = 0.0;    % レートFB
            obj.Gains.K_linear   = 1.0;    % 線形化FF有効
            
            obj.Gains.kp_phi = 0.0;
            obj.Gains.kp_psi = 1.2; 
            obj.Gains.ki_psi = 0.0; 
            obj.Gains.Kp_alt = 1.0;        % 高度制御ゲイン
            
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
        % 公開メソッド
        % =========================================================
        function set_longitudinal_control(obj, enable_flag)
            obj.EnableLongitudinalControl = logical(enable_flag);
            fprintf('Autopilot: Longitudinal Control set to %d\n', obj.EnableLongitudinalControl);
        end
        
        % ★追加: 実行スクリプトのエラー修正用メソッド
        function set_nominal_glide_ratio(obj, val)
            obj.NominalGlideRatio = val;
            fprintf('Autopilot: Nominal L/D set to %.2f\n', val);
        end
        
        function set_target_path(obj, path_matrix)
            if size(path_matrix, 2) < 2, return; end
            obj.ReferencePathMatrix = path_matrix;
            obj.LastIndex = 1;
            obj.NeedGlobalSearch = true;
            fprintf('Autopilot: Path Set (Points: %d)\n', size(path_matrix, 1));
        end
        
        function import_mission_data(obj, planner)
            % Plannerから全フェーズのデータを結合して取り込む
            if isprop(planner, 'ResultData') && ~isempty(planner.ResultData)
                d = planner.ResultData;
                
                fields = {'runup', 'entry', 'loiter', 'dubins', 'final'};
                dx=[]; dy=[]; dz=[]; dt=[]; dpsi=[]; dk=[];
                
                for i=1:length(fields)
                    f = fields{i};
                    if isfield(d, f) && ~isempty(d.(f).x)
                        dx = [dx; d.(f).x(:)];
                        dy = [dy; d.(f).y(:)];
                        dz = [dz; d.(f).z(:)];
                        dt = [dt; d.(f).t(:)];
                        dpsi = [dpsi; d.(f).psi(:)];
                        
                        if isfield(d.(f), 'kappa')
                            dk = [dk; d.(f).kappa(:)];
                        else
                            dk = [dk; zeros(size(d.(f).x(:)))];
                        end
                    end
                end
                
                if isempty(dx)
                    warning('Autopilot: Imported path is empty.');
                    return;
                end
                
                dv = zeros(size(dx));
                path_mat = [dx, dy, dz, dt, dpsi, dk, dv];
                
                [~, unique_idx] = unique(path_mat(:,4), 'stable');
                path_mat = path_mat(unique_idx, :);
                
                obj.set_target_path(path_mat);
            end
            
            if isprop(planner, 'WindVector')
                w = planner.WindVector;
                obj.BaseWindVector = [w(1); w(2)];
                fprintf('Autopilot: Imported Base Wind = [%.1f, %.1f]\n', w(1), w(2));
            end
        end
        
        % =========================================================
        % ★★★ メイン更新ループ (Update) ★★★
        % =========================================================
        function [delta_R, delta_L, log_struct] = update(obj, t, current_state, wind_vec_est, ~)
            
            if ~obj.IsReady || isempty(obj.ReferencePathMatrix)
                delta_R=0; delta_L=0; log_struct=[]; return;
            end
            
            dt = obj.UpdateInterval;
            
            % --- 状態量の展開 ---
            u=current_state(1); v=current_state(2); w=current_state(3);
            phi=current_state(7); theta=current_state(8); psi=current_state(9);
            pos_ned = current_state(10:12);
            z_curr = pos_ned(3); 
            
            V_tas = sqrt(u^2+v^2+w^2);
            V_horiz = max(0.1, V_tas * cos(theta)); 
            
            R_b2n = obj.rotation_matrix(phi, theta, psi);
            Vel_ned = R_b2n * [u;v;w];
            Vg_n = Vel_ned(1); Vg_e = Vel_ned(2);
            Vg_horiz = sqrt(Vg_n^2 + Vg_e^2);
            chi_curr = atan2(Vg_e, Vg_n);
            
            Wx_est = wind_vec_est(1); Wy_est = wind_vec_est(2);
            Wx_base = obj.BaseWindVector(1); Wy_base = obj.BaseWindVector(2);
            wind_disturbance = [Wx_est - Wx_base; Wy_est - Wy_base];
            
            % =====================================================
            % Step 1: Geometry & Kinematics (統合誘導則)
            % =====================================================
            
            % 1. 現在位置の検索
            idx = obj.find_closest_index(pos_ned(1), pos_ned(2));
            P = obj.ReferencePathMatrix;
            
            % 現在位置に対応する参照値 (位置誤差計算用)
            path_x   = P(idx, 1);
            path_y   = P(idx, 2);
            z_ref    = P(idx, 3);
            path_psi = P(idx, 5);
            
            % ▼▼▼▼▼▼ 先読みロジック (Lookahead) に差し替え ▼▼▼▼▼▼
            if obj.LookAheadTime > 0
                % 未来の距離を計算 (速度 x 先読み時間)
                dist_lookahead = Vg_horiz * obj.LookAheadTime;
                
                idx_future = idx;
                dist_acc = 0;
                
                % 経路に沿って未来のインデックスを探す
                while dist_acc < dist_lookahead && idx_future < size(P, 1) - 1
                    d_step = norm(P(idx_future+1, 1:2) - P(idx_future, 1:2));
                    if d_step < 1e-3, d_step = 0.5; end % 安全策
                    dist_acc = dist_acc + d_step;
                    idx_future = idx_future + 1;
                end
                
                % ★未来の点の曲率を採用する
                kappa_ref = P(idx_future, 6); 
            else
                % 先読みオフ
                kappa_ref = P(idx, 6); 
            end
            % ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
            v_ref    = P(idx, 7);
            
            dx = pos_ned(1) - path_x; 
            dy = pos_ned(2) - path_y;
            e_cross = -sin(path_psi)*dx + cos(path_psi)*dy;
            
            chi_rel = obj.normalize_angle(chi_curr - path_psi);
            
            % (A) 経路追従項
            denom = 1.0 - kappa_ref * e_cross;
            if abs(denom) < 0.1, denom = sign(denom) * 0.1; end
            scale_factor = 1.0 / denom;
            
            dot_chi_path = Vg_horiz * kappa_ref * scale_factor * cos(chi_rel);
            
            % (B) 誤差収束項 (VF Term)
            k_vf = obj.Gains.k_vf;
            chi_inf = obj.Gains.chi_inf;
            dot_e = Vg_horiz * sin(chi_rel);
            
            numerator = -(chi_inf * 2 / pi) * (k_vf * dot_e);
            denominator = 1 + (k_vf * e_cross)^2;
            dot_chi_app = numerator / denominator;
            
            dot_chi_cmd = dot_chi_path + dot_chi_app;
            
            chi_app_val = -(chi_inf * 2 / pi) * atan(k_vf * e_cross);
            current_chi_target = path_psi + chi_app_val;
            
            [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                dot_chi_cmd, Vg_horiz, wind_disturbance, psi, current_chi_target, phi, theta, V_tas);
            
            % =====================================================
            % Step 2: Dynamics & Control (制御則)
            % =====================================================
            
            g = 9.81;
            phi_ref = atan( (V_horiz * r_req) / g );
            phi_ref = max(-0.8, min(0.8, phi_ref)); 
            
            ref_state.V = V_tas;
            ref_state.phi = phi_ref;
            ref_state.theta = theta;
            ref_state.r_req = r_req;
            if abs(u) > 0.1, ref_state.alpha = atan2(w, u); else, ref_state.alpha = 0; end
            
            [delta_delta_a, da_ref] = obj.calc_linear_correction(current_state, ref_state);
            dda_term = obj.Gains.K_linear * delta_delta_a;
            
            err_psi = obj.normalize_angle(cmd_psi - psi);
            obj.PsiErrIntegral = max(-0.3, min(0.3, obj.PsiErrIntegral + err_psi * dt));
            
            phi_pid = obj.Gains.kp_psi * err_psi + obj.Gains.ki_psi * obj.PsiErrIntegral;
            phi_cmd_final = phi_ref + phi_pid;
            
            da_fb = obj.Gains.kp_phi * (phi_cmd_final - phi);
            
            da_total = da_ref + dda_term + da_fb;
            % ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
            % ★★★ 追加: 操舵反転ロジック (ここが抜けていました) ★★★
            % ▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼▼
            if obj.ReverseControlInput
                da_total = -da_total;
            end
            % ▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲▲
            % =====================================================
            % Step 3: Longitudinal Control (縦制御)
            % =====================================================
            delta_s = 0;
            
            if obj.EnableLongitudinalControl
                [delta_s, ~, ~] = obj.calc_strict_height_control(...
                    current_state, z_ref, z_curr, phi_ref, v_ref, obj.Gains.Kp_alt);
            end     
            
            max_ds = 1.0 - abs(da_total);
            delta_s = max(0.0, min(max_ds, delta_s));
            
            [delta_R, delta_L] = obj.apply_mixing(da_total, delta_s);
            
            obj.LastUpdateTime = t;
            obj.LastDeltaR = delta_R; 
            obj.LastDeltaL = delta_L;
            
            log_struct.chi_cmd = cmd_chi;
            log_struct.psi_cmd = cmd_psi;
            log_struct.r_req   = r_req;
            log_struct.phi_ref = phi_ref;
            log_struct.da_total= da_total;
            log_struct.z_ref   = z_ref;
            log_struct.idx     = idx;
        end
    end
    
    methods (Access = private)
        
        function [r_req, chi_cmd_out, psi_cmd_out] = apply_wind_and_feedback(obj, ...
                dot_chi_cmd, Vg, wind_vec, psi, chi_target_now, current_phi, current_theta, V_tas)
            
            Wx = wind_vec(1); Wy = wind_vec(2);
            W_mag = norm([Wx, Wy]);
            chi_wind = atan2(Wy, Wx);
            
            chi_diff = chi_target_now - chi_wind;
            sin_eta = -(W_mag / max(1.0, Vg)) * sin(chi_diff);
            sin_eta = max(-0.95, min(0.95, sin_eta));
            eta = asin(sin_eta);
            
            Va_horiz = max(1.0, V_tas * cos(current_theta));
            term_drift = (Vg / (Va_horiz * cos(eta))) - 1.0;
            dot_eta_pred = term_drift * dot_chi_cmd;
            
            safe_phi = max(-1.0, min(1.0, current_phi)); 
            coord_factor = cos(current_theta) / cos(safe_phi);
            
            r_ff = coord_factor * (dot_chi_cmd - dot_eta_pred);
            
            psi_cmd = chi_target_now + eta;
            err_psi = obj.normalize_angle(psi_cmd - psi);
            r_fb = obj.Gains.K_guidance * err_psi;
            
            r_req = r_ff + r_fb;
            chi_cmd_out = chi_target_now;
            psi_cmd_out = psi_cmd;
        end
        
        function [delta_delta_a, da_ref] = calc_linear_correction(obj, current_state, ref_state)
            p = obj.Params;
            if isfield(p, 'params'), par=p.params; else, par=p; end
            if isfield(p, 'prop'), b=p.prop.b; else, b=p.b; end
            d=par.d; Cnr=par.C_n_r; Cnda=par.C_n_delta_a;
            Ixx=par.I_xx; Iyy=par.I_yy; S=p.S_c;
            
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
            
            if ~isempty(obj.AtmoModel), rho=obj.AtmoModel.get_density(-z_curr/1000); else, rho=1.225; end
            
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
        
        function [ds_cmd, G_brake, Zw_total] = calc_strict_height_control(obj, state, z_ref, z_curr, phi_trim, v_ref, Kp)
            u = state(1); v = state(2); w = state(3);
            phi = state(7); theta = state(8); z_ned = state(12);
            p = obj.Params;
            if isfield(p, 'prop'), S=p.prop.S; m=p.prop.TotalMass; else, S=p.S_c; m=1.0; end
            if isfield(p, 'params'), par=p.params; else, par=p; end
            g = 9.81;
            
            if ~isempty(obj.AtmoModel), rho=obj.AtmoModel.get_density(-z_ned/1000); else, rho=1.225; end
            
            V_sq = u^2 + v^2 + w^2;
            if abs(u) < 0.1, u = 0.1; end
            alpha = atan2(w, u); c_a = cos(alpha); s_a = sin(alpha);
            
            CL = par.C_L_0 + par.C_L_alpha * alpha;
            dCL_da = par.C_L_alpha;
            CD0 = par.C_D_0;
            if isfield(par, 'C_D_alpha_sq')
                CD = CD0 + par.C_D_alpha_sq * alpha^2;
                dCD_da = 2 * par.C_D_alpha_sq * alpha; 
            else
                CD = CD0 + par.C_D_alpha * alpha;
                dCD_da = par.C_D_alpha;
            end
            CL_ds = par.C_L_delta_s; CD_ds = par.C_D_delta_s;
            
            Zw1 = - (rho * w * S / m) * (CL * c_a + CD * s_a);
            term_sideslip = 1.0 + (v^2) / (u^2 + w^2);
            dCz_da_part = (dCL_da + CD) * c_a + (dCD_da - CL) * s_a;
            Zw2 = - (rho * u * S / (2*m)) * term_sideslip * dCz_da_part;
            Zw_total = Zw1 + Zw2; if abs(Zw_total) < 1e-6, Zw_total = -1e-6; end
            
            Z_ds = - (rho * V_sq * S / (2*m)) * (CL_ds * c_a + CD_ds * s_a);
            
            cp = cos(phi); ct = cos(theta); sp = sin(phi);
            G_brake = (cp * ct) * (Z_ds / Zw_total);
            if abs(G_brake) < 1e-3, G_brake = 1e-3; end
            
            term_gravity = (m * g * cp * ct) / abs(Zw_total);
            G_bank = sp * ct * (term_gravity + w);
            
            delta_z = z_curr - z_ref;
            delta_phi = phi - phi_trim;
            
            u_input = Kp * delta_z + G_bank * delta_phi;
            ds_cmd = u_input / G_brake;
            ds_cmd = max(0.0, min(1.0, ds_cmd));
        end
        
        function idx = find_closest_index(obj, x, y)
            P = obj.ReferencePathMatrix;
            if isempty(P), idx=1; return; end
            if obj.NeedGlobalSearch
                d_sq = (P(:,1) - x).^2 + (P(:,2) - y).^2;
                [~, idx] = min(d_sq);
                obj.NeedGlobalSearch = false;
            else
                center = obj.LastIndex;
                start_i = max(1, center - 20);
                end_i   = min(size(P,1), center + obj.WindowSize);
                range_idx = start_i : end_i;
                sub_P = P(range_idx, 1:2);
                d_sq = (sub_P(:,1) - x).^2 + (sub_P(:,2) - y).^2;
                [~, loc] = min(d_sq);
                idx = range_idx(loc);
            end
            obj.LastIndex = idx;
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