classdef ParafoilAutopilot < handle
    % PARAFOILAUTOPILOT
    % ベクトル場誘導 (Vector Field Guidance) と 解析的線形化補正 (Analytical Linearization)
    % を融合した、風外乱に強い高精度パラフォイル自律制御クラス。
    
    properties
        % --- 機体パラメータ ---
        Params       % 質量、空力係数、幾何形状などを含む構造体
        Yaw_Factor   % [FF用] バンク角 -> トグル操作量 への変換係数 (簡易モデル)
        AtmoModel    % ★追加: 大気モデル (AtmoTempPressRho)
        
        % --- 制御ゲイン (要調整) ---
        Gains
        % .k_vf_loiter, .k_vf_mission, .chi_inf
        % .K_linear, .kp_psi, .ki_psi, .kd_psi, .kp_phi
        
        % --- 状態・軌道管理 ---
        CurrentMode = 'Loiter';
        LoiterParams
        MissionPath
        
        LastIndex = 1;      % Mission用インデックス
        WindowSize = 100;
        
        % --- 内部状態 ---
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
        function obj = ParafoilAutopilot(params)
            obj.Params = params;
            
            % --- ★追加: 大気モデルの初期化 ---
            % ParafoilControlMapper_linear_old と同様に外部クラスを使用
            try
                obj.AtmoModel = AtmoTempPressRho();
            catch
                warning('AtmoTempPressRho class not found. Falling back to simple model.');
                obj.AtmoModel = [];
            end
            
            % --- ゲイン設定 ---
            obj.Gains.k_vf_loiter  = 0.04;
            obj.Gains.k_vf_mission = 0.1;
            obj.Gains.chi_inf      = pi/2;
            obj.Gains.K_linear     = 1.0; 
            
            obj.Gains.kp_psi = 1.0;
            obj.Gains.ki_psi = 0.1;
            obj.Gains.kd_psi = 0.5;
            obj.Gains.kp_phi = 1.0;
            
            % --- Yaw_Factor計算 ---
            if isfield(params, 'prop'), b=params.prop.b; else, b=params.b; end
            if isfield(params, 'params')
                d=params.params.d; Cnr=params.params.C_n_r; Cnda=params.params.C_n_delta_a;
            else
                d=params.d; Cnr=params.C_n_r; Cnda=params.C_n_delta_a;
            end
            obj.Yaw_Factor = -d * b * Cnr / (2 * Cnda);
            
            obj.IsReady = true;
        end
        
        % =========================================================
        % 軌道データの取り込み
        % =========================================================
        function import_mission_data(obj, planner)
            d = planner.ResultData;
            if isfield(d, 'loiter') && length(d.loiter.x) > 5
                lx=d.loiter.x; ly=d.loiter.y;
                xc = (min(lx)+max(lx))/2; yc = (min(ly)+max(ly))/2;
                R = (max(lx)-min(lx))/2;
                ang1 = atan2(ly(1)-yc, lx(1)-xc); 
                ang2 = atan2(ly(5)-yc, lx(5)-xc);
                lambda = sign(obj.normalize_angle(ang2 - ang1));
                obj.LoiterParams = struct('xc',xc, 'yc',yc, 'R',R, 'lambda',lambda);
            end
            if isfield(d, 'final')
                obj.MissionPath = d.final; 
                obj.MissionPath.num_points = length(obj.MissionPath.x);
            end
            obj.LastIndex = 1;
        end
        
        function set_mode(obj, mode_str)
            if any(strcmpi(mode_str, {'Loiter', 'Mission'}))
                obj.CurrentMode = mode_str;
                obj.PsiErrIntegral = 0;
            end
        end

        % =========================================================
        % ★★★ メイン更新ループ ★★★
        % =========================================================
        function [delta_R, delta_L, log_struct] = update(obj, t, current_state, wind_vec, V_horiz, delta_s_bias)
            
            if (obj.LastUpdateTime >= 0) && (t - obj.LastUpdateTime < obj.UpdateInterval)
                delta_R = obj.LastDeltaR; delta_L = obj.LastDeltaL; log_struct=[]; return;
            end
            
            dt = obj.UpdateInterval;
            pos_ned = current_state(10:12);
            psi = current_state(9);
            
            if strcmpi(obj.CurrentMode, 'Mission') && ~isempty(obj.MissionPath)
                obj.LastIndex = obj.find_closest_index(pos_ned(1), pos_ned(2));
            end
            
            % --- Step 1 & 2: VF誘導と空間微分 ---
            cmd_now = obj.run_vf_guidance(pos_ned, V_horiz, wind_vec);
            
            dt_lookahead = 0.5;
            u_g = V_horiz * cos(psi); v_g = V_horiz * sin(psi);
            pos_next = pos_ned + [u_g; v_g; 0] * dt_lookahead;
            cmd_next = obj.run_vf_guidance(pos_next, V_horiz, wind_vec);
            
            d_psi = obj.normalize_angle(cmd_next.psi_cmd - cmd_now.psi_cmd);
            r_req = d_psi / dt_lookahead;
            
            % --- Step 3: 動的参照生成 ---
            g = 9.81;
            phi_ref = atan( (V_horiz * r_req) / g );
            phi_ref = max(-0.8, min(0.8, phi_ref));
            
            ref_state.V = sqrt(sum(current_state(1:3).^2));
            ref_state.phi = phi_ref;
            ref_state.theta = current_state(8);
            ref_state.r_req = r_req;
            if current_state(1) > 0.1
                ref_state.alpha = atan2(current_state(3), current_state(1));
            else
                ref_state.alpha = 0;
            end
            
            % --- Step 4: 解析的線形化補正 ---
            [delta_delta_a, da_ref] = obj.calc_linear_correction(current_state, ref_state);
            dda_term = obj.Gains.K_linear * delta_delta_a;
            
            % --- Step 5: 統合とPID補正 ---
            psi_err = obj.normalize_angle(cmd_now.psi_cmd - psi);
            obj.PsiErrIntegral = max(-0.3, min(0.3, obj.PsiErrIntegral + psi_err * dt));
            
            phi_pid = obj.Gains.kp_psi * psi_err ...
                    + obj.Gains.ki_psi * obj.PsiErrIntegral ...
                    - obj.Gains.kd_psi * (current_state(6) - r_req);
            
            phi_cmd_final = phi_ref + phi_pid;
            
            phi_curr = current_state(7);
            da_fb = obj.Gains.kp_phi * (phi_cmd_final - phi_curr);
            
            da_total = da_ref + dda_term + da_fb;
            [delta_R, delta_L] = obj.apply_mixing(da_total, delta_s_bias);
            
            % ログ
            obj.LastUpdateTime = t;
            obj.LastDeltaR = delta_R; obj.LastDeltaL = delta_L;
            %log_struct.chi_cmd = cmd_now.chi_cmd;
            log_struct.psi_cmd = cmd_now.psi_cmd;
            log_struct.r_req   = r_req;
            log_struct.phi_ref = phi_ref;
            log_struct.da_ref  = da_ref;
            log_struct.dda     = dda_term;
            log_struct.da_total= da_total;
        end
    end
    
    methods (Access = private)
        
        % =========================================================
        % VF誘導 (V_horiz風補正)
        % =========================================================
        function cmd = run_vf_guidance(obj, pos_ned, V_horiz, wind_vec)
            cmd.chi_cmd = 0; cmd.psi_cmd = 0;
            x = pos_ned(1); y = pos_ned(2);
            
            if strcmpi(obj.CurrentMode, 'Loiter') && ~isempty(obj.LoiterParams)
                p = obj.LoiterParams;
                dx = x - p.xc; dy = y - p.yc;
                dist = sqrt(dx^2+dy^2);
                phi_pos = atan2(dy, dx);
                err = dist - p.R;
                chi_raw = phi_pos + p.lambda*(pi/2) + p.lambda * atan(obj.Gains.k_vf_loiter * err);
            elseif strcmpi(obj.CurrentMode, 'Mission') && ~isempty(obj.MissionPath)
                idx = obj.find_closest_index(x, y);
                P = obj.MissionPath;
                path_x = P.x(idx); path_y = P.y(idx); chi_path = P.chi(idx);
                dx = x - path_x; dy = y - path_y;
                e_cross = -sin(chi_path)*dx + cos(chi_path)*dy;
                chi_raw = chi_path - obj.Gains.chi_inf * (2/pi) * atan(obj.Gains.k_vf_mission * e_cross);
            else
                chi_raw = 0;
            end
            
            chi_cmd = obj.normalize_angle(chi_raw);
            cmd.chi_cmd = chi_cmd;
            
            if V_horiz > 1.0
                Wx = wind_vec(1); Wy = wind_vec(2);
                chi_wind = atan2(Wy, Wx);
                W_mag = norm([Wx, Wy]);
                chi_diff = chi_cmd - chi_wind;
                sin_eta = -(W_mag / V_horiz) * sin(chi_diff);
                sin_eta = max(-0.95, min(0.95, sin_eta));
                eta = asin(sin_eta);
                psi_cmd = chi_cmd + eta;
            else
                psi_cmd = chi_cmd;
            end
            cmd.psi_cmd = obj.normalize_angle(psi_cmd);
        end
        
        % =========================================================
        % 解析的線形化補正
        % =========================================================
        function [delta_delta_a, da_ref] = calc_linear_correction(obj, current_state, ref_state)
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
            
            % --- ★修正: 高度から密度への変換 (Old Mapper準拠) ---
            % NED座標系 z_curr (Down) -> 高度(km) = -z_curr / 1000
            if ~isempty(obj.AtmoModel)
                rho = obj.AtmoModel.get_density(-z_curr / 1000);
            else
                % クラスが無い場合のフォールバック (高度zはm単位)
                rho = 1.225 * exp(-(-z_curr)/8500); 
            end
            
            K_damp = 0.25 * rho * S * b^2 * Cnr;
            K_ctl  = 0.5  * rho * S * b / d * Cnda;
            
            dN_dV = K_damp * r_ref_body + 2 * K_ctl * V_ref * da_ref;
            N_u = dN_dV * (u_ref / V_ref);
            N_w = dN_dV * (w_ref / V_ref);
            
            N_r = K_damp * V_ref;
            N_q = (Iyy - Ixx) * p_ref_body;
            N_phi = 0;
            S_phi = N_r * (-q_ref_body) + N_q * (r_ref_body) + N_phi;
            
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
        
        % --- ヘルパー ---
        function idx = find_closest_index(obj, x, y)
             P = obj.MissionPath;
             center = obj.LastIndex;
             range = max(1, center-10) : min(P.num_points, center+obj.WindowSize);
             d_sq = (P.x(range)-x).^2 + (P.y(range)-y).^2;
             [~, loc] = min(d_sq);
             idx = range(loc);
        end
        function a = normalize_angle(~, a), a = mod(a + pi, 2*pi) - pi; end
        function [dR, dL] = apply_mixing(~, da, ds)
            da = max(-1, min(1, da));
            if da>0, dR=da+ds; dL=ds; else, dR=ds; dL=abs(da)+ds; end
            dR=max(0,min(1,dR)); dL=max(0,min(1,dL));
        end
    end
end