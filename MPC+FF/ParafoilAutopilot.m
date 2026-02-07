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
        % ★追加: 縦制御のON/OFFフラグ
        EnableLongitudinalControl = true;
        % ★追加: Plannerで計算された基準滑空比 (トリムL/D)
        NominalGlideRatio = 3.0; % 初期値 (適当な値、あとで上書きされる)
        % ★追加: 物理モデルのインスタンス (厳密な迎角取得用)
        Dynamics
        % ★追加: 計画時に想定されていた風 (North, East)
        BaseWindVector = [0; 0];

    end
    
    methods
        % ★★★ ここに追加してください (public methods) ★★★
        function set_longitudinal_control(obj, enable_flag)
            obj.EnableLongitudinalControl = logical(enable_flag);
            fprintf('Autopilot: Longitudinal Control set to %d\n', obj.EnableLongitudinalControl);
        end
        % ★★★ ここまで ★★★

        % =========================================================
        % コンストラクタ
        % =========================================================
        function obj = ParafoilAutopilot(params,dynamics_instance, linearizer)
            obj.Params = params;
            % ★追加: Dynamicsを受け取る
            if nargin > 1
                obj.Dynamics = dynamics_instance;
            end
            if nargin > 2, obj.Linearizer = linearizer; end
            
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
            % ★追加: 高度制御ゲイン (フゴイド時定数相当: g/V ~ 0.5 - 1.0)
            obj.Gains.Kp_alt = 1.0;
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

                    % ★追加: 目標横滑り速度 v_ref (なければ 0)
                    % Plannerが 'v' または 'v_ref' を出力していると仮定
                    if isfield(d.dubins, 'v_ref')
                        dv = [d.dubins.v_ref(:); d.final.v_ref(:)];
                    elseif isfield(d.dubins, 'v')
                        dv = [d.dubins.v(:); d.final.v(:)];
                    else
                        dv = zeros(size(dx));
                    end
                    % 6列構成でセット
                    obj.set_target_path([dx, dy, dz, dt, dpsi, dk, dv]);
                end
            end

            % ★追加: Plannerが持っている風速情報(計画風)を取り込む
            if isprop(planner, 'WindVector')
                w = planner.WindVector;
                obj.BaseWindVector = [w(1); w(2)];
                fprintf('Autopilot: Imported Base Wind = [%.1f, %.1f]\n', w(1), w(2));
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

        % ★追加: 基準滑空比をセットするメソッド
        function set_nominal_glide_ratio(obj, val)
            obj.NominalGlideRatio = val;
            fprintf('Autopilot: Nominal L/D set to %.2f\n', val);
        end
        
        % =========================================================
        % ★★★ メイン更新ループ (Update) ★★★
        % =========================================================
        function [delta_R, delta_L, log_struct] = update(obj, t, current_state, wind_vec_est, ~)
            
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
            % ★★★ 【修正】ここを追加してください ★★★
            z_curr = pos_ned(3); % または z_curr = current_state(12);
            % ★★★★★★★★★★★★★★★★★★★★★★★
            V_tas = sqrt(u^2+v^2+w^2);
            V_horiz = max(0.1, V_tas * cos(theta)); % 対気水平速度
            
            % 現在の対地速度・コース角の計算
            R_b2n = obj.rotation_matrix(phi, theta, psi);
            Vel_ned = R_b2n * [u;v;w];
            Vg_n = Vel_ned(1);
            Vg_e = Vel_ned(2) ;
            Vg_horiz = sqrt(Vg_n^2 + Vg_e^2);
            %Vg_vert = Vel_ned(3);
            chi_curr = atan2(Vg_e, Vg_n);
            
            % ★変更: デフォルト値の設定
            z_ref = z_curr; % 値がない場合は現在高度維持
            v_ref = 0;      % 値がない場合は横滑りなし
            r_req = 0; cmd_chi = 0; cmd_psi = 0; % 初期化
            % =====================================================
            % ★修正箇所: 風の差分計算 (Disturbanceのみ抽出)
            % =====================================================
            Wx_est = wind_vec_est(1); Wy_est = wind_vec_est(2);
            Wx_base = obj.BaseWindVector(1); Wy_base = obj.BaseWindVector(2);
            % 誘導則に渡す風速ベクトル = 「現在の推定風」 - 「計画時の風」
            % 計画通りならゼロになり、カニ歩き補正は行われない（軌道形状で対応済みのため）
            wind_for_guidance = [Wx_est - Wx_base; Wy_est - Wy_base];
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

                % ★追加: Loiterパラメータから目標高度を取得
                if isfield(p, 'z_ref')
                    z_ref = p.z_ref;
                else
                    % パラメータにない場合は、Loiter開始時の高度などを保持するか、現在高度を使用
                    z_ref = z_curr; 
                end

                % 2. 目標横滑り v_ref
                if isfield(p, 'v_ref')
                    v_ref = p.v_ref;
                else
                    % 指定がなければゼロ（Coordinated Turn）
                    v_ref = 0;
                end
                
                % 厳密風補正メソッドへ
                % Loiter中はターゲット方位が刻々と変わるため、簡便に chi_curr をターゲットとみなすか
                % または phi_pos + 90deg を使う
                current_chi_target = chi_curr; 
                % ★修正: wind_vec ではなく wind_for_guidance を渡す
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_for_guidance, psi, current_chi_target, phi, theta, V_tas);
                
                
            elseif strcmpi(obj.CurrentMode, 'Mission') && ~isempty(obj.ReferencePathMatrix)
                % === Missionモード (厳密版) ===
                
                idx = obj.find_closest_index(pos_ned(1), pos_ned(2));
                P = obj.ReferencePathMatrix;
                path_x = P(idx, 1); path_y = P(idx, 2); path_psi = P(idx, 5);
                % ★目標高度 z_ref の取得 (3列目)
                z_ref  = P(idx, 3);

                % ★修正: 曲率 kappa の取得
                if size(P, 2) >= 6
                    kappa_ref = P(idx, 6);
                else
                    kappa_ref = 0;
                end

                % ★追加: 目標横滑り v_ref の取得 (7列目)
                if size(P, 2) >= 7
                    v_ref = P(idx, 7);
                else
                    v_ref = 0;
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
                % ★修正: wind_vec ではなく wind_for_guidance を渡す
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_for_guidance, psi, current_chi_target, phi, theta, V_tas);
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

            % -----------------------------------------------------
            % Step 3: 縦制御 (Longitudinal Control) - 厳密版実装
            % -----------------------------------------------------
            delta_s = 0;
            log_G_brake = 0; log_Zw = 0;
            
            if obj.EnableLongitudinalControl
                % phi_trim は、旋回に必要なバンク角 (phi_ref) を使用
                phi_trim = phi_ref;
                
                % ★変更: v_ref を固定の0ではなく、上で取得した変数を使用
                % 厳密メソッドの呼び出し
                [delta_s, log_G_brake, log_Zw] = obj.calc_strict_height_control(...
                    current_state, z_ref, z_curr, phi_trim, v_ref, obj.Gains.Kp_alt);
            end     
                    
            
            % --- 3. 優先順位付きミキシング (修正) ---
            % 横操作で使った残りを上限とする
            max_ds = 1.0 - abs(da_total);
            delta_s = max(0.0, min(max_ds, delta_s));
            
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
            %log_struct.G_ground = G_req;
            %log_struct.Target_CL_CD = Target_CL_CD;
            % ログへの追加
            log_struct.z_ref = z_ref; % ログに残す
            log_struct.v_ref = v_ref; % ログに残す
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
        
        %{ 
        ★完全新規追加
        function ds = calc_longitudinal_control(obj, state, target_CL_CD, da_cmd, wind_vec)
            % 引数: target_CL_CD は補正済みの目標揚抗比
            
            % 1. 厳密な迎角の取得 (Dynamicsを利用)
            if ~isempty(obj.Dynamics)
                if isfield(obj.Params, 'ang_psi_rigging')
                    GAMMA = obj.Params.ang_psi_rigging;
                else
                    GAMMA = 0;
                end
                
                % ★Dynamicsのメソッド呼び出し
                [alpha, ~, ~] = obj.Dynamics.get_aerodynamic_state(state, wind_vec, GAMMA);
            else
                % フォールバック (簡易計算)
                u=state(1); w=state(3);
                alpha = atan2(w, u);
            end
            
            p = obj.Params;
    
            % ベース空力係数 (A, E) の計算
            % ★ポイント: 旋回操作 da_cmd による抗力増・揚力変化を「ベース」に含める
            A = p.C_L_0 + p.C_L_alpha * alpha + p.C_L_delta_a * abs(da_cmd);
            
            % 抗力係数 (2乗則モデルなどを適用)
            if isfield(p, 'C_D_alpha_sq')
                CD_alpha_term = p.C_D_alpha_sq * alpha^2;
            else
                CD_alpha_term = p.C_D_alpha * alpha^2;
            end
            E = p.C_D_0 + CD_alpha_term + p.C_D_delta_a * abs(da_cmd);
    
            % 制御感度 (B, F)
            B = p.C_L_delta_s;
            F = p.C_D_delta_s;
            target_G =target_CL_CD;
            % 解析解: ds = (A - G*E) / (G*F - B)
            numerator = A - target_G * E;
            denominator = target_G * F - B;
    
            if abs(denominator) < 1e-6
                ds = 0;
            else
                ds = numerator / denominator;
            end
            
            ds = max(0.0, ds);
        end
        %}

        % =========================================================
        % ★★★ 厳密な高度制御メソッド (paramsエラー修正版) ★★★
        % =========================================================
        function [ds_cmd, G_brake, Zw_total] = calc_strict_height_control(obj, state, z_ref, z_curr, phi_trim, v_ref, Kp)
            
            % --- 1. 状態量とパラメータ ---
            u = state(1); v = state(2); w = state(3);
            phi = state(7); theta = state(8);
            z_ned = state(12);
            
            % ★ここが重要: obj.Params を ローカル変数 p に入れる
            p = obj.Params; 
            
            % 物理パラメータ取得 (構造体の階層チェック)
            if isfield(p, 'prop')
                S = p.prop.S; m = p.prop.TotalMass;
            else
                S = p.S_c; m = 1.0; 
            end
            g = 9.81;
            
            % 大気密度
            if ~isempty(obj.AtmoModel)
                rho = obj.AtmoModel.get_density(-z_ned/1000);
            else
                rho = 1.225;
            end
            
            % 速度・迎角
            V_sq = u^2 + v^2 + w^2;
            if abs(u) < 0.1, u = 0.1; end
            alpha = atan2(w, u);
            c_a = cos(alpha); s_a = sin(alpha);
            
            % --- 2. 空力係数 (階層構造に対応) ---
            % p.params がある場合と、p直下にある場合の両方に対応
            if isfield(p, 'params') && isfield(p.params, 'C_L_alpha')
                par = p.params;
            else
                par = p;
            end
            
            CL = par.C_L_0 + par.C_L_alpha * alpha;
            dCL_da = par.C_L_alpha;
            
            CD0 = par.C_D_0;
            if isfield(par, 'C_D_alpha_sq')
                CDa2 = par.C_D_alpha_sq;
                CD = CD0 + CDa2 * alpha^2;
                dCD_da = 2 * CDa2 * alpha; 
            else
                CDa = par.C_D_alpha;
                CD = CD0 + CDa * alpha;
                dCD_da = CDa;
            end
            
            CL_ds = par.C_L_delta_s;
            CD_ds = par.C_D_delta_s;
            
            % --- 3. 厳密な微係数の導出 ---
            
            % Z_w 第一項 (動圧変化項)
            Zw1 = - (rho * w * S / m) * (CL * c_a + CD * s_a);
            
            % Z_w 第二項 (迎角変化項)
            term_sideslip = 1.0 + (v^2) / (u^2 + w^2);
            dCz_da_part = (dCL_da + CD) * c_a + (dCD_da - CL) * s_a;
            
            Zw2 = - (rho * u * S / (2*m)) * term_sideslip * dCz_da_part;
            
            Zw_total = Zw1 + Zw2;
            if abs(Zw_total) < 1e-6, Zw_total = -1e-6; end
            
            % Z_delta_s
            Z_ds = - (rho * V_sq * S / (2*m)) * (CL_ds * c_a + CD_ds * s_a);
            
            % --- 4. ゲイン計算 ---
            cp = cos(phi); ct = cos(theta); sp = sin(phi);
            
            % G_brake
            G_brake = (cp * ct) * (Z_ds / Zw_total);
            if abs(G_brake) < 1e-3, G_brake = 1e-3; end
            
            % G_bank
            term_gravity = (m * g * cp * ct) / abs(Zw_total);
            G_bank = sp * ct * (term_gravity + w);
            
            % G_slip
            G_slip = sp * ct;
            
            % --- 5. 制御則 ---
            delta_z = z_curr - z_ref;
            delta_phi = phi - phi_trim;
            delta_v   = v - v_ref;
            
            %u_input = Kp * delta_z + G_bank * delta_phi + G_slip * delta_v;
            u_input = Kp * delta_z + G_bank * delta_phi;
            ds_cmd = u_input / G_brake;
            ds_cmd = max(0.0, min(1.0, ds_cmd));
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