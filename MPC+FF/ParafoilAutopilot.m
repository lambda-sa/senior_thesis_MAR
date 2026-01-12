classdef ParafoilAutopilot < handle
    % PARAFOILAUTOPILOT (Analytical VF Guidance Version)
    %
    % 概要:
    %   パラフォイルの自律誘導制御を行うクラス。
    %   Loiter/Mission双方において「解析的微分（Analytical Differentiation）」
    %   を用いたベクトル場誘導を実装し、遅れのない高精度な追従を実現する。
    
    properties
        % =========================================================
        % 1. 機体・環境パラメータ
        % =========================================================
        Params       % 質量、翼面積、慣性モーメント、空力係数などを含む構造体
        Yaw_Factor   % [FF用] バンク角(phi) -> 左右トグル差(delta_a) への変換係数
        AtmoModel    % 大気モデル (高度に応じた密度計算用)
        Linearizer   % (オプション) 外部線形化クラス
        
        % =========================================================
        % 2. 制御ゲイン
        % =========================================================
        Gains
        % .k_vf_loiter  : Loiter時の経路収束ゲイン
        % .k_vf_mission : Mission時の経路収束ゲイン
        % .chi_inf      : 最大進入角度 [rad] (通常 pi/2)
        % .K_guidance   : 方位誤差フィードバックゲイン (外側ループ)
        % .K_linear     : 線形化補正の適用率 (0.0=OFF, 1.0=Full)
        
        % PIDゲイン (バンク角制御 - 内側ループ)
        % .kp_phi (必須), .kp_psi, .ki_psi (補助)
        
        % =========================================================
        % 3. 誘導・軌道管理
        % =========================================================
        CurrentMode = 'Loiter'; 
        LoiterParams        % .xc, .yc, .R, .lambda
        ReferencePathMatrix = []; % [x, y, z, t, psi, (kappa)]
        
        LastIndex = 1;      % Mission経路探索用
        WindowSize = 100;
        NeedGlobalSearch = true;
        
        % =========================================================
        % 4. 内部状態変数
        % =========================================================
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
            
            % --- ゲイン初期設定 ---
            obj.Gains.k_vf_loiter  = 0.08;
            obj.Gains.k_vf_mission = 0.15;
            obj.Gains.chi_inf      = pi/2;
            
            obj.Gains.K_guidance   = 0.0;    % 方位角FB (レート指令へ加算)
            obj.Gains.K_linear     = 1.0;    % 線形化FF
            
            % 内側ループ (バンク角FB)
            obj.Gains.kp_phi = 0.0; % これが重要 (指令バンク角へ追従させる)
            
            % 補助PID (レートFBを入れたので、こちらのP項は弱めるかゼロに)
            obj.Gains.kp_psi = 0.0; 
            obj.Gains.ki_psi = 0.0; % I項は定常偏差消去用に少し残す
            obj.Gains.kd_psi = 0.0;
            
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
        % データセット用メソッド
        % =========================================================
        function set_target_path(obj, path_matrix)
            if size(path_matrix, 2) < 2, return; end
            obj.ReferencePathMatrix = path_matrix;
            obj.LastIndex = 1;
            obj.NeedGlobalSearch = true;
        end
        
        function import_mission_data(obj, planner)
            % Mission軌道の取り込み
            if ismethod(planner, 'get_unified_trajectory')
                [path_mat, ~] = planner.get_unified_trajectory();
                if ~isempty(path_mat), obj.set_target_path(path_mat); end
            elseif isprop(planner, 'ResultData') && isfield(planner.ResultData, 'dubins')
                d = planner.ResultData;
                if isfield(d, 'dubins') && ~isempty(d.dubins.x)
                    dx = [d.dubins.x(:); d.final.x(:)];
                    dy = [d.dubins.y(:); d.final.y(:)];
                    dz = [d.dubins.z(:); d.final.z(:)];
                    dt = [d.dubins.t(:); d.final.t(:)];
                    dpsi = [d.dubins.psi(:); d.final.psi(:)];
                    
                    % もしPlannerがkappa(6列目)を持っていればそれも取り込む
                    % ここでは標準的な5列構成と仮定
                    obj.set_target_path([dx, dy, dz, dt, dpsi]);
                end
            end
            
            % Loiterパラメータは Scheduler 側で設定されるためここでは最低限
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
            
            % --- 状態量の展開 ---
            u=current_state(1); v=current_state(2); w=current_state(3);
            phi=current_state(7); theta=current_state(8); psi=current_state(9);
            pos_ned = current_state(10:12);
            
            V_tas = sqrt(u^2+v^2+w^2);
            V_horiz = max(0.1, V_tas * cos(theta)); % 対地水平速度の目安
            
            % -----------------------------------------------------
            % Step 1: ベクトル場誘導 & レート要求 (解析的微分)
            % -----------------------------------------------------
            r_req = 0;
            cmd_chi = 0; % ログ用
            cmd_psi = 0; % ログ用
            
            if strcmpi(obj.CurrentMode, 'Loiter') && ~isempty(obj.LoiterParams)
                % === Loiterモード ===
                p = obj.LoiterParams;
                dx = pos_ned(1) - p.xc; dy = pos_ned(2) - p.yc;
                d = sqrt(dx^2 + dy^2);
                phi_pos = atan2(dy, dx);
                
                % (A) 現在の対地コース角(chi)の厳密計算
                % (風がある場合、対気速度ベクトル+風ベクトルから求めるべきだが、
                %  ここでは簡易的に慣性速度から計算するか、前回のコマンドを利用)
                %  今回は機体速度(u,v,w)からNED速度への変換で計算
                R_b2n = obj.rotation_matrix(phi, theta, psi);
                Vel_ned = R_b2n * [u;v;w]; % 風成分を含まない対気速度(NED)
                % ここに風を足して対地速度にする
                Vg_n = Vel_ned(1) + wind_vec(1);
                Vg_e = Vel_ned(2) + wind_vec(2);
                chi_curr = atan2(Vg_e, Vg_n);
                Vg_horiz = sqrt(Vg_n^2 + Vg_e^2); % 真の対地速度
                
                % (B) 解析的レート計算
                ang_diff = chi_curr - phi_pos;
                dot_phi_pos = Vg_horiz * sin(ang_diff) / d; % 幾何学項
                
                k = obj.Gains.k_vf_loiter;
                err = d - p.R;
                dot_d = Vg_horiz * cos(ang_diff);
                dot_chi_app = (p.lambda * k * dot_d) / (1 + (k * err)^2); % 収束項
                
                dot_chi_cmd = dot_phi_pos + dot_chi_app;
                
                % (C) 風補正 & FB
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_vec, psi, chi_curr);
                
            elseif strcmpi(obj.CurrentMode, 'Mission') && ~isempty(obj.ReferencePathMatrix)
                % === Missionモード (解析的微分へ変更) ===
                
                % 1. 経路情報の取得
                idx = obj.find_closest_index(pos_ned(1), pos_ned(2));
                P = obj.ReferencePathMatrix;
                path_x = P(idx, 1);
                path_y = P(idx, 2);
                path_psi = P(idx, 5); % 経路の方位角
                
                % 2. 曲率(kappa)の局所推定
                % ★対策2: 曲率のリミット処理
                kappa = 0;
                if idx < size(P, 1)
                    next_x = P(idx+1, 1); next_y = P(idx+1, 2);
                    dist_seg = sqrt((next_x - path_x)^2 + (next_y - path_y)^2);
                    
                    % 0.5m以下の微小移動はノイズとみなして無視
                    if dist_seg > 0.5
                        d_psi_path = obj.normalize_angle(P(idx+1, 5) - path_psi);
                        raw_k = d_psi_path / dist_seg;
                        
                        % 最大曲率制限 (R_min = 30m)
                        limit_k = 1.0 / 30.0;
                        kappa = max(-limit_k, min(limit_k, raw_k));
                    end
                end
                
                % 3. 誤差計算
                dx = pos_ned(1) - path_x; dy = pos_ned(2) - path_y;
                % クロストラック誤差 (進行方向左を正とする定義)
                e_cross = -sin(path_psi)*dx + cos(path_psi)*dy;
                
                % 現在の対地コース (簡易計算)
                R_b2n = obj.rotation_matrix(phi, theta, psi);
                Vel_ned = R_b2n * [u;v;w];
                Vg_n = Vel_ned(1) + wind_vec(1);
                Vg_e = Vel_ned(2) + wind_vec(2);
                chi_curr = atan2(Vg_e, Vg_n);
                Vg_horiz = sqrt(Vg_n^2 + Vg_e^2);
                
                % 4. レート計算 (解析解)
                % (A) パス追従項
                dot_chi_path = kappa * Vg_horiz;
                
                % (B) 誤差収束項 (d_chi_app/dt)
                k_vf = obj.Gains.k_vf_mission;
                chi_inf = obj.Gains.chi_inf;
                
                % 誤差の変化率 de/dt = Vg * sin(chi_curr - path_psi)
                chi_rel = obj.normalize_angle(chi_curr - path_psi);
                dot_e = Vg_horiz * sin(chi_rel);
                
                % 微分式: -(2/pi)*chi_inf * (k * dot_e) / (1 + (ke)^2)
                numerator = -(chi_inf * 2 / pi) * (k_vf * dot_e);
                denominator = 1 + (k_vf * e_cross)^2;
                dot_chi_app = numerator / denominator;
                
                dot_chi_cmd = dot_chi_path + dot_chi_app;
                
                % (C) 風補正 & FB
                % ここで目標対地コース chi_cmd も再構築して渡す
                % chi_cmd = path_psi + chi_app(e)
                chi_app_val = -(chi_inf * 2 / pi) * atan(k_vf * e_cross);
                current_chi_target = path_psi + chi_app_val;
                
                [r_req, cmd_chi, cmd_psi] = obj.apply_wind_and_feedback(...
                    dot_chi_cmd, Vg_horiz, wind_vec, psi, current_chi_target);
            end
            
            % -----------------------------------------------------
            % Step 2: バンク角指令 (phi_ref) の生成
            % -----------------------------------------------------
            g = 9.81;
            phi_ref = atan( (V_horiz * r_req) / g );
            phi_ref = max(-0.8, min(0.8, phi_ref)); % リミッタ
            
            % -----------------------------------------------------
            % Step 3: 線形化フィードフォワード (Linearization)
            % -----------------------------------------------------
            ref_state.V = V_tas;
            ref_state.phi = phi_ref;
            ref_state.theta = theta;
            ref_state.r_req = r_req;
            if abs(u) > 0.1, ref_state.alpha = atan2(w, u); else, ref_state.alpha = 0; end
            
            [delta_delta_a, da_ref] = obj.calc_linear_correction(current_state, ref_state);
            dda_term = obj.Gains.K_linear * delta_delta_a;
            
            % -----------------------------------------------------
            % Step 4: 残留誤差補正 (PID) & ミキシング
            % -----------------------------------------------------
            % 目標方位との誤差 (I項用)
            err_psi = obj.normalize_angle(cmd_psi - psi);
            obj.PsiErrIntegral = max(-0.3, min(0.3, obj.PsiErrIntegral + err_psi * dt));
            
            % バンク角補正 (PID)
            % ※ K_guidance(Step1)で主な方位補正は済んでいるので、ここは積分項メイン
            phi_pid = obj.Gains.kp_psi * err_psi ...
                    + obj.Gains.ki_psi * obj.PsiErrIntegral;
            
            phi_cmd_final = phi_ref + phi_pid;
            
            % バンク角FB (最重要: 指令通りに傾ける)
            da_fb = obj.Gains.kp_phi * (phi_cmd_final - phi);
            
            % 合計操作量
            da_total = da_ref + dda_term + da_fb;
            [delta_R, delta_L] = obj.apply_mixing(da_total, 0);
            
            % -----------------------------------------------------
            % Step 5: 更新 & ログ
            % -----------------------------------------------------
            obj.LastUpdateTime = t;
            obj.LastDeltaR = delta_R; 
            obj.LastDeltaL = delta_L;
            
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
        % 共通処理: 風補正とフィードバックの適用
        % =========================================================
        function [r_req, chi_cmd_out, psi_cmd_out] = apply_wind_and_feedback(obj, dot_chi_cmd, Vg, wind_vec, psi, chi_target_now)
            
            % 1. 風係数 (Factor) の計算
            if Vg > 1.0
                Wx = wind_vec(1); Wy = wind_vec(2);
                W_mag = norm([Wx, Wy]);
                chi_wind = atan2(Wy, Wx);
                
                % クラブ角 eta の計算
                % eta = asin( - (W/Vg) * sin(chi - chi_wind) )
                chi_diff = chi_target_now - chi_wind;
                sin_eta = -(W_mag / Vg) * sin(chi_diff);
                sin_eta = max(-0.95, min(0.95, sin_eta));
                eta = asin(sin_eta);
                
                % 補正係数
                numer = W_mag * cos(chi_target_now - chi_wind);
                denom = Vg * cos(eta);
                if abs(denom) > 0.1
                    factor = 1 - (numer / denom);
                else
                    factor = 1.0;
                end
                
                r_ff = dot_chi_cmd * factor;
                psi_cmd = chi_target_now + eta;
            else
                r_ff = dot_chi_cmd;
                psi_cmd = chi_target_now;
            end
            
            % 2. フィードバック項 (Rate FB)
            err_psi = obj.normalize_angle(psi_cmd - psi);
            r_fb = obj.Gains.K_guidance * err_psi;
            
            r_req = r_ff + r_fb;
            
            chi_cmd_out = chi_target_now;
            psi_cmd_out = psi_cmd;
        end
        
        % =========================================================
        % ヘルパー: 最近傍点探索
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
        
        % =========================================================
        % ヘルパー: 線形化補正 (前回と同じ)
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
        
        % =========================================================
        % ユーティリティ
        % =========================================================
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