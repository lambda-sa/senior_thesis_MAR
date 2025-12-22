classdef ParafoilDynamics3DOF_Aero_DeltaA_yaw < handle
    
    properties
        params       
        atmo_model   
        
        % 新しいモデル用の定数・係数
        Yaw_Control_Factor  % b / (2 * Cn_delta_a)
        C_n_r
        C_n_p
        
        max_bank_rad
        tau_delta_a
    end
    
    methods
        function obj = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, atmo_model)
            obj.params = params;
            obj.atmo_model = atmo_model;
            
            % --- 制御則用パラメータの抽出 (ヨー軸主導モデル) ---
            % 式: delta_a = (b * g / (2 * V^2 * Cn_da)) * [ ... ]
            % 定数部分 factor = b / (2 * Cn_da) を事前計算
            
            b = params.b;
            d = params.d;

            Cn_delta_a = params.C_n_delta_a; % ※paramsに追加されている必要あり
            
            if abs(Cn_delta_a) < 1e-9
                error('Cn_delta_a is zero. Cannot use Yaw-based control model.'); 
            end
           
            obj.C_n_r = params.C_n_r;
            %obj.Cn_p = params.Cn_p;
            obj.Yaw_Control_Factor = -d*b *obj.C_n_r/ (2 * Cn_delta_a);

            % 制限値と時定数
            obj.max_bank_rad = deg2rad(params.guidance.max_bank_deg);
            
            if isfield(params, 'tau_delta_a')
                obj.tau_delta_a = params.tau_delta_a;
            else
                obj.tau_delta_a = params.tau_sigma; 
            end
        end
        
        function [dy, L, D] = get_derivatives(obj, t, y, control_inputs)
            V     = y(1); 
            gamma = y(2); 
            psi   = y(3); 
            h     = y(6); 
            delta_a_actual = y(7); 
            delta_s_actual = y(8); 
            mu_actual      = y(9); % ★追加: 現在のリギング角 (rad)

            delta_a_cmd = control_inputs.delta_a_cmd; 
            delta_s_cmd = control_inputs.delta_s_cmd; 
            w_x = control_inputs.wind_I(1);
            w_y = control_inputs.wind_I(2);
            mu_cmd = control_inputs.mu_cmd; % ★追加: リギング角指令値
            
            %時間遅れを無視した場合    
            %delta_a_actual = control_inputs.delta_a_cmd;
            %delta_s_actual = control_inputs.delta_s_cmd;
            %mu_actual      = control_inputs.mu_cmd;

            h_km = h / 1000.0; if h_km < 0, h_km = 0; end
            rho_current = obj.atmo_model.get_density(h_km);
            g_current = obj.atmo_model.get_gravity(h);
            
            temp_params = obj.params;
            temp_params.rho = rho_current;
            temp_params.g = g_current;
            
            try
                % alphaを計算 (sigma=0と仮定してもCmには影響しないモデル構成)
                alpha_rad = obj.solve_trim_alpha(temp_params, delta_a_actual,delta_s_actual, mu_actual);
        
            catch
                alpha_rad = 0; 
            end
            % --- 2. 実バンク角 sigma (phi) の計算 ---
            % 現在の速度 V とピッチ角 theta を用いて、釣り合うバンク角を逆算
            
            theta_current = alpha_rad + gamma;
            %fprintf(' (α,θ:%.3f, %.3f)\n', rad2deg(alpha_rad),rad2deg(theta_current));
            sigma_actual = obj.calc_bank_angle_from_delta_a(...
                delta_a_actual, ...
                V, ...
                theta_current, ...
                temp_params, ...
                g_current ...
            );
            
            % --- 空力計算 ---
            [CL_total, C_Dc,C_Ds] = obj.calculate_cl_cd_total(temp_params, alpha_rad,delta_a_actual, delta_s_actual);
            CD_total = C_Dc + C_Ds;
            
            q = 0.5 * rho_current * V^2; 
            S = obj.params.S_c; 
            L = q * S * CL_total; 
            D = q * S * CD_total; 
            W = temp_params.m_total * g_current; 
            
            % --- 運動方程式 ---
            V_dot = -(D + W * sin(gamma)) / temp_params.m_total;
            
            % sigma_actual はここでは物理的なバンク角 phi として扱われます
            gamma_dot = (L * cos(sigma_actual) - W * cos(gamma)) / (temp_params.m_total * V);
            
            cos_gamma = cos(gamma);
            if abs(cos_gamma) < 1e-6, psi_dot = 0.0; 
            else, psi_dot = (L * sin(sigma_actual)) / (temp_params.m_total * V * cos_gamma); end
            
            x_dot = V * cos_gamma * cos(psi) + w_x;
            y_dot = V * cos_gamma * sin(psi) + w_y;
            h_dot = V * sin(gamma);
            
            delta_a_dot = (delta_a_cmd - delta_a_actual) / obj.tau_delta_a;
            delta_s_dot = (delta_s_cmd - delta_s_actual) / obj.params.tau_delta_s;
            % params.tau_mu (時定数) が必要
            mu_dot = (mu_cmd - mu_actual) / obj.params.tau_mu;
            
            % --- 微分ベクトルの結合 ---
            dy = [V_dot; gamma_dot; psi_dot; x_dot; y_dot; h_dot; delta_a_dot; delta_s_dot; mu_dot];
            %dy = [V_dot; gamma_dot; psi_dot; x_dot; y_dot; h_dot; 0; 0; 0];
        end
        
        %% 初期トリム計算 (ロジック変更に伴い修正)
       function [alpha_trim, gamma_trim, V_trim, theta_trim, sigma_trim] = ...
                calculate_initial_trim(obj, h_init, delta_a_init, delta_s_init,mu_val)
            
            h_km = max(0, h_init / 1000.0);
            rho = obj.atmo_model.get_density(h_km);
            g = obj.atmo_model.get_gravity(h_init);
            
            p = obj.params;
            p.rho = rho;
            p.g = g;
            
            % 1. alpha (縦トリム) の決定
            alpha_trim = obj.solve_trim_alpha(p,delta_a_init, delta_s_init,mu_val);
            
            % 2. V, gamma, sigma の連立計算 (反復法)
            %    V <-> Sigma の相互依存を解消して整合性を取る
            
            % 空力係数 (alpha固定)
            %[CL, CDc, CDs] = obj.calculate_cl_cd_total(p, alpha_trim, delta_s_init, delta_a_init);
            [CL, CDc, CDs] = obj.calculate_cl_cd_total(p, alpha_trim,delta_a_init, delta_s_init);
            CD = CDc + CDs;
            
            % 初期推定
            sigma_curr = 0.0;
            V_curr = 15.0; % 仮の初期速度
            
            W = p.m_total * g;
            S = p.S_c;
            
            % 反復ループ
            for iter = 1:20
                % A. 経路角 gamma (釣り合い)
                % tan(gamma) = -CD / (CL * cos(sigma))
                cos_sigma = cos(sigma_curr);
                if abs(cos_sigma) < 0.01, cos_sigma = 0.01; end
                gamma_curr = atan( -CD / (CL * cos_sigma) );
                
                % B. ピッチ角 theta
                theta_curr = alpha_trim + gamma_curr;
                
                % C. 速度 V (揚力釣り合い)
                % L * cos(sigma) = W * cos(gamma)
                lift_req = W * cos(gamma_curr) / cos_sigma;
                V_new = sqrt( (2 * lift_req) / (rho * S * CL) );
                
                % D. バンク角 sigma (制御則逆算)
                % 新しい V と theta を使って再計算
                sigma_new = obj.calc_bank_angle_from_delta_a(...
                    delta_a_init, V_new, theta_curr, p, g);
                
                % バンクリミット
                sigma_new = max(-obj.max_bank_rad, min(obj.max_bank_rad, sigma_new));
                
                % 収束判定
                if abs(sigma_new - sigma_curr) < 1e-4 && abs(V_new - V_curr) < 1e-3
                    sigma_curr = sigma_new;
                    V_curr = V_new;
                    break;
                end
                
                sigma_curr = sigma_new;
                V_curr = V_new;
            end
            
            gamma_trim = gamma_curr;
            V_trim     = V_curr;
            theta_trim = alpha_trim + gamma_trim;
            sigma_trim = sigma_curr;
        end
    
    
            
        %% --- ヘルパー: 制御則逆算 (delta_a -> sigma) ---
        function sigma = calc_bank_angle_from_delta_a(obj, delta_a, V, theta, params, g)
            if V < 0.1 || abs(params.C_n_delta_a) < 1e-9
                sigma = 0.0; return;
            end
            
            % 1. 共通係数 (b*d / 2*Cn_da) * (g / V^2 *cos(theta))
            factor_geo = obj.Yaw_Control_Factor;
            factor_dyn = g / (V^2)* cos(theta);
            
           
            
            % 3. 全体ゲイン K (delta_a = K * sin(sigma))
            K_total = factor_geo * factor_dyn;
            
            % 4. 逆算 (asin)
            if abs(K_total) < 1e-9
                sigma = 0.0;
            else
                sin_val = delta_a / K_total;
                sin_val = max(-1, min(1, sin_val)); % クランプ
                sigma = asin(sin_val);
            end
            
            % リミット適用
            %sigma = max(-obj.max_bank_rad, min(obj.max_bank_rad, sigma));
        end

        %% --- ヘルパー: トリム迎角の計算 ---
        function alpha = solve_trim_alpha(obj, p, delta_a,delta_s, mu_val)
            %disp(['Mu_Cmd Size(初期トリム計算内): ' num2str(size(mu_val))]);
            C_m_func = @(a) obj.calculate_cm_total(p, a, delta_a,delta_s, mu_val);
            alpha_guess = deg2rad(7.0);
            try
                alpha = fzero(C_m_func, alpha_guess, optimset('Display','off'));
                %fprintf(' alpha_trim_ソルバー内": %.2f deg\n', rad2deg(alpha));
            catch ME
                % fzero が失敗した場合 (初期値で符号が同じ、など)
                % ★ 修正: '%' を '%%' にエスケープする安全な形式
                warning_message = sprintf('fzero が解を見つけられませんでした (epsilon=%.2f)。NaNを返します。\nエラー詳細: %s', ...
                                          delta_s, ME.message);
                escaped_message = strrep(warning_message, '%', '%%');
                warning(escaped_message); 
                
                return;
            end

            if ~isnan(alpha)
                delta_alpha = deg2rad(0.001);
                Cm_plus = C_m_func(alpha + delta_alpha);
                Cm_minus = C_m_func(alpha - delta_alpha);
                dCm_dAlpha = (Cm_plus - Cm_minus) / (2 * delta_alpha);
                
                if dCm_dAlpha >= 0
                    % 静的に不安定 (dCm/dAlpha が正またはゼロ)
                    %warning('静的安定性警告: epsilon=%.2f, alpha=%.2f deg で dCm/dAlpha = %.4f (>= 0) です。', ...
                           %delta_s, rad2deg(alpha), dCm_dAlpha);
                end
            end
        end

        %{
        % --- 初期トリム計算用のヘルパー (空力固定版) ---
        function resid = calculate_sigma_residual_with_fixed_aero(obj, sigma, delta_a, p, alpha, CL, CD)
            % 1. sigma から gamma, V を計算
            [gamma, V] = obj.calc_gamma_V_from_fixed_aero(p, sigma, CL, CD);
            theta = alpha + gamma;
            
            % 2. 提案式で delta_a を計算
            common_factor = obj.Yaw_Control_Factor * (p.g / (V^2));
            term_main = -obj.Cn_r * cos(theta) * sin(sigma);
            term_sub = obj.Cn_p * sin(theta) * tan(sigma);
            
            delta_a_calc = common_factor * (term_main + term_sub);
            
            resid = delta_a_calc - delta_a;
        end

        
        % --- 既存のヘルパー関数群 (変更なし) ---
        function [gamma, V] = calc_gamma_V_from_fixed_aero(~, p, sigma, CL, CD)
            if CL <= 1e-3, CL = 1e-3; end
            gamma = atan( -CD / (CL * cos(sigma)) );
            W = p.m_total * p.g;
            aero_norm = sqrt( (CL * cos(sigma))^2 + CD^2 );
            V = sqrt( (2 * W) / (p.rho * p.S_c * aero_norm) );
        end
        
        function [alpha_trim_rad, gamma_trim_rad, V_trim, V_h, V_v] = ...
                calculate_trim_alpha_gamma_multibody(obj, params, epsilon, sigma_rad)
            C_m_total = @(alpha) obj.calculate_cm_total(params, alpha, epsilon);
            try
                alpha_guess = deg2rad(5.0);
                alpha_trim_rad = fzero(C_m_total, alpha_guess, optimset('Display', 'off'));
            catch ME, rethrow(ME); end
            
            [CL_total, C_Dc,C_Ds] = obj.calculate_cl_cd_total(params, alpha_trim_rad, epsilon);
            CD_total = C_Dc + C_Ds;
            if CL_total <= 0, gamma_trim_rad = NaN; V_trim = NaN; V_h = NaN; V_v = NaN; return; end
            gamma_trim_rad = atan(-CD_total / (CL_total*cos(sigma_rad)));
            W = params.m_total * params.g;
            rho = params.rho;
            S_c = params.S_c;
            numerator = 2 * W * cos(gamma_trim_rad);
            denominator = rho * S_c * CL_total * cos(sigma_rad); 
            if denominator <= 0, V_trim = NaN; V_h = NaN; V_v = NaN; return; end
            V_trim = sqrt(numerator / denominator);
            V_h = V_trim * cos(gamma_trim_rad); V_v = V_trim * sin(gamma_trim_rad); 
        end
        
        function [C_Lc,C_Dc,C_Ds] = calculate_cl_cd_total(~,params, alpha, epsilon)
            C_Lc = params.C_L_0 + params.C_L_alpha * alpha + params.C_L_delta_s * epsilon;
            C_Dc = params.C_D_0 + params.C_D_alpha * (alpha^2) + params.C_D_delta_s * epsilon;
            C_Ds = params.C_D_s * (params.S_p / params.S_c);
        end
        function Cm_total = calculate_cm_total(obj, params, alpha, epsilon)
            C_mc = params.C_m_0; 
            C_m_alpha = params.C_m_alpha;
            [C_Lc,C_Dc,C_Ds] = obj.calculate_cl_cd_total(params, alpha, epsilon);
            Rcg_x = params.Rcg_x; Rcg_z = params.Rcg_z;
            Rpg_x = params.Rpg_x; Rpg_z = params.Rpg_z;
            c_ref = params.c;
            r_cg = [Rcg_x;0;Rcg_z]; r_pg = [Rpg_x;0;Rpg_z];
            r_cg_cross = obj.vec_to_skew(r_cg);
            r_pg_cross = obj.vec_to_skew(r_pg);
            turn_alpha_rigging = obj.get_T_AC(alpha + params.ang_psi_rigging+pi);
            Cm_canopy_moment = C_mc +C_m_alpha *alpha; 
            Cm_canopy_forces = r_cg_cross*turn_alpha_rigging*[C_Dc;0;C_Lc]/c_ref;
            Cm_payload_forces = r_pg_cross*turn_alpha_rigging*[C_Ds;0;0]/c_ref;
            Cm_total = Cm_canopy_moment + Cm_canopy_forces(2,1) + Cm_payload_forces(2,1);
        end
        function S = vec_to_skew(~,v)
            S = [ 0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0 ];
        end
        function T_AC = get_T_AC(~,alpha)
            c_a = cos(alpha); s_a = sin(alpha);
            T_AC = [c_a, 0, -s_a; 0, 1, 0; s_a, 0, c_a];
        end
        %}
         %% --- ヘルパー関数: C_L_total, C_D_total の計算 ---
        function [C_Lc,C_Dc,C_Ds] = calculate_cl_cd_total(~,params, alpha,delta_a, epsilon)
            % alpha (rad), epsilon [0-1] を入力
            % システム全体の空力係数を計算 (S_c 基準で正規化)
            
            % 1. キャノピー空力係数
            %    (Excelから読み込んだパラメータ名 C_L_0, C_L_alpha 等を使用)
            C_Lc = params.C_L_0 + params.C_L_alpha * alpha + params.C_L_delta_a * delta_a+ params.C_L_delta_s * epsilon;
            C_Dc = params.C_D_0 + params.C_D_alpha * (alpha^2) +params.C_D_delta_a * delta_a+  params.C_D_delta_s * epsilon;
            
            % 2. ペイロード抗力係数 (S_c 基準に正規化)
            C_Ds = params.C_D_s * (params.S_p / params.S_c);

            %{
            % 3. 機体軸での力係数 (風軸 -> 機体軸)
            % キャノピー
            F_x_c_coeff = -C_Dc * cos(alpha) + C_Lc * sin(alpha);
            F_z_c_coeff = -C_Dc * sin(alpha) - C_Lc * cos(alpha);
            % ペイロード (抗力が -X機体軸方向と仮定)
            F_x_p_coeff = -C_Ds_norm;
            F_z_p_coeff = 0;
            
            % 4. システム全体の機体軸力係数
            F_x_total = F_x_c_coeff + F_x_p_coeff;
            F_z_total = F_z_c_coeff + F_z_p_coeff;
            
            % 5. システム全体の揚力・抗力係数 (機体軸 -> 風軸)
            %    (風軸における L, D が正になるように定義)
            CL_total = F_x_total * sin(alpha) - F_z_total * cos(alpha);
            CD_total = -F_x_total * cos(alpha) - F_z_total * sin(alpha);
            %}
            
        end

        %% --- ヘルパー関数: C_m_total の計算 ---
        function Cm_total = calculate_cm_total(obj, params, alpha, delta_a,epsilon, mu_rigging)
            % alpha (rad), epsilon [0-1] を入力
            
            % 1. キャノピー単体の空力係数 (C_Lc, C_Dc, C_mc) を計算
            %    (Excelから読み込んだパラメータ名 C_L_0, C_L_alpha 等を使用)
            
            C_mc = params.C_m_0; % (Cmc/4 は Cm0 と同じとする)
            C_m_alpha = params.C_m_alpha;
            [C_Lc,C_Dc,C_Ds] = obj.calculate_cl_cd_total(params, alpha,delta_a, epsilon);
            % 2. ペイロードの抗力係数 (C_Ds) を S_c 基準に正規化
            

            % 3. モーメントアーム (Rcg, Rpg) と 基準長 (c)
            Rcg_x = params.Rcg_x;
            Rcg_z = params.Rcg_z;
            Rpg_x = params.Rpg_x;
            Rpg_z = params.Rpg_z;
            c_ref = params.c;

            % 4. モーメント計算 
            r_cg = [Rcg_x;0;Rcg_z];
            r_pg = [Rpg_x;0;Rpg_z];
            
            r_cg_cross = obj.vec_to_skew(r_cg);
            r_pg_cross = obj.vec_to_skew(r_pg);
            turn_alpha_rigging = obj.get_T_AC(alpha + mu_rigging+pi);

            % M_canopy_moment (キャノピー自体のモーメント)
            Cm_canopy_moment = C_mc +C_m_alpha *alpha; 
            
            % M_canopy_forces (キャノピー空力がCG周りに作るモーメント)
            % 機体軸 (Body Axis) での力 F_x, F_z
            %{
            F_x_c_coeff = -C_Dc * cos(alpha) + C_Lc * sin(alpha);
            F_z_c_coeff = -C_Dc * sin(alpha) - C_Lc * cos(alpha);
            % Cm = (R_x * F_z - R_z * F_x) / c_ref
            Cm_canopy_forces = (Rcg_x * F_z_c_coeff - Rcg_z * F_x_c_coeff) / c_ref;
            
            % M_payload_forces (ペイロード抗力がCG周りに作るモーメント)
            F_x_p_coeff = -C_Ds* cos(alpha); % (ペイロード抗力 Ds は -X 機体軸方向 と仮定)
            F_z_p_coeff = 0;
            Cm_payload_forces = (Rpg_x * F_z_p_coeff - Rpg_z * F_x_p_coeff) / c_ref;
            %}
            Cm_canopy_forces = r_cg_cross*turn_alpha_rigging*[C_Dc;0;C_Lc]/c_ref;
            Cm_payload_forces = r_pg_cross*turn_alpha_rigging*[C_Ds;0;0]/c_ref;
            % 5. 合計
            Cm_total = Cm_canopy_moment + Cm_canopy_forces(2,1) + Cm_payload_forces(2,1);
        end


        
        %%
        function S = vec_to_skew(~,v)
        % VEC_TO_SKEW 3次元ベクトルから対応する歪対称行列を作成します。
        %
        %   S = VEC_TO_SKEW(v)
        %
        %   入力:
        %     v - 3つの要素を持つベクトル (例: [v1; v2; v3])
        %
        %   出力:
        %     S - 対応する3x3の歪対称行列
        %
        %   歪対称行列 S は、任意のベクトル r に対して v x r = S * r となるように定義されます。
        
            % 入力ベクトルのサイズをチェック
            if numel(v) ~= 3
                error('入力ベクトルは3つの要素を持つ必要があります。');
            end
        
            % 歪対称行列の構成
            % S = [ 0, -v(3),  v(2);
            %      v(3),   0, -v(1);
            %     -v(2),  v(1),   0 ]
        
            S = [ 0,    -v(3),  v(2);
                 v(3),  0,    -v(1);
                -v(2),  v(1),   0 ];
        end
        %% 迎え角 alpha による座標変換 (風軸 -> 揚力・抗力軸)
        function T_AC = get_T_AC(~,alpha)
            
            c_a = cos(alpha); s_a = sin(alpha);
            T_AC = [c_a, 0, -s_a; 
                    0,   1,    0; 
                    s_a, 0,  c_a];
        end
    end
end