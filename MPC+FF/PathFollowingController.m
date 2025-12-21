classdef PathFollowingController
    properties
        model      % ParafoilDynamics
        params     % パラメータ
        
        % ゲイン
        Kp_guidance = 0.5; 
        Kp_rate     = 2.0; 
        Kd_rate     = 4.0; 
    end
    
    methods
        function obj = PathFollowingController(model_instance)
            obj.model = model_instance;
            obj.params = model_instance.params;
        end
        
        function [ctrl_input, debug_info] = compute_control(obj, t, y, current_path_segment)
            % --- 1. 状態量の抽出 ---
            u=y(1); v=y(2); w=y(3); 
            p=y(4); q=y(5); r=y(6);
            psi=y(9); % 現在の方位角
            pos_x=y(10); pos_y=y(11); h=-y(12);
            V_A = sqrt(u^2 + v^2 + w^2);
            
            % --- 2. ガイダンス則 (経路タイプによる分岐) ---
            if strcmp(current_path_segment.type, 'line')
                % === 直線追従ロジック ===
                wp_start = current_path_segment.start;
                wp_end   = current_path_segment.end;
                
                path_vec = wp_end - wp_start;
                psi_path = atan2(path_vec(2), path_vec(1));
                
                % クロストラックエラー計算
                d_vec = [pos_x; pos_y] - wp_start;
                t_vec = path_vec / norm(path_vec);
                e = d_vec(1)*(-t_vec(2)) + d_vec(2)*t_vec(1); 
                
                % 直線用ガイダンス (Vector Field)
                chi_infty = deg2rad(60); 
                psi_cmd = psi_path - 2/pi * chi_infty * atan(obj.Kp_guidance * e);
                
            elseif strcmp(current_path_segment.type, 'circle')
                % === 円追従ロジック ===
                center = current_path_segment.center;
                radius = current_path_segment.radius;
                direction = current_path_segment.direction; % 1: CW, -1: CCW
                
                % 中心からの距離ベクトル
                d_vec = [pos_x; pos_y] - center;
                dist_to_center = norm(d_vec);
                
                % 半径誤差 (正: 円の外側, 負: 円の内側)
                e = dist_to_center - radius;
                
                % 接線方向の方位角 (現在位置における円の接線)
                % atan2(dy, dx) は中心から外に向かう角度。
                % ここに +90度(CCW) または -90度(CW) すると接線になる。
                psi_tangent = atan2(d_vec(2), d_vec(1)) + direction * (-pi/2);
                
                % 補正角 (外側にいるときは内側へ、内側にいるときは外側へ)
                % errorが正(外側)なら、CW時は右(-)、CCW時は左(+)に向けたい
                % つまり direction * (-atan(Kp * e))
                chi_infty = deg2rad(60); % 最大修正角
                correction = direction * (-2/pi * chi_infty * atan(obj.Kp_guidance * e));
                
                psi_cmd = psi_tangent + correction;
            else
                psi_cmd = psi; % 未定義なら直進維持
                e = 0;
            end
            
            % --- 共通処理: 角度差の計算 ---
            psi_err = atan2(sin(psi_cmd - psi), cos(psi_cmd - psi));
            
            % --- 3. 運動制御 (逆モデル制御) ---
            % ここは直線でも円でも全く同じ！
            r_cmd = obj.Kp_rate * psi_err;
            r_dot_req = obj.Kd_rate * (r_cmd - r);
            
            % --- 4. 逆モデル演算 ---
            I = obj.model.inertia_tensor;
            N_inertial_coupling = (I(1,1) - I(2,2)) * p * q; 
            N_req = I(3,3) * r_dot_req - N_inertial_coupling;
            
            rho = obj.model.atmo_model.get_density(h/1000);
            q_bar = 0.5 * rho * V_A^2;
            S = obj.params.S_c;
            b = obj.model.prop.b;
            
            if V_A > 0.1, beta = asin(v/V_A); else, beta=0; end
            r_hat = (r * b) / (2 * V_A);
            
            if isfield(obj.params, 'C_n_beta'), Cn_beta_val = obj.params.C_n_beta; else, Cn_beta_val = 0; end
            
            Cn_base = Cn_beta_val * beta + obj.params.C_n_r * r_hat;
            Cn_ctrl = obj.params.C_n_delta_a;
            
            denom = q_bar * S * b * Cn_ctrl;
            
            if abs(denom) > 1e-6
                term1 = N_req / (q_bar * S * b);
                delta_a_raw = (term1 - Cn_base) / Cn_ctrl;
            else
                delta_a_raw = 0;
            end
            
            % --- 5. 入力整形 ---
            delta_s_nominal = 0.1;
            delta_a_clamped = max(min(delta_a_raw, 1.0), -1.0);
            
            if delta_a_clamped > 0
                delta_R = delta_s_nominal + delta_a_clamped;
                delta_L = delta_s_nominal;
            else
                delta_R = delta_s_nominal;
                delta_L = delta_s_nominal - delta_a_clamped;
            end
            
            ctrl_input.delta_R = max(min(delta_R, 1), 0);
            ctrl_input.delta_L = max(min(delta_L, 1), 0);
            ctrl_input.GAMMA = 0;
            ctrl_input.wind_I = [0;0;0];
            
            debug_info.e = e;
            debug_info.psi_cmd = psi_cmd;
            debug_info.r_req = r_cmd;
            debug_info.da_raw = delta_a_raw;
        end
    end
end