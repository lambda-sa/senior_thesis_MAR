classdef ParafoilDynamics
    properties
        prop
        inertia_tensor
        params
        atmo_model

         % --- 追加: 線形化用スイッチ ---
        use_soft_min = false; % デフォルトは false (厳密な min を使用)
        soft_min_k = 50;      % Linearizerから上書きされる値

    end
    
    methods
        function obj = ParafoilDynamics(params)
            obj.prop = params.prop;
            obj.inertia_tensor = params.I_total_body;
            obj.params = params;
            obj.atmo_model = AtmoTempPressRho();

           
        end
        
        function dy  = get_derivatives(obj, t, y, control_inputs)
            % 状態変数の抽出
            u=y(1); v=y(2); w=y(3); p=y(4); q=y(5); r=y(6);
            phi=y(7); theta=y(8); psi=y(9); x=y(10); y_pos=y(11); 
            z_pos = y(12); % 慣性系 Z 座標 (下向き正)
            h_altitude_m = -z_pos; % 高度 (上向き正, [m])
            
            % 外部引数から入力を受け取る
            delta_R = control_inputs.delta_R; % 右ブレーキ
            delta_L = control_inputs.delta_L; % 左ブレーキ
            GAMMA   = control_inputs.GAMMA;   % ガンマ角
            wind_I  = control_inputs.wind_I; % 風速を構造体から取得

            % --- 修正: モードによる関数の切り替え ---
            if obj.use_soft_min
                % 線形化用: SoftMinを使用 (滑らか)
                delta_s = obj.calculate_soft_min(delta_R, delta_L, obj.soft_min_k);
            else
                % シミュレーション用: 厳密な min を使用 (カドあり)
                delta_s = min(delta_R, delta_L);
            end
            
            delta_a = delta_R - delta_L;
            
            %rho = obj.params.rho_sea_level;
            %rho_calc = AtmoTempPressRho();
            rho = obj.atmo_model.get_density(h_altitude_m/1000);
            
            % ベクトル定義
            V_cg_B = [u; v; w];
            omega_B = [p; q; r];
            V_wind_I = wind_I;

            % 座標変換行列の取得
            T_IB = RotationMatrices.get_inertial_to_body_matrix(phi, theta, psi);
            T_BC = RotationMatrices.get_T_BC(GAMMA);
            V_wind_B = T_IB * V_wind_I;
            
            % 各部の速度計算
            V_canopy_C = obj.calculate_canopy_velocity_C(V_cg_B, omega_B, V_wind_B, T_BC);
            u_c = V_canopy_C(1); v_c = V_canopy_C(2); w_c = V_canopy_C(3);
            V_payload_B = obj.calculate_payload_velocity_B(V_cg_B, omega_B, V_wind_B);
            
            % 空力角度の計算
            V_A = norm(V_canopy_C);
            if V_A < 1e-6, alpha=0; beta=0; 
            else
                alpha = atan2(w_c, u_c);
                
                % ★★★ 保護コードを追加 ★★★
                arg_for_asin = v_c / V_A;
                if arg_for_asin > 1.0
                    arg_for_asin = 1.0;
                elseif arg_for_asin < -1.0
                    arg_for_asin = -1.0;
                end
                beta = asin(arg_for_asin);
                % ★★★ ここまで ★★★
            end
            % 現在高度における重力加速度を取得
            g_current = obj.atmo_model.get_gravity(h_altitude_m);

            % 力とモーメントの計算
            F_A = obj.force_F_A(alpha, V_A, GAMMA, delta_a, delta_s, beta, rho);
            F_S = obj.force_F_S(V_payload_B, rho);
            F_W = obj.force_F_W(phi, theta,g_current);
            %F_total = F_A + F_S + F_W;
            F_total =  F_A + F_S+F_W;

            M_A = obj.moment_M_A(p, q, r, V_A, alpha, delta_a, rho, phi, GAMMA);
            
            %空力座標系の原点の、機体重心空の位置ベクトルを求める
            r_aero = obj.prop.r_total_cm_to_canopy_origin_B + T_BC* obj.prop.r_canopy_origin_to_ac_B;
            M_aero = obj.calculate_moment_from_force(r_aero, F_A);


            r_pay = obj.prop.r_total_cm_to_payload_cm_B;
            M_PayloadDrag = obj.calculate_moment_from_force(r_pay, F_S);

            M_total = M_A +M_aero+ M_PayloadDrag;
            %M_total = M_PayloadDrag;
            % 運動方程式
            total_mass = obj.prop.TotalMass;
            acce_vec = F_total / total_mass - cross(omega_B, V_cg_B);
            u_dot = acce_vec(1); v_dot = acce_vec(2); w_dot = acce_vec(3);
            
            ang_acce_vec = obj.inertia_tensor \ (M_total - cross(omega_B, obj.inertia_tensor * omega_B));
            p_dot = ang_acce_vec(1); q_dot = ang_acce_vec(2); r_dot = ang_acce_vec(3);
            
            attitude_dot_mat = RotationMatrices.get_attitude_dot_matrix(phi, theta);
            attitude_dot_vec = attitude_dot_mat * omega_B;
            phi_dot = attitude_dot_vec(1); theta_dot = attitude_dot_vec(2); psi_dot = attitude_dot_vec(3);
            
            inertial_vel = T_IB' * V_cg_B;
            
            x_dot = inertial_vel(1); y_dot = inertial_vel(2); z_dot = inertial_vel(3);
            

            dy = [u_dot; v_dot; w_dot; p_dot; q_dot; r_dot; phi_dot; theta_dot; psi_dot; x_dot; y_dot; z_dot];
        end

        % ★3. デバッグ用の新しいメソッドを追加★
        function [F, M] = get_forces_and_moments(obj, t, y_vec)
             % このメソッドは get_derivatives の計算を呼び出して、
             % 力とモーメントだけを返すラッパーです。
             control_inputs.delta_R = 0;
             control_inputs.delta_L = 0;
             control_inputs.GAMMA = 0;
             wind_I = [0;0;0];
             [~, F, M] = obj.get_derivatives(t, y_vec, control_inputs, wind_I);
        end

        % --- 誘導則クラス用のヘルパーメソッド ---
        function vel_I = get_inertial_velocity(~, y)
            % 状態ベクトルyから慣性座標系の速度ベクトル [Vx; Vy; Vz] を計算して返す
            u=y(1); v=y(2); w=y(3);
            phi=y(7); theta=y(8); psi=y(9);
            V_cg_B = [u; v; w];
            T_IB = RotationMatrices.get_inertial_to_body_matrix(phi, theta, psi);
            vel_I = T_IB' * V_cg_B;
        end

        function V_xy = get_horizontal_speed(obj, y)
            % 状態ベクトルyから慣性座標系の水平速度の大きさ V_xy を計算して返す
            vel_I = obj.get_inertial_velocity(y);
            V_xy = sqrt(vel_I(1)^2 + vel_I(2)^2);
        end

        function V_down = get_downward_speed(obj, y)
            % 状態ベクトルyから慣性座標系の降下速度 V_down (正が下向き) を計算して返す
            vel_I = obj.get_inertial_velocity(y);
            V_down = -vel_I(3);
        end

        function z = get_altitude(~, y)
            % 状態ベクトルyから高度z (正が上向き) を返す
            z = y(12);
        end

       function eul_rates = get_euler_rates(obj, y)
            % 状態ベクトルyから必要な値を取得
            p = y(4); q = y(5); r = y(6);
            phi = y(7); theta = y(8);
        
            omega_B = [p; q; r];
        
            % get_derivatives と同じ計算を実行
            attitude_dot_mat = RotationMatrices.get_attitude_dot_matrix(phi, theta);
            eul_rates = attitude_dot_mat * omega_B; % [phi_dot; theta_dot; psi_dot]
        end

    end
    
    methods (Access = private)
        function V_canopy_C = calculate_canopy_velocity_C(obj, V_cg_B, omega_B, V_wind_B, T_BC)
            % (全体CG -> キャノピー原点) + (キャノピー原点 -> AC)
            r_cm_to_canopy_ac_B = obj.prop.r_total_cm_to_canopy_origin_B + obj.prop.r_canopy_origin_to_ac_B;
            
            V_canopy_ac_B = V_cg_B + cross(omega_B, r_cm_to_canopy_ac_B);
            V_relative_ac_B = V_canopy_ac_B - V_wind_B;
            V_canopy_C = T_BC * V_relative_ac_B;
        end

        function V_payload_B = calculate_payload_velocity_B(obj, V_cg_B, omega_B, V_wind_B)
            r_cm_to_payload_cm_B = obj.prop.r_total_cm_to_payload_cm_B;
            V_payload_cm_B = V_cg_B + cross(omega_B, r_cm_to_payload_cm_B);
            V_payload_B = V_payload_cm_B - V_wind_B;
        end

        function F_W = force_F_W(obj, phi, theta,g_current)
            F_W = obj.prop.TotalMass *g_current * [-sin(theta); sin(phi)*cos(theta); cos(phi)*cos(theta)];
        end

        function f_A = force_F_A(obj, alpha, V_A, GAMMA, delta_a, delta_s, beta, rho)
            q_barS = 0.5 * rho * V_A^2 * obj.params.S_c;
            C_L = obj.params.C_L_0 + obj.params.C_L_alpha * alpha + obj.params.C_L_delta_s * delta_s + obj.params.C_L_delta_a * delta_a;
            C_D = obj.params.C_D_0 + obj.params.C_D_alpha * alpha^2 + obj.params.C_D_delta_a * delta_a + obj.params.C_D_delta_s * delta_s;
            C_Y = obj.params.C_Y_beta * beta;
            F_aerodynamic_axis = q_barS * [-C_D; C_Y; -C_L];
            T_AC = RotationMatrices.get_T_AC(alpha);
            T_BC = RotationMatrices.get_T_BC(GAMMA);
            f_A = T_BC' * T_AC * F_aerodynamic_axis;
        end
       
      function M_A = moment_M_A(obj, p, q, r, V_A, alpha, delta_a, rho, phi, GAMMA)
            q_barS = 0.5 * rho * V_A^2 * obj.params.S_c;
            b = obj.prop.b; 
            c = obj.prop.c;
            d = obj.params.d;
            
            % ★★★ 追加: 角速度ベクトルを機体軸(B)から安定軸(C)へ変換 ★★★
            omega_B = [p; q; r];
            T_BC = RotationMatrices.get_T_BC(GAMMA);
            omega_C = T_BC * omega_B;
            p_c = omega_C(1);
            q_c = omega_C(2);
            r_c = omega_C(3);
            % ★★★ 追加ここまで ★★★
            
            if V_A < 1e-6
                Roll = 0; Pitch = 0; Yaw = 0;
            else
                % ★★★ 変更: モーメント係数の計算に p_c, q_c, r_c を使用 ★★★
                % 1. ロールモーメント係数 (Roll)
                Roll = obj.params.C_l_phi * phi ...
                     + (obj.params.C_l_p * p_c * b) / (2 * V_A) ...
                     + (obj.params.C_l_delta_a * delta_a) / d;
                
                % 2. ピッチモーメント係数 (Pitch)
                Pitch = obj.params.C_m_0 ...
                      + obj.params.C_m_alpha * alpha ...
                      + (obj.params.C_m_q * q_c * c) / (2 * V_A);
                      
                % 3. ヨーモーメント係数 (Yaw)
                Yaw = (obj.params.C_n_r * r_c * b) / (2 * V_A) ...
                    + (obj.params.C_n_delta_a * delta_a) / d;
            end
            
            % 各モーメントをベクトルにまとめる
            moment_vec = q_barS * [b * Roll; c * Pitch; b * Yaw];
            
            % 最終的な座標変換 (C -> B) を適用
            M_A = T_BC' * moment_vec;
    end

        function f_S = force_F_S(obj, V_payload_B, rho)
            V_s = norm(V_payload_B);
            if V_s < 1e-6, f_S = [0;0;0];
            else, f_S = -0.5 * rho * V_s * obj.params.S_p * obj.params.C_D_s * V_payload_B;
            end
        end

        function M = calculate_moment_from_force(~, r, F)
            M = cross(r, F);
        end

        % SoftMin 計算メソッド
        function val = calculate_soft_min(~, x, y, k)
            % Log-Sum-Exp Trick
            val = - (1/k) * log(exp(-k*x) + exp(-k*y));
        end
        
    end
end
