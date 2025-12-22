classdef ParafoilLinearizerAnalytic
    % ParafoilLinearizerAnalytic
    % 物理パラメータから安定微係数を解析的に計算し、線形モデルを構築するクラス
    
    properties
        model % ParafoilDynamics のインスタンス
    end
    
    methods
        function obj = ParafoilLinearizerAnalytic(parafoil_dynamics_instance)
            obj.model = parafoil_dynamics_instance;
        end
        
        function [A, B, E] = get_linear_model(obj, t, y_trim, u_trim_vec, wind_trim_vec)
            % 解析解を用いて線形化行列 A, B, E を導出する
            % 注意: このメソッドは簡単のため、トリム状態での微係数を計算します。
            % 厳密な非線形性の考慮には、y_trim から alpha0, V0 等を逆算する必要があります。
            
            p = obj.model.params;
            prop = obj.model.prop;
            
            % --- 1. トリム状態量の展開 ---
            u0 = y_trim(1); w0 = y_trim(3); q0 = y_trim(5); theta0 = y_trim(8);
            V0 = sqrt(u0^2 + w0^2);
            rho = obj.model.atmo_model.get_density(-y_trim(12)/1000);
            Q0 = 0.5 * rho * V0^2;
            S = p.S_c; b = prop.b; c = prop.c;
            
            % 重心から空力中心までの距離 (機体座標系)
            % ※ ParafoilDynamics の実装に合わせてベクトルを合成
            r_ac = prop.r_total_cm_to_canopy_origin_B + prop.r_canopy_origin_to_ac_B; 
            dx = r_ac(1); dz = r_ac(3); % 重心から見たAC位置 (通常 dz < 0)

            % --- 2. 基本的な空力係数の準備 (簡易版) ---
            % 本来は alpha0, Gamma0 を考慮した回転が必要ですが、
            % ここでは主要項である C_L_alpha, C_D_0 等を抽出します。
            
            % 質量・慣性特性
            m = prop.TotalMass;
            I = obj.model.inertia_tensor;
            Ixx = I(1,1); Iyy = I(2,2); Izz = I(3,3); Ixz = I(1,3);
            Gamma_in = Ixx*Izz - Ixz^2; % 慣性行列式
            
            % --- 3. 次元付き微係数の計算 (Dimensional Derivatives) ---
            
            % 縦系 (X, Z, M)
            % Xu: 速度減衰 (抗力 + ペイロード抗力)
            Xu = -rho * V0 * (S * p.C_D_0 + p.S_p * p.C_D_s);
            % Xw: 誘導抗力・揚力傾斜の前方成分 (簡易化: CL0 - CD_alpha)
            Xw = (Q0 * S / V0) * (p.C_L_0 - 0); % CD_alpha は通常小さいとして省略
            
            % Zu: 揚力 (速度比例)
            Zu = -rho * V0 * S * p.C_L_0;
            % Zw: 揚力傾斜 (ダンピング)
            Zw = -(Q0 * S / V0) * (p.C_L_alpha + p.C_D_0);
            
            % Mu, Mw, Mq (モーメントアーム dx, dz を考慮)
            % Mu = M_aero_u + Fx_u * dz + Fz_u * dx
            Mu = (rho * V0 * S * c / 2) * p.C_m_0 + (Xu * (-dz) + Zu * dx); % dzは負なので -dz は正の腕
            
            % Mw = M_aero_w + Fx_w * dz + Fz_w * dx
            Mw = (Q0 * S * c / (2*V0)) * p.C_m_alpha + (Xw * (-dz) + Zw * dx);
            
            % Mq: ピッチダンピング
            Mq = (Q0 * S * c^2 / (4*V0)) * p.C_m_q;

            % 横・方向系 (Y, L, N)
            % Yv: 横滑り力
            Yv = (Q0 * S / V0) * p.C_Y_beta;
            
            % Lv: ロールモーメント (空力 + 横力×高さ) -> 復元モーメント
            % Cl_beta がパラメータにない場合は Yv*(-dz) が支配的
            Lv = (Q0 * S * b / V0) * 0 + Yv * (-dz); 
            
            % Nv: ヨーモーメント (空力 + 横力×前後位置)
            Nv = (Q0 * S * b / V0) * 0 + Yv * dx;
            
            % Lp, Np, Lr, Nr (ダンピング項)
            Lp = (Q0 * S * b^2 / (2*V0)) * p.C_l_p;
            Lr = (Q0 * S * b^2 / (2*V0)) * 0; % C_l_r がない場合
            Np = (Q0 * S * b^2 / (2*V0)) * 0; % C_n_p がない場合
            Nr = (Q0 * S * b^2 / (2*V0)) * p.C_n_r;

            % --- 4. 行列要素の構築 (正規化 & 慣性連成解消) ---
            
            % 縦系 (Longitudinal)
            A_lon = zeros(4,4);
            A_lon(1,1) = Xu/m;  A_lon(1,2) = Xw/m;  A_lon(1,3) = -w0;          A_lon(1,4) = -9.81*cos(theta0);
            A_lon(2,1) = Zu/m;  A_lon(2,2) = Zw/m;  A_lon(2,3) = u0;           A_lon(2,4) = -9.81*sin(theta0);
            A_lon(3,1) = Mu/Iyy; A_lon(3,2) = Mw/Iyy; A_lon(3,3) = Mq/Iyy;
            A_lon(4,3) = 1;

            % 横系 (Lateral) - Primed Coefficients
            % L' = (Izz*L + Ixz*N)/Gamma_in
            % N' = (Ixz*L + Ixx*N)/Gamma_in
            L_prime = @(L, N) (Izz*L + Ixz*N) / Gamma_in;
            N_prime = @(L, N) (Ixz*L + Ixx*N) / Gamma_in;
            
            A_lat = zeros(4,4);
            A_lat(1,1) = Yv/m;        A_lat(1,2) = w0;          A_lat(1,3) = -u0;         A_lat(1,4) = 9.81*cos(theta0);
            A_lat(2,1) = L_prime(Lv, Nv); A_lat(2,2) = L_prime(Lp, Np); A_lat(2,3) = L_prime(Lr, Nr);
            A_lat(3,1) = N_prime(Lv, Nv); A_lat(3,2) = N_prime(Lp, Np); A_lat(3,3) = N_prime(Lr, Nr);
            A_lat(4,2) = 1; A_lat(4,3) = tan(theta0);

            % 全体行列 A (12x12) へのマッピング
            A = zeros(12,12);
            % 縦系: u, w, q, theta -> indices 1, 3, 5, 8
            idx_lon = [1, 3, 5, 8];
            A(idx_lon, idx_lon) = A_lon;
            % 横系: v, p, r, phi -> indices 2, 4, 6, 7
            idx_lat = [2, 4, 6, 7];
            A(idx_lat, idx_lat) = A_lat;
            % 位置 (x, y, z) の運動学 (簡易)
            A(10:12, 1:3) = eye(3); % dot_pos = vel (近似)

            % --- 5. 入力行列 B の構築 ---
            % 制御入力 u_vec = [delta_a; delta_s; Gamma]
            B = zeros(12, 3);
            
            % 非対称ブレーキ (delta_a) -> 横系モーメント
            L_da = Q0 * S * b * (p.C_l_delta_a / p.d);
            N_da = Q0 * S * b * (p.C_n_delta_a / p.d);
            
            B(4, 1) = L_prime(L_da, N_da); % p_dot
            B(6, 1) = N_prime(L_da, N_da); % r_dot
            
            % 対称ブレーキ (delta_s) -> 縦系 力・モーメント
            X_ds = -Q0 * S * p.C_D_delta_s;
            Z_ds = -Q0 * S * p.C_L_delta_s; % 通常0
            M_ds = Q0 * S * c * 0 + (X_ds*(-dz) + Z_ds*dx); % C_m_delta_s があれば追加
            
            B(1, 2) = X_ds / m;
            B(3, 2) = Z_ds / m;
            B(5, 2) = M_ds / Iyy;
            
            % 入射角 (Gamma) -> 力ベクトルの回転効果 (簡易化)
            % Fx_Gamma = Fz0, Fz_Gamma = -Fx0
            Fx0 = -rho*V0*S*p.C_D_0; % トリム抗力概算
            Fz0 = -rho*V0*S*p.C_L_0; % トリム揚力概算
            
            X_gam = Fz0; 
            Z_gam = -Fx0;
            M_gam = (X_gam*(-dz) + Z_gam*dx);
            
            B(1, 3) = X_gam / m;
            B(3, 3) = Z_gam / m;
            B(5, 3) = M_gam / Iyy;

            % --- 6. 外乱行列 E の構築 (風速入力) ---
            % w_vec = [uw; vw; ww] (機体座標系)
            E = zeros(12, 3);
            
            % 空力項の符号反転 (慣性項は含めない)
            % 縦系 (u, w) -> uw, ww
            E(1, 1) = -Xu/m;  E(1, 3) = -Xw/m;
            E(3, 1) = -Zu/m;  E(3, 3) = -Zw/m;
            E(5, 1) = -Mu/Iyy; E(5, 3) = -Mw/Iyy;
            
            % 横系 (v) -> vw
            E(2, 2) = -Yv/m;
            E(4, 2) = -L_prime(Lv, Nv);
            E(6, 2) = -N_prime(Lv, Nv);
            
        end
    end
end