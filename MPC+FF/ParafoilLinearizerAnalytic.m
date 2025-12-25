classdef ParafoilLinearizerAnalytic
    % ParafoilLinearizerAnalytic
    % 川口康太氏の修士論文モデルに基づき、任意の運動状態における
    % 6自由度パラフォイル運動方程式の線形化行列(A, B)を解析的に導出するクラス。
    
    properties
        model % ParafoilDynamics のインスタンス (パラメータ保持用)
        k_softmin = 20; % Softminの平滑化係数
    end
    
    methods
        function obj = ParafoilLinearizerAnalytic(parafoil_dynamics_instance, k_softmin)
            obj.model = parafoil_dynamics_instance;
            if nargin > 1
                obj.k_softmin = k_softmin;
            end
        end
        
        function [A, B] = get_linear_model(obj, t, y, u_vec, ~)
            % get_linear_model
            % 指定された状態 y と入力 u_vec における線形化行列 A, B を返す
            %
            % Inputs:
            %   y: 状態ベクトル [u; v; w; p; q; r; phi; theta; psi; x; y; z] (12x1)
            %   u_vec: 入力ベクトル [delta_R; delta_L] (2x1)
            % Outputs:
            %   A: システム行列 (12x12)
            %   B: 入力行列 (12x2)
            
            p = obj.model.params;   % 空力係数などが含まれる構造体
            prop = obj.model.prop;  % 質量、慣性、幾何形状
            
            % --- 1. 状態量の展開 ---
            u = y(1); v = y(2); w = y(3);
            pp = y(4); q = y(5); r = y(6); % pp = roll rate
            phi = y(7); theta = y(8); psi = y(9);
            z = y(12);
            
            delta_R = u_vec(1);
            delta_L = u_vec(2);
            
            % --- 2. 基礎物理量の計算 ---
            V_sq = u^2 + v^2 + w^2;
            V = sqrt(V_sq);
            if V < 0.1, V = 0.1; end % ゼロ除算防止
            
            rho = obj.model.atmo_model.get_density(-z/1000); % 高度は -z
            Q = 0.5 * rho * V_sq;
            S = p.S_c; % 基準面積
            b = prop.b;
            c = prop.c;
            
            alpha = atan2(w, u);
            beta = asin(v / V);
            
            % Softmin 計算
            k = obj.k_softmin;
            exp_nR = exp(-k * delta_R);
            exp_nL = exp(-k * delta_L);
            sum_exp = exp_nR + exp_nL;
            
            delta_s = -(1/k) * log(sum_exp);
            delta_a = delta_R - delta_L;
            
            % Softmin 勾配 (d_ds_d_dR, d_ds_d_dL)
            sigma_R = exp_nR / sum_exp;
            sigma_L = exp_nL / sum_exp;
            
            % 質量・慣性
            m = prop.TotalMass;
            I = obj.model.inertia_tensor;
            I_inv = inv(I);
            
            % --- 3. 偏微分係数の準備 (Partial Derivatives) ---
            
            % (3.1) 動圧・迎角・横滑り角の速度微分
            % dQ/du = rho*u, dQ/dv = rho*v, dQ/dw = rho*w
            dQ_du = rho * u; dQ_dv = rho * v; dQ_dw = rho * w;
            
            % d_alpha/du, d_alpha/dw
            den_a = u^2 + w^2;
            da_du = -w / den_a;
            da_dw =  u / den_a;
            
            % d_beta/dv (近似: 1/V)
            db_dv = 1 / V; 
            
            % (3.2) 空力係数の計算と微分
            % ※ パラメータ構造体 p の定義に依存しますが、ここでは一般的なモデルを想定
            
            % 基本係数 (線形近似: C = C0 + Ca*alpha + ...)
            CL = p.C_L_0 + p.C_L_alpha * alpha;
            CD = p.C_D_0 + p.C_D_alpha * alpha + p.C_D_delta_s * delta_s; % 対称ブレーキ考慮
            CY = p.C_Y_beta * beta;
            
            % 機体軸への変換 (簡易化: Gamma=0, beta小)
            % F_Ax = Q * S * (-CD*cos(a) + CL*sin(a))  (正確には T_AC 変換)
            % F_Az = Q * S * (-CD*sin(a) - CL*cos(a))
            
            sa = sin(alpha); ca = cos(alpha);
            CX_body = -CD * ca + CL * sa;
            CZ_body = -CD * sa - CL * ca;
            CY_body = CY;
            
            % 空力係数のalpha微分
            dCL_da = p.C_L_alpha;
            dCD_da = p.C_D_alpha; % 通常は小さいが考慮
            
            dCX_da = -(dCD_da * ca - CD * sa) + (dCL_da * sa + CL * ca);
            dCZ_da = -(dCD_da * sa + CD * ca) - (dCL_da * ca - CL * sa);
            
            % (3.3) 力の速度微分 (dX/du, dX/dw 等)
            % Fx = Q * S * CX_body
            % dFx/du = (dQ/du)*S*CX + Q*S*(dCX/da * da/du)
            
            QS = Q * S;
            Fx_u = dQ_du * S * CX_body + QS * dCX_da * da_du;
            Fx_w = dQ_dw * S * CX_body + QS * dCX_da * da_dw;
            % Fx_v = dQ_dv * S * CX_body; (横滑りによる動圧変化のみ)
            Fx_v = dQ_dv * S * CX_body;

            Fz_u = dQ_du * S * CZ_body + QS * dCZ_da * da_du;
            Fz_w = dQ_dw * S * CZ_body + QS * dCZ_da * da_dw;
            Fz_v = dQ_dv * S * CZ_body;
            
            % Fy (横力)
            % Fy = Q * S * CY_body
            % dFy/dv = dQ_dv * S * CY + QS * p.C_Y_beta * db_dv
            Fy_v = dQ_dv * S * CY_body + QS * p.C_Y_beta * db_dv;
            Fy_u = dQ_du * S * CY_body;
            Fy_w = dQ_dw * S * CY_body;
            
            % (3.4) モーメントの微分
            % 簡単のため、CM_alpha, Cl_beta, Cn_beta 等の主要項のみ展開
            
            % Cl (Roll)
            %dCl_db = p.C_l_beta; 
            dCl_db = 0;
            %Mx_v = (dQ_dv * S * b * (p.C_l_0 + dCl_db*beta) + QS * b * dCl_db * db_dv); 
            %c_l_0=0なので上の式は不適

            Mx_v = (dQ_dv * S * b * (dCl_db*beta) + QS * b * dCl_db * db_dv);
            % Cm (Pitch)
            Cm = p.C_m_0 + p.C_m_alpha * alpha;
            dCm_da = p.C_m_alpha;
            My_u = (dQ_du * S * c * Cm + QS * c * dCm_da * da_du);
            My_w = (dQ_dw * S * c * Cm + QS * c * dCm_da * da_dw);
            
            % Cn (Yaw)
            %dCn_db = p.C_n_beta;
            dCn_db =0;

            %Mz_v = (dQ_dv * S * b * (p.C_n_0 + dCn_db*beta) + QS * b * dCn_db * db_dv);
            %c_n_0=0なので上の式は不適
            Mz_v = (dQ_dv * S * b * ( dCn_db*beta) + QS * b * dCn_db * db_dv);
            % --- 4. 行列 A の構築 (12x12) ---
            A = zeros(12, 12);
            
            % [Block 1,1] dv/dv (並進速度)
            % (1/m)*dF/dv - Omega x
            A(1,1) = Fx_u/m;      A(1,2) = Fx_v/m + r;  A(1,3) = Fx_w/m - q;
            A(2,1) = Fy_u/m - r;  A(2,2) = Fy_v/m;      A(2,3) = Fy_w/m + pp;
            A(3,1) = Fz_u/m + q;  A(3,2) = Fz_v/m - pp; A(3,3) = Fz_w/m;
            
            % [Block 1,2] dv/domega (角速度依存: コリオリ + 空力ダンピングがあれば)
            A(1,5) = -w; A(1,6) = v;
            A(2,4) = w;  A(2,6) = -u;
            A(3,4) = -v; A(3,5) = u;
            
            % [Block 1,3] dv/dTheta (重力項)
            g = 9.81;
            st = sin(theta); ct = cos(theta);
            sp = sin(phi);   cp = cos(phi);
            
            A(1,8) = -g * ct;
            A(2,7) =  g * ct * cp; A(2,8) = -g * st * sp;
            A(3,7) = -g * ct * sp; A(3,8) = -g * st * cp;
            
            % [Block 2,1] domega/dv (静安定)
            % I_inv * dM/dv
            dM_dv_mat = zeros(3,3);
            dM_dv_mat(1,2) = Mx_v; % Roll due to sideslip
            dM_dv_mat(2,1) = My_u; dM_dv_mat(2,3) = My_w; % Pitch terms
            dM_dv_mat(3,2) = Mz_v; % Yaw due to sideslip
            
            A(4:6, 1:3) = I_inv * dM_dv_mat;
            
            % [Block 2,2] domega/domega (ダンピング + ジャイロ)
            % 空力ダンピング行列 D
            D = zeros(3,3);
            coeff_p = 0.25 * rho * V * S * b^2; % 簡易化 (本来は u 依存等あるが V で代表)
            coeff_q = 0.25 * rho * V * S * c^2;
            
            D(1,1) = coeff_p * p.C_l_p; %D(1,3) = coeff_p * p.C_l_r;
            D(1,3) = 0;%C_l_r=0
            D(2,2) = coeff_q * p.C_m_q;
            %D(3,1) = coeff_p * p.C_n_p; 
            D(3,1) = 0;
            D(3,3) = coeff_p * p.C_n_r;
            
            % ジャイロ項 G = d/dw (w x Iw)
            % I = [Ixx 0 -Ixz; 0 Iyy 0; -Ixz 0 Izz] を仮定
            Ixx = I(1,1); Iyy = I(2,2); Izz = I(3,3); Ixz = I(1,3);
            G = zeros(3,3);
            G(1,2) = -Iyy*q + Izz*r + Ixz*pp; G(1,3) = (Iyy-Izz)*q;
            G(2,1) = (Ixx-Izz)*r - 2*Ixz*pp;  G(2,3) = (Ixx-Izz)*pp + 2*Ixz*r;
            G(3,1) = (Iyy-Ixx)*q;             G(3,2) = (Ixx-Iyy)*pp - Ixz*r;
            
            A(4:6, 4:6) = I_inv * (D - G);
            
            % [Block 3,2] dTheta/domega (キネマティクス)
            tt = tan(theta); sect = sec(theta);
            R_kin = [1, sp*tt, cp*tt;
                     0, cp,    -sp;
                     0, sp*sect, cp*sect];
            A(7:9, 4:6) = R_kin;
            
            % [Block 3,3] dTheta/dTheta (角速度がある場合のみ非ゼロ)
            if norm([pp,q,r]) > 1e-6
               % 複雑になるため、ここではq,r項の一部のみ実装例として示す
               % (完全な実装には前述のLaTeX式の通り記述が必要)
               dq_phi = q*cp - r*sp;
               dq_psi = q*sp + r*cp;
               
               A(7,7) = dq_phi * tt;      A(7,8) = dq_psi * (sect^2);
               A(8,7) = -dq_psi;
               A(9,7) = dq_phi * sect;    A(9,8) = dq_psi * sect * tt;
            end
            
            % [Block 4,1] dPos/dv (座標変換)
            % T_IB'
            T_IB_T = [ct*cos(psi), sp*st*cos(psi)-cp*sin(psi), cp*st*cos(psi)+sp*sin(psi);
                      ct*sin(psi), sp*st*sin(psi)+cp*cos(psi), cp*st*sin(psi)-sp*cos(psi);
                      -st,         sp*ct,                      cp*ct];
            A(10:12, 1:3) = T_IB_T;
            
            % --- 5. 行列 B の構築 (12x2) ---
            B = zeros(12, 2);
            
            % 操舵微係数の準備 (dC/d_delta_a, dC/d_delta_s)
            % 1. 非対称 (delta_a) -> Roll, Yaw, Y-force
            %    F = Q*S * ...
            dCl_da_coef = p.C_l_delta_a / p.d; % 無次元化注意
            dCn_da_coef = p.C_n_delta_a / p.d;
            dCy_da_coef = 0; % 必要なら追加
            
            % 2. 対称 (delta_s) -> Drag, Lift, Pitch
            dCD_ds_coef = p.C_D_delta_s;
            dCL_ds_coef = p.C_L_delta_s; 
            dCm_ds_coef = 0; % 必要なら p.C_m_delta_s
            
            % 機体軸への投影 (対称ブレーキ)
            dCX_ds = -dCD_ds_coef * ca + dCL_ds_coef * sa;
            dCZ_ds = -dCD_ds_coef * sa - dCL_ds_coef * ca;
            
            % 物理入力への変換行列 T_soft (2x2)
            % [d_da/d_dR, d_da/d_dL;
            %  d_ds/d_dR, d_ds/d_dL]
            T_soft = [1, -1; 
                      sigma_R, sigma_L];
            
            % 並進力入力 B_force_raw (3x2) -> [Fx_da, Fx_ds; Fy_da...; Fz_da...]
            B_force_raw = zeros(3, 2);
            B_force_raw(1, 2) = QS * dCX_ds;
            B_force_raw(3, 2) = QS * dCZ_ds;
            % B_force_raw(2, 1) = QS * dCy_da_coef; % もしあれば
            
            % モーメント入力 B_mom_raw (3x2)
            B_mom_raw = zeros(3, 2);
            B_mom_raw(1, 1) = QS * b * dCl_da_coef;
            B_mom_raw(3, 1) = QS * b * dCn_da_coef;
            B_mom_raw(2, 2) = QS * c * dCm_ds_coef;
            
            % 状態空間行列へマッピング
            B_force = (B_force_raw * T_soft) / m;
            B_mom   = I_inv * (B_mom_raw * T_soft);
            
            B(1:3, :) = B_force;
            B(4:6, :) = B_mom;
            
        end
    end
end