classdef ParafoilControlMapper_linear < handle
    % PARAFOILCONTROLMAPPER
    % 軌道計画の結果(V, theta, phi)からフィードフォワード操作量を計算し、
    % さらに現在の状態誤差に基づいて線形化モデルによる補正(delta_delta_a)を加えるクラス
    %
    % ※ 変更点: Ixz=0, Izzでの除算なし (モーメント釣り合いベース)
    
    properties
        Params       % 全パラメータ構造体
        Yaw_Factor   % ベースFF用係数
    end
    
    methods
        function obj = ParafoilControlMapper_linear(params)
            % コンストラクタ
            obj.Params = params;
            
            % --- パラメータ展開 ---
            if isfield(params, 'prop')
                b = params.prop.b;
            else
                b = params.b;
            end
            
            if isfield(params, 'params')
                d = params.params.d;
                Cn_r = params.params.C_n_r;
                Cn_da = params.params.C_n_delta_a;
            else
                d = params.d;
                Cn_r = params.C_n_r;
                Cn_da = params.C_n_delta_a;
            end
            
            if abs(Cn_da) < 1e-9
                error('C_n_delta_a is too small.');
            end
            
            % ベースFF係数 (モーメント係数の比率なので慣性は不要)
            obj.Yaw_Factor = -d * b * Cn_r / (2 * Cn_da);
        end
        
        % --- 既存メソッド (変更なし) ---
        function [delta_R, delta_L, delta_a] = compute_feedforward_input(obj, state_vec, phi_cmd, delta_s_bias)
            if nargin < 4, delta_s_bias = 0; end
            u = state_vec(1); v = state_vec(2); w = state_vec(3);
            theta = state_vec(8);
            V_sq = u^2 + v^2 + w^2;
            g = 9.81;
            if V_sq < 1.0, delta_a = 0;
            else
                K_dyn = (g / V_sq) * cos(theta);
                delta_a = obj.Yaw_Factor * K_dyn * sin(phi_cmd);
            end
            [delta_R, delta_L] = obj.apply_mixing(delta_a, delta_s_bias);
        end
        
        function [delta_R, delta_L, delta_a] = compute_input_from_reference(obj, phi_ref, V_ref, theta_ref, delta_s_bias)
            if nargin < 5, delta_s_bias = 0; end
            g = 9.81; 
            if V_ref < 1.0, V_ref = 1.0; end
            K_dyn = (g / (V_ref^2)) * cos(theta_ref);
            delta_a = obj.Yaw_Factor * K_dyn * sin(phi_ref);
            [delta_R, delta_L] = obj.apply_mixing(delta_a, delta_s_bias);
        end
        
        % ★★★ モーメント釣り合い(Izz不要)版 補正計算メソッド ★★★
        function [delta_R, delta_L, delta_a_total, delta_delta_a] = compute_corrected_input(obj, current_state, ref_state, delta_s_bias)
            if nargin < 4, delta_s_bias = 0; end
            
            % 1. 参照状態の展開
            if isstruct(ref_state)
                V_ref = ref_state.V;
                phi_ref = ref_state.phi;
                theta_ref = ref_state.theta;
                
                % ジャイロ項計算用の角速度
                if isfield(ref_state, 'q') && isfield(ref_state, 'r')
                    q_ref = ref_state.q;
                    r_ref = ref_state.r;
                else
                    g = 9.81;
                    psi_dot = (g / V_ref) * tan(phi_ref);
                    q_ref = psi_dot * sin(phi_ref) * cos(theta_ref);
                    r_ref = psi_dot * cos(phi_ref) * cos(theta_ref);
                end
                
                % 参照操作量 (Trim)
                if isfield(ref_state, 'delta_a')
                    da_ref = ref_state.delta_a;
                else
                    [~, ~, da_ref] = obj.compute_input_from_reference(phi_ref, V_ref, theta_ref, 0);
                end
            else
                % 簡易ベクトルの場合
                V_ref = ref_state(1); theta_ref = ref_state(2); phi_ref = ref_state(3);
                g = 9.81;
                psi_dot = (g / V_ref) * tan(phi_ref);
                q_ref = psi_dot * sin(phi_ref) * cos(theta_ref);
                r_ref = psi_dot * cos(phi_ref) * cos(theta_ref);
                [~, ~, da_ref] = obj.compute_input_from_reference(phi_ref, V_ref, theta_ref, 0);
            end
            
            % 2. 状態偏差 (Delta)
            u_curr = current_state(1); w_curr = current_state(3);
            phi_curr = current_state(7);
            
            % 参照速度成分 (alpha=0, beta=0 仮定)
            u_ref = V_ref; 
            w_ref = 0;
            
            du = u_curr - u_ref;
            dw = w_curr - w_ref;
            dphi = phi_curr - phi_ref;
            
            % 3. モーメント感度係数の計算 (Izで割らず、Nm次元で計算)
            p = obj.Params;
            S = p.S_c; b = p.prop.b;
            rho = 1.225; % 必要に応じて高度補正
            
            % 定数定義 (PDFの k1, k2, k3 相当)
            % N = (1/4 rho S b^2 Cn_r) * V * r  + (1/2 rho S b Cn_da) * V^2 * da
            
            K_damp = 0.25 * rho * S * b^2 * p.C_n_r;       % V*r の係数
            K_ctl  = 0.5  * rho * S * b   * p.C_n_delta_a; % V^2*da の係数
            
            % (A) 速度感度項 dN/du, dN/dw
            % dN/dV = K_damp * r + 2 * K_ctl * V * da
            dN_dV = K_damp * r_ref + 2 * K_ctl * V_ref * da_ref;
            
            % u, w 方向への分解 (Izz不要)
            N_u = dN_dV * (u_ref / V_ref);
            N_w = dN_dV * (w_ref / V_ref);
            
            % (B) バンク角・ジャイロ補正項 (N_r, N_q, N_phi)
            % PDF式: delta_r = -q * dphi, delta_q = r * dphi
            % 項: (N_r * (-q) + N_q * r + N_phi) * dphi
            
            % N_r (Yaw Damping term): dN/dr
            % N_r = K_damp * V
            N_r = K_damp * V_ref;
            
            % N_q (Gyro term): dN/dq
            % ジャイロモーメント N_gyro approx (Ixx - Iyy) * p * q
            % dN_gyro/dq = (Ixx - Iyy) * p_ref
            % ※ Ixz=0 としたので Ixx, Iyy の差分項のみ考慮
            Ixx = p.I_xx; Iyy = p.I_yy;
            p_ref = -(9.81 / V_ref * tan(phi_ref)) * sin(theta_ref); % Kinematics
            
            N_q = (Ixx - Iyy) * p_ref;
            
            % N_phi (Gravity term)
            % Ixz=0 の場合、重力はYawモーメントに直接影響しないためゼロ
            N_phi = 0;
            
            % 感度総和
            S_phi = N_r * (-q_ref) + N_q * (r_ref) + N_phi;
            
            % (C) 制御効力係数 N_da
            % dN/dda = K_ctl * V^2
            N_da = K_ctl * V_ref^2;
            
            % 4. 補正量 delta_delta_a の計算
            % dda = -1/N_da * [ (N_u*du + N_w*dw) + S_phi*dphi ]
            % (モーメントの釣り合い式 0 = dN + N_da * dda より)
            
            if abs(N_da) < 1e-9
                delta_delta_a = 0;
            else
                term_vel = N_u * du + N_w * dw;
                term_phi = S_phi * dphi;
                delta_delta_a = - (1 / N_da) * (term_vel + term_phi);
            end
            
            % 5. 合成
            delta_a_total = da_ref + delta_delta_a;
            
            [delta_R, delta_L] = obj.apply_mixing(delta_a_total, delta_s_bias);
        end
    end
    
    methods (Access = private)
        function [delta_R, delta_L] = apply_mixing(~, delta_a, delta_s_bias)
            delta_a = max(-1.0, min(1.0, delta_a));
            if delta_a > 0
                delta_R = delta_a + delta_s_bias;
                delta_L = delta_s_bias;
            else
                delta_R = delta_s_bias;
                delta_L = abs(delta_a) + delta_s_bias;
            end
            delta_R = max(0, min(1, delta_R));
            delta_L = max(0, min(1, delta_L));
        end
    end
end