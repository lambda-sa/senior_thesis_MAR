classdef RotationMatrices
    % 回転行列と座標変換に関する静的メソッドを提供するユーティリティクラス
    
    methods (Static)
        function T_IB = get_inertial_to_body_matrix(phi, theta, psi)
            % 慣性座標系(I)からボディ座標系(B)への座標変換行列
            c_phi = cos(phi);   s_phi = sin(phi);
            c_th  = cos(theta); s_th  = sin(theta);
            c_psi = cos(psi);   s_psi = sin(psi);
            
            L_BI = [c_psi*c_th, c_psi*s_th*s_phi - s_psi*c_phi, c_psi*s_th*c_phi + s_psi*s_phi;
                    s_psi*c_th, s_psi*s_th*s_phi + c_psi*c_phi, s_psi*s_th*c_phi - c_psi*s_phi;
                    -s_th,      c_th*s_phi,                     c_th*c_phi];
            T_IB = L_BI';
        end
        
        function T = get_attitude_dot_matrix(phi, theta)
            % ボディ角速度をオイラー角微分に変換する行列
            c_phi = cos(phi);   s_phi = sin(phi);
            t_th  = tan(theta); c_th  = cos(theta);
            if abs(c_th) < 1e-9, warning('Pitch angle is at a singularity.'); end
            T = [1, s_phi * t_th, c_phi * t_th;
                 0, c_phi,       -s_phi;
                 0, s_phi / c_th, c_phi / c_th];
        end

        function T_AC = get_T_AC(alpha)
            % 迎え角 alpha による座標変換 (風軸 -> 揚力・抗力軸)
            c_a = cos(alpha); s_a = sin(alpha);
            T_AC = [c_a, 0, -s_a; 
                    0,   1,    0; 
                    s_a, 0,  c_a];
        end
        
        function T_BC = get_T_BC(GAMMA)
            % GAMMA角による座標変換 (ボディ軸 -> 風軸)
            c_g = cos(GAMMA); s_g = sin(GAMMA);
            T_BC = [c_g, 0, -s_g; 
                    0,   1,    0; 
                    s_g, 0,  c_g];
        end
    end
end