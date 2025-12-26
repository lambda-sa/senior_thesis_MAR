classdef ParafoilProperties
    % パラフォイルの物理的・幾何学的特性と、それに基づく重心位置を定義するクラス
    
    properties (SetAccess = private)
        % 基本物理パラメータ
        M_payload, l_x, l_y, l_z
        M_canopy, b, c, t
        r_inner, r_center
        delta_x_offset, delta_y_offset, delta_z_offset
        n = 10000

        % 計算により導出されるプロパティ (ペイロード原点基準)
        TotalMass           
        r_payload_cm_abs
        r_canopy_cm_abs
        r_total_cm_abs
        r_canopy_origin_abs
        
        % 機体全体重心(CG)を基準とするベクトル
        r_total_cm_to_payload_cm_B
        r_total_cm_to_canopy_cm_B
        r_total_cm_to_canopy_origin_B
        
        % キャノピーのローカルな幾何ベクトル
        r_canopy_origin_to_ac_B
        
        % キャノピー翼弦上の幾何学的定義
        x_ac_from_te
        x_origin_from_te
        
        % キャノピー形状パラメータ
        catenary_a, catenary_C
    end
    
    methods
        function obj = ParafoilProperties(params)
            % コンストラクタ: パラメータを読み込み、重心などを自動計算
            fields = fieldnames(params);
            for i = 1:length(fields)
                prop_name = fields{i};
                if strcmp(prop_name, 'delta_x_pay'), prop_name = 'delta_x_offset'; end
                if strcmp(prop_name, 'delta_y_pay'), prop_name = 'delta_y_offset'; end
                if strcmp(prop_name, 'delta_z_pay'), prop_name = 'delta_z_offset'; end
                
                if isprop(obj, prop_name)
                    obj.(prop_name) = params.(fields{i});
                end
            end
            
            % 幾何学的定義をプロパティに格納
            obj.x_ac_from_te = 3 * obj.c / 4;
            obj.x_origin_from_te = obj.c / 4;

            % キャノピーの式及び重心を求めるに当たって、ペイロード原点からのオフセットを定義　キャノピーの重心を求める上での ペイロード原点(not重心)からの距離なのでobjには格納しない
            
            z_p1_magnitude = (obj.l_z / 2) + obj.r_center + obj.t / 2;
            r_canopy_offset = [0; 0; -z_p1_magnitude];
            

            %キャノピー座標系のペイロード原点からの位置ベクトル
            canopy_core_x = obj.c/2- obj.x_ac_from_te; 
            obj.r_canopy_origin_abs = [canopy_core_x; 0; -z_p1_magnitude];


            % キャノピー原点からACへのベクトルを計算
            ac_offset_x = obj.x_ac_from_te - obj.x_origin_from_te;
            obj.r_canopy_origin_to_ac_B = [ac_offset_x; 0; 0];
            
            % カテナリー形状を定義
            [obj.catenary_a, obj.catenary_C] = obj.calculate_catenary_params_internal(r_canopy_offset);
            
            % 各コンポーネントと全体の重心を計算
            obj.TotalMass = obj.M_payload + obj.M_canopy;
            obj.r_payload_cm_abs = [obj.delta_x_offset; obj.delta_y_offset; obj.delta_z_offset];
            obj.r_canopy_cm_abs = obj.calculate_canopy_cm_internal(r_canopy_offset);
            
            mass_moment = obj.M_payload * obj.r_payload_cm_abs + obj.M_canopy * obj.r_canopy_cm_abs;
            obj.r_total_cm_abs = mass_moment / obj.TotalMass;
            
            % 機体全体重心を基準とする相対ベクトルを計算
            obj.r_total_cm_to_payload_cm_B = obj.r_payload_cm_abs - obj.r_total_cm_abs;
            obj.r_total_cm_to_canopy_cm_B = obj.r_canopy_cm_abs - obj.r_total_cm_abs;
            obj.r_total_cm_to_canopy_origin_B = obj.r_canopy_origin_abs - obj.r_total_cm_abs;
        end
        % の end の直前に追加
 
        function visualize_geometry(obj)
            % ParafoilProperties で計算されたベクトルと重心を可視化するメソッド
            
            figure('Name', 'Parafoil Geometry Visualization (Payload Origin Ref.)');
            hold on;
            grid on;
            axis equal;
            xlabel('X [m] (Forward)');
            ylabel('Y [m] (Right)');
            zlabel('Z [m] (Up)'); % Z軸は上向きを正としてプロットします
            title('Parafoil Geometry (Payload Origin [0,0,0] Reference)');
            view(3);

            % --- 1. 主要な「点」をプロット (ペイロード原点 [0,0,0] 基準) ---
            
            % ペイロード原点 (基準)
            p_origin = [0; 0; 0];
            plot3(p_origin(1), p_origin(2), p_origin(3), 'k*', 'MarkerSize', 10, 'DisplayName', 'Payload Origin [0,0,0]');
            
            % ペイロード重心 (P_CM) [cite: 1881, 1929]
            p_cm = obj.r_payload_cm_abs;
            plot3(p_cm(1), p_cm(2), p_cm(3), 'bo', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Payload CM (abs)');

            % キャノピー原点 (C_Origin) [cite: 1884, 1921]
            c_origin = obj.r_canopy_origin_abs;
            plot3(c_origin(1), c_origin(2), c_origin(3), 'rs', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Canopy Origin (abs)');
            
            % キャノピー重心 (C_CM) [cite: 1882, 1930]
            c_cm = obj.r_canopy_cm_abs;
            plot3(c_cm(1), c_cm(2), c_cm(3), 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Canopy CM (abs)');

            % 機体全体重心 (Overall CG / B_Origin) [cite: 1883, 1933]
            b_origin = obj.r_total_cm_abs;
            plot3(b_origin(1), b_origin(2), b_origin(3), 'mx', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'OVERALL CG (B-Frame Origin)');

            % 空力中心 (AC) [cite: 1890, 1924]
            % (C_Origin から r_canopy_origin_to_ac_B だけ移動した点)
            % (注: GAMMA=0 の場合の位置関係として描画)
            ac_vec_b = obj.r_canopy_origin_to_ac_B;
            ac_abs = c_origin + ac_vec_b; % T_BCは単位行列と仮定
            plot3(ac_abs(1), ac_abs(2), ac_abs(3), 'gd', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Aerodynamic Center (abs, \Gamma=0)');
            
            
            % --- 2. 運動方程式で使われる「ベクトル」をプロット (Overall CG 基準) ---
            
            % (Overall CG) -> (Payload CM) [cite: 1886, 1935]
            vec_b_p = obj.r_total_cm_to_payload_cm_B;
            quiver3(b_origin(1), b_origin(2), b_origin(3), ...
                    vec_b_p(1), vec_b_p(2), vec_b_p(3), ...
                    'b', 'LineWidth', 1.5, 'DisplayName', 'r_{CG \rightarrow P_{CM}}');

            % (Overall CG) -> (Canopy Origin) [cite: 1888, 1939]
            vec_b_c = obj.r_total_cm_to_canopy_origin_B;
            quiver3(b_origin(1), b_origin(2), b_origin(3), ...
                    vec_b_c(1), vec_b_c(2), vec_b_c(3), ...
                    'r', 'LineWidth', 1.5, 'DisplayName', 'r_{CG \rightarrow C_{Origin}}');
            
            % (Canopy Origin) -> (AC) [cite: 1890, 1924]
            % (GAMMA=0と仮定し、C_Originから描画)
            quiver3(c_origin(1), c_origin(2), c_origin(3), ...
                    ac_vec_b(1), ac_vec_b(2), ac_vec_b(3), ...
                    'g', 'LineWidth', 1.5, 'DisplayName', 'r_{C_{Origin} \rightarrow AC} (\Gamma=0)');

            legend('show');
            hold off;
        end
    end % この 'methods' ブロックを追加


    
    methods (Access = private)
        function r_canopy_cm = calculate_canopy_cm_internal(obj, r_canopy_origin_abs)
            sum_mass_times_position = zeros(3, 1);
            m_element = obj.M_canopy / obj.n;
            for k = 1:obj.n
                y_element = -obj.b / 2 + (k - 0.5) * obj.b / obj.n;
                z_rel = obj.catenary_a * cosh(y_element / obj.catenary_a) + obj.catenary_C;
                r_element_vec = r_canopy_origin_abs + [0; y_element; z_rel];
                sum_mass_times_position = sum_mass_times_position + m_element * r_element_vec;
            end
            r_canopy_cm = sum_mass_times_position / obj.M_canopy;
        end

        function [a, C] = calculate_catenary_params_internal(obj, r_canopy_origin_abs)
            p1 = [0, 0];
            payload_half_diag = sqrt(obj.l_y^2 + obj.l_z^2) / 2;
          
    
            r_at_end_ref_magnitude = obj.r_inner + payload_half_diag + obj.t / 2;
            b_half = obj.b / 2;
            if r_at_end_ref_magnitude^2 < b_half^2, z_p2_abs_val = 0; else, z_p2_abs_val = -sqrt(r_at_end_ref_magnitude^2 - b_half^2); end
            p2_abs = [r_canopy_origin_abs(1); b_half; z_p2_abs_val];
            p2_rel = p2_abs - r_canopy_origin_abs;
            p2 = [p2_rel(2), p2_rel(3)];
            y1 = p1(1); z1 = p1(2); y2 = p2(1); z2 = p2(2);
            f = @(a_param) a_param * (cosh(y2/a_param) - cosh(y1/a_param)) - (z2 - z1);
            a0 = abs(y2-y1)/2; if a0==0, a0=1; end
            options = optimoptions('fsolve', 'Display','off');
            a = fsolve(f, a0, options);
            C = z1 - a*cosh(y1/a);
        end
    
    end
end