classdef ParafoilInertiaCalculator
    % ParafoilPropertiesに基づき、慣性テンソルのみを計算するクラス
    
    properties (SetAccess = private)
        I_payload_total_cm  % 機体全体重心周りの、ペイロード部分の慣性テンソル
        I_total_cm          % 機体全体重心周りの、機体全体の慣性テンソル
        I_canopy_transformed % 機体全体重心周りの、キャノピー部分の慣性テンソル

        I_canopy_own_cm     % 【New】キャノピー単体重心周りの、キャノピーの慣性テンソル
    end
    
    methods
        function obj = ParafoilInertiaCalculator(props)
            % コンストラクタ: ParafoilPropertiesオブジェクトを受け取り、慣性計算を実行
            if ~isa(props, 'ParafoilProperties')
                error('Input must be an object of the ParafoilProperties class.');
            end
            
            % 1. ペイロード部分の慣性テンソルを計算
            I_payload_local = obj.get_payload_inertia_local(props);
            r_vec_p = props.r_total_cm_to_payload_cm_B;
            obj.I_payload_total_cm = obj.apply_parallel_axis_theorem(I_payload_local, props.M_payload, r_vec_p);

            % 2. キャノピー部分の慣性テンソルを計算
            obj.I_canopy_transformed = obj.get_canopy_inertia_transformed(props);

            % 3. 全体の慣性テンソルを合算
            obj.I_total_cm = obj.I_payload_total_cm + obj.I_canopy_transformed;
            obj.I_total_cm = diag(diag(obj.I_total_cm));

            % 4. 【New】キャノピー単体重心周りの慣性テンソルを計算
            obj.I_canopy_own_cm = obj.get_canopy_inertia_about_own_cm(props);
        end
    end
    
    methods (Access = private)
        function I_local = get_payload_inertia_local(~, props)
            I_xx = (1/12) * props.M_payload * (props.l_y^2 + props.l_z^2);
            I_yy = (1/12) * props.M_payload * (props.l_x^2 + props.l_z^2);
            I_zz = (1/12) * props.M_payload * (props.l_x^2 + props.l_y^2);
            I_local = diag([I_xx, I_yy, I_zz]);
        end
        
        function I_canopy_transformed = get_canopy_inertia_transformed(obj, props)
            I_canopy_transformed = zeros(3,3);
            m_element = props.M_canopy / props.n;
            
            a_h = props.c/2; b_h = props.b/(2*props.n); c_h = props.t/2;
            I_element_local = (1/12) * m_element * diag([(2*b_h)^2+(2*c_h)^2, (2*a_h)^2+(2*c_h)^2, (2*a_h)^2+(2*b_h)^2]);
            
            % キャノピー原点の絶対位置を取得
            z_p1_magnitude = (props.l_z / 2) + props.r_center + props.t / 2;
            r_canopy_origin_abs = [0; 0; -z_p1_magnitude];

            for k = 1:props.n
                y_element = -props.b/2 + (k-0.5)*props.b/props.n;
                z_rel = props.catenary_a * cosh(y_element / props.catenary_a) + props.catenary_C;
                r_element_abs = r_canopy_origin_abs + [0; y_element; z_rel];
                r_vec_c = r_element_abs - props.r_total_cm_abs;
                
                I_element_transformed = obj.apply_parallel_axis_theorem(I_element_local, m_element, r_vec_c);
                I_canopy_transformed = I_canopy_transformed + I_element_transformed;
            end
        end

        function I_transformed = apply_parallel_axis_theorem(~, I_local_cm, mass, r_vec)
            r_sq = r_vec' * r_vec;
            r_outer_prod = r_vec * r_vec';
            I_transformed = I_local_cm + mass * (r_sq * eye(3) - r_outer_prod);
        end

        function I_canopy_own = get_canopy_inertia_about_own_cm(obj, props)
            % 【修正版】キャノピー単体重心周りの慣性テンソル
            % ペイロードの位置情報(l_zなど)には一切依存せず、キャノピー形状のみで計算する
            
            I_canopy_own = zeros(3,3);
            m_element = props.M_canopy / props.n;
            
            % 1. 要素単体の慣性モーメント (形状のみに依存)
            a_h = props.c/2; b_h = props.b/(2*props.n); c_h = props.t/2;
            I_element_local = (1/12) * m_element * diag([(2*b_h)^2+(2*c_h)^2, (2*a_h)^2+(2*c_h)^2, (2*a_h)^2+(2*b_h)^2]);
            
            % 2. キャノピーの「ローカル重心」を再計算
            % (props.r_canopy_cm_abs はペイロード基準のオフセットを含むため、ここでは使わない)
            r_cm_local = obj.calculate_local_canopy_cm(props);
            
            % 3. 各要素の寄与を合算
            for k = 1:props.n
                % キャノピー基準(Local)での要素位置を取得
                r_element_local = obj.calculate_element_pos_local(props, k);
                
                % ベクトル: ローカル要素位置 - ローカル重心
                % (純粋にキャノピー内部の距離関係のみになる)
                r_vec_c = r_element_local - r_cm_local;
                
                % 平行軸の定理
                r_sq = r_vec_c' * r_vec_c;
                r_outer_prod = r_vec_c * r_vec_c';
                I_element_shifted = I_element_local + m_element * (r_sq * eye(3) - r_outer_prod);
                
                I_canopy_own = I_canopy_own + I_element_shifted;
            end
        end
    
        function r_element_local = calculate_element_pos_local(~, props, k)
            % キャノピーのアーチ頂点(または基準点)を原点としたローカル座標
            % ペイロードからの距離オフセット(z_p1_magnitude)は加算しない
            y_element = -props.b/2 + (k-0.5)*props.b/props.n;
            z_rel = props.catenary_a * cosh(y_element / props.catenary_a) + props.catenary_C;
            
            % ローカル座標: [0; y; z_rel]
            r_element_local = [0; y_element; z_rel];
        end
    
        function r_cm_local = calculate_local_canopy_cm(obj, props)
            % キャノピー単体のローカル座標系における重心位置を計算
            sum_pos = [0; 0; 0];
            for k = 1:props.n
                sum_pos = sum_pos + obj.calculate_element_pos_local(props, k);
            end
            r_cm_local = sum_pos / props.n; % 質量は一様なので個数で割るのと同義
        end

    end
end