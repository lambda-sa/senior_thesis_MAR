classdef ParafoilLinearizerNumeric
    properties
        model       % ParafoilDynamics のインスタンス
        epsilon     % 摂動の微小量
        
        % --- 追加: 線形化時に適用したいパラメータ ---
        soft_min_k  % 平滑化係数 (例: 20~50)
    end
    
    methods
        function obj = ParafoilLinearizerNumeric(parafoil_dynamics_instance)
            obj.model = parafoil_dynamics_instance;
            obj.epsilon = 1e-5;
            obj.soft_min_k = 50; % デフォルト値
        end
        
        function [A, B] = get_linear_model(obj, t, y_trim, u_trim_vec, extra_params)
            
            % 1. 現在のモデル設定をバックアップ (退避)
            original_use_soft = obj.model.use_soft_min;
            original_k = obj.model.soft_min_k;
            
            % 2. 線形化用に設定を上書き (SoftMinモードをON)
            obj.model.use_soft_min = true;
            obj.model.soft_min_k = obj.soft_min_k; % Linearizerのk値を適用
            
            try
                % 3. 数値微分の実行 (この間だけ SoftMin が使われる)
                [A, B] = obj.calculate_jacobian(t, y_trim, u_trim_vec, extra_params);
                
            catch ME
                % エラーが起きても設定を必ず元に戻すための try-catch
                obj.model.use_soft_min = original_use_soft;
                obj.model.soft_min_k = original_k;
                rethrow(ME);
            end
            
            % 4. 設定を元に戻す (復元)
            % これにより、シミュレーション時は元の min 関数に戻る
            obj.model.use_soft_min = original_use_soft;
            obj.model.soft_min_k = original_k;
        end

        function input_struct = vector_to_struct(~, u_vec, extra_params)
            input_struct.delta_L = u_vec(1);
            input_struct.delta_R = u_vec(2);
            input_struct.GAMMA = extra_params.GAMMA;
            input_struct.wind_I = extra_params.wind_I;
        end
        
    end

    
    
    methods (Access = private)
        function [A, B] = calculate_jacobian(obj, t, y_trim, u_trim_vec, extra_params)
            % (前回提示した数値微分のループ処理をここに記述)
            n_states = 12;
            n_inputs = 2;
            
            A = zeros(n_states, n_states);
            B = zeros(n_states, n_inputs);
            
            % --- A行列の計算 ---
            % (前回のコードと同じ)
            for i = 1:n_states
                y_plus = y_trim; y_minus = y_trim;
                y_plus(i) = y_trim(i) + obj.epsilon;
                y_minus(i) = y_trim(i) - obj.epsilon;
                
                input_struct = obj.vector_to_struct(u_trim_vec, extra_params);
                
                dy_plus = obj.model.get_derivatives(t, y_plus, input_struct);
                dy_minus = obj.model.get_derivatives(t, y_minus, input_struct);
                A(:, i) = (dy_plus - dy_minus) / (2 * obj.epsilon);
            end
            
            % --- B行列の計算 ---
            for j = 1:n_inputs
                u_plus = u_trim_vec; u_minus = u_trim_vec;
                u_plus(j) = u_trim_vec(j) + obj.epsilon;
                u_minus(j) = u_trim_vec(j) - obj.epsilon;
                
                struct_plus = obj.vector_to_struct(u_plus, extra_params);
                struct_minus = obj.vector_to_struct(u_minus, extra_params);
                
                dy_plus = obj.model.get_derivatives(t, y_trim, struct_plus);
                dy_minus = obj.model.get_derivatives(t, y_trim, struct_minus);
                B(:, j) = (dy_plus - dy_minus) / (2 * obj.epsilon);
            end
        end
        
        
    end
end