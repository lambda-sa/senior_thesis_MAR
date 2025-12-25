classdef PlannerTrackingScheduler < handle
    % PLANNERTRACKINGSCHEDULER
    % 修正版: バンク角だけでなく、速度とピッチ角も計画値を使用する
    
    properties
        Mapper
        RefTime
        RefBankAngle
        
        % ★★★ 追加 ★★★
        RefTrueAirspeed % 計画されたTAS [m/s]
        RefPitchAngle   % 計画されたピッチ角 [rad]
        
        WindVector
        phase = 1
    end
    
    methods
        % ★★★ コンストラクタ引数を追加 ★★★
        function obj = PlannerTrackingScheduler(mapper, t_ref, phi_ref, V_ref, theta_ref, wind_vec)
            obj.Mapper = mapper;
            obj.RefTime = t_ref;
            obj.RefBankAngle = phi_ref;
            
            % 追加データの保存
            obj.RefTrueAirspeed = V_ref;
            obj.RefPitchAngle = theta_ref;
            
            if length(wind_vec) == 2
                obj.WindVector = [wind_vec(1); wind_vec(2); 0];
            else
                obj.WindVector = wind_vec;
            end
        end
        
        function inputs = get_inputs(obj, t, ~, ~) % 第2引数(y)は使いません
            % 1. 各参照値を線形補間
            if t > obj.RefTime(end)
                phi_cmd   = obj.RefBankAngle(end);
                V_cmd     = obj.RefTrueAirspeed(end);
                theta_cmd = obj.RefPitchAngle(end);
            else
                % まとめて補間 (interp1は行列も扱えるがわかりやすく個別に記述)
                phi_cmd   = interp1(obj.RefTime, obj.RefBankAngle, t, 'linear', 'extrap');
                V_cmd     = interp1(obj.RefTime, obj.RefTrueAirspeed, t, 'linear', 'extrap');
                theta_cmd = interp1(obj.RefTime, obj.RefPitchAngle, t, 'linear', 'extrap');
            end
            
            % 2. Mapperに「計画値」を渡して計算
            [dR, dL, ~] = obj.Mapper.compute_input_from_reference(phi_cmd, V_cmd, theta_cmd);
            
            % 3. パッキング
            inputs.delta_R = dR;
            inputs.delta_L = dL;
            inputs.delta_a_cmd = dR - dL;
            inputs.delta_s_cmd = min(dR, dL);
            inputs.mu_cmd = 0;
            inputs.GAMMA = 0;
            inputs.wind_I = obj.WindVector;
        end
    end
end