classdef ConstantControlScheduler < handle
    % Excel等で指定された固定値を常に出力するスケジューラ
    
    properties
        InitialInputs % 構造体: delta_R, delta_L, GAMMA, wind_I
        phase = 1;    % 固定フェーズID
    end
    
    methods
        function obj = ConstantControlScheduler(params)
            % コンストラクタ: params構造体から初期値を読み取る
            
            % デフォルト値の設定 (Excelに無い場合用)
            defaults.delta_R_init = 0;
            defaults.delta_L_init = 0;
            defaults.ang_gamma_rigging = 0; % GAMMAのこと (rad単位である前提)
            defaults.wind_x = 0;
            defaults.wind_y = 0;
            defaults.wind_z = 0;
            
            % paramsの内容で上書き
            fnames = fieldnames(defaults);
            for i = 1:length(fnames)
                fn = fnames{i};
                if isfield(params, fn)
                    val = params.(fn);
                else
                    val = defaults.(fn);
                end
                % 内部プロパティに保存
                obj.InitialInputs.(fn) = val;
            end
        end
        
        function inputs = get_inputs(obj, ~, ~, ~)
            % 常に初期値を返す (時刻tや状態yによらない)
            
            inputs.delta_R = obj.InitialInputs.delta_R_init;
            inputs.delta_L = obj.InitialInputs.delta_L_init;
            inputs.GAMMA   = obj.InitialInputs.ang_gamma_rigging;
            
            % 風速ベクトル
            inputs.wind_I  = [obj.InitialInputs.wind_x;
                              obj.InitialInputs.wind_y;
                              obj.InitialInputs.wind_z];
        end
    end
end