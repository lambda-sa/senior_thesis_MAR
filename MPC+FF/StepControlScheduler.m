classdef StepControlScheduler < handle
    % STEPCONTROLSCHEDULER 指定時刻にステップ入力を入れるスケジューラ
    %
    % 動作:
    %   - リギング角 (GAMMA): 初期設定値で固定
    %   - 右ブレーキ (delta_R): 常に 0
    %   - 左ブレーキ (delta_L): t < 50s までは 0, t >= 50s で 0.2
    %   - 風 (Wind): 初期設定値で固定
    
    properties
        InitialGamma  % リギング角 [rad]
        WindVector    % 風速ベクトル [m/s]
        phase = 1;    % 現在のフェーズ (1:待機, 2:操舵)
    end
    
    methods
        function obj = StepControlScheduler(sim_settings)
            % コンストラクタ: sim_settings から初期値を読み込む
            
            % リギング角 (なければ0)
            if isfield(sim_settings, 'ang_gamma_rigging')
                obj.InitialGamma = sim_settings.ang_gamma_rigging;
            else
                obj.InitialGamma = 0;
            end
            
            % 風速ベクトル
            wx = 0; wy = 0; wz = 0;
            if isfield(sim_settings, 'wind_x'), wx = sim_settings.wind_x; end
            if isfield(sim_settings, 'wind_y'), wy = sim_settings.wind_y; end
            if isfield(sim_settings, 'wind_z'), wz = sim_settings.wind_z; end
            obj.WindVector = [wx; wy; wz];
        end
        
        function inputs = get_inputs(obj, t, ~, ~)
            % シミュレーション時刻 t に応じた制御入力を返す
            
            % --- 固定値の設定 ---
            inputs.delta_R = 0;
            inputs.GAMMA   = obj.InitialGamma;
            inputs.wind_I  = obj.WindVector;
            
            % --- ステップ入力のロジック ---
            if t < 50
                % 50秒未満: 操作なし
                inputs.delta_L = 0;
                obj.phase = 1; % Phase 1: 直進/待機
            else
                % 50秒以降: 左旋回入力
                inputs.delta_L = 0.2;
                obj.phase = 2; % Phase 2: 左旋回
            end
        end
    end
end