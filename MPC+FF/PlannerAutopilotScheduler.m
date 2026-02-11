classdef PlannerAutopilotScheduler < handle
    % PLANNERAUTOPILOTSCHEDULER (Unified Mode Version)
    %
    % 概要:
    %   モード切替(Loiter/Mission)ロジックを廃止。
    %   一本化された軌道を追従する Autopilot に対して、
    %   シンプルに時刻と状態を渡して制御量を受け取るブリッジクラス。
    
    properties
        Autopilot
        WindVector_Truth % 物理モデルに渡す「真の風」 (外乱込み)
        WindVector_Est   % 制御器に渡す「推定風」 (計画値など)
        phase=1;
    end
    
    methods
        function obj = PlannerAutopilotScheduler(mission, autopilot_instance)
            obj.Autopilot = autopilot_instance;
            
            % 1. 風設定 (デフォルト値)
            % 実行スクリプト側で run_simulation_task 内で上書きされることを想定
            w_init = [0; 0; 0];
            if isprop(mission, 'Planner') && ~isempty(mission.Planner.WindVector)
                w_vec = mission.Planner.WindVector;
                w_init = [w_vec(1); w_vec(2); 0];
            end
            obj.WindVector_Truth = w_init;
            obj.WindVector_Est   = w_init;
            
            % 2. Autopilotにデータをロード
            % Plannerが生成した連続軌道データ(RunUp->Loiter->Dubins->Final)を一括で渡す
            if isprop(mission, 'Planner')
                obj.Autopilot.import_mission_data(mission.Planner);
            end
            
            fprintf('[Scheduler] Initialized for Unified Trajectory Following.\n');
        end
        
        function inputs = get_inputs(obj, t, y, ~)
            % --- 制御量計算 ---
            % Autopilotは現在位置に基づいて自動的にパス上の目標点を探索するため、
            % 外部からのモード切替指示は不要です。
            
            % Autopilotには「推定風 (Est)」を渡す
            % (外乱シミュレーション時: Estは計画風、Truthは実風)
            [dR, dL, ~] = obj.Autopilot.update(t, y, obj.WindVector_Est, []);
            
            % --- Dynamicsへの入力作成 ---
            inputs.delta_R = dR;
            inputs.delta_L = dL;
            
            % 物理モデル (Dynamics) には「真の風 (Truth)」を渡す
            % これにより、機体は「実風」に流される挙動をする
            inputs.wind_I  = obj.WindVector_Truth; 
            
            % 補助入力 (必要に応じて)
            inputs.GAMMA   = 0;
            inputs.delta_a_cmd = dR - dL;
            inputs.delta_s_cmd = min(dR, dL);
            inputs.mu_cmd = 0;
        end
    end
end