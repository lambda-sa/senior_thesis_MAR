classdef PlannerAutopilotScheduler < handle
    % PLANNERTRACKINGSCHEDULER (Autopilot Adapter Version)
    %
    % 役割:
    % シミュレーションの各ステップで、ParafoilAutopilot クラスを呼び出し、
    % 現在の状態に応じた操作量 (delta_R, delta_L) を計算して返すアダプター。
    %
    % ※ 以前のような「時刻tに応じた参照値の補間」は行わず、
    %    「現在位置」に基づくVF誘導を Autopilot に委譲する。
    
    properties
        Autopilot       % ParafoilAutopilot インスタンス (旧 Mapper)

        % ★★★ 追加: 風ベクトルを2つに分離 ★★★
        WindVector_Truth  % [Wx; Wy; 0] (環境の真の風 -> 物理モデルへ)
        WindVector_Est    % [Wx; Wy; 0] (推定された風 -> Autopilotへ)
        % ---------------------------------------
        phase =1;
        
        % フェーズ管理用パラメータ
        MissionSwitchAlt = -300; % [m] Missionモードへ切り替える高度(負の値)
                                 % 例: -300m (高度300m) を切ったら着陸進入
    end
    
    methods
        function obj = PlannerAutopilotScheduler(mission, autopilot_instance)
            obj.Autopilot = autopilot_instance;
            
            % ★★★ 追加: 風プロパティの初期化 ★★★
            % デフォルト値として、Plannerが持っている風(予測値)を入れておく
            w_init = [0; 0; 0];
            if isprop(mission, 'Planner') && ~isempty(mission.Planner.WindVector)
                w_vec = mission.Planner.WindVector;
                if length(w_vec) == 2
                    w_init = [w_vec(1); w_vec(2); 0];
                else
                    w_init = w_vec; % 既に3次元ならそのまま
                end
            end
            
            obj.WindVector_Truth = w_init;
            obj.WindVector_Est   = w_init;
            % ---------------------------------------
            
            obj.Autopilot.import_mission_data(mission.Planner);
            obj.Autopilot.set_mode('Loiter');
        end
        
        function inputs = get_inputs(obj, ~, y, ~)
            % (中略: フェーズ切り替えロジック等はそのまま)
            
            % ★★★ 変更点 1: Autopilotには「推定風 (Est)」を渡す ★★★
            % これにより、外乱込みの風をセットしていればクラブ角補正が効く
            [dR, dL, ~] = obj.Autopilot.update(y, obj.WindVector_Est, 0);
            
            % ★★★ 変更点 2: 物理モデルには「真の風 (Truth)」を渡す ★★★
            inputs.delta_R = dR;
            inputs.delta_L = dL;
            inputs.wind_I  = obj.WindVector_Truth; % 機体を流すのはこの風
            inputs.GAMMA   = 0;
            
            inputs.delta_a_cmd = dR - dL;
            inputs.delta_s_cmd = min(dR, dL);
        end
    end
end