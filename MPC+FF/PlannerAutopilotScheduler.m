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
        
        function inputs = get_inputs(obj, t, y, ~)
            % 1. モード切替判定
            current_z = y(12);
            if strcmp(obj.Autopilot.CurrentMode, 'Loiter') && current_z > obj.MissionSwitchAlt
                obj.Autopilot.set_mode('Mission');
            end
            
            
            % ====================================================
            % ★★★ キネマティクスによる対気状態推定 ★★★
            % ====================================================
            
            % (1) センサ値の取得
            V_g_Body = y(1:3); % 対地速度(機体軸) [u; v; w]
            phi=y(7); theta=y(8); psi=y(9);
            
            % (2) 回転行列 (NED -> Body)
            % ※ RotationMatricesクラスがない場合は自前で定義してください
            T_IB = RotationMatrices.get_inertial_to_body_matrix(phi, theta, psi);
            
            % (3) 風ベクトルの座標変換 (NED -> Body)
            % Schedulerが持っている「推定風」を使う
            wind_Body = T_IB * obj.WindVector_Est;
            
            % (4) 対気速度ベクトルの計算 (Body Frame)
            % V_air = V_ground - Wind
            V_air_Body = V_g_Body - wind_Body;
            
            u_a = V_air_Body(1);
            w_a = V_air_Body(3);
            V_tas_est = norm(V_air_Body);
            
            % (5) 迎角 alpha の推定 (幾何学的定義)
            alpha_est = atan2(w_a, u_a);
            
            % (6) 経路角 gamma と 水平速度 V_horiz の推定
            % 理論: theta = gamma + alpha
            gamma_est = theta - alpha_est;
            
            % これが求めたかった「真の水平対気速度」
            V_horiz_est = V_tas_est * cos(gamma_est);
            
            % ゼロ割防止
            if V_horiz_est < 0.1, V_horiz_est = 0.1; end
            
            % ====================================================
            % Autopilotへ推定値を渡す
            % ====================================================
            
            % ★修正: 推定した V_horiz_est を渡す
            [dR, dL, ~] = obj.Autopilot.update(t, y, obj.WindVector_Est, V_horiz_est);
            
            % ★★★ 変更点 2: 物理モデルには「真の風 (Truth)」を渡す ★★★
            inputs.delta_R = dR;
            inputs.delta_L = dL;
            inputs.wind_I  = obj.WindVector_Truth; % 機体を流すのはこの風
            inputs.GAMMA   = 0;
            
            inputs.delta_a_cmd = dR - dL;
            inputs.delta_s_cmd = min(dR, dL);

            inputs.mu_cmd = 0;
        end
    end
end