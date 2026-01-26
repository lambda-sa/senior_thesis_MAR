classdef PlannerAutopilotScheduler < handle
    % PLANNERAUTOPILOTSCHEDULER (Debugged Version)
    
    properties
        Autopilot
        WindVector_Truth
        WindVector_Est
        MissionSwitchAlt  % NED座標系 (負の値)
        phase=1;
    end
    
    methods
        function obj = PlannerAutopilotScheduler(mission, autopilot_instance)
            obj.Autopilot = autopilot_instance;
            
            % 1. 風設定
            w_init = [0; 0; 0];
            if isprop(mission, 'Planner') && ~isempty(mission.Planner.WindVector)
                w_vec = mission.Planner.WindVector;
                w_init = [w_vec(1); w_vec(2); 0];
            end
            obj.WindVector_Truth = w_init;
            obj.WindVector_Est   = w_init;
            
            % 2. Autopilot初期化
            if isprop(mission, 'Planner')
                obj.Autopilot.import_mission_data(mission.Planner);
            end
            obj.Autopilot.set_mode('Loiter');
            
            % 3. デフォルト高度設定 (万が一 Dubins設定が失敗したとき用)
            obj.MissionSwitchAlt = -300; 
            
            % 4. Dubinsから設定上書き
            obj.setup_loiter_from_dubins(mission);
        end
        
        function inputs = get_inputs(obj, t, y, ~)
            current_z = y(12); 
            
            % --- フェーズ管理 (変更なし) ---
            if strcmp(obj.Autopilot.CurrentMode, 'Loiter')
                if current_z > obj.MissionSwitchAlt
                    obj.Autopilot.set_mode('Mission');
                    fprintf('>>> PHASE CHANGE: Loiter -> MISSION (Alt: %.1f m) <<<\n', -current_z);
                end
            end
            
            % --- 制御量計算 ---
            % ★修正1: Autopilotには「推定風 (Est)」を渡す
            % (run_simulation_task で Est=[0;0] に設定されていれば、カニ歩きせず軌道誤差のみで制御する)
            [dR, dL, ~] = obj.Autopilot.update(t, y, obj.WindVector_Est, []);
            
            inputs.delta_R = dR;
            inputs.delta_L = dL;
            
            % ★修正2: 物理モデル (Dynamics) には「真の風 (Truth)」を渡す
            % ParafoilDynamics は control_inputs.wind_I を見て流される計算を行う
            inputs.wind_I  = obj.WindVector_Truth; 
            
            inputs.GAMMA   = 0;
            inputs.delta_a_cmd = dR - dL;
            inputs.delta_s_cmd = min(dR, dL);
            inputs.mu_cmd = 0;
        end
        
        function setup_loiter_from_dubins(obj, mission)
            if ~isprop(mission, 'Planner') || ~isfield(mission.Planner.ResultData, 'dubins')
                return;
            end
            d_res = mission.Planner.ResultData.dubins;
            if isempty(d_res.x), return; end
            
            fprintf('[Scheduler] Setting up Loiter from Dubins...\n');
            
            % --- Loiterパラメータ ---
            xs = d_res.x(1); ys = d_res.y(1); psis = d_res.psi(1);
            
            if isprop(mission.Planner, 'R_fixed') && ~isempty(mission.Planner.R_fixed)
                R_turn = mission.Planner.R_fixed;
            elseif isprop(mission.Planner, 'R_min')
                R_turn = mission.Planner.R_min;
            else
                R_turn = 50.0; 
            end
            
            lambda = 1;
            if isfield(d_res, 'modes') && ~isempty(d_res.modes)
                if strcmp(d_res.modes{1}, 'L'), lambda = -1; end
            else
                idx_chk = min(5, length(d_res.psi));
                d_psi = mod(d_res.psi(idx_chk) - psis + pi, 2*pi) - pi;
                if d_psi < -1e-3, lambda = -1; end
            end
            
            xc = xs + R_turn * cos(psis + lambda * pi/2);
            yc = ys + R_turn * sin(psis + lambda * pi/2);
            
            obj.Autopilot.LoiterParams = struct('xc',xc, 'yc',yc, 'R',R_turn, 'lambda',lambda);
            
            % --- 切り替え高度 ---
            if isfield(d_res, 'segment_end_idx') && ~isempty(d_res.segment_end_idx) && d_res.segment_end_idx(1) > 0
                idx_switch = d_res.segment_end_idx(1);
            else
                % フォールバック: 変化率がゼロになる点
                dp = diff(unwrap(d_res.psi));
                idx_switch = find(abs(dp) < 1e-4, 1, 'first');
                if isempty(idx_switch), idx_switch = floor(length(d_res.z)*0.3); end
            end
            
            idx_switch = max(1, min(idx_switch, length(d_res.z)));
            raw_z = d_res.z(idx_switch);
            
            % ★最重要修正: 必ず負の値にする (NED座標対応)
            obj.MissionSwitchAlt = -abs(raw_z);
            
            fprintf('   -> Loiter R=%.1f, Dir=%d, Switch Alt=%.1f m (NED z > %.1f)\n', ...
                    R_turn, lambda, abs(raw_z), obj.MissionSwitchAlt);
        end
    end
end