classdef SimulationEngine
    % SIMULATIONENGINE 6自由度パラフォイルシミュレーション実行エンジン
    % 修正版: 制御入力のログ記録機能を追加
    
    properties (SetAccess = private)
        DynamicsModel
        ControlScheduler
        TimeStep
        MaxTime
        InitialConditions
        
        % 結果格納用
        TimeVector
        ResultsMatrix
        PhaseVector
        InertialVelMatrix
        EulRateMatrix
        
        % ★★★ 追加: 制御入力ログ用 ★★★
        % 列構成: [delta_R, delta_L, delta_s, delta_a, GAMMA, Wind_x, Wind_y, Wind_z]
        ControlInputMatrix 
    end
    
    methods
        function obj = SimulationEngine(dynamics_model, control_scheduler, h, t_max, y0)
            obj.DynamicsModel = dynamics_model;
            obj.ControlScheduler = control_scheduler;
            obj.TimeStep = h;
            obj.MaxTime = t_max;
            obj.InitialConditions = y0;
        end
    
        function obj = run(obj)
            % 時間ベクトルの生成
            t_vec = (0:obj.TimeStep:obj.MaxTime)';
            num_steps = length(t_vec);
            num_states = length(obj.InitialConditions);
    
            % 配列の事前確保
            y_res = zeros(num_steps, num_states);
            phase_res = zeros(num_steps, 1);
            vel_res = zeros(num_steps, 3);
            eul_rate_res = zeros(num_steps, 3);
            
            % ★★★ 追加: 制御入力ログ配列の確保 (8列) ★★★
            ctrl_res = zeros(num_steps, 8); 
    
            % 初期条件のセット
            y_res(1, :) = obj.InitialConditions';
            phase_res(1) = obj.ControlScheduler.phase;
            
            % 初期速度記録
            try
                initial_vel_I = obj.DynamicsModel.get_inertial_velocity(obj.InitialConditions);
                vel_res(1, :) = initial_vel_I';
            catch
                vel_res(1, :) = [NaN, NaN, NaN];
            end
            
            current_y = obj.InitialConditions;
            last_valid_step = 1;
    
            % --- メインループ ---
            try
                for i = 1:(num_steps - 1)
                    current_t = t_vec(i);
                    
                    % 1. スケジューラから制御入力を取得
                    % (ControlScheduler は delta_R, delta_L, GAMMA, wind_I を持つ構造体を返す想定)
                    control_inputs = obj.ControlScheduler.get_inputs(current_t, current_y, obj.TimeStep);
                    
                    % ★★★ 追加: 制御入力をログに記録 ★★★
                    % 構造体から値を抽出して計算
                    dR = control_inputs.delta_R;
                    dL = control_inputs.delta_L;
                    ds = min(dR, dL);     % 対称
                    da = dR - dL;         % 非対称
                    gam = control_inputs.GAMMA;
                    w = control_inputs.wind_I;
                    
                    % i番目の行に保存 (現在の入力)
                    ctrl_res(i, :) = [dR, dL, ds, da, gam, w(1), w(2), w(3)];
                    
                    
                    % 2. フェーズ記録
                    phase_res(i) = obj.ControlScheduler.phase; % 入力決定時のフェーズを記録
                    
                    % 3. ダイナミクス計算 (ルンゲクッタ)
                    func = @(t, y) obj.DynamicsModel.get_derivatives(t, y, control_inputs);
                    current_y = obj.runge_kutta_step(func, current_t, current_y, obj.TimeStep);
                    
                    % 4. 結果保存 (次のステップの状態)
                    y_res(i+1, :) = current_y';
                    last_valid_step = i + 1;
                    
                    % 5. 速度・角速度記録
                    try
                        vel_res(i+1, :) = (obj.DynamicsModel.get_inertial_velocity(current_y))';
                        eul_rate_res(i+1, :) = (obj.DynamicsModel.get_euler_rates(current_y))';
                    catch
                        vel_res(i+1, :) = [NaN, NaN, NaN];
                        eul_rate_res(i+1, :) = [NaN, NaN, NaN];
                    end
    
                    % 6. 発散・高度チェック
                    if any(isnan(current_y))
                        error('シミュレーションが発散しました (NaNを検出)。');
                    end
                    
                    % 高度チェック (-z_pos)
                    if -current_y(12) <= 1e-5
                        fprintf('地面に到達しました (t=%.2f s)。\n', t_vec(i+1));
                        % 最後のステップの入力も記録しておく（便宜上前ステップと同じ値）
                        ctrl_res(i+1, :) = ctrl_res(i, :);
                        phase_res(i+1) = obj.ControlScheduler.phase;
                        break; 
                    end
                end
                
                % 最後のステップの制御ログ埋め合わせ（ループを抜けた場合用）
                if ctrl_res(last_valid_step, 1) == 0 && last_valid_step > 1
                     ctrl_res(last_valid_step, :) = ctrl_res(last_valid_step-1, :);
                     phase_res(last_valid_step) = phase_res(last_valid_step-1);
                end
                
                disp('シミュレーション完了。');
    
            catch ME
                fprintf('\n!! エラーにより中断 (t=%.2f s) !!\n%s\n', t_vec(last_valid_step), ME.message);
                %rethrow(ME);
            end
            
            % --- 結果のトリミングと保存 ---
            valid_idx = 1:last_valid_step;
            
            obj.TimeVector = t_vec(valid_idx);
            obj.ResultsMatrix = y_res(valid_idx, :);
            obj.PhaseVector = phase_res(valid_idx);
            obj.InertialVelMatrix = vel_res(valid_idx, :);
            obj.EulRateMatrix = eul_rate_res(valid_idx, :);
            obj.ControlInputMatrix = ctrl_res(valid_idx, :); % ★★★ 格納 ★★★
        end
        
        % ★★★ 修正: create_results_table に制御入力を追加 ★★★
        function simDataTable = create_results_table(obj)
            t = obj.TimeVector;
            y = obj.ResultsMatrix;
            ph = obj.PhaseVector;
            vel = obj.InertialVelMatrix;
            eul = obj.EulRateMatrix;
            ctrl = obj.ControlInputMatrix;
            
            simDataTable = table(t, ...
                y(:, 10), y(:, 11), y(:, 12), ...           % 位置
                y(:, 1), y(:, 2), y(:, 3), ...             % 機体系速度
                y(:, 4), y(:, 5), y(:, 6), ...             % 角速度(p,q,r)
                rad2deg(y(:, 7)), rad2deg(y(:, 8)), rad2deg(y(:, 9)), ... % 姿勢角
                ph, ...                                     % フェーズ
                vel(:, 1), vel(:, 2), vel(:, 3), ...       % 慣性系速度
                eul(:, 1), eul(:, 2), eul(:, 3), ...       % オイラー角変化率
                ctrl(:, 1), ctrl(:, 2), ctrl(:, 3), ctrl(:, 4), rad2deg(ctrl(:, 5)), ... % 制御入力
                ctrl(:, 6), ctrl(:, 7), ctrl(:, 8), ...    % 風速
                'VariableNames', { ...
                'Time', ...
                'Inertial_X_Position', 'Inertial_Y_Position', 'Inertial_Z_Position', ...
                'Body_U_Vel', 'Body_V_Vel', 'Body_W_Vel', ...
                'Body_P_Vel', 'Body_Q_Vel', 'Body_R_Vel', ...
                'Roll_Angle', 'Pitch_Angle', 'Yaw_Angle', ...
                'Phase', ...
                'Inertial_Vx', 'Inertial_Vy', 'Inertial_Vz', ...
                'Eul_Phi_Dot', 'Eul_Theta_Dot', 'Eul_Psi_Dot', ...
                'delta_R', 'delta_L', 'delta_s', 'delta_a', 'GAMMA_deg', ...
                'Wind_X', 'Wind_Y', 'Wind_Z'});
        end
    end
    
    methods (Access = private)
        function y_next = runge_kutta_step(~, func, t, y, h)
            k1 = func(t, y);
            k2 = func(t + h/2, y + (h/2)*k1);
            k3 = func(t + h/2, y + (h/2)*k2);
            k4 = func(t + h, y + h*k3);
            y_next = y + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end
    end
end