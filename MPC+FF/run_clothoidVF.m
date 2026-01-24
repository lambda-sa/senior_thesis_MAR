% run_clothoidVF.m
% 軌道計画 -> FF制御 -> 6自由度シミュレーション 統合実行スクリプト
clear; clc; close all;

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_ref.xlsx'; % パラメータファイル
wind_planned = [0;0];  % 風速 (North, East) [m/s]
wind_disturbance = [1;1];  % 予期せぬ外乱（未知の風）
wind_actual = wind_planned + wind_disturbance; % 実際に吹いている風

target_pos = [0, 600, 2743]; % 目標 [N, E, Alt]
L_final = 500;

% パラメータ読み込みとモデル構築
if exist('load_params_from_excel', 'file')
    [params, sim_settings] = load_params_from_excel(excelFileName);
else
    warning('load_params_from_excelが見つかりません。デフォルト値を使用します。');
    params = ParafoilParams(); 
    sim_settings = struct('X_initial', 0, 'Y_initial', 50, 'h_init', 1000, 't_max', 60, 'psi_initial_deg', 0);
end
% --- 6-DOFモデル用のパラメータ整形 ---
params.I_total_body = [
    params.I_xx,    0,          params.I_xz;
    0,              params.I_yy, 0;
    params.I_xz,    0,          params.I_zz
];

prop = struct();
prop.TotalMass = params.m_total;
prop.b = params.b;
prop.c = params.c;
prop.d = params.d;
prop.r_total_cm_to_canopy_origin_B = [params.Rcg_x; 0; params.Rcg_z];
prop.r_total_cm_to_payload_cm_B   = [params.Rpg_x; 0; params.Rpg_z];
prop.r_canopy_origin_to_ac_B      = [0; 0; 0]; 
params.prop = prop;

%% --- 2. Phase A: 軌道計画 (Mission & Planner) ---
fprintf('=== Phase A: Path Planning ===\n');
mission = ParafoilMissionClothoid(excelFileName);

% クロゾイド設定 (ロールレート10deg/s, 先行係数0.5, 助走30m)
mission.set_clothoid_options(0.5, 0.5, 0.0);

% シミュレーション実行 (内部で風補正計算が行われる)
mission.run_wind_simulation(target_pos, L_final, wind_planned);

% ★重要: 計画された「対気」軌道データを取得
% 6DOFは空気力を受けて飛ぶため、対気的なバンク角プロファイルが必要です。
trajPlanAir = mission.export_detailed_trajectory('Air');
trajPlanGnd = mission.export_detailed_trajectory('Ground'); % 比較用

t_plan = trajPlanAir.Time;
phi_plan = trajPlanAir.Euler_RPY(:, 1); % 目標バンク角

fprintf('  -> Planning Complete. Final Time: %.2f s\n', t_plan(end));

%% --- 3. Phase B & C: シミュレーション実行 (2ケース) ---
fprintf('\n=== Starting Comparison Simulation ===\n');

% --- Case 1: 無風 (No Disturbance) ---
% Truth: 計画風速, Est: 計画風速
fprintf('Running Case 1: No Disturbance...\n');
simData_NoWind = run_simulation_task(mission, params, trajPlanAir, target_pos, ...
                                     wind_planned, wind_disturbance);

% --- Case 2: 外乱あり (With Disturbance) ---
% Truth: 実際の風(外乱込み), Est: 計画風速(外乱を知らない)
wind_actual = wind_planned + wind_disturbance;
fprintf('Running Case 2: With Disturbance (Wind=[%.1f, %.1f])...\n', wind_actual(1), wind_actual(2));
simData_Wind = run_simulation_task(mission, params, trajPlanAir, target_pos, ...
                                   wind_actual, wind_disturbance);

%% --- 4. Phase D: 結果の比較・可視化 ---
fprintf('\n=== Phase D: Visualization ===\n');

% (1) フェーズ遷移時刻の抽出
plan_trans_times = [];
% データソースの取得
if isprop(mission.Planner, 'ResultDataWind') && ~isempty(mission.Planner.ResultDataWind)
    pd = mission.Planner.ResultDataWind;
else
    pd = mission.Planner.ResultData;
end

if ~isempty(pd)
    % (A) 基本フェーズ (Runup, Entry, Loiter)
    if isfield(pd, 'runup') && ~isempty(pd.runup.t), plan_trans_times(end+1) = pd.runup.t(end); end
    if isfield(pd, 'entry') && ~isempty(pd.entry.t), plan_trans_times(end+1) = pd.entry.t(end); end
    if isfield(pd, 'loiter') && ~isempty(pd.loiter.t), plan_trans_times(end+1) = pd.loiter.t(end); end
    
    % (B) Dubins区間 & 内部セグメント
    if isfield(pd, 'dubins') && ~isempty(pd.dubins.t)
        % Dubins終了 (=Final開始)
        plan_trans_times(end+1) = pd.dubins.t(end);
        
        % Dubins内部の切り替わり (Turn->Straight等)
        if isfield(pd.dubins, 'segment_end_idx')
            seg_idxs = pd.dubins.segment_end_idx;
            valid_idxs = seg_idxs(seg_idxs > 0 & seg_idxs < length(pd.dubins.t));
            
            if ~isempty(valid_idxs)
                internal_times = pd.dubins.t(valid_idxs);
                
                % ★ここをご提示のコード（安全版）にします
                plan_trans_times = [plan_trans_times, internal_times(:)'];
            end
        end
    end
    
    plan_trans_times = unique(plan_trans_times);
end
fprintf('Planned Transition Times: %s s\n', mat2str(plan_trans_times, 2));

% (2) 比較クラスのインスタンス化 (WindTrajectoryComparatorを使用)
comparator = WindTrajectoryComparator(...
    simData_NoWind, ...   % Case 1 (青)
    simData_Wind, ...     % Case 2 (マゼンタ)
    trajPlanGnd, ...
    trajPlanAir, ...
    target_pos, ...
    plan_trans_times);

% (3) 描画
comparator.plotAll();

fprintf('Comparison Done.\n');


function simData = run_simulation_task(mission, params, trajPlanAir, target_pos, wind_truth, wind_disturbance)
    % =========================================================
    % 相対運動シミュレーション (Moving Frame Approach)
    % =========================================================
    
    % 1. シミュレーション空間の風 = 外乱のみ
    wind_sim_env = wind_disturbance; 
    
    % 2. Plant & Autopilot 設定
    plant = ParafoilDynamics(params);
    if isprop(plant, 'WindVector')
        plant.WindVector = [wind_sim_env(1); wind_sim_env(2); 0];
    end
    
    autopilot = ParafoilAutopilot(mission.Params);
    autopilot.import_mission_data(mission.Planner);
    
    trim_ld = mission.PhysicsParams.glide_ratio;
    autopilot.set_longitudinal_control(true);
    autopilot.set_nominal_glide_ratio(trim_ld);
    
    % 3. Scheduler 設定
    scheduler = PlannerAutopilotScheduler(mission, autopilot);
    scheduler.WindVector_Truth = [wind_sim_env; 0]; 
    scheduler.WindVector_Est   = [0; 0];            
    
    % 4. 初期条件
    start_pos_ned = trajPlanAir.Position(1, :); 
    start_vel_tas = trajPlanAir.V_Air(1, :);   
    start_euler   = trajPlanAir.Euler_RPY(1, :); 
    
    y0 = [start_vel_tas(1); start_vel_tas(2); start_vel_tas(3); ...
          0; 0; 0; ...
          start_euler(1); start_euler(2); start_euler(3); ...
          start_pos_ned(1); start_pos_ned(2); -start_pos_ned(3)];
      
    % 5. エンジン設定
    dt_6dof = 0.05;
    t_max_6dof = trajPlanAir.Time(end) + 30; % ターゲット通過を見越して長めに
    engine = SimulationEngine(plant, scheduler, dt_6dof, t_max_6dof, y0);
    
    % ★重要: エンジン内部の自動停止機能はOFFにする
    % (理由: エンジン内は対気座標なので、地上の固定ターゲットとの距離が正しく測れないため)
    engine.TargetNED = []; 
    
    % 6. 実行 (最後まで計算させる)
    engine = engine.run();
    raw_simData = engine.create_results_table();
    
    % =========================================================
    % 7. 後処理: 「地上のターゲット」との距離判定 & トリミング
    % =========================================================
    
    % (A) 計画風成分の計算
    wind_planned = wind_truth - wind_disturbance;
    
    % (B) 生データ(対気)を対地座標に変換するための準備
    time_vec = raw_simData.Time;
    pos_air_N = raw_simData.Inertial_X_Position;
    pos_air_E = raw_simData.Inertial_Y_Position;
    pos_D     = raw_simData.Inertial_Z_Position;
    
    drift_N = wind_planned(1) * time_vec;
    drift_E = wind_planned(2) * time_vec;
    
    % (C) 対地位置の計算
    pos_gnd_N = pos_air_N + drift_N;
    pos_gnd_E = pos_air_E + drift_E;
    
    % (D) 「地上のターゲット」との距離を計算
    % target_pos = [North, East, Alt] -> NED=[N, E, -Alt]
    tg_N = target_pos(1);
    tg_E = target_pos(2);
    tg_D = -target_pos(3); 
    
    dist_sq = (pos_gnd_N - tg_N).^2 + ...
              (pos_gnd_E - tg_E).^2 + ...
              (pos_D     - tg_D).^2;
          
    dist_history = sqrt(dist_sq);
    
    % (E) 最接近点を探す
    [min_dist, idx_min] = min(dist_history);
    
    % (F) データを最接近点でカット (Trim)
    simData = raw_simData(1:idx_min, :);
    
    % (G) テーブル内の値を「対地座標系」に更新して確定させる
    simData.Inertial_X_Position = pos_gnd_N(1:idx_min);
    simData.Inertial_Y_Position = pos_gnd_E(1:idx_min);
    
    % 速度も対地速度に直す
    simData.Inertial_Vx = simData.Inertial_Vx + wind_planned(1);
    simData.Inertial_Vy = simData.Inertial_Vy + wind_planned(2);
    
    % 風速ログも全風速に直す
    if ismember('Wind_X', simData.Properties.VariableNames)
        simData.Wind_X = simData.Wind_X + wind_planned(1);
        simData.Wind_Y = simData.Wind_Y + wind_planned(2);
    end
    
    % 結果表示
    closest_t = simData.Time(end);
    fprintf('\n=== Ground Target Reached (Trimmed) ===\n');
    fprintf('  Time        : %.2f s\n', closest_t);
    fprintf('  Min Distance: %.2f m\n', min_dist);
    fprintf('  Ground Pos  : N=%.1f, E=%.1f, Alt=%.1f m\n', ...
        simData.Inertial_X_Position(end), ...
        simData.Inertial_Y_Position(end), ...
        -simData.Inertial_Z_Position(end));
    fprintf('=======================================\n');
end