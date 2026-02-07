% run_clothoidVF.m
% 軌道計画 -> FF制御 -> 6自由度シミュレーション 統合実行スクリプト
clear; clc; close all;

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_ref.xlsx'; % パラメータファイル
wind_planned = [0;0];  % 風速 (North, East) [m/s]
wind_disturbance = [0;0];  % 予期せぬ外乱（未知の風）
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

% --- Case 2: リプランニングあり (With Replanning) ---
wind_actual = wind_planned + wind_disturbance;
replan_time = 30.0; % ★追加: 再計画を行う時刻 [s] (ロイター中など)

fprintf('\nRunning Case 2: Replanning Mode (Wind Actual=[%.1f, %.1f])...\n', wind_actual(1), wind_actual(2));

% ★変更: 戻り値を2つ受け取る
[simData_Replanned, trajPlan_Replanned] = run_simulation_with_replanning(...
    mission, params, target_pos, ...
    wind_planned, ...
    wind_actual, ...
    replan_time);


%% --- 結果の可視化 (リプランニング効果の確認) ---
figure('Name', 'Replanning Result', 'Color', 'w');

% 1. 3D軌跡プロット
subplot(1, 2, 1); hold on; grid on; axis equal;
xlabel('North [m]'); ylabel('East [m]'); zlabel('Altitude [m]');
title('3D Trajectory: Replanning vs Actual');

% (A) 初期計画 (点線・薄い青)
plot3(trajPlanGnd.Position(:,2), trajPlanGnd.Position(:,1), trajPlanGnd.Position(:,3), ...
    'b--', 'LineWidth', 1, 'DisplayName', 'Initial Plan');

% (B) リプラン後の計画 (実線・緑)
% ※リプラン時刻以降のみ有効な軌道
plot3(trajPlan_Replanned.Position(:,2), trajPlan_Replanned.Position(:,1), trajPlan_Replanned.Position(:,3), ...
    'g-', 'LineWidth', 2, 'DisplayName', 'Replanned Path');

% (C) 実際の飛行軌跡 (マゼンタ)
plot3(simData_Replanned.Inertial_Y_Position, simData_Replanned.Inertial_X_Position, -simData_Replanned.Inertial_Z_Position, ...
    'm-', 'LineWidth', 1.5, 'DisplayName', 'Actual Flight');

% (D) ターゲットとリプラン地点
plot3(target_pos(2), target_pos(1), target_pos(3), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'Target');
    
% リプラン地点をマーク (時刻で検索)
[~, idx_replan] = min(abs(simData_Replanned.Time - replan_time));
pos_replan = [simData_Replanned.Inertial_Y_Position(idx_replan), ...
              simData_Replanned.Inertial_X_Position(idx_replan), ...
              -simData_Replanned.Inertial_Z_Position(idx_replan)];
plot3(pos_replan(1), pos_replan(2), pos_replan(3), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'y', 'DisplayName', 'Replan Point');

legend('Location', 'best');
view(-45, 30);


% 2. 2D平面図 (上空視点)
subplot(1, 2, 2); hold on; grid on; axis equal;
xlabel('East [m]'); ylabel('North [m]');
title('Top View (Ground Path)');

% (A) 初期計画
plot(trajPlanGnd.Position(:,2), trajPlanGnd.Position(:,1), 'b--', 'DisplayName', 'Initial Plan');

% (B) リプラン後の計画
plot(trajPlan_Replanned.Position(:,2), trajPlan_Replanned.Position(:,1), 'g-', 'LineWidth', 2, 'DisplayName', 'Replanned Path');

% (C) 実際の飛行軌跡
plot(simData_Replanned.Inertial_Y_Position, simData_Replanned.Inertial_X_Position, 'm-', 'LineWidth', 1.5, 'DisplayName', 'Actual Flight');

% (D) ターゲット
plot(target_pos(2), target_pos(1), 'rp', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
plot(pos_replan(1), pos_replan(2), 'ko', 'MarkerFaceColor', 'y');

legend('Location', 'best');
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
    simData_Replanned, ...     % Case 2 (マゼンタ)
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

% =========================================================================
%  以下、新規追加する関数群 (リプランニング用)
% =========================================================================

function [simData, trajPlan_Replanned] = run_simulation_with_replanning(mission, params, initial_target_pos, wind_planned, wind_actual, replan_time)
    % RUN_SIMULATION_WITH_REPLANNING (Clean Slate / SimSettings Rebuild)
    
    fprintf('\n=======================================================\n');
    fprintf('   Starting Simulation with Replanning Sequence\n');
    fprintf('=======================================================\n');

    %% 1. Phase 1: 初期計画に基づく飛行
    fprintf('[Phase 1] Flying with Initial Plan...\n');
    
    disturbance_1 = wind_actual - wind_planned;
    y0_phase1 = get_initial_state_from_planner(mission.Planner);
    
    plant1 = ParafoilDynamics(params);
    if isprop(plant1, 'WindVector'), plant1.WindVector = [disturbance_1; 0]; end
    
    autopilot1 = ParafoilAutopilot(mission.Params);
    autopilot1.import_mission_data(mission.Planner);
    autopilot1.set_longitudinal_control(true);
    autopilot1.set_nominal_glide_ratio(mission.PhysicsParams.glide_ratio);
    
    scheduler1 = PlannerAutopilotScheduler(mission, autopilot1);
    scheduler1.WindVector_Truth = [disturbance_1; 0]; 
    scheduler1.WindVector_Est   = [0; 0];             
    
    dt = 0.05;
    engine1 = SimulationEngine(plant1, scheduler1, dt, replan_time, y0_phase1);
    engine1.TargetNED = []; 
    
    engine1 = engine1.run();
    res1 = engine1.create_results_table();

    % ---------------------------------------------------------
    % 1. Phase 1 終了時の全状態ベクトルを取得
    % ---------------------------------------------------------
    % engine1.ResultsMatrix の最後の行が、リプラン時刻(t=replan_time)の状態です
    % y_end_raw = [u; v; w; p; q; r; phi; theta; psi; N_sim; E_sim; D_sim]
    y_end_raw = engine1.ResultsMatrix(end, :)';
    
    % ---------------------------------------------------------
    % 2. 状態変数の分解と整理
    % ---------------------------------------------------------
    % (A) 速度・角速度・姿勢角
    % これらは「機体固定座標」や「姿勢」なので、座標系の移動に関係なくそのまま使えます
    u_now     = y_end_raw(1);
    v_now     = y_end_raw(2);
    w_now     = y_end_raw(3);
    p_now     = y_end_raw(4);
    q_now     = y_end_raw(5);
    r_now     = y_end_raw(6);
    phi_now   = y_end_raw(7);
    theta_now = y_end_raw(8);
    psi_now   = y_end_raw(9); % 方位角 (Yaw)
    
    % 対気速度の大きさ (PlannerのV0設定用)
    V_tas_now = norm([u_now, v_now, w_now]);
    
    % (B) 位置座標の変換 (Sim座標 -> 対地座標)
    % Phase 1 は「計画風(wind_planned)」で流される座標系だったので、
    % その「流された分」を足して、現在の本当の対地位置を計算します。
    pos_sim_N = y_end_raw(10);
    pos_sim_E = y_end_raw(11);
    pos_sim_D = y_end_raw(12);
    
    % ドリフト量の計算
    drift_N = wind_planned(1) * replan_time;
    drift_E = wind_planned(2) * replan_time;
    
    % 再計画のスタート地点 (対地座標)
    start_N_gnd = pos_sim_N + drift_N;
    start_E_gnd = pos_sim_E + drift_E;
    start_Alt_gnd = -pos_sim_D; % NEDのDは下向き正なので、高度はマイナス
    
    fprintf('[Replanning Init] Pos: N=%.1f, E=%.1f, Alt=%.1f, Heading=%.1f deg\n', ...
        start_N_gnd, start_E_gnd, start_Alt_gnd, rad2deg(psi_now));
    
    % ---------------------------------------------------------
    % 3. 次のステップへの適用
    % ---------------------------------------------------------
    
    % [用途1] Planner (SimSettings) 用の初期条件
    % Excel設定を上書きするために使う値
    new_settings_struct = struct();
    new_settings_struct.X_initial = start_N_gnd;
    new_settings_struct.Y_initial = start_E_gnd;
    new_settings_struct.h_init    = start_Alt_gnd;
    new_settings_struct.psi_initial_deg = rad2deg(psi_now);
    
    % [用途2] Phase 2 シミュレーション用 (Dynamics) の初期状態ベクトル y0
    % 位置(N,E,D)だけ対地座標に差し替え、他はそのままコピーする
    y0_phase2 = y_end_raw;
    y0_phase2(10) = start_N_gnd;
    y0_phase2(11) = start_E_gnd;
    y0_phase2(12) = pos_sim_D; % Phase 2の高度基準をどうするかによりますが、
                               % 通常は対地高度(-start_Alt_gnd)を入れます。
                               % ただし、Phase 2も相対座標でやるなら0リセット等の考慮が必要。
                               % 今回の「SimSettings再構築版」なら対地座標でOKです。
    y0_phase2(12) = -start_Alt_gnd;
    
    %% 2. Estimation & Handoff (現在位置の取得)
    fprintf('[Handoff] Calculating Start State for Replanning...\n');
    
    % 風推定
    wind_estimated = estimate_wind_ideal(res1, wind_planned, disturbance_1);
    
    % 現在の状態 (Sim座標系)
    y_end_p1 = engine1.ResultsMatrix(end, :)';
    pos_sim_end = y_end_p1(10:12); 
    
    % 対地位置に変換 (ここが新しいスタート地点)
    pos_gnd_N = pos_sim_end(1) + wind_planned(1) * replan_time;
    pos_gnd_E = pos_sim_end(2) + wind_planned(2) * replan_time;
    pos_gnd_Alt = -pos_sim_end(3); % 高度
    
    V_tas_now = norm(y_end_p1(1:3));
    psi_now   = y_end_p1(9);
    
    fprintf('  New Start Pos : N=%.1f, E=%.1f, Alt=%.1f\n', pos_gnd_N, pos_gnd_E, pos_gnd_Alt);

    %% 3. ★★★ Replanning: SimSettingsの再構築 ★★★
    fprintf('[Replanning] Rebuilding SimSettings from scratch...\n');
    
    % (A) 元のExcelファイルからパラメータをロード (デフォルト値を取得)
    if isprop(mission, 'ExcelFileName')
        excelFile = mission.ExcelFileName;
    else
        excelFile = 'parafoil_parameters_ref.xlsx'; 
    end
    
    % 設定読み込み関数がある前提 (なければデフォルトで作る)
    if exist('load_params_from_excel', 'file')
        [~, new_sim_settings] = load_params_from_excel(excelFile);
    else
        % 万が一関数がない場合、既存のmissionからコピー
        if isprop(mission, 'SimSettings')
            new_sim_settings = mission.SimSettings;
        else
            error('Cannot load SimSettings. load_params_from_excel missing.');
        end
    end
    
    % (B) 設定構造体を「現在の状態」で上書きする (Excelを書き換えるのと同じ効果)
    new_sim_settings.X_initial = pos_gnd_N;
    new_sim_settings.Y_initial = pos_gnd_E;
    new_sim_settings.h_init    = pos_gnd_Alt;
    new_sim_settings.psi_initial_deg = rad2deg(psi_now);
    
    % (C) 新しいMissionオブジェクトを生成
    mission_replan = ParafoilMissionClothoid(excelFile);
    
    % (D) 書き換えた設定を注入 (これが "Excelを作り直す" に相当)
    if isprop(mission_replan, 'SimSettings')
        mission_replan.SimSettings = new_sim_settings;
    elseif isprop(mission_replan, 'Settings')
        mission_replan.Settings = new_sim_settings;
    end
    
    % 念のため PhysicsParams も更新 (ダブルチェック)
    pp = mission_replan.PhysicsParams;
    pp.start_pos = [pos_gnd_N, pos_gnd_E, pos_gnd_Alt, psi_now];
    if isprop(pp, 'V0'), pp.V0 = V_tas_now; end
    if isprop(pp, 'h_init'), pp.h_init = pos_gnd_Alt; end
    mission_replan.PhysicsParams = pp;
    
    % (E) RunUp無効化 & 再計算
    mission_replan.set_clothoid_options(0.5, 0.5, 0.0);
    
    try
        % これで Planner は「今注入された SimSettings」を見て計算するはず
        mission_replan.run_wind_simulation(initial_target_pos, 500, wind_estimated);
        fprintf('  -> Success.\n');
    catch ME
        warning('Replanning Failed: %s', message);
        simData = apply_wind_drift(res1, wind_planned); 
        trajPlan_Replanned = [];
        return; 
    end
    
    % (F) 軌道取得 & 時刻調整
    trajPlan_Replanned = mission_replan.export_detailed_trajectory('Ground');
    trajPlan_Replanned.Time = trajPlan_Replanned.Time + replan_time; 
    
    % デバッグ: 本当に更新されたかチェック
    dist_err = norm(trajPlan_Replanned.Position(1, 1:2) - [pos_gnd_N, pos_gnd_E]);
    if dist_err > 10
        warning('Start position still mismatches! Planner might be ignoring injected Settings.');
    else
        fprintf('  [Check] Start position matched correctly (Diff: %.2f m)\n', dist_err);
    end

    %% 4. Phase 2: 修正計画に基づく飛行
    fprintf('[Phase 2] Flying Updated Plan...\n');
    
    disturbance_2 = wind_actual - wind_estimated;
    
    % 初期条件: 位置はそのまま対地位置を使う
    y0_p2 = y_end_p1; 
    y0_p2(10) = pos_gnd_N;
    y0_p2(11) = pos_gnd_E;
    y0_p2(12) = -pos_gnd_Alt;
    
    plant2 = ParafoilDynamics(params);
    if isprop(plant2, 'WindVector'), plant2.WindVector = [disturbance_2; 0]; end
    
    % ★重要: 新しい mission_replan を使う
    autopilot2 = ParafoilAutopilot(mission_replan.Params); 
    autopilot2.import_mission_data(mission_replan.Planner); 
    
    autopilot2.set_longitudinal_control(true);
    autopilot2.set_nominal_glide_ratio(mission_replan.PhysicsParams.glide_ratio);
    
    scheduler2 = PlannerAutopilotScheduler(mission_replan, autopilot2);
    scheduler2.WindVector_Truth = [disturbance_2; 0];
    scheduler2.WindVector_Est   = [0; 0];
    
    % 終了条件
    trajPlan_Air_2 = mission_replan.export_detailed_trajectory('Air');
    t_max_2 = trajPlan_Air_2.Time(end) + 60;
    air_tgt = trajPlan_Air_2.Position(end, :);
    tgt_ned_2 = [air_tgt(1); air_tgt(2); -air_tgt(3)];
    
    % 実行
    engine2 = SimulationEngine(plant2, scheduler2, dt, t_max_2, y0_p2);
    engine2.TargetNED = tgt_ned_2; 
    
    engine2 = engine2.run();
    res2 = engine2.create_results_table(); 
    
    %% 5. 結合
    res1 = apply_wind_drift(res1, wind_planned);
    
    if ~isempty(res2)
        res2.Time = res2.Time + replan_time;
        res2 = apply_wind_drift(res2, wind_estimated);
        simData = [res1; res2(2:end, :)];
    else
        simData = res1;
    end
    
    fprintf('=== Simulation Complete ===\n');
end

% --- Helpers (変更なし) ---
function wind_est_total = estimate_wind_ideal(simData, wind_planned, disturbance_truth)
    idx = height(simData);
    u_g = simData.Body_U_Vel(idx);
    v_g = simData.Body_V_Vel(idx);
    w_g = simData.Body_W_Vel(idx);
    phi = deg2rad(simData.Roll_Angle(idx));
    theta = deg2rad(simData.Pitch_Angle(idx));
    psi = deg2rad(simData.Yaw_Angle(idx));
    R_ib = rotation_matrix(phi, theta, psi);
    V_gnd = R_ib * [u_g; v_g; w_g];
    V_air_true = V_gnd - [disturbance_truth; 0];
    pitot_val = norm(V_air_true);
    V_air_meas_body = [pitot_val; 0; 0];
    V_air_meas_inertial = R_ib * V_air_meas_body;
    w_dist_est = V_gnd(1:2) - V_air_meas_inertial(1:2);
    wind_est_total = wind_planned + w_dist_est;
end

function T = apply_wind_drift(T, w_vec)
    drift_N = w_vec(1) * T.Time;
    drift_E = w_vec(2) * T.Time;
    T.Inertial_X_Position = T.Inertial_X_Position + drift_N;
    T.Inertial_Y_Position = T.Inertial_Y_Position + drift_E;
    T.Inertial_Vx = T.Inertial_Vx + w_vec(1);
    T.Inertial_Vy = T.Inertial_Vy + w_vec(2);
    if ismember('Wind_X', T.Properties.VariableNames)
        T.Wind_X = T.Wind_X + w_vec(1);
        T.Wind_Y = T.Wind_Y + w_vec(2);
    end
end

function y0 = get_initial_state_from_planner(planner)
    d = planner.ResultDataNaive;
    x=d.runup.x(1); y=d.runup.y(1); z=d.runup.z(1); psi=d.runup.psi(1);
    if isprop(planner, 'V0'), V=planner.V0; else, V=15; end
    y0 = [V; 0; 0; 0; 0; 0; 0; 0; psi; x; y; -z];
end

function R = rotation_matrix(phi, theta, psi)
    cph=cos(phi); sph=sin(phi);
    cth=cos(theta); sth=sin(theta);
    cps=cos(psi); sps=sin(psi);
    R = [ ...
        cth*cps,   sph*sth*cps - cph*sps,   cph*sth*cps + sph*sps;
        cth*sps,   sph*sth*sps + cph*cps,   cph*sth*sps - sph*cps;
       -sth,       sph*cth,                 cph*cth ...
    ];
end
