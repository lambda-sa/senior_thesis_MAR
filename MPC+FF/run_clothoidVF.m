% run_clothoidVF.m
% 軌道計画 -> FF制御 -> 6自由度シミュレーション 統合実行スクリプト
clear; clc; close all;

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_ref.xlsx'; % パラメータファイル
wind_planned = [0; 0];  % 風速 (North, East) [m/s]
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

%% --- 3. Phase B: 6DOFシミュレーション準備 ---
fprintf('\n=== Phase B: 6DOF Simulation Setup ===\n');
%% --- 2. モデル & 初期条件の準備 ---
fprintf('=== Initializing Models ===\n');



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

% (A) 6-DOF Dynamics Model

%% --- 3. Phase B: 6DOFシミュレーション準備 ---
fprintf('\n=== Phase B: 6DOF Simulation Setup (Pure Feedforward) ===\n');

% 3-1. モデルとコントローラのインスタンス化
plant = ParafoilDynamics(params);
%mapper = ParafoilControlMapper_linear_old(params);




% (1) オートパイロットの初期化
autopilot = ParafoilAutopilot(mission.Params);

% (2) 計画データのインポート
% Plannerが作った軌道(Entry->Dubins->Final)とLoiter情報をAutopilotに渡す
autopilot.import_mission_data(mission.Planner);

% (3) ゲインの調整 (必要に応じて)
%autopilot.Gains.k_vf_loiter = 0.04;  % 円旋回は滑らかに
%autopilot.Gains.k_vf_mission = 0.12; % 着陸進入は強めに誘導

% (4) スケジューラ(Adapter)の作成
% ダイナミクスとオートパイロットの仲介役
scheduler = PlannerAutopilotScheduler(mission, autopilot);

% ★★★ ここを追加: Schedulerが持つ風ベクトルを「実際の風」で上書きする ★★★
% SchedulerはデフォルトでPlannerの風(wind_planned)を読み込んでしまうため、
% ここで強制的に wind_actual (NED 3次元) に書き換えます。
% =========================================================
% ★★★ Schedulerへの風設定 (クラブ角補正を有効化) ★★★
% =========================================================

% 1. 物理モデル用 (Truth): 「実際の風 (Actual)」
% (3) 風情報の更新 (Plannerの予測風ではなく、Simulation用の"正解"風を与える)
% Scheduler作成時にデフォルト値が入っていますが、ここでシミュレーション条件(wind_actual)で上書きします
scheduler.WindVector_Truth = [wind_actual; 0]; % 物理環境の風
scheduler.WindVector_Est   = [wind_actual; 0]; % 制御器が知っている風

fprintf('Setting: Scheduler configured with Dubins-Loiter & Actual Wind.\n');
fprintf('Simulation Setup: Wind Truth=[%.1f, %.1f], Est=[%.1f, %.1f]\n', ...
    wind_actual(1), wind_actual(2), wind_actual(1), wind_actual(2));

% --- ★削除した部分★ ---
% MissionSwitchAlt の手動計算ロジックは全て削除しました。
% Scheduler が自動計算した値が使われます。
fprintf('Mission Switch Altitude is automatically set to: %.1f m (NED)\n', ...
    scheduler.MissionSwitchAlt);


% 3-3. 初期条件の抽出 (Plannerの開始状態に合わせる)
start_pos_ned = trajPlanGnd.Position(1, :); 
start_vel_tas = trajPlanAir.V_Air(1, :);   
start_euler   = trajPlanAir.Euler_RPY(1, :); 

V_abs = norm(start_vel_tas);
u_init = V_abs; v_init = 0; w_init = 0; 
x_init = start_pos_ned(1);
y_init = start_pos_ned(2);
z_init = -start_pos_ned(3); % Height -> Down

y0 = [u_init; v_init; w_init; ... 
      0; 0; 0; ...                
      start_euler(1); start_euler(2); start_euler(3); ... 
      x_init; y_init; z_init];
% 3-4. エンジンの初期化
dt_6dof = 0.05;
t_max_6dof = t_plan(end); % 計画時間より少し長く
engine = SimulationEngine(plant, scheduler, dt_6dof, t_max_6dof, y0);

%% --- 4. Phase C: 実行 ---
fprintf('\n=== Phase C: Running 6DOF Simulation ===\n');
% =========================================================
% ★追加: 会合点(Target)への接近停止条件
% =========================================================
% target_pos は [North, East, Altitude] なので NED (Down) に変換
target_pos_ned = [target_pos(1); target_pos(2); -target_pos(3)];

engine.TargetNED = target_pos_ned;  % 目標座標 [N; E; D]
%engine.StopThreshold = 50.0;       % 停止距離 [m]

%fprintf('Termination Condition: Stop if distance to target < %.1f m\n', engine.StopThreshold);

% 実行
engine = engine.run();

% 結果テーブル取得
simData = engine.create_results_table();
% ゴール地点（最接近点）の座標を取得（プロット用）
final_pos_ned = [simData.Inertial_X_Position(end), simData.Inertial_Y_Position(end), simData.Inertial_Z_Position(end)];
final_pos_alt = [final_pos_ned(1), final_pos_ned(2), -final_pos_ned(3)];

%% --- 5. Phase D: 結果の比較・可視化 ---
fprintf('\n=== Phase D: Verification ===\n');
simData = engine.create_results_table();
PayloadPlotter.plotResults(simData, engine.ControlScheduler); % [cite: 2848]


% =========================================================
% ★設定: フェーズ区切りの縦線を表示するかどうかのフラグ
% =========================================================
SHOW_PHASE_LINES = false;  % true: 表示する, false: 表示しない

% 1. Plannerからフェーズ遷移時刻を抽出 (ログ表示用に計算は常に行う)
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
                plan_trans_times = [plan_trans_times, internal_times(:)'];
            end
        end
    end
    
    plan_trans_times = unique(plan_trans_times);
end

fprintf('Planned Transition Times: %s s\n', mat2str(plan_trans_times, 2));

% 2. 比較クラスへの引数をフラグで切り替え
if SHOW_PHASE_LINES
    times_arg = plan_trans_times; % 時刻リストを渡す
else
    times_arg = [];               % 空配列を渡す（線は引かれない）
end
% 2. 比較クラスのインスタンス化
comparator = TrajectoryComparator(simData, trajPlanGnd, trajPlanAir, target_pos, wind_planned, plan_trans_times);

% 3. 描画
comparator.plotAll();

fprintf('Done.\n')