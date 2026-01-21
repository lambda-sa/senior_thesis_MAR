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
                                     wind_planned, wind_planned);

% --- Case 2: 外乱あり (With Disturbance) ---
% Truth: 実際の風(外乱込み), Est: 計画風速(外乱を知らない)
wind_actual = wind_planned + wind_disturbance;
fprintf('Running Case 2: With Disturbance (Wind=[%.1f, %.1f])...\n', wind_actual(1), wind_actual(2));
simData_Wind = run_simulation_task(mission, params, trajPlanAir, target_pos, ...
                                   wind_actual, wind_planned);

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


%% ========================================================
%  共通シミュレーション実行関数 (提示コードのPhase B/Cを移植)
% =========================================================
function simData = run_simulation_task(mission, params, trajPlanAir, target_pos, wind_truth, wind_est)
    
    % 1. モデルとコントローラのインスタンス化
    % (paramsはメインスクリプトで整形済みなのでそのまま渡す)
    plant = ParafoilDynamics(params);
    
    autopilot = ParafoilAutopilot(mission.Params);
    autopilot.import_mission_data(mission.Planner);

    % 2. 計算された L/D を取得
    trim_ld = mission.PhysicsParams.glide_ratio;
    autopilot.set_longitudinal_control(true);

    % ★ここでPlannerの計算値を渡す
    autopilot.set_nominal_glide_ratio(trim_ld);
    % 2. スケジューラの作成
    scheduler = PlannerAutopilotScheduler(mission, autopilot);
    
    % ★重要: 風の設定 (Truth vs Est)
    scheduler.WindVector_Truth = [wind_truth; 0]; % 物理環境
    scheduler.WindVector_Est   = [wind_est; 0];   % 制御器の認識
    
    % 3. 初期条件の抽出 (Plannerの開始状態に合わせる)
    start_pos_ned = trajPlanAir.Position(1, :); 
    start_vel_tas = trajPlanAir.V_Air(1, :);   
    start_euler   = trajPlanAir.Euler_RPY(1, :); 
    
    % 初期状態ベクトルの作成
    y0 = [start_vel_tas(1); start_vel_tas(2); start_vel_tas(3); ... % uvw
          0; 0; 0; ...                                               % pqr
          start_euler(1); start_euler(2); start_euler(3); ...        % phi, theta, psi
          start_pos_ned(1); start_pos_ned(2); -start_pos_ned(3)];    % N, E, D (Alt->Down変換)
      
    % ※念のため trajPlanAir の第3成分が高度(正)なら Down(負)にする
    if start_pos_ned(3) > 0
        y0(12) = -start_pos_ned(3);
    end

    % 4. エンジンの初期化
    dt_6dof = 0.05;
    t_max_6dof = trajPlanAir.Time(end) + 10; % マージン
    engine = SimulationEngine(plant, scheduler, dt_6dof, t_max_6dof, y0);
    
    % ターゲット設定 (NED)
    target_pos_ned = [target_pos(1); target_pos(2); -target_pos(3)];
    engine.TargetNED = target_pos_ned;
    
    % 5. 実行
    engine = engine.run();
    
    % 6. 結果テーブル取得
    simData = engine.create_results_table();
end