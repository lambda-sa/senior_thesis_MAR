%% 
% main_reproduce_paper.m
% 論文 [IAC-21-D2-3-2-x66968] Figure 19, 20 再現用スクリプト
% ベース: ユーザー提供の main.m 構造

clear; clc; close all;

%% --- 1. 設定 & パラメータ読み込み ---
fprintf('=== Phase 1: Initialization ===\n');

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_ref.xlsx'; % パラメータファイル
wind_vector_2d = [0; 0];  % 風速 (North, East) [m/s]
target_pos = [0, 600, 1000]; % 目標 [N, E, Alt]
L_final = 500;

% パラメータ読み込みとモデル構築
if exist('load_params_from_excel', 'file')
    [params, sim_settings] = load_params_from_excel(excelFileName);
else
    warning('load_params_from_excelが見つかりません。デフォルト値を使用します。');
    params = ParafoilParams(); 
    sim_settings = struct('X_initial', 0, 'Y_initial', 50, 'h_init', 1000, 't_max', 60, 'psi_initial_deg', 0);
end 

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
% 1-2. シミュレーション設定
%t_max = 400; % 論文の着陸までの時間に合わせて設定
%dt = 0.05;   % タイムステップ

% 1-3. 物理プロパティの整合性確保 (ParafoilDynamics用)
% ParamParafoilPaper内で既に prop 構造体は作成されていますが、
% ParafoilDynamics が要求する形であることを確認
if ~isfield(params, 'prop')
    error('params に prop 構造体が含まれていません。ParamParafoilPaperを確認してください。');
end


%% --- 2. Phase A: 環境と軌道のセットアップ (Mission Setup) ---
fprintf('\n=== Phase A: Environment & Guidance Setup ===\n');

% シミュレーション時間とステップ
t_max = 600; 
dt = 0.01; % 制御周期に合わせて少し細かく設定

%% 2. 初期条件の設定
% 位置: 高度 1000m
pos_0 = [0; 0; -3300];

% 速度: 対気速度 13m/s (安定飛行速度)
u_init = 13; v_init = 0; w_init = 0;
vel_B_0 = [u_init; v_init; w_init];

% 姿勢・角速度: 初期値ゼロ
att_0 = [0; 0; 0];
omega_0 = [0; 0; 0];

% 状態ベクトルの結合
y0 = [vel_B_0; omega_0; att_0; pos_0];

%% 3. 環境設定 (一定風のみ)
% 例: 東向き 2m/s の一定風 (乱気流なし)
mean_wind = [0; 0; 0]; 

%% --- 4. スケジューラとエンジンの準備 ---
% ★重要: ここで新しい物理ベーススケジューラを作成
fprintf('Using Physics-Based Guidance Scheduler...\n');

% 物理パラメータを抽出してスケジューラに渡す
phys_params.b    = params.b;
phys_params.Cnr  = params.C_n_r;   % ※ paramsから取得できるなら params.Cnr
phys_params.Cnda = params.C_n_delta_a;   % ※ paramsから取得できるなら params.Cnda

scheduler =VFPathFollowingScheduler_Cascade(dt, phys_params);

% SimulationEngineが使うプロパティを設定
% (クラス側でこれらのプロパティを定義したのでエラーになりません)
scheduler.Waypoints = [0, 0; 800, 800];   
scheduler.OrbitCenter = [600, 600];       
scheduler.OrbitRadius = 300;              
scheduler.OrbitDirection = 1;             
scheduler.d_boundary = 100; % この距離に入ると自動でOrbitモードになります

% ダイナミクス設定
plant = ParafoilDynamics(params);
plant.use_soft_min = false; 

% シミュレーションエンジン作成 (既存のまま)
engine = SimulationEngine(plant, scheduler, dt, t_max, y0);

%% 5. 実行
fprintf('=== Phase 2: Running Simulation ===\n');
try
    % engineは内部で scheduler.do_control(t, y, mean_wind) を呼び出します
    engine = engine.run(); 
    simData = engine.create_results_table();
    fprintf('  -> Simulation Complete.\n');
catch ME
    fprintf('  -> Simulation Failed: %s\n', ME.message);
    % エラー箇所の詳細表示
    for k = 1:length(ME.stack)
        fprintf('     File: %s, Line: %d, Name: %s\n', ...
            ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
    end
    return;
end

%% 6. 可視化
fprintf('=== Phase 3: Visualization ===\n');

% --- XY軌跡 ---
figure('Name', 'Trajectory XY (Physics-Based)', 'Color', 'w');
plot(simData.Inertial_Y_Position, simData.Inertial_X_Position, 'b-', 'LineWidth', 1.5); hold on;
line([0, 800], [0, 800], 'Color', 'r', 'LineStyle', '--');
viscircles(scheduler.OrbitCenter([2,1]), scheduler.OrbitRadius, 'Color', 'r', 'LineStyle', ':');
axis equal; grid on;
xlabel('East [m]'); ylabel('North [m]');
title('Trajectory XY (Physics-Based Cascade)');
legend('Flight Path', 'Target Path');

% --- 制御応答 ---
figure('Name', 'Control Response', 'Color', 'w');
subplot(3,1,1);
plot(simData.Time, rad2deg(simData.Yaw_Angle), 'b');
ylabel('Yaw [deg]'); grid on; title('Yaw Angle (Target: \psi_{ref})');

subplot(3,1,2);
plot(simData.Time, rad2deg(simData.Body_R_Vel), 'm');
ylabel('Rate [deg/s]'); grid on; title('Yaw Rate (Target: r_{cmd})');

subplot(3,1,3);
plot(simData.Time, simData.delta_a, 'r');
ylabel('\delta_a'); grid on; title('Control Input');
xlabel('Time [s]');


engine = SimulationEngine(plant, scheduler, dt, t_max, y0);
simData = engine.create_results_table();
PayloadPlotter.plotResults(simData, engine.ControlScheduler); % [cite: 2848]