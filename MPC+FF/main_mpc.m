%% =========================================================================
%% main_mpc_with_engine.m (Corrected Version)
%% SimulationEngine + MPC Scheduler による実行スクリプト
%% =========================================================================
clear; clc; close all;

% ★修正1: MPCの不要な警告("...が空です")を消す
mpcverbosity('off'); 

%% --- 1. 経路生成 (Reference Generation) ---
fprintf('=== Path Planning ===\n');
excelFileName = 'parafoil_parameters_SRC.xlsx'; 
target_pos = [0, 2000, 3000];
L_final = 1000;
wind_vector_2d = [3.0; 3.0]; % 北東の風

mission = ParafoilMissionWind(excelFileName);
mission.run_wind_simulation(target_pos, L_final, wind_vector_2d);

% MPC用参照軌道
ref_table = mission.export_detailed_trajectory('Ground');


%% --- 2. モデル & 初期条件の準備 ---
fprintf('=== Initializing Models ===\n');

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

% (A) 6-DOF Dynamics Model
plant = ParafoilDynamics(params); 

% (B) 初期条件の設定 (Missionから高精度な初期値を取得)
%     export_6dof_state() を使って、位置・姿勢・速度・角速度をすべて取得します
initState = mission.export_6dof_state(); 
% 最初の行 (t=0) を取得
initRow = initState(1, :);

% 状態ベクトル y0 (12x1) の構築
% y = [u, v, w, p, q, r, phi, theta, psi, N, E, D]
y0 = zeros(12, 1);
y0(1) = initRow.u;
y0(2) = initRow.v;
y0(3) = initRow.w;
y0(4) = initRow.p;   % ★重要: 旋回開始ならp,q,rはゼロではない
y0(5) = initRow.q;
y0(6) = initRow.r;
y0(7) = initRow.phi;
y0(8) = initRow.theta;
y0(9) = initRow.psi;
y0(10)= initRow.N;
y0(11)= initRow.E;
y0(12)= initRow.D;   % Mission側ですでに「高度の反転(Down)」になっているのでそのまま代入


%% --- 3. MPC Scheduler の構築 (★ここが修正の核心) ---
fprintf('=== Configuring MPC Scheduler ===\n');
Ts_mpc = 0.1;    % MPC制御周期 [s]
Horizon = 20;    % 予測ホライゾン (長すぎると重くなるので20~50推奨)

% ★修正3: 新しい統合クラス「ParafoilMPC_Integrated」を使用する
% (前回作成したクラスファイル名が ParafoilMPC_Integrated.m である前提)
mpc_ctrl = ParafoilMPC_6DOF(Ts_mpc, Horizon, params);

% ★修正4: 参照軌道を事前にセットする (これが抜けていたためエラーになった)
% 位置データは [North, East, Down] に変換して渡す
%ref_time = ref_table.Time;
%ref_pos_ned = [ref_table.Position(:,1), ref_table.Position(:,2), -ref_table.Position(:,3)];
%alpha_trim = 0.1; % トリム迎え角 [rad] (適当な初期値)

%mpc_ctrl.set_reference_trajectory(ref_time, ref_pos_ned, alpha_trim);

% ★修正: Missionオブジェクトを丸ごと渡す
% これにより、内部で「参照軌道データのロード」と「初期線形化(A, B行列計算)」が自動実行されます
mpc_ctrl.load_reference(mission);

% スケジューラの作成
wind_3d = [wind_vector_2d; 0];
gamma_rigging = 0; % リギング角
scheduler = MPCScheduler(mpc_ctrl, ref_table, wind_3d, gamma_rigging);


%% --- 4. SimulationEngine のセットアップと実行 ---
fprintf('=== Running Simulation Engine ===\n');
dt_sim = 0.01; % 物理計算ステップ [s]
T_max = ref_table.Time(end) + 10.0;

% エンジン生成
engine = SimulationEngine(plant, scheduler, dt_sim, T_max, y0);

% 実行 (try-catchでエラー内容を見やすく表示)
try
    engine = engine.run();
catch ME
    fprintf(2, 'エラーが発生しました:\n%s\n', ME.message);
    % スタックトレースの表示
    for k = 1:length(ME.stack)
        fprintf('  File: %s, Line: %d, Name: %s\n', ...
            ME.stack(k).file, ME.stack(k).line, ME.stack(k).name);
    end
    return;
end


%% --- 5. 結果の確認 ---
fprintf('=== Verification & Plotting ===\n');

% 結果テーブルの取得
res = engine.create_results_table();

% シミュレーションが即死していないかチェック
if height(res) < 5
    warning('シミュレーションがすぐに終了しました。エラーログを確認してください。');
else
    % 軌跡プロット
    figure('Name', 'Trajectory Result', 'Color', 'w');
    plot(res.East, res.North, 'b-', 'LineWidth', 1.5); hold on;
    
    % 参照軌道のプロット
    plot(ref_table.Position(:,2), ref_table.Position(:,1), 'g--', 'LineWidth', 1.5);
    
    xlabel('East [m]'); ylabel('North [m]');
    legend('MPC Controlled', 'Reference Path');
    grid on; axis equal;
    title('SimulationEngine with MPC');
    
    % 入力とヨー角の確認
    figure('Name', 'Control Inputs', 'Color', 'w');
    subplot(2,1,1);
    plot(res.Time, res.delta_L, 'r', 'DisplayName', 'Left'); hold on;
    plot(res.Time, res.delta_R, 'b', 'DisplayName', 'Right');
    ylabel('Brake Input'); legend; grid on;
    
    subplot(2,1,2);
    plot(res.Time, rad2deg(res.psi), 'b'); hold on;
    % 参照ヨー角の補間 (表示用)
    ref_psi = interp1(ref_table.Time, ref_table.Euler_RPY(:,3), res.Time, 'linear', 'extrap');
    plot(res.Time, rad2deg(ref_psi), 'g--');
    ylabel('Yaw [deg]'); legend('Actual', 'Target'); grid on;
    
    % 高度チェックのプロット
    figure('Name', 'Altitude Check', 'Color', 'w');
    plot(res.Time, -res.Down, 'b'); % Downを反転して高度表示
    hold on;
    plot(ref_table.Time, ref_table.Position(:,3), 'g--');
    xlabel('Time [s]'); ylabel('Altitude [m]');
    legend('Actual', 'Reference');
    grid on; title('Flight Altitude');
end