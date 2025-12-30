% run_full_simulation.m
% 軌道計画 -> FF制御 -> 6自由度シミュレーション 統合実行スクリプト
clear; clc; close all;

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_SRC.xlsx'; % パラメータファイル
wind_vector_2d = [0; 0];  % 風速 (North, East) [m/s]
target_pos = [0, 2000, 3000]; % 目標 [N, E, Alt]
L_final = 2000;

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
mission = ParafoilMissionWind(excelFileName);

% シミュレーション実行 (内部で風補正計算が行われる)
mission.run_wind_simulation(target_pos, L_final, wind_vector_2d);

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
mapper = ParafoilControlMapper_linear(params);

% ★★★ 3-2. 参照データの抽出 ★★★
t_plan = trajPlanAir.Time;
phi_plan = trajPlanAir.Euler_RPY(:, 1);   % バンク角
theta_plan = trajPlanAir.Euler_RPY(:, 2); % ピッチ角 (ここ重要)

% 速度 (V_Air は [Vx, Vy, Vz] なのでノルムを計算)
V_plan = sqrt(sum(trajPlanAir.V_Air.^2, 2));

% ★★★ 3-3. スケジューラの作成 (VとThetaも渡す) ★★★
%scheduler = PlannerTrackingScheduler(mapper, t_plan, phi_plan, V_plan, theta_plan, wind_vector_2d);

scheduler = PlannerTrackingScheduler(mission,mapper);
% 3-3. 初期条件の抽出 (Plannerの開始状態に合わせる)
% Plannerの初期状態:
start_pos_ned = trajPlanGnd.Position(1, :); % [N, E, Alt] (Altは正)
start_vel_tas = trajPlanAir.V_Air(1, :);    % [Vx, Vy, Vz] (Wind Frame)
start_euler   = trajPlanAir.Euler_RPY(1, :); % [phi, theta, psi]

% 6DOF状態ベクトル: [u, v, w, p, q, r, phi, theta, psi, x, y, z]
% ※ 注意: 6DOFの z は "Down" (下向き正) なので、高度の符号反転が必要
% ※ 注意: u,v,w は機体座標系。start_vel_tas は対気慣性系に近い。
%   ここでは簡易的に、初期はトリム状態で u=V_tas, v=0, w=0 と仮定して近似するか、
%   あるいは回転行列で厳密に変換する。
%   ParafoilDynamics内部で整合性をとるため、ここでは近似値を与え、
%   最初の数ステップで物理的に落ち着くのを待つのが一般的です。

V_abs = norm(start_vel_tas);
u_init = V_abs; v_init = 0; w_init = 0; % Body frame approx
x_init = start_pos_ned(1);
y_init = start_pos_ned(2);
z_init = -start_pos_ned(3); % Height -> Down

y0 = [u_init; v_init; w_init; ... % Velocity (Body)
      0; 0; 0; ...                % Angular Rates
      start_euler(1); start_euler(2); start_euler(3); ... % Euler
      x_init; y_init; z_init];    % Position (NED)

% 3-4. エンジンの初期化
dt_6dof = 0.05;
t_max_6dof = t_plan(end) + 10; % 計画時間より少し長く
engine = SimulationEngine(plant, scheduler, dt_6dof, t_max_6dof, y0);

%% --- 4. Phase C: 実行 ---
fprintf('\n=== Phase C: Running 6DOF Simulation ===\n');
engine = engine.run();

% 結果テーブル取得
simData = engine.create_results_table();

%% --- 5. Phase D: 結果の比較・可視化 ---
fprintf('\n=== Phase D: Verification ===\n');

% レイアウト調整のためウィンドウサイズを縦に少し大きくしました
figure('Name', 'Planner vs 6DOF', 'Color', 'w', 'Position', [100, 50, 1200, 900]);

% 1. 3D軌跡比較 (左側の列すべてを使用: [1 3 5])
subplot(3, 2, [1 3 5]);
plot3(trajPlanGnd.Position(:,2), trajPlanGnd.Position(:,1), trajPlanGnd.Position(:,3), ...
      'g--', 'LineWidth', 1.5, 'DisplayName', 'Planner (Ideal)');
hold on;
plot3(simData.East, simData.North, -simData.Down, ...
      'b-', 'LineWidth', 2.0, 'DisplayName', '6DOF (Actual)');
plot3(target_pos(2), target_pos(1), target_pos(3), 'rx', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Target');
grid on; axis equal; view(3);
xlabel('East [m]'); ylabel('North [m]'); zlabel('Altitude [m]');
legend('Location', 'best');
title(sprintf('Trajectory Comparison\nWind: N=%.1f, E=%.1f m/s', wind_vector_2d(1), wind_vector_2d(2)));

% 2. バンク角追従確認 (右上: 2)
subplot(3, 2, 2);
plot(t_plan, rad2deg(phi_plan), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Command (Planner)');
hold on;
plot(simData.Time, simData.phi, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Response (6DOF)');
grid on; legend('Location','best');
xlabel('Time [s]'); ylabel('Bank Angle [deg]');
title('Bank Angle Tracking');
xlim([0, t_plan(end)]);

% 3. 【新規追加】高度 vs 時間 (右中: 4)
subplot(3, 2, 4);
% Plannerの高度: trajPlanGnd.Position(:, 3)
plot(t_plan, trajPlanGnd.Position(:, 3), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Planner)');
hold on;
% 6DOFの高度: -simData.Down (Downは下向き正なので反転)
plot(simData.Time, -simData.Down, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Act (6DOF)');
grid on; legend('Location','best');
xlabel('Time [s]'); ylabel('Altitude [m]');
title('Altitude History');
xlim([0, t_plan(end)]);

% 4. 制御入力履歴 (右下: 6)
subplot(3, 2, 6);
plot(simData.Time, simData.delta_R, 'r-', 'LineWidth', 1, 'DisplayName', '\delta_R');
hold on;
plot(simData.Time, simData.delta_L, 'b-', 'LineWidth', 1, 'DisplayName', '\delta_L');
plot(simData.Time, simData.delta_a, 'k--', 'LineWidth', 1, 'DisplayName', '\delta_a (Net)');
grid on; legend('Location','best');
xlabel('Time [s]'); ylabel('Control Input');
title('Control Inputs Generated by Mapper');
xlim([0, t_plan(end)]);

% 最終誤差の表示
final_pos_6dof = [simData.North(end), simData.East(end)];
err_dist = norm(final_pos_6dof - target_pos(1:2));
fprintf('Final Landing Error (6DOF): %.2f m\n', err_dist);

%% --- Figure 2: Velocity Comparison (NED Ground Frame) ---
fprintf('Plotting Ground Velocity Comparison...\n');

figure('Name', 'Ground Velocity Comparison (NED)', 'Color', 'w', 'Position', [150, 100, 800, 800]);

% --- 1. 実際の軌道 (6DOF) の対地速度計算 ---
% 機体座標系速度 [u, v, w] を NED座標系 [Vn, Ve, Vd] に変換
V_ned_6dof = zeros(height(simData), 3);

for i = 1:height(simData)
    % オイラー角 (rad)
    phi   = simData.phi(i);
    theta = simData.theta(i);
    psi   = simData.psi(i);
    
    % 回転行列 (Body -> NED)
    % c: cos, s: sin
    c_th = cos(theta); s_th = sin(theta);
    c_ps = cos(psi);   s_ps = sin(psi);
    c_ph = cos(phi);   s_ph = sin(phi);
    
    R_b2n = [ ...
        c_th*c_ps,    s_ph*s_th*c_ps - c_ph*s_ps,    c_ph*s_th*c_ps + s_ph*s_ps;
        c_th*s_ps,    s_ph*s_th*s_ps + c_ph*c_ps,    c_ph*s_th*s_ps - s_ph*c_ps;
       -s_th,         s_ph*c_th,                     c_ph*c_th
    ];
    
    % 機体固定座標系の速度ベクトル
    V_body = [simData.u(i); simData.v(i); simData.w(i)];
    
    % 座標変換
    V_ned_6dof(i, :) = (R_b2n * V_body)';
end

% --- 2. 参照軌道 (Planner) の対地速度計算 ---
% 対気速度 [V_air] を変換し、風ベクトルを加算して対地速度にする
V_ned_plan = zeros(length(t_plan), 3);
wind_3d = [wind_vector_2d; 0]; % 風ベクトル [North, East, Down] (Downは通常0)

for i = 1:length(t_plan)
    % Plannerの姿勢角
    phi_p   = trajPlanAir.Euler_RPY(i, 1);
    theta_p = trajPlanAir.Euler_RPY(i, 2);
    psi_p   = trajPlanAir.Euler_RPY(i, 3);
    
    % 回転行列 (Body -> NED)
    c_th = cos(theta_p); s_th = sin(theta_p);
    c_ps = cos(psi_p);   s_ps = sin(psi_p);
    c_ph = cos(phi_p);   s_ph = sin(phi_p);
    
    R_b2n_p = [ ...
        c_th*c_ps,    s_ph*s_th*c_ps - c_ph*s_ps,    c_ph*s_th*c_ps + s_ph*s_ps;
        c_th*s_ps,    s_ph*s_th*s_ps + c_ph*c_ps,    c_ph*s_th*s_ps - s_ph*c_ps;
       -s_th,         s_ph*c_th,                     c_ph*c_th
    ];
    
    % Plannerの対気速度 (Body Frame)
    V_air_body = trajPlanAir.V_Air(i, :)'; 
    
    % 対地速度 = (回転行列 * 対気速度) + 風ベクトル
    V_ned_plan(i, :) = (R_b2n_p * V_air_body + wind_3d)'; 
end

% --- 3. プロット描画 ---

% (1) North Velocity (Vx)
subplot(3, 1, 1);
plot(t_plan, V_ned_plan(:, 1), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Planner)');
hold on;
plot(simData.Time, V_ned_6dof(:, 1), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Act (6DOF)');
grid on; legend('Location', 'best');
ylabel('V_{North} (x) [m/s]');
title('Ground Velocity Comparison (NED Frame)');
xlim([0, t_plan(end)]);

% (2) East Velocity (Vy)
subplot(3, 1, 2);
plot(t_plan, V_ned_plan(:, 2), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Planner)');
hold on;
plot(simData.Time, V_ned_6dof(:, 2), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Act (6DOF)');
grid on; legend('Location', 'best');
ylabel('V_{East} (y) [m/s]');
xlim([0, t_plan(end)]);

% (3) Down Velocity (Vz) - 降下速度
subplot(3, 1, 3);
plot(t_plan, V_ned_plan(:, 3), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Planner)');
hold on;
plot(simData.Time, V_ned_6dof(:, 3), 'b-', 'LineWidth', 1.5, 'DisplayName', 'Act (6DOF)');
grid on; legend('Location', 'best');
ylabel('V_{Down} (z) [m/s]');
xlabel('Time [s]');
xlim([0, t_plan(end)]);

fprintf('  -> Ground velocity plot created.\n');