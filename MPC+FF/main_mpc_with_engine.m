%% =========================================================================
%% FILE: main_mpc_with_engine.m
%% =========================================================================
% SimulationEngine + MPCScheduler を使用した統合シミュレーション

clearvars; close all; clc;

%% --- 1. 初期化 & パラメータ設定 ---
fprintf('=== 1. Initialization ===\n');
excelFileName = 'parafoil_parameters_ref.xlsx'; 
[params, sim_settings] = load_params_from_excel(excelFileName);

% 6-DOF用パラメータ拡張 (main_6dof_adapter.m と同様)
params.I_total_body = [params.I_xx, 0, params.I_xz; 0, params.I_yy, 0; params.I_xz, 0, params.I_zz];
prop = struct();
prop.TotalMass = params.m_total;
prop.b = params.b; prop.c = params.c; prop.d = params.d;
prop.r_total_cm_to_canopy_origin_B = [params.Rcg_x; 0; params.Rcg_z];
prop.r_total_cm_to_payload_cm_B = [params.Rpg_x; 0; params.Rpg_z];
prop.r_canopy_origin_to_ac_B = [0; 0; 0];
params.prop = prop;

% 6-DOF プラントモデル
real_plant = ParafoilDynamics(params);


%% --- 2. 参照軌道生成 (Guidance) ---
fprintf('=== 2. Guidance Generation ===\n');
target_pos = [0, 0, 0];
L_final = 100;
wind_vector = [5.0; -3.0]; % [Nx, Ey]

mission = ParafoilMissionWind(excelFileName);
mission.set_turn_input_ratio(0.8);
mission.run_wind_simulation(target_pos, L_final, wind_vector);

tabGround = mission.export_detailed_trajectory('Ground');
T_ref = tabGround.Time;
% 参照軌道行列の作成 [12 x N]
X_ref_all = zeros(12, length(T_ref));
for k = 1:length(T_ref)
    pos = tabGround.Position(k,:); eul = tabGround.Euler_RPY(k,:); vel = tabGround.V_Air(k,:);
    % NED変換: z = -Alt
    X_ref_all(:,k) = [vel(1); vel(2); vel(3); 0; 0; 0; eul(1); eul(2); eul(3); pos(1); pos(2); -pos(3)];
end


%% --- 3. MPCコントローラ & スケジューラ準備 ---
fprintf('=== 3. Controller Setup ===\n');
% 線形化
linearizer = ParafoilLinearizer(real_plant);
h_trim = 1000; V_trim = 15;
[A_lin, B_lin] = linearizer.get_linear_model([V_trim;0;0;0;0;0;0;0;0;0;0;-h_trim], 0);

% MPC設定
dt_control = 0.1;
mpc_ctl = ParafoilMPC(A_lin, B_lin, dt_control, 20); % Horizon=20
mpc_ctl.Q = diag([0 0 0, 0 0 0, 10 0 100, 50 50 0]); % 重み
mpc_ctl.R = 50;

% ★★★ ここが変更点: Scheduler の作成 ★★★
% エンジンに渡すためのスケジューラインスタンスを作成
wind_3d = [wind_vector(1); wind_vector(2); 0];
scheduler = MPCScheduler(mpc_ctl, X_ref_all, T_ref, wind_3d);


%% --- 4. SimulationEngine による実行 ---
fprintf('=== 4. Running Simulation Engine ===\n');

% 初期条件 (参照軌道始点 + 位置誤差)
y0 = X_ref_all(:, 1);
y0(11) = y0(11) + 20; % 東に20mずらす

sim_time = T_ref(end) + 10;

% ★★★ エンジンの生成と実行 ★★★
engine = SimulationEngine(real_plant, scheduler, dt_control, sim_time, y0);
engine = engine.run();

% 結果テーブルの取得
results = engine.create_results_table();


%% --- 5. 結果の可視化 ---
fprintf('=== 5. Visualization ===\n');
figure('Name', 'Engine MPC Result', 'Color', 'w', 'Position', [100, 100, 1000, 600]);

% 軌跡
subplot(2,2,1); hold on; grid on; axis equal;
plot(X_ref_all(11,:), X_ref_all(10,:), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref');
plot(results.East, results.North, 'b-', 'LineWidth', 2, 'DisplayName', 'MPC');
plot(0, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('East'); ylabel('North'); legend; title('Trajectory');

% 入力 (Engineが記録した delta_a を使用)
subplot(2,2,2); hold on; grid on;
plot(results.Time, results.delta_a, 'r-', 'LineWidth', 1.5);
yline([1, -1], 'k:');
xlabel('Time [s]'); ylabel('\delta_a'); title('Control Input');

% 高度
subplot(2,2,3); hold on; grid on;
plot(results.Time, -results.Down, 'b-');
xlabel('Time [s]'); ylabel('Alt [m]');

% 方位角
subplot(2,2,4); hold on; grid on;
plot(results.Time, rad2deg(X_ref_all(9, 1:length(results.Time))), 'g--'); % 簡易比較
plot(results.Time, results.psi, 'b-');
xlabel('Time [s]'); ylabel('Heading [deg]');

fprintf('Final Error: %.2f m\n', norm([results.North(end), results.East(end)]));