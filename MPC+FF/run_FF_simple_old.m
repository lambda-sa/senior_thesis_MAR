% run_simple_test.m
% 単純軌道（直線/旋回）を用いた6DOFシミュレーションテスト
clear; clc; close all;

%% --- 1. 設定 ---
excelFileName = 'parafoil_parameters_SRC.xlsx';
wind_vector_2d = [0; 0];  % テスト用に風は無風推奨
% wind_vector_2d = [3; -3]; % 必要なら風を入れる

% パラメータ読み込みとモデル構築
if exist('load_params_from_excel', 'file')
    [params, sim_settings] = load_params_from_excel(excelFileName);
else
    warning('load_params_from_excelが見つかりません。デフォルト値を使用します。');
    params = ParafoilParams(); 
    sim_settings = struct('X_initial', 0, 'Y_initial', 50, 'h_init', 1000, 't_max', 60, 'psi_initial_deg', 0);
end
% テストモード設定
% シミュレーション条件
test_type = 'Turn';      % 'Straight' or 'Turn'
bank_cmd  = 7.0;        % 旋回時のバンク角 [deg]
sim_duration = 600;       % 秒
wind_vec = [0; 0];       % 検証用は無風推奨（風を入れても計算はされます）

%% 2. ミッション実行 (Physics-Based)
fprintf('=== Initializing Mission with Physics ===\n');
mission = ParafoilMissionWind(excelFileName);

% 必要なら風を設定 (Planner内部にセットされる)
mission.Planner.WindVector = [wind_vec; 0];

% ★ここが核心：Excelパラメータからトリムを計算し、軌道生成
mission.run_simple_simulation(test_type, sim_duration, bank_cmd);

% 参照軌道の取得
trajRef = mission.export_detailed_trajectory('Air');

%% --- 3. Phase B & C: 6DOFシミュレーション ---
fprintf('\n=== Phase B: 6DOF Simulation Setup ===\n');

% モデル準備
%params = mission.Params;
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
plant = ParafoilDynamics(params);
mapper = ParafoilControlMapper_linear(params);

% トリム計算結果を初期値に使う（これがやりたかったこと）
% Missionが計算した物理パラメータを取得
pp = mission.PhysicsParams;

% 初期条件ベクトル (NED)
% [u, v, w, p, q, r, phi, theta, psi, N, E, D]
% ※ theta_trim は trajectoryの最初の値を使うのが最も整合性が取れる
phi_init   = trajRef.Euler_RPY(1, 1);
theta_init = trajRef.Euler_RPY(1, 2);
psi_init   = trajRef.Euler_RPY(1, 3);
V_tas_init = norm(trajRef.V_Air(1, :));

y0 = [V_tas_init; 0; 0; ...       % u, v, w (Body)
      0; 0; 0; ...                % p, q, r
      phi_init; theta_init; psi_init; ... % Euler
      trajRef.Position(1,1); trajRef.Position(1,2); -trajRef.Position(1,3)]; % NED Pos

% スケジューラ作成
scheduler = PlannerTrackingScheduler(mission,mapper);

% 実行
engine = SimulationEngine(plant, scheduler, 0.05, sim_duration, y0);
engine = engine.run();
simData = engine.create_results_table();

V_down_act = zeros(height(simData), 1);

for i = 1:height(simData)
    u = simData.u(i);
    v = simData.v(i);
    w = simData.w(i);
    phi = simData.phi(i);
    theta = simData.theta(i);
    
    % 回転行列の計算 (Body -> NED の Z成分のみ)
    % V_down = (RotMatrix * [u;v;w]) の 3要素目
    V_down_act(i) = -u * sin(theta) + ...
                     v * sin(phi) * cos(theta) + ...
                     w * cos(phi) * cos(theta);
end

%% 4. 比較プロット
figure('Name', 'Physics-Based Simple Path Verification', 'Color', 'w', 'Position', [100,100,1000,600]);

subplot(2,2,[1 3]);
plot3(trajRef.Position(:,2), trajRef.Position(:,1), trajRef.Position(:,3), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Trim)');
hold on;
plot3(simData.East, simData.North, -simData.Down, 'b-', 'LineWidth', 2, 'DisplayName', '6DOF');
grid on; axis equal; view(3);
xlabel('East'); ylabel('North'); zlabel('Alt');
legend; title(['Trajectory: ' test_type]);

subplot(2,2,2);
plot(trajRef.Time, rad2deg(trajRef.Euler_RPY(:,1)), 'g--', 'LineWidth',1.5, 'DisplayName','Ref');
hold on;
plot(simData.Time, simData.phi, 'b-', 'LineWidth',1.5, 'DisplayName','6DOF');
ylabel('Bank [deg]'); grid on; title('Bank Angle');

subplot(2,2,4);
plot(trajRef.Time, -trajRef.V_Air(:,3), 'g--', 'LineWidth',1.5, 'DisplayName','Ref');
hold on;
plot(simData.Time, simData.w, 'b-', 'LineWidth',1.5, 'DisplayName','6DOF w');
ylabel('Sink Speed [m/s]'); grid on; title('Descent Rate (Check Trim)');
legend;

% 2. プロット作成
figure('Name', 'Descent Rate Verification', 'Color', 'w');

% 参照軌道の降下速度 (V_AirのZ成分 = V_down)
% ※ trajRef.V_Air は [North, East, Down] のベクトルなので、3列目が降下速度
plot(trajRef.Time, trajRef.V_Air(:, 3), 'g--', 'LineWidth', 2.0, 'DisplayName', 'Ref (Planner)');

hold on;

% 6DOFの降下速度 (計算した V_down_act)
plot(simData.Time, V_down_act, 'b-', 'LineWidth', 1.5, 'DisplayName', '6DOF (Actual V_{down})');

% (参考) これまで表示していた w を点線で表示 -> 値が全然違うことがわかるはず
plot(simData.Time, simData.w, 'r:', 'LineWidth', 1.0, 'DisplayName', '6DOF (Body w - Wrong!)');

grid on; legend('Location', 'best');
xlabel('Time [s]'); ylabel('Vertical Velocity [m/s] (Positive = Down)');
title('Correct Descent Rate Comparison');
% --- 結果の取得 ---
t_rk4_result = engine.TimeVector;
y_rk4_result = engine.ResultsMatrix;

% --- 結果の処理 (6_simulatorの形式) ---
simData = table(t_rk4_result, ...
                y_rk4_result(:, 10), y_rk4_result(:, 11), y_rk4_result(:, 12), ...
                y_rk4_result(:, 1), y_rk4_result(:, 2), y_rk4_result(:, 3), ...
                y_rk4_result(:, 4), y_rk4_result(:, 5), y_rk4_result(:, 6), ...
                rad2deg(y_rk4_result(:, 7)), rad2deg(y_rk4_result(:, 8)), rad2deg(y_rk4_result(:, 9)), ...
                engine.PhaseVector, ...
                engine.InertialVelMatrix(:, 1), engine.InertialVelMatrix(:, 2), engine.InertialVelMatrix(:, 3), ...
                engine.EulRateMatrix(:, 1), engine.EulRateMatrix(:, 2), engine.EulRateMatrix(:, 3), ...
                'VariableNames', {'Time', ...
                                 'Inertial_X_Position', 'Inertial_Y_Position', 'Inertial_Z_Position', ...
                                 'Body_U_Vel', 'Body_V_Vel', 'Body_W_Vel', ...
                                 'Body_P_Vel', 'Body_Q_Vel', 'Body_R_Vel', ...
                                 'Roll_Angle', 'Pitch_Angle', 'Yaw_Angle', ...
                                 'Phase', ...
                                 'Inertial_Vx', 'Inertial_Vy', 'Inertial_Vz', ...
                                 'Eul_Phi_Dot', 'Eul_Theta_Dot', 'Eul_Psi_Dot'}); % [cite: 2829-2847]

% --- 結果のプロットと保存 (6_simulatorのプロッターを使用) ---
PayloadPlotter.plotResults(simData, engine.ControlScheduler); % [cite: 2848]
outputFileName = 'parafoil_simulation_results_KawaguchiParams.csv';
writetable(simData, outputFileName, 'Encoding', 'UTF-8'); % [cite: 2849]
disp(['シミュレーションデータを ', outputFileName, ' に保存しました。']);