% run_simple_test.m
% 単純軌道（直線/旋回）を用いた6DOFシミュレーションテスト
clear; clc; close all;
%% --- 1. 設定と参照軌道の生成 (Phase A) ---
excelFileName = 'parafoil_parameters_SRC.xlsx';
test_duration = 600;   % 秒
bank_cmd_deg  = 5;   % 旋回テスト用のバンク角 [deg]
test_type = "turn";
fprintf('=== Phase A: Reference Generation (TestMission) ===\n');
% 既存の MissionWind の機能を継承したテストクラスを使用
mission = ParafoilTestMission(excelFileName);

% 12変数の参照テーブルを取得 (u,v,w, p,q,r, phi,theta,psi, N,E,D)
% 内部でトリム計算と12変数変換 (export_6dof_state) が行われる
trajRef = mission.run_test_maneuver(test_duration, bank_cmd_deg);

%% --- 3. Phase B & C: 6DOFシミュレーション ---
fprintf('\n=== Phase B: 6DOF Simulation Setup ===\n');

% モデル準備
%params = mission.Params;
params = mission.Params; % Missionが読み込んだパラメータを共有
atmo   = mission.AtmoModel;

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
mapper = ParafoilControlMapper(params);

% --- 初期状態の設定 (参照軌道の最初の1行に完全一致させる) ---
% trajRef(1, :) から 12変数を抽出
u0 = trajRef.u(1); v0 = trajRef.v(1); w0 = trajRef.w(1);
p0 = trajRef.p(1); q0 = trajRef.q(1); r0 = trajRef.r(1);
ph0 = trajRef.phi(1); th0 = trajRef.theta(1); ps0 = trajRef.psi(1);
N0 = trajRef.x(1); E0 = trajRef.y(1); D0 = -trajRef.z(1);

y0 = [u0; v0; w0; p0; q0; r0; ph0; th0; ps0; N0; E0; D0];

% --- スケジューラ用参照データの整理 ---

% 1. 時間軸の取得 (trajRefにTimeがない場合は、元のミッションデータから取得)
if ismember('Time', trajRef.Properties.VariableNames)
    t_plan = trajRef.Time;
else
    % export_6dof_stateの内部で使われていた時間軸を代入
    t_plan = mission.export_detailed_trajectory('Ground').Time;
    trajRef.Time = t_plan; % 後続の処理のためにテーブルに追加しておく
end

% 2. 各状態量の抽出
phi_plan   = trajRef.phi;   % バンク角 [rad]
theta_plan = trajRef.theta; % ピッチ角 [rad]

% 3. 対気速度スカラ (V) の計算
% export_6dof_state は u, v, w を出すため、そのノルムを計算します
V_plan = sqrt(trajRef.u.^2 + trajRef.v.^2 + trajRef.w.^2);
wind_vec = [0; 0];    % Step 1: 制御アルゴリズムの純粋な検証（まずはこれ！）
% 4. スケジューラのインスタンス化
% 定義した変数名を使用して、可読性を高めます
scheduler = PlannerTrackingScheduler(mapper, ...
    t_plan, ...
    phi_plan, ...
    V_plan, ...
    theta_plan, ...
    wind_vec);

% 実行
engine = SimulationEngine(plant, scheduler, 0.05, test_duration, y0);
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
%plot3(trajRef.Position(:,2), trajRef.Position(:,1), trajRef.Position(:,3), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Trim)');
% テーブルの列 (y=East, x=North, z=Alt) を指定
plot3(trajRef.E, trajRef.N, trajRef.D, 'g--', 'LineWidth', 1.5, 'DisplayName', 'Ref (Mission)');
hold on;
plot3(simData.East, simData.North, -simData.Down, 'b-', 'LineWidth', 2, 'DisplayName', '6DOF');
grid on; axis equal; view(3);
xlabel('East'); ylabel('North'); zlabel('Alt');
legend; title(['Trajectory: ' test_type]);

subplot(2,2,2);
plot(trajRef.Time, rad2deg(trajRef.phi), 'g--', 'LineWidth',1.5, 'DisplayName','Ref');
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