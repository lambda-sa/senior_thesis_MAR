% run_simulation.m
% パラフォイル 誘導・制御シミュレーション実行スクリプト
% (TAS旋回半径 & 目標点風補正 & 物理整合性確認版)

clear; clc; close all;

%% --- 1. シミュレーション条件設定 ---
excelFileName = 'parafoil_parameters_SRC.xlsx'; 

% 目標着地位置 (North, East, Altitude) [m]
target_pos = [0, 5000, 3000];       

% ファイナルレグ長
L_final = 1000;                

% 風速ベクトル (North, East) [m/s]
% 例: 北へ5m/s, 西へ3m/s (かなり強い風を設定して補正能力を見る)
wind_vector_2d = [3; -3]; 

%% --- 2. ミッション初期化 & 実行 ---
fprintf('=== Simulation Start ===\n');

% インスタンス生成
mission = ParafoilMissionWind(excelFileName);

% 設定 (任意)
mission.set_turn_input_ratio(0.8);
mission.set_minimum_orbit(1);

% ★実行
% 内部で以下が行われます:
% 1. 初期高度でのTAS計算 (Rの決定)
% 2. 沈下時間予測と総ドリフト量の計算
% 3. 目標点のオフセット (Air Target = Target - Drift)
% 4. 経路計画 & 軌道生成
mission.run_wind_simulation(target_pos, L_final, wind_vector_2d);

%% --- 3. パラメータと精度の確認 ---
fprintf('\n=== Verification ===\n');

% A. パラメータ確認
params = mission.Params;
fprintf('>> Parameters Used:\n');
fprintf('   Mass: %.2f kg\n', params.m_total);
if isfield(params, 'S_c'), fprintf('   Wing Area: %.2f m^2\n', params.S_c); end

% B. 着地精度確認
tabGround = mission.export_detailed_trajectory('Ground');
final_pos = tabGround.Position(end, :);
dist_err = norm(final_pos(1:2) - target_pos(1:2));

fprintf('\n>> Landing Accuracy (Ground Frame):\n');
fprintf('   Target:      [%.1f, %.1f]\n', target_pos(1), target_pos(2));
fprintf('   Actual Land: [%.1f, %.1f]\n', final_pos(1), final_pos(2));
fprintf('   Error Dist:  %.2f m\n', dist_err);

if dist_err < 5.0
    fprintf('   [SUCCESS] 風補正は正常に機能しています。\n');
else
    fprintf('   [WARNING] 着地誤差が大きいです。風速と滑空比のバランスを確認してください。\n');
end

%% --- 4. 結果の可視化 ---

% データ取得
tabAir = mission.export_detailed_trajectory('Air');

figure('Name', 'Simulation Result', 'Color', 'w', 'Position', [100, 100, 1200, 800]);

% 1. 高度プロファイル
subplot(2, 3, 1);
plot(tabGround.Time, tabGround.Position(:,3), 'k-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Altitude [m]');
title('Altitude Profile');

% 2. 速度検証 (TAS vs Ground Speed)
subplot(2, 3, 2);
V_TAS = sqrt(sum(tabGround.V_Air.^2, 2));
V_GS  = sqrt(sum(tabGround.V_Ground.^2, 2));
plot(tabGround.Time, V_GS, 'b-', 'DisplayName', 'Ground Speed'); hold on;
plot(tabGround.Time, V_TAS, 'r-', 'LineWidth', 2.0, 'DisplayName', 'True Airspeed (TAS)');
grid on; legend('Location','best'); xlabel('Time [s]'); ylabel('Speed [m/s]');
title('Speed Verification');
% 解説: 赤線(TAS)は高度が高いほど速く、着地(高度0)で遅くなるカーブを描くはずです。

% 3. 迎角検証 (旋回時の増加確認)
subplot(2, 3, 3);
plot(tabGround.Time, rad2deg(tabGround.Alpha), 'm-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Alpha [deg]');
title('Angle of Attack (\alpha)');
% 解説: 旋回中(Bankが入っている間)だけ値が増加していれば、ダイナミクス計算は成功です。

% 4. バンク角
subplot(2, 3, 4);
plot(tabGround.Time, rad2deg(tabGround.Euler_RPY(:,1)), 'b-', 'LineWidth', 1.5);
grid on; xlabel('Time [s]'); ylabel('Bank [deg]');
title('Bank Angle');

% 5. 【重要】軌道の比較 (Top View)
subplot(2, 3, [5 6]);
% 対気軌道 (Air Path)
plot(tabAir.Position(:,2), tabAir.Position(:,1), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Air Path (Wind Frame)');
hold on;
% 対地軌道 (Ground Path)
plot(tabGround.Position(:,2), tabGround.Position(:,1), 'k-', 'LineWidth', 2.0, 'DisplayName', 'Ground Path (Actual)');
% 目標点
plot(target_pos(2), target_pos(1), 'rx', 'MarkerSize', 15, 'LineWidth', 3, 'DisplayName', 'Target (0,0)');
% 対気目標点 (オフセット先)
plot(tabAir.Position(end,2), tabAir.Position(end,1), 'go', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'Offset Aim Point');

axis equal; grid on; xlabel('East [m]'); ylabel('North [m]'); legend('Location','bestoutside');
title(sprintf('Wind Correction Check\nWind: North %.1fm/s, East %.1fm/s', wind_vector_2d(1), wind_vector_2d(2)));

% 解説:
% 緑の点線(Air Path)の最後は、ターゲット(赤×)から風上側に大きくズレた位置(緑〇)で終わっているはずです。
% しかし、黒の実線(Ground Path)は、風に流された結果、見事にターゲット(赤×)に着地するはずです。

%% --- 5. アニメーション (任意) ---
% fprintf('Press any key to start animation...\n');
% pause;
% mission.animate_trajectory(20.0, 'Ground');

%% 2. アニメーション開始
fprintf('\n=== Starting Dual View Animation ===\n');
fprintf('新しいウィンドウが開きます。再生速度は speed_factor で調整可能です。\n');

% 再生速度 (例: 20倍速)
speed_factor =100.0;

% ★新メソッド呼び出し
mission.animate_dual_views(speed_factor);

% フェーズ別の静止画グラフを表示
mission.plot_phases_static();

%% --- 6. 12変数データの出力 ---
fprintf('\n=== Exporting 12-State Variables ===\n');

% 12変数テーブルを取得
state12 = mission.export_6dof_state();

% 内容の確認 (最初の5行)
disp(head(state12, 5));

% ExcelやCSVに保存する場合
% writetable(state12, 'simulation_12states.csv');
% fprintf('Saved to simulation_12states.csv\n');

% グラフで確認 (角速度 p, q, r)
figure('Name', 'Body Rates', 'Color', 'w');
plot(mission.TrajectoryLogGround.Time, state12.p, 'r', 'DisplayName', 'p (Roll Rate)'); hold on;
plot(mission.TrajectoryLogGround.Time, state12.q, 'g', 'DisplayName', 'q (Pitch Rate)');
plot(mission.TrajectoryLogGround.Time, state12.r, 'b', 'DisplayName', 'r (Yaw Rate)');
grid on; legend; xlabel('Time [s]'); ylabel('Rad/s');
title('Body Angular Rates');