% run_clothoidVF.m
% 軌道計画 -> FF制御 -> 6自由度シミュレーション 統合実行スクリプト
clear; clc; close all;

%% --- 1. 共通設定 ---
excelFileName = 'parafoil_parameters_ref.xlsx'; % パラメータファイル
wind_planned = [0; 0];  % 風速 (North, East) [m/s]
wind_disturbance = [0;0];  % 予期せぬ外乱（未知の風）
wind_actual = wind_planned + wind_disturbance; % 実際に吹いている風

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
plot3(simData.Inertial_Y_Position, simData.Inertial_X_Position, -simData.Inertial_Z_Position, ...
      'b-', 'LineWidth', 2.0, 'DisplayName', '6DOF (Actual)');
plot3(target_pos(2), target_pos(1), target_pos(3), 'rx', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Target');
grid on; axis equal; view(3);
xlabel('East [m]'); ylabel('North [m]'); zlabel('Altitude [m]');
legend('Location', 'best');
title(sprintf('Trajectory Comparison\nWind: N=%.1f, E=%.1f m/s', wind_planned(1), wind_planned(2)));

% 2. バンク角追従確認 (右上: 2)
subplot(3, 2, 2);
plot(t_plan, rad2deg(phi_plan), 'g--', 'LineWidth', 1.5, 'DisplayName', 'Command (Planner)');
hold on;
plot(simData.Time, simData.Roll_Angle, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Response (6DOF)');
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
plot(simData.Time, -simData.Inertial_Z_Position, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Act (6DOF)');
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
final_pos_6dof = [simData.Inertial_Y_Position(end), simData.Inertial_Z_Position(end)];
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
    phi   = simData.Roll_Angle(i);
    theta = simData.Pitch_Angle(i);
    psi   = simData.Yaw_Angle(i);
    
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
    V_body = [simData.Body_U_Vel(i); simData.Body_V_Vel(i); simData.Body_W_Vel(i)];
    
    % 座標変換
    V_ned_6dof(i, :) = (R_b2n * V_body)';
end

% --- 2. 参照軌道 (Planner) の対地速度計算 ---
% 対気速度 [V_air] を変換し、風ベクトルを加算して対地速度にする
V_ned_plan = zeros(length(t_plan), 3);
wind_3d = [wind_planned; 0]; % 風ベクトル [North, East, Down] (Downは通常0)

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


simData = engine.create_results_table();
PayloadPlotter.plotResults(simData, engine.ControlScheduler); % [cite: 2848]
%{
%% --- Debug: Loiter Vector Field Visualization ---
if strcmpi(autopilot.CurrentMode, 'Loiter') || ~isempty(autopilot.LoiterParams)
    fprintf('\n=== Debugging Loiter Guidance ===\n');
    
    % 1. Loiterパラメータの取得
    p = autopilot.LoiterParams;
    if isempty(p)
        % もしAutopilot内に残ってなければ、Plannerの結果から再構築
        if isfield(mission.Planner.ResultData, 'loiter')
             lx = mission.Planner.ResultData.loiter.x;
             ly = mission.Planner.ResultData.loiter.y;
             p.xc = (min(lx)+max(lx))/2; p.yc = (min(ly)+max(ly))/2; 
             p.R = (max(lx)-min(lx))/2;
             % 方向はとりあえず CW(1) と仮定して確認
             p.lambda = 1; 
        end
    end

    figure('Name', 'Loiter Guidance Debug', 'Color', 'w', 'Position', [100, 100, 800, 800]);
    hold on; grid on; axis equal;
    
    % 2. 目標円の描画
    th = 0:0.01:2*pi;
    cx = p.xc + p.R * cos(th);
    cy = p.yc + p.R * sin(th);
    plot(cy, cx, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Target Circle'); % Plotは (East, North)
    
    % 3. 軌跡の描画
    plot(simData.Inertial_Y_Position, simData.Inertial_X_Position, 'b-', 'LineWidth', 1.2, 'DisplayName', 'Trajectory');
    
    % 4. ベクトル（矢印）の描画
    % 全データだと重すぎるので間引く
    step_idx = 1:50:height(simData); 
    
    for i = step_idx
        % 位置 (NED -> Plot: Y, X)
        pos_n = simData.Inertial_X_Position(i);
        pos_e = simData.Inertial_Y_Position(i);
        
        % 実際のヘディング (青矢印)
        psi = simData.Yaw_Angle(i);
        scale = 40; % 矢印の長さ
        quiver(pos_e, pos_n, scale*sin(psi), scale*cos(psi), ...
            'Color', 'b', 'LineWidth', 1, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
        
        % --- ここで当時の誘導計算を再現 ---
        % (ログに psi_cmd が残っていればそれを使うのがベストですが、再現計算します)
        
        % 1. 基本コース角 (Chi)
        dx = pos_n - p.xc; dy = pos_e - p.yc;
        dist = sqrt(dx^2 + dy^2);
        phi_pos = atan2(dy, dx);
        tilde_d = (dist - p.R) / p.R;
        k = autopilot.Gains.k_vf_loiter;
        lambda = p.lambda;
        
        % 収束角
        correction = (pi/2) + atan(k * tilde_d);
        chi_cmd = phi_pos + lambda * correction;
        
        % 2. 風補正 (Crab)
        % 風データは Scheduler の Est を使う想定 (定常風と仮定)
        if exist('wind_actual', 'var')
            Wx = wind_actual(1); Wy = wind_actual(2);
        else
            Wx = 0; Wy = 0;
        end
        W_mag = norm([Wx, Wy]);
        chi_wind = atan2(Wy, Wx);
        
        % V_horiz 推定 (ログになければ近似)
        u=simData.Body_U_Vel(i); v=simData.Body_V_Vel(i); w=simData.Body_W_Vel(i);
        theta = simData.Pitch_Angle(i);
        % キネマティクス推定を簡易再現
        % V_air ~ V_ground - Wind
        % ここでは簡易的に V_horiz = V_ground_h (風が弱い場合) 
        % または 前回の議論の V_tas * cos(theta-alpha)
        V_tas = sqrt(u^2+v^2+w^2);
        V_horiz = V_tas * cos(theta); % 簡易版
        
        if V_horiz > 1
            W_cross = W_mag * sin(chi_cmd - chi_wind);
            sin_eta = max(-0.95, min(0.95, W_cross/V_horiz));
            eta = asin(sin_eta);
            psi_cmd = chi_cmd - eta;
        else
            psi_cmd = chi_cmd;
        end
        
        % 目標ヘディング (緑矢印)
        quiver(pos_e, pos_n, scale*sin(psi_cmd), scale*cos(psi_cmd), ...
            'Color', 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5, 'HandleVisibility', 'off');
            
        % 風向 (赤矢印, 短く)
        if i == step_idx(1) % 凡例用に1回だけ描画
             quiver(pos_e, pos_n, scale*0.5*sin(chi_wind), scale*0.5*cos(chi_wind), ...
            'Color', 'r', 'LineWidth', 1, 'DisplayName', 'Wind Dir');
        else
             quiver(pos_e, pos_n, scale*0.5*sin(chi_wind), scale*0.5*cos(chi_wind), ...
            'Color', 'r', 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end
    
    % 凡例用ダミープロット
    plot(nan, nan, 'b-', 'DisplayName', 'Actual Heading');
    plot(nan, nan, 'g-', 'DisplayName', 'Command Heading');
    
    xlabel('East [m]'); ylabel('North [m]');
    title('Loiter Debug: Blue=Actual, Green=Command, Red=Wind');
    legend('Location', 'bestoutside');
    xlim([p.yc-p.R*1.5, p.yc+p.R*1.5]);
    ylim([p.xc-p.R*1.5, p.xc+p.R*1.5]);
end

%}