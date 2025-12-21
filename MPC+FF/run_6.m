% --- main_6dof_adapter.m ---
clearvars; close all; clc;

% 1. Excelからパラメータを読み込む
excelFileName = 'parafoil_parameters_ref.xlsx'; % ファイル名は適切に変更してください
[params, sim_settings] = load_params_from_excel(excelFileName);

% 2. 6自由度モデル用に params を拡張・整形する
%    (ParafoilDynamicsクラスは params.prop や params.I_total_body を要求するため)

% --- 慣性行列の構築 ---
params.I_total_body = [
    params.I_xx,    0,          params.I_xz;
    0,              params.I_yy, 0;
    params.I_xz,    0,          params.I_zz
];

% --- prop 構造体の構築 ---
% 6-DOFモデル内の関数が参照するプロパティをマッピング
prop = struct();
prop.TotalMass = params.m_total;
prop.b = params.b;
prop.c = params.c;
prop.d = params.d;

% 位置ベクトル (Excelの重心位置情報を使用)
% ※ 3-DOFでは (Rcg_x, Rcg_z) しか使いませんでしたが、
%    6-DOFでは以下の3つのベクトルが必要です。Excelの値から構成します。
%    ここでは「全重心(CG)を原点」とした相対位置ベクトルとして定義します。

% 例: 全重心からキャノピー空力中心(AC)へのベクトル
%     (Excelの Rcg_x, Rcg_z は CG -> Canopy の距離として定義されていると仮定)
%     符号に注意してください: 
%     通常、キャノピーはCGより「上(-Z)」かつ「後ろ(-X)または前」にあります。
prop.r_total_cm_to_canopy_origin_B = [
    params.Rcg_x;  % X方向距離
    0;             % Y方向 (対称なら0)
    params.Rcg_z   % Z方向 (通常は負の値: 上方向)
];

% ペイロード重心位置 (CGからの距離)
prop.r_total_cm_to_payload_cm_B = [
    params.Rpg_x; 
    0; 
    params.Rpg_z
];

% キャノピーの回転中心から空力中心へのオフセット (なければ0)
prop.r_canopy_origin_to_ac_B = [0; 0; 0]; 

% params に prop を統合
params.prop = prop;


% 3. シミュレーションの準備
% --------------------------------------------------------
% ダイナミクスモデルの初期化
dynamics_model = ParafoilDynamics(params);

% 初期条件 (Excelの sim_settings から設定)
% 状態ベクトル: [u, v, w, p, q, r, phi, theta, psi, x, y, z]
y0 = zeros(12, 1);

% 速度 (対気速度 V_init から u へ。簡単のため u=V, v=w=0 と仮定)
y0(1) = 6; 
y0(2) = 0.0;   % v 
y0(3) = 3.0;   % w 
y0(4) = 0.0;   % p 
y0(5) = 0.0;   % q 
y0(6) = 0.0;   % r 
% 姿勢角
y0(7) = 0; % phi
y0(8) = 0; % theta
y0(9) = deg2rad(sim_settings.psi_initial_deg); % psi

% 位置
y0(10) = sim_settings.X_initial;
y0(11) = sim_settings.Y_initial;
y0(12) = -sim_settings.h_init; % NED座標系なので高度はマイナス

% シミュレーション設定
h_step = 0.05;
t_max  = sim_settings.t_max;

% スケジューラ (制御なし)
%scheduler = ZeroScheduler(); 
scheduler = WindStepScheduler(sim_settings);
% エンジン初期化と実行
fprintf('--- 6-DOF シミュレーション開始 ---\n');
engine = SimulationEngine(dynamics_model, scheduler, h_step, t_max, y0);
engine = engine.run();

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