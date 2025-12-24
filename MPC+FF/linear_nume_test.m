%% VerifyLinearization.m
% ParafoilDynamics (非線形) と Linear Model (線形近似) の比較検証

clear; clc; close all;

% --- 1. セットアップ ---
% パラメータとモデルの準備
excelFileName = 'parafoil_parameters_SRC.xlsx';
[params, sim_settings] = load_params_from_excel(excelFileName);


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
%params = ParafoilParams(); % 既存のパラメータクラスを使用
model = ParafoilDynamics(params);

% [重要] 線形化モデルと挙動を一致させるため、非線形モデルもSoftminモードにする
% (ParafoilDynamicsの実装に依存しますが、プロパティがあれば設定)
if isprop(model, 'use_soft_min')
    model.use_soft_min = true;
    model.soft_min_k = 30; % 生成した関数と同じk値
end

%% 2. パラメータベクトルの作成 (get_parafoil_jacobians用)
% get_parafoil_jacobians が要求する引数の順番:
% [m, Ixx, Iyy, Izz, Ixz, rho, g, 
%  Sp, Ss, b, c, 
%  xc, yc, zc, xp, yp, zp, Gamma, 
%  CD0, CD_a2, CD_da, CD_ds, 
%  CL0, CL_a, CL_da, 
%  Cm0, Cm_a, Cm_q, 
%  Cl_phi, Cl_p, Cl_da, 
%  Cn_r, Cn_da, Cy_beta, CDS, k_soft]

% --- Softminパラメータ (Excelにない場合はここで定義) ---
val_k_soft = 30;

% --- 構造体から変数へのマッピング ---
% 環境・機体
val_m   = params.m_total;
val_Ixx = params.I_xx;
val_Iyy = params.I_yy;
val_Izz = params.I_zz;
val_Ixz = params.I_xz;
val_rho = params.rho;
val_g   = params.g;

% 形状 (load_params_from_excel の定義に基づく)
val_Sp = params.S_c; % キャノピー面積 (Script内の S_c)
val_Ss = params.S_p; % ペイロード面積 (Script内の S_p)
val_b  = params.b;
val_c  = params.c;

% 位置ベクトル (重心からのオフセット)
% load_params_from_excel で計算された Rcg, Rpg を使用
val_xc = params.Rcg_x; val_yc = 0; val_zc = params.Rcg_z; 
val_xp = params.Rpg_x; val_yp = 0; val_zp = params.Rpg_z; 

% 取り付け角 (sim_settings に読み込まれている Gamma)
val_Gamma = sim_settings.ang_gamma_rigging; 

% 空力係数 (Aero_Coeffs シートから読み込まれた値)
val_CD0   = params.C_D_0;
val_CD_a2 = params.C_D_alpha; % ※モデルが alpha^2 依存か線形依存か要確認
val_CD_da = params.C_D_delta_a;
val_CD_ds = params.C_D_delta_s;

val_CL0   = params.C_L_0;
val_CL_a  = params.C_L_alpha;
val_CL_da = params.C_L_delta_a;

val_Cm0   = params.C_m_0;
val_Cm_a  = params.C_m_alpha;
val_Cm_q  = params.C_m_q;

val_Cl_phi= params.C_l_phi;
val_Cl_p  = params.C_l_p;
val_Cl_da = params.C_l_delta_a;

val_Cn_r  = params.C_n_r;
val_Cn_da = params.C_n_delta_a;

val_Cy_beta = params.C_Y_beta;

% ペイロード空力
val_CDS   = params.C_D_s;

% --- ベクトル結合 ---
params_vec = [ ...
    val_m, val_Ixx, val_Iyy, val_Izz, val_Ixz, val_rho, val_g, ...
    val_Sp, val_Ss, val_b, val_c, ...
    val_xc, val_yc, val_zc, val_xp, val_yp, val_zp, val_Gamma, ...
    val_CD0, val_CD_a2, val_CD_da, val_CD_ds, ...
    val_CL0, val_CL_a, val_CL_da, ...
    val_Cm0, val_Cm_a, val_Cm_q, ...
    val_Cl_phi, val_Cl_p, val_Cl_da, ...
    val_Cn_r, val_Cn_da, val_Cy_beta, val_CDS, val_k_soft ...
];

%% 3. トリム状態の定義
% Excelの Simulation_Settings から初期値を推定、または手動設定
if isfield(sim_settings, 'V_init_guess')
    V_trim = sim_settings.V_init_guess;
else
    V_trim = 20.0; % デフォルト値
end

theta_trim = -10 * (pi/180); % [rad] 滑空角の目安
h_trim = 1000; % [m]

% トリム速度成分 (簡易近似: 迎角0, 横滑り0)
u0 = V_trim * cos(0); 
w0 = V_trim * sin(0); 

% 状態ベクトル x = [u; v; w; p; q; r; phi; theta; psi; x; y; z]
x_trim = zeros(12, 1);
x_trim(1) = u0; 
x_trim(3) = w0;
x_trim(8) = theta_trim;
x_trim(12) = -h_trim; 

% トリム入力 [dL; dR]
% Excelに初期値があればそれを使用
dL_trim = 0; dR_trim = 0;
if isfield(sim_settings, 'delta_L_init'), dL_trim = sim_settings.delta_L_init; end
if isfield(sim_settings, 'delta_R_init'), dR_trim = sim_settings.delta_R_init; end
u_input_trim = [dL_trim; dR_trim];

%% 4. 線形化行列 (A, B) の計算
fprintf('Calculating Linear Model using get_parafoil_jacobians...\n');

% 生成された関数を呼び出し
[A, B] = get_parafoil_jacobians(x_trim, u_input_trim, params_vec);

% トリム点でのドリフト成分 f0 (微係数の残り) の計算
input_struct_trim.delta_L = u_input_trim(1);
input_struct_trim.delta_R = u_input_trim(2);
% 非線形モデル側の実装に合わせてパラメータを渡す
% (ParafoilDynamicsの実装によっては、Gammaや風の渡し方が異なる場合があります)
input_struct_trim.GAMMA = val_Gamma; 
input_struct_trim.wind_I = [0;0;0];

% 非線形モデルで f(x0, u0) を計算
f0 = model.get_derivatives(0, x_trim, input_struct_trim);

%% 5. 検証シミュレーション
T_sim = 10.0;
dt = 0.01;
time = 0:dt:T_sim;
N_steps = length(time);

% 入力: 1秒後に非対称操作 (右ターン)
step_time = 1.0;
delta_add = 0.1; % 入力変化量

Y_nonlinear = zeros(12, N_steps);
Y_linear    = zeros(12, N_steps);

Y_nonlinear(:, 1) = x_trim;
Y_linear(:, 1)    = x_trim;

fprintf('=== Running Simulation (Step Input +%.2f) ===\n', delta_add);

for k = 1:N_steps-1
    t_curr = time(k);
    
    % 入力作成 (ステップ入力)
    if t_curr < step_time
        u_curr = u_input_trim;
    else
        % 右のみ引く操作などを想定
        u_curr = u_input_trim + [0; delta_add]; 
    end
    
    % 偏差入力 delta_u
    du = u_curr - u_input_trim;
    
    % --- A. 非線形モデル (RK4積分) ---
    input_struct_nl.delta_L = u_curr(1);
    input_struct_nl.delta_R = u_curr(2);
    input_struct_nl.GAMMA = val_Gamma;
    input_struct_nl.wind_I = [0;0;0];
    
    state_nl = Y_nonlinear(:,k);
    
    % RK4 Step
    k1 = model.get_derivatives(t_curr, state_nl, input_struct_nl);
    k2 = model.get_derivatives(t_curr + dt/2, state_nl + k1*dt/2, input_struct_nl);
    k3 = model.get_derivatives(t_curr + dt/2, state_nl + k2*dt/2, input_struct_nl);
    k4 = model.get_derivatives(t_curr + dt, state_nl + k3*dt, input_struct_nl);
    
    Y_nonlinear(:, k+1) = state_nl + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    
    % --- B. 線形モデル (Euler積分) ---
    % dx/dt = f(x0,u0) + A*(x-x0) + B*(u-u0)
    dx = Y_linear(:,k) - x_trim;
    dy_lin = f0 + A * dx + B * du;
    
    Y_linear(:, k+1) = Y_linear(:,k) + dy_lin * dt;
end

%% 6. 結果の可視化
fig = figure('Name', 'Linear vs Nonlinear Verification', 'Color', 'w', 'Position', [100 100 1200 800]);

% Yaw Rate r
subplot(2, 3, 1);
plot(time, rad2deg(Y_nonlinear(6,:)), 'b', 'LineWidth', 1.5); hold on;
plot(time, rad2deg(Y_linear(6,:)), 'r--', 'LineWidth', 1.5);
ylabel('Yaw Rate [deg/s]'); title('Yaw Rate (r)'); grid on; legend('Nonlinear', 'Linear');

% Roll Angle phi
subplot(2, 3, 2);
plot(time, rad2deg(Y_nonlinear(7,:)), 'b', 'LineWidth', 1.5); hold on;
plot(time, rad2deg(Y_linear(7,:)), 'r--', 'LineWidth', 1.5);
ylabel('Roll [deg]'); title('Roll Angle (\phi)'); grid on;

% Velocity u
subplot(2, 3, 3);
plot(time, Y_nonlinear(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, Y_linear(1,:), 'r--', 'LineWidth', 1.5);
ylabel('u [m/s]'); title('Forward Velocity u'); grid on;

% Forward Position (North-East)
subplot(2, 3, [4 5]);
plot(Y_nonlinear(11,:), Y_nonlinear(10,:), 'b', 'LineWidth', 1.5); hold on;
plot(Y_linear(11,:), Y_linear(10,:), 'r--', 'LineWidth', 1.5);
xlabel('East [m]'); ylabel('North [m]'); 
title('Ground Track'); grid on; axis equal;

% 誤差ノルム
diff_norm = norm(Y_nonlinear(:,end) - Y_linear(:,end));
fprintf('Final State Difference Norm: %.4f\n', diff_norm);