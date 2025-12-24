%% generate_parafoil_model_final.m
% 参考文献: 東京理科大学 修士論文「複数の風情報を用いたパラフォイル回収システムの誘導性能評価」
% 目的: 6DoFモデルのシンボリック微分を行い、線形化行列(A,B)と状態微分(dx)の関数を生成する
% 特徴:
%  1. Softmin関数によりブレーキ入力の微分可能性を担保
%  2. 慣性行列の逆行列をシンボリック的に明示計算し、処理落ちを回避
%  3. トリム探索用と線形化用の2つの関数を出力

clear; clc;
fprintf('=== 6-DOF Parafoil Model Linearization (Final Version) ===\n');

%% 1. シンボリック変数の定義
% 状態変数 x (12状態)
syms u v w p q r phi theta psi x_pos y_pos z_pos real
% 入力変数 u_in (2入力)
syms dL dR real 
% Softmin用パラメータ (線形化時の滑らかさを調整)
syms k_soft real 
assume(k_soft > 0);

% --- パラメータ定義 (35個) ---
syms m Ixx Iyy Izz Ixz real         % 機体特性
syms rho g real                     % 環境
syms Sp Ss b c real                 % 幾何形状
syms x_c y_c z_c x_p y_p z_p real   % 重心位置
syms Gamma real                     % 取り付け角
syms CD0 CD_a2 CD_da CD_ds real     % キャノピー抗力
syms CL0 CL_a CL_da real            % キャノピー揚力
syms Cm0 Cm_a Cm_q real             % ピッチモーメント
syms Cl_phi Cl_p Cl_da real         % ロールモーメント
syms Cn_r Cn_da real                % ヨーモーメント
syms Cy_beta real                   % 横力
syms CDS real                       % ペイロード抗力

% ベクトル化
x_state = [u; v; w; p; q; r; phi; theta; psi; x_pos; y_pos; z_pos];
u_input = [dL; dR];

% パラメータリスト (順序重要)
params_sym = [m, Ixx, Iyy, Izz, Ixz, rho, g, ...
              Sp, Ss, b, c, ...
              x_c, y_c, z_c, x_p, y_p, z_p, Gamma, ...
              CD0, CD_a2, CD_da, CD_ds, ...
              CL0, CL_a, CL_da, ...
              Cm0, Cm_a, Cm_q, ...
              Cl_phi, Cl_p, Cl_da, ...
              Cn_r, Cn_da, Cy_beta, CDS, k_soft];

%% 2. Softminによる操作量の定義
% 論文における delta_s = min(dL, dR) を滑らかに近似
delta_a = dR - dL; % 非対称操作量 (ターン)
term_L = exp(-k_soft * dL);
term_R = exp(-k_soft * dR);
delta_s = -(1/k_soft) * log(term_L + term_R); % 対称操作量 (ブレーキ)

%% 3. 座標変換行列
s_phi = sin(phi); c_phi = cos(phi);
s_tht = sin(theta); c_tht = cos(theta);
s_psi = sin(psi); c_psi = cos(psi);

% 慣性系 -> 機体座標系 (T_IB) の転置として定義されることが多いが、
% 重力項で使うため [3x3] を定義
T_IB = [c_tht*c_psi, c_tht*s_psi, -s_tht;
        s_phi*s_tht*c_psi - c_phi*s_psi, s_phi*s_tht*s_psi + c_phi*c_psi, s_phi*c_tht;
        c_phi*s_tht*c_psi + s_phi*s_psi, c_phi*s_tht*s_psi - s_phi*c_psi, c_phi*c_tht];

% 機体座標系 -> キャノピー座標系 (Incidence Angle Gamma)
T_BC = [cos(Gamma), 0, -sin(Gamma);
        0,          1,  0;
        sin(Gamma), 0,  cos(Gamma)];

%% 4. 速度・迎角・横滑り角の計算
V_b = [u; v; w];
Omega_b = [p; q; r];
S_omega_B = [0, -r, q; r, 0, -p; -q, p, 0]; % 外積行列

% キャノピー位置での速度
Delta_c_vec = [x_c; y_c; z_c];
V_canopy_body = V_b + S_omega_B * Delta_c_vec;
V_c_vec = T_BC * V_canopy_body; % キャノピー座標系へ
uc = V_c_vec(1); vc = V_c_vec(2); wc = V_c_vec(3);
Va_c = sqrt(uc^2 + vc^2 + wc^2);

alpha_c = atan2(wc, uc);
beta_c  = asin(vc / Va_c);

% ペイロード位置での速度
Delta_p_vec = [x_p; y_p; z_p];
V_payload_body = V_b + S_omega_B * Delta_p_vec;
us = V_payload_body(1); vs = V_payload_body(2); ws = V_payload_body(3);
Va_s = sqrt(us^2 + vs^2 + ws^2);

%% 5. 力とモーメントの定式化
Q_c = 0.5 * rho * Va_c^2;

% 空力係数
CL = CL0 + CL_a * alpha_c + CL_da * delta_a;
CD = CD0 + CD_a2 * alpha_c^2 + CD_da * delta_a + CD_ds * delta_s;
CY = Cy_beta * beta_c;

% キャノピー空力 (機体座標系へ変換)
T_AC = [cos(alpha_c), 0, -sin(alpha_c);
        0,            1,  0;
        sin(alpha_c), 0,  cos(alpha_c)];
F_aero_frame = [-CD; CY; -CL]; 
F_A_vec = Q_c * Sp * T_BC.' * T_AC * F_aero_frame;

% ペイロード空力 (速度と逆方向)
F_S_vec = -0.5 * rho * Ss * CDS * Va_s * [us; vs; ws];

% 重力 (機体座標系)
F_W_vec = T_IB * [0; 0; m * g];

% 合力
F_total = F_W_vec + F_A_vec + F_S_vec;

% モーメント
Omega_c = T_BC * Omega_b;
pc = Omega_c(1); qc = Omega_c(2); rc = Omega_c(3);
term_p = (pc * b) / (2 * Va_c);
term_q = (qc * c) / (2 * Va_c);
term_r = (rc * b) / (2 * Va_c);

Cm = Cm0 + Cm_a * alpha_c + Cm_q * term_q;
Cl = Cl_phi * phi + Cl_p * term_p + Cl_da * delta_a / b;
Cn = Cn_r * term_r + Cn_da * delta_a / b;

M_A_canopy = Q_c * Sp * [b * Cl; c * Cm; b * Cn];
M_A_vec = T_BC.' * M_A_canopy; % 機体座標系へ

% モーメントアームによる外積項
M_arm_A = cross(Delta_c_vec, F_A_vec);
M_arm_S = cross(Delta_p_vec, F_S_vec);
M_total = M_A_vec + M_arm_A + M_arm_S;

%% 6. 運動方程式 (Equations of Motion)
% 並進: m * (dot_V + omega x V) = F
dot_V_b = (1/m) * F_total - cross(Omega_b, V_b);

% 回転: I * dot_omega + omega x (I * omega) = M
I_mat = [Ixx, 0, Ixz; 0, Iyy, 0; Ixz, 0, Izz];
I_inv = inv(I_mat); % 3x3逆行列を明示的に計算 (高速化の鍵)
dot_Omega_b = I_inv * (M_total - cross(Omega_b, I_mat * Omega_b));

% 姿勢 (Euler Rate)
Mat_Euler = [1, s_phi*s_tht/c_tht, c_phi*s_tht/c_tht;
             0, c_phi,             -s_phi;
             0, s_phi/c_tht,       c_phi/c_tht];
dot_Euler = Mat_Euler * Omega_b;

% 位置 (Navigation)
dot_Pos = T_IB.' * V_b;

% 全システム微分方程式 f(x, u)
f_sys = [dot_V_b; dot_Omega_b; dot_Euler; dot_Pos];

%% 7. ヤコビアンの計算 (Linearization)
fprintf('Calculating Jacobian A (df/dx)...\n');
% simplifyは計算時間がかかりすぎるため省略
A_sym = jacobian(f_sys, x_state);

fprintf('Calculating Jacobian B (df/du)...\n');
B_sym = jacobian(f_sys, u_input);

%% 8. 関数ファイルの生成
fprintf('Generating optimized function files...\n');

% 1. 線形化用: [A, B] = get_parafoil_jacobians(x, u, params)
% Optimize=falseにより、式が長くても確実に生成させる
matlabFunction(A_sym, B_sym, ...
    'File', 'get_parafoil_jacobians', ...
    'Vars', {x_state, u_input, params_sym}, ...
    'Optimize', true);

% 2. トリム探索用: dx = get_parafoil_derivatives(x, u, params)
% 非線形シミュレーションやfminsearchでの釣合い計算に使用
matlabFunction(f_sys, ...
    'File', 'get_parafoil_derivatives', ...
    'Vars', {x_state, u_input, params_sym}, ...
    'Optimize', false);

fprintf('[SUCCESS] Linearization complete.\n');
fprintf('Generated files:\n');
fprintf(' - get_parafoil_jacobians.m (Returns A, B matrices)\n');
fprintf(' - get_parafoil_derivatives.m (Returns dx vector)\n');