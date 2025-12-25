%% generate_and_visualize_model.m
% 1. 6自由度パラフォイルモデルをシンボリックに定義
% 2. ヤコビアン(A, B)を解析的に導出
% 3. 数式を人間が読める形式で表示 (Visualization)
% 4. MPC用の高速な関数ファイルを生成 (Generation)

clear; clc;
fprintf('=== 6-DOF Model Linearization & Generation ===\n');

%% 1. シンボリック変数の定義
syms u v w p q r phi theta psi N E D real % 状態 x (12)
syms dL dR real                           % 入力 u (2)
syms rho g real                           % 環境
syms m S b real                           % 機体形状
syms Ixx Iyy Izz Ixz real                 % 慣性

% 空力係数 (線形係数 + 操舵係数)
syms CL0 CL_a CD0 CD_a2 real
syms Clp Cnr Cld Cnd real 

% 変数ベクトル化
x_sym = [u; v; w; p; q; r; phi; theta; psi; N; E; D];
u_sym = [dL; dR];

% パラメータリスト (関数生成時に引数として渡すもの)
params_sym = [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CL_a, CD0, CD_a2, Clp, Cnr, Cld, Cnd];

%% ========================================================================
%% 2. 物理モデルの記述 (ParafoilDynamicsの内容を再現)
%% ========================================================================

% --- 補助変数 ---
V_sq = u^2 + v^2 + w^2;
Va   = sqrt(V_sq);           % 対気速度
alpha = atan2(w, u);         % 迎角 (簡易)
Q_dyn = 0.5 * rho * V_sq;    % 動圧

% 操舵入力
delta_turn = dR - dL;        % 差動 (Roll/Yaw生成)

% --- 力 (Forces) ---
% 揚力・抗力係数
CL = CL0 + CL_a * alpha;
CD = CD0 + CD_a2 * alpha^2;  % ブレーク抵抗項を入れても良い (+ C_D_brake * (dL+dR))

% Body軸上の力 (Wind -> Body変換含む近似)
% Fx ~ -Drag, Fz ~ -Lift
Fx_aero = -Q_dyn * S * CD; 
Fy_aero = 0; % 横滑り力省略
Fz_aero = -Q_dyn * S * CL;

% 重力 (Body Frame)
Fg_x = -m * g * sin(theta);
Fg_y =  m * g * cos(theta) * sin(phi);
Fg_z =  m * g * cos(theta) * cos(phi);

Fx = Fx_aero + Fg_x;
Fy = Fy_aero + Fg_y;
Fz = Fz_aero + Fg_z;

% --- モーメント (Moments) ---
% L (Roll): 減衰(p) + 操舵(dR-dL)
% N (Yaw):  減衰(r) + 操舵(dR-dL)
% ※ 無次元化: p_hat = pb/(2V)
term_p = (p * b) / (2 * Va);
term_r = (r * b) / (2 * Va);

L_aero = Q_dyn * S * b * (Clp * term_p + Cld * delta_turn);
M_aero = 0; % Pitch省略
N_aero = Q_dyn * S * b * (Cnr * term_r + Cnd * delta_turn);

Mx = L_aero;
My = M_aero;
Mz = N_aero;

%% 3. 運動方程式 (Equations of Motion)

% --- 並進 (Newton) ---
dot_u = Fx/m - (q*w - r*v);
dot_v = Fy/m - (r*u - p*w);
dot_w = Fz/m - (p*v - q*u);

% --- 回転 (Euler, xz面対称近似) ---
% 慣性テンソルの逆行列等は複雑になるため、簡易的に主要項のみ記述
% dot_p * Ixx = Mx - (Izz-Iyy)qr + Ixz(dot_r + pq) ...
% ここでは Ixz=0, 対角項のみの簡易形として解く (可視化の見やすさ優先)
% ※ 厳密にやるなら solve() を使いますが、式が爆発的に長くなります。
dot_p = (Mx - (Izz - Iyy)*q*r) / Ixx;
dot_q = (My - (Ixx - Izz)*r*p) / Iyy;
dot_r = (Mz - (Iyy - Ixx)*p*q) / Izz;

% --- 姿勢 (Kinematics) ---
dot_phi   = p + (q*sin(phi) + r*cos(phi))*tan(theta);
dot_theta = q*cos(phi) - r*sin(phi);
dot_psi   = (q*sin(phi) + r*cos(phi))/cos(theta);

% --- 位置 (Navigation) ---
% 簡易的に省略せず記述 (MPCでは使わない場合が多いがモデルとしては正しい)
dot_N = u*cos(theta)*cos(psi); 
dot_E = u*cos(theta)*sin(psi);
dot_D = -u*sin(theta); % 近似

% 全システム f(x, u)
f_sys = [dot_u; dot_v; dot_w; dot_p; dot_q; dot_r; dot_phi; dot_theta; dot_psi; dot_N; dot_E; dot_D];

%% ========================================================================
%% 4. ヤコビアンの導出 (Symbolic Differentiation)
%% ========================================================================
fprintf('Calculating Jacobian A (df/dx)...\n');
A_sym = jacobian(f_sys, x_sym);

fprintf('Calculating Jacobian B (df/du)...\n');
B_sym = jacobian(f_sys, u_sym);

%% ========================================================================
%% 5. 数式の可視化 (Visualization)
%% ========================================================================
fprintf('\n------------------------------------------------------------\n');
fprintf(' Visualization of Linearized Model (at Trim Condition)\n');
fprintf('------------------------------------------------------------\n');

% 数式をきれいにするために、トリム状態(平衡点)の仮定を代入します
% (u=V, v=w=p=q=r=phi=theta=0)
% これをしないと、sin(phi)などの項が残って式が巨大になります
A_disp = subs(A_sym, [v, w, p, q, r, phi, theta], [0, 0, 0, 0, 0, 0, 0]);
B_disp = subs(B_sym, [v, w, p, q, r, phi, theta], [0, 0, 0, 0, 0, 0, 0]);

% さらに Va -> u (トリム速度) として整理
A_disp = simplify(subs(A_disp, Va, u));
B_disp = simplify(subs(B_disp, Va, u));

% 必要な部分 (MPCで使う部分: p, r, phi, psi) だけ抽出して表示
% indices: 4(p), 6(r), 7(phi), 9(psi)
idx_mpc = [4, 6, 7, 9];
labels = {'p', 'r', 'phi', 'psi'};

fprintf('--- Matrix A (Subset for MPC) ---\n');
% A_sub = A_disp(idx_mpc, idx_mpc)
% pretty(A_sub)
% 見やすくするため成分ごとに表示
for i=1:4
    for j=1:4
        term = A_disp(idx_mpc(i), idx_mpc(j));
        if term ~= 0
            fprintf('d(%s)/d(%s) = \n', labels{i}, labels{j});
            pretty(term);
        end
    end
end

fprintf('\n--- Matrix B (Subset for MPC) ---\n');
% B_sub = B_disp(idx_mpc, :)
for i=1:4
    % dL
    termL = B_disp(idx_mpc(i), 1);
    if termL ~= 0
        fprintf('d(%s)/d(dL) = \n', labels{i});
        pretty(termL);
    end
    % dR
    termR = B_disp(idx_mpc(i), 2);
    if termR ~= 0
        fprintf('d(%s)/d(dR) = \n', labels{i});
        pretty(termR);
    end
end

%% ========================================================================
%% 6. 関数ファイルの生成 (Generation)
%% ========================================================================
fprintf('\n------------------------------------------------------------\n');
fprintf(' Generating Optimized Function File...\n');
fprintf('------------------------------------------------------------\n');

% ※可視化用の A_disp ではなく、元の完全な A_sym を保存します
% これにより、トリム速度 V が変わっても計算できるようになります。

matlabFunction(A_sym, B_sym, ...
    'File', 'get_jacobians_sym', ...
    'Vars', {x_sym, u_sym, params_sym}, ...
    'Optimize', true);

fprintf('[SUCCESS] "get_jacobians_sym.m" has been created.\n');
fprintf('You can now run "ParafoilMPC_6DOF" using this file.\n');