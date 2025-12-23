classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_6DOF
    % ParafoilMissionWindの生成データを参照軌道として追従する
    % "Scratch-pad" MPC implementation using quadprog
    
    properties
        Ts                % 制御周期
        N                 % 予測ホライゾン
        ParamsList        % 物理パラメータ
        
        % 線形モデル (離散時間)
        Ad
        Bd
        
        % 重み行列
        Q, R, P
        
        % 状態保持
        LastU
        RefFullTable      % 参照軌道テーブル

        % ★追加: 線形化計算用オブジェクト
        Linearizer 
        DynamicsModel % 線形化に必要な物理モデル

        % ★追加: 線形化モード設定
        % 'Fixed'    : 参照軌道の初期値で1回だけ線形化 (計算高速・安定)
        % 'Adaptive' : 毎ステップ現在の状態で線形化 (精度高い・計算重い)
        LinearizationMode = 'Fixed';
    end
    
    methods
        %% コンストラクタ
        function obj = ParafoilMPC_6DOF(ts, horizon, params)
            if nargin < 1, ts = 0.1; end
            if nargin < 2, horizon = 20; end
            
            obj.Ts = ts;
            obj.N  = horizon;
            obj.LastU = [0; 0]; % 初期入力 (左右ブレークなし)
            
            % --- パラメータ展開 ---
            % (最低限必要なものだけ展開)
            obj.ParamsList = params; 
            
            % --- 線形モデル構築 (初期化時に固定) ---
            % ※本来は逐次線形化(LTV)が良いですが、ここではLTIで実装
            %[plant_d, ~] = obj.get_linear_model_discrete(ts);
            %obj.Ad = plant_d.A;
            %obj.Bd = plant_d.B;
            
           % --- 重み行列の設定 (ブライソンの法則適用) ---
            % Q_ii = 1 / (最大許容誤差)^2
            % これにより、「許容誤差」に達したときのコストが等しく「1」になります。
            
            % 1. 最大許容誤差の定義 (Maximum Allowable Deviations)
            % --------------------------------------------------------
            % [速度 (Body Velocity)]  単位: m/s
            err_u = 2.0;   % 対気速度のズレ: ±2 m/s くらいまでは許容
            err_v = 1.0;   % 横滑り速度: 横風で ±1 m/s 程度はずれる
            err_w = 1.0;   % 降下速度: ±1 m/s の変動
            
            % [角速度 (Angular Rates)] 単位: rad/s
            % ※ 10 deg/s = 0.17 rad/s 程度の揺れは許容範囲
            err_p = deg2rad(10); 
            err_q = deg2rad(10);
            err_r = deg2rad(10);
            
            % [姿勢角 (Euler Angles)] 単位: rad
            % ※ ロール・ピッチは ±15度ずれると怖いが、制御上は許容範囲とする
            err_phi   = deg2rad(15); % ロール
            err_theta = deg2rad(15); % ピッチ
            err_psi   = deg2rad(5);  % ヨー (★重要: 方位は ±5度以内に抑えたい)
            
            % [位置 (Position NED)] 単位: m
            % ※ 経路から ±5m 以内なら「オンコース」とみなす
            err_N = 5.0; 
            err_E = 5.0; 
            err_D = 10.0; % 高度は XY より多少甘くても良い
            
            % 2. 重みベクトルの計算 (1/x^2)
            % --------------------------------------------------------
            w_vec = zeros(1, 12);
            
            w_vec(1) = 1 / err_u^2;
            w_vec(2) = 1 / err_v^2;
            w_vec(3) = 1 / err_w^2;
            
            w_vec(4) = 1 / err_p^2;
            w_vec(5) = 1 / err_q^2;
            w_vec(6) = 1 / err_r^2;
            
            w_vec(7) = 1 / err_phi^2;
            w_vec(8) = 1 / err_theta^2;
            w_vec(9) = 1 / err_psi^2;
            
            w_vec(10)= 1 / err_N^2;
            w_vec(11)= 1 / err_E^2;
            w_vec(12)= 1 / err_D^2;
            
            % 3. 追加の調整係数 (Priority Factor)
            % ブライソンの法則はあくまで「正規化」なので、
            % さらに「絶対に守りたい項目」には倍率を掛けます
            obj.Q = diag(w_vec);
            obj.P = obj.Q;           % 終端コスト
            % --- 入力重み R もブライソンの法則で設定 ---
            % 入力 u は 0.0 ~ 1.0 の範囲。最大変動幅をどう見るか。
            % 例: 1ステップで 0.2 (20%) 以上の急操作はコストが高いとする
            max_du = 0.2; 
            r_val  = 1 / max_du^2;
            obj.R = diag([r_val, r_val]);

            %★追加: 線形化のための準備
            % 物理パラメータ構造体から ParafoilDynamics インスタンスを生成
            obj.DynamicsModel = ParafoilDynamics(params); 
            obj.Linearizer = ParafoilLinearizerAnalytic(obj.DynamicsModel);
            
            % 一旦ダミーで線形化しておく（エラー防止のため）
            % 後で load_reference が呼ばれたら上書きされる
            obj.Ad = eye(12); 
            obj.Bd = zeros(12, 2);
        end
        
        %% ★重要修正: Missionクラスからデータをロードするメソッド
        function load_reference(obj, missionObj)
            % Missionクラスで計算済みの「高精度な物理整合データ」をロード
            
            % 1. 12変数状態量の取得
            stateData = missionObj.export_6dof_state(); % table形式
            
            % 2. 時間軸の取得 (Ground Trajectoryから)
            trajGround = missionObj.export_detailed_trajectory('Ground');
            timeVec = trajGround.Time;
            
            % 3. テーブル結合して保持
            % MPC内部で扱いやすいように Time 列を先頭に付与
            obj.RefFullTable = [table(timeVec, 'VariableNames', {'Time'}), stateData];
            
            fprintf('MPC: Reference trajectory loaded (%d steps).\n', height(obj.RefFullTable));

            % 2. ★ここで参照軌道の初期値を使ってモデルを再線形化
            %obj.update_linear_model_from_initial_ref();

            % ★重要: モードに関わらず、最初は必ず「参照軌道の初期値」で線形モデルを作る
            obj.update_linear_model(obj.RefFullTable(1, :));

        end

    
        %% Stepメソッド (QPソルバー実行)
        function [u_cmd, info] = step(obj, x_curr, current_time)
            
            if isempty(obj.RefFullTable)
                error('Reference trajectory not loaded. Call load_reference() first.');
            end
            % --- 1. 線形モデルの更新 (Adaptiveモードの場合) ---
            if strcmpi(obj.LinearizationMode, 'Adaptive')
                % 現在の状態に基づいてモデルを更新 (u_prevを使用)
                obj.update_linear_model_from_state(x_curr, obj.LastU);
            end
            % 1. ホライゾン分の参照軌道を取得
            xref_horizon = obj.get_horizon_ref(current_time, x_curr(9));
            xref_vec = xref_horizon(:); % ベクトル化
            
            % 2. QP行列構築 (Condensing)
            [H, f] = obj.construct_qp_matrices(x_curr, xref_vec);
            
            % 3. 制約条件 (0 <= u <= 1)
            lb = zeros(2 * obj.N, 1);
            ub = ones(2 * obj.N, 1);
            
            % --- 5. QPソルバー実行 ---
            options = optimoptions('quadprog', 'Display', 'off', 'Algorithm', 'interior-point-convex');
            
            % 前回値を初期値として使うと少し速い場合がある(quadprogのアルゴリズムによる)
           [U_seq, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, u0_qp, options);
            
            if exitflag < 0
                warning('QP Solver failed (Flag: %d). Using previous input.', exitflag);
                u_cmd = obj.LastU';
            else
                u_cmd = U_seq(1:2)';
                obj.LastU = u_cmd(:);
            end
            
            info.ExitFlag = exitflag;
            info.RefState = xref_horizon(:, 1); % 直近の目標値
        end
        
        %% 内部: 線形モデル更新 (テーブル行から)
        function update_linear_model(obj, ref_row)
            % ref_row: tableの1行 (Time, u, v, ..., D)
            y_trim = zeros(12, 1);
            y_trim(1) = ref_row.u;
            y_trim(2) = ref_row.v;
            y_trim(3) = ref_row.w;
            y_trim(4) = ref_row.p;
            y_trim(5) = ref_row.q;
            y_trim(6) = ref_row.r;
            y_trim(7) = ref_row.phi;
            y_trim(8) = ref_row.theta; 
            y_trim(9) = ref_row.psi;
            y_trim(10)= ref_row.N;
            y_trim(11)= ref_row.E;
            y_trim(12)= ref_row.D;
            
            u_trim = [0; 0]; 
            wind = [0;0;0]; % 必要に応じてMissionから取得した風を入れる
            
            % 解析的線形化
            [A_c, B_c] = obj.Linearizer.get_linear_model(0, y_trim, u_trim, wind);
            
            % 離散化
            sys_d = c2d(ss(A_c, B_c, eye(12), zeros(12,2)), obj.Ts);
            obj.Ad = sys_d.A;
            obj.Bd = sys_d.B;
            
            if strcmpi(obj.LinearizationMode, 'Fixed')
                fprintf('MPC: Linear Model Fixed at V=%.1fm/s, Theta=%.1fdeg\n', y_trim(1), rad2deg(y_trim(8)));
            end
        end
        
        %% 内部: 線形モデル更新 (状態ベクトルから) - Adaptive用
        function update_linear_model_from_state(obj, x_vec, u_prev)
            y_trim = x_vec;
            u_trim = u_prev;
            wind = [0;0;0];
            
            [A_c, B_c] = obj.Linearizer.get_linear_model(0, y_trim, u_trim, wind);
            sys_d = c2d(ss(A_c, B_c, eye(12), zeros(12,2)), obj.Ts);
            obj.Ad = sys_d.A;
            obj.Bd = sys_d.B;
        end

        %% 内部: QP行列構築
        function [H, f] = construct_qp_matrices(obj, x0, xref_vec)
            n = 12; m = 2; Np = obj.N;
            
            % 事前割り当て
            F = zeros(n*Np, n);
            Phi = zeros(n*Np, m*Np);
            
            % 重み行列の拡張 (Kronecker積)
            Q_bar = kron(eye(Np), obj.Q);
            R_bar = kron(eye(Np), obj.R);
            
            % システム行列の累乗計算
            A_pow = eye(n);
            for i = 1:Np
                % Phiの各列を埋める (畳み込み積分的な処理)
                for j = 1:i
                    if i == j
                        term = obj.Bd;
                    else
                        term = (obj.Ad^(i-j)) * obj.Bd;
                    end
                    % 行: (i-1)n+1 ~ in, 列: (j-1)m+1 ~ jm
                    Phi((i-1)*n+1 : i*n, (j-1)*m+1 : j*m) = term;
                end
                
                A_pow = obj.Ad * A_pow;
                F((i-1)*n+1 : i*n, :) = A_pow;
            end
            
            % QP形式へ変換: min 1/2 U'HU + f'U
            H = 2 * (Phi' * Q_bar * Phi + R_bar);
            H = (H + H') / 2; % 対称化
            
            Y_free = F * x0;       % 自由応答 (入力ゼロ時の未来)
            E = Y_free - xref_vec; % 誤差ベクトル
            
            f = (2 * E' * Q_bar * Phi)';
        end
        
        %% 内部: ホライゾン参照値取得 (補間 + Psi Unwrap)
        function xref_mat = get_horizon_ref(obj, t_curr, psi_curr)
            % 予測ホライゾンの時刻列
            t_query = t_curr + (1:obj.N)' * obj.Ts;
            
            % データ抽出
            ref_time = obj.RefFullTable.Time;
            ref_data = obj.RefFullTable{:, 2:13}; % Time列を除いた12変数
            
            % 線形補間
            % 'extrap' で範囲外(着陸後)は最後の値を維持
            xref_mat = interp1(ref_time, ref_data, t_query, 'linear', 'extrap');
            
            xref_mat = xref_mat'; % [12 x N]
            
            % Psi (9番目の要素) の不連続性(Unwrap)処理
            idx_psi = 9;
            ref_psi = xref_mat(idx_psi, :);
            
            % 現在のPsiと、参照軌道の最初のPsiの差分を計算し、2piのズレを補正
            d_psi = ref_psi(1) - psi_curr;
            d_psi_wrapped = atan2(sin(d_psi), cos(d_psi));
            offset = psi_curr + d_psi_wrapped - ref_psi(1);
            
            xref_mat(idx_psi, :) = ref_psi + offset;
        end
       
    end
end