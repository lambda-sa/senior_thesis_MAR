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
            [plant_d, ~] = obj.get_linear_model_discrete(ts);
            obj.Ad = plant_d.A;
            obj.Bd = plant_d.B;
            
            % --- 重み行列の設定 ---
            % x = [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            w_vec = zeros(1, 12);
            w_vec(1:3) = 1.0;     % 速度追従
            w_vec(7)   = 10.0;    % phi (ロール安定)
            w_vec(9)   = 50.0;    % psi (方位維持を重視)
            w_vec(10:11) = 20.0;  % N, E (位置ずれ補正)
            w_vec(12)  = 0.0;     % D (高度は成り行きでも良いが、着陸点合わせるなら重み付け)
            
            obj.Q = diag(w_vec);
            obj.P = obj.Q;           % 終端コスト
            obj.R = diag([5.0, 5.0]); % 入力変化抑制

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
            obj.update_linear_model_from_initial_ref();

        end

        %% ★新規: 初期参照値から線形モデルを更新するメソッド
        function update_linear_model_from_initial_ref(obj)
            if isempty(obj.RefFullTable)
                warning('参照軌道がないため線形化できません。');
                return;
            end
            
            % --- A. 参照軌道の1行目(初期値)を取得 ---
            ref_init = obj.RefFullTable(1, :);
            
            % y_trim ベクトル (12x1) の構築
            % 順序: [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            y_trim = zeros(12, 1);
            y_trim(1) = ref_init.u;
            y_trim(2) = ref_init.v;
            y_trim(3) = ref_init.w;
            y_trim(4) = ref_init.p;
            y_trim(5) = ref_init.q;
            y_trim(6) = ref_init.r;
            y_trim(7) = ref_init.phi;
            y_trim(8) = ref_init.theta;
            y_trim(9) = ref_init.psi;
            y_trim(10)= ref_init.N;
            y_trim(11)= ref_init.E;
            y_trim(12)= ref_init.D; % ※RefFullTableのDが高度反転値になっているか確認
            
            % --- B. 入力トリムの設定 ---
            % 通常、参照軌道生成時の入力は左右0（滑空）か、トリム位置
            u_trim = [0; 0]; 
            
            % --- C. 風情報の取得 ---
            % 本来は外部からセットすべきですが、ここではゼロまたはDynamicsModel内の情報を使う
            % ※必要なら set_wind メソッド等で事前にセットしておく
            wind_trim = [0; 0; 0]; 
            
            fprintf('MPC: Updating Linear Model using Initial Reference State...\n');
            fprintf('     V_trim=%.1f m/s, Theta=%.1f deg\n', y_trim(1), rad2deg(y_trim(8)));
            
            % --- D. 解析的線形化の実行 ---
            [A_cont, B_cont, ~] = obj.Linearizer.get_linear_model(0, y_trim, u_trim, wind_trim);
            
            % --- E. 離散化 (c2d) ---
            sys_c = ss(A_cont, B_cont, eye(12), zeros(12,2));
            sys_d = c2d(sys_c, obj.Ts);
            
            obj.Ad = sys_d.A;
            obj.Bd = sys_d.B;
            
            fprintf('MPC: Linear Model Updated successfully.\n');
        end
        
        %% Stepメソッド (QPソルバー実行)
        function [u_cmd, info] = step(obj, x_curr, current_time)
            
            if isempty(obj.RefFullTable)
                error('Reference trajectory not loaded. Call load_reference() first.');
            end
            
            % 1. ホライゾン分の参照軌道を取得
            xref_horizon = obj.get_horizon_ref(current_time, x_curr(9));
            xref_vec = xref_horizon(:); % ベクトル化
            
            % 2. QP行列構築 (Condensing)
            [H, f] = obj.construct_qp_matrices(x_curr, xref_vec);
            
            % 3. 制約条件 (0 <= u <= 1)
            lb = zeros(2 * obj.N, 1);
            ub = ones(2 * obj.N, 1);
            
            % 4. ソルバー実行
            options = optimoptions('quadprog', 'Display', 'off');
            
            % 前回値を初期値として使うと少し速い場合がある(quadprogのアルゴリズムによる)
            [U_seq, ~, exitflag] = quadprog(H, f, [], [], [], [], lb, ub, [], options);
            
            if exitflag < 0
                warning('QP Solver failed (Flag: %d). Using previous input.', exitflag);
                u_cmd = obj.LastU';
            else
                u_cmd = U_seq(1:2)';
                obj.LastU = u_cmd(:);
            end
            
            info.ExitFlag = exitflag;
            info.PredictedPath = []; % 必要ならここに予測軌道を格納
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