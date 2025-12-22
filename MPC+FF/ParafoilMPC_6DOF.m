classdef ParafoilMPC_6DOF < handle
    % PARAFOILMPC_6DOF (Complete Version)
    % 6自由度シンボリック線形化モデルに基づくMPCコントローラ
    %
    % 依存ファイル: get_jacobians_sym.m (generate_linear_model.m で生成)
    
    properties
        MpcObj            % MPCオブジェクト (Model Predictive Control Toolbox)
        MpcState          % MPC内部状態推定器
        Ts                % 制御周期 [s]
        PredictionHorizon % 予測ホライゾン [steps]
        ParamsList        % get_jacobians_sym に渡すパラメータ配列
        CurrentU          % 直前の入力 [dL, dR]
    end
    
    methods
        %% コンストラクタ
        function obj = ParafoilMPC_6DOF(ts, prediction_horizon, params)
            % ts: 制御周期 (例: 0.1)
            % prediction_horizon: 予測ステップ数 (例: 20)
            % params: 機体諸元構造体 (Excelから読み込んだもの)
            
            if nargin < 1, ts = 0.1; end
            if nargin < 2, prediction_horizon = 20; end
            
            obj.Ts = ts;
            obj.PredictionHorizon = prediction_horizon;
            obj.CurrentU = [0, 0];
            
            % 0. 依存ファイルのチェック
            if exist('get_jacobians_sym', 'file') ~= 2
                error('エラー: "get_jacobians_sym.m" が見つかりません。先に "generate_linear_model.m" を実行してください。');
            end
            
            % 1. パラメータリストの構築
            % generate_linear_model.m の定義順序と厳密に合わせる必要があります
            % [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CL_a, CD0, CD_a2, Clp, Cnr, Cld, Cnd]
            
            % ※ 構造体のフィールド名は、お使いの環境に合わせて修正してください
            m = params.m_total;
            S = params.S_c;
            b = params.b;
            
            % 慣性モーメント (なければ概算)
            if isfield(params, 'Ixx'), Ixx = params.Ixx; else, Ixx = m*(b/4)^2; end
            if isfield(params, 'Iyy'), Iyy = params.Iyy; else, Iyy = m*(b/5)^2; end % 仮
            if isfield(params, 'Izz'), Izz = params.Izz; else, Izz = m*(b/3)^2; end
            if isfield(params, 'Ixz'), Ixz = params.Ixz; else, Ixz = 0; end
            
            % 空力係数 (なければ一般的な初期値)
            rho = 1.225; g = 9.81;
            % 空力係数 (パラフォイル固有)
            CL0=params.C_L_0; CLa=params.C_L_alpha;
            CD0=params.C_D_0; CDa2=params.C_D_alpha;
            Clp=params.C_l_p; Cnr=params.C_n_r; Cld=params.C_l_delta_a; Cnd=params.C_n_delta_a;
            
            % リスト化
            obj.ParamsList = [m, Ixx, Iyy, Izz, Ixz, rho, S, b, g, CL0, CLa, CD0, CDa2, Clp, Cnr, Cld, Cnd];
            
            % 2. 線形モデルの取得 (シンボリック関数呼び出し)
            plant = obj.get_symbolic_linear_model(ts);
            
            % 3. MPCオブジェクト作成
            obj.MpcObj = mpc(plant, ts);
            
            % --- MPCパラメータ設定 ---
            obj.MpcObj.PredictionHorizon = prediction_horizon;
            obj.MpcObj.ControlHorizon = 3; 
            
            % 制約: 入力は 0.0 ~ 1.0 (ブレーク量)
            obj.MpcObj.MV(1).Min = 0; obj.MpcObj.MV(1).Max = 1.0;
            obj.MpcObj.MV(2).Min = 0; obj.MpcObj.MV(2).Max = 1.0;
            
            % 制約: 急激な操作の禁止 (レートリミット)
            obj.MpcObj.MV(1).RateMin = -2.0; obj.MpcObj.MV(1).RateMax = 2.0;
            obj.MpcObj.MV(2).RateMin = -2.0; obj.MpcObj.MV(2).RateMax = 2.0;
            
            % 重みづけ (チューニングポイント)
            % [入力] 0.0~1.0の範囲なので、最大振れ幅は1.0とする
            max_u_val  = 1.0; 
            
            % [入力変化率] 物理制約が2.0/sなら、1ステップ(Ts)あたりの最大変化量は 2.0*Ts
            max_du_val = 2.0 * ts; 
            
            % [出力: ロール角(phi)] パラフォイルは揺れやすいので少し広めに許容 (例: 10度)
            max_phi_err = deg2rad(10); 
            
            % [出力: 方位角(psi)] 軌道追従の要なので厳しく設定 (例: 5度)
            max_psi_err = deg2rad(5);

            % 2. 重み係数の計算 (1 / max^2)
            % ---------------------------------------------------------
            w_u   = 1 / (max_u_val^2);
            w_du  = 1 / (max_du_val^2);
            w_phi = 1 / (max_phi_err^2);
            w_psi = 1 / (max_psi_err^2);

            % 入力(MV): 小さくして動きやすくする
            obj.MpcObj.Weights.ManipulatedVariables = [w_u, w_u];
            % 入力変化率(MVRate): ガタつきを抑える
            obj.MpcObj.Weights.ManipulatedVariablesRate = [1.0, 1.0];
            % 出力(OV): [phi, psi] の追従重要度 (psi優先)
            obj.MpcObj.Weights.OutputVariables = [w_phi, w_psi];
            
            % 4. 状態推定器の初期化
            obj.MpcState = mpcstate(obj.MpcObj);
        end
        
        %% Step実行メソッド
        function [u_cmd, info] = step(obj, current_x_full, current_time, ref_table)
            % current_x_full: 12次元状態ベクトル [u,v,w, p,q,r, phi,theta,psi, N,E,D]
            % ref_table: 参照軌道テーブル
            
            % 1. 状態量の抽出
            % MPC内部モデルは [p, r, phi, psi] の4状態
            % ※ 12次元ベクトルの定義順序: 4:p, 6:r, 7:phi, 9:psi
            p_curr   = current_x_full(4);
            r_curr   = current_x_full(6);
            phi_curr = current_x_full(7);
            psi_curr = current_x_full(9);
            
            % 観測値 (Measured Output): [phi, psi]
            y_measure = [phi_curr; psi_curr];
            
            % 2. 参照信号 (Lookahead Reference) の作成
            horizon_times = current_time + (0 : obj.PredictionHorizon)' * obj.Ts;
            
            % 線形補間で未来の目標値を取得
            % Euler_RPY = [roll, pitch, yaw]
            ref_phi = interp1(ref_table.Time, ref_table.Euler_RPY(:,1), horizon_times, 'linear', 'extrap');
            ref_psi = interp1(ref_table.Time, ref_table.Euler_RPY(:,3), horizon_times, 'linear', 'extrap');
            
            % --- 方位角の Unwrap (2pi問題の解決) ---
            % 現在の方位 psi_curr に対して、目標値 ref_psi が近くなるように 2pi ずらす
            % 例: 現在 350度(-10度)、目標 10度 -> 差分20度として扱う (差分340度ではなく)
            
            % 基準偏差
            diff_psi = ref_psi(1) - psi_curr;
            diff_psi = atan2(sin(diff_psi), cos(diff_psi)); % -pi ~ pi
            
            % オフセット後の開始点
            start_ref_psi = psi_curr + diff_psi;
            
            % 未来の軌道形状(増分)を維持したままシフト
            ref_psi_relative = ref_psi - ref_psi(1);
            ref_psi_adjusted = start_ref_psi + ref_psi_relative;
            
            ref_signal = [ref_phi, ref_psi_adjusted];
            
            % 3. 内部状態の更新
            % 非線形シミュレータ(Plant)の値をMPCの推定器に教える
            x_mpc = [p_curr; r_curr; phi_curr; psi_curr];
            obj.MpcState.Plant = x_mpc;
            
            % 4. 最適化計算
            [u_opt, info] = mpcmove(obj.MpcObj, obj.MpcState, y_measure, ref_signal);
            
            obj.CurrentU = u_opt;
            u_cmd = u_opt; % [delta_L, delta_R]
        end
        
        %% 線形モデル取得 (内部メソッド)
        function plant = get_symbolic_linear_model(obj, Ts)
            % トリム条件 (動作点) の設定
            % この条件周りでのみ A, B 行列が正確になります
            
            % 平衡飛行状態 (Level Flight)
            V_trim = 15.0; % [m/s]
            
            % x_trim = [u; v; w; p; q; r; phi; theta; psi; N; E; D]
            x_trim = zeros(12, 1);
            x_trim(1) = V_trim; % u = V
            % 必要であれば theta_trim = alpha_trim を計算して入れる
            
            u_trim = [0; 0]; % 入力ゼロ付近
            
            % ★ 生成された関数 get_jacobians_sym を呼び出す
            % 引数: (x, u, params)
            [A_full, B_full] = get_jacobians_sym(x_trim, u_trim, obj.ParamsList);
            
            % MPC用サブシステムへの縮退
            % 全12状態から [p, r, phi, psi] (Indices: 4, 6, 7, 9) を抽出
            idx_mpc = [4, 6, 7, 9];
            
            A_mpc = A_full(idx_mpc, idx_mpc);
            B_mpc = B_full(idx_mpc, :);
            
            % C行列: 出力は phi(3番目) と psi(4番目)
            C_mpc = [0, 0, 1, 0;
                     0, 0, 0, 1];
                 
            D_mpc = [0, 0;
                     0, 0];
            
            % 状態空間モデル作成
            plant_c = ss(A_mpc, B_mpc, C_mpc, D_mpc);
            plant_c.StateName = {'p', 'r', 'phi', 'psi'};
            plant_c.InputName = {'delta_L', 'delta_R'};
            plant_c.OutputName = {'phi', 'psi'};
            
            % 離散化
            plant = c2d(plant_c, Ts);
            
            % 可制御性チェック (デバッグ用)
            % if rank(ctrb(plant)) < 4, warning('Linearized model is uncontrollable!'); end
        end
    end
end