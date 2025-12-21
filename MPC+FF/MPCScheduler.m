classdef MPCScheduler < handle
    % MPCScheduler (Complete Version)
    % SimulationEngine (高周波物理演算) と ParafoilMPC_6DOF (低周波制御演算) の仲介役
    % Zero-Order Hold (0次ホールド) 機能を提供します。
    
    properties
        MpcCtrl      % MPCコントローラ本体 (ParafoilMPC_6DOF)
        RefTable     % 参照軌道テーブル
        WindVec      % 風ベクトル [Wx; Wy; Wz] (Inertial Frame)
        GammaRigging % 固定リギング角 (トリム用)
        
        LastInput    % 直前の制御入力キャッシュ [delta_L, delta_R]
        NextMpcTime  % 次にMPC計算を行うべき時刻
        
        phase = 1;   % SimulationEngine用フェーズID (常に1:Guidedとして扱う)
    end
    
    methods
        %% コンストラクタ
        function obj = MPCScheduler(mpc_ctrl, ref_table, wind_vec, gamma)
            % mpc_ctrl: 初期化済みの ParafoilMPC_6DOF インスタンス
            % ref_table: 参照軌道データ
            % wind_vec: 風ベクトル
            
            obj.MpcCtrl = mpc_ctrl;
            obj.RefTable = ref_table;
            obj.WindVec = wind_vec;
            obj.GammaRigging = gamma;
            
            obj.LastInput = [0, 0]; % 初期入力ゼロ
            obj.NextMpcTime = 0;    % 最初は即時実行
        end
        
        %% 入力取得メソッド (SimulationEngineから呼ばれる)
        function inputs = get_inputs(obj, current_time, current_state, ~)
            % GET_INPUTS
            % current_state: 12次元状態ベクトル [u, v, w, p, q, r, phi, theta, psi, N, E, D]
            
            % --- 1. MPC実行判定 (Zero-Order Hold Logic) ---
            % 現在時刻が「次の実行予定時刻」を過ぎていれば計算更新
            if current_time >= obj.NextMpcTime
                
                % ★重要: 12次元ベクトルをそのまま渡す (内部で必要な要素を抽出させる)
                % 以前のエラー原因(抽出ミス)を回避
                try
                    [u_cmd, ~] = obj.MpcCtrl.step(current_state, current_time, obj.RefTable);
                    
                    % 計算成功なら値を更新
                    obj.LastInput = u_cmd;
                    
                catch ME
                    % 万が一MPCがエラーを吐いた場合、前回の値を維持して警告
                    warning('MPC computation failed at t=%.2f: %s', current_time, ME.message);
                end
                
                % 次の実行時刻を更新 (現在時刻 + 周期)
                % ※ "obj.NextMpcTime + Ts" だと遅れが蓄積する場合があるため current_time 基準にする
                obj.NextMpcTime = current_time + obj.MpcCtrl.Ts;
            end
            
            % --- 2. 制御入力構造体の作成 ---
            % SimulationEngine が期待するフォーマットに変換
            
            % 現在有効な入力 (ホールドされている値)
            u_current = obj.LastInput;
            dL = u_current(1);
            dR = u_current(2);
            
            % 必須フィールド
            inputs.delta_L = dL;
            inputs.delta_R = dR;
            inputs.GAMMA   = obj.GammaRigging;
            inputs.wind_I  = obj.WindVec;
            
            % ログ用・互換性用のフィールド (非対称・対称成分)
            inputs.delta_s = min(dR, dL);     % 対称成分 (ブレーキ)
            inputs.delta_a = dR - dL;         % 非対称成分 (ターン)
            
            % 必要に応じて追加 (SimulationEngineの実装依存)
            % inputs.debug_info = ...;
        end
    end
end