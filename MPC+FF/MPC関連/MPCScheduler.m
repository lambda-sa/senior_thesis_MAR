classdef MPCScheduler < handle
    % MPCScheduler (Fixed Version)
    % 修正点: MPC.step メソッドへの引数を (State, Time) の2つに変更
    
    properties
        MpcCtrl      % MPCコントローラ本体 (ParafoilMPC_Integrated)
        RefTable     % 参照軌道 (※旧互換のため保持するが、stepには渡さない)
        WindVec      % 風ベクトル
        GammaRigging % 固定リギング角
        
        LastInput    % 直前の制御入力キャッシュ
        NextMpcTime  % 次にMPC計算を行うべき時刻
        
        phase = 1;   % SimulationEngine用フェーズID
    end
    
    methods
        %% コンストラクタ
        function obj = MPCScheduler(mpc_ctrl, ref_table, wind_vec, gamma)
            obj.MpcCtrl = mpc_ctrl;
            obj.RefTable = ref_table; % (注) 新MPCクラスは内部でこれを持っています
            obj.WindVec = wind_vec;
            obj.GammaRigging = gamma;
            
            obj.LastInput = [0, 0]; 
            obj.NextMpcTime = 0;    
        end
        
        %% 入力取得メソッド
        function inputs = get_inputs(obj, current_time, current_state, ~)
            
            % --- 1. MPC実行判定 ---
            if current_time >= obj.NextMpcTime
                
                try
                    % =========================================================
                    % ★修正箇所: 第3引数 (obj.RefTable) を削除しました
                    % ParafoilMPC_Integrated.step(x, t) に合わせるためです
                    % =========================================================
                    [u_cmd, ~] = obj.MpcCtrl.step(current_state, current_time);
                    
                    obj.LastInput = u_cmd;
                    
                catch ME
                    % エラー時は前回の入力を維持し、警告を表示
                    warning('MPC computation failed at t=%.2f: %s', current_time, ME.message);
                end
                
                obj.NextMpcTime = current_time + obj.MpcCtrl.Ts;
            end
            
            % --- 2. 制御入力構造体の作成 ---
            u_current = obj.LastInput;
            dR = u_current(1);
            dL= u_current(2);
            
            inputs.delta_L = dL;
            inputs.delta_R = dR;
            inputs.GAMMA   = obj.GammaRigging;
            inputs.wind_I  = obj.WindVec;
            
            inputs.delta_s = min(dR, dL);
            inputs.delta_a = dR - dL;
        end
    end
end