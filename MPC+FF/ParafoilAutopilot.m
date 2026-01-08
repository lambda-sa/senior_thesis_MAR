classdef ParafoilAutopilot < handle
    % PARAFOILAUTOPILOT
    % 誘導(Guidance)と制御(Control)を統合した自律飛行クラス
    %
    % [概要]
    % Plannerが作成した軌道データを読み込み、現在の機体状態(位置・速度・風)に基づいて
    % 最適なトグル操作量(delta_R, delta_L)を計算して出力します。
    %
    % [座標系の定義]
    % - NED座標系 (North-East-Down): X=北, Y=東, Z=下
    % - 方位角(Psi): 北を0として時計回り(CW)が正
    %
    % [誘導ロジックの特徴]
    % 1. Loiter (待機): 位相角ベースの解析的ベクトル場 (指定論文の式)
    % 2. Mission (着陸): 点群に対するウィンドウ探索 + General Path VF
    
    properties
        % --- 機体・空力パラメータ ---
        Params       % 全パラメータ構造体 (質量, 翼面積, 空力係数など)
        Yaw_Factor   % トグル操作量 -> ヨーモーメントへの変換係数 (FF制御用)
        
        % --- ゲイン設定 (要調整) ---
        Gains
        % .k_vf_loiter : Loiter時の収束ゲイン (小さい=マイルド, 大きい=急激)
        % .k_vf_mission: Mission時の収束ゲイン (着陸精度に直結)
        % .chi_inf     : 無限遠での最大進入角 [rad] (通常 pi/2 = 90度)
        % .kp_psi      : [制御] ヘディング誤差 -> 目標バンク角 の比例ゲイン
        % .kp_phi      : [制御] バンク角誤差 -> トグル操作量 の比例ゲイン
        
        % --- 軌道データ (import_mission_dataでセット) ---
        CurrentMode     % 現在のモード: 'Loiter' or 'Mission'
        
        % Loiter用: 円軌道の幾何学情報
        LoiterParams    % .xc, .yc (中心), .R (半径), .lambda (+1:右回り, -1:左回り)
        
        % Mission用: 点群データ (Entry -> Dubins -> Final)
        MissionPath     % .x, .y, .z, .chi(接線), .kappa(曲率), .num_points
        LastIndex       % 前回の最近傍点インデックス (探索高速化 & 逆走防止用)
        WindowSize      % 前方探索ウィンドウサイズ (点の個数)
        
        % --- 状態管理 ---
        IsReady         % 初期化完了フラグ
    end
    
    methods
        function obj = ParafoilAutopilot(params)
            % コンストラクタ: パラメータの展開と初期設定
            obj.Params = params;
            obj.IsReady = false;
            
            % -----------------------------------------------------
            % 1. トグル操作量のFF係数 (Yaw_Factor) の計算
            % -----------------------------------------------------
            % 「ある旋回レート(r)を出したい時、どれくらいトグルを引けばいいか？」
            % を、線形化モデルのモーメント釣り合い式から逆算するための係数です。
            % 式: delta_a_req = Yaw_Factor * (g / V^2) * sin(phi)
            
            if isfield(params, 'prop'), b = params.prop.b; else, b = params.b; end
            if isfield(params, 'params')
                d = params.params.d; % 吊り下げ長
                Cn_r = params.params.C_n_r; % ヨー減衰係数
                Cn_da = params.params.C_n_delta_a; % 舵効き係数
            else
                d = params.d; Cn_r = params.C_n_r; Cn_da = params.C_n_delta_a;
            end
            
            % 物理的意味: (減衰項) / (操舵効力項)
            obj.Yaw_Factor = -d * b * Cn_r / (2 * Cn_da);
            
            % -----------------------------------------------------
            % 2. 制御ゲインのデフォルト設定
            % -----------------------------------------------------
            obj.Gains.k_vf_loiter  = 0.04;   % Loiterは滑らかさ重視 (低め)
            obj.Gains.chi_inf      = pi/3;   % 最大90度でコースに戻る
            obj.Gains.k_vf_mission = 0.10;   % Missionは追従精度重視 (高め)
            
            obj.Gains.kp_psi       = 0;    % 方位ズレ1radに対し、バンク1.5rad指令
            obj.Gains.kp_phi       = 0;    % バンクズレ1radに対し、トグル0.5操作
            
            % 内部変数の初期化
            obj.CurrentMode = 'Loiter';
            obj.WindowSize = 100; % 探索範囲は前方100点まで
            obj.LastIndex = 1;
        end
        
        % =========================================================
        % データ取り込み: Plannerの出力結果を解釈する
        % =========================================================
        function import_mission_data(obj, planner)
            % ParafoilPathPlannerClothoid の結果(ResultData)を読み込み、
            % 誘導則が扱いやすい形に変換して保持します。
            
            % 風を考慮した計画データがあればそれを優先
            if isprop(planner, 'ResultDataWind') && ~isempty(planner.ResultDataWind)
                d = planner.ResultDataWind;
            else
                d = planner.ResultData;
            end
            if isempty(d), error('Planner has no data.'); end
            
            % --- A. Loiterデータの解析 (点群 -> 円パラメータ) ---
            % Loiterは点群を追うより、中心と半径から解析的に計算する方が
            % 計算負荷が低く、かつ真円に近い滑らかな旋回ができます。
            if isfield(d, 'loiter') && length(d.loiter.x) > 10
                lx = d.loiter.x; ly = d.loiter.y;
                
                % 最大値・最小値から中心と半径を推定
                xc = (min(lx) + max(lx)) / 2;
                yc = (min(ly) + max(ly)) / 2;
                R  = (max(lx) - min(lx)) / 2;
                
                % 回転方向(Lambda)の自動判定
                % 始点付近の2つのベクトル外積を用いて判定
                % NED系(X:北, Y:東)における外積 (x1*y2 - x2*y1)
                % 時計回り(CW)なら外積はプラスになる傾向があるが、場所による。
                % 確実な方法: 角度の変化を見る
                ang1 = atan2(ly(1)-yc, lx(1)-xc);
                ang2 = atan2(ly(5)-yc, lx(5)-xc);
                diff = ang2 - ang1;
                % 角度補正 (-pi~piのまたぎ処理)
                if diff > pi, diff = diff - 2*pi; end
                if diff < -pi, diff = diff + 2*pi; end
                
                % NED系では時計回りが「角度増加」なので diff > 0 ならCW
                if diff > 0
                    lambda = 1;  % 右旋回 (Clockwise)
                else
                    lambda = -1; % 左旋回 (Counter-Clockwise)
                end
                
                obj.LoiterParams = struct('xc',xc, 'yc',yc, 'R',R, 'lambda',lambda);
            else
                warning('No Loiter data found.');
                obj.LoiterParams = [];
            end
            
            % --- B. Missionデータの結合 (Entry + Dubins + Final) ---
            % 複雑な経路は「一本の点群」として扱います。
            raw_x=[]; raw_y=[]; raw_z=[]; raw_psi=[];
            phases = {'entry', 'dubins', 'final'}; 
            for i = 1:length(phases)
                p = phases{i};
                if isfield(d, p) && ~isempty(d.(p).x)
                    raw_x = [raw_x, d.(p).x(:)'];
                    raw_y = [raw_y, d.(p).y(:)'];
                    raw_z = [raw_z, d.(p).z(:)'];
                    raw_psi = [raw_psi, d.(p).psi(:)'];
                end
            end
            
            if isempty(raw_x)
                warning('No Mission data found.');
                obj.MissionPath = [];
            else
                % 曲率 (Kappa) の数値再計算
                % Plannerが出力する座標から、数値微分で曲率を求めます。
                % これにより、バンク角のフィードフォワード制御が可能になります。
                dx = gradient(raw_x); dy = gradient(raw_y);
                ds = sqrt(dx.^2 + dy.^2); ds(ds<0.1) = 0.1; % ゼロ割防止
                dpsi = gradient(unwrap(raw_psi));
                
                kappa = dpsi ./ ds; % 曲率定義: d(角度)/d(距離)
                kappa = smoothdata(kappa, 'movmean', 5); % 微分ノイズを除去
                
                obj.MissionPath.x = raw_x;
                obj.MissionPath.y = raw_y;
                obj.MissionPath.z = raw_z;
                obj.MissionPath.chi = unwrap(raw_psi);
                obj.MissionPath.kappa = kappa;
                obj.MissionPath.num_points = length(raw_x);
            end
            
            obj.IsReady = true;
            obj.CurrentMode = 'Loiter'; % 初期モード設定
            obj.LastIndex = 1;
            
            fprintf('[Autopilot] Initialized. Loiter: %s (R=%.1fm), Mission: %d pts.\n', ...
                calc_rot_str(obj.LoiterParams.lambda), obj.LoiterParams.R, obj.MissionPath.num_points);
        end
        
        % モード手動切替 (Schedulerから呼ばれる)
        function set_mode(obj, mode_str)
            if strcmpi(mode_str, 'Loiter') || strcmpi(mode_str, 'Mission')
                obj.CurrentMode = mode_str;
            else
                error('Mode must be "Loiter" or "Mission"');
            end
        end
        
        % =========================================================
        % ★★★ メイン更新メソッド (シミュレーション毎ステップ実行) ★★★
        % =========================================================
        function [delta_R, delta_L, log_struct] = update(obj, current_state, wind_vec, delta_s_bias)
            % 入力:
            %   current_state : [u, v, w, p, q, r, phi, theta, psi, x, y, z]
            %   wind_vec      : [Wx, Wy, 0] (NED系)
            %   delta_s_bias  : 対称操作量 (高度制御用)
            
            if nargin < 4, delta_s_bias = 0; end
            if ~obj.IsReady, error('Autopilot not initialized.'); end
            
            % --- 1. 状態量の展開 ---
            u=current_state(1); v=current_state(2); w=current_state(3);
            phi=current_state(7); theta=current_state(8); psi=current_state(9);
            pos_ned = current_state(10:12);
            
            % 対気速度 (TAS) と 対地速度の概算
            V_tas = sqrt(u^2+v^2+w^2); if V_tas<1, V_tas=1; end
            % バンク角FF計算には対地速度(Ground Speed)を使うのが物理的に正確
            V_g_sq = (V_tas*cos(theta))^2; 
            
            % -----------------------------------------------------
            % Step 1: 誘導 (Guidance)
            % -----------------------------------------------------
            % 現在位置に基づき、目指すべき「対地コース(Chi)」と「曲率(Kappa)」を計算
            cmd = obj.run_vf_guidance(pos_ned);
            
            % -----------------------------------------------------
            % Step 2: 制御目標の計算 (Control Target)
            % -----------------------------------------------------
            % (A) 風補正 (Crab Angle Compensation)
            % VFが出した「対地コース」通りに進むために、機首を何度ずらすべきか？
            chi_cmd = cmd.chi_cmd;
            
            % 風向・風速から横風成分を計算
            W_mag = norm(wind_vec(1:2));
            chi_wind = atan2(wind_vec(2), wind_vec(1));
            W_cross = W_mag * sin(chi_cmd - chi_wind); % 進行方向に対する横風
            
            % クラブ角 eta の計算 (sin(eta) = W_cross / V_tas)
            eta = asin(max(-0.9, min(0.9, W_cross/V_tas)));
            
            % 最終的な目標ヘディング角
            psi_cmd = chi_cmd - eta;
            
            % (B) 目標バンク角の決定
            g = 9.81;
            
            % FF項: カーブを曲がるための遠心力釣り合いバンク角
            % tan(phi) = (V^2 / g) * kappa
            phi_ff = atan( (V_g_sq / g) * cmd.kappa_ref ); 
            
            % FB項: ヘディング誤差を埋めるための修正バンク角
            % (PID制御のP項に相当)
            psi_err = obj.normalize_angle(psi_cmd - psi);
            phi_fb = obj.Gains.kp_psi * psi_err;
            
            % 合成とリミッター (パラフォイルは45度以上傾けると危険)
            phi_cmd = max(-0.8, min(0.8, phi_ff + phi_fb));
            
            % -----------------------------------------------------
            % Step 3: アクチュエータ操作 (Actuation)
            % -----------------------------------------------------
            % (A) FF項: 定常旋回操作量
            % 空力モデルに基づき、そのバンク角を維持するのに必要なトグル量を予測
            K_dyn = (g / V_tas^2) * cos(theta);
            da_ff = obj.Yaw_Factor * K_dyn * sin(phi_cmd);
            
            % (B) FB項: バンク角追従誤差補正
            % モデル化誤差や外乱でバンク角がズレた分をPID(P)で戻す
            phi_err = phi_cmd - phi;
            da_fb = obj.Gains.kp_phi * phi_err;
            
            % 合計操作量
            da_total = da_ff + da_fb;
            
            % (C) ミキシング (左右トグルへの配分)
            % ブレーキ量(delta_s)をベースに、旋回成分(da_total)を差動で加える
            [delta_R, delta_L] = obj.apply_mixing(da_total, delta_s_bias);
            
            % -----------------------------------------------------
            % ログ出力
            % -----------------------------------------------------
            log_struct.mode = obj.CurrentMode;
            log_struct.chi_cmd = chi_cmd;
            log_struct.psi_cmd = psi_cmd;
            log_struct.phi_cmd = phi_cmd;
            log_struct.kappa_ref = cmd.kappa_ref;
            log_struct.error = cmd.error;
            log_struct.delta_a = da_total;
        end
    end
    
    methods (Access = private)
        % =========================================================
        % ★★★ 誘導ロジック (Hybrid Vector Field) ★★★
        % =========================================================
        function cmd = run_vf_guidance(obj, pos)
            x = pos(1); y = pos(2);
            
            if strcmpi(obj.CurrentMode, 'Loiter')
                % =================================================
                % 【Loiterモード】: 指定された位相角ベースVF則
                % Formula: chi_d = phi + lambda * (pi/2 + atan(k * ~d))
                % =================================================
                if isempty(obj.LoiterParams)
                    cmd=struct('chi_cmd',0,'kappa_ref',0,'error',0); return; 
                end
                p = obj.LoiterParams;
                
                % 1. 幾何学計算
                dx = x - p.xc; dy = y - p.yc;
                dist = sqrt(dx^2 + dy^2);
                
                % 位相角 phi (NED系: 北=0, 東=90)
                phi = atan2(dy, dx);
                
                % 2. 正規化誤差 (~d)
                % 外側(dist>R)でプラス、内側でマイナス
                tilde_d = (dist - p.R) / p.R;
                
                % 3. 回転方向の符号定義
                % 式において「右旋回(CW)」のときに係数がプラスになるように設定
                % NED系では、lambda=1 が右旋回(CW)
                lambda_k = p.lambda; 
                
                % 4. 誘導則計算
                % (pi/2): 接線成分
                % atan(k*d): 収束成分
                k = obj.Gains.k_vf_loiter;
                
                correction_term = (pi/2) + atan(k * tilde_d);
                chi_cmd_val = phi + lambda_k * correction_term;
                
                % 5. 結果格納
                cmd.chi_cmd = obj.normalize_angle(chi_cmd_val);
                cmd.kappa_ref = lambda_k * (1/p.R); % 曲率 (CW正, CCW負)
                cmd.error = dist - p.R;             % 誤差 (m)
                
            else
                % =================================================
                % 【Missionモード】: 点群追従VF (General Path)
                % =================================================
                if isempty(obj.MissionPath)
                    cmd=struct('chi_cmd',0,'kappa_ref',0,'error',0); return; 
                end
                P = obj.MissionPath;
                
                % 1. ウィンドウ探索 (Window Search)
                % 前回の探索点(LastIndex)から、少し先(WindowSize)までを探す。
                % これにより計算を高速化し、かつコース交差時の誤認を防ぐ。
                idx_start = obj.LastIndex;
                idx_end = min(idx_start + obj.WindowSize, P.num_points);
                range = idx_start:idx_end;
                
                % 距離の2乗計算
                d_sq = (P.x(range) - x).^2 + (P.y(range) - y).^2;
                [~, min_loc] = min(d_sq);
                
                idx = range(min_loc);
                obj.LastIndex = idx; % 次回のために更新
                
                % 2. 参照値の取得
                chi_ref = P.chi(idx);
                kappa_ref = P.kappa(idx);
                
                % 3. クロストラックエラー (e) の計算
                % 参照点での接線ベクトルに対し、機体が右にいるか左にいるか？
                % 定義: 進行方向に向かって右ズレを正とする
                dx = x - P.x(idx); 
                dy = y - P.y(idx);
                
                % 外積 (Cross Product) in 2D
                e_cross = -sin(chi_ref)*dx + cos(chi_ref)*dy;
                
                % 4. 誘導則計算 (General Path VF)
                % chi_cmd = chi_ref - chi_inf * (2/pi) * atan(k * e)
                % 右ズレ(e>0)なら、左へ(-方向)戻る必要があるためマイナス
                
                k = obj.Gains.k_vf_mission;
                chi_inf = obj.Gains.chi_inf;
                
                correction = chi_inf * (2/pi) * atan(k * e_cross);
                
                % 5. 結果格納
                cmd.chi_cmd = obj.normalize_angle(chi_ref - correction);
                cmd.kappa_ref = kappa_ref;
                cmd.error = e_cross;
            end
        end
        
        % --- ユーティリティ関数 ---
        
        function [dR, dL] = apply_mixing(~, da, ds)
            % 操舵ミキシングロジック
            % da (差動成分): 正なら右ターン(右引き)、負なら左ターン(左引き)
            % ds (同相成分): ブレーキ(両方引く)
            
            da = max(-1, min(1, da)); % -1~1に制限
            
            if da > 0
                dR = da + ds; % 右を引く
                dL = ds;
            else
                dR = ds;
                dL = abs(da) + ds; % 左を引く
            end
            
            % 物理的限界 (0=全開, 1=全閉)
            dR = max(0, min(1, dR)); 
            dL = max(0, min(1, dL));
        end
        
        function a = normalize_angle(~, a)
            % 角度を -pi ~ pi の範囲に正規化
            a = mod(a + pi, 2*pi) - pi;
        end
    end
end

% ヘルパー関数 (ファイル外またはメソッド外)
function s = calc_rot_str(lam)
    if lam > 0, s='CW (Right)'; else, s='CCW (Left)'; end
end