classdef VFPathFollowingScheduler_Cascade < handle
    % VFPATHFOLLOWINGSCHEDULER_CASCADE
    % 修正版: 物理モデルベースの誘導・制御則
    
    properties
        % --- Simulation & Environment ---
        dt              % Time step
        
        % --- Path Definition ---
        Waypoints       % [Start_N, Start_E; End_N, End_E]
        OrbitCenter     % [N, E]
        OrbitRadius     % [m]
        OrbitDirection  % +1 (CW) or -1 (CCW)
        d_boundary      % Switching distance [m]
        
        % --- Guidance Gains ---
        LookAheadTime   % [s] LOS Look-ahead time
        
        % --- Control Gains ---
        Tau_yaw         % [s] Yaw tracking time constant (Outer Loop)
        Kp_rate         % P gain (Inner Loop)
        Ki_rate         % I gain (Inner Loop)
        Kd_rate         % D gain (Inner Loop)
        
        % --- Physical Parameters ---
        b               % Span [m]
        Cnr             % Yaw damping derivative
        Cnda            % Control derivative
        
        % --- Internal State ---
        current_mode    % 1: Line, 2: Orbit
        current_idx     % ★追加: 現在のウェイポイントインデックス
        
        integ_rate_err  % 積分器の状態
        prev_rate_err   % 微分用の前回値
        
        % --- Limits ---
        MaxYawRate      % [rad/s]
        Limit_delta_a   % [0-1]
        base_brake      % [0-1]
        
        % --- Logging ---
        LastChi_d
        LastPsi_ref
        LastR_cmd
        LastFF_term
        LastFB_term

        phase=1;
    end
    
    methods
        function obj = VFPathFollowingScheduler_Cascade(dt_init, params)
            % コンストラクタ
            obj.dt = dt_init;
            
            % --- Path Settings (Default) ---
            obj.Waypoints = [0, 0; 800, 800];
            obj.OrbitCenter = [600, 600];
            obj.OrbitRadius = 200;
            obj.OrbitDirection = 1; 
            obj.d_boundary = 50; % 切り替え距離
            
            % --- Guidance Parameters ---
            obj.LookAheadTime = 0.01;
            
            % --- Control Gains ---
            obj.Tau_yaw = 1.0 / 0.5; % Bandwidth ~ 0.8 rad/s
            obj.Kp_rate = 0;
            obj.Ki_rate = 0.07;
            obj.Kd_rate = 0;
            
            % --- Physical Parameters ---
            if isfield(params, 'b'), obj.b = params.b; else, obj.b = 10; end
            if isfield(params, 'Cnr'), obj.Cnr = params.Cnr; else, obj.Cnr = -0.27; end
            if isfield(params, 'Cnda'), obj.Cnda = params.Cnda; else, obj.Cnda = 0.066; end
            
            % Limits
            obj.MaxYawRate = deg2rad(25);
            obj.Limit_delta_a = 0.8;
            obj.base_brake = 0.0;
            
            % Init
            obj.reset();
        end
        
        function reset(obj)
            obj.current_mode = 1;
            obj.current_idx = 1; % ★初期化
            obj.integ_rate_err = 0;
            obj.prev_rate_err = 0;
        end
        
        % ★入力引数に Wind を追加
        function inputs = get_inputs(obj, t, y, dt_sim, wind_vector)
            obj.dt = dt_sim;
            
            % 1. 状態量の抽出
            vel   = y(1:3);   % u, v, w (Body)
            omega = y(4:6);   % p, q, r (Body)
            att   = y(7:9);   % phi, theta, psi
            pos   = y(10:12); % N, E, D
            Va    = norm(vel);
            
            % 風ベクトルが未指定の場合はゼロとする
            if nargin < 5 || isempty(wind_vector)
                wind_vector = [0;0];
            end
            
            % 2. コアロジックの実行
            [delta_a, log_data] = obj.step_logic(pos, vel, att, omega, Va, wind_vector);
            
            % 3. Control Allocation
            inputs.delta_a = delta_a;
            if delta_a >= 0
                inputs.delta_R = obj.base_brake + delta_a;
                inputs.delta_L = obj.base_brake;
            else
                inputs.delta_R = obj.base_brake;
                inputs.delta_L = obj.base_brake + abs(delta_a);
            end
            
            % リミット処理
            inputs.delta_R = max(min(inputs.delta_R, 1), 0);
            inputs.delta_L = max(min(inputs.delta_L, 1), 0);
            inputs.GAMMA = 0;
            inputs.wind_I = [wind_vector(1); wind_vector(2); 0];
            
            % ログ保存
            obj.LastChi_d   = log_data.chi_d;
            obj.LastPsi_ref = log_data.psi_ref;
            obj.LastR_cmd   = log_data.r_cmd;
            obj.LastFF_term = log_data.FF_term;
            obj.LastFB_term = log_data.FB_term;
        end
        
        % ==========================================================
        % ★ Physics-Based Cascade Control Logic
        % ==========================================================
        function [delta_a, log] = step_logic(obj, pos, vel, att, omega, Va, Wind)
            
            x = pos(1);
            y = pos(2);
            psi = att(3);
            
            % --- 対地速度の計算 (簡易版) ---
            % 本来はDCMで計算すべきですが、ここでは近似計算
            % V_ground = V_body_proj + Wind
            x_dot = vel(1)*cos(psi) + Wind(1);
            y_dot = vel(1)*sin(psi) + Wind(2);
            V_ground = norm([x_dot, y_dot]);
            
            % =============================================================
            % Step 1: Guidance (Line-of-Sight + Wind Correction)
            % =============================================================
            
            % --- モード判定 (Orbit切り替え) ---
            if obj.current_mode == 1
                % 直線モード: ウェイポイント管理
                if obj.current_idx <= size(obj.Waypoints, 2)
                    wp_next = obj.Waypoints(1:2, obj.current_idx);
                    
                    dist = norm([x;y] - wp_next);
                    if dist < 20.0 % 到達半径
                        obj.current_idx = obj.current_idx + 1;
                    end
                end
                
                % 全部終わったらOrbitへ
                if obj.current_idx > size(obj.Waypoints, 2)
                    obj.current_mode = 2;
                end
            end
            
            % --- 目標方位角 psi_ref の計算 ---
            if obj.current_mode == 1
                % === 直線追従 (Line Following) ===
                if obj.current_idx > size(obj.Waypoints, 2)
                    % ガード: 万が一インデックスが溢れていたら最後のWPを使う
                    wp_next = obj.Waypoints(1:2, end);
                    wp_prev = obj.Waypoints(1:2, end); 
                else
                    wp_next = obj.Waypoints(1:2, obj.current_idx);
                    if obj.current_idx == 1
                        wp_prev = [0; 0];
                    else
                        wp_prev = obj.Waypoints(1:2, obj.current_idx-1);
                    end
                end
                
                path_vec = wp_next - wp_prev;
                curr_vec = [x; y] - wp_prev;
                
                psi_path = atan2(path_vec(2), path_vec(1));
                e_ct = -(curr_vec(1) * sin(psi_path)) + (curr_vec(2) * cos(psi_path));
                
                % Look-ahead
                Delta = max(V_ground * obj.LookAheadTime, 15.0);
                chi_d = psi_path - atan(e_ct / Delta);
                %chi_d = psi_path ;
                
            else
                % === 円旋回 (Orbit Following) ===
                d_vec = [x; y] - obj.OrbitCenter(:);
                dist_center = norm(d_vec);
                psi_tangent = atan2(d_vec(2), d_vec(1)) + obj.OrbitDirection * (pi/2);
                
                % 距離誤差補正
                k_orbit = 0.4;
                correction = atan(k_orbit * (dist_center - obj.OrbitRadius) / obj.OrbitRadius);
                chi_d = psi_tangent + obj.OrbitDirection * correction * 1.5;
            end
            
            % --- 風補正 (Wind Correction) ---
            W_cross = -Wind(1)*sin(chi_d) + Wind(2)*cos(chi_d);
            ratio = W_cross / max(Va, 0.1);
            if abs(ratio) > 1.0, ratio = sign(ratio); end
            wca = asin(ratio);
            
            psi_ref_target = chi_d - wca;
            
            % 角度の連続化
            err_psi_raw = psi_ref_target - psi;
            err_psi_raw = mod(err_psi_raw + pi, 2*pi) - pi;
            psi_ref = psi + err_psi_raw;
            
            
            % =============================================================
            % Step 2: Kinematics (Heading Error -> Yaw Rate Command)
            % ★ここが抜けていました
            % =============================================================
            
            % 目標ヨーレート r_cmd の計算 (P制御)
            % 時定数 Tau_yaw で追従させる
            r_cmd = (1.0 / obj.Tau_yaw) * err_psi_raw;
            
            % ヨーレートリミッタ
            r_cmd = max(min(r_cmd, obj.MaxYawRate), -obj.MaxYawRate);
            
            
            % =============================================================
            % Step 3: Dynamics (Rate Tracking -> delta_a)
            % =============================================================
            
            r_curr = omega(3);
            Va_safe = max(Va, 5.0);
            
            % Feedforward Term
            FF_term = - (obj.Cnr / obj.Cnda) * (obj.b / (2 * Va_safe)) * r_cmd;
            
            % Feedback Term (PID)
            err_r = r_cmd - r_curr;
            
            % Anti-windup (簡易版: 誤差が大きいときは積分しない)
            if abs(err_r) < 0.2
                obj.integ_rate_err = obj.integ_rate_err + err_r * obj.dt;
            end
            
            d_err_r = (err_r - obj.prev_rate_err) / obj.dt;
            obj.prev_rate_err = err_r;
            
            FB_term = obj.Kp_rate * err_r + ...
                      obj.Ki_rate * obj.integ_rate_err + ...
                      obj.Kd_rate * d_err_r;
            
            % 合計
            delta_a = FF_term + FB_term;
            
            % Control Saturation
            delta_a = max(min(delta_a, obj.Limit_delta_a), -obj.Limit_delta_a);
            
            % --- Log Data Output ---
            log.chi_d = chi_d;
            log.psi_ref = psi_ref;
            log.r_cmd = r_cmd;
            log.FF_term = FF_term;
            log.FB_term = FB_term;
        end
    end
end