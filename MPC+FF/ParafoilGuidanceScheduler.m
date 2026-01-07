classdef ParafoilGuidanceScheduler < handle
    % PARAFOILGUIDANCESCHEDULER
    % 物理モデルベースのカスケード誘導・制御則
    % SimulationEngine対応版: get_inputs メソッドを実装

    properties
        % --- 設定パラメータ ---
        dt              % 制御周期 [s]
        
        % 1. Guidance (Vector Field)
        VF_chi_inf      % 最大進入角 [rad]
        VF_k_line       % 直線収束ゲイン
        VF_k_orbit      % 旋回収束ゲイン
        
        % 2. Wind Rejection (Course Control)
        Kp_course, Ki_course
        
        % 3. Outer Loop (Attitude)
        Tau_yaw         % 時定数
        
        % 4. Inner Loop (Rate Control)
        Kp_rate, Ki_rate, Kd_rate
        
        % --- 物理パラメータ ---
        b, Cnr, Cnda
        
        % --- 内部状態 ---
        integ_course_err
        integ_rate_err
        prev_rate_err
        
        % --- ミッション設定 (SimulationEngineからセットされる) ---
        Waypoints       % [Start_N, Start_E; End_N, End_E]
        OrbitCenter     % [N, E]
        OrbitRadius     % [m]
        OrbitDirection  % +1 or -1
        d_boundary      % 切り替え距離 [m]
        
        % --- 現在のモード ---
        current_mode    % 'line' or 'orbit'
        
        % --- 制約 ---
        max_delta_a
        base_brake

        phase=1;
    end
    
    methods
        function obj = ParafoilGuidanceScheduler(dt, phys_params)
            % コンストラクタ
            obj.dt = dt;
            
            % 物理パラメータ
            obj.b = phys_params.b;
            obj.Cnr = phys_params.Cnr;
            obj.Cnda = phys_params.Cnda;
            
            % --- デフォルトゲイン設定 ---
            obj.VF_chi_inf = deg2rad(80);
            obj.VF_k_line  = 0.02; 
            obj.VF_k_orbit = 0.04;
            
            obj.Kp_course  = 2.0;
            obj.Ki_course  = 0.5; 
            
            obj.Tau_yaw    = 1.0 / 0.5; 
            
            obj.Kp_rate    = 0.1; 
            obj.Ki_rate    = 0.5;
            obj.Kd_rate    = 0.01;
            
            obj.max_delta_a = 0.8;
            obj.base_brake  = 0.0;
            
            obj.current_mode = 'line';
            obj.reset();
        end
        
        function reset(obj)
            obj.integ_course_err = 0;
            obj.integ_rate_err   = 0;
            obj.prev_rate_err    = 0;
            obj.current_mode     = 'line';
        end
        
        % ==========================================================
        % ★ SimulationEngine用のインターフェース (ここが重要)
        % ==========================================================
        function inputs = get_inputs(obj, t, y, wind_local)
            % GET_INPUTS
            % SimulationEngineから各ステップで呼ばれる関数
            
            % 1. 状態ベクトルの展開
            state.vel   = y(1:3);   % u, v, w
            state.omega = y(4:6);   % p, q, r
            state.att   = y(7:9);   % phi, theta, psi
            state.pos   = y(10:12); % N, E, D
            state.Va    = norm(state.vel);
            
            % 2. ミッションロジック (Line -> Orbit 自動切り替え)
            dist_to_center = norm(state.pos(1:2) - obj.OrbitCenter(:));
            
            if strcmp(obj.current_mode, 'line') && (dist_to_center < obj.d_boundary)
                obj.current_mode = 'orbit';
                % fprintf('Time %.2f: Guidance switched to ORBIT mode.\n', t);
            end
            
            % ターゲット構造体の作成
            target.type = obj.current_mode;
            if strcmp(target.type, 'line')
                target.wpt1 = obj.Waypoints(1, :)';
                target.wpt2 = obj.Waypoints(2, :)';
            else
                target.center = obj.OrbitCenter(:);
                target.radius = obj.OrbitRadius;
                target.dir    = obj.OrbitDirection;
            end
            
            % 3. 制御計算の実行 (Core Logic)
            [inputs, ~] = obj.step(state, target, wind_local);
            
            % ※ SimulationEngineによっては inputs に余計なフィールドがあると
            % エラーになる場合があるため、必要なものだけ返す設計にしています。
            % ログが必要な場合は SimulationEngine 側で保存する必要があります。
        end
        
        % ==========================================================
        % コアロジック (Physics-Based Cascade Control)
        % ==========================================================
        function [ctrl_input, log] = step(obj, state, target, wind_inertial)
            %% Step 1: Guidance (Vector Field)
            chi_d = obj.calc_vector_field(state.pos(1:2), target);
            
            %% Step 2: Wind Rejection (Course -> Yaw)
            % 対地速度ベクトルの方向から現在のコース角を計算
            % (簡易版: 水平飛行近似)
            V_N = state.vel(1)*cos(state.att(3)) - state.vel(2)*sin(state.att(3)); 
            V_E = state.vel(1)*sin(state.att(3)) + state.vel(2)*cos(state.att(3));
            chi_curr = atan2(V_E, V_N);
            
            err_chi = angdiff(chi_curr, chi_d);
            obj.integ_course_err = obj.integ_course_err + err_chi * obj.dt;
            
            psi_ref = chi_d + obj.Kp_course * err_chi + obj.Ki_course * obj.integ_course_err;
            
            %% Step 3: Kinematics (Outer Loop)
            err_psi = angdiff(state.att(3), psi_ref);
            psi_dot_req = (1 / obj.Tau_yaw) * err_psi;
            
            % 座標変換
            phi = state.att(1); theta = state.att(2); q = state.omega(2);
            cos_phi_safe = sign(cos(phi)) * max(abs(cos(phi)), 0.1);
            r_cmd = (psi_dot_req * cos(theta) - q * sin(phi)) / cos_phi_safe;
            
            % リミッター
            r_cmd = max(min(r_cmd, deg2rad(30)), -deg2rad(30));
            
            %% Step 4: Dynamics (Inner Loop: FF + PID)
            r_curr = state.omega(3);
            Va_safe = max(state.Va, 5.0);
            
            % ★物理フィードフォワード
            FF_term = - (obj.Cnr / obj.Cnda) * (obj.b / (2 * Va_safe)) * r_cmd;
            
            % PID
            err_r = r_cmd - r_curr;
            obj.integ_rate_err = obj.integ_rate_err + err_r * obj.dt;
            d_err_r = (err_r - obj.prev_rate_err) / obj.dt;
            obj.prev_rate_err = err_r;
            
            FB_term = obj.Kp_rate * err_r + obj.Ki_rate * obj.integ_rate_err + obj.Kd_rate * d_err_r;
            
            delta_a = FF_term + FB_term;
            delta_a = max(min(delta_a, obj.max_delta_a), -obj.max_delta_a);
            
            %% Step 5: Allocation
            if delta_a >= 0
                delta_R = obj.base_brake + delta_a;
                delta_L = obj.base_brake;
            else
                delta_R = obj.base_brake;
                delta_L = obj.base_brake + abs(delta_a);
            end
            
            ctrl_input.delta_R = max(min(delta_R, 1), 0);
            ctrl_input.delta_L = max(min(delta_L, 1), 0);
            ctrl_input.delta_a = delta_a;
            
            % ログ用データ
            log.chi_d = chi_d;
            log.psi_ref = psi_ref;
            log.r_cmd = r_cmd;
        end
        
        function chi_d = calc_vector_field(obj, pos, target)
            px = pos(1); py = pos(2);
            if strcmp(target.type, 'line')
                p1 = target.wpt1; p2 = target.wpt2;
                path_vec = p2 - p1;
                chi_q = atan2(path_vec(2), path_vec(1));
                e = -(px - p1(1))*sin(chi_q) + (py - p1(2))*cos(chi_q);
                chi_d = chi_q - obj.VF_chi_inf * (2/pi) * atan(obj.VF_k_line * e);
            else
                center = target.center; rho = target.radius; lam = target.dir;
                dx = px - center(1); dy = py - center(2);
                phi_ang = atan2(dy, dx);
                d_tilde = sqrt(dx^2 + dy^2) - rho;
                chi_d = phi_ang + lam * (pi/2 + atan(obj.VF_k_orbit * d_tilde / rho));
            end
        end
    end
end

function d = angdiff(x, y)
    d = angle(exp(1i*(y - x)));
end