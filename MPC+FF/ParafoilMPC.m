classdef ParafoilMPC < handle
    properties
        A_d, B_d
        N, dt
        Q, R
    end
    
    methods
        function obj = ParafoilMPC(A_c, B_c, dt, N)
            obj.dt = dt;
            obj.N = N;
            
            nx = size(A_c, 1);
            
            % オイラー離散化
            obj.A_d = eye(nx) + A_c * dt;
            obj.B_d = B_c * dt;
            
            % 重み行列 (デフォルト)
            % 状態: [u v w p q r phi theta psi x y z]
            % psi(9) と y(11) に重みを置く
            obj.Q = diag([1 1 1, 1 1 1, 10 10 5, 2 2 1]); 
            
            % 入力: [delta_R; delta_L] の変化率抑制
            obj.R = diag([20, 20]); 
        end
        
        % ParafoilMPC.m の solve メソッドを差し替えてください

        function [delta_u_opt, predicted_traj] = solve(obj, x_current, ref_traj, x_trim, u_trim_vec)
            nx = 12;
            nu = 2;
            N_p = obj.N;
            
            % 1. 現在の状態偏差 dx = x_curr - x_trim
            dx = x_current - x_trim;
            % 方位角の偏差は角度差分関数を使う (重要)
            dx(9) = obj.angdiff(x_trim(9), x_current(9)); 
            
            % 2. 参照軌道の偏差 dref = x_ref - x_trim
            % 行列演算用に reshape する前に、角度差分を正しく計算する
            dref = ref_traj - repmat(x_trim, 1, N_p);
            
            % ★追加修正: Referenceの角度変化も angdiff で計算し直す
            for i = 1:N_p
                % ref_traj(9, i) と x_trim(9) の差分を計算して上書き
                dref(9, i) = obj.angdiff(x_trim(9), ref_traj(9, i));
            end
            
            % 3. QP行列の構築
            [H, f] = obj.construct_qp_matrices(dx, dref);
            
            % 4. 制約条件 (u_min <= u <= u_max)
            lb = zeros(nu * N_p, 1);
            ub = zeros(nu * N_p, 1);
            
            for i = 1:N_p
                idx = (i-1)*nu + (1:nu);
                lb(idx) = -u_trim_vec;       % Delta u >= -u_trim
                ub(idx) = 1.0 - u_trim_vec;  % Delta u <= 1 - u_trim
            end
            
            % 5. QPソルバー実行
            options = optimoptions('quadprog','Display','off');
            dU_seq = quadprog(H, f, [], [], [], [], lb, ub, [], options);
            
            if isempty(dU_seq)
                delta_u_opt = [0; 0];
                predicted_traj = [];
            else
                delta_u_opt = dU_seq(1:2);
                
                % デバッグ用予測軌道
                predicted_traj = zeros(nx, N_p);
                xt = dx;
                for k=1:N_p
                    ut = dU_seq((k-1)*nu+(1:nu));
                    xt = obj.A_d * xt + obj.B_d * ut;
                    predicted_traj(:,k) = xt + x_trim; 
                    % 注意: 線形近似上の予測なので、大きな角度変化時には実際の軌道とズレます
                end
            end
        end
        
        function [H, f] = construct_qp_matrices(obj, dx0, dref)
            nx = size(obj.A_d, 1);
            nu = size(obj.B_d, 2);
            N_p = obj.N;
            
            Q_bar = kron(eye(N_p), obj.Q);
            R_bar = kron(eye(N_p), obj.R);
            
            Sx = zeros(nx * N_p, nx);
            Su = zeros(nx * N_p, nu * N_p);
            
            A_pow = eye(nx);
            for i = 1:N_p
                A_pow = obj.A_d * A_pow;
                Sx((i-1)*nx+1 : i*nx, :) = A_pow;
                for j = 1:i
                    if i==j, term = obj.B_d;
                    else, term = (obj.A_d^(i-j)) * obj.B_d; end
                    Su((i-1)*nx+1 : i*nx, (j-1)*nu+1 : j*nu) = term;
                end
            end
            
            ref_vec = reshape(dref, [], 1);
            d_err = Sx * dx0 - ref_vec;
            
            H = 2 * (Su' * Q_bar * Su + R_bar);
            H = (H + H') / 2;
            f = 2 * Su' * Q_bar * d_err;
        end
        
        function d = angdiff(~, a, b)
            d = b - a;
            d = mod(d + pi, 2*pi) - pi;
        end
    end
end