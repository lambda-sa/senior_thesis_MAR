classdef PathFollowingScheduler < handle
    % PATHFOLLOWINGSCHEDULERLINEARFF 線形化モデルに基づくFF制御
    
    properties
        RefTime
        RefPsi
        RefPsiDot
        
        K_FF       % ★ 線形化により求めたゲイン
        K_P = 0.05; % フィードバックゲイン
        
        WindVector
        phase = 2;
    end
    
    methods
        function obj = PathFollowingScheduler(sim_data_3dof, wind_vec, params, y0)
            % コンストラクタ
            % y0: 6-DOFの初期状態ベクトル (線形化の中心点として使用)
            
            % 1. 経路データの処理
            raw_time = sim_data_3dof.time;
            [unique_time, unique_idx] = unique(raw_time, 'stable');
            obj.RefTime = unique_time;
            
            ref_x = sim_data_3dof.position(unique_idx, 1);
            ref_y = sim_data_3dof.position(unique_idx, 2);
            
            vx = gradient(ref_x, unique_time);
            vy = gradient(ref_y, unique_time);
            obj.RefPsi = atan2(vy, vx); 
            
            psi_unwrap = unwrap(obj.RefPsi);
            obj.RefPsiDot = gradient(psi_unwrap, unique_time);
            
            obj.WindVector = wind_vec(:);
            
            % 2. ★★★ 線形化によるFFゲインの自動計算 ★★★
            % 線形化用の一時的なモデルを作成
            temp_model = ParafoilDynamics(params);
            
            % 外部関数(または下部メソッド)で計算
            obj.K_FF = calculate_linearized_gain(temp_model, y0, params);
        end
        
        function inputs = get_inputs(obj, t, y, ~)
            psi_current = y(9);
            
            t_ref = min(t, obj.RefTime(end));
            target_psi = interp1(obj.RefTime, obj.RefPsi, t_ref, 'linear', 'extrap');
            target_psi_dot = interp1(obj.RefTime, obj.RefPsiDot, t_ref, 'linear', 'extrap');
            
            % FF制御 (線形化ゲインを使用)
            u_ff = obj.K_FF * target_psi_dot;
            
            % FB制御
            err_psi = target_psi - psi_current;
            while err_psi > pi,  err_psi = err_psi - 2*pi; end
            while err_psi < -pi, err_psi = err_psi + 2*pi; end
            
            u_fb = obj.K_P * err_psi;
            
            u_total = u_ff + u_fb;
            
            % リミッタ
            if u_total > 1.0, u_total = 1.0; end
            if u_total < -1.0, u_total = -1.0; end
            
            inputs.delta_a = u_total;
            if u_total >= 0
                inputs.delta_R = u_total; inputs.delta_L = 0;
            else
                inputs.delta_R = 0; inputs.delta_L = -u_total;
            end
            inputs.GAMMA = 0;
            inputs.wind_I = obj.WindVector;
        end
    end
end