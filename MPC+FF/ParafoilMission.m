classdef ParafoilMission < handle
    % PARAFOILMISSION (決定版: sim_settings対応 & TAS旋回半径)
    
    properties
        ExcelFileName
        Params        % 機体・空力パラメータ (固定値)
        SimSettings   % シミュレーション設定 (初期位置など)
        AtmoModel     % 大気モデル
        Planner       % 経路計画クラス
        PhysicsParams % 計算された物理パラメータ
    end
    
    methods
        function obj = ParafoilMission(excelFileName)
            obj.ExcelFileName = excelFileName;
            
            % 大気モデルの初期化
            try
                obj.AtmoModel = AtmoTempPressRho();
            catch
                warning('AtmoTempPressRho initialization failed. Using Standard Atmosphere.');
                obj.AtmoModel = AtmoModel_Standard(); 
            end
            
            obj.Planner = ParafoilPathPlanner(obj.AtmoModel);
        end
        
        function compute_physics_parameters(obj)
            % 1. パラメータ読み込み (params と sim_settings を分離)
            [params, sim_settings] = load_params_from_excel(obj.ExcelFileName);
            obj.Params = params;
            obj.SimSettings = sim_settings;
            
            % 2. 初期高度の取得
            if isfield(sim_settings, 'h_init')
                h_start = sim_settings.h_init;
            elseif isfield(sim_settings, 'z_initial')
                h_start = -sim_settings.z_initial; % NED座標系
            else
                h_start = 1000; 
            end
            
            % 3. ダイナミクスモデルでTAS計算
            dyn = ParafoilDynamics3DOF_Aero_DeltaA_yaw(params, obj.AtmoModel);
            
            % Rigging角の取得 (sim_settings優先)
            if isfield(sim_settings, 'ang_gamma_rigging')
                mu_rigging = sim_settings.ang_gamma_rigging;
            elseif isfield(params, 'ang_psi_rigging')
                mu_rigging = params.ang_psi_rigging; 
            else
                mu_rigging = 0;
            end
            
            try
                % 初期高度での釣り合い計算
                [alpha_trim, gamma_trim, V_trim_TAS, ~, ~] = ...
                    dyn.calculate_initial_trim(h_start, 0, 0, mu_rigging);
            catch
                warning('Trim calculation failed. Using defaults.');
                V_trim_TAS = 15.0; gamma_trim = deg2rad(-15); alpha_trim = 0;
            end
            
            % 4. 速度変換 (プランナー入力用 EAS)
            rho_start = obj.AtmoModel.get_density(h_start/1000);
            rho_0 = 1.225;
            V_EAS = V_trim_TAS * sqrt(rho_start / rho_0);
            
            % 滑空比
            gr = -1.0 / tan(gamma_trim);
            
            % 5. 旋回半径 R の計算 (TAS基準)
            % 高高度ではTASが速いため、EAS一定計算よりも旋回半径は大きくなる
            g = obj.AtmoModel.get_gravity(h_start);
            if isfield(params.guidance, 'max_bank_deg')
                phi_max = deg2rad(params.guidance.max_bank_deg);
            else
                phi_max = deg2rad(30);
            end
            
            % R = V_TAS^2 / (g * tan(phi))
            R_min = V_trim_TAS^2 / (g * tan(phi_max));
            
            % 6. 開始位置等の設定
            if isfield(sim_settings, 'X_initial'), sx = sim_settings.X_initial; else, sx = -500; end
            if isfield(sim_settings, 'Y_initial'), sy = sim_settings.Y_initial; else, sy = 0; end
            if isfield(sim_settings, 'psi_initial_deg'), syaw = deg2rad(sim_settings.psi_initial_deg); else, syaw = 0; end
            if isfield(sim_settings, 'landing_direction'), land_dir = sim_settings.landing_direction; else, land_dir = 180; end

            pp.start_pos = [sx, sy, h_start, syaw];
            pp.landing_direction = land_dir;
            
            pp.V_air = V_EAS;          
            pp.V_trim_actual = V_trim_TAS; 
            pp.glide_ratio = gr;
            pp.min_turn_radius = R_min;
            pp.alpha_trim = alpha_trim; 
            
            obj.PhysicsParams = pp;
            
            fprintf('Physics Parameters Computed:\n');
            fprintf('  Start Alt: %.0fm\n', h_start);
            fprintf('  V_TAS: %.2f m/s (at Start), V_EAS: %.2f m/s\n', V_trim_TAS, V_EAS);
            fprintf('  L/D: %.2f, Turn Radius: %.1f m\n', gr, R_min);
        end
    end
end