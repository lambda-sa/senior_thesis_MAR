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
           if isfield(params.guidance, 'nominal_bank_deg')
                % ガイダンス設定に「標準バンク角」があればそれを使う
                phi_plan = deg2rad(params.guidance.nominal_bank_deg);
            elseif isfield(params.guidance, 'max_bank_deg')
                % なければ最大バンク角の8割程度などを目安にする（あるいはそのまま使う）
                % ここでは単純化のため max_bank_deg をそのまま使用
                phi_plan = deg2rad(params.guidance.max_bank_deg);
            else
                phi_plan = deg2rad(5); % デフォルト値
            end
            
            % R = V_TAS^2 / (g * tan(phi))
            R_calc = V_trim_TAS^2 / (g * tan(phi_plan));
            
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
            pp.min_turn_radius = R_calc;
            pp.alpha_trim = alpha_trim; 
            
            obj.PhysicsParams = pp;
            
            fprintf('Physics Parameters Computed:\n');
            fprintf('  Start Alt: %.0fm\n', h_start);
            fprintf('  V_TAS: %.2f m/s (at Start), V_EAS: %.2f m/s\n', V_trim_TAS, V_EAS);
            fprintf('  L/D: %.2f, Turn Radius: %.1f m\n', gr, R_calc);
        end

        function run_simple_simulation(obj, type, duration, bank_angle_input)
            % RUN_SIMPLE_SIMULATION
            % Excelパラメータからトリム計算を行い、単純軌道を生成する
            % type: 'Straight' or 'Turn'
            
            % 1. Excelからパラメータ読み込み & 初期高度設定
            % (load_params_from_excel は compute_physics_parameters 内で呼ばれる)
            obj.compute_physics_parameters();
            
            pp = obj.PhysicsParams;
            
            % Start Positionの取得 (ExcelのSimSettingsから)
            % x, y, h, psi
            start_vec = pp.start_pos; 
            
            % 2. 制御入力の決定
            if strcmpi(type, 'Straight')
                bank_deg = 0;
            else
                bank_deg = bank_angle_input;
            end
            
            % 3. トリム計算結果の確認
            % compute_physics_parameters で計算された V_air (EAS) と glide_ratio を使う
            fprintf('Using Trim Conditions from Excel:\n');
            fprintf('  Mass: %.2f kg, Area: %.2f m^2\n', obj.Params.m_total, obj.Params.S_c);
            fprintf('  Trim V_EAS: %.2f m/s\n', pp.V_air);
            fprintf('  Trim L/D:   %.2f\n', pp.glide_ratio);
            fprintf('  Start Alt:  %.1f m\n', start_vec(3));
            
            % 4. Plannerへ委譲
            obj.Planner.plan_simple_physics_based(...
                type, ...
                duration, ...
                start_vec, ...
                pp.V_air, ...       % トリムEAS
                pp.glide_ratio, ... % トリムL/D
                bank_deg ...        % バンク角
            );
            
            % 5. ログデータの生成 (Air/Ground)
            % Planner.WindVector がセットされていれば風も適用される
            obj.Planner.save_naive_path();
            
            % Missionクラスが持つ風情報を使う
            if isa(obj.Planner, 'ParafoilPathPlannerWind')
                 wx = obj.Planner.WindVector(1);
                 wy = obj.Planner.WindVector(2);
                 obj.Planner.apply_wind_effect(wx, wy);
            end
            
            % 詳細テーブル計算 (6DOF検証用に必須)
            obj.compute_dual_trajectories();
        end
    end
end