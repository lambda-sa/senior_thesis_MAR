% load_params_from_excel.m
%
% Excelファイルからパラメータを読み込み、
% 'params' 構造体と 'sim_settings' 構造体を構築する。
% (★ $\delta_a$ 関連の空力係数・慣性・制御リミットを読み込むよう改良)

function [params, sim_settings] = load_params_from_excel(excelFileName)
    disp('--- Excelファイルからパラメータを読み込み中 ...');
    
    % --- 1. シートの読み込み ---
    try
        tbl_aircraft = readtable(excelFileName, 'Sheet', 'Aircraft_Model');
        tbl_payload = readtable(excelFileName, 'Sheet', 'Payload_Model');
        tbl_aero = readtable(excelFileName, 'Sheet', 'Aero_Coeffs');
        tbl_sim = readtable(excelFileName, 'Sheet', 'Simulation_Settings');
        tbl_inertia = readtable(excelFileName, 'Sheet', 'Inertia');
        disp('Excelシートの読み込み完了。');
    catch ME
        error('Excelファイルまたはシートの読み込みに失敗しました: %s', ME.message);
    end
    
    % --- 2. ヘルパー関数 (find_val) の定義 ---
    find_val = @(tbl, name) find_val_numeric_helper(tbl, name);
    
    % --- 3. params 構造体の構築 ---
    params = struct();
    
    try
        % (A) 物理定数
        params.m_canopy = find_val(tbl_aircraft, 'M_canopy');
        params.m_payload = find_val(tbl_payload, 'M_payload');
        params.m_total = params.m_canopy + params.m_payload;
        params.g = find_val(tbl_aircraft, 'g');
        params.rho = find_val(tbl_aircraft, 'rho_sea_level');
        
        % (B) 基準面積・基準長
        params.S_c = find_val(tbl_aircraft, 'S_c'); 
        params.S_p = find_val(tbl_payload, 'S_p'); 
        params.c = find_val(tbl_aircraft, 'c');   
        %params.ang_psi_rigging = find_val(tbl_aircraft, 'ang_psi_rigging');
        %params.ang_psi_rigging = deg2rad(params.ang_psi_rigging);
        params.b = find_val(tbl_aircraft, 'b'); % ★追加: 翼幅
        params.d = find_val(tbl_aircraft, 'd'); % ★追加: コントロールライン長

        % --- 慣性モーメント (Aircraft_Model) ---
        params.I_xx = find_val(tbl_inertia, 'I_xx');
        params.I_yy = find_val(tbl_inertia, 'I_yy');
        params.I_zz = find_val(tbl_inertia, 'I_zz');
        params.I_xz = find_val(tbl_inertia, 'I_xz');

        % (C) ペイロード空力
        params.C_D_s = find_val(tbl_payload, 'C_D_s'); 
        
        % (D) 空力係数 (Aero_Coeffs)
        params.C_L_0 = find_val(tbl_aero, 'C_L_0');
        params.C_L_alpha = find_val(tbl_aero, 'C_L_alpha');
        params.C_L_delta_a = find_val(tbl_aero, 'C_L_delta_a');
        params.C_L_delta_s = find_val(tbl_aero, 'C_L_delta_s');
        params.C_D_0 = find_val(tbl_aero, 'C_D_0');
        params.C_D_alpha = find_val(tbl_aero, 'C_D_alpha');
        params.C_D_delta_a = find_val(tbl_aero, 'C_D_delta_a');
        params.C_D_delta_s = find_val(tbl_aero, 'C_D_delta_s');

        params.C_m_0 = find_val(tbl_aero, 'C_m_0');
        params.C_m_alpha = find_val(tbl_aero, 'C_m_alpha');
        params.C_m_q        = find_val(tbl_aero, 'C_m_q');
        %params.C_n_delta_a = find_val(tbl_aero, 'C_n_delta_a')

        params.C_l_phi = find_val(tbl_aero, 'Cl_phi');
        params.C_l_p = find_val(tbl_aero, 'Cl_p');
        params.C_l_delta_a = find_val(tbl_aero, 'Cl_delta_a');

        params.C_n_r = find_val(tbl_aero, 'Cn_r');
        params.C_n_delta_a = find_val(tbl_aero, 'Cn_delta_a');
        % --- 横・方向・回転の空力係数 (Aero_Coeffs) ---
        params.C_Y_beta     = find_val(tbl_aero, 'C_Y_beta');
        

        % (E) モーメントアーム
        Rcg_x = find_val(tbl_aircraft, 'r_total_cm_to_canopy_origin_B_X') + ...
                find_val(tbl_aircraft, 'r_canopy_origin_to_ac_B_X');
        Rcg_z = find_val(tbl_aircraft, 'r_total_cm_to_canopy_origin_B_Z') + ...
                find_val(tbl_aircraft, 'r_canopy_origin_to_ac_B_Z');
        Rpg_x = find_val(tbl_payload, 'delta_x_pay');
        Rpg_z = find_val(tbl_payload, 'delta_z_pay');
        params.Rcg_x = Rcg_x;
        params.Rcg_z = Rcg_z;
        params.Rpg_x = Rpg_x;
        params.Rpg_z = Rpg_z;
        
        % (F) 誘導則パラメータ
        params.guidance = struct();
        params.guidance.turn_start_time = find_val(tbl_sim, 'turn_start_time');
        params.guidance.turn_end_time = find_val(tbl_sim, 'turn_end_time');
        params.guidance.max_bank_deg = find_val(tbl_sim, 'max_bank_deg');
        params.guidance.brake_start_time = find_val(tbl_sim, 'brake_start_time');
        params.guidance.brake_end_time = find_val(tbl_sim, 'brake_end_time');
        
        
        % ★変更: $\dot{\psi}$ (psi_dot) の代わりに $\delta_a$ を読み込む
        params.guidance.target_delta_a = find_val(tbl_sim, 'target_delta_a');
        params.guidance.target_delta_s = find_val(tbl_sim, 'target_delta_s');
        % (G) 時定数
        params.tau_delta_a = find_val(tbl_sim, 'tau_delta_a');
        params.tau_delta_s = find_val(tbl_sim, 'tau_delta_s');
        
        % ★追加: (H) 慣性テンソル (Aero_Coeffs シートから)
        %params.I_XXI = find_val(tbl_aero, 'I_XXI');
        %params.I_XZI = find_val(tbl_aero, 'I_XZI');
        
        disp('params 構造体の構築完了。');
        
    catch ME
        disp('params 構造体の構築中にエラーが発生しました。');
        disp('Excelの Parameter 名が正しいか確認してください。');
        rethrow(ME);
    end
    
    % --- 4. sim_settings 構造体の構築 ---
    try
        sim_settings = struct();
        sim_settings.h_init = find_val(tbl_sim, 'h_init');
        %sim_settings.delta_a_init = find_val(tbl_sim, 'delta_a_init');
        %sim_settings.delta_s_init = find_val(tbl_sim, 'delta_s_init');
        try sim_settings.delta_R_init = find_val(tbl_sim, 'delta_R_init'); catch, sim_settings.delta_R_init = 0; end
        try sim_settings.delta_L_init = find_val(tbl_sim, 'delta_L_init'); catch, sim_settings.delta_L_init = 0; end
        sim_settings.t_max = find_val(tbl_sim, 't_max');
        sim_settings.h_step = find_val(tbl_sim, 'h_step');
        sim_settings.ang_gamma_rigging = deg2rad(find_val(tbl_aircraft, 'ang_psi_rigging'));
        
        % ★追加: 初期状態
        %sim_settings.V_initial = find_val(tbl_sim, 'V_initial');
        %sim_settings.gamma_initial_deg = find_val(tbl_sim, 'gamma_initial_deg');
        sim_settings.psi_initial_deg = find_val(tbl_sim, 'psi_initial_deg');
        sim_settings.X_initial = find_val(tbl_sim, 'X_initial');
        sim_settings.Y_initial = find_val(tbl_sim, 'Y_initial');
        sim_settings.V_init_guess = find_val(tbl_sim, 'V_init_guess');
        %sim_settings.ang_psi_rigging = deg2rad(params.ang_psi_rigging);
        % ★追加: 制御リミット
        params.delta_a_min = find_val(tbl_sim, 'delta_a_min');
        params.delta_a_max = find_val(tbl_sim, 'delta_a_max');
        params.delta_s_min = find_val(tbl_sim, 'delta_s_min');
        params.delta_s_max = find_val(tbl_sim, 'delta_s_max');

        % ★追加: リギング角制御用パラメータ (Simulation_Settings シートなどを想定)
        % 単位変換: deg -> rad
        params.tau_mu     = find_val(tbl_sim, 'tau_mu');         % 時定数
        params.guidance.mu_min     = deg2rad(find_val(tbl_sim, 'mu_min_deg')); 
        params.guidance.mu_max     = deg2rad(find_val(tbl_sim, 'mu_max_deg'));
        params.guidance.target_mu  = deg2rad(find_val(tbl_sim, 'target_mu_deg')); % 目標値

        % 風速
        try sim_settings.wind_x = find_val(tbl_sim, 'wind_x'); catch, sim_settings.wind_x = 0; end
        try sim_settings.wind_y = find_val(tbl_sim, 'wind_y'); catch, sim_settings.wind_y = 0; end
        try sim_settings.wind_z = find_val(tbl_sim, 'wind_z'); catch, sim_settings.wind_z = 0; end
                
    catch ME
        disp('sim_settings 構造体の構築中にエラーが発生しました。');
        disp('Excelの Parameter 名 (h_init, t_max など) が正しいか確認してください。');
        rethrow(ME);
    end
    
    % --- 5. CONTROL_TYPE の読み込み ---
    try
        opts_str = detectImportOptions(excelFileName, 'Sheet', 'Simulation_Settings');
        opts_str = setvartype(opts_str, 'Value', 'string');
        tbl_sim_str = readtable(excelFileName, opts_str);
        str_val = tbl_sim_str.Value(strcmp(tbl_sim_str.Parameter, 'con_type'));
        
        if isempty(str_val)
            error('パラメータ「con_type」が見つかりません。');
        end
        
        params.con_type = char(str_val);
        CONTROL_TYPE = params.con_type;
        
    catch ME
        disp('CONTROL_TYPE の読み込みに失敗しました。');
        rethrow(ME);
    end
    
    fprintf('--- Excelパラメータ読み込み完了 (制御モード: %s) ---\n', CONTROL_TYPE);
end

% --- ヘルパー関数 (数値専用) ---
function val = find_val_numeric_helper(tbl, name)
    try
        val = tbl.Value(strcmp(tbl.Parameter, name));
        if isempty(val)
             error('パラメータ「%s」が見つかりません。', name);
        end
        if isnan(val)
             error('パラメータ「%s」の値がNaNです。Excelファイルを確認してください。', name);
        end
    catch ME
        sheet_name = 'Unknown';
        if exist('tbl', 'var') && ~isempty(tbl.Properties.Description)
            sheet_name = tbl.Properties.Description;
        end
        error('Excelシート "%s" からパラメータ「%s」の検索中にエラー: %s', ...
              sheet_name, name, ME.message);
    end
end