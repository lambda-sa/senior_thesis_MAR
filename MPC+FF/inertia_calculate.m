% --- main.m ---
%clearvars
clear all
close all

% --- パラメータをExcelファイルから読み込む ---
excelFileName = 'parameter_inertia_ref.xlsx';
try
    params_aircraft = readtable(excelFileName, 'Sheet', '機体パラメータ');
    params_payload = readtable(excelFileName, 'Sheet', 'ペイロード');
    params_aero_coeffs = readtable(excelFileName, 'Sheet', '空力係数');
    initial_conditions_table = readtable(excelFileName, 'Sheet', '初期条件');
    mission_params_table = readtable(excelFileName, 'Sheet', '誘導目標');

    disp('Excelからパラメータを読み込みました。');
    
catch ME
    error('Excelファイルの読み込みに失敗しました。ファイル名とシート名を確認してください: %s', ME.message);
end

% 全パラメータを1つの構造体に統合
params = struct();
params = add_to_struct(params, params_aircraft);
params = add_to_struct(params, params_payload);
params = add_to_struct(params, params_aero_coeffs);

% 初期条件を別の構造体として保持
initial_conditions_struct = struct();
initial_conditions_struct = add_to_struct(initial_conditions_struct, initial_conditions_table);
mission_params = struct();
mission_params = add_to_struct(mission_params, mission_params_table);

% 新しいクラス群を使って慣性モーメントと物理プロパティを計算
fprintf("\n--- 物理プロパティと慣性モーメントの計算 ---\n");
% 1. ParafoilPropertiesオブジェクトを作成（質量、寸法、重心などが計算される）
prop = ParafoilProperties(params);
prop.visualize_geometry(); % ★★★ この行を追加 ★★★

% 2. ParafoilInertiaCalculatorオブジェクトにpropを渡して慣性モーメントを計算
calc = ParafoilInertiaCalculator(prop);

% --- ヘルパー関数 ---
function combined_struct = add_to_struct(target_struct, source_table)
    for k = 1:height(source_table)
        param_name = source_table.Parameter{k};
        param_value = source_table.Value(k);
        if ~isempty(param_name) && ischar(param_name) && isvarname(param_name)
            target_struct.(param_name) = param_value;
        else, warning('無効なパラメータ名 "%s" はスキップされます。', param_name);
        end
    end
    combined_struct = target_struct;
end