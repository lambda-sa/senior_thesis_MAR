%% inspect_equations.m
% ワークスペースにある A_sym, B_sym を解析して表示します
% 事前に generate_parafoil_model_softmin.m を実行しておいてください

clc;
fprintf('=== Linearized Model Inspector ===\n');

if ~exist('A_sym', 'var') || ~exist('B_sym', 'var')
    error('ワークスペースに A_sym, B_sym が見つかりません。先に生成スクリプトを実行してください。');
end

while true
    fprintf('\n----------------------------------------\n');
    fprintf('どの行列を見ますか？ (終了するには q を入力)\n');
    mat_choice = input('A または B を入力 [A]: ', 's');
    
    if isempty(mat_choice), mat_choice = 'A'; end
    if strcmpi(mat_choice, 'q'), break; end
    
    if strcmpi(mat_choice, 'A')
        TargetMat = A_sym;
        str_mat = 'A';
    elseif strcmpi(mat_choice, 'B')
        TargetMat = B_sym;
        str_mat = 'B';
    else
        fprintf('無効な入力です。\n');
        continue;
    end
    
    % 行と列の入力
    r = input(sprintf('%s行列の 行番号 (Row) : ', str_mat));
    c = input(sprintf('%s行列の 列番号 (Col) : ', str_mat));
    
    if isempty(r) || isempty(c) || r<=0 || c<=0
        fprintf('インデックスが不正です。\n');
        continue;
    end
    
    try
        fprintf('\n計算中... (数式を簡略化しています)\n');
        
        % 1. 対象の要素を抽出
        element = TargetMat(r, c);
        
        % 2. 簡略化 (これに時間がかかる場合があります)
        % simplifyStepsを増やすとより綺麗になりますが時間がかかります
        element_simple = simplify(element, 'Steps', 50);
        
        if element_simple == 0
            fprintf('\n結果: 0 (この項は影響しません)\n');
        else
            % 3. コンソールに見やすく表示 (ASCIIアート形式)
            fprintf('\n--- MATLAB Pretty Print ---\n');
            pretty(element_simple);
            
            % 4. LaTeX形式で表示 (論文・レポート用)
            fprintf('\n--- LaTeX Code ---\n');
            latex_code = latex(element_simple);
            disp(latex_code);
        end
        
    catch ME
        fprintf('エラーが発生しました: %s\n', ME.message);
    end
end
fprintf('終了します。\n');