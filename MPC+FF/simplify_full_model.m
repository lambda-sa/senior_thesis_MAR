%% simplify_full_model.m
% 目的: A_sym, B_sym の全要素に対して simplify を適用し、数式を最適化する
% 注意: 計算には数分〜数十分かかる場合があります。
clc;
fprintf('=== Full Model Simplification (A & B) ===\n');

if ~exist('A_sym', 'var') || ~exist('B_sym', 'var')
    error('ワークスペースに A_sym, B_sym が見つかりません。先に生成スクリプトを実行してください。');
end

% 結果格納用の変数を初期化
A_simple = A_sym; % コピーを作成
B_simple = B_sym; % コピーを作成

% タスク総数の計算
total_elems_A = numel(A_sym);
total_elems_B = numel(B_sym);
total_tasks = total_elems_A + total_elems_B;
current_task = 0;

% プログレスバーの表示
h = waitbar(0, 'Initializing simplification...');
tic; % 時間計測開始

%% 1. A行列の Simplify
[rows_A, cols_A] = size(A_sym);
fprintf('Processing A Matrix (%d x %d)...\n', rows_A, cols_A);

for r = 1:rows_A
    for c = 1:cols_A
        current_task = current_task + 1;
        
        % 進捗更新
        progress = current_task / total_tasks;
        msg = sprintf('Simplifying A(%d,%d)... [%d/%d]', r, c, current_task, total_tasks);
        waitbar(progress, h, msg);
        
        % simplify実行 (Stepsオプションで深さを調整可能)
        % 計算が止まるのを防ぐため try-catch は使いませんが、
        % 重すぎる項はそのままにしておく等の処理が必要なら追加可能です。
        A_simple(r, c) = simplify(A_sym(r, c));
    end
end

%% 2. B行列の Simplify
[rows_B, cols_B] = size(B_sym);
fprintf('Processing B Matrix (%d x %d)...\n', rows_B, cols_B);

for r = 1:rows_B
    for c = 1:cols_B
        current_task = current_task + 1;
        
        % 進捗更新
        progress = current_task / total_tasks;
        msg = sprintf('Simplifying B(%d,%d)... [%d/%d]', r, c, current_task, total_tasks);
        waitbar(progress, h, msg);
        
        % simplify実行
        B_simple(r, c) = simplify(B_sym(r, c));
    end
end

% 完了処理
close(h);
total_time = toc;
fprintf('\n=== Simplification Complete ===\n');
fprintf('総所要時間: %.1f 秒\n', total_time);

%% 3. 結果の比較とファイル生成
fprintf('\n--- 最適化結果の確認 ---\n');
% 適当な要素 (例えば A(1,1)) で文字数を比較
len_raw = length(char(A_sym(1,1)));
len_sim = length(char(A_simple(1,1)));
fprintf('例 A(1,1)の文字数: %d -> %d (%.1f%% 削減)\n', ...
    len_raw, len_sim, (1 - len_sim/len_raw)*100);

% ファイル生成
fprintf('\n軽量化された関数ファイル "get_parafoil_jacobians_simple.m" を生成します...\n');

% 数式が綺麗になったので、matlabFunctionの最適化('Optimize', true)を
% ONにしても成功する可能性が高いですが、安全のため今回は false のままで
% 「読みやすいコード」を出力させます。
matlabFunction(A_simple, B_simple, ...
    'File', 'get_parafoil_jacobians_simple', ...
    'Vars', {x_state, u_input, params_sym}, ...
    'Optimize', false); 

fprintf('[SUCCESS] ファイル生成完了: get_parafoil_jacobians_simple.m\n');
fprintf('以前のファイルと入れ替えて使用できます。\n');