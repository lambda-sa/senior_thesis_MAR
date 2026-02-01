classdef WindTrajectoryComparatorZoomed < WindTrajectoryComparator
    % WINDTRAJECTORYCOMPARATORZOOMED
    % シミュレーションの最後尾から指定秒数(duration)分だけデータを切り取り、
    % 拡大して比較描画するためのクラス。
    
    properties
        ZoomDuration % タイトル表示用に秒数を保存するプロパティ
    end
    
    methods
        function obj = WindTrajectoryComparatorZoomed(simNoWind, simWind, trajPlanGnd, trajPlanAir, targetPos, planTimes, duration)
            % duration: 切り出す秒数 (省略時は70秒)
            if nargin < 7, duration = 70; end
            
            % 1. 終了時刻の特定 (両ケースの最大値)
            t_end = max(simNoWind.Time(end), simWind.Time(end));
            t_start = max(0, t_end - duration);
            
            fprintf('Zooming in on the last %d seconds (Time: %.1f s to %.1f s)...\n', duration, t_start, t_end);
            
            % 2. データのスライス処理
            simNoWind_sliced = sliceTable(simNoWind, t_start);
            simWind_sliced   = sliceTable(simWind, t_start);
            
            trajPlanGnd_sliced = sliceStructOrTable(trajPlanGnd, t_start);
            trajPlanAir_sliced = sliceStructOrTable(trajPlanAir, t_start);
            
            % 3. フェーズ遷移時刻のフィルタリング
            if ~isempty(planTimes)
                planTimes_sliced = planTimes(planTimes >= t_start & planTimes <= t_end);
            else
                planTimes_sliced = [];
            end
            
            % 4. 親クラスのコンストラクタ呼び出し
            obj@WindTrajectoryComparator(simNoWind_sliced, simWind_sliced, ...
                                         trajPlanGnd_sliced, trajPlanAir_sliced, ...
                                         targetPos, planTimes_sliced);
            
            % 5. 秒数をプロパティに保存 (あとで plotAll で使うため)
            obj.ZoomDuration = duration;
        end
        
        function plotAll(obj)
            % 1. 親クラス(WindTrajectoryComparator)の描画を実行
            %    ここで Figure が生成され、Handles.Figure が有効になります
            plotAll@WindTrajectoryComparator(obj);
            
            % 2. Figureの名前をここで上書き変更
            if isfield(obj.Handles, 'Figure') && isvalid(obj.Handles.Figure)
                set(obj.Handles.Figure, 'Name', sprintf('Comparison: Last %d sec (Zoomed)', obj.ZoomDuration));
            end
        end
    end
end

%% --- Local Helper Functions for Slicing ---

function T_out = sliceTable(T_in, t_start)
    % Table型のデータを Time 列に基づいてスライス
    if isempty(T_in), T_out = T_in; return; end
    T_out = T_in(T_in.Time >= t_start, :);
end

function S_out = sliceStructOrTable(S_in, t_start)
    % StructまたはTable型の軌道データをスライス
    % (Timeフィールドを持つことを前提とする)
    if isempty(S_in)
        S_out = S_in; 
        return; 
    end
    
    if istable(S_in)
        % Tableの場合
        S_out = sliceTable(S_in, t_start);
    elseif isstruct(S_in)
        % Structの場合 (各フィールドを行方向にスライス)
        if ~isfield(S_in, 'Time')
            warning('Structure does not have "Time" field. Returning original.');
            S_out = S_in;
            return;
        end
        
        % フィルタ用インデックス作成
        mask = S_in.Time >= t_start;
        
        fields = fieldnames(S_in);
        S_out = S_in;
        for i = 1:length(fields)
            val = S_in.(fields{i});
            % 行数がTimeと一致するものだけスライス対象とする
            if size(val, 1) == length(mask)
                S_out.(fields{i}) = val(mask, :);
            end
        end
    else
        % その他
        S_out = S_in;
    end
end