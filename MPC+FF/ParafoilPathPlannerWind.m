classdef ParafoilPathPlannerWind < ParafoilPathPlanner
    % PARAFOILPATHPLANNERWIND (決定版: 風適用 & 可視化)
    
    properties
        ResultDataWind   % 風適用後の対地経路データ
        ResultDataNaive  % 補正なしデータ
    end
    
    methods
        function obj = ParafoilPathPlannerWind(atmo_model)
            obj@ParafoilPathPlanner(atmo_model);
            obj.WindVector = [0; 0; 0];
        end
        
        function save_naive_path(obj)
            if ~isempty(obj.ResultData)
                obj.ResultDataNaive = obj.ResultData;
            end
        end
        
        function apply_wind_effect(obj, wx, wy)
            % 親クラスの WindVector を更新
            obj.WindVector = [wx; wy; 0];
            
            if isempty(obj.ResultData), return; end
            d = obj.ResultData;
            w_data = d; 
            
            % 各セグメントに風の変位(drift = W * t)を加算
            if ~isempty(d.loiter.t)
                w_data.loiter.x = d.loiter.x + wx * d.loiter.t;
                w_data.loiter.y = d.loiter.y + wy * d.loiter.t;
            end
            if ~isempty(d.dubins.t)
                w_data.dubins.x = d.dubins.x + wx * d.dubins.t;
                w_data.dubins.y = d.dubins.y + wy * d.dubins.t;
            end
            if ~isempty(d.final.t)
                w_data.final.x = d.final.x + wx * d.final.t;
                w_data.final.y = d.final.y + wy * d.final.t;
            end
            
            obj.ResultDataWind = w_data;
            
            if ~isempty(w_data.final.x)
                fx = w_data.final.x(end); fy = w_data.final.y(end);
                fprintf('  -> 風適用後の着地点: [%.1f, %.1f] (Target: 0,0)\n', fx, fy);
            end
        end
        
        function plot_wind_comparison(obj)
            obj.plot_comparison();
        end

        function plot_comparison(obj)
            if isempty(obj.ResultDataWind), warning('風データがありません。'); return; end
            
            figure('Color', 'w', 'Name', 'Wind Effect Comparison');
            hold on; grid on; axis equal; view(3);
            xlabel('X [m]'); ylabel('Y [m]'); zlabel('Altitude [m]');
            
            % Naive Path
            if ~isempty(obj.ResultDataNaive)
                n = obj.ResultDataNaive;
                plot3(n.loiter.x, n.loiter.y, n.loiter.z, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'DisplayName', 'No Wind (Baseline)');
                if ~isempty(n.dubins.x), plot3(n.dubins.x, n.dubins.y, n.dubins.z, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility','off'); end
                if ~isempty(n.final.x), plot3(n.final.x, n.final.y, n.final.z, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility','off'); end
            end
            
            % Wind Path
            w = obj.ResultDataWind;
            if ~isempty(w.loiter.x)
                plot3(w.loiter.x, w.loiter.y, w.loiter.z, 'c-', 'LineWidth', 1.5, 'DisplayName', 'Wind: Loiter');
                plot3(w.loiter.x(1), w.loiter.y(1), w.loiter.z(1), 'bs', 'MarkerFaceColor','b', 'HandleVisibility','off');
            end
            if ~isempty(w.dubins.x)
                plot3(w.dubins.x, w.dubins.y, w.dubins.z, 'b-', 'LineWidth', 2.0, 'DisplayName', 'Wind: Dubins');
            end
            if ~isempty(w.final.x)
                plot3(w.final.x, w.final.y, w.final.z, 'r-', 'LineWidth', 3.0, 'DisplayName', 'Wind: Final');
                plot3(w.final.x(end), w.final.y(end), w.final.z(end), 'ro', 'MarkerFaceColor','r', 'HandleVisibility','off');
            end
            
            plot3(0, 0, 0, 'gx', 'MarkerSize', 15, 'LineWidth', 2, 'DisplayName', 'Target (0,0)');
            legend('Location', 'bestoutside');
            
            wx = obj.WindVector(1); wy = obj.WindVector(2);
            title(sprintf('Wind Flight Path\nWind: [%.1f, %.1f] m/s', wx, wy));
            view([-45, 30]);
        end
    end
end