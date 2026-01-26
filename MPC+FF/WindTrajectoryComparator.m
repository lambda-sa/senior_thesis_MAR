classdef WindTrajectoryComparator < TrajectoryComparator
    % WINDTRAJECTORYCOMPARATOR (Subclass)
    % 有風データを赤色で重ね書きするクラス
    
    properties
        SimDataWind % 有風シミュレーション結果
        % 鮮やかな黄みの赤 (#FF4B00)
        ColorWind = [255, 75, 0] / 255; 
    end
    
    methods
        function obj = WindTrajectoryComparator(simNoWind, simWind, trajPlanGnd, trajPlanAir, targetPos, planTimes)
            obj@TrajectoryComparator(simNoWind, trajPlanGnd, trajPlanAir, targetPos, [], planTimes);
            obj.SimDataWind = simWind;
        end
        
        function plotAll(obj)
            plotAll@TrajectoryComparator(obj);
            set(obj.Handles.Figure, 'Name', 'Comparison: No Wind vs Wind');
            obj.overlayWindData();
        end
    end
    
    methods (Access = protected)
        % ★制御入力比較 (Override)
        function plotControlInputs(obj, ax)
            d1 = obj.SimData;     % No Wind
            d2 = obj.SimDataWind; % Wind
            C = obj.Colors;
            CW = obj.ColorWind;
            
            % 定義: δa = δL - δR, δs = min(δL, δR)
            da1 = d1.delta_L - d1.delta_R; ds1 = min(d1.delta_L, d1.delta_R);
            da2 = d2.delta_L - d2.delta_R; ds2 = min(d2.delta_L, d2.delta_R);
            
            % No Wind (Blue)
            plot(ax, d1.Time, da1, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', '\delta_{\ita} (without wind)');
            hold(ax, 'on');
            plot(ax, d1.Time, ds1, ':', 'Color', C.NoWind, 'LineWidth', 2.0, 'DisplayName', '\delta_{\its} (without wind)');
            
            % Wind (Red)
            plot(ax, d2.Time, da2, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', '\delta_{\ita} (with wind)');
            plot(ax, d2.Time, ds2, ':', 'Color', CW, 'LineWidth', 2.0, 'DisplayName', '\delta_{\its} (with wind)');
            
            xlabel(ax, 'Time [s]'); ylabel(ax, 'Control input [m]');
            legend(ax, 'show', 'Location', 'bestoutside');
        end
        
        function overlayWindData(obj)
            d1 = obj.SimData; d2 = obj.SimDataWind;
            CW = obj.ColorWind; H = obj.Handles;
            
            % 3D View
            ax = H.Ax3D; hold(ax, 'on');
            plot3(ax, d2.Inertial_X_Position, d2.Inertial_Y_Position, -d2.Inertial_Z_Position, ...
                '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            % Top View
            ax = H.AxTop; hold(ax, 'on');
            plot(ax, d2.Inertial_X_Position, d2.Inertial_Y_Position, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            % Side View
            ax = H.AxSide; hold(ax, 'on');
            plot(ax, d2.Inertial_X_Position, -d2.Inertial_Z_Position, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            % Error (Recalculate and Rescale)
            [he1, ve1] = obj.calcErrors(d1);
            [he2, ve2] = obj.calcErrors(d2);
            
            axH = H.AxErrH; hold(axH, 'on');
            plot(axH, d2.Time, he2, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            axV = H.AxErrV; hold(axV, 'on');
            plot(axV, d2.Time, ve2, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            % Scale unification using all data
            obj.adjustAxesLimits([axH, axV], [he1; ve1; he2; ve2]);
            
            % Dynamics
            ax = H.AxBank; hold(ax, 'on');
            plot(ax, d2.Time, d2.Roll_Angle, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
            
            ax = H.AxAlt; hold(ax, 'on');
            plot(ax, d2.Time, -d2.Inertial_Z_Position, '-', 'Color', CW, 'LineWidth', 1.5, 'DisplayName', 'With wind');
        end
    end
end