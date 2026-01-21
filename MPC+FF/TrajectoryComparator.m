classdef TrajectoryComparator < handle
    % TRAJECTORYCOMPARATOR (Base Class)
    % 論文用CUD配色・スタイル適用版
    % Green(Plan) / Blue(No Wind) / Red(Wind)
    
    properties
        SimData         % 6DOFシミュレーション結果 (No Wind)
        TrajPlanGnd     % 計画軌道 (Ground)
        TrajPlanAir     % 計画軌道 (Air)
        TargetPos       % 目標位置
        WindPlanned     % 計画風速
        PlanTransitionTimes = [] 
        
        Handles = struct(); 
        
        % ★CUD推奨配色定義 (RGB 0-1正規化)
        Colors = struct( ...
            'Plan',    [3, 175, 122] / 255, ...   % #03AF7A (Ideal: Green)
            'NoWind',  [0, 90, 255] / 255, ...    % #005AFF (No Wind: Blue)
            'Wind',    [255, 75, 0] / 255, ...    % #FF4B00 (Wind: Red)
            'Target',  [0.6, 0.2, 0.2], ...       % Marker (Dark Red)
            'Phase',   [0.4, 0.4, 0.4]  ...       % Phase Lines (Dark Gray)
        );
    end
    
    methods
        function obj = TrajectoryComparator(simData, trajPlanGnd, trajPlanAir, targetPos, windPlanned, planTimes)
            obj.SimData = simData;
            obj.TrajPlanGnd = trajPlanGnd;
            obj.TrajPlanAir = trajPlanAir;
            obj.TargetPos = targetPos;
            obj.WindPlanned = windPlanned;
            if nargin > 5
                obj.PlanTransitionTimes = planTimes;
            end
        end
        
        function plotAll(obj)
            fig = figure('Name', 'Simulation Analysis', 'Color', 'w', ...
                'Position', [50, 50, 1200, 1000]);
            obj.Handles.Figure = fig;
            
            tabGroup = uitabgroup(fig);
            
            t1 = uitab(tabGroup, 'Title', '3D View', 'BackgroundColor', 'w');
            obj.plotTab1_3D(t1);
            
            t2 = uitab(tabGroup, 'Title', 'Analysis (2x2)', 'BackgroundColor', 'w');
            obj.plotTab2_Analysis(t2);
            
            t3 = uitab(tabGroup, 'Title', 'Dynamics', 'BackgroundColor', 'w');
            obj.plotTab3_Dynamics(t3);
        end
    end
    
    methods (Access = protected)
        
        %% --- Tab 1: 3D View ---
        function plotTab1_3D(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            t = tiledlayout(bgPanel, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            ax = nexttile(t);
            obj.Handles.Ax3D = ax; 
            
            simD = obj.SimData; planG = obj.TrajPlanGnd; tgt = obj.TargetPos; C = obj.Colors;
            
            % 計画 (緑・点線)
            plot3(ax, planG.Position(:,1), planG.Position(:,2), planG.Position(:,3), ...
                  ':', 'Color', C.Plan, 'LineWidth', 2.0, 'DisplayName', 'Refference');
            hold(ax, 'on');
            % 無風 (青・実線)
            plot3(ax, simD.Inertial_X_Position, simD.Inertial_Y_Position, -simD.Inertial_Z_Position, ...
                  '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            % ターゲット
            plot3(ax, tgt(1), tgt(2), tgt(3), ...
                  '^', 'Color', C.Target, 'MarkerSize', 8, 'MarkerFaceColor', C.Target, 'DisplayName', 'Target');
            
            view(ax, 3); axis(ax, 'equal');
            xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]'); zlabel(ax, 'Altitude [m]');
            legend(ax, 'show'); obj.applyStyle(ax);
        end
        
        %% --- Tab 2: Analysis (Errors) ---
        function plotTab2_Analysis(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            t = tiledlayout(bgPanel, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            simD = obj.SimData; planG = obj.TrajPlanGnd; tgt = obj.TargetPos; C = obj.Colors;
            [horiz_err, vert_err] = obj.calcErrors(simD);

            % 1. Top View
            ax1 = nexttile(t); obj.Handles.AxTop = ax1;
            plot(ax1, planG.Position(:,1), planG.Position(:,2), ':', 'Color', C.Plan, 'LineWidth', 2.0, 'DisplayName', 'Refference');
            hold(ax1, 'on');
            plot(ax1, simD.Inertial_X_Position, simD.Inertial_Y_Position, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            plot(ax1, tgt(1), tgt(2), '^', 'Color', C.Target, 'MarkerSize', 8, 'MarkerFaceColor', C.Target, 'DisplayName', 'Target');
            axis(ax1, 'equal'); xlabel(ax1, 'x [m]'); ylabel(ax1, 'y [m]');
            legend(ax1, 'show'); obj.applyStyle(ax1);
            
            % 2. Side View
            ax2 = nexttile(t); obj.Handles.AxSide = ax2;
            plot(ax2, planG.Position(:,1), planG.Position(:,3), ':', 'Color', C.Plan, 'LineWidth', 2.0, 'DisplayName', 'Refference');
            hold(ax2, 'on');
            plot(ax2, simD.Inertial_X_Position, -simD.Inertial_Z_Position, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            plot(ax2, tgt(1), tgt(3), '^', 'Color', C.Target, 'MarkerSize', 8, 'MarkerFaceColor', C.Target, 'DisplayName', 'Target');
            xlabel(ax2, 'x [m]'); ylabel(ax2, 'Altitude [m]');
            legend(ax2, 'show'); obj.applyStyle(ax2);
            
            % 3. Horizontal Error
            ax3 = nexttile(t); obj.Handles.AxErrH = ax3;
            plot(ax3, simD.Time, horiz_err, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            xlabel(ax3, 'Time [s]'); ylabel(ax3, 'Horizontal error [m]'); 
            legend(ax3, 'show'); 
            obj.addPhaseLines(ax3); obj.applyStyle(ax3);
            obj.adjustAxesLimits(ax3, horiz_err);
            
            % 4. Vertical Error
            ax4 = nexttile(t); obj.Handles.AxErrV = ax4;
            plot(ax4, simD.Time, vert_err, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            xlabel(ax4, 'Time [s]'); ylabel(ax4, 'Vertical error [m]'); 
            legend(ax4, 'show'); 
            obj.addPhaseLines(ax4); obj.applyStyle(ax4);
            obj.adjustAxesLimits(ax4, vert_err);
        end
        
        %% --- Tab 3: Dynamics ---
        function plotTab3_Dynamics(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            t = tiledlayout(bgPanel, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            simD = obj.SimData; planA = obj.TrajPlanAir; planG = obj.TrajPlanGnd; C = obj.Colors;
            t_plan = planA.Time; t_max = max(t_plan(end), simD.Time(end));
            
            % Bank
            ax1 = nexttile(t); obj.Handles.AxBank = ax1;
            plot(ax1, t_plan, -rad2deg(planA.Euler_RPY(:, 1)), ':', 'Color', C.Plan, 'LineWidth', 2.0, 'DisplayName', 'Refference');
            hold(ax1, 'on');
            plot(ax1, simD.Time, simD.Roll_Angle, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            xlabel(ax1, 'Time [s]'); ylabel(ax1, 'Bank angle [deg]');
            xlim(ax1, [0, t_max]);
            legend(ax1, 'show'); obj.addPhaseLines(ax1); obj.applyStyle(ax1);
            
            % Alt
            ax2 = nexttile(t); obj.Handles.AxAlt = ax2;
            plot(ax2, t_plan, planG.Position(:, 3), ':', 'Color', C.Plan, 'LineWidth', 2.0, 'DisplayName', 'Refference');
            hold(ax2, 'on');
            plot(ax2, simD.Time, -simD.Inertial_Z_Position, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', 'Without wind');
            xlabel(ax2, 'Time [s]'); ylabel(ax2, 'Altitude [m]');
            xlim(ax2, [0, t_max]);
            legend(ax2, 'show'); obj.addPhaseLines(ax2); obj.applyStyle(ax2);
            
            % Control
            ax3 = nexttile(t); obj.Handles.AxCtrl = ax3;
            obj.plotControlInputs(ax3);
            xlim(ax3, [0, t_max]);
            obj.addPhaseLines(ax3); obj.applyStyle(ax3);
        end
        
        % ★制御入力描画 (Base: No Wind Only)
        % 定義: δa = δL - δR (非対称), δs = min(δL, δR) (対称)
        function plotControlInputs(obj, ax)
            simD = obj.SimData; C = obj.Colors;
            
            delta_a = simD.delta_L - simD.delta_R;
            delta_s = min(simD.delta_L, simD.delta_R);
            
            % δa: 実線
            plot(ax, simD.Time, delta_a, '-', 'Color', C.NoWind, 'LineWidth', 1.5, 'DisplayName', '\delta_{\ita} (without wind)');
            hold(ax, 'on');
            % δs: 点線 (太め)
            plot(ax, simD.Time, delta_s, ':', 'Color', C.NoWind, 'LineWidth', 2.0, 'DisplayName', '\delta_{\its} (without wind)');
            
            xlabel(ax, 'Time [s]'); ylabel(ax, 'Control input [m]');
            legend(ax, 'show');
        end

        % --- Helpers ---
        function [horiz_err, vert_err] = calcErrors(obj, simD)
            ref_pts = obj.TrajPlanGnd.Position;
            act_pts = [simD.Inertial_X_Position, simD.Inertial_Y_Position, -simD.Inertial_Z_Position];
            n_act = size(act_pts, 1);
            err_vec = zeros(n_act, 3);
            for i = 1:n_act
                [~, idx] = min(sum((ref_pts - act_pts(i, :)).^2, 2));
                err_vec(i, :) = act_pts(i, :) - ref_pts(idx, :);
            end
            horiz_err = sqrt(err_vec(:,1).^2 + err_vec(:,2).^2);
            vert_err = abs(err_vec(:,3));
        end

        function adjustAxesLimits(~, axes_list, data_array)
            max_val = max(data_array, [], 'all');
            if max_val > 0
                common_ylim = [0, max_val * 1.1];
                for i = 1:length(axes_list)
                    set(axes_list(i), 'YLim', common_ylim);
                end
            end
        end

        % ★xlineでフェーズ線を上下端まで描画
        function addPhaseLines(obj, ax)
            times = obj.PlanTransitionTimes;
            if isempty(times), return; end
            
            for k = 1:length(times)
                t = times(k);
                xl = get(ax, 'XLim');
                if t > xl(1) && t < xl(2)
                    xline(ax, t, ':', 'Color', obj.Colors.Phase, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                end
            end
        end
        
        function applyStyle(~, ax)
            set(ax, 'FontName', 'Arial', 'FontSize', 14); 
            set(ax, 'TickLabelInterpreter', 'tex');
            set(ax.XLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k');
            set(ax.YLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k');
            if isprop(ax, 'ZLabel'), set(ax.ZLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k'); end
            title(ax, ''); grid(ax, 'off');
            set(ax, 'Color', 'w', 'XColor','k', 'YColor','k', 'ZColor','k');
            set(ax, 'LineWidth', 1.0, 'TickDir', 'in', 'TickLength', [0.015 0.025], 'Box', 'on');
            lgd = legend(ax); 
            if ~isempty(lgd), set(lgd, 'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', 12, 'Location', 'best', 'Box', 'off', 'TextColor', 'k'); end
        end
    end
end