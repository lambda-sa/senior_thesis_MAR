classdef TrajectoryComparator
    % TRAJECTORYCOMPARATOR 参照軌道と実行軌道を比較・可視化するクラス
    %
    % [仕様]
    % - Tab 1: Single large 3D plot
    % - Tab 2: 2x2 Grid (Top View, Side View, Horiz Error, Vert Error)
    % - Tab 3: Dynamics (Bank, Alt, Control)
    % - Layout: TiledLayout (Tight spacing)
    % - Colors: Blue (6DoF) & Orange (Ideal)
    % - Phase Lines: Uses "Planned Transition Times"
    % - Style: Arial, TeX, White Background, No Grid/Title
    
    properties
        SimData         % 6DOFシミュレーション結果 (table)
        TrajPlanGnd     % 計画軌道 (Ground: table)
        TrajPlanAir     % 計画軌道 (Air: table - バンク角参照用)
        TargetPos       % 目標位置 [N, E, Alt]
        WindPlanned     % 計画風速 [N, E]
        
        % 計画上のフェーズ遷移時刻リスト
        PlanTransitionTimes = [] 
        
        % カラーパレット
        Colors = struct( ...
            'Blue',   [0.1216, 0.4667, 0.7059], ... % 6DoF (Actual)
            'Orange', [1.0000, 0.4980, 0.0549], ... % Ideal (Plan)
            'Red',    [0.8392, 0.1529, 0.1569], ... % Target / Control R
            'Green',  [0.1725, 0.6275, 0.1725], ... % Control L
            'Gray',   [0.4980, 0.4980, 0.4980], ... % Control a
            'Phase',  [0.5000, 0.5000, 0.5000]  ... % Phase Lines (Dark Gray)
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
            % ウィンドウ作成
            fig = figure('Name', 'Simulation Analysis', 'Color', 'w', ...
                'Position', [50, 50, 1200, 1000]); % 正方形に近い比率に変更
            
            % タブグループ
            tabGroup = uitabgroup(fig);
            
            % --- Tab 1: 3D View (Single Plot) ---
            t1 = uitab(tabGroup, 'Title', '3D View', 'BackgroundColor', 'w');
            obj.plotTab1_3D(t1);
            
            % --- Tab 2: Analysis (2x2: Top/Side/Horiz/Vert) ---
            t2 = uitab(tabGroup, 'Title', 'Analysis (2x2)', 'BackgroundColor', 'w');
            obj.plotTab2_Analysis(t2);
            
            % --- Tab 3: Dynamics (Controls etc.) ---
            t3 = uitab(tabGroup, 'Title', 'Dynamics', 'BackgroundColor', 'w');
            obj.plotTab3_Dynamics(t3);
        end
        
        %% --- Tab 1: 3D View (Single) ---
        function plotTab1_3D(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            
            % 1x1 全画面レイアウト
            t = tiledlayout(bgPanel, 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            simD = obj.SimData;
            planG = obj.TrajPlanGnd;
            tgt = obj.TargetPos;
            C = obj.Colors;
            
            ax = nexttile(t);
            plot3(ax, planG.Position(:,1), planG.Position(:,2), planG.Position(:,3), ...
                  '--', 'Color', C.Orange, 'LineWidth', 1.5, 'DisplayName', 'Ideal Trajectory');
            hold(ax, 'on');
            plot3(ax, simD.Inertial_X_Position, simD.Inertial_Y_Position, -simD.Inertial_Z_Position, ...
                  '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', '6DoF Simulation');
            plot3(ax, tgt(1), tgt(2), tgt(3), ...
                  '^', 'Color', C.Red, 'MarkerSize', 8, 'MarkerFaceColor', C.Red, 'DisplayName', 'Target');
            
            view(ax, 3); axis(ax, 'equal');
            xlabel(ax, 'North: {\itx} [m]'); ylabel(ax, 'East: {\ity} [m]'); zlabel(ax, 'Altitude: {\itz} [m]');
            legend(ax, 'show'); 
            obj.applyStyle(ax);
        end
        
        %% --- Tab 2: Analysis (2x2 Grid) ---
        function plotTab2_Analysis(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            
            % 2行2列レイアウト
            t = tiledlayout(bgPanel, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            simD = obj.SimData;
            planG = obj.TrajPlanGnd;
            tgt = obj.TargetPos;
            C = obj.Colors;
            
            % 誤差計算
            ref_pts = planG.Position; 
            act_pts = [simD.Inertial_X_Position, simD.Inertial_Y_Position, -simD.Inertial_Z_Position];
            n_act = size(act_pts, 1);
            err_vec = zeros(n_act, 3); 
            for i = 1:n_act
                [~, idx] = min(sum((ref_pts - act_pts(i, :)).^2, 2));
                err_vec(i, :) = act_pts(i, :) - ref_pts(idx, :); 
            end
            horiz_err = sqrt(err_vec(:,1).^2 + err_vec(:,2).^2);
            vert_err = abs(err_vec(:,3));

            % 1. Top View (左上)
            ax1 = nexttile(t);
            plot(ax1, planG.Position(:,1), planG.Position(:,2), '--', 'Color', C.Orange, 'LineWidth', 1.5, 'DisplayName', 'Ideal Trajectory');
            hold(ax1, 'on');
            plot(ax1, simD.Inertial_X_Position, simD.Inertial_Y_Position, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', '6DoF Simulation');
            plot(ax1, tgt(1), tgt(2), '^', 'Color', C.Red, 'MarkerSize', 8, 'MarkerFaceColor', C.Red, 'DisplayName', 'Target');
            axis(ax1, 'equal');
            xlabel(ax1, 'North: {\itx} [m]'); ylabel(ax1, 'East: {\ity} [m]');
            legend(ax1, 'show'); obj.applyStyle(ax1);
            
            % 2. Side View (右上)
            ax2 = nexttile(t);
            plot(ax2, planG.Position(:,1), planG.Position(:,3), '--', 'Color', C.Orange, 'LineWidth', 1.5, 'DisplayName', 'Ideal Trajectory');
            hold(ax2, 'on');
            plot(ax2, simD.Inertial_X_Position, -simD.Inertial_Z_Position, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', '6DoF Simulation');
            plot(ax2, tgt(1), tgt(3), '^', 'Color', C.Red, 'MarkerSize', 8, 'MarkerFaceColor', C.Red, 'DisplayName', 'Target');
            xlabel(ax2, 'North: {\itx} [m]'); ylabel(ax2, 'Altitude: {\itz} [m]');
            legend(ax2, 'show'); obj.applyStyle(ax2);
            
            % 3. Horizontal Error (左下)
            ax3 = nexttile(t);
            plot(ax3, simD.Time, horiz_err, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', 'Horiz. Error');
            xlabel(ax3, 'Time: {\itt} [s]'); ylabel(ax3, 'Horiz. Error [m]');
            legend(ax3, 'show'); obj.addPhaseLines(ax3); obj.applyStyle(ax3);
            
            % 4. Vertical Error (右下)
            ax4 = nexttile(t);
            plot(ax4, simD.Time, vert_err, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', 'Vert. Error');
            xlabel(ax4, 'Time: {\itt} [s]'); ylabel(ax4, 'Vert. Error [m]');
            legend(ax4, 'show'); obj.addPhaseLines(ax4); obj.applyStyle(ax4);
        end
        
        %% --- Tab 3: Dynamics (Bank, Alt, Control) ---
        function plotTab3_Dynamics(obj, parentTab)
            bgPanel = uipanel('Parent', parentTab, 'BackgroundColor', 'w', 'BorderType', 'none');
            t = tiledlayout(bgPanel, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            simD = obj.SimData;
            planA = obj.TrajPlanAir;
            planG = obj.TrajPlanGnd;
            C = obj.Colors;
            t_plan = planA.Time;
            t_max = max(t_plan(end), simD.Time(end));
            
            % Bank Angle
            ax1 = nexttile(t);
            plot(ax1, t_plan, -rad2deg(planA.Euler_RPY(:, 1)), '--', 'Color', C.Orange, 'LineWidth', 1.5, 'DisplayName', 'Ideal Trajectory');
            hold(ax1, 'on');
            plot(ax1, simD.Time, simD.Roll_Angle, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', '6DoF Simulation');
            xlabel(ax1, 'Time: {\itt} [s]'); ylabel(ax1, 'Bank Angle: {\it\phi} [deg]');
            xlim(ax1, [0, t_max]);
            legend(ax1, 'show'); obj.addPhaseLines(ax1); obj.applyStyle(ax1);
            
            % Altitude (Time Series)
            ax2 = nexttile(t);
            plot(ax2, t_plan, planG.Position(:, 3), '--', 'Color', C.Orange, 'LineWidth', 1.5, 'DisplayName', 'Ideal Trajectory');
            hold(ax2, 'on');
            plot(ax2, simD.Time, -simD.Inertial_Z_Position, '-', 'Color', C.Blue, 'LineWidth', 1.5, 'DisplayName', '6DoF Simulation');
            xlabel(ax2, 'Time: {\itt} [s]'); ylabel(ax2, 'Altitude: {\itz} [m]');
            xlim(ax2, [0, t_max]);
            legend(ax2, 'show'); obj.addPhaseLines(ax2); obj.applyStyle(ax2);
            
            % Control Input
            ax3 = nexttile(t);
            plot(ax3, simD.Time, simD.delta_R, '-', 'Color', C.Red,   'LineWidth', 1.0, 'DisplayName', 'Right: \delta_{\itR}');
            hold(ax3, 'on');
            plot(ax3, simD.Time, simD.delta_L, '-', 'Color', C.Green, 'LineWidth', 1.0, 'DisplayName', 'Left: \delta_{\itL}');
            plot(ax3, simD.Time, simD.delta_a, '--', 'Color', C.Gray, 'LineWidth', 1.0, 'DisplayName', 'Net: \delta_{\ita}');
            xlabel(ax3, 'Time: {\itt} [s]'); ylabel(ax3, 'Control Input (\delta)');
            xlim(ax3, [0, t_max]);
            legend(ax3, 'show'); obj.addPhaseLines(ax3); obj.applyStyle(ax3);
        end
    end
    
    methods (Access = private)
        % Phase Lines
        function addPhaseLines(obj, ax)
            times = obj.PlanTransitionTimes;
            if isempty(times), return; end
            y_lim = get(ax, 'YLim');
            hold(ax, 'on');
            for k = 1:length(times)
                t_switch = times(k);
                curr_xlim = get(ax, 'XLim');
                if t_switch < curr_xlim(1) || t_switch > curr_xlim(2), continue; end
                line(ax, [t_switch, t_switch], y_lim, ...
                    'Color', obj.Colors.Phase, 'LineStyle', ':', 'LineWidth', 1.2, 'HandleVisibility', 'off');
            end
        end

        % Style Application
        function applyStyle(~, ax)
            set(ax, 'FontName', 'Arial', 'FontSize', 14); % 2x2で見やすいようフォントサイズ調整(18->14)
            set(ax, 'TickLabelInterpreter', 'tex');
            set(ax.XLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k');
            set(ax.YLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k');
            if isprop(ax, 'ZLabel'), set(ax.ZLabel, 'Interpreter', 'tex', 'FontName', 'Arial', 'Color', 'k'); end
            title(ax, '');
            grid(ax, 'off');
            set(ax, 'Color', 'w'); set(ax, 'XColor', 'k'); set(ax, 'YColor', 'k'); set(ax, 'ZColor', 'k');
            set(ax, 'LineWidth', 1.0); set(ax, 'TickDir', 'in'); set(ax, 'TickLength', [0.015 0.025]); set(ax, 'Box', 'on');
            lgd = legend(ax); 
            if ~isempty(lgd)
                set(lgd, 'Interpreter', 'tex', 'FontName', 'Arial', 'FontSize', 12, ...
                    'Location', 'best', 'Box', 'off', 'TextColor', 'k');
            end
        end
    end
end