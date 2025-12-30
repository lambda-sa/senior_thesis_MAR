classdef PayloadPlotter
    methods (Static)
        function plotResults(simData, scheduler)
            if nargin < 2, scheduler = []; end
            
            disp('=== PayloadPlotter Debug Info ===');
            t = simData.Time;
            
            % --- データの準備 ---
            % NED座標: X=North, Y=East, Z=Down
            % プロット用: X=East, Y=North, Z=Alt(-Z)
            Data_X = simData.Inertial_Y_Position; % East
            Data_Y = simData.Inertial_X_Position; % North
            Data_Z = -simData.Inertial_Z_Position; % Altitude
            
            % ★デバッグ: データ範囲の確認 (コンソールに出力されます)
            fprintf('  East  (X) Range: %8.1f ~ %8.1f [m] (Diff: %.1f)\n', min(Data_X), max(Data_X), range(Data_X));
            fprintf('  North (Y) Range: %8.1f ~ %8.1f [m] (Diff: %.1f)\n', min(Data_Y), max(Data_Y), range(Data_Y));
            fprintf('  Alt   (Z) Range: %8.1f ~ %8.1f [m] (Diff: %.1f)\n', min(Data_Z), max(Data_Z), range(Data_Z));
            
            % 異常値警告
            if max(Data_Z) > 20000
                warning('高度データに 20,000m を超える値が含まれています。初期条件や発散を確認してください。');
            end
            
            % 1. メインウィンドウ作成
            fig = figure('Name', '6-DOF Simulation Results', ...
                         'NumberTitle', 'off', 'Color', 'w', ...
                         'Position', [50, 50, 1400, 900]);
            
            % 2. タブグループ作成
            tabGroup = uitabgroup(fig);
            
            % =========================================================
            % Tab 1: 3D Trajectory
            % =========================================================
            tab3D = uitab(tabGroup, 'Title', '3D Trajectory');
            ax3D = axes('Parent', tab3D);
            hold(ax3D, 'on');
            
            % 1. 全体パス (黒線)
            plot3(ax3D, Data_X, Data_Y, Data_Z, 'k-', 'LineWidth', 0.5, 'HandleVisibility', 'off');
            
            % 2. Phaseごとの色分け
            if ismember('Phase', simData.Properties.VariableNames)
                phases = unique(simData.Phase);
                colors = lines(length(phases));
                legendEntries = {};
                for k = 1:length(phases)
                    idx = (simData.Phase == phases(k));
                    if sum(idx) < 2, continue; end
                    
                    % 異常値(極端に遠い点)を除外してプロットする場合の簡易ロジック
                    % (ここではそのままプロットしますが、表示範囲制限で対応します)
                    plot3(ax3D, Data_X(idx), Data_Y(idx), Data_Z(idx), ...
                          'Color', colors(k,:), 'LineWidth', 2.0);
                    legendEntries{end+1} = ['Phase ' num2str(phases(k))];
                end
                if ~isempty(legendEntries)
                    legend(ax3D, legendEntries, 'Location', 'best');
                end
            else
                plot3(ax3D, Data_X, Data_Y, Data_Z, 'b-', 'LineWidth', 2.0);
            end
            
            % 始点・終点
            plot3(ax3D, Data_X(1), Data_Y(1), Data_Z(1), 'bs', 'MarkerFaceColor','b', 'MarkerSize', 10, 'DisplayName','Start');
            plot3(ax3D, Data_X(end), Data_Y(end), Data_Z(end), 'ro', 'MarkerFaceColor','r', 'MarkerSize', 10, 'DisplayName','End');
            
            grid(ax3D, 'on'); 
            box(ax3D, 'on');
            
            % --- ★最重要修正: 見た目の強制補正 ---
            
            % 1. 軸の数値を「物理等倍(equal)」ではなく「枠に合わせる(normal)」にする
            axis(ax3D, 'normal'); 
            
            % 2. プロットボックスのアスペクト比を強制的に [1 1 0.7] (横:奥:高) にする
            %    これにより、データがどれだけ偏っていても、見た目はサイコロ状の箱に収まる
            pbaspect(ax3D, [1 1 0.7]); 
            
            % 3. 表示範囲の制限 (異常な飛び値がある場合に備えてズーム)
            %    Z軸は 0m ～ 初期高度の1.2倍まで に制限
            z_max_limit = max(max(Data_Z(1)*1.2, 100), max(Data_Z)); 
            if z_max_limit > 15000 && Data_Z(1) < 5000 
                % もし初期高度が低いのにデータだけ異常に高い場合は、初期高度基準で切る
                z_max_limit = Data_Z(1) * 1.5; 
            end
            ylim(ax3D, [min(Data_Y)-50, max(Data_Y)+50]);
            xlim(ax3D, [min(Data_X)-50, max(Data_X)+50]);
            zlim(ax3D, [min(0, min(Data_Z)), z_max_limit]);
            
            view(ax3D, [-30, 30]); 
            rotate3d(ax3D, 'on'); 
            
            xlabel(ax3D, 'East [m]'); 
            ylabel(ax3D, 'North [m]'); 
            zlabel(ax3D, 'Altitude [m]');
            title(ax3D, '3D Flight Trajectory (Scale Adjusted)');
            hold(ax3D, 'off');
            
            
            % =========================================================
            % Tab 2: Flight Data Analysis
            % =========================================================
            tabAnalysis = uitab(tabGroup, 'Title', 'Flight Data Analysis');
            
            % --- 1. 2D Ground Track ---
            ax1 = subplot(3,3,1, 'Parent', tabAnalysis);
            plot(ax1, Data_Y, Data_X, 'b-', 'LineWidth', 1.5); hold(ax1, 'on');
            plot(ax1, Data_Y(end), Data_X(end), 'rx', 'MarkerSize', 8, 'LineWidth', 2);
            xlabel(ax1, 'East [m]'); ylabel(ax1, 'North [m]');
            %title(ax1, '2D Ground Track');
            axis(ax1, 'equal'); %grid(ax1, 'on');

            % --- 2. Altitude ---
            ax2 = subplot(3,3,2, 'Parent', tabAnalysis);
            plot(ax2, t, Data_Z, 'k-', 'LineWidth', 1.5);
            xlabel(ax2, 'Time [s]'); ylabel(ax2, 'Altitude [m]');
            %title(ax2, 'Altitude'); grid(ax2, 'on');

            % --- 3. Body Velocities ---
            ax3 = subplot(3,3,3, 'Parent', tabAnalysis);
            plot(ax3, t, simData.Body_U_Vel, 'b-', 'DisplayName', 'u'); hold(ax3, 'on');
            plot(ax3, t, simData.Body_V_Vel, 'g-', 'DisplayName', 'v');
            plot(ax3, t, simData.Body_W_Vel, 'r-', 'DisplayName', 'w');
            xlabel(ax3, 'Time [s]'); ylabel(ax3, 'velocity [m/s]');
            %title(ax3, 'Body Velocities'); 
            legend(ax3,'Location','best'); %grid(ax3, 'on');

            % --- 4. Yaw Angle & Rate ---
            ax4 = subplot(3,3,4, 'Parent', tabAnalysis);
            yyaxis(ax4, 'left');
            plot(ax4, t, simData.Yaw_Angle, 'b-', 'LineWidth', 1.5); ylabel(ax4, 'Yaw [deg]');
            yyaxis(ax4, 'right');
            plot(ax4, t, rad2deg(simData.Eul_Psi_Dot), 'r:', 'LineWidth', 1.5); ylabel(ax4, 'Rate [deg/s]');
            xlabel(ax4, 'Time [s]'); 
            %title(ax4, 'Yaw'); grid(ax4, 'on');

            % --- 5. Attitude ---
            ax5 = subplot(3,3,5, 'Parent', tabAnalysis);
            plot(ax5, t, simData.Roll_Angle, 'b-', 'DisplayName','Roll'); hold(ax5, 'on');
            plot(ax5, t, simData.Pitch_Angle, 'r-', 'DisplayName','Pitch');
            xlabel(ax5, 'Time [s]'); ylabel(ax5, 'Angle [deg]');
            legend(ax5,'Location','best'); 
            %title(ax5, 'Attitude'); grid(ax5, 'on');

            % --- 6. Body Rates ---
            ax6 = subplot(3,3,6, 'Parent', tabAnalysis);
            plot(ax6, t, rad2deg(simData.Body_P_Vel), 'b-', 'DisplayName','p'); hold(ax6, 'on');
            plot(ax6, t, rad2deg(simData.Body_Q_Vel), 'r-', 'DisplayName','q');
            plot(ax6, t, rad2deg(simData.Body_R_Vel), 'g-', 'DisplayName','r');
            xlabel(ax6, 'Time [s]'); ylabel(ax6, 'Angle [deg]');
            legend(ax6,'Location','best'); 
            %title(ax6, 'Body Rates'); grid(ax6, 'on');

            % --- 7. Control Inputs ---
            ax7 = subplot(3,3,7, 'Parent', tabAnalysis);
            if ismember('delta_R', simData.Properties.VariableNames)
                yyaxis(ax7, 'left');
                plot(ax7, t, simData.delta_R, 'r-'); hold(ax7, 'on');
                plot(ax7, t, simData.delta_L, 'g-'); ylabel(ax7, 'Brake');
                yyaxis(ax7, 'right');
                plot(ax7, t, simData.delta_a, 'b-'); ylabel(ax7, 'Asym');
                legend(ax7,'Location','best');
            else
                text(ax7, 0.5, 0.5, 'No Data');
            end
            title(ax7, 'Control'); grid(ax7, 'on');

            % --- 8. Inertial Vel ---
            ax8 = subplot(3,3,8, 'Parent', tabAnalysis);
            plot(ax8, t, simData.Inertial_Vx, 'b-', 'DisplayName', 'V_x'); hold(ax8, 'on');
            plot(ax8, t, simData.Inertial_Vy, 'g-', 'DisplayName', 'V_y');
            plot(ax8, t, simData.Inertial_Vz, 'r-', 'DisplayName', 'V_z');
            xlabel(ax8, 'Time [s]'); ylabel(ax8, 'm/s');
            %title(ax8, 'Inertial Vel'); 
            legend(ax8,'Location','best'); %grid(ax8, 'on');

            % --- 9. Phase ---
            ax9 = subplot(3,3,9, 'Parent', tabAnalysis);
            if ismember('Phase', simData.Properties.VariableNames)
                plot(ax9, t, simData.Phase, 'k-');
                %title(ax9, 'Phase');
            end
            %grid(ax9, 'on');

            % ★追加項目 Tab 3: Attitude History (3軸独立表示)
            % =========================================================
            tabAttitude = uitab(tabGroup, 'Title', 'Attitude History');
            
            % --- Roll ---
            axA1 = subplot(3,1,1, 'Parent', tabAttitude);
            plot(axA1, t, simData.Roll_Angle, 'b-', 'LineWidth', 1.5);
            ylabel(axA1, 'Roll (phi) [deg]');
            grid(axA1, 'on'); title(axA1, 'Attitude Angles vs Time');
            
            % --- Pitch ---
            axA2 = subplot(3,1,2, 'Parent', tabAttitude);
            plot(axA2, t, simData.Pitch_Angle, 'r-', 'LineWidth', 1.5);
            ylabel(axA2, 'Pitch (theta) [deg]');
            grid(axA2, 'on');
            
            % --- Yaw ---
            axA3 = subplot(3,1,3, 'Parent', tabAttitude);
            plot(axA3, t, simData.Yaw_Angle, 'g-', 'LineWidth', 1.5);
            ylabel(axA3, 'Yaw (psi) [deg]');
            xlabel(axA3, 'Time [s]');
            grid(axA3, 'on');

            % 3つのグラフの横軸を同期（ズームした時に連動するようになります）
            linkaxes([axA1, axA2, axA3], 'x');

        end
    end
end