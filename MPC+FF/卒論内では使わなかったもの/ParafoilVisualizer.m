classdef ParafoilVisualizer < handle
    properties
        VideoFileName = 'parafoil_flight.mp4';
        FrameRate = 30;
        SpeedMultiplier = 5; % 再生速度 (1なら等倍, 5なら5倍速)
        WingSpan = 10;       % 描画用の翼幅 [m]
        WingChord = 3;       % 描画用の翼弦長 [m]
        AxisRange = 20;      % 追従カメラの表示範囲 [m]
    end
    
    methods
        function obj = ParafoilVisualizer(fileName)
            if nargin > 0
                obj.VideoFileName = fileName;
            end
        end
        
        function animate(obj, time, state)
            % animate 3D軌道を描画し動画を保存する
            % Input:
            %   time: 時間ベクトル [N x 1]
            %   state: 状態行列 [N x 12] (u,v,w, p,q,r, phi,theta,psi, N,E,D)
            
            % --- データの抽出と座標変換 ---
            % MATLABのplot3は通常 (X, Y, Z) = (East, North, Up) で描くと見やすいですが、
            % ここでは航空力学慣例 (North, East, Down) を (X, Y, -Z) にマッピングします。
            
            % NED座標系
            pos_N = state(:, 10);
            pos_E = state(:, 11);
            pos_D = state(:, 12);
            
            % 姿勢角 (Euler Angles)
            phi   = state(:, 7);
            theta = state(:, 8);
            psi   = state(:, 9);
            
            % 描画用座標 (X:North, Y:East, Z:Altitude)
            X_plot = pos_N;
            Y_plot = pos_E;
            Z_plot = -pos_D; % 高度はマイナスD
            
            % --- 動画保存の準備 ---
            v = VideoWriter(obj.VideoFileName, 'MPEG-4');
            v.FrameRate = obj.FrameRate;
            open(v);
            
            % --- 図の準備 ---
            fig = figure('Name', 'Parafoil 3D Animation', 'Color', 'w', 'Position', [100 100 800 600]);
            ax = axes('Parent', fig);
            grid(ax, 'on');
            axis(ax, 'equal');
            hold(ax, 'on');
            xlabel(ax, 'North [m]');
            ylabel(ax, 'East [m]');
            zlabel(ax, 'Altitude [m]');
            view(ax, 3); % 3Dビュー
            
            % 軌跡のライン（初期化）
            trail = plot3(ax, X_plot(1), Y_plot(1), Z_plot(1), 'b-', 'LineWidth', 1);
            
            % 機体モデル（初期化）
            % パッチ（多角形）で翼を描画
            wing_patch = patch(ax, 'XData', [], 'YData', [], 'ZData', [], 'FaceColor', 'r', 'FaceAlpha', 0.8);
            
            % 経過時間の表示
            time_text = text(ax, 0, 0, 0, '', 'Units', 'normalized', 'Position', [0.05 0.95], 'FontSize', 12);

            % --- アニメーションループ ---
            dt_sim = time(2) - time(1);
            skip_step = round(1 / (dt_sim * obj.FrameRate) * obj.SpeedMultiplier);
            if skip_step < 1, skip_step = 1; end
            
            total_frames = length(time);
            
            fprintf('Creating Animation: %s...\n', obj.VideoFileName);
            
            for i = 1:skip_step:total_frames
                if ~isvalid(fig), break; end % ウィンドウが閉じられたら終了
                
                % 現在の位置と姿勢
                current_pos = [X_plot(i); Y_plot(i); Z_plot(i)];
                current_att = [phi(i); theta(i); psi(i)];
                
                % 1. 軌跡の更新
                set(trail, 'XData', X_plot(1:i), 'YData', Y_plot(1:i), 'ZData', Z_plot(1:i));
                
                % 2. 機体(翼)の描画更新
                obj.update_wing_shape(wing_patch, current_pos, current_att);
                
                % 3. カメラと軸の更新 (機体中心に追従)
                range = obj.AxisRange;
                xlim(ax, [current_pos(1)-range, current_pos(1)+range]);
                ylim(ax, [current_pos(2)-range, current_pos(2)+range]);
                zlim(ax, [current_pos(3)-range, current_pos(3)+range]);
                
                % 視点を少し機体後方上空に固定する場合（お好みでコメントアウト解除）
                % cam_offset = [-30; -10; 10]; % 北・東・上へのオフセット
                % campos(ax, current_pos + cam_offset);
                % camtarget(ax, current_pos);
                
                % 4. テキスト更新
                set(time_text, 'String', sprintf('Time: %.2fs', time(i)));
                
                % 5. フレーム保存
                frame = getframe(fig);
                writeVideo(v, frame);
                
                drawnow limitrate;
            end
            
            close(v);
            fprintf('Animation saved successfully.\n');
        end
    end
    
    methods (Access = private)
        function update_wing_shape(obj, patch_handle, pos, att)
            % 翼の形状定義 (Body座標系: x=前方, y=右, z=下)
            % MATLABプロット用には z=上 に変換して扱う
            b = obj.WingSpan / 2;
            c = obj.WingChord;
            
            % 翼の4頂点 (Body系) [x; y; z]
            % 前縁左, 前縁右, 後縁右, 後縁左
            vertices_B = [
                 c/2, -b,  0;
                 c/2,  b,  0;
                -c/2,  b,  0;
                -c/2, -b,  0;
            ]'; 
        
            % 回転行列 (NED -> Body の逆 = Body -> NED)
            phi = att(1); theta = att(2); psi = att(3);
            
            % MATLABのplot3座標系(X=N, Y=E, Z=Up)に合わせるための変換
            % 航空力学(NED)での回転行列 R_nb (Body to Nav)
            % R_nb = R_z(psi) * R_y(theta) * R_x(phi)
            
            cph = cos(phi); sph = sin(phi);
            cth = cos(theta); sth = sin(theta);
            cps = cos(psi); sps = sin(psi);
            
            % Body(NED) -> Inertial(NED)
            R_nb = [
                cth*cps, sin(phi)*sth*cps - cos(phi)*sps, cos(phi)*sth*cps + sin(phi)*sps;
                cth*sps, sin(phi)*sth*sps + cos(phi)*cps, cos(phi)*sth*sps - sin(phi)*cps;
                -sth,    sin(phi)*cth,                    cos(phi)*cth
            ];
            
            % 頂点の回転
            vertices_NED = R_nb * vertices_B;
            
            % NED -> Plot座標 (Z_plot = -Z_ned) への変換
            % X_plot = X_ned
            % Y_plot = Y_ned
            % Z_plot = -Z_ned
            vertices_Plot = [
                vertices_NED(1, :);
                vertices_NED(2, :);
                -vertices_NED(3, :)
            ];
            
            % 並進移動
            vertices_Final = vertices_Plot + pos;
            
            % パッチの更新
            set(patch_handle, ...
                'XData', vertices_Final(1, :), ...
                'YData', vertices_Final(2, :), ...
                'ZData', vertices_Final(3, :));
        end
    end
end