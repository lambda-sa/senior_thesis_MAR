% --- AtmosphericTemperatureEstimator.m ---
classdef AtmosphericTemperatureEstimator < handle
    properties (Access = protected)
        file_path;
        altitude_data;
        temperature_data;
    end
    
    methods
        function obj = AtmosphericTemperatureEstimator(file_path)
            if nargin > 0
                obj.file_path = file_path;
                obj.load_temperature_data();
            end
        end
        
        function T = get_temperature(obj, altitude)
            if altitude < 0 || altitude > 150
                error('高度は0 以上で、150 以下でなければなりません。');
             %if  altitude > 150
                %error('高度は150 以下でなければなりません。');
            end
            
            % 86km以下の場合は近似式を使用
            if altitude <= 86
                T = obj.approximate_below_86_temp(altitude);
            else
                % CSVデータからの線形補間
                [x0, y0, x1, y1] = obj.find_bracketing_points(obj.file_path, altitude);
                if isnan(x0)
                    error('補間対象の範囲がCSV 内に見つかりませんでした。');
                end
                T = y0 + (altitude - x0) * (y1 - y0) / (x1 - x0); % 手動線形補間
            end
        end
    end
    
    methods (Access = protected)
        function load_temperature_data(obj)
            data = readtable(obj.file_path);
            obj.altitude_data = data{:, 1}; % 1列目が高度
            obj.temperature_data = data{:, 2}; % 2列目が温度
        end
        
        function [x0, y0, x1, y1] = find_bracketing_points(~, filePath, target_altitude)
            % CSVファイルを読み込む
            data = readtable(filePath);
            altitudes = data{:,1}; % 高度データ
            values = data{:,2};    % 対応する値 (温度, 圧力, 密度など)

            % 補間区間を検索
            idx = find(altitudes <= target_altitude, 1, 'last');

            if isempty(idx) || idx == length(altitudes)
                % ターゲット高度がデータの範囲外の場合、NaNを返す
                x0 = NaN; y0 = NaN; x1 = NaN; y1 = NaN;
                return;
            end

            x0 = altitudes(idx);
            y0 = values(idx);
            x1 = altitudes(idx+1);
            y1 = values(idx+1);
        end
        
        function T = approximate_below_86_temp(~, h_m)
            r_0 = 6356.756; % 地球の平均半径 (km)
            h = r_0 * h_m / (r_0 + h_m); % 地上からの高度 h (km)

            if h >= 0 && h <= 11
                T = 288.15 - 6.5 * h;
            elseif h > 11 && h <= 20
                T = 216.65;
            elseif h > 20 && h <= 32
                T = 216.65 + 1.0 * (h - 20);
            elseif h > 32 && h <= 47
                T = 228.65 + 2.8 * (h - 32);
            elseif h > 47 && h <= 51
                T = 270.65;
            elseif h > 51 && h <= 71
                T = 270.65 - 2.8 * (h - 51);
            elseif h > 71 && h <= 84.852
                T = 214.65 - 2.0 * (h - 71);
            else
                T = NaN; % 範囲外
            end
        end
    end
end