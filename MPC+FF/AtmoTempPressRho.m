% --- AtmoTempPressRho.m ---
classdef AtmoTempPressRho < AtmosphericTemperatureEstimator
    properties (GetAccess = protected)
        
        file_pre = "data_pre.csv";
        file_rho = "data_rho.csv";
    end
    properties (Constant)
        g0 = 9.80665;    % 標準重力加速度 [m/s^2]
        Re = 6371000;    % 地球の平均半径 [m]
        m_0 = 28.9644*10^(-3); % kg/mol
        r_gus = 8.31432; % Nm/(mol/K)
    end

    methods
        function obj = AtmoTempPressRho()
            obj@AtmosphericTemperatureEstimator('P4_input.csv');
        end

        function p = get_press(obj, altitude)
            %if  altitude > 150
                %error('高度は150 以下でなければなりません。');
            if altitude < 0 || altitude > 150
                error('高度は0 以上で、150 以下でなければなりません。');
                
            end
            if altitude <= 86
                T = obj.get_temperature(altitude);
                p = obj.approximate_below_86_press(altitude, T);
            else
                [x0, y0, x1, y1] = obj.find_bracketing_points(obj.file_pre, altitude);
                if isnan(x0)
                    error('補間対象の範囲がCSV 内に見つかりませんでした。');
                end
                p = y0 + (altitude - x0) * (y1 - y0) / (x1 - x0); % 手動線形補間
            end
        end

        function rho = get_density(obj, altitude)
            if altitude < 0 || altitude > 150
                error('高度は0 以上で、150 以下でなければなりません。');
            %if  altitude > 150
                %error('高度は150 以下でなければなりません。');
            end
            if altitude <= 86
                T = obj.get_temperature(altitude);
                p = obj.get_press(altitude);
                rho = p * obj.m_0 / (obj.r_gus * T);
            else
                [x0, y0, x1, y1] = obj.find_bracketing_points(obj.file_rho, altitude);
                if isnan(x0)
                    error('補間対象の範囲がCSV 内に見つかりませんでした。');
                end
                rho = y0 + (altitude - x0) * (y1 - y0) / (x1 - x0); % 手動線形補間
            end
        end
        function g = get_gravity(obj, altitude_m)
            % 高度h[m]における重力加速度を計算するメソッド
            g = obj.g0 * (obj.Re / (obj.Re + altitude_m))^2;
        end
    end

    methods (Access = protected)
        function p = approximate_below_86_press(~, z, T)
            r_0 = 6356.756;
            h = r_0 * z / (r_0 + z);
            if h >= 0 && h <= 11
                p = 101325 * (288.15/T)^(-5.256);
            elseif h > 11 && h <= 20
                p = 22632.264 * exp(-0.1577 * (h-11));
            elseif h > 20 && h <= 32
                p = 5474.889 * (216.65/T)^(34.163);
            elseif h > 32 && h <= 47
                p = 868.019 * (228.65/T)^(12.201);
            elseif h > 47 && h <= 51
                p = 110.906 * exp(-0.1262 * (h-47));
            elseif h > 51 && h <= 71
                p = 66.939 * (270.65/T)^(-12.201);
            elseif h > 71 && h <= 84.852
                p = 3.956 * (214.65/T)^(-17.082);
            else
                p = NaN;
            end
        end
    end
end