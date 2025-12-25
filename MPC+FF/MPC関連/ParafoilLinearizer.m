classdef ParafoilLinearizer < handle
    properties
        Plant
    end
    
    methods
        function obj = ParafoilLinearizer(plant_instance)
            obj.Plant = plant_instance;
        end
        
        function [A, B] = get_linear_model(obj, trim_state, trim_input_struct, h_step)
            % trim_input_struct: {delta_R, delta_L, GAMMA, wind_I}
            if nargin < 4, h_step = 1e-5; end
            
            nx = 12;
            nu = 2; % [delta_R, delta_L]
            
            A = zeros(nx, nx);
            B = zeros(nx, nu);
            
            % ベースライン
            y0 = obj.Plant.get_derivatives(0, trim_state, trim_input_struct);
            
            % --- A行列 (df/dx) ---
            for i = 1:nx
                x_p = trim_state;
                x_p(i) = x_p(i) + h_step;
                y_p = obj.Plant.get_derivatives(0, x_p, trim_input_struct);
                A(:, i) = (y_p - y0) / h_step;
            end
            
            % --- B行列 (df/du) ---
            % 1. delta_R を摂動
            u_R_p = trim_input_struct;
            u_R_p.delta_R = u_R_p.delta_R + h_step;
            y_R = obj.Plant.get_derivatives(0, trim_state, u_R_p);
            B(:, 1) = (y_R - y0) / h_step;
            
            % 2. delta_L を摂動
            u_L_p = trim_input_struct;
            u_L_p.delta_L = u_L_p.delta_L + h_step;
            y_L = obj.Plant.get_derivatives(0, trim_state, u_L_p);
            B(:, 2) = (y_L - y0) / h_step;
        end
    end
end