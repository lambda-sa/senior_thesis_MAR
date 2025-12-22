classdef ParafoilControlMapper < handle
    % PARAFOILCONTROLMAPPER
    % 軌道計画の結果(V, theta, phi)からフィードフォワード操作量を計算するクラス
    
    properties
        Yaw_Factor   % 事前計算する係数: -d * b * Cn_r / (2 * Cn_da)
    end
    
    methods
        function obj = ParafoilControlMapper(params)
            % コンストラクタ: 係数の事前計算
            
            % パラメータ構造体の階層に対応
            if isfield(params, 'prop')
                b = params.prop.b;
            else
                b = params.b;
            end
            
            if isfield(params, 'params')
                d = params.params.d;
                Cn_r = params.params.C_n_r;
                Cn_da = params.params.C_n_delta_a;
            else
                d = params.d;
                Cn_r = params.C_n_r;
                Cn_da = params.C_n_delta_a;
            end
            
            if abs(Cn_da) < 1e-9
                error('C_n_delta_a is too small. Cannot invert control logic.');
            end
            
            % 3DOFクラスと同じ係数計算
            obj.Yaw_Factor = -d * b * Cn_r / (2 * Cn_da);
        end
        
        % --- 既存のメソッド (状態ベクトルを使用する場合) ---
        function [delta_R, delta_L, delta_a] = compute_feedforward_input(obj, state_vec, phi_cmd, delta_s_bias)
            if nargin < 4, delta_s_bias = 0; end
            
            u = state_vec(1); v = state_vec(2); w = state_vec(3);
            theta = state_vec(8);
            V_sq = u^2 + v^2 + w^2;
            
            % 内部で下のメソッドを呼ぶ形にリファクタリングも可能ですが、
            % ここでは既存コードを壊さないよう独立させておきます
            g = 9.81;
            if V_sq < 1.0, delta_a = 0;
            else
                K_dyn = (g / V_sq) * cos(theta);
                delta_a = obj.Yaw_Factor * K_dyn * sin(phi_cmd);
            end
            
            [delta_R, delta_L] = obj.apply_mixing(delta_a, delta_s_bias);
        end
        
        % ★★★ 今回追加が必要なメソッド (参照値 V_ref, theta_ref を使用) ★★★
        function [delta_R, delta_L, delta_a] = compute_input_from_reference(obj, phi_ref, V_ref, theta_ref, delta_s_bias)
            % compute_input_from_reference
            % 計画された V_ref, theta_ref を用いて操作量を計算します。
            
            if nargin < 5, delta_s_bias = 0; end
            
            g = 9.81; 
            
            % ゼロ除算防止
            if V_ref < 1.0, V_ref = 1.0; end
            
            % 逆算ロジック (参照値を使用)
            % 式: delta_a = Factor * (g / V^2 * cos(theta)) * sin(phi)
            K_dyn = (g / (V_ref^2)) * cos(theta_ref);
            delta_a = obj.Yaw_Factor * K_dyn * sin(phi_ref);
            
            % ミキシング処理
            [delta_R, delta_L] = obj.apply_mixing(delta_a, delta_s_bias);
        end
    end
    
    methods (Access = private)
        % 共通ミキシング処理
        function [delta_R, delta_L] = apply_mixing(~, delta_a, delta_s_bias)
            % クランプ (物理的な限界 -1 ~ 1)
            delta_a = max(-1.0, min(1.0, delta_a));
            
            if delta_a > 0
                % 右旋回 (右ブレーキを引く)
                delta_R = delta_a + delta_s_bias;
                delta_L = delta_s_bias;
            else
                % 左旋回 (左ブレーキを引く)
                delta_R = delta_s_bias;
                delta_L = abs(delta_a) + delta_s_bias;
            end
            
            % 最終リミッター (0~1)
            delta_R = max(0, min(1, delta_R));
            delta_L = max(0, min(1, delta_L));
        end
    end
end