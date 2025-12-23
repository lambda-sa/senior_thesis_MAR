function [px, py, pyaw, mode, lengths, bcost] = dubins_path_planning(s, e, c, v, step, start_constraint)
% DUBINS_PATH_PLANNING (With Start Constraint)
% start_constraint: 'L', 'R', or []

    if nargin < 6, start_constraint = []; end

    sx = s(1); sy = s(2); syaw = s(3);
    ex = e(1); ey = e(2); eyaw = e(3);

    ex_local = ex - sx; ey_local = ey - sy;
    lex = cos(syaw) * ex_local + sin(syaw) * ey_local;
    ley = -sin(syaw) * ex_local + cos(syaw) * ey_local;
    leyaw = eyaw - syaw;

    [lpx, lpy, lpyaw, mode, lengths, bcost] = dubins_path_planning_from_origin(lex, ley, leyaw, c, step, start_constraint);

    if ~isempty(lpx)
        px = cos(-syaw) .* lpx + sin(-syaw) .* lpy + sx;
        py = -sin(-syaw) .* lpx + cos(-syaw) .* lpy + sy;
        pyaw = mod2pi(lpyaw + syaw);
    else
        px=[]; py=[]; pyaw=[];
    end
end

function [px, py, pyaw, mode, lengths, bcost] = dubins_path_planning_from_origin(ex, ey, eyaw, c, step, start_cons)
    dx = ex; dy = ey; D = sqrt(dx^2 + dy^2); d = D * c;
    theta = mod2pi(atan2(dy, dx)); alpha = mod2pi(-theta); beta  = mod2pi(eyaw - theta);

    % Mode definitions
    all_planners = {'LSL', @LSL; 'RSR', @RSR; 'LSR', @LSR; 'RSL', @RSL; 'RLR', @RLR; 'LRL', @LRL};

    bcost = inf; bt=[]; bp=[]; bq=[]; bmode=[];

    for i = 1:size(all_planners, 1)
        mode_name = all_planners{i, 1};
        planner_func = all_planners{i, 2};
        
        % ★ 制約チェック: 指定された方向で始まらないモードはスキップ
        if ~isempty(start_cons)
            if mode_name(1) ~= start_cons
                continue; 
            end
        end

        [t, p, q, mode] = planner_func(alpha, beta, d);
        if isempty(t), continue; end
        
        cost = (abs(t) + abs(p) + abs(q));
        if bcost > cost
            bt = t; bp = p; bq = q; bmode = mode; bcost = cost;
        end
    end
    
    lengths = [bt, bp, bq] ./ c;
    [px, py, pyaw] = generate_course(lengths, bmode, c, step);
    mode = bmode;
end

% --- Helper Functions (Same as before) ---
function [px, py, pyaw] = generate_course(lengths, mode, c, step_size)
    % GENERATE_COURSE (厳密解版)
    % オイラー積分(近似)ではなく、円の幾何学式を用いて正確な座標を計算します。
    
    px = [0.0]; py = [0.0]; pyaw = [0.0];
    
    for i = 1:length(mode)
        m = mode{i}; 
        l = lengths(i);
        
        if l < 1e-6, continue; end
        
        % ステップ数の計算
        n_steps = floor(l / step_size);
        
        % メインループ
        for k = 1:n_steps
            cx = px(end); cy = py(end); cyaw = pyaw(end);
            
            if strcmp(m, 'S')
                % 直線: そのまま進む
                nx = cx + step_size * cos(cyaw);
                ny = cy + step_size * sin(cyaw);
                nyaw = cyaw;
            else
                % 旋回: 厳密な円弧計算
                % 弦長近似ではなく、積分形を使用:
                % x = integral(cos(theta)) = sin(theta)
                % y = integral(sin(theta)) = -cos(theta)
                
                d_ang = step_size * c; % 角度変化量
                
                if strcmp(m, 'L')
                    % 左旋回 (角度が増える)
                    % nx = cx + (sin(cyaw + d_ang) - sin(cyaw)) / c;
                    % ny = cy + (cos(cyaw) - cos(cyaw + d_ang)) / c;
                    % nyaw = cyaw + d_ang;
                    
                    % 安定性のために相対移動量として計算
                    nx = cx + (sin(cyaw + d_ang) - sin(cyaw)) / c;
                    ny = cy - (cos(cyaw + d_ang) - cos(cyaw)) / c;
                    nyaw = cyaw + d_ang;
                    
                elseif strcmp(m, 'R')
                    % 右旋回 (角度が減る)
                    % 右旋回は曲率が -c 扱いと同義
                    % nx = cx + (sin(cyaw - d_ang) - sin(cyaw)) / -c;
                    % ny = cy + (cos(cyaw) - cos(cyaw - d_ang)) / -c;
                    
                    nx = cx - (sin(cyaw - d_ang) - sin(cyaw)) / c;
                    ny = cy + (cos(cyaw - d_ang) - cos(cyaw)) / c;
                    nyaw = cyaw - d_ang;
                end
            end
            px(end+1) = nx; py(end+1) = ny; pyaw(end+1) = nyaw;
        end
        
        % 端数処理 (Rem) - 最後の隙間を埋める
        rem = l - n_steps * step_size;
        if rem > 1e-4
            cx = px(end); cy = py(end); cyaw = pyaw(end);
            if strcmp(m, 'S')
                px(end+1) = cx + rem * cos(cyaw);
                py(end+1) = cy + rem * sin(cyaw);
                pyaw(end+1) = cyaw;
            else
                d_rem = rem * c;
                if strcmp(m, 'L')
                    px(end+1) = cx + (sin(cyaw + d_rem) - sin(cyaw)) / c;
                    py(end+1) = cy - (cos(cyaw + d_rem) - cos(cyaw)) / c;
                    pyaw(end+1) = cyaw + d_rem;
                elseif strcmp(m, 'R')
                    px(end+1) = cx - (sin(cyaw - d_rem) - sin(cyaw)) / c;
                    py(end+1) = cy + (cos(cyaw - d_rem) - cos(cyaw)) / c;
                    pyaw(end+1) = cyaw - d_rem;
                end
            end
        end
    end
end
% (LSL, RSR等の関数定義は省略なしで実装してください。以前のコードと同じです)
function [t, p, q, mode] = LSL(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    tmp0 = d + sa - sb; p_squared = 2 + d*d - 2*c_ab + 2*d*(sa - sb);
    if p_squared < 0, t=[]; p=[]; q=[]; mode={'L','S','L'}; return; end
    tmp1 = atan2((cb - ca), tmp0); t = mod2pi(-alpha + tmp1); p = sqrt(p_squared); q = mod2pi(beta - tmp1); mode = {'L', 'S', 'L'};
end
function [t, p, q, mode] = RSR(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    tmp0 = d - sa + sb; p_squared = 2 + d*d - 2*c_ab + 2*d*(sb - sa);
    if p_squared < 0, t=[]; p=[]; q=[]; mode={'R','S','R'}; return; end
    tmp1 = atan2((ca - cb), tmp0); t = mod2pi(alpha - tmp1); p = sqrt(p_squared); q = mod2pi(-beta + tmp1); mode = {'R', 'S', 'R'};
end
function [t, p, q, mode] = LSR(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    p_squared = -2 + d*d + 2*c_ab + 2*d*(sa + sb);
    if p_squared < 0, t=[]; p=[]; q=[]; mode={'L','S','R'}; return; end
    p = sqrt(p_squared); tmp2 = atan2((-ca - cb), (d + sa + sb)) - atan2(-2.0, p); t = mod2pi(-alpha + tmp2); q = mod2pi(-mod2pi(beta) + tmp2); mode = {'L', 'S', 'R'};
end
function [t, p, q, mode] = RSL(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    p_squared = d*d - 2 + 2*c_ab - 2*d*(sa + sb);
    if p_squared < 0, t=[]; p=[]; q=[]; mode={'R','S','L'}; return; end
    p = sqrt(p_squared); tmp2 = atan2((ca + cb), (d - sa - sb)) - atan2(2.0, p); t = mod2pi(alpha - tmp2); q = mod2pi(beta - tmp2); mode = {'R', 'S', 'L'};
end
function [t, p, q, mode] = RLR(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    tmp_rlr = (6.0 - d*d + 2.0*c_ab + 2.0*d*(sa - sb)) / 8.0;
    if abs(tmp_rlr) > 1.0, t=[]; p=[]; q=[]; mode={'R','L','R'}; return; end
    p = mod2pi(2*pi - acos(tmp_rlr)); t = mod2pi(alpha - atan2(ca - cb, d - sa + sb) + mod2pi(p/2.0)); q = mod2pi(alpha - beta - t + mod2pi(p)); mode = {'R', 'L', 'R'};
end
function [t, p, q, mode] = LRL(alpha, beta, d)
    sa = sin(alpha); sb = sin(beta); ca = cos(alpha); cb = cos(beta); c_ab = cos(alpha - beta);
    tmp_lrl = (6.0 - d*d + 2.0*c_ab + 2.0*d*(-sa + sb)) / 8.0;
    if abs(tmp_lrl) > 1.0, t=[]; p=[]; q=[]; mode={'L','R','L'}; return; end
    p = mod2pi(2*pi - acos(tmp_lrl)); t = mod2pi(-alpha - atan2(ca - cb, d + sa - sb) + p/2.0); q = mod2pi(mod2pi(beta) - alpha - t + mod2pi(p)); mode = {'L', 'R', 'L'};
end
function angle = mod2pi(theta)
    angle = theta - 2.0 * pi * floor(theta / (2.0 * pi));
end