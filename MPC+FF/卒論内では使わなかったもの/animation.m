% --- シミュレーション実行後に追加 ---

fprintf('\n=== Starting Dual-View Animation ===\n');

% アニメーションの再生速度倍率 (例: 10倍速)
% 数値を大きくすると速く、小さくすると遅くなります。
speed_factor = 40.0; 

% 対地(Actual) と 対気(Reference) の比較アニメーションを実行
% mission は ParafoilMissionWind のインスタンスである必要があります
mission.animate_dual_views(speed_factor);

fprintf('Animation finished.\n');