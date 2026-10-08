%DEMO_SIMULATION  Simulate a pulsed-exposure image with known motion and estimate it.
%   Run from anywhere in MATLAB or Octave:
%     run('matlab/examples/demo_simulation.m')

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

L = 4;                      % 2L+1 = 9 exposures
v_true = [5 -3];            % px per exposure step, x right, y down
n = 256;                    % analysed window, n x n pixels

% Moving window over a real surface recorded with the Basler camera.
% Replace 'Texture', tex by 'CorrLength', 1.5 for a random texture.
tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
[bp, frames] = simulate_pulsed_image(n, v_true, L, 'Texture', tex, 'Noise', 0.1);

est = estimate_velocity(bp, 'L', L, 'Diagnostics', true);

true_angle = mod(atan2d(-v_true(2), v_true(1)), 180);
fprintf('true:      direction %6.2f deg, speed %6.3f px/step\n', true_angle, norm(v_true));
fprintf('estimate:  direction %6.2f deg, speed %6.3f px/step (quality %.2f, GEBA %.2f)\n', ...
    est.angle_deg, est.speed, est.quality, est.geba);

plot_estimate(bp, est, frames(:, :, L + 1));
