%DEMO_SIMULATION  Simulate the mean of 9 frames with known motion and estimate it.
%   The simulated image is equally one exposure with 9 short flashes.
%   Run from anywhere in MATLAB or Octave:
%     run('matlab/examples/demo_simulation.m')

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

K = 9;                      % number of frames (or flashes)
L = (K - 1)/2;              % L = 4, so K = 2L+1 = 9 (simulate_pulsed_image takes L)
v_true = [5 -3];            % px per step, x right, y down
n = 256;                    % analysed window, n x n pixels

% Moving window over a real surface recorded with the Basler camera.
% Replace 'Texture', tex by 'CorrLength', 1.5 for a random texture.
tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
[bp, frames] = fftvel.simulate_pulsed_image(n, v_true, L, 'Texture', tex, 'Noise', 0.1);

est = fftvel.estimate_velocity(bp, 'Exposures', K, 'Diagnostics', true);

true_angle = mod(atan2d(-v_true(2), v_true(1)), 180);
fprintf('true:      direction %6.2f deg, speed %6.3f px per frame\n', true_angle, norm(v_true));
fprintf('estimate:  direction %6.2f deg, speed %6.3f px per frame (quality %.2f, GEBA %.2f)\n', ...
    est.angle_deg, est.speed, est.quality, est.geba);

fftvel.plot_estimate(bp, est, frames(:, :, L + 1));    % the middle frame
