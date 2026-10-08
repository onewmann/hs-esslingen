%DEMO_TESTSEQUENCES  Sliding 9-frame windows over the test sequences.
%   The sequences in data/testsequences.mat move by a known amount per frame.
%   Each frame is the previous one shifted with zero fill, so every window
%   is cropped to the area that holds image content in all of its frames.
%   The reference is the mean frame-to-frame shift from phase correlation.
%   Run in MATLAB or Octave:
%     run('matlab/examples/demo_testsequences.m')

here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));
d = load(fullfile(here, '..', '..', 'data', 'testsequences.mat'));
L = 4;
names = {'model', 'tyre', 'accelerating'};

figure('Name', 'Test sequences', 'Color', 'w');
for k = 1:numel(names)
    s = double(d.(names{k}))/255;
    nf = size(s, 3);
    shift = zeros(nf - 1, 2);
    for f = 1:nf - 1
        [dx, dy] = phase_correlation(s(:, :, f), s(:, :, f + 1));
        shift(f, :) = [dx dy];
    end
    starts = 1:(nf - 2*L);
    est_speed = nan(size(starts));
    ref_speed = nan(size(starts));
    valid = false(size(starts));
    for w = 1:numel(starts)
        idx = starts(w) + (0:2*L);
        border = ceil(max(abs(sum(shift(1:idx(end) - 1, :), 1))));   % zero fill
        bp = mean(s(border + 1:end, border + 1:end, idx), 3);
        e = estimate_velocity(bp, 'L', L);
        est_speed(w) = e.speed;
        valid(w) = e.valid;
        ref_speed(w) = norm(mean(shift(idx(1:end-1), :), 1));
    end
    subplot(1, numel(names), k);
    plot(starts, ref_speed, 'k-', 'LineWidth', 1.5); hold on;
    plot(starts(valid), est_speed(valid), 'o', 'Color', [0.16 0.47 0.84], 'MarkerFaceColor', [0.16 0.47 0.84]);
    plot(starts(~valid), est_speed(~valid), 'o', 'Color', [0.92 0.41 0.20]);
    hold off;
    ylim([0 1.5*max(ref_speed)]);
    xlabel('first frame of window'); ylabel('speed (px/frame)');
    title(names{k});
    if k == 1
        legend('phase correlation', 'estimate (valid)', 'estimate (not valid)', 'Location', 'south');
    end
    fprintf('%-13s %d windows, %d valid, median |error| of valid %.3f px/frame\n', names{k}, ...
        numel(starts), sum(valid), median(abs(est_speed(valid) - ref_speed(valid))));
end
