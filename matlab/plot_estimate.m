function fig = plot_estimate(bp, est, frame)
%PLOT_ESTIMATE  Show every stage of the estimate in one figure.
%   PLOT_ESTIMATE(BP, EST) needs EST from ESTIMATE_VELOCITY(..., 'Diagnostics',
%   true). PLOT_ESTIMATE(BP, EST, FRAME) also shows a single exposure.
%   Works in MATLAB and Octave.

    if ~isfield(est, 'diag')
        error('plot_estimate:diag', 'Call estimate_velocity with ''Diagnostics'', true.');
    end
    d = est.diag;
    has_frame = nargin >= 3 && ~isempty(frame);
    fig = figure('Name', 'Pulsed-exposure velocity estimate', 'Color', 'w');
    colormap(gray);
    if has_frame
        slot = struct('frame', 1, 'bp', 2, 'spec', 3, 'wr', 4, 'prof', 5, 'text', 6);
    else
        slot = struct('frame', 0, 'bp', 1, 'spec', 2, 'wr', 4, 'prof', [5 6], 'text', 3);
    end

    if has_frame
        subplot(2, 3, slot.frame);
        imagesc(crop_square(double(frame))); axis image; axis off;
        title('Single exposure');
    end

    subplot(2, 3, slot.bp);
    imagesc(crop_square(double(bp))); axis image; axis off;
    title(sprintf('Pulsed image (%d exposures)', 2*est.L + 1));

    subplot(2, 3, slot.spec);
    S = d.spectrum;
    npad = size(S, 1);
    c = floor(npad/2) + 1;
    ax = ((1:npad) - c)/npad;
    imagesc(ax, -ax, S); axis image; axis xy;
    hold on;
    t = est.angle_deg*pi/180;
    plot(0.5*[-cos(t) cos(t)], 0.5*[-sin(t) sin(t)], 'r-', 'LineWidth', 1);
    hold off;
    xlabel('f_x (cycles/px)'); ylabel('f_y (cycles/px)');
    title('Log spectrum, direction of v');

    subplot(2, 3, slot.wr);
    plot(d.theta, d.wr, 'k-', 'LineWidth', 1);
    hold on;
    yl = get(gca, 'YLim');
    plot([est.angle_deg est.angle_deg], yl, 'r--');
    hold off;
    xlim([0 180]);
    xlabel('\theta (deg)'); ylabel('WR(\theta)');
    title(sprintf('Ripple profile, GEBA = %.2f', est.geba));

    subplot(2, 3, slot.prof);
    plot(d.profile_j, d.profile, 'k-');
    hold on;
    plot(d.profile_j, d.model, 'r-', 'LineWidth', 1);
    hold off;
    xlabel('distance from centre (bins)'); ylabel('mean log spectrum');
    legend('measured', 'Dirichlet model', 'Location', 'best');
    title(sprintf('Projection along v, fit quality %.2f', est.quality));

    subplot(2, 3, slot.text);
    axis off;
    if est.valid
        state = 'valid';
    else
        state = 'not valid';
    end
    txt = {sprintf('direction  %.1f deg', est.angle_deg), ...
           sprintf('speed      %.2f px/step', est.speed), ...
           sprintf('v (image)  [%.2f, %.2f]', est.vx, est.vy), ...
           sprintf('GEBA       %.2f', est.geba), ...
           sprintf('quality    %.2f (%s)', est.quality, state)};
    for k = 1:numel(txt)
        text(0, 1 - 0.18*k, txt{k}, 'FontName', 'Courier', 'FontSize', 10);
    end
end
