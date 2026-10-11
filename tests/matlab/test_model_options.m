function test_model_options()
%TEST_MODEL_OPTIONS  Exposures, Duty, ROI, ClipOutliers, isotropy check, frame stacks.
%   The same checks are in tests/python/test_model_options.py. Random
%   inputs are used only where every realisation passes; the rest is
%   deterministic or loaded from tests/parity/fixed_inputs.mat.
    seed_rng(7);
    here = fileparts(mfilename('fullpath'));
    tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
    fx = load(fullfile(here, '..', 'parity', 'fixed_inputs.mat'));

    % even number of exposures; Exposures wins over L
    for K = [8 10]
        for c = [0.6 25; 5 110].'
            bp = simulate_pulsed_image(256, c(1)*[cosd(c(2)), -sind(c(2))], 4, 'Exposures', K);
            e = estimate_velocity(bp, 'Exposures', K);
            check(e, c(1), c(2), sprintf('Exposures %d', K));
            assert(e.exposures == K && e.L == 4, 'exposures field');
            e1 = estimate_velocity(bp, 'L', 1, 'Exposures', K);
            assert(e1.speed == e.speed, 'Exposures wins over L');
        end
    end
    [bp, fr] = simulate_pulsed_image(64, [2.5 -1], 4, 'Exposures', 6, 'Duty', 0.3);
    assert(size(fr, 3) == 6 && max(abs(bp(:) - reshape(mean(fr, 3), [], 1))) < 1e-12, ...
        'simulated image is the mean of the exposures');
    seed_rng(2);
    a = simulate_pulsed_image(64, [2.5 -1], 4);
    seed_rng(2);
    assert(isequal(a, simulate_pulsed_image(64, [2.5 -1], 4, 'Exposures', 9)), 'Exposures 2L+1');

    % motion blur within each exposure, real surface texture
    for duty = [0.2 0.5]
        for c = [34.6 10.1; 12 120].'
            bp = simulate_pulsed_image(256, c(1)*[cosd(c(2)), -sind(c(2))], 4, 'Texture', tex, 'Duty', duty);
            check(estimate_velocity(bp, 'Duty', duty), c(1), c(2), sprintf('Duty %g', duty));
        end
    end
    bp = simulate_pulsed_image(128, [3 -2], 4);
    assert(isequal(estimate_velocity(bp, 'Duty', 0), estimate_velocity(bp)), 'Duty 0');

    % analysed square
    bp = simulate_pulsed_image(300, [3 -1], 4);
    bp = bp(:, 1:280);
    e = estimate_velocity(bp, 'ROI', 128);
    assert(e.n == 128 && e.npad == 256, 'ROI size');
    check(e, hypot(3, 1), atan2d(1, 3), 'ROI');
    assert(isequal(e, estimate_velocity(crop_square(bp, 128))), 'ROI is the central square');
    assert(isequal(estimate_velocity(bp, 'ROI', 1000), estimate_velocity(bp)), 'ROI larger than the image');

    % two static hot pixels: reported as motion at their separation unless clipped
    bp = simulate_pulsed_image(96, [1.5 -0.5], 4);
    clean = min(max(round(100 + 15*bp/std(bp(:))), 0), 255);
    x = clean;
    x(31, 31) = 255;
    x(51, 61) = 255;
    e = estimate_velocity(x);
    assert(e.valid && abs(e.speed/hypot(20, 30) - 1) < 0.01, 'hot pixels without ClipOutliers');
    check(estimate_velocity(x, 'ClipOutliers', true), hypot(1.5, 0.5), atan2d(0.5, 1.5), 'ClipOutliers');
    assert(isequal(estimate_velocity(clean, 'ClipOutliers', true), estimate_velocity(clean)), ...
        'ClipOutliers changes nothing on a clean image');

    % isotropic spectra: defocused static texture and a square repeating texture
    [xx, yy] = meshgrid(-6:6);
    disk = double(xx.^2 + yy.^2 <= 2.5^2);
    img = conv2(tex(301:568, 301:568), disk/sum(disk(:)), 'valid');
    e = estimate_velocity(img);
    assert(e.quality >= e.threshold && e.in_range, 'defocus: passed before the isotropy check');
    assert(e.isotropic && ~e.valid, 'defocused static texture marked valid');
    e = estimate_velocity(double(fx.tile16));
    assert(e.quality >= e.threshold && e.in_range && abs(e.speed - 16) < 0.1, 'tile: fit');
    assert(e.isotropic && ~e.valid, 'repeating texture marked valid');
    e = estimate_velocity(simulate_pulsed_image(128, [3 -2], 4));
    assert(~e.isotropic && e.quality - e.quality_perp > 0.15, 'motion marked isotropic');

    % zero padding does not lower the threshold
    for n = [32 64]
        e2 = estimate_velocity(double(fx.(sprintf('noise%d_1', n))));
        for pad = [6 8]
            for k = 1:3
                e = estimate_velocity(double(fx.(sprintf('noise%d_%d', n, k))), 'PadFactor', pad);
                assert(~e.valid && e.threshold > e2.threshold - 0.02, ...
                    sprintf('white noise %d px, PadFactor %d marked valid', n, pad));
            end
        end
    end
    % beyond the window (PadFactor > 2) the model is flat and must not score rounding errors
    e = estimate_velocity(double(fx.noise32_1), 'PadFactor', 8, 'Diagnostics', true);
    flat = e.diag.speed_grid*max(abs(cosd(e.angle_deg)), abs(sind(e.angle_deg))) >= e.n - 1;
    assert(any(flat) && all(e.diag.speed_score(flat) == 0), 'flat model scores');

    % frame stacks
    [~, fr] = simulate_pulsed_image(64, [3 -2], 5);           % 11 frames
    u8 = uint8(round(255*(fr - min(fr(:)))/(max(fr(:)) - min(fr(:)))));
    [bp, K] = pulsed_from_stack(u8);
    assert(K == 11 && isa(bp, 'double') && isequal(size(bp), [64 64]), 'stack of 11');
    assert(isequal(bp, mean(double(u8), 3)), 'stack mean');
    camera = reshape(u8(:, :, 1:10), 64, 64, 1, 10);         % H x W x 1 x F, as from getdata
    [bp, K] = pulsed_from_stack(camera, 'L', 4);
    assert(K == 9 && isequal(bp, mean(double(u8(:, :, 1:9)), 3)), 'camera stack, L');
    [bp, K] = pulsed_from_stack(camera, 'L', 2, 'Exposures', 8);
    assert(K == 8 && isequal(bp, mean(double(u8(:, :, 1:8)), 3)), 'camera stack, Exposures');
    expect_error(@() pulsed_from_stack(camera), 'pulsed_from_stack:even');
    expect_error(@() pulsed_from_stack(camera, 'L', 5), 'pulsed_from_stack:frames');
    expect_error(@() pulsed_from_stack(u8(:, :, 1)), 'pulsed_from_stack:dims');
    expect_error(@() pulsed_from_stack(u8, 'L', '2'), 'pulsed_from_stack:L');
    expect_error(@() pulsed_from_stack(u8, 'Exposures', 1), 'pulsed_from_stack:Exposures');

    % colour texture for the simulation
    rgb = cat(3, tex, 0.5*tex, tex.^2);
    grey = 0.2989*rgb(:, :, 1) + 0.5870*rgb(:, :, 2) + 0.1140*rgb(:, :, 3);
    a = simulate_pulsed_image(64, [2 1], 2, 'Texture', rgb);
    assert(isequal(a, simulate_pulsed_image(64, [2 1], 2, 'Texture', grey)), 'RGB texture');
    assert(isequal(a, simulate_pulsed_image(64, [2 1], 2, 'Texture', cat(3, rgb, ones(size(tex))))), ...
        'RGBA texture');
    expect_error(@() simulate_pulsed_image(64, [2 1], 2, 'Texture', rgb(:, :, 1:2)), ...
        'simulate_pulsed_image:texture');

    % the simulator checks its options as the estimator does
    expect_error(@() simulate_pulsed_image(64, [2 1], '2'), 'simulate_pulsed_image:L');
    expect_error(@() simulate_pulsed_image(64, [2 1], 0), 'simulate_pulsed_image:L');
    expect_error(@() simulate_pulsed_image(64, [2 1], 2, 'Exposures', 1), 'simulate_pulsed_image:Exposures');
    expect_error(@() simulate_pulsed_image(64, [2 1], 2, 'Exposures', 6.5), 'simulate_pulsed_image:Exposures');
    expect_error(@() simulate_pulsed_image(64, [2 1], 2, 'Duty', 1), 'simulate_pulsed_image:Duty');
    expect_error(@() simulate_pulsed_image(64, [2 1], 2, 'Duty', '0.2'), 'simulate_pulsed_image:Duty');
end

function check(e, speed, angle, what)
    assert(abs(angle_error(e.angle_deg, angle)) < 0.5, ...
        sprintf('%s: angle %.2f, expected %g', what, e.angle_deg, angle));
    assert(abs(e.speed/speed - 1) < 0.01, sprintf('%s: speed %.4f, expected %g', what, e.speed, speed));
    assert(e.valid, sprintf('%s: not valid', what));
end

function expect_error(f, id)
    try
        f();
    catch err
        assert(strcmp(err.identifier, id), sprintf('expected %s, got %s: %s', id, err.identifier, err.message));
        return
    end
    error('test_model_options:noerror', 'expected the error %s', id);
end
