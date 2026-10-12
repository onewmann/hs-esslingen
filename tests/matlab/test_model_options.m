function test_model_options()
%TEST_MODEL_OPTIONS  Exposures, Duty, ROI, ClipOutliers, isotropy and blur checks, saturation,
%   frame stacks.
%   The same checks are in tests/python/test_model_options.py. MATLAB and
%   Octave draw different random numbers from the same seed, so inputs that
%   fail the checks for some random draws are loaded from
%   tests/parity/fixed_inputs.mat (the same images as in the Python tests);
%   the rest is deterministic or passed for each of 200 seeds in Octave.
    seed_rng(7);
    here = fileparts(mfilename('fullpath'));
    tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
    fx = load(fullfile(here, '..', 'parity', 'fixed_inputs.mat'));

    % even number of exposures; Exposures wins over L. The image of 8 exposures at 0.6 px per step
    % is stored: about 1 in 200 random textures gives an angle more than 0.5 deg off
    for K = [8 10]
        for c = [0.6 25; 5 110].'
            if K == 8 && c(1) == 0.6
                bp = double(fx.exposures8);
            else
                bp = fftvel.simulate_pulsed_image(256, c(1)*[cosd(c(2)), -sind(c(2))], 4, 'Exposures', K);
            end
            e = fftvel.estimate_velocity(bp, 'Exposures', K);
            check(e, c(1), c(2), sprintf('Exposures %d', K));
            assert(e.exposures == K && e.L == 4, 'exposures field');
            e1 = fftvel.estimate_velocity(bp, 'L', 1, 'Exposures', K);
            assert(e1.speed == e.speed, 'Exposures wins over L');
        end
    end
    [bp, fr] = fftvel.simulate_pulsed_image(64, [2.5 -1], 4, 'Exposures', 6, 'Duty', 0.3);
    assert(size(fr, 3) == 6 && max(abs(bp(:) - reshape(mean(fr, 3), [], 1))) < 1e-12, ...
        'simulated image is the mean of the exposures');
    seed_rng(2);
    a = fftvel.simulate_pulsed_image(64, [2.5 -1], 4);
    seed_rng(2);
    assert(isequal(a, fftvel.simulate_pulsed_image(64, [2.5 -1], 4, 'Exposures', 9)), 'Exposures 2L+1');

    % motion blur within each exposure, real surface texture
    for duty = [0.2 0.5]
        for c = [34.6 10.1; 12 120].'
            bp = fftvel.simulate_pulsed_image(256, c(1)*[cosd(c(2)), -sind(c(2))], 4, ...
                'Texture', tex, 'Duty', duty);
            check(fftvel.estimate_velocity(bp, 'Duty', duty), c(1), c(2), sprintf('Duty %g', duty));
        end
    end
    bp = fftvel.simulate_pulsed_image(128, [3 -2], 4);
    assert(isequal(fftvel.estimate_velocity(bp, 'Duty', 0), fftvel.estimate_velocity(bp)), 'Duty 0');

    % analysed square
    bp = fftvel.simulate_pulsed_image(300, [3 -1], 4);
    bp = bp(:, 1:280);
    e = fftvel.estimate_velocity(bp, 'ROI', 128);
    assert(e.n == 128 && e.npad == 256, 'ROI size');
    check(e, hypot(3, 1), atan2d(1, 3), 'ROI');
    assert(isequal(e, fftvel.estimate_velocity(fftvel.internal.crop_square(bp, 128))), ...
        'ROI is the central square');
    assert(isequal(fftvel.estimate_velocity(bp, 'ROI', 1000), fftvel.estimate_velocity(bp)), ...
        'ROI larger than the image');

    % two static hot pixels: reported as motion at their separation unless clipped. The 8-bit
    % texture moving by [1.5 -0.5] is stored: about 4 % of the random textures are off by more
    % than 1 % or 0.5 deg after clipping
    clean = double(fx.hot96);
    x = clean;
    x(31, 31) = 255;
    x(51, 61) = 255;
    e = fftvel.estimate_velocity(x);
    assert(e.valid && abs(e.speed/hypot(20, 30) - 1) < 0.01, 'hot pixels without ClipOutliers');
    check(fftvel.estimate_velocity(x, 'ClipOutliers', true), hypot(1.5, 0.5), atan2d(0.5, 1.5), ...
        'ClipOutliers');
    assert(isequal(fftvel.estimate_velocity(clean, 'ClipOutliers', true), ...
        fftvel.estimate_velocity(clean)), ...
        'ClipOutliers changes nothing on a clean image');
    % the clipped pixels and all pixels within CLIP_RADIUS of them are set to the median
    cst = fftvel.internal.constants();
    y = x;
    [cols, rows] = meshgrid(1:96, 1:96);
    for rc = [31 31; 51 61].'
        y((rows - rc(1)).^2 + (cols - rc(2)).^2 <= cst.CLIP_RADIUS^2) = median(x(:));
    end
    assert(isequal(fftvel.internal.log_spectrum(x, 2, [], true), fftvel.internal.log_spectrum(y)), ...
        'ClipOutliers takes the rim along');
    m = false(7, 9);
    m(4, 5) = true;
    d = fftvel.internal.dilate(m, 2);
    assert(nnz(d) == 13 && d(4, 3) && d(2, 5) && d(3, 4) && ~d(2, 4), 'dilate: disc of radius 2');
    m(4, 5) = false;
    m(1, 1) = true;
    assert(nnz(fftvel.internal.dilate(m, 2)) == 6, 'dilate: corner, nothing wraps');
    assert(isequal(fftvel.internal.dilate(m, 0), m) && ~any(any(fftvel.internal.dilate(false(5), 3))), ...
        'dilate: radius 0 and empty mask');

    % blur check: the box template of a blur of K*BLUR_STEP px is the Dirichlet template of K
    % exposures BLUR_STEP apart
    j = (2:229).';
    for K = [2 9 40]
        a = fftvel.internal.box_template(j, K*cst.BLUR_STEP, 512, 256, 30);
        b = fftvel.internal.dirichlet_template(j, cst.BLUR_STEP, 4, 512, 256, 30, K);
        assert(max(abs(a - b)) < 1e-12, sprintf('box template, K = %d', K));
    end
    % one continuous exposure looks like K flashes d/K apart. The inputs are fixed (parity cases 16
    % and 17: a single exposure blurred over 6 px with 5 % noise, no flashes), because a few noise
    % realisations give a wrong direction and then a wrong blur length
    pc = load(fullfile(here, '..', 'parity', 'cases.mat'), 'bp16', 'bp17');
    for bp = {double(pc.bp16), double(pc.bp17)}
        e = fftvel.estimate_velocity(bp{1});
        assert(e.valid && e.blur_like && abs(e.blur_length/6 - 1) < 0.05 && abs(9*e.speed/6 - 1) < 0.05, ...
            'blur of 6 px: valid and blur_like');
        r = fftvel.estimate_velocity(bp{1}, 'RejectBlur', true);
        assert(~r.valid, 'blur of 6 px: valid with RejectBlur');
        r.valid = e.valid;
        assert(isequal(r, e), 'RejectBlur changes only valid');
    end
    bp = fftvel.simulate_pulsed_image(128, [7 3], 4);
    e = fftvel.estimate_velocity(bp);
    assert(~e.blur_like && e.quality - e.quality_blur > 0.3, 'fast flashes: not blur_like');
    assert(isequal(fftvel.estimate_velocity(bp, 'RejectBlur', true), e), 'fast flashes: RejectBlur');

    % saturated: share of pixels at the maximum (same rules as tests/python/test_model_options.py)
    x8 = repmat(uint8(7), 20, 20);
    x8(1, 1:2) = 255;
    assert(fftvel.internal.saturated_share(x8) == 2/400, 'saturated: uint8 maximum');
    x8(1, 1:2) = 254;
    assert(fftvel.internal.saturated_share(x8) == 2/400, 'saturated: below the maximum, 0.5 %');
    x8(1, 2) = 7;
    assert(fftvel.internal.saturated_share(x8) == 0 && fftvel.internal.saturated_share(double(x8)) == 0, ...
        'saturated: 0.25 % below the maximum');
    x8(1, 2) = 254;
    assert(fftvel.internal.saturated_share(single(x8)) == 2/400, 'saturated: single');
    big = repmat(uint16(50), 100, 100);
    big(1, 1) = 65535;
    assert(fftvel.internal.saturated_share(big) == 1e-4 && ...
        fftvel.internal.saturated_share(double(big)) == 0, 'saturated: uint16');
    assert(fftvel.internal.saturated_share(repmat(uint8(255), 8, 8)) == 1 && ...
        fftvel.internal.saturated_share(repmat(255, 8, 8)) == 0 && ...
        fftvel.internal.saturated_share(repmat(int16(7), 8, 8)) == 0, 'saturated: constant images');
    i16 = zeros(10, 'int16');
    i16(1, 1) = 32767;
    lg = false(10);
    lg(1, 1) = true;
    assert(fftvel.internal.saturated_share(i16) == 0.01 && fftvel.internal.saturated_share(lg) == 0.01, ...
        'saturated: int16 and logical');
    bp = fftvel.simulate_pulsed_image(160, [3 -1], 4);
    x8 = uint8(min(max(round(100 + 15*bp/std(bp(:))), 0), 255));
    x8(1:16, 1:16) = 255;                       % outside the central 128 px
    e1 = fftvel.estimate_velocity(x8, 'ROI', 128);
    e2 = fftvel.estimate_velocity(x8);
    assert(e1.saturated == 0 && e2.saturated == 256/160^2, 'saturated: analysed square only');
    x8(61:70, 61:70) = 255;                     % inside
    e = fftvel.estimate_velocity(x8, 'ROI', 128, 'ClipOutliers', true);
    e1 = fftvel.estimate_velocity(x8, 'ROI', 128);
    e2 = fftvel.estimate_velocity(double(x8), 'ROI', 128);
    assert(e.saturated == 100/128^2 && e1.saturated == e.saturated, 'saturated: before clipping');
    assert(e2.saturated == 100/128^2, 'saturated: double');

    % isotropic spectra: defocused static texture and a square repeating texture
    [xx, yy] = meshgrid(-6:6);
    disk = double(xx.^2 + yy.^2 <= 2.5^2);
    img = conv2(tex(301:568, 301:568), disk/sum(disk(:)), 'valid');
    e = fftvel.estimate_velocity(img);
    assert(e.quality >= e.threshold && e.in_range, 'defocus: passed before the isotropy check');
    assert(e.isotropic && ~e.valid, 'defocused static texture marked valid');
    e = fftvel.estimate_velocity(double(fx.tile16));
    assert(e.quality >= e.threshold && e.in_range && abs(e.speed - 16) < 0.1, 'tile: fit');
    assert(e.isotropic && ~e.valid, 'repeating texture marked valid');
    e = fftvel.estimate_velocity(fftvel.simulate_pulsed_image(128, [3 -2], 4));
    assert(~e.isotropic && e.quality - e.quality_perp > 0.15, 'motion marked isotropic');

    % zero padding does not lower the threshold
    for n = [32 64]
        e2 = fftvel.estimate_velocity(double(fx.(sprintf('noise%d_1', n))));
        for pad = [6 8]
            for k = 1:3
                e = fftvel.estimate_velocity(double(fx.(sprintf('noise%d_%d', n, k))), 'PadFactor', pad);
                assert(~e.valid && e.threshold > e2.threshold - 0.02, ...
                    sprintf('white noise %d px, PadFactor %d marked valid', n, pad));
            end
        end
    end
    % beyond the window (PadFactor > 2) the model is flat and must not score rounding errors
    e = fftvel.estimate_velocity(double(fx.noise32_1), 'PadFactor', 8, 'Diagnostics', true);
    flat = e.diag.speed_grid*max(abs(cosd(e.angle_deg)), abs(sind(e.angle_deg))) >= e.n - 1;
    assert(any(flat) && all(e.diag.speed_score(flat) == 0), 'flat model scores');

    % frame stacks
    [~, fr] = fftvel.simulate_pulsed_image(64, [3 -2], 5);           % 11 frames
    u8 = uint8(round(255*(fr - min(fr(:)))/(max(fr(:)) - min(fr(:)))));
    [bp, K] = fftvel.pulsed_from_stack(u8);
    assert(K == 11 && isa(bp, 'double') && isequal(size(bp), [64 64]), 'stack of 11');
    assert(isequal(bp, mean(double(u8), 3)), 'stack mean');
    camera = reshape(u8(:, :, 1:10), 64, 64, 1, 10);         % H x W x 1 x F, as from getdata
    [bp, K] = fftvel.pulsed_from_stack(camera, 'L', 4);
    assert(K == 9 && isequal(bp, mean(double(u8(:, :, 1:9)), 3)), 'camera stack, L');
    [bp, K] = fftvel.pulsed_from_stack(camera, 'L', 2, 'Exposures', 8);
    assert(K == 8 && isequal(bp, mean(double(u8(:, :, 1:8)), 3)), 'camera stack, Exposures');
    expect_error(@() fftvel.pulsed_from_stack(camera), 'pulsed_from_stack:even');
    expect_error(@() fftvel.pulsed_from_stack(camera, 'L', 5), 'pulsed_from_stack:frames');
    expect_error(@() fftvel.pulsed_from_stack(u8(:, :, 1)), 'pulsed_from_stack:dims');
    expect_error(@() fftvel.pulsed_from_stack(u8, 'L', '2'), 'pulsed_from_stack:L');
    expect_error(@() fftvel.pulsed_from_stack(u8, 'Exposures', 1), 'pulsed_from_stack:Exposures');

    % colour texture for the simulation
    rgb = cat(3, tex, 0.5*tex, tex.^2);
    grey = 0.2989*rgb(:, :, 1) + 0.5870*rgb(:, :, 2) + 0.1140*rgb(:, :, 3);
    a = fftvel.simulate_pulsed_image(64, [2 1], 2, 'Texture', rgb);
    assert(isequal(a, fftvel.simulate_pulsed_image(64, [2 1], 2, 'Texture', grey)), 'RGB texture');
    assert(isequal(a, fftvel.simulate_pulsed_image(64, [2 1], 2, ...
        'Texture', cat(3, rgb, ones(size(tex))))), ...
        'RGBA texture');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 2, 'Texture', rgb(:, :, 1:2)), ...
        'simulate_pulsed_image:texture');

    % the simulator checks its options as the estimator does
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], '2'), 'simulate_pulsed_image:L');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 0), 'simulate_pulsed_image:L');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 2, 'Exposures', 1), ...
        'simulate_pulsed_image:Exposures');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 2, 'Exposures', 6.5), ...
        'simulate_pulsed_image:Exposures');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 2, 'Duty', 1), 'simulate_pulsed_image:Duty');
    expect_error(@() fftvel.simulate_pulsed_image(64, [2 1], 2, 'Duty', '0.2'), 'simulate_pulsed_image:Duty');
end

function check(e, speed, angle, what)
    assert(abs(fftvel.angle_difference(e.angle_deg, angle)) < 0.5, ...
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
