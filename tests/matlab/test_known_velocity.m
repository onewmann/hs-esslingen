function test_known_velocity()
%TEST_KNOWN_VELOCITY  Simulated motion with known v must be recovered.
    seed_rng(42);
    speeds = [0.7 3 12 30];
    angles = [15 80 125];
    for s = speeds
        for a = angles
            v = s*[cosd(a), -sind(a)];
            bp = fftvel.simulate_pulsed_image(256, v, 4);
            e = fftvel.estimate_velocity(bp, 'L', 4);
            assert(abs(fftvel.angle_difference(e.angle_deg, a)) < 0.5, ...
                sprintf('angle %.2f, expected %g (speed %g)', e.angle_deg, a, s));
            assert(abs(e.speed/s - 1) < 0.01, ...
                sprintf('speed %.3f, expected %g (angle %g)', e.speed, s, a));
            assert(e.valid, sprintf('not valid for speed %g, angle %g', s, a));
        end
    end

    % more exposures resolve slower motion
    bp = fftvel.simulate_pulsed_image(256, [0.3 0], 20);
    e = fftvel.estimate_velocity(bp, 'L', 20);
    assert(abs(e.speed/0.3 - 1) < 0.02, sprintf('L=20: speed %.3f, expected 0.3', e.speed));

    % real surface texture recorded with the Basler camera
    here = fileparts(mfilename('fullpath'));
    tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
    bp = fftvel.simulate_pulsed_image(256, [5 -3], 4, 'Texture', tex);
    e = fftvel.estimate_velocity(bp);
    assert(abs(fftvel.angle_difference(e.angle_deg, atan2d(3, 5))) < 0.5, 'real texture: angle');
    assert(abs(e.speed/hypot(5, 3) - 1) < 0.01, 'real texture: speed');

    % motion below the measurable range must not produce a valid estimate
    for c = [1 0.9; 1 1.0; 3 0.45; 4 0.3].'
        bp = fftvel.simulate_pulsed_image(256, c(2)*[cosd(37), -sind(37)], c(1));
        e = fftvel.estimate_velocity(bp, 'L', c(1));
        assert(~e.in_range && ~e.valid, sprintf('L=%d, speed %g marked valid', c(1), c(2)));
    end

    % no motion must not produce a valid estimate, also on small images.
    % At 32 px about 1 in 1000 such images passes, so these inputs are fixed
    % (saved by tests/parity/make_cases.py, the same as in the Python tests)
    fx = load(fullfile(here, '..', 'parity', 'fixed_inputs.mat'));
    for k = 1:3
        for name = {sprintf('static32_%d', k), sprintf('noise32_%d', k)}
            e = fftvel.estimate_velocity(double(fx.(name{1})));
            assert(~e.valid && e.quality < e.threshold - 0.05, ...
                sprintf('%s marked valid (quality %.2f)', name{1}, e.quality));
        end
    end
    for n = [64 128 256]
        for k = 1:3
            e = fftvel.estimate_velocity(fftvel.simulate_pulsed_image(n, [0 0], 4));
            assert(~e.valid, sprintf('static scene %d px marked valid (quality %.2f)', n, e.quality));
            e = fftvel.estimate_velocity(randn(n));
            assert(~e.valid, sprintf('white noise %d px marked valid (quality %.2f)', n, e.quality));
        end
    end
end
