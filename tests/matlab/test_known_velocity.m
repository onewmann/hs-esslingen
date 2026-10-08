function test_known_velocity()
%TEST_KNOWN_VELOCITY  Simulated motion with known v must be recovered.
    seed_rng(42);
    speeds = [0.7 3 12 30];
    angles = [15 80 125];
    for s = speeds
        for a = angles
            v = s*[cosd(a), -sind(a)];
            bp = simulate_pulsed_image(256, v, 4);
            e = estimate_velocity(bp, 'L', 4);
            assert(abs(angle_error(e.angle_deg, a)) < 0.5, ...
                sprintf('angle %.2f, expected %g (speed %g)', e.angle_deg, a, s));
            assert(abs(e.speed/s - 1) < 0.01, ...
                sprintf('speed %.3f, expected %g (angle %g)', e.speed, s, a));
            assert(e.valid, sprintf('not valid for speed %g, angle %g', s, a));
        end
    end

    % more exposures resolve slower motion
    bp = simulate_pulsed_image(256, [0.3 0], 20);
    e = estimate_velocity(bp, 'L', 20);
    assert(abs(e.speed/0.3 - 1) < 0.02, sprintf('L=20: speed %.3f, expected 0.3', e.speed));

    % real surface texture recorded with the Basler camera
    here = fileparts(mfilename('fullpath'));
    tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
    bp = simulate_pulsed_image(256, [5 -3], 4, 'Texture', tex);
    e = estimate_velocity(bp);
    assert(abs(angle_error(e.angle_deg, atan2d(3, 5))) < 0.5, 'real texture: angle');
    assert(abs(e.speed/hypot(5, 3) - 1) < 0.01, 'real texture: speed');

    % no motion must not produce a valid estimate
    for k = 1:3
        e = estimate_velocity(simulate_pulsed_image(256, [0 0], 4));
        assert(~e.valid, sprintf('static scene marked valid (quality %.2f)', e.quality));
    end
end
