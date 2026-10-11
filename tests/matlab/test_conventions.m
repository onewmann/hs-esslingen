function test_conventions()
%TEST_CONVENTIONS  Angle convention, image coordinates and helper functions.
    seed_rng(11);
    % v in image coordinates (x right, y down) -> expected angle (y up)
    cases = [ 2  0    0;
              0 -2   90;
              2 -2   45;
              2  2  135;
             -3  1   18.434948822922];   % -v = (3, -1): up and to the right
    for k = 1:size(cases, 1)
        v = cases(k, 1:2);
        bp = fftvel.simulate_pulsed_image(128, v, 4);
        e = fftvel.estimate_velocity(bp);
        assert(abs(fftvel.angle_difference(e.angle_deg, cases(k, 3))) < 0.5, ...
            sprintf('angle %.2f for v=[%g %g], expected %.2f', e.angle_deg, v, cases(k, 3)));
        assert(e.angle_deg >= 0 && e.angle_deg < 180, 'angle outside [0, 180)');
        vv = [e.vx e.vy];
        assert(min(norm(vv - v), norm(vv + v)) < 0.03*norm(v) + 0.02, ...
            sprintf('[vx vy] = [%.3f %.3f] for v=[%g %g]', vv, v));
    end

    assert(max(abs(fftvel.internal.hann_window(5) - [0; 0.5; 1; 0.5; 0])) < 1e-15, 'Hann window, 5 points');
    assert(isequal(fftvel.internal.hann_window(1), 1), 'Hann window, 1 point');

    img = reshape(1:24, 4, 6);
    assert(isequal(fftvel.internal.crop_square(img), img(:, 2:5)), 'crop_square wide image');
    assert(isequal(fftvel.internal.crop_square(img.'), img(:, 2:5).'), 'crop_square tall image');

    S = 3*ones(64);
    [M, j] = fftvel.internal.mean_projection(S, [0 33 90], 31);
    assert(max(abs(M(:) - 3)) < 1e-12, 'mean projection of a constant');
    assert(isequal(j, (-31:31).'), 'bin positions');

    y = -(((1:9) - 4.3).^2);
    [~, i] = max(y);
    assert(abs(i + fftvel.internal.parabolic_offset(y, i, false) - 4.3) < 1e-12, 'parabolic vertex');
    assert(fftvel.internal.parabolic_offset(y, 1, false) == 0, 'parabolic offset at the edge');

    % directions modulo 180 degrees, wrapped to [-90, 90) as in Python
    d = fftvel.angle_difference([179 1 10 100 90 0 270 -90], [1 179 0 0 0 90 0 0]);
    assert(isequal(d, [-2 2 10 -80 -90 -90 -90 -90]), 'angle_difference');
    assert(isequal(fftvel.angle_difference([30; 150], 135), [75; 15]), ...
        'angle_difference with a scalar');
    % integer inputs are converted to double: no saturation, no rounding
    d = [fftvel.angle_difference(uint8(10), uint8(20)), fftvel.angle_difference(30.5, int32(0))];
    assert(isa(d, 'double') && isequal(d, [-10 30.5]), 'angle_difference with integer input');
end
