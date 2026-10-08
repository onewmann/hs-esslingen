function test_options()
%TEST_OPTIONS  Input checks and option handling.
    seed_rng(3);
    bp = simulate_pulsed_image(64, [2 1], 4);
    expect_error(@() estimate_velocity(bp, 'L', 0), 'L = 0');
    expect_error(@() estimate_velocity(bp, 'L', 2.5), 'non-integer L');
    expect_error(@() estimate_velocity(bp, 'PadFactor', 1.5), 'non-integer PadFactor');
    expect_error(@() estimate_velocity(bp, 'Bogus', 1), 'unknown option');
    expect_error(@() estimate_velocity(bp, 'L'), 'missing value');
    expect_error(@() estimate_velocity(zeros(8)), 'too small');
    expect_error(@() estimate_velocity(zeros(32, 32, 3)), 'not 2-D');

    e = estimate_velocity(bp, 'l', 4, 'minquality', 0.99);
    assert(e.valid == (e.quality >= 0.99), 'MinQuality and case-insensitive names');
    assert(~isfield(e, 'diag'), 'no diagnostics by default');
    e = estimate_velocity(bp, 'Diagnostics', true);
    assert(isfield(e, 'diag') && size(e.diag.spectrum, 1) == 128, 'diagnostics');
    e = estimate_velocity(ones(64));
    assert(~e.valid && e.quality == 0 && e.geba == 1, 'blank image');
    e1 = estimate_velocity(bp);
    e2 = estimate_velocity(uint8(255*(bp - min(bp(:)))/(max(bp(:)) - min(bp(:)))));
    assert(abs(angle_error(e1.angle_deg, e2.angle_deg)) < 1, 'integer input');
end

function expect_error(f, what)
    try
        f();
    catch
        return
    end
    error('test_options:noerror', 'expected an error for %s', what);
end
