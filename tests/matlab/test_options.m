function test_options()
%TEST_OPTIONS  Input checks and option handling.
    seed_rng(3);
    bp = fftvel.simulate_pulsed_image(64, [2 1], 4);
    expect_error(@() fftvel.estimate_velocity(bp, 'Bogus', 1), 'unknown option', 'parse_options:unknown');
    expect_error(@() fftvel.estimate_velocity(bp, 'L'), 'missing value', 'parse_options:pairs');
    expect_error(@() fftvel.estimate_velocity(bp, 5, 3), 'option name not text', 'parse_options:name');
    expect_error(@() fftvel.estimate_velocity(zeros(8)), 'too small', 'log_spectrum:size');
    expect_error(@() fftvel.estimate_velocity(zeros(32, 32, 3)), 'not 2-D', 'crop_square:dims');
    bad = bp;
    bad(3, 4) = NaN;
    expect_error(@() fftvel.estimate_velocity(bad), 'NaN pixel', 'log_spectrum:finite');

    % every option value is checked; the same rules as in tests/python/test_core.py
    noise = randn(32);
    bad_options = {'L', 0; 'L', 2.5; 'L', '4'; 'L', true; 'L', Inf; 'L', NaN; 'L', [4 4]; ...
        'L', []; 'L', complex(4, 1); ...
        'PadFactor', 1.5; 'PadFactor', '2'; 'PadFactor', 0; 'PadFactor', false; ...
        'MinQuality', NaN; 'MinQuality', []; 'MinQuality', '0.8'; 'MinQuality', 1.5; ...
        'MinQuality', -0.1; 'MinQuality', Inf; ...
        'Exposures', 1; 'Exposures', 8.5; 'Exposures', '8'; 'Exposures', true; ...
        'Duty', 1; 'Duty', -0.1; 'Duty', NaN; 'Duty', '0.2'; 'Duty', []; ...
        'ROI', 8; 'ROI', 100.5; 'ROI', '128'; 'ROI', Inf; ...
        'ClipOutliers', 2; 'ClipOutliers', 'yes'; 'ClipOutliers', []; ...
        'RejectBlur', 2; 'RejectBlur', 'yes'; 'RejectBlur', []; ...
        'Diagnostics', 'no'; 'Diagnostics', 2};
    for k = 1:size(bad_options, 1)
        name = bad_options{k, 1};
        value = bad_options{k, 2};
        expect_error(@() fftvel.estimate_velocity(noise, name, value), ...
            sprintf('%s, bad value %d', name, k), ['estimate_velocity:' name]);
    end
    msg = '';
    try
        fftvel.estimate_velocity(noise, 'L', '4');
    catch err
        msg = err.message;
    end
    assert(~isempty(strfind(msg, 'text ''4''')), 'message for text input');
    expect_error(@() fftvel.internal.log_spectrum(noise, '2'), 'log_spectrum PadFactor', 'log_spectrum:pad');
    expect_error(@() fftvel.internal.log_spectrum(noise, 2, '64'), 'log_spectrum ROI', 'log_spectrum:roi');
    expect_error(@() fftvel.internal.log_spectrum(noise, 2, [], 2), 'log_spectrum ClipOutliers', ...
        'log_spectrum:clip');

    e1 = fftvel.estimate_velocity(bp, 'L', int32(4), 'PadFactor', uint8(2), 'MinQuality', single(0.5), ...
        'Exposures', int64(9), 'Duty', single(0), 'ROI', uint16(64), 'ClipOutliers', false, ...
        'RejectBlur', false);
    e2 = fftvel.estimate_velocity(bp, 'L', 4, 'MinQuality', 0.5);
    assert(isequal(e1, e2) && e1.exposures == 9 && e1.L == 4, 'integer-class options');
    e3 = fftvel.estimate_velocity(bp, 'Exposures', [], 'ROI', [], 'ClipOutliers', 0, 'Duty', 0);
    assert(isequal(e3, fftvel.estimate_velocity(bp)), 'empty means default');

    e = fftvel.estimate_velocity(bp, 'l', 4, 'minquality', 0.99);
    assert(e.valid == (e.quality >= max(0.99, e.threshold) && e.in_range && ~e.isotropic), ...
        'MinQuality and case-insensitive names');
    assert(~isfield(e, 'diag'), 'no diagnostics by default');
    e = fftvel.estimate_velocity(bp, 'Diagnostics', true);
    assert(isfield(e, 'diag') && size(e.diag.spectrum, 1) == 128, 'diagnostics');
    e = fftvel.estimate_velocity(ones(64));
    assert(~e.valid && e.quality == 0 && e.geba == 1, 'blank image');
    e1 = fftvel.estimate_velocity(bp);
    e2 = fftvel.estimate_velocity(uint8(255*(bp - min(bp(:)))/(max(bp(:)) - min(bp(:)))));
    assert(abs(fftvel.angle_difference(e1.angle_deg, e2.angle_deg)) < 1, 'integer input');
end

function expect_error(f, what, id)
    try
        f();
    catch err
        if ~strcmp(err.identifier, id)
            error('test_options:id', 'wrong error for %s: %s (%s)', what, err.identifier, err.message);
        end
        return
    end
    error('test_options:noerror', 'expected an error for %s', what);
end
