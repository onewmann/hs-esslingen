function test_parity()
%TEST_PARITY  Same input, same result as the Python implementation.
%   tests/parity/cases.mat holds the inputs and tests/parity/expected.csv
%   the results of python/src/fftvel (see tests/parity/make_cases.py).
%   Columns 9-13 give the options of each case (0: default): Exposures,
%   Duty, ROI, ClipOutliers and stack (the input is an H x W x 1 x F uint8
%   stack for fftvel.pulsed_from_stack); columns 14-16 the results
%   threshold, quality_perp and isotropic; column 17 the option RejectBlur,
%   column 18 whether the image is passed as uint8 as stored; columns
%   19-22 the results quality_blur, blur_length, blur_like and saturated.
    here = fileparts(mfilename('fullpath'));
    pdir = fullfile(here, '..', 'parity');
    d = load(fullfile(pdir, 'cases.mat'));
    ex = dlmread(fullfile(pdir, 'expected.csv'), ',', 1, 0);
    assert(isequal(size(ex), [20 22]), 'expected.csv: unexpected size');
    for i = 1:size(ex, 1)
        L = ex(i, 2);
        K = [];
        if ex(i, 9) > 0
            K = ex(i, 9);
        end
        roi = [];
        if ex(i, 11) > 0
            roi = ex(i, 11);
        end
        x = d.(sprintf('bp%02d', i));
        if ex(i, 13)
            [bp, Ks] = fftvel.pulsed_from_stack(x, 'L', L, 'Exposures', K);
            if isempty(K)
                assert(Ks == 2*L + 1, sprintf('case %d: stack depth', i));
            else
                assert(Ks == K, sprintf('case %d: stack depth', i));
            end
        elseif ex(i, 18)
            assert(isa(x, 'uint8'), sprintf('case %d: expected a uint8 image', i));
            bp = x;              % integer input: saturated counts the 255 pixels
        else
            bp = double(x);
        end
        e = fftvel.estimate_velocity(bp, 'L', L, 'PadFactor', ex(i, 3), 'Exposures', K, ...
            'Duty', ex(i, 10), 'ROI', roi, 'ClipOutliers', ex(i, 12) ~= 0, 'RejectBlur', ex(i, 17) ~= 0);
        assert(abs(fftvel.angle_difference(e.angle_deg, ex(i, 4))) < 1e-9, sprintf('case %d: angle', i));
        assert(abs(e.speed/ex(i, 5) - 1) < 1e-9, sprintf('case %d: speed', i));
        assert(abs(e.geba/ex(i, 6) - 1) < 1e-9, sprintf('case %d: geba', i));
        assert(abs(e.quality - ex(i, 7)) < 1e-9, sprintf('case %d: quality', i));
        assert(e.valid == ex(i, 8), sprintf('case %d: valid', i));
        assert(abs(e.threshold - ex(i, 14)) < 1e-12, sprintf('case %d: threshold', i));
        assert(abs(e.quality_perp - ex(i, 15)) < 1e-9, sprintf('case %d: quality_perp', i));
        assert(e.isotropic == ex(i, 16), sprintf('case %d: isotropic', i));
        assert(abs(e.quality_blur - ex(i, 19)) < 1e-9, sprintf('case %d: quality_blur', i));
        assert(abs(e.blur_length/ex(i, 20) - 1) < 1e-9, sprintf('case %d: blur_length', i));
        assert(e.blur_like == ex(i, 21), sprintf('case %d: blur_like', i));
        assert(abs(e.saturated - ex(i, 22)) < 1e-12, sprintf('case %d: saturated', i));
    end
end
