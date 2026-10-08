function test_parity()
%TEST_PARITY  Same input, same result as the Python implementation.
%   tests/parity/cases.mat holds the inputs and tests/parity/expected.csv
%   the results of python/src/fftvel (see tests/parity/make_cases.py).
    here = fileparts(mfilename('fullpath'));
    pdir = fullfile(here, '..', 'parity');
    d = load(fullfile(pdir, 'cases.mat'));
    ex = dlmread(fullfile(pdir, 'expected.csv'), ',', 1, 0);
    for i = 1:size(ex, 1)
        bp = double(d.(sprintf('bp%02d', i)));
        e = estimate_velocity(bp, 'L', ex(i, 2), 'PadFactor', ex(i, 3));
        assert(abs(angle_error(e.angle_deg, ex(i, 4))) < 1e-6, sprintf('case %d: angle', i));
        assert(abs(e.speed/ex(i, 5) - 1) < 1e-7, sprintf('case %d: speed', i));
        assert(abs(e.geba/ex(i, 6) - 1) < 1e-7, sprintf('case %d: geba', i));
        assert(abs(e.quality - ex(i, 7)) < 1e-7, sprintf('case %d: quality', i));
        assert(e.valid == ex(i, 8), sprintf('case %d: valid', i));
    end
end
