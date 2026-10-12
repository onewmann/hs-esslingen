function test_phase_correlation()
%TEST_PHASE_CORRELATION  Reference shift between two frames.
%   tests/parity/phase_cases.mat holds image pairs and the results of the
%   Python implementation (see tests/parity/make_phase_cases.py);
%   tests/python/test_phase_correlation.py checks Python against the same
%   numbers. The random texture of the sign, near-half and single-row checks
%   is the one the Python tests draw, loaded from
%   tests/parity/fixed_inputs.mat: MATLAB and Octave draw different numbers
%   from the same seed, and about 4 % of the textures miss a tolerance.
    here = fileparts(mfilename('fullpath'));

    % simulated frames one and two steps apart
    seed_rng(5);
    [~, fr] = fftvel.simulate_pulsed_image(128, [3 -2], 1);
    [dx, dy] = fftvel.phase_correlation(fr(:, :, 1), fr(:, :, 2));
    assert(abs(dx - 3) < 0.02 && abs(dy + 2) < 0.02, sprintf('integer shift: %.4f %.4f', dx, dy));
    [~, fr] = fftvel.simulate_pulsed_image(128, [-1.5 0.75], 1);
    [dx, dy] = fftvel.phase_correlation(fr(:, :, 1), fr(:, :, 3));
    assert(abs(dx + 3) < 0.02 && abs(dy - 1.5) < 0.02, sprintf('sub-pixel shift: %.4f %.4f', dx, dy));

    % parity with Python
    d = load(fullfile(here, '..', 'parity', 'phase_cases.mat'));
    assert(isequal(size(d.expected), [10 3]) && isequal(size(d.shift), [10 2]), ...
        'phase_cases.mat: unexpected size');
    for k = 1:size(d.expected, 1)
        [dx, dy, peak] = fftvel.phase_correlation(d.(sprintf('a%02d', k)), d.(sprintf('b%02d', k)));
        assert(max(abs([dx dy peak] - d.expected(k, :))) < 1e-9, ...
            sprintf('case %d: %.12g %.12g %.12g, Python %.12g %.12g %.12g', k, dx, dy, peak, ...
            d.expected(k, :)));
        if k ~= 8 && ~any(isnan(d.shift(k, :)))     % not the case with strong noise or the blank one
            assert(max(abs([dx dy] - d.shift(k, :))) < 0.06, sprintf('case %d: true shift', k));
        end
    end

    % sign and axis order: content moves 5.3 columns right and 2.6 rows up
    seed_rng(11);
    fx = load(fullfile(here, '..', 'parity', 'fixed_inputs.mat'), 'texture256');
    tex = double(fx.texture256);        % periodic, so shifts wrap cleanly
    rows = 41:104;                      % 64 rows x 96 columns
    cols = 21:116;
    a = tex(rows, cols);
    s = shifted(tex, 5.3, -2.6);
    b = s(rows, cols);
    [dx, dy, peak] = fftvel.phase_correlation(a, b);
    assert(max(abs([dx dy] - [5.3 -2.6])) < 0.02, sprintf('sign: %.4f %.4f', dx, dy));
    assert(peak > 0.5 && peak <= 1, sprintf('peak %.4f', peak));
    [dx, dy] = fftvel.phase_correlation(b, a);
    assert(max(abs([dx dy] - [-5.3 2.6])) < 0.02, sprintf('swapped: %.4f %.4f', dx, dy));
    [dx, dy] = fftvel.phase_correlation(a.', b.');
    assert(max(abs([dx dy] - [-2.6 5.3])) < 0.02, sprintf('transposed: %.4f %.4f', dx, dy));

    % dx in (-w/2, w/2], dy in (-h/2, h/2]: a periodic shift by 60 of 64 columns is -4
    a = fftvel.make_texture(64, 1.5);
    [dx, dy] = fftvel.phase_correlation(a, circshift(a, [-5 60]));
    assert(max(abs([dx dy] - [-4 -5])) < 0.02, sprintf('wrapped: %.4f %.4f', dx, dy));

    % sub-pixel shifts just inside +-w/2 and +-h/2 keep their sign through the
    % refinement and the wrapping; whole periodic frames, so only the window
    % limits the overlap
    a = tex(1:4:end, 1:2:end);          % 64 rows x 128 columns, periodic
    for d = [63.6 0.3; -63.6 -0.3; 0.3 31.6; -0.3 -31.6].'
        [dx, dy] = fftvel.phase_correlation(a, shifted(a, d(1), d(2)));
        assert(max(abs([dx dy] - d.')) < 0.05, sprintf('near half the size: %.4f %.4f', dx, dy));
    end

    % identical images; uint8 frames give the same result as their double copies
    a = randn(48, 40);
    [dx, dy, peak] = fftvel.phase_correlation(a, a);
    assert(abs(dx) < 1e-12 && abs(dy) < 1e-12 && abs(peak - 1) < 1e-12, 'identical images');
    u = uint8(100 + 20*a);
    v = circshift(u, [0 3]);
    [dx, dy, peak] = fftvel.phase_correlation(u, v);
    [dx2, dy2, peak2] = fftvel.phase_correlation(double(u), double(v));
    assert(isequal([dx dy peak], [dx2 dy2 peak2]), 'uint8 input');

    % blank images: no shift and peak 0; a single row or column: the shift
    % along it, exactly 0 across it
    [dx, dy, peak] = fftvel.phase_correlation(ones(32), 7*ones(32));
    assert(isequal([dx dy peak], [0 0 0]), 'blank images');
    % NaN in an image: no shift and a NaN peak, not a shift to the corner of
    % the fine grid
    a = tex(1:32, 1:32);
    a(6, 8) = NaN;
    [dx, dy, peak] = fftvel.phase_correlation(a, tex(1:32, 1:32));
    assert(dx == 0 && dy == 0 && isnan(peak), 'NaN in an image');
    s = shifted(tex(1, :), -5.7, 0);
    a = tex(1, 41:168);
    b = s(41:168);
    [dx, dy] = fftvel.phase_correlation(a, b);
    assert(abs(dx + 5.7) < 0.03 && dy == 0, sprintf('single row: %.4f %.4f', dx, dy));
    [dx2, dy2] = fftvel.phase_correlation(a.', b.');
    assert(dx2 == 0 && abs(dy2 - dx) < 1e-12, 'single column');

    % sub-pixel accuracy on Fourier-shifted crops of the camera texture: the
    % parabolic fit used before was off by up to 0.17 px between pixels
    tex = double(imread(fullfile(here, '..', '..', 'data', 'surface_texture.png')));
    tex = tex(1:384, 1:384);
    whole = [3 -5; -6 2; 0 0; 4 4; -2 6; 5 -1; -4 -3; 1 -6; 6 3; -1 1];
    err = zeros(10, 2);
    for k = 1:10
        frac = (k - 1)/10;
        sh = whole(k, :) + [frac, 1 - frac];
        s = shifted(tex, sh(1), sh(2));
        [dx, dy] = fftvel.phase_correlation(tex(129:256, 129:256), s(129:256, 129:256));
        err(k, :) = [dx dy] - sh;
    end
    assert(max(abs(err(:))) < 0.02, sprintf('sub-pixel error up to %.4f px', max(abs(err(:)))));
    assert(max(abs(mean(err, 1))) < 0.005, sprintf('mean error %.4f px', max(abs(mean(err, 1)))));
end

function s = shifted(tex, dx, dy)
%SHIFTED  TEX moved by (DX, DY) px (x along columns, y down), periodic band-limited shift.
    [h, w] = size(tex);
    fy = (mod((0:h-1) + floor(h/2), h) - floor(h/2))/h;
    fx = (mod((0:w-1) + floor(w/2), w) - floor(w/2))/w;
    s = real(ifft2(fft2(tex) .* exp(-2i*pi*(fy(:)*dy + fx*dx))));
end
