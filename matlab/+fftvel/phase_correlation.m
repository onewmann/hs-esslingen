function [dx, dy, peak] = phase_correlation(a, b)
%PHASE_CORRELATION  Translation between two images.
%   [DX, DY] = fftvel.phase_correlation(A, B) returns the shift that moves
%   the content of A onto B, in pixels (DX along columns, DY along rows
%   pointing down), each in (-w/2, w/2] and (-h/2, h/2]. Both images are
%   made mean-free and Hann windowed; the lag map is the inverse FFT of
%   their normalised cross-power spectrum, and its largest sample gives the
%   shift to the nearest pixel. It serves as an independent reference for
%   frame sequences.
%
%   Sub-pixel refinement: the lag map is a trigonometric polynomial in the
%   lag, so it can be evaluated between the samples directly from the
%   cross-power spectrum (band-limited interpolation). It is evaluated on a
%   grid of 1/UPSAMPLE px (UPSAMPLE = 20) within 1 px of the largest
%   sample, as two small matrix products with complex exponentials instead
%   of a large zero-padded FFT, and the largest value on that grid is
%   refined by a parabola along each axis. A parabola through three samples
%   of the full-resolution lag map instead, whose peak is only about one
%   pixel wide, is off by up to about 0.17 px for shifts between pixels. On
%   Fourier-shifted 256 px crops of data/surface_texture.png without noise
%   the refined shift is off by 0.002 px (median; at most 0.005 px); the
%   remaining error comes from content moving under the fixed window.
%
%   [DX, DY, PEAK] = ... also returns the height of the interpolated peak,
%   the largest value on the fine grid (1 for two identical images, lower
%   when content moves in or out of the window, changes or is noisy).
%   Blank images give DX = DY = PEAK = 0; NaN or Inf in either image gives
%   DX = DY = 0 and PEAK = NaN.
%
%   Same as fftvel.phase_correlation in Python (python/src/fftvel/phasecorr.py).

    a = double(a);
    b = double(b);
    [h, w] = size(a);
    win = fftvel.internal.hann_window(h)*fftvel.internal.hann_window(w).';
    A = fft2((a - mean(a(:))) .* win);
    B = fft2((b - mean(b(:))) .* win);
    X = B .* conj(A);
    X = X ./ max(abs(X), eps);
    r = real(ifft2(X));
    [~, idx] = max(r(:));
    [i, j] = ind2sub(size(r), idx);
    [oy, ox, peak] = upsampled_peak(X, i - 1, j - 1);
    dy = i - 1 + oy;
    dx = j - 1 + ox;
    if dy > h/2, dy = dy - h; end
    if dx > w/2, dx = dx - w; end
end

function [oy, ox, peak] = upsampled_peak(X, i0, j0)
%UPSAMPLED_PEAK  Maximum of the lag map real(ifft2(X)) next to the lag (I0, J0).
%   Returns the offsets of the maximum from (I0, J0) and the largest value
%   on the grid. Row P and column Q of R hold the lag map at
%   (I0 + UY(P), J0 + UX(Q)); lags count from 0.

    UPSAMPLE = 20;          % same as UPSAMPLE in python/src/fftvel/phasecorr.py
    [h, w] = size(X);
    u = (-UPSAMPLE:UPSAMPLE)/UPSAMPLE;
    uy = u;
    ux = u;
    if h == 1, uy = 0; end  % an axis of length 1 has no sub-pixel lag
    if w == 1, ux = 0; end
    ky = mod((0:h-1) + floor(h/2), h) - floor(h/2);     % signed frequency indices, as in the FFT
    kx = mod((0:w-1) + floor(w/2), w) - floor(w/2);
    ey = exp(2i*pi*((i0 + uy(:))*ky)/h);
    ex = exp(2i*pi*(kx(:)*(j0 + ux))/w);
    R = real(ey*X*ex)/(h*w);
    [~, idx] = max(R(:));
    [p, q] = ind2sub(size(R), idx);
    cy = (numel(uy) + 1)/2;
    cx = (numel(ux) + 1)/2;
    if ~(R(p, q) > R(cy, cx))   % nothing higher than the sample itself: blank images, NaN
        p = cy;
        q = cx;
    end
    peak = R(p, q);
    oy = uy(p) + fftvel.internal.parabolic_offset(R(:, q), p, false)/UPSAMPLE;
    ox = ux(q) + fftvel.internal.parabolic_offset(R(p, :), q, false)/UPSAMPLE;
end
