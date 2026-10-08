function [dx, dy, peak] = phase_correlation(a, b)
%PHASE_CORRELATION  Translation between two images.
%   [DX, DY] = PHASE_CORRELATION(A, B) returns the shift that moves the
%   content of A onto B, in pixels (DX along columns, DY along rows pointing
%   down), with parabolic sub-pixel refinement. Both images are Hann
%   windowed. It serves as an independent reference for frame sequences.
%   [DX, DY, PEAK] = ... also returns the correlation peak height.

    a = double(a);
    b = double(b);
    [h, w] = size(a);
    win = hann_window(h)*hann_window(w).';
    A = fft2((a - mean(a(:))) .* win);
    B = fft2((b - mean(b(:))) .* win);
    X = B .* conj(A);
    r = real(ifft2(X ./ max(abs(X), eps)));
    [peak, idx] = max(r(:));
    [i, j] = ind2sub(size(r), idx);
    di = parabolic_offset(r(:, j), i, true);
    dj = parabolic_offset(r(i, :), j, true);
    dy = i - 1 + di;
    dx = j - 1 + dj;
    if dy > h/2, dy = dy - h; end
    if dx > w/2, dx = dx - w; end
end
