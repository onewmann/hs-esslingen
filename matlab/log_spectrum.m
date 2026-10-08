function [S, n, npad] = log_spectrum(bp, pad_factor)
%LOG_SPECTRUM  Centred log-magnitude spectrum of a windowed image.
%   [S, N, NPAD] = LOG_SPECTRUM(BP, PAD_FACTOR) crops BP to its central
%   N-by-N square, removes the mean, applies a 2-D Hann window and returns
%   log(|B| + 0.05*median|B|) of the NPAD-by-NPAD FFT, NPAD = PAD_FACTOR*N,
%   with the zero frequency at index floor(NPAD/2)+1 (fftshift layout).
%   The floor relative to the median keeps the result independent of the
%   image scale.

    if nargin < 2
        pad_factor = 2;
    end
    if pad_factor < 1 || pad_factor ~= round(pad_factor)
        error('log_spectrum:pad', 'PadFactor must be a positive integer.');
    end
    x = crop_square(double(bp));
    n = size(x, 1);
    if n < 16
        error('log_spectrum:size', 'Image must be at least 16x16 pixels.');
    end
    npad = pad_factor*n;
    x = x - mean(x(:));
    w = hann_window(n);
    B = fftshift(fft2(x .* (w*w.'), npad, npad));
    A = abs(B);
    fl = 0.05*median(A(:));
    if fl <= 0
        fl = 0.05*mean(A(:));
    end
    if fl <= 0
        fl = 1;          % blank image: flat spectrum, the estimate comes out not valid
    end
    S = log(A + fl);
end
