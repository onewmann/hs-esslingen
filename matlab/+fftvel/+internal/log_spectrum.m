function [S, n, npad] = log_spectrum(bp, pad_factor, roi, clip_outliers)
%LOG_SPECTRUM  Centred log-magnitude spectrum of a windowed image.
%   [S, N, NPAD] = fftvel.internal.log_spectrum(BP, PAD_FACTOR) crops BP
%   to its central N-by-N square, removes the mean, applies a 2-D Hann
%   window and returns log(|B| + 0.05*median|B|) of the NPAD-by-NPAD FFT,
%   NPAD = PAD_FACTOR*N, with the zero frequency at index floor(NPAD/2)+1
%   (fftshift layout). The floor relative to the median keeps the result
%   independent of the image scale.
%
%   [S, N, NPAD] = fftvel.internal.log_spectrum(..., ROI) crops to the central
%   ROI-by-ROI square instead (ROI empty, or larger than the image: the
%   largest central square).
%
%   [S, N, NPAD] = fftvel.internal.log_spectrum(..., ROI, CLIP_OUTLIERS) with
%   CLIP_OUTLIERS true first sets the pixels more than 6 robust standard
%   deviations (1.4826*MAD) from the median, and all pixels within 3 px of
%   them (fftvel.internal.dilate), to the median.

    c = fftvel.internal.constants();
    if nargin < 2
        pad_factor = c.DEFAULT_PAD_FACTOR;
    end
    if nargin < 3
        roi = [];
    end
    if nargin < 4
        clip_outliers = false;
    end
    pad_factor = fftvel.internal.check_option(pad_factor, 'PadFactor', 'positive', 'log_spectrum:pad');
    side = [];
    if ~isempty(roi)
        roi = fftvel.internal.check_option(roi, 'ROI', 'roi', 'log_spectrum:roi');
        if ndims(bp) == 2
            side = min([roi, size(bp)]);
        end
    end
    clip_outliers = fftvel.internal.check_option(clip_outliers, 'ClipOutliers', 'flag', ...
        'log_spectrum:clip');
    x = double(fftvel.internal.crop_square(bp, side));
    n = size(x, 1);
    if n < c.MIN_SIZE
        error('log_spectrum:size', 'Image must be at least %dx%d pixels.', c.MIN_SIZE, c.MIN_SIZE);
    end
    if any(~isfinite(x(:)))
        error('log_spectrum:finite', 'Image contains NaN or Inf values.');
    end
    if clip_outliers
        med = median(x(:));
        s = c.MAD_SCALE*median(abs(x(:) - med));
        if s > 0
            x(fftvel.internal.dilate(abs(x - med) > c.CLIP_SIGMA*s, c.CLIP_RADIUS)) = med;
        end
    end
    npad = pad_factor*n;
    x = x - mean(x(:));
    w = fftvel.internal.hann_window(n);
    B = fftshift(fft2(x .* (w*w.'), npad, npad));
    A = abs(B);
    fl = c.LOG_FLOOR*median(A(:));
    if fl <= 0
        fl = c.LOG_FLOOR*mean(A(:));
    end
    if fl <= 0
        fl = c.BLANK_FLOOR;      % blank image: flat spectrum, the estimate comes out not valid
    end
    S = log(A + fl);
end
