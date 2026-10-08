function w = hann_window(n)
%HANN_WINDOW  Symmetric Hann window as a column vector.
%   W = HANN_WINDOW(N) equals hann(N) from the Signal Processing Toolbox
%   and numpy.hanning(N), without needing either.

    if n == 1
        w = 1;
    else
        w = 0.5 - 0.5*cos(2*pi*(0:n-1).'/(n - 1));
    end
end
