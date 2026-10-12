function r = window_autocorr(tau, n)
%WINDOW_AUTOCORR  Normalised autocorrelation of the N-point Hann window.
%   R = fftvel.internal.window_autocorr(TAU, N) evaluates the closed form
%   for a continuous lag TAU (any size); it matches the discrete
%   autocorrelation of fftvel.internal.hann_window(N) at integer lags to
%   better than 1e-7 for N >= 64 (3e-5 at N = 16). Used by
%   fftvel.internal.dirichlet_template and fftvel.internal.box_template.
    x = abs(tau)/(n - 1);
    r = zeros(size(x));
    in = x < 1;
    x = x(in);
    r(in) = (1 - x).*(2 + cos(2*pi*x))/3 + sin(2*pi*x)/(2*pi);
end
