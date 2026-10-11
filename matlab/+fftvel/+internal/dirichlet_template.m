function T = dirichlet_template(jpos, v, L, npad, n, angle_deg, K, duty)
%DIRICHLET_TEMPLATE  Expected log spectrum along v for K exposures.
%   T = fftvel.internal.dirichlet_template(JPOS, V, L, NPAD, N, ANGLE_DEG)
%   returns one column per speed in V, evaluated at the bins JPOS of the
%   projection along ANGLE_DEG, for K = 2L+1 exposures. For a random
%   texture, the expected power of the windowed pulsed image along v is
%
%     P(j) = [ K + 2*sum_{m=1}^{K-1} (K-m) * rho(m*v) * cos(2*pi*j*m*|v|/NPAD) ] / K^2
%
%   where rho(m*v) = rho_w(m*vx)*rho_w(m*vy) is the autocorrelation of the
%   separable N-point Hann window at the offset between two exposures m
%   steps apart. Without the window this is the squared Dirichlet kernel;
%   the window blurs it exactly as the FFT of the windowed image does, and
%   the closed form avoids any aliasing of the fine side lobes. T is
%   0.5*log(P + 1e-3).
%
%   T = fftvel.internal.dirichlet_template(..., K) uses K exposures
%   instead of 2L+1 (K empty: 2L+1).
%   T = fftvel.internal.dirichlet_template(..., K, DUTY) with DUTY > 0
%   (exposure time / time step) adds the motion blur within each exposure:
%   P is multiplied by sinc(j*|v|*DUTY/NPAD)^2.

    if nargin < 7 || isempty(K)
        K = 2*L + 1;
    end
    if nargin < 8
        duty = 0;
    end
    jpos = jpos(:);
    v = v(:).';
    c = abs(cosd(angle_deg));
    s = abs(sind(angle_deg));
    P = K*ones(numel(jpos), numel(v));
    for m = 1:K-1
        w = (K - m) * window_autocorr(m*v*c, n) .* window_autocorr(m*v*s, n);
        P = P + 2*cos(2*pi*jpos*(m*v)/npad) .* repmat(w, numel(jpos), 1);
    end
    P = max(P/K^2, 0);
    if duty > 0
        x = jpos*v*duty/npad;
        sc = ones(size(x));          % sinc without the Signal Processing Toolbox
        nz = x ~= 0;
        sc(nz) = sin(pi*x(nz))./(pi*x(nz));
        P = P .* sc.^2;
    end
    T = 0.5*log(P + 1e-3);
end

function r = window_autocorr(tau, n)
%WINDOW_AUTOCORR  Normalised autocorrelation of the N-point Hann window.
%   Closed form for a continuous lag; matches the discrete autocorrelation
%   of fftvel.internal.hann_window(N) at integer lags to better than 1e-7
%   for N >= 64 (3e-5 at N = 16).
    x = abs(tau)/(n - 1);
    r = zeros(size(x));
    in = x < 1;
    x = x(in);
    r(in) = (1 - x).*(2 + cos(2*pi*x))/3 + sin(2*pi*x)/(2*pi);
end
