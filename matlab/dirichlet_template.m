function T = dirichlet_template(jpos, v, L, npad, n, angle_deg)
%DIRICHLET_TEMPLATE  Expected log spectrum along v for 2L+1 exposures.
%   T = DIRICHLET_TEMPLATE(JPOS, V, L, NPAD, N, ANGLE_DEG) returns one
%   column per speed in V, evaluated at the bins JPOS of the projection
%   along ANGLE_DEG. For a random texture, the expected power of the
%   windowed pulsed image along v is
%
%     P(j) = [ (2L+1) + 2*sum_{m=1}^{2L} (2L+1-m) * rho(m*v) * cos(2*pi*j*m*|v|/NPAD) ]
%            / (2L+1)^2
%
%   where rho(m*v) = rho_w(m*vx)*rho_w(m*vy) is the autocorrelation of the
%   separable N-point Hann window at the offset between two exposures m
%   steps apart. Without the window this is the squared Dirichlet kernel;
%   the window blurs it exactly as the FFT of the windowed image does, and
%   the closed form avoids any aliasing of the fine side lobes. T is
%   0.5*log(P + 1e-3).

    jpos = jpos(:);
    v = v(:).';
    c = abs(cosd(angle_deg));
    s = abs(sind(angle_deg));
    K = 2*L + 1;
    P = K*ones(numel(jpos), numel(v));
    for m = 1:2*L
        w = (K - m) * window_autocorr(m*v*c, n) .* window_autocorr(m*v*s, n);
        P = P + 2*cos(2*pi*jpos*(m*v)/npad) .* repmat(w, numel(jpos), 1);
    end
    P = max(P/K^2, 0);
    T = 0.5*log(P + 1e-3);
end

function r = window_autocorr(tau, n)
%WINDOW_AUTOCORR  Normalised autocorrelation of the N-point Hann window.
%   Closed form for a continuous lag; matches the discrete autocorrelation
%   of HANN_WINDOW(N) at integer lags to better than 1e-7 for N >= 64
%   (3e-5 at N = 16).
    x = abs(tau)/(n - 1);
    r = zeros(size(x));
    in = x < 1;
    x = x(in);
    r(in) = (1 - x).*(2 + cos(2*pi*x))/3 + sin(2*pi*x)/(2*pi);
end
