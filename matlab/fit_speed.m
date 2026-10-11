function [speed, quality, fit] = fit_speed(g, jpos, L, npad, n, angle_deg, K, duty)
%FIT_SPEED  Speed from the projection along the direction of motion.
%   [SPEED, QUALITY] = FIT_SPEED(G, JPOS, L, NPAD, N, ANGLE_DEG) compares the
%   measured half profile G(JPOS) along ANGLE_DEG with the Dirichlet model
%   of 2L+1 exposures (DIRICHLET_TEMPLATE) for speeds on a logarithmic grid.
%   A cubic polynomial in JPOS absorbs the smooth texture spectrum. SPEED
%   maximises the correlation between the detrended profile and the
%   detrended model; QUALITY is that correlation.
%
%   FIT_SPEED(..., K, DUTY) uses K exposures instead of 2L+1 (K empty:
%   2L+1) and the motion blur DUTY within each exposure (default 0), see
%   DIRICHLET_TEMPLATE.
%
%   The measurable range runs from VMIN = 1.5*NPAD/(K*max(JPOS)), where
%   the first spectral zero still lies well inside the profile, to NPAD/4,
%   where the motion lines are four bins apart. The grid extends down to
%   VMIN/4 as a guard band: slower motion then lands there instead of on a
%   wrong in-range speed, and FIT.in_range is false. With PadFactor > 2
%   the grid reaches speeds at which exposures one step apart no longer
%   overlap within the window; the model is then flat (DUTY = 0), and its
%   score is set to 0 instead of a ratio of rounding errors.
%
%   FIT holds the grid (FIT.v), the scores (FIT.score), VMIN, the flag
%   FIT.in_range and the fitted model curve (FIT.model) for plotting.

    if nargin < 7 || isempty(K)
        K = 2*L + 1;
    end
    if nargin < 8
        duty = 0;
    end
    step = 0.005;
    jpos = jpos(:);
    g = g(:);
    vmin = 1.5*npad/(K*max(jpos));
    vmax = npad/4;
    nv = floor((log(vmax) - log(vmin))/step) + 1;
    kext = floor(log(4)/step);
    v = exp(log(vmin) + (-kext:nv-1)*step);

    T = dirichlet_template(jpos, v, L, npad, n, angle_deg, K, duty);
    x = jpos/max(jpos);
    X = [x.^3, x.^2, x, ones(size(x))];
    [Q, R] = qr(X, 0); %#ok<ASGLU>
    gp = g - Q*(Q.'*g);
    Tp = T - Q*(Q.'*T);
    tn = sqrt(sum(Tp.^2, 1));
    score = (gp.'*Tp) ./ (norm(gp)*tn + realmin);
    score(tn <= 1e-6*sqrt(sum(T.^2, 1))) = 0;       % flat model: no information

    [quality, k] = max(score);
    dk = parabolic_offset(score, k, false);
    speed = exp(log(v(k)) + dk*step);

    fit = struct();
    fit.v = v;
    fit.score = score;
    fit.vmin = vmin;
    fit.in_range = k > kext && k < numel(v);
    tk = Tp(:, k);
    fit.model = (g - gp) + (gp.'*tk)/(tk.'*tk)*tk;
end
