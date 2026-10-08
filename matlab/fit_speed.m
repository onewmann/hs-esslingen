function [speed, quality, fit] = fit_speed(g, jpos, L, npad, n, angle_deg)
%FIT_SPEED  Speed from the projection along the direction of motion.
%   [SPEED, QUALITY] = FIT_SPEED(G, JPOS, L, NPAD, N, ANGLE_DEG) compares the
%   measured half profile G(JPOS) along ANGLE_DEG with the Dirichlet model
%   of 2L+1 exposures (DIRICHLET_TEMPLATE) for speeds on a logarithmic grid.
%   A cubic polynomial in JPOS absorbs the smooth texture spectrum. SPEED
%   maximises the correlation between the detrended profile and the
%   detrended model; QUALITY is that correlation.
%
%   The measurable range runs from VMIN = 1.5*NPAD/((2L+1)*max(JPOS)), where
%   the first spectral zero still lies well inside the profile, to NPAD/4,
%   where the motion lines are four bins apart. The grid extends down to
%   VMIN/4 as a guard band: slower motion then lands there instead of on a
%   wrong in-range speed, and FIT.in_range is false.
%
%   FIT holds the grid (FIT.v), the scores (FIT.score), VMIN, the flag
%   FIT.in_range and the fitted model curve (FIT.model) for plotting.

    step = 0.005;
    jpos = jpos(:);
    g = g(:);
    vmin = 1.5*npad/((2*L + 1)*max(jpos));
    vmax = npad/4;
    nv = floor((log(vmax) - log(vmin))/step) + 1;
    kext = floor(log(4)/step);
    v = exp(log(vmin) + (-kext:nv-1)*step);

    T = dirichlet_template(jpos, v, L, npad, n, angle_deg);
    x = jpos/max(jpos);
    X = [x.^3, x.^2, x, ones(size(x))];
    [Q, R] = qr(X, 0); %#ok<ASGLU>
    gp = g - Q*(Q.'*g);
    Tp = T - Q*(Q.'*T);
    tn = sqrt(sum(Tp.^2, 1));
    score = (gp.'*Tp) ./ (norm(gp)*tn + realmin);

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
