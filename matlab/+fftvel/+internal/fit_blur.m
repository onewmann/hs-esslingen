function [blur_length, quality, fit] = fit_blur(g, jpos, npad, n, angle_deg)
%FIT_BLUR  Length of a continuous motion blur that explains the profile.
%   [BLUR_LENGTH, QUALITY] = fftvel.internal.fit_blur(G, JPOS, NPAD, N,
%   ANGLE_DEG) is the fit of fftvel.internal.fit_speed with
%   fftvel.internal.box_template (one long exposure, no flashes) in place
%   of the Dirichlet model. A blur of D px has the main lobe of K flashes
%   D/K apart; only the periodic side lines of the flashes tell them
%   apart. The grid of blur lengths is K times the speed grid of
%   fftvel.internal.fit_speed (it does not depend on K), from
%   1.5*NPAD/max(JPOS)/4 up to NPAD/4, where the zeros of the sinc are four
%   bins apart. FIT holds the grid (FIT.d) and the scores (FIT.score).

    c = fftvel.internal.constants();
    step = c.SPEED_STEP;
    jpos = jpos(:);
    g = g(:);
    dref = c.VMIN_FACTOR*npad/max(jpos);          % K*vmin of fit_speed
    dmax = npad/c.MIN_LINE_SPACING;
    nd = floor((log(dmax) - log(dref))/step) + 1;
    kext = floor(log(c.GUARD_FACTOR)/step);
    d = exp(log(dref) + (-kext:nd-1)*step);
    score = fftvel.internal.template_scores(g, jpos, ...
        fftvel.internal.box_template(jpos, d, npad, n, angle_deg));
    [quality, k] = max(score);
    dk = fftvel.internal.parabolic_offset(score, k, false);
    blur_length = exp(log(d(k)) + dk*step);
    fit = struct();
    fit.d = d;
    fit.score = score;
end
