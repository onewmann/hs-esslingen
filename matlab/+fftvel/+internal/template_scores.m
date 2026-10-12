function [score, gp, Tp] = template_scores(g, jpos, T)
%TEMPLATE_SCORES  Correlation of a detrended profile with detrended templates.
%   [SCORE, GP, TP] = fftvel.internal.template_scores(G, JPOS, T) removes a
%   cubic polynomial in JPOS from the profile G (a column) and from each
%   column of T, and returns the correlation SCORE of the two for each
%   column, with the detrended profile GP and templates TP. Columns whose
%   detrended norm is below FLAT_MODEL (1e-6) of their norm (a flat model)
%   score 0. Shared by fftvel.internal.fit_speed and
%   fftvel.internal.fit_blur.

    c = fftvel.internal.constants();
    x = jpos/max(jpos);
    X = zeros(numel(x), c.TREND_DEGREE + 1);     % cubic: [x.^3, x.^2, x, 1]
    for p = 0:c.TREND_DEGREE
        X(:, c.TREND_DEGREE + 1 - p) = x.^p;
    end
    [Q, R] = qr(X, 0); %#ok<ASGLU>
    gp = g - Q*(Q.'*g);
    Tp = T - Q*(Q.'*T);
    tn = sqrt(sum(Tp.^2, 1));
    score = (gp.'*Tp) ./ (norm(gp)*tn + realmin);
    score(tn <= c.FLAT_MODEL*sqrt(sum(T.^2, 1))) = 0;      % flat model: no information
end
