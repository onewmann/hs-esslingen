function [M, j, C] = mean_projection(S, theta_deg, rmax)
%MEAN_PROJECTION  Mean of a centred spectrum along parallel lines (Radon).
%   [M, J] = MEAN_PROJECTION(S, THETA_DEG, RMAX) projects the square array S
%   (zero frequency at floor(size/2)+1) onto the axis at angle THETA_DEG,
%   measured counter-clockwise from +x with y pointing up. M(:,k) holds the
%   mean of S along the lines at signed distance J = -RMAX..RMAX from the
%   centre, for every angle in THETA_DEG. Only samples inside the disc of
%   radius RMAX are used and every sample is split linearly between its two
%   neighbouring bins.
%
%   Using the mean instead of the sum (as in RADON) removes the dependence
%   on the chord length of the square support, which otherwise biases the
%   ripple profile towards 45 and 135 degrees.
%
%   [M, J, C] = MEAN_PROJECTION(...) also returns the bin weights. They
%   depend only on the geometry; the weights of the last call with more
%   than one angle are cached.

    persistent cache_key cache_C

    npad = size(S, 1);
    c = floor(npad/2) + 1;
    k = (1:npad) - c;
    [U, V] = meshgrid(k, -k);          % U: x to the right, V: y up
    inside = U.^2 + V.^2 <= rmax^2;
    u = U(inside);
    v = V(inside);
    s = S(inside);
    nb = 2*rmax + 3;                   % bins -rmax-1 .. rmax+1
    na = numel(theta_deg);

    key = [npad, rmax, theta_deg(:).'];
    have_C = isequal(key, cache_key);
    if have_C
        C = cache_C;
    else
        C = zeros(2*rmax + 1, na);
    end
    M = zeros(2*rmax + 1, na);
    for a = 1:na
        t = theta_deg(a)*pi/180;
        r = u*cos(t) + v*sin(t);
        j0 = floor(r);
        f = r - j0;
        b = j0 + rmax + 2;
        w = s.*f;
        acc = accumarray(b, s - w, [nb 1]) + accumarray(b + 1, w, [nb 1]);
        if ~have_C
            cnt = accumarray(b, 1 - f, [nb 1]) + accumarray(b + 1, f, [nb 1]);
            C(:, a) = cnt(2:end-1);
        end
        M(:, a) = acc(2:end-1) ./ max(C(:, a), 1e-12);
    end
    if ~have_C && na > 1          % keep the full angle grid, not single-angle calls
        cache_key = key;
        cache_C = C;
    end
    j = (-rmax:rmax).';
end
