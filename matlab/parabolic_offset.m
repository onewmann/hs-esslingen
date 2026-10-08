function d = parabolic_offset(y, i, periodic)
%PARABOLIC_OFFSET  Sub-sample position of a maximum by parabolic fit.
%   D = PARABOLIC_OFFSET(Y, I, PERIODIC) fits a parabola through Y(I-1),
%   Y(I), Y(I+1) and returns the offset of its vertex from I, in samples
%   (between -0.5 and 0.5). With PERIODIC true the ends of Y wrap around;
%   otherwise D is 0 at the ends. D is 0 if the three points are not
%   concave.

    n = numel(y);
    if periodic
        ym = y(mod(i - 2, n) + 1);
        yp = y(mod(i, n) + 1);
    else
        if i <= 1 || i >= n
            d = 0;
            return
        end
        ym = y(i - 1);
        yp = y(i + 1);
    end
    den = ym - 2*y(i) + yp;
    if den < 0
        d = 0.5*(ym - yp)/den;
    else
        d = 0;
    end
end
