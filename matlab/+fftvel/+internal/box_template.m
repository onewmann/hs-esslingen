function T = box_template(jpos, d, npad, n, angle_deg)
%BOX_TEMPLATE  Expected log spectrum along the motion for one continuous exposure.
%   T = fftvel.internal.box_template(JPOS, D, NPAD, N, ANGLE_DEG) returns
%   one column per blur length in D, evaluated at the bins JPOS of the
%   projection along ANGLE_DEG. A single exposure during which the texture
%   moves by D px (motion blur, no flashes) multiplies the spectrum by
%   sinc(j*D/NPAD)^2. As in fftvel.internal.dirichlet_template the window
%   is included exactly: the blur is sampled every H = BLUR_STEP = 0.5 px,
%
%     P(j) = (H/D) * [ 1 + 2*sum_{q>=1} max(1 - q*H/D, 0) * rho(q*H) * cos(2*pi*j*q*H/NPAD) ]
%
%   with rho the autocorrelation of the window
%   (fftvel.internal.window_autocorr). For D = K*H this is the template of
%   K exposures H apart. T is 0.5*log(P + 1e-3).

    cst = fftvel.internal.constants();
    h = cst.BLUR_STEP;
    jpos = jpos(:);
    d = d(:).';
    c = abs(cosd(angle_deg));
    s = abs(sind(angle_deg));
    u = h*(1:ceil(max(d)/h));                % offsets 0 < u <= max(d)
    rho = fftvel.internal.window_autocorr(u*c, n) .* fftvel.internal.window_autocorr(u*s, n);
    C = cos(2*pi*jpos*u/npad) .* repmat(rho, numel(jpos), 1);
    W = max(1 - repmat(u.', 1, numel(d)) ./ repmat(d, numel(u), 1), 0);
    P = max(repmat(h./d, numel(jpos), 1) .* (1 + 2*(C*W)), 0);
    T = 0.5*log(P + cst.TEMPLATE_FLOOR);
end
