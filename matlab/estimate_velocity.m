function est = estimate_velocity(bp, varargin)
%ESTIMATE_VELOCITY  Velocity of a moving texture from one pulsed-exposure image.
%   EST = ESTIMATE_VELOCITY(BP) estimates direction and speed of the motion
%   recorded in the pulsed-exposure image BP. BP is the sum or mean of 2L+1
%   exposures taken at equal time steps while the texture moves by v per
%   step (a strobe within one camera frame, or the average of 2L+1
%   consecutive video frames).
%
%   EST = ESTIMATE_VELOCITY(BP, 'Name', Value, ...) sets options:
%     'L'            half number of exposures, 2L+1 in total (default 4)
%     'PadFactor'    zero padding factor of the FFT (default 2)
%     'MinQuality'   lowest fit quality accepted for EST.valid (default 0.6);
%                    small images need more, see EST.threshold
%     'Diagnostics'  also return intermediate results in EST.diag
%                    (default false)
%
%   Fields of EST:
%     angle_deg  direction of motion in degrees in [0, 180), counter-
%                clockwise from the +x axis as seen on screen (y up).
%                A single image cannot tell v from -v.
%     speed      |v| in pixels per exposure step
%     vx, vy     speed*[cosd(angle_deg), -sind(angle_deg)], i.e. the
%                velocity in image coordinates (x along columns, y along
%                rows pointing down); -[vx vy] is equally possible
%     geba       max(WR)/mean(WR) of the ripple profile
%     quality    correlation between the measured profile and the model
%     in_range   false if the best speed lies below the measurable range
%     threshold  quality needed for this image size: the larger of MinQuality
%                and tanh(9/sqrt(m-5)), m = number of profile points
%     valid      quality >= threshold and in_range
%
%   Method: the motion multiplies the image spectrum by a Dirichlet kernel,
%   which shows up as straight lines perpendicular to v, spaced
%   Npad/|v| bins apart. The direction is the angle whose mean Radon
%   projection of the log spectrum has the largest total variation (ripple
%   profile WR). The speed is found by fitting the Dirichlet model for
%   2L+1 exposures to the projection along that direction.
%
%   The model, the ripple profile WR and GEBA follow Chapter 3 of the
%   doctoral thesis cited in README.md; the speed fit and the validity
%   threshold were added in this implementation.
%
%   Non-square images are cropped to the central square. They are never
%   resized, so the speed stays in camera pixels.
%
%   See also SIMULATE_PULSED_IMAGE, PLOT_ESTIMATE.

    opts = parse_options(varargin, struct('L', 4, 'PadFactor', 2, ...
        'MinQuality', 0.6, 'Diagnostics', false));
    L = double(opts.L);
    if ~isscalar(L) || ~isreal(L) || L < 1 || L ~= round(L)
        error('estimate_velocity:L', 'L must be a positive integer.');
    end

    [S, n, npad] = log_spectrum(bp, opts.PadFactor);
    rmax = floor(npad/2) - 1;
    rfrac = 0.9;

    % Direction: maximum of the ripple profile WR(theta)
    theta = 0:179;
    [M, j] = mean_projection(S, theta, rmax);
    wr = ripple_profile(M, j, rfrac);
    [wrmax, i] = max(wr);
    angle = mod(theta(i) + parabolic_offset(wr, i, true), 180);
    if mean(wr) > 0
        geba = wrmax / mean(wr);
    else
        geba = 1;        % no ripple in any direction
    end

    % Speed: fit the Dirichlet model to the projection along that direction
    m = mean_projection(S, angle, rmax);
    pos = j >= 2 & j <= rfrac*rmax;
    mf = flipud(m);
    g = 0.5*(m(pos) + mf(pos));
    jpos = j(pos);
    [speed, quality, fit] = fit_speed(g, jpos, L, npad, n, angle);

    est = struct();
    est.angle_deg = angle;
    est.speed = speed;
    est.vx = speed*cosd(angle);
    est.vy = -speed*sind(angle);
    est.geba = geba;
    est.quality = quality;
    est.in_range = fit.in_range;
    est.threshold = max(opts.MinQuality, null_quality(numel(jpos)));
    est.valid = quality >= est.threshold && fit.in_range;
    est.L = L;
    est.n = n;
    est.npad = npad;
    if opts.Diagnostics
        d = struct();
        d.spectrum = S;
        d.theta = theta;
        d.wr = wr;
        d.profile_j = jpos;
        d.profile = g;
        d.model = fit.model;
        d.speed_grid = fit.v;
        d.speed_score = fit.score;
        est.diag = d;
    end
end

function q = null_quality(m)
%NULL_QUALITY  Fit quality that images without motion stay below.
%   Under the null hypothesis the Fisher transform atanh(q) of the best
%   correlation scales with 1/sqrt(m-5) for a profile of m points (four
%   are taken by the cubic, one by the fit). On 1008 images without motion
%   (white noise, random and camera textures, 32 to 256 px, L = 1 to 10)
%   the largest value was tanh(8.6/sqrt(m-5)); 9 adds a margin.
    if m <= 5
        q = Inf;
    else
        q = tanh(9/sqrt(m - 5));
    end
end
