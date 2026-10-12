function est = estimate_velocity(bp, varargin)
%ESTIMATE_VELOCITY  Velocity of a moving texture from one pulsed-exposure image.
%   EST = fftvel.estimate_velocity(BP) estimates direction and speed of the
%   motion recorded in the pulsed-exposure image BP. BP is the sum or mean
%   of 2L+1 exposures taken at equal time steps while the texture moves by
%   v per step (a strobe within one camera frame, or the average of 2L+1
%   consecutive video frames, see fftvel.pulsed_from_stack).
%
%   EST = fftvel.estimate_velocity(BP, 'Name', Value, ...) sets options:
%     'L'             half number of exposures, 2L+1 in total (default 4)
%     'Exposures'     exact number of exposures (flashes) K >= 2, also
%                     even; overrides L for the model (default [] = 2L+1).
%                     For slow motion a wrong count scales the speed by
%                     K_true/K without lowering the quality.
%     'Duty'          duration of one exposure (flash pulse, or camera
%                     exposure time) divided by the time step, in [0, 1)
%                     (default 0 = short flashes). The motion blur within
%                     each exposure, Duty*|v| pixels, must either be given
%                     here or stay well below about 3 px.
%     'ROI'           side in pixels of the central square that is analysed,
%                     at least 16 (default [] = the largest central square;
%                     larger values also give the largest square). The run
%                     time grows with ROI^2.
%     'PadFactor'     zero padding factor of the FFT (default 2)
%     'MinQuality'    lowest fit quality accepted for EST.valid, in [0, 1]
%                     (default 0.6); small images need more, see
%                     EST.threshold
%     'ClipOutliers'  set pixels more than 6 robust standard deviations
%                     (1.4826*MAD) from the median, and all pixels within
%                     3 px of them, to the median before the FFT, against
%                     static bright spots such as glints or hot pixels
%                     (default false)
%     'RejectBlur'    also require EST.blur_like to be false for EST.valid
%                     (default false), against single exposures of a moving
%                     surface (motion blur, no flashes). Most multi-flash
%                     images of slow motion (below about 2.2 px per step)
%                     with noise are blur-like too, so this rejects most
%                     valid results there.
%     'Diagnostics'   also return intermediate results in EST.diag
%                     (default false)
%   All option values are checked; text such as '4' is an error.
%
%   Fields of EST:
%     angle_deg     direction of motion in degrees in [0, 180), counter-
%                   clockwise from the +x axis as seen on screen (y up).
%                   A single image cannot tell v from -v.
%     speed         |v| in pixels per exposure step
%     vx, vy        speed*[cosd(angle_deg), -sind(angle_deg)], i.e. the
%                   velocity in image coordinates (x along columns, y along
%                   rows pointing down); -[vx vy] is equally possible
%     geba          max(WR)/mean(WR) of the ripple profile
%     quality       correlation between the measured profile and the model
%     quality_perp  quality of the same fit perpendicular to angle_deg
%     in_range      false if the best speed lies in the guard band below
%                   the measurable range or on the last point of the speed
%                   grid, its top end at about npad/4
%     isotropic     true if the perpendicular fit is about as good
%                   (quality - quality_perp < 0.15) and gives about the
%                   same speed (|log(speed_perp/speed)| < 0.12). Motion
%                   makes the spectrum anisotropic; defocus blur and square
%                   repeating textures do not.
%     threshold     quality needed for this image size: the larger of
%                   MinQuality and tanh(9/sqrt(m-5)), with m the number of
%                   independent profile points, min(points,
%                   floor(2*n*points/npad)); padding beyond PadFactor 2
%                   does not lower it
%     valid         quality >= threshold and in_range and not isotropic
%                   (and not blur_like with RejectBlur)
%     L             L as passed
%     exposures     number of exposures K used in the model
%     n             side of the analysed central square in pixels
%     npad          FFT length
%     quality_blur  quality of the same fit along angle_deg with a single
%                   continuous exposure (motion blur of a free length, no
%                   flashes) in place of the K flashes
%     blur_length   blur length in pixels of that fit; a blur of d px
%                   looks like K flashes about d/K apart
%     blur_like     true if the blur model fits better (quality_blur >
%                   quality) and the speed is below npad/max(profile_j),
%                   about 2.2 px per step, so that the repeated lines of
%                   the flashes (npad/speed bins apart), which alone tell
%                   flashes from a blur, lie beyond the fitted profile.
%                   Most single exposures of a moving surface are
%                   blur-like, but below that speed so are most
%                   multi-flash images with noise: there the flag says
%                   that the image cannot tell flashes from a blur, not
%                   that it is a blur. Only RejectBlur makes it count for
%                   valid
%     saturated     share of the analysed pixels at the maximum of the
%                   square, before any processing: for an integer image
%                   whose maximum is INTMAX of its class, the share at that
%                   value; otherwise the share at the maximum if it is at
%                   least 0.5 % and the square is not constant, else 0.
%                   Saturated areas and glints can give wrong results that
%                   are marked valid
%     diag          only with 'Diagnostics', true: a struct of intermediate
%                   results, spectrum (log spectrum, npad x npad), theta
%                   and wr (angle grid in degrees and ripple profile WR),
%                   profile_j, profile and model (bins of the half profile
%                   along angle_deg, the profile and the fitted model),
%                   speed_grid and speed_score (speed grid and fit quality
%                   per speed), speed_perp (speed of the perpendicular
%                   fit), blur_grid and blur_score (blur lengths and
%                   quality of the blur fit)
%
%   Method: the motion multiplies the image spectrum by a Dirichlet kernel,
%   which shows up as straight lines perpendicular to v, spaced
%   Npad/|v| bins apart. The direction is the angle whose mean Radon
%   projection of the log spectrum has the largest total variation (ripple
%   profile WR). The speed is found by fitting the Dirichlet model for
%   2L+1 exposures to the projection along that direction.
%
%   The model, the ripple profile WR and GEBA follow Chapter 3 of J. Horn's
%   doctoral thesis (Universitaet Karlsruhe (TH), 2006; full reference in
%   README.md). The speed fit, the validity threshold, the isotropy check
%   and the blur check were added in this implementation.
%
%   Non-square images are cropped to the central square. They are never
%   resized, so the speed stays in camera pixels.
%
%   fftvel.to_physical converts EST to length units per second. The
%   defaults and all other numeric constants of the method are listed in
%   fftvel.internal.constants (and in tests/constants.json).
%
%   See also fftvel.to_physical, fftvel.pulsed_from_stack,
%   fftvel.simulate_pulsed_image, fftvel.plot_estimate,
%   fftvel.angle_difference.

    c = fftvel.internal.constants();
    opts = fftvel.internal.parse_options(varargin, struct('L', c.DEFAULT_L, 'Exposures', [], ...
        'Duty', 0, 'ROI', [], 'PadFactor', c.DEFAULT_PAD_FACTOR, 'MinQuality', ...
        c.DEFAULT_MIN_QUALITY, 'ClipOutliers', false, 'RejectBlur', false, 'Diagnostics', false));
    L = fftvel.internal.check_option(opts.L, 'L', 'positive', 'estimate_velocity:L');
    if isempty(opts.Exposures)
        K = 2*L + 1;
    else
        K = fftvel.internal.check_option(opts.Exposures, 'Exposures', 'exposures', ...
            'estimate_velocity:Exposures');
    end
    duty = fftvel.internal.check_option(opts.Duty, 'Duty', 'duty', 'estimate_velocity:Duty');
    roi = opts.ROI;
    if ~isempty(roi)
        roi = fftvel.internal.check_option(roi, 'ROI', 'roi', 'estimate_velocity:ROI');
    end
    pad = fftvel.internal.check_option(opts.PadFactor, 'PadFactor', 'positive', ...
        'estimate_velocity:PadFactor');
    min_quality = fftvel.internal.check_option(opts.MinQuality, 'MinQuality', 'unit', ...
        'estimate_velocity:MinQuality');
    clip = fftvel.internal.check_option(opts.ClipOutliers, 'ClipOutliers', 'flag', ...
        'estimate_velocity:ClipOutliers');
    reject_blur = fftvel.internal.check_option(opts.RejectBlur, 'RejectBlur', 'flag', ...
        'estimate_velocity:RejectBlur');
    diagnostics = fftvel.internal.check_option(opts.Diagnostics, 'Diagnostics', 'flag', ...
        'estimate_velocity:Diagnostics');

    [S, n, npad] = fftvel.internal.log_spectrum(bp, pad, roi, clip);
    saturated = fftvel.internal.saturated_share(fftvel.internal.crop_square(bp, n));  % as given
    rmax = floor(npad/2) - 1;
    rfrac = c.RFRAC;

    % Direction: maximum of the ripple profile WR(theta)
    theta = 0:c.ANGLE_STEP:180 - c.ANGLE_STEP;
    [M, j] = fftvel.internal.mean_projection(S, theta, rmax);
    wr = fftvel.internal.ripple_profile(M, j, rfrac);
    [wrmax, i] = max(wr);
    angle = mod(theta(i) + c.ANGLE_STEP*fftvel.internal.parabolic_offset(wr, i, true), 180);
    if mean(wr) > 0
        geba = wrmax / mean(wr);
    else
        geba = 1;        % no ripple in any direction
    end

    % Speed: fit the Dirichlet model to the projection along that direction
    m = fftvel.internal.mean_projection(S, angle, rmax);
    pos = j >= c.PROFILE_START & j <= rfrac*rmax;
    mf = flipud(m);
    g = 0.5*(m(pos) + mf(pos));
    jpos = j(pos);
    [speed, quality, fit] = fftvel.internal.fit_speed(g, jpos, L, npad, n, angle, K, duty);

    % Isotropy check: the same fit perpendicular to the direction found
    ap = mod(angle + 90, 180);
    mp = fftvel.internal.mean_projection(S, ap, rmax);
    mpf = flipud(mp);
    [speed_perp, quality_perp] = fftvel.internal.fit_speed(0.5*(mp(pos) + mpf(pos)), jpos, L, ...
        npad, n, ap, K, duty);
    isotropic = quality - quality_perp < c.ISO_GAP && abs(log(speed_perp/speed)) < c.ISO_LOG_SPEED;

    % Blur check: one continuous exposure (motion blur, no flashes) along the same direction. Only the
    % repeated lines of the flashes, npad/speed bins apart, tell the two apart; below npad/max(jpos)
    % they lie beyond the profile
    [blur_length, quality_blur, fitb] = fftvel.internal.fit_blur(g, jpos, npad, n, angle);
    blur_like = quality - quality_blur < c.BLUR_GAP && speed*max(jpos) < npad;

    % Independent profile points: zero padding beyond 2 only interpolates
    npts = numel(jpos);
    m_eff = min(npts, floor(c.INDEPENDENT_PAD*n*npts/npad));

    est = struct();
    est.angle_deg = angle;
    est.speed = speed;
    est.vx = speed*cosd(angle);
    est.vy = -speed*sind(angle);
    est.geba = geba;
    est.quality = quality;
    est.quality_perp = quality_perp;
    est.in_range = fit.in_range;
    est.isotropic = isotropic;
    est.threshold = max(min_quality, null_quality(m_eff, c));
    est.valid = quality >= est.threshold && fit.in_range && ~isotropic && ~(reject_blur && blur_like);
    est.L = L;
    est.exposures = K;
    est.n = n;
    est.npad = npad;
    est.quality_blur = quality_blur;
    est.blur_length = blur_length;
    est.blur_like = blur_like;
    est.saturated = saturated;
    if diagnostics
        d = struct();
        d.spectrum = S;
        d.theta = theta;
        d.wr = wr;
        d.profile_j = jpos;
        d.profile = g;
        d.model = fit.model;
        d.speed_grid = fit.v;
        d.speed_score = fit.score;
        d.speed_perp = speed_perp;
        d.blur_grid = fitb.d;
        d.blur_score = fitb.score;
        est.diag = d;
    end
end

function q = null_quality(m, c)
%NULL_QUALITY  Fit quality that images without motion stay below.
%   Q = NULL_QUALITY(M, C) with C from fftvel.internal.constants is
%   tanh(C.NULL_K/sqrt(M - C.NULL_OFFSET)), that is tanh(9/sqrt(M-5)).
%   Under the null hypothesis the Fisher transform atanh(q) of the best
%   correlation scales with 1/sqrt(m-5) for a profile of m independent
%   points (four are taken by the cubic, one by the fit). On the 1008
%   calibration images without motion of experiments/calibrate_threshold.py
%   (white noise, random and camera textures, 32 to 256 px, L = 1 to 10,
%   PadFactor 2) the largest value was tanh(8.41/sqrt(m-5)); on all 2268
%   simulated images without motion, held-out textures included, it was
%   tanh(8.64/sqrt(m-5)), and none of them was valid with 9.
    if m <= c.NULL_OFFSET
        q = Inf;
    else
        q = tanh(c.NULL_K/sqrt(m - c.NULL_OFFSET));
    end
end
