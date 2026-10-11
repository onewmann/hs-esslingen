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
%                     (1.4826*MAD) from the median to the median before the
%                     FFT, against static bright spots such as glints or hot
%                     pixels (default false)
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
%     in_range      false if the best speed lies below the measurable range
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
%     L             L as passed
%     exposures     number of exposures K used in the model
%     n             side of the analysed central square in pixels
%     npad          FFT length
%
%   Method: the motion multiplies the image spectrum by a Dirichlet kernel,
%   which shows up as straight lines perpendicular to v, spaced
%   Npad/|v| bins apart. The direction is the angle whose mean Radon
%   projection of the log spectrum has the largest total variation (ripple
%   profile WR). The speed is found by fitting the Dirichlet model for
%   2L+1 exposures to the projection along that direction.
%
%   The model, the ripple profile WR and GEBA follow Chapter 3 of J. Horn's
%   doctoral thesis (Universitaet Karlsruhe, 2007; full reference in
%   README.md). The speed fit, the validity threshold and the isotropy check
%   were added in this implementation.
%
%   Non-square images are cropped to the central square. They are never
%   resized, so the speed stays in camera pixels.
%
%   See also fftvel.pulsed_from_stack, fftvel.simulate_pulsed_image,
%   fftvel.plot_estimate, fftvel.angle_difference.

    opts = fftvel.internal.parse_options(varargin, struct('L', 4, 'Exposures', [], 'Duty', 0, ...
        'ROI', [], 'PadFactor', 2, 'MinQuality', 0.6, 'ClipOutliers', false, ...
        'Diagnostics', false));
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
    diagnostics = fftvel.internal.check_option(opts.Diagnostics, 'Diagnostics', 'flag', ...
        'estimate_velocity:Diagnostics');

    [S, n, npad] = fftvel.internal.log_spectrum(bp, pad, roi, clip);
    rmax = floor(npad/2) - 1;
    rfrac = 0.9;

    % Direction: maximum of the ripple profile WR(theta)
    theta = 0:179;
    [M, j] = fftvel.internal.mean_projection(S, theta, rmax);
    wr = fftvel.internal.ripple_profile(M, j, rfrac);
    [wrmax, i] = max(wr);
    angle = mod(theta(i) + fftvel.internal.parabolic_offset(wr, i, true), 180);
    if mean(wr) > 0
        geba = wrmax / mean(wr);
    else
        geba = 1;        % no ripple in any direction
    end

    % Speed: fit the Dirichlet model to the projection along that direction
    m = fftvel.internal.mean_projection(S, angle, rmax);
    pos = j >= 2 & j <= rfrac*rmax;
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
    isotropic = quality - quality_perp < 0.15 && abs(log(speed_perp/speed)) < 0.12;

    % Independent profile points: zero padding beyond 2 only interpolates
    npts = numel(jpos);
    m_eff = min(npts, floor(2*n*npts/npad));

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
    est.threshold = max(min_quality, null_quality(m_eff));
    est.valid = quality >= est.threshold && fit.in_range && ~isotropic;
    est.L = L;
    est.exposures = K;
    est.n = n;
    est.npad = npad;
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
        est.diag = d;
    end
end

function q = null_quality(m)
%NULL_QUALITY  Fit quality that images without motion stay below.
%   Under the null hypothesis the Fisher transform atanh(q) of the best
%   correlation scales with 1/sqrt(m-5) for a profile of m independent
%   points (four are taken by the cubic, one by the fit). On 1008 images
%   without motion (white noise, random and camera textures, 32 to 256 px,
%   L = 1 to 10, PadFactor 2) the largest value was tanh(8.6/sqrt(m-5));
%   9 adds a margin.
    if m <= 5
        q = Inf;
    else
        q = tanh(9/sqrt(m - 5));
    end
end
