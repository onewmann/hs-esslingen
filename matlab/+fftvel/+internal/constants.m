function c = constants()
%CONSTANTS  Numeric constants of the estimator and the simulator.
%   C = fftvel.internal.constants() returns a struct with one field per
%   constant. The names and values are those of the module constants in
%   python/src/fftvel/core.py and of tests/constants.json, which also
%   gives the meaning of each; tests/matlab/test_constants.m and
%   tests/python/test_constants.py check all three against each other.
%   The struct is built on the first call and kept, so a function fetches
%   it once at its start, never inside a loop.

    persistent cache
    if isempty(cache)
        c = struct();

        % Defaults of fftvel.estimate_velocity
        c.DEFAULT_L = 4;                 % half the number of exposures, K = 2L+1 = 9
        c.DEFAULT_PAD_FACTOR = 2;        % zero padding factor of the FFT
        c.DEFAULT_MIN_QUALITY = 0.6;     % lowest fit quality accepted for valid

        % Input checks
        c.MIN_SIZE = 16;                 % smallest analysed square (and ROI) in px
        c.MIN_EXPOSURES = 2;             % smallest number of exposures K

        % Log spectrum
        c.LOG_FLOOR = 0.05;              % S = log(|B| + LOG_FLOOR*median|B|), mean|B| if the median is 0
        c.BLANK_FLOOR = 1.0;             % floor instead when |B| is 0 everywhere (blank image)
        c.CLIP_SIGMA = 6;                % ClipOutliers: robust standard deviations from the median ...
        c.MAD_SCALE = 1.4826;            % ... with robust standard deviation = MAD_SCALE*MAD
        c.CLIP_RADIUS = 3;               % ... and the pixels within this distance (px) of such a pixel

        % Direction
        c.ANGLE_STEP = 1.0;              % angle grid 0:ANGLE_STEP:180-ANGLE_STEP degrees
        c.RFRAC = 0.9;                   % fraction of the spectrum radius used for WR and the fit
        c.WEIGHT_FLOOR = 1e-12;          % mean_projection: lower bound of the bin weights

        % Speed
        c.PROFILE_START = 2;             % first bin of the half profile used for the speed fit
        c.SPEED_STEP = 0.005;            % step of the logarithmic speed grid
        c.VMIN_FACTOR = 1.5;             % vmin = VMIN_FACTOR*npad/(K*max(jpos))
        c.MIN_LINE_SPACING = 4;          % vmax = npad/MIN_LINE_SPACING
        c.GUARD_FACTOR = 4;              % the speed grid extends down to vmin/GUARD_FACTOR
        c.TREND_DEGREE = 3;              % degree of the polynomial that absorbs the texture spectrum
        c.TEMPLATE_FLOOR = 1e-3;         % dirichlet_template and box_template return 0.5*log(P + TEMPLATE_FLOOR)
        c.FLAT_MODEL = 1e-6;             % fit_speed, fit_blur: score 0 where the detrended model is flat

        % Validity
        c.NULL_K = 9;                    % null quality tanh(NULL_K/sqrt(m - NULL_OFFSET)) ...
        c.NULL_OFFSET = 5;               % ... TREND_DEGREE + 1 points for the polynomial, 1 for the fit
        c.INDEPENDENT_PAD = 2;           % m = min(points, floor(INDEPENDENT_PAD*n*points/npad))
        c.ISO_GAP = 0.15;                % isotropic if quality - quality_perp is below this ...
        c.ISO_LOG_SPEED = 0.12;          % ... and |log(speed_perp/speed)| is below this
        c.BLUR_STEP = 0.5;               % px between the samples of the continuous blur in box_template
        c.BLUR_GAP = 0.0;                % blur_like if quality - quality_blur is below this (and slow)

        % Saturation
        c.SAT_SHARE = 0.005;             % saturated: least share of pixels at the maximum (unless intmax)
        c.SAT_WARN = 0.01;               % front ends warn above this share

        % Simulation and grey conversion of colour input
        c.SIM_MARGIN = 32;               % extra texture size in px beyond N + (K-1+Duty)*max|V|
        c.DEFAULT_CORR_LEN = 1.5;        % correlation length of the random texture in px
        c.GREY_WEIGHTS = [0.2989, 0.5870, 0.1140];   % R, G, B weights of rgb2gray

        cache = c;
    end
    c = cache;
end
