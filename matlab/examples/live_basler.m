function live_basler(varargin)
%LIVE_BASLER  Live velocity estimate from a Basler camera (MATLAB only).
%   LIVE_BASLER('Mode', MODE, 'Name', Value, ...) streams from a GenTL
%   camera and estimates direction and speed of the moving surface. MODE
%   must be given:
%     'flash'    every camera frame is one pulsed image: a strobe fires
%                K = 2L+1 times (or 'Exposures') during one exposure. The
%                newest frame is estimated on its own. The speed is in px
%                per flash step, and in px/s if 'FlashPeriod' is given.
%     'average'  continuous light: the newest K consecutive frames are
%                averaged into a pulsed image (do not use this for
%                multi-flash frames). The speed is in px per frame and px/s.
%   Needs the Image Acquisition Toolbox with the GenICam/GenTL support
%   package (Octave has no camera interface; use the Python version
%   python/examples/live_camera.py there).
%
%   Options:
%     'Mode'          'flash' or 'average' (required)
%     'L'             half number of flashes or frames, K = 2L+1 (4)
%     'Exposures'     exact number K >= 2 of flashes or frames, overrides
%                     L ([])
%     'Duty'          flash mode: flash pulse length / flash period
%                     ([] = 0). Average mode: exposure time / frame step
%                     ([] = read the exposure time from the camera)
%     'FlashPeriod'   flash mode: time between flashes in s, for px/s ([])
%     'Exposure'      exposure time in us, switches auto exposure off ([]
%                     keeps the camera setting)
%     'FrameRate'     camera frame rate in Hz ([] keeps the camera setting)
%     'Tol'           average mode: tolerance of the frame-interval check
%                     described below (0.1)
%     'ClipOutliers'  clip static bright spots before the FFT (false)
%     'ROI'           side of the central square read from the sensor, px
%                     (512; one estimate then takes several seconds)
%     'Period'        seconds between estimates (0.5)
%     'Device'        index into the GenTL device list (1)
%
%   Frame times (average mode): GETDATA returns toolbox times, taken on the
%   host when a frame arrives, not camera timestamps (those would need
%   GenICam chunk data, which this script does not use). They carry
%   transfer jitter. If the camera reports its frame rate
%   (ResultingFrameRate or ResultingFrameRateAbs), a window is accepted
%   when its mean frame interval agrees with that period within 'Tol' and
%   no interval is off by half a period or more (a dropped frame); px/s
%   then uses the camera's period. This assumes that the camera runs freely
%   at that rate, not on a hardware trigger. Without a reported rate, or if
%   the mean interval does not agree with it, the spread of the intervals
%   (max - min) must be at most 'Tol' times their median. Other windows are
%   skipped. In flash mode the frame times only give the time axis; the
%   flash timing cannot be checked from them.
%
%   Exposure (average mode): each frame is blurred over its exposure time.
%   The script passes Duty = exposure / frame step to the estimator and
%   warns when it lies between 0.4 and 0.6, where the speed fit depends
%   strongly on the exact value, or above 0.9. Keep the exposure below
%   about 0.3 of the frame step.
%
%   Note: written against the documented Image Acquisition Toolbox API; it
%   has not been run in MATLAB or on camera hardware for this repository.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..'));
    opts = parse_options(varargin, struct('Mode', '', 'L', 4, 'Exposures', [], 'Duty', [], ...
        'FlashPeriod', [], 'Exposure', [], 'FrameRate', [], 'Tol', 0.1, 'ClipOutliers', false, ...
        'ROI', 512, 'Period', 0.5, 'Device', 1));
    if ~(ischar(opts.Mode) || isstring(opts.Mode)) || ~any(strcmpi(char(opts.Mode), {'flash', 'average'}))
        error('live_basler:Mode', ['Choose ''Mode'', ''flash'' (every frame holds several strobe ' ...
            'flashes) or ''Mode'', ''average'' (consecutive frames are averaged).']);
    end
    flash = strcmpi(char(opts.Mode), 'flash');
    L = check_option(opts.L, 'L', 'positive', 'live_basler:L');
    if isempty(opts.Exposures)
        K = 2*L + 1;
    else
        K = check_option(opts.Exposures, 'Exposures', 'exposures', 'live_basler:Exposures');
    end
    if ~isempty(opts.Duty)
        check_option(opts.Duty, 'Duty', 'duty', 'live_basler:Duty');
    end
    tol = check_option(opts.Tol, 'Tol', 'unit', 'live_basler:Tol');
    clip = check_option(opts.ClipOutliers, 'ClipOutliers', 'flag', 'live_basler:ClipOutliers');
    check_option(opts.ROI, 'ROI', 'roi', 'live_basler:ROI');
    check_option(opts.Device, 'Device', 'positive', 'live_basler:Device');
    check_positive(opts.FlashPeriod, 'FlashPeriod');
    check_positive(opts.Exposure, 'Exposure');
    check_positive(opts.FrameRate, 'FrameRate');
    check_positive(opts.Period, 'Period');
    if isempty(opts.Period)
        error('live_basler:Period', 'Period must be a positive number.');
    end
    if ~flash && ~isempty(opts.FlashPeriod)
        error('live_basler:FlashPeriod', ['FlashPeriod is for ''Mode'', ''flash''; in average ' ...
            'mode the frame times give px/s.']);
    end
    if flash
        nwin = 1;
        unit = 'px/flash';
    else
        nwin = K;
        unit = 'px/frame';
    end

    % --- camera -----------------------------------------------------------
    info = imaqhwinfo('gentl');
    if isempty(info.DeviceIDs)
        error('live_basler:nocamera', 'No GenTL camera found.');
    end
    dev = info.DeviceIDs{opts.Device};
    formats = info.DeviceInfo(opts.Device).SupportedFormats;
    k = find(~cellfun(@isempty, strfind(formats, 'Mono8')), 1);
    if isempty(k)
        error('live_basler:format', 'Camera offers no Mono8 format.');
    end
    vid = videoinput('gentl', dev, formats{k});
    vid.ReturnedColorspace = 'grayscale';
    vid.FramesPerTrigger = Inf;
    triggerconfig(vid, 'immediate');
    res = vid.VideoResolution;                         % [width height]
    side = min([opts.ROI, res]);
    vid.ROIPosition = [floor((res(1) - side)/2), floor((res(2) - side)/2), side, side];
    src = getselectedsource(vid);
    if ~isempty(opts.FrameRate)
        set_first(src, {'AcquisitionFrameRateEnable'}, 'True');
        if isempty(set_first(src, {'AcquisitionFrameRate', 'AcquisitionFrameRateAbs'}, opts.FrameRate))
            warning('live_basler:framerate', ['Could not set the frame rate (neither ' ...
                'AcquisitionFrameRate nor AcquisitionFrameRateAbs).']);
        end
    end
    if ~isempty(opts.Exposure)
        set_first(src, {'ExposureAuto'}, 'Off');
        if isempty(set_first(src, {'ExposureTime', 'ExposureTimeAbs'}, opts.Exposure))
            warning('live_basler:exposure', ['Could not set the exposure time (neither ' ...
                'ExposureTime nor ExposureTimeAbs).']);
        end
    elseif ~any(strcmpi(char(get_first(src, {'ExposureAuto'})), {'', 'Off'}))
        warning('live_basler:autoexposure', ['Auto exposure is on, so the exposure time can ' ...
            'change; set ''Exposure''.']);
    end
    cam_period = [];                                   % nominal frame period from the camera, s

    % --- results and figure -------------------------------------------------
    t_est = [];
    angle_est = [];
    speed_est = [];
    speed_pxs = [];
    valid_est = false(0);
    t_start = [];
    warned = {};

    fig = figure('Name', 'Live velocity (Basler)', 'NumberTitle', 'off', ...
        'Color', 'w', 'CloseRequestFcn', @on_close);
    subplot(2, 2, 1); h_frame = imagesc(zeros(side)); axis image off; colormap(gray);
    title('Current frame');
    subplot(2, 2, 2); h_bp = imagesc(zeros(side)); axis image off;
    if flash
        title(sprintf('Log spectrum (%d flashes)', K));
    else
        title(sprintf('Pulsed image (%d frames)', K));
    end
    ax_wr = subplot(2, 2, [3 4]);
    h_wr = plot(ax_wr, 0:179, zeros(1, 180), 'k-');
    xlim(ax_wr, [0 180]); xlabel(ax_wr, '\theta (deg)'); ylabel(ax_wr, 'WR');
    h_btn = uicontrol('Style', 'togglebutton', 'String', 'Start', ...
        'Units', 'pixels', 'Position', [10 10 70 28], 'Callback', @on_toggle);

    tmr = timer('ExecutionMode', 'fixedSpacing', 'Period', opts.Period, ...
        'BusyMode', 'drop', 'TimerFcn', @on_timer, 'ErrorFcn', @on_timer_error);

    % --- callbacks ------------------------------------------------------------
    function on_toggle(btn, ~)
        if get(btn, 'Value')
            set(btn, 'String', 'Stop');
            t_est = []; angle_est = []; speed_est = []; speed_pxs = []; valid_est = false(0);
            t_start = [];
            start(vid);
            rate = get_first(src, {'ResultingFrameRate', 'ResultingFrameRateAbs'});
            if isnumeric(rate) && isscalar(rate) && rate > 0
                cam_period = 1/rate;
            else
                cam_period = [];
                rate = NaN;
            end
            fprintf('Camera: exposure %.0f us, frame rate %.2f Hz\n', 1e6*exposure_seconds(src), rate);
            if flash && ~isempty(opts.FlashPeriod)
                need = (K - 1 + duty_or_zero())*opts.FlashPeriod;
                if exposure_seconds(src) < need
                    warning('live_basler:flashes', ['The exposure (%.0f us) is shorter than %d ' ...
                        'flashes %.0f us apart need (%.0f us).'], 1e6*exposure_seconds(src), K, ...
                        1e6*opts.FlashPeriod, 1e6*need);
                end
            end
            start(tmr);
        else
            stop_acquisition();
            plot_results();
        end
    end

    function on_timer(~, ~)
        n = vid.FramesAvailable;
        if n < nwin
            return
        end
        [frames, ts] = getdata(vid, n);
        frames = double(frames(:, :, 1, end-nwin+1:end));
        frames = reshape(frames, size(frames, 1), size(frames, 2), nwin);
        ts = ts(end-nwin+1:end);
        if isempty(t_start)
            t_start = ts(1);
        end
        set(h_frame, 'CData', frames(:, :, end));
        if flash
            bp = frames(:, :, end);
            duty = duty_or_zero();
            if isempty(opts.FlashPeriod)
                step = NaN;
            else
                step = opts.FlashPeriod;
            end
        else
            [step, regular] = frame_step(ts, cam_period, tol);
            if ~regular
                dt = diff(ts);
                if isempty(cam_period)
                    nominal = 'unknown';
                else
                    nominal = sprintf('%.2f ms', 1e3*cam_period);
                end
                title(ax_wr, sprintf(['t = %.1f s: frame intervals %.2f-%.2f ms (camera period ' ...
                    '%s), skipped'], ts(end) - t_start, 1e3*min(dt), 1e3*max(dt), nominal));
                drawnow limitrate;
                return
            end
            if isempty(opts.Duty)
                duty = exposure_duty(exposure_seconds(src)/step);
            else
                duty = opts.Duty;
            end
            bp = mean(frames, 3);
        end
        est = estimate_velocity(bp, 'L', L, 'Exposures', opts.Exposures, 'Duty', duty, ...
            'ClipOutliers', clip, 'Diagnostics', true);
        if flash
            set(h_bp, 'CData', est.diag.spectrum);
        else
            set(h_bp, 'CData', bp);
        end
        set(h_wr, 'YData', est.diag.wr);
        if est.valid
            state = '';
        elseif est.isotropic
            state = ' (not valid: isotropic)';
        else
            state = ' (not valid)';
        end
        pxs = est.speed/step;
        if isnan(pxs)
            per_s = '';
        else
            per_s = sprintf(' = %.0f px/s', pxs);
        end
        title(ax_wr, sprintf('direction %.1f deg, speed %.2f %s%s, quality %.2f%s', ...
            est.angle_deg, est.speed, unit, per_s, est.quality, state));
        drawnow limitrate;
        t_est(end + 1) = ts(end) - t_start;
        angle_est(end + 1) = est.angle_deg;
        speed_est(end + 1) = est.speed;
        speed_pxs(end + 1) = pxs;
        valid_est(end + 1) = est.valid;
    end

    function duty = duty_or_zero()
        if isempty(opts.Duty)
            duty = 0;
        else
            duty = opts.Duty;
        end
    end

    function duty = exposure_duty(duty)
        % Duty = exposure / frame step for the estimator; warns once per risky range.
        if ~(duty > 0)
            duty = 0;
            return
        end
        if duty >= 0.4 && duty <= 0.6 && ~any(strcmp(warned, 'half'))
            warned{end + 1} = 'half';
            warning('live_basler:duty', ['The exposure is %.2f of the frame step. Between 0.4 ' ...
                'and 0.6 the speed fit depends strongly on the exact exposure (at 0.5 every ' ...
                'second spectral line vanishes); use an exposure below about 0.3 of the step.'], duty);
        elseif duty > 0.9 && ~any(strcmp(warned, 'full'))
            warned{end + 1} = 'full';
            warning('live_basler:duty', ['The exposure is %.2f of the frame step, so the frames ' ...
                'are strongly blurred and most estimates will not be valid; shorten the ' ...
                'exposure.'], duty);
        end
        duty = min(duty, 0.99);                        % the estimator needs Duty < 1
    end

    function on_timer_error(~, evt)
        stop_acquisition();
        set(h_btn, 'Value', 0, 'String', 'Start');
        warning('live_basler:timer', 'Acquisition stopped: %s', evt.Data.message);
    end

    function stop_acquisition()
        if isvalid(tmr) && strcmp(tmr.Running, 'on')
            stop(tmr);
        end
        if isvalid(vid) && isrunning(vid)
            stop(vid);
            flushdata(vid);
        end
    end

    function plot_results()
        if isempty(t_est)
            disp('No estimates recorded.');
            return
        end
        figure('Name', 'Velocity over time', 'Color', 'w');
        subplot(2, 1, 1);
        plot(t_est(valid_est), angle_est(valid_est), 'o', 'MarkerFaceColor', 'auto'); hold on;
        plot(t_est(~valid_est), angle_est(~valid_est), 'x'); hold off;
        ylim([0 180]); ylabel('direction (deg)'); grid on;
        legend('valid', 'not valid', 'Location', 'best');
        subplot(2, 1, 2);
        if all(isnan(speed_pxs))
            speed = speed_est;
            label = sprintf('speed (%s)', unit);
        else
            speed = speed_pxs;
            label = 'speed (px/s)';
        end
        plot(t_est(valid_est), speed(valid_est), 'o', 'MarkerFaceColor', 'auto'); hold on;
        plot(t_est(~valid_est), speed(~valid_est), 'x'); hold off;
        xlabel('time (s)'); ylabel(label); grid on;
    end

    function on_close(~, ~)
        stop_acquisition();
        if isvalid(tmr)
            delete(tmr);
        end
        if isvalid(vid)
            delete(vid);
        end
        delete(fig);
    end
end

function [step, regular] = frame_step(ts, period, tol)
%FRAME_STEP  Time step of a window of frame times and whether it is regular.
%   With the nominal camera PERIOD (s), host times are accepted when their
%   mean interval agrees with PERIOD within TOL and no interval is off by
%   half a period or more (a dropped frame); STEP is then PERIOD. Otherwise
%   STEP is the median interval and the spread of the intervals must be at
%   most TOL*STEP. Same rules as frame_step in python/src/fftvel/stream.py.
    dt = diff(ts(:));
    if isempty(dt)
        step = NaN;
        regular = false;
    elseif ~isempty(period) && abs(mean(dt)/period - 1) <= tol
        step = period;
        regular = all(abs(dt/period - 1) < 0.5);
    else
        step = median(dt);
        regular = step > 0 && (max(dt) - min(dt)) <= tol*step;
    end
end

function s = exposure_seconds(src)
%EXPOSURE_SECONDS  Exposure time of the camera in s (NaN if unknown).
    us = get_first(src, {'ExposureTime', 'ExposureTimeAbs'});  % USB3 and newer GigE, GigE ace classic
    if isnumeric(us) && isscalar(us) && us > 0
        s = double(us)*1e-6;
    else
        s = NaN;
    end
end

function value = get_first(src, names)
%GET_FIRST  Value of the first readable camera property in NAMES ([] if none).
    value = [];
    for k = 1:numel(names)
        try
            value = get(src, names{k});
            return
        catch
        end
    end
end

function name = set_first(src, names, value)
%SET_FIRST  Set the first writable camera property in NAMES; return its name ('' if none).
    name = '';
    for k = 1:numel(names)
        try
            set(src, names{k}, value);
            name = names{k};
            return
        catch
        end
    end
end

function check_positive(x, name)
%CHECK_POSITIVE  Error unless X is empty or a positive, finite real number.
    if ~isempty(x) && ~(isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0)
        error(['live_basler:' name], '%s must be a positive number or [].', name);
    end
end
