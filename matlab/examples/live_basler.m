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
%     'SaveDir'       folder to save every estimated window in, for
%                     experiments/real_measurements.py ([] = do not save)
%     'Note'          with 'SaveDir': text for the note column of meta.csv
%                     ('')
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
%   Saturation: the script warns once when more than 1 % of the analysed
%   pixels are saturated (EST.saturated), because glints and saturated
%   areas can give wrong results that are marked valid.
%
%   Saving ('SaveDir'): every window taken for an estimate (K frames in
%   average mode, also those skipped for uneven frame intervals; one frame
%   in flash mode) is saved as burst_0001.mat or frame_0001.mat, numbered on
%   from the files already in the folder, with the variables frames (as
%   GETDATA returns them, H x W x 1 x K, in the class of the camera), t
%   (the toolbox times of the frames in s), exposure_s (NaN if unknown),
%   frame_rate (Hz, if the camera reports it) and clock ('host'). A line
%   for each file is appended to meta.csv in the folder (written with a
%   header first if it does not exist), with the mode, the number of
%   flashes or frames, the exposure, the periods known here and 'Note'
%   filled in; in average mode frame_period_s is the camera's period,
%   since the toolbox times jitter. Fill in gsd_mm_per_px, set_speed_mm_s,
%   set_angle_deg and static before running
%   python experiments/real_measurements.py FOLDER/meta.csv
%   (data/real/README.md). An existing meta.csv keeps its header and
%   separator; with semicolons (German Excel) the numbers get decimal
%   commas. If meta.csv cannot be written (open in Excel), the line is
%   printed instead.
%
%   Note: written against the documented Image Acquisition Toolbox API; it
%   has not been run in MATLAB or on camera hardware for this repository.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..'));
    cst = fftvel.internal.constants();
    opts =fftvel.internal.parse_options(varargin, struct('Mode', '', 'L', 4, 'Exposures', [], ...
        'Duty', [], 'FlashPeriod', [], 'Exposure', [], 'FrameRate', [], 'Tol', 0.1, ...
        'ClipOutliers', false, 'ROI', 512, 'Period', 0.5, 'Device', 1, 'SaveDir', [], 'Note', ''));
    if ~(ischar(opts.Mode) || isstring(opts.Mode)) || ~any(strcmpi(char(opts.Mode), {'flash', 'average'}))
        error('live_basler:Mode', ['Choose ''Mode'', ''flash'' (every frame holds several strobe ' ...
            'flashes) or ''Mode'', ''average'' (consecutive frames are averaged).']);
    end
    flash = strcmpi(char(opts.Mode), 'flash');
    L = fftvel.internal.check_option(opts.L, 'L', 'positive', 'live_basler:L');
    if isempty(opts.Exposures)
        K = 2*L + 1;
    else
        K = fftvel.internal.check_option(opts.Exposures, 'Exposures', 'exposures', ...
            'live_basler:Exposures');
    end
    if ~isempty(opts.Duty)
        fftvel.internal.check_option(opts.Duty, 'Duty', 'duty', 'live_basler:Duty');
    end
    tol = fftvel.internal.check_option(opts.Tol, 'Tol', 'unit', 'live_basler:Tol');
    clip = fftvel.internal.check_option(opts.ClipOutliers, 'ClipOutliers', 'flag', ...
        'live_basler:ClipOutliers');
    fftvel.internal.check_option(opts.ROI, 'ROI', 'roi', 'live_basler:ROI');
    fftvel.internal.check_option(opts.Device, 'Device', 'positive', 'live_basler:Device');
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
        stem = 'frame';
    else
        nwin = K;
        unit = 'px/frame';
        stem = 'burst';
    end
    if ~(ischar(opts.Note) || (isstring(opts.Note) && isscalar(opts.Note)))
        error('live_basler:Note', 'Note must be text.');
    end
    note = char(opts.Note);
    save_dir = '';
    if ~isempty(opts.SaveDir)
        if ~((ischar(opts.SaveDir) && size(opts.SaveDir, 1) == 1) || ...
                (isstring(opts.SaveDir) && isscalar(opts.SaveDir)))
            error('live_basler:SaveDir', 'SaveDir must be the name of a folder, or [].');
        end
        save_dir = char(opts.SaveDir);
        if exist(save_dir, 'file') && ~exist(save_dir, 'dir')
            error('live_basler:SaveDir', 'SaveDir ''%s'' is a file, not a folder.', save_dir);
        end
    elseif ~isempty(note)
        error('live_basler:Note', 'Note needs SaveDir.');
    end

    n_saved = 0;                                       % number of the last file in save_dir
    saved_now = 0;                                     % files saved since the last Start
    if ~isempty(save_dir)
        if ~exist(save_dir, 'dir') && ~mkdir(save_dir)
            error('live_basler:SaveDir', 'Cannot create the folder ''%s''.', save_dir);
        end
        n_saved = last_number(save_dir, stem);
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
            saved_now = 0;
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
            if saved_now > 0
                fprintf(['Saved %d files to %s. Fill in gsd_mm_per_px, set_speed_mm_s, set_angle_deg ' ...
                    'and static in its meta.csv, then run\n  python experiments/real_measurements.py %s\n'], ...
                    saved_now, save_dir, fullfile(save_dir, 'meta.csv'));
            end
        end
    end

    function on_timer(~, ~)
        n = vid.FramesAvailable;
        if n < nwin
            return
        end
        [frames, ts] = getdata(vid, n);
        frames = frames(:, :, :, end-nwin+1:end);       % H x W x 1 x nwin, as getdata returns them
        ts = ts(end-nwin+1:end);
        saved = '';
        if ~isempty(save_dir)
            saved = [' [' save_window(frames, ts) ']'];
        end
        frames = double(frames);
        frames = reshape(frames, size(frames, 1), size(frames, 2), nwin);
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
                    '%s), skipped%s'], ts(end) - t_start, 1e3*min(dt), 1e3*max(dt), nominal, saved));
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
        est = fftvel.estimate_velocity(bp, 'L', L, 'Exposures', opts.Exposures, 'Duty', duty, ...
            'ClipOutliers', clip, 'Diagnostics', true);
        if est.saturated > cst.SAT_WARN && ~any(strcmp(warned, 'saturated'))
            warned{end + 1} = 'saturated';
            warning('live_basler:saturated', ['%.1f %% of the analysed pixels are saturated ' ...
                '(more than %g %%); saturated areas and glints can give wrong results that are ' ...
                'marked valid. Lower the exposure or the flash energy.'], 100*est.saturated, ...
                100*cst.SAT_WARN);
        end
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
        title(ax_wr, sprintf('direction %.1f deg, speed %.2f %s%s, quality %.2f%s%s', ...
            est.angle_deg, est.speed, unit, per_s, est.quality, state, saved));
        drawnow limitrate;
        t_est(end + 1) = ts(end) - t_start;
        angle_est(end + 1) = est.angle_deg;
        speed_est(end + 1) = est.speed;
        speed_pxs(end + 1) = pxs;
        valid_est(end + 1) = est.valid;
    end

    function name = save_window(frames, ts)
        % Save one window as a .mat file that real_measurements.py reads, and add its line to meta.csv.
        n_saved = n_saved + 1;
        saved_now = saved_now + 1;
        name = sprintf('%s_%04d.mat', stem, n_saved);
        s = struct('frames', frames, 't', ts, 'exposure_s', exposure_seconds(src), 'clock', 'host');
        if ~isempty(cam_period)
            s.frame_rate = 1/cam_period;
        end
        save(fullfile(save_dir, name), '-struct', 's', '-v7');
        row = struct('file', name, 'mode', 'flash', 'exposures', '', 'frames', '', 'flash_period_s', '', ...
            'frame_period_s', '', 'duty', '', 'exposure_s', '', 'note', note);
        if flash
            row.exposures = sprintf('%d', K);
            if ~isempty(opts.FlashPeriod)
                row.flash_period_s = sprintf('%.6g', opts.FlashPeriod);
            end
        else
            row.mode = 'average';
            row.frames = sprintf('%d', K);
            if ~isempty(cam_period)                    % toolbox times jitter: the camera's period is better
                row.frame_period_s = sprintf('%.6g', cam_period);
            end
        end
        if ~isempty(opts.Duty)
            row.duty = sprintf('%.6g', opts.Duty);
        end
        if s.exposure_s > 0
            row.exposure_s = sprintf('%.6g', s.exposure_s);
        end
        append_meta(save_dir, row);
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

function k = last_number(folder, stem)
%LAST_NUMBER  Highest number k of the files STEM_k.mat in FOLDER (0 if there are none).
    k = 0;
    files = dir(fullfile(folder, [stem '_*.mat']));
    for j = 1:numel(files)
        tok = regexp(files(j).name, ['^' stem '_(\d+)\.mat$'], 'tokens', 'once');
        if ~isempty(tok)
            k = max(k, str2double(tok{1}));
        end
    end
end

function append_meta(folder, row)
%APPEND_META  Append ROW (a struct of text) as a line to FOLDER/meta.csv, in the order of its header.
%   A new meta.csv gets the header of the columns below. An existing one
%   keeps its header and separator (comma, semicolon or tab); with
%   semicolons or tabs the numbers get decimal commas, as German Excel
%   writes them. If its header lacks a column that ROW fills, or the file
%   cannot be written, the line is printed instead.
    columns = {'file', 'mode', 'exposures', 'frames', 'flash_period_s', 'frame_period_s', 'duty', ...
        'exposure_s', 'gsd_mm_per_px', 'set_speed_mm_s', 'set_angle_deg', 'static', 'roi', 'note'};
    file = fullfile(folder, 'meta.csv');
    header = columns;
    delim = ',';
    info = dir(file);
    fresh = isempty(info) || info.bytes == 0;
    ends_with_newline = true;
    if ~fresh
        fid = fopen(file, 'r');
        if fid < 0
            print_line(row, columns, sprintf('cannot read %s', file));
            return
        end
        bytes = fread(fid, Inf, 'uint8=>char').';    % any encoding (Excel: Windows-1252), byte by byte
        fclose(fid);
        ends_with_newline = bytes(end) == char(10);
        breaks = [0, find(bytes == char(10) | bytes == char(13)), numel(bytes) + 1];
        first = '';                                    % first line that is not blank
        for j = 1:numel(breaks) - 1
            first = bytes(breaks(j) + 1:breaks(j + 1) - 1);
            if any(~isspace(first))
                break
            end
            first = '';
        end
        fresh = isempty(first);                        % only blank lines: write the header
    end
    if ~fresh
        if numel(first) >= 3 && isequal(double(first(1:3)), [239 187 191])
            first = first(4:end);                      % UTF-8 byte order mark
        end
        delims = {',', ';', char(9)};
        [~, k] = max(cellfun(@(d) sum(first == d), delims));
        delim = delims{k};
        header = lower(strtrim(strsplit(first, delim)));
        names = fieldnames(row);
        filled = names(~cellfun(@isempty, struct2cell(row)));
        missing = filled(~ismember(filled, header));
        if ~isempty(missing)
            print_line(row, columns, sprintf('the header of %s has no column %s', file, ...
                strjoin(missing.', ', ')));
            return
        end
        if ~strcmp(delim, ',')                         % semicolons or tabs: decimal commas, as Excel writes them
            for name = {'exposures', 'frames', 'flash_period_s', 'frame_period_s', 'duty', 'exposure_s'}
                row.(name{1}) = strrep(row.(name{1}), '.', ',');
            end
        end
    end
    fid = fopen(file, 'a');
    if fid < 0
        print_line(row, columns, sprintf('cannot write %s (open in Excel?)', file));
        return
    end
    if ~ends_with_newline
        fprintf(fid, '\n');
    end
    if fresh
        fprintf(fid, '%s\n', strjoin(columns, ','));
    end
    fprintf(fid, '%s\n', csv_line(row, header, delim));
    fclose(fid);
end

function print_line(row, columns, why)
%PRINT_LINE  Print the meta.csv line of ROW (with the header) when it cannot be appended.
    warning('live_basler:meta', '%s; add this line to meta.csv yourself:\n%s\n%s', why, ...
        strjoin(columns, ','), csv_line(row, columns, ','));
end

function line = csv_line(row, header, delim)
%CSV_LINE  Values of ROW in the order of HEADER, separated by DELIM; quoted where needed.
    cells = repmat({''}, 1, numel(header));
    for k = 1:numel(header)
        if isfield(row, header{k})
            c = row.(header{k});
            if any(c == delim) || any(c == '"') || any(c == char(10))
                c = ['"' strrep(c, '"', '""') '"'];
            end
            cells{k} = c;
        end
    end
    line = strjoin(cells, delim);
end

function check_positive(x, name)
%CHECK_POSITIVE  Error unless X is empty or a positive, finite real number.
    if ~isempty(x) && ~(isnumeric(x) && isscalar(x) && isreal(x) && isfinite(x) && x > 0)
        error(['live_basler:' name], '%s must be a positive number or [].', name);
    end
end
