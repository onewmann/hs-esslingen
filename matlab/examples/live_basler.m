function live_basler(varargin)
%LIVE_BASLER  Live velocity estimate from a Basler camera (MATLAB only).
%   LIVE_BASLER() streams from the first GenTL camera, averages every 2L+1
%   consecutive frames into a pulsed image and estimates direction and speed.
%   Needs the Image Acquisition Toolbox with the GenICam/GenTL support
%   package (Octave has no camera interface; use the Python version
%   python/examples/live_camera.py there).
%
%   LIVE_BASLER('Name', Value, ...) options:
%     'L'          half number of frames per estimate, 2L+1 in total (4)
%     'ROI'        side of the central square read from the sensor, px (512)
%     'Period'     seconds between estimates (0.5)
%     'FrameRate'  camera frame rate in Hz, [] keeps the camera setting ([])
%     'Device'     index into the GenTL device list (1)
%
%   Frames come from the continuous camera stream with their timestamps, so
%   the time step between exposures is known and constant. Windows whose
%   frame intervals differ by more than 20 % are skipped, because the
%   pulsed-exposure model needs equal steps. The speed is reported in
%   px/frame and in px/s.
%
%   Note: written for and tested against the documented Image Acquisition
%   Toolbox API; it has not been run on hardware in this repository's CI.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..'));
    opts = parse_options(varargin, struct('L', 4, 'ROI', 512, 'Period', 0.5, ...
        'FrameRate', [], 'Device', 1));
    L = opts.L;
    nwin = 2*L + 1;

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
    if ~isempty(opts.FrameRate)
        src = getselectedsource(vid);
        try
            src.AcquisitionFrameRateEnable = 'True';
            src.AcquisitionFrameRate = opts.FrameRate;
        catch err
            warning('live_basler:framerate', 'Could not set the frame rate: %s', err.message);
        end
    end

    % --- results and figure -------------------------------------------------
    t_est = [];
    angle_est = [];
    speed_est = [];
    speed_pxs = [];
    valid_est = false(0);
    t_start = [];

    fig = figure('Name', 'Live velocity (Basler)', 'NumberTitle', 'off', ...
        'Color', 'w', 'CloseRequestFcn', @on_close);
    subplot(2, 2, 1); h_frame = imagesc(zeros(side)); axis image off; colormap(gray);
    title('Current frame');
    subplot(2, 2, 2); h_bp = imagesc(zeros(side)); axis image off;
    title(sprintf('Pulsed image (%d frames)', nwin));
    ax_wr = subplot(2, 2, [3 4]);
    h_wr = plot(ax_wr, 0:179, zeros(1, 180), 'k-');
    xlim(ax_wr, [0 180]); xlabel(ax_wr, '\theta (deg)'); ylabel(ax_wr, 'WR');
    h_btn = uicontrol('Style', 'togglebutton', 'String', 'Start', ...
        'Units', 'pixels', 'Position', [10 10 70 28], 'Callback', @on_toggle);

    tmr = timer('ExecutionMode', 'fixedSpacing', 'Period', opts.Period, ...
        'BusyMode', 'drop', 'TimerFcn', @on_timer, 'ErrorFcn', @on_timer_error);

    % --- callbacks ------------------------------------------------------------
    function on_toggle(src, ~)
        if get(src, 'Value')
            set(src, 'String', 'Stop');
            t_est = []; angle_est = []; speed_est = []; speed_pxs = []; valid_est = false(0);
            t_start = [];
            start(vid);
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
        frames = double(squeeze(frames(:, :, 1, end-nwin+1:end)));
        ts = ts(end-nwin+1:end);
        if isempty(t_start)
            t_start = ts(1);
        end
        dt = diff(ts);
        step = median(dt);
        set(h_frame, 'CData', frames(:, :, end));
        if max(abs(dt - step)) > 0.2*step
            title(ax_wr, sprintf('t = %.1f s: frame intervals not equal (%.1f-%.1f ms), skipped', ...
                ts(end) - t_start, 1e3*min(dt), 1e3*max(dt)));
            drawnow limitrate;
            return
        end
        bp = mean(frames, 3);
        est = estimate_velocity(bp, 'L', L, 'Diagnostics', true);
        set(h_bp, 'CData', bp);
        set(h_wr, 'YData', est.diag.wr);
        if est.valid
            state = '';
        else
            state = ' (not valid)';
        end
        title(ax_wr, sprintf('direction %.1f deg, speed %.2f px/frame = %.0f px/s, quality %.2f%s', ...
            est.angle_deg, est.speed, est.speed/step, est.quality, state));
        drawnow limitrate;
        t_est(end + 1) = ts(end) - t_start;
        angle_est(end + 1) = est.angle_deg;
        speed_est(end + 1) = est.speed;
        speed_pxs(end + 1) = est.speed/step;
        valid_est(end + 1) = est.valid;
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
        plot(t_est(valid_est), speed_pxs(valid_est), 'o', 'MarkerFaceColor', 'auto'); hold on;
        plot(t_est(~valid_est), speed_pxs(~valid_est), 'x'); hold off;
        xlabel('time (s)'); ylabel('speed (px/s)'); grid on;
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
