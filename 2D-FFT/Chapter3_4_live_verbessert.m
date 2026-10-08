function Chapter3_4_live()
% Chapter3_4_live: Live-Erfassung von Kamera 3 in Graustufen und Selbstbeurteilung (Kap. 3.4)

    %% 1. Konfiguration
    numFrames = 20;
    interval  = 0.2;
    dsFactor  = 1;

    %% 2. Kamera auswählen (robust)
    camList = webcamlist;
    nCams   = numel(camList);
    if nCams == 0
        error('Keine Kamera erkannt.');
    end
    camIdx = min(2,nCams);
    cam    = webcam(camIdx);
    if any(strcmp(cam.AvailableResolutions,'640x480'))
        cam.Resolution = '640x480';
    end
    dims = sscanf(cam.Resolution,'%dx%d');  w = dims(1); h = dims(2);

    %% 3. GUI aufbauen
    hFig = figure('Name','Kapitel 3.4 Live & GEBA','NumberTitle','off', ...
        'MenuBar','none','ToolBar','none','CloseRequestFcn',@onClose);
    axLive = axes('Parent',hFig,'Position',[0.05 0.2 0.9 0.75]);
    hImg   = imshow(zeros(h,w,'uint8'),'Parent',axLive);
    title(axLive,'Live-View (Graustufen)');

    uicontrol(hFig,'Style','pushbutton','String','Start', ...
        'Units','normalized','Position',[0.30 0.05 0.15 0.08],'Callback',@onStart);
    uicontrol(hFig,'Style','pushbutton','String','Stop', ...
        'Units','normalized','Position',[0.55 0.05 0.15 0.08],'Callback',@onStop);

    %% 4. Initialisierung
    samples   = zeros(h,w,numFrames,'uint8');
    frameCnt  = 0;
    isRunning = false;

    %% Start-Callback
    function onStart(~,~)
        frameCnt = 0; isRunning = true;
        samples  = zeros(h,w,numFrames,'uint8');

        while isRunning && frameCnt < numFrames
            frameCnt = frameCnt + 1;
            rgbFrame = snapshot(cam);
            grayFrame = rgb2gray(rgbFrame);
            samples(:,:,frameCnt) = grayFrame;
            set(hImg,'CData',grayFrame); drawnow;
            pause(interval);
        end

        isRunning = false;
        if frameCnt >= numFrames
            close(hFig);
            analyzeGrayFrames(samples, dsFactor);
        end
    end

    %% Stop-Callback
    function onStop(~,~)
        isRunning = false;
    end

    %% Close-Callback
    function onClose(~,~)
        isRunning = false;
        clear cam;
        delete(hFig);
    end

    %% Analyse-Funktion (Graustufen)
    function analyzeGrayFrames(frames, dsF)
        [h0,w0,K] = size(frames);
        if dsF > 1
            frames = frames(1:dsF:end,1:dsF:end,:);
            [h0,w0,~] = size(frames);
        end

        bp = mean(double(frames),3);                   % Mittelbild
        bp = bp - mean(bp(:));                         % DC-Subtraktion
        win2d = hann(h0) * hann(w0)';
        bp_win = bp .* win2d;

        P = 2^nextpow2(max(h0,w0));
        bp_pad = padarray(bp_win, [P-h0, P-w0], 0, 'post');
        B      = fftshift(fft2(bp_pad));
        logB   = log(1 + abs(B));

        theta = 0:179;
        [R, xp] = radon(logB, theta);
        WR = sum(abs(diff(R,1,1)),1);

        [WRmax, idxMax] = max(WR);
        v_angle_est = theta(idxMax);
        GEBA        = WRmax / mean(WR);

        figure('Name','Kapitel 3.4 Analyse','NumberTitle','off');
        subplot(2,3,1); imshow(mat2gray(bp)); title('Summiertes Graubild');
        subplot(2,3,2); imshow(mat2gray(logB)); title('Log-FFT Spektrum');
        subplot(2,3,3); imagesc(theta,xp,R); colormap(hot); colorbar;
        xlabel('Winkel (°)'); ylabel('r'); title('Radon\{logB\}');
        subplot(2,3,4); plot(theta,WR,'b-','LineWidth',1.5); grid on;
        xlabel('Winkel (°)'); ylabel('WR');
        title(sprintf('WR-Profil (θ*=%.1f°, GEBA=%.2f)', v_angle_est, GEBA));
        subplot(2,3,5); axis off;
        text(0.1,0.6,sprintf('θ* = %.1f°', v_angle_est),'FontSize',12);
        text(0.1,0.3,sprintf('GEBA = %.2f',GEBA),'FontSize',12);
        subplot(2,3,6); axis off;
    end
end
