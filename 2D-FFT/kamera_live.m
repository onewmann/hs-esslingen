function kamera_live()
% KAMERA_LIVE  Zeichnet 2L+1 Frames von Kamera 3 alle 10 ms auf und 
%              führt danach eine Kapitel 3.3-Analyse (Mittelung, 
%              Pixelreduktion, FFT) durch.
%
% Aufruf:
%   >> kamera_live
%
% Drücke im GUI auf „Start“, um die Erfassung zu beginnen, auf „Stop“, 
% um sie abzubrechen.

    %% 1. Einstellungen
    L         = 3;                  % Halbe Frame-Anzahl → 2L+1
    numFrames = 2*L + 1;            % z.B. 7 Bilder
    interval  = 0.01;               % 10 ms Pause zwischen Snapshots

    %% 2. Webcam initialisieren
    camIdx = 3;                     
    cam    = webcam(camIdx);
    % Optional: Auflösung setzen
    if any(strcmp(cam.AvailableResolutions,'640x480'))
        cam.Resolution = '640x480';
    end
    dims = sscanf(cam.Resolution,'%dx%d');
    w    = dims(1);                 
    h    = dims(2);

    %% 3. GUI erstellen
    hFig = figure('Name','Live-Erfassung Kapitel 3.3', ...
                  'NumberTitle','off', ...
                  'MenuBar','none','ToolBar','none', ...
                  'CloseRequestFcn',@onClose);

    % Live-View
    axCam = subplot(2,3,1,'Parent',hFig);
    hImg  = imshow(zeros(h,w,3,'uint8'),'Parent',axCam);
    title(axCam,'Live-View Kamera 3');

    % Platzhalter für Analyse-Resultate
    subplot(2,3,2,'Parent',hFig); title('Summiertes Bild');
    subplot(2,3,3,'Parent',hFig); title('Pyramid-Down');
    subplot(2,3,4,'Parent',hFig); title('FFT Full');
    subplot(2,3,5,'Parent',hFig); title('FFT Pyramid');
    subplot(2,3,6,'Parent',hFig); title('FFT Crop');

    % Start/Stop Buttons
    uicontrol(hFig,'Style','pushbutton','String','Start', ...
              'Units','normalized','Position',[0.30 0.01 0.15 0.07], ...
              'Callback',@onStart);
    uicontrol(hFig,'Style','pushbutton','String','Stop',  ...
              'Units','normalized','Position',[0.55 0.01 0.15 0.07], ...
              'Callback',@onStop);

    %% 4. Initialisierung
    samples   = zeros(h,w,3,numFrames,'uint8');
    frameCnt  = 0;
    isRunning = false;

    %% 5. Callback: Start
    function onStart(~,~)
        if isRunning, return; end
        frameCnt  = 0;
        samples   = zeros(h,w,3,numFrames,'uint8');
        isRunning = true;

        % Aufnahme-Schleife
        while isRunning && frameCnt < numFrames
            frameCnt = frameCnt + 1;
            img = snapshot(cam);
            set(hImg,'CData',img);      % Live-View aktualisieren
            drawnow;                    % GUI-Events verarbeiten
            samples(:,:,:,frameCnt) = img;
            pause(interval);
        end

        % Beenden und Analyse starten
        if frameCnt >= numFrames
            isRunning = false;
            close(hFig);
            analyzeSequence(samples, L);
        end
    end

    %% 6. Callback: Stop
    function onStop(~,~)
        isRunning = false;
    end

    %% 7. Callback: Fenster schließen
    function onClose(~,~)
        isRunning = false;
        clear cam;
        delete(hFig);
    end

    %% 8. Unterfunktion: Kapitel 3.3-Analyse
    function analyzeSequence(seqRGB, L)
        % seqRGB: h×w×3×(2L+1)
        [h0,w0,~,K] = size(seqRGB);
        % 1) In Graustufen umwandeln
        seqGray = zeros(h0,w0,K,'double');
        for k = 1:K
            seqGray(:,:,k) = rgb2gray(seqRGB(:,:,:,k));
        end
        % 2) Mittelung
        bp = mean(seqGray,3);

        % 3) Pixelreduktion
        downFactor = 0.5;
        bp_pyr     = imresize(bp, downFactor, 'bilinear');
        cropSize   = round(h0*downFactor);
        startIdx   = floor((h0-cropSize)/2)+1;
        bp_crop    = bp(startIdx:startIdx+cropSize-1, ...
                        startIdx:startIdx+cropSize-1);

        % 4) FFT-Spektren
        B_full   = fftshift(fft2(bp));      logB_full  = log(1+abs(B_full));
        B_pyr    = fftshift(fft2(bp_pyr));  logB_pyr   = log(1+abs(B_pyr));
        B_crop   = fftshift(fft2(bp_crop)); logB_crop  = log(1+abs(B_crop));

        % 5) Darstellung
        figure('Name','Analyse Kapitel 3.3','NumberTitle','off');

        subplot(2,3,1); imshow(mat2gray(bp));       title('Summiertes Bild');
        subplot(2,3,2); imshow(mat2gray(bp_pyr));   title('Pyramid-Down');
        subplot(2,3,3); imshow(mat2gray(bp_crop));  title('Zentraler Crop');
        subplot(2,3,4); imshow(mat2gray(logB_full));title('FFT Full');
        subplot(2,3,5); imshow(mat2gray(logB_pyr)); title('FFT Pyramid');
        subplot(2,3,6); imshow(mat2gray(logB_crop));title('FFT Crop');
    end
end
