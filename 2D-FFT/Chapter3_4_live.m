function Chapter3_4_live()
% Chapter3_4_live  Live-Aufnahme von Kamera 3 und Selbstbeurteilung (Kap. 3.4)
%
%  - Klick auf Start: Nimmt numFrames Bilder im 10 ms-Rhythmus auf
%  - Anschließend: 2D-FFT, log-Spektrum, Radon, WR und GEBA berechnen
%  - Klick auf Stop oder Schließen: Abbruch der Aufnahme

    %% 1. Parameter
    numFrames = 1000;       % Anzahl Frames (z.B. 50)
    interval  = 0.005;    % 10 ms Pause zwischen den Snapshots

    %% 2. Kamera initialisieren
    camIdx = 3;
    cam    = webcam(camIdx);

    % Optional: auf 640×480 setzen, wenn unterstützt
    if any(strcmp(cam.AvailableResolutions,'640x480'))
        cam.Resolution = '640x480';
    end
    dims = sscanf(cam.Resolution, '%dx%d');
    w    = dims(1);
    h    = dims(2);

    %% 3. Live-GUI mit Start/Stop
    hFig = figure( ...
        'Name',           'Kapitel 3.4 Live-Erfassung', ...
        'NumberTitle',    'off', ...
        'MenuBar',        'none', ...
        'ToolBar',        'none', ...
        'CloseRequestFcn',@onClose ...
    );

    axLive = axes('Parent',hFig, 'Position',[0.05 0.2 0.9 0.75]);
    hImg   = imshow(zeros(h,w,3,'uint8'),'Parent',axLive);
    title(axLive,'Live-View Kamera 3');

    uicontrol( ...
        'Parent', hFig, ...
        'Style',  'pushbutton', ...
        'String', 'Start', ...
        'Units',  'normalized', ...
        'Position',[0.30 0.05 0.15 0.1], ...
        'Callback',@onStart ...
    );

    uicontrol( ...
        'Parent', hFig, ...
        'Style',  'pushbutton', ...
        'String', 'Stop', ...
        'Units',  'normalized', ...
        'Position',[0.55 0.05 0.15 0.1], ...
        'Callback',@onStop ...
    );

    %% 4. Puffer und Steuerflag
    samples   = zeros(h, w, 3, numFrames, 'uint8');
    frameCnt  = 0;
    isRunning = false;

    %% 5. Callback: Start
    function onStart(~,~)
        if isRunning, return; end
        frameCnt  = 0;
        samples   = zeros(h, w, 3, numFrames, 'uint8');
        isRunning = true;

        while isRunning && frameCnt < numFrames
            frameCnt = frameCnt + 1;
            img = snapshot(cam);
            set(hImg, 'CData', img);
            drawnow;                         % UI-Events abarbeiten
            samples(:,:,:,frameCnt) = img;   % speichern
            pause(interval);
        end

        % Aufnahme fertig oder gestoppt
        isRunning = false;
        if frameCnt >= numFrames
            close(hFig);
            analyzeChapter3_4(samples);
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
end


%% Unterfunktion: Kapitel 3.4 Analyse
function analyzeChapter3_4(seqRGB)
    % seqRGB: h×w×3×K uint8
    [h, w, ~, K] = size(seqRGB);

    % 1) In Graustufen konvertieren
    seqGray = zeros(h, w, K, 'double');
    for k = 1:K
        seqGray(:,:,k) = rgb2gray(seqRGB(:,:,:,k));
    end

    % 2) Summiertes Bild (arithm. Mittel)
    bp = mean(seqGray, 3);

    % 3) 2D-FFT & log-Spektrum
    B    = fftshift(fft2(bp));
    logB = log(1 + abs(B));

    % 4) Radon-Transformation
    theta = 0:179;
    [R, xp] = radon(logB, theta);

    % 5) WR-Profil berechnen
    WR = zeros(1, numel(theta));
    for i = 1:numel(theta)
        gradProj = abs(diff(R(:,i)));
        WR(i)    = sum(gradProj);
    end

    % 6) Winkel & GEBA
    [WR_max, idx] = max(WR);
    v_angle_est   = theta(idx);
    GEBA          = WR_max / mean(WR);

    % 7) Ergebnisse plotten
    figure('Name','Kapitel 3.4 Analyse','NumberTitle','off','Position',[200 200 900 400]);

    subplot(2,3,1);
    imshow(mat2gray(bp));
    title('Summiertes Bild bp(x;v)');

    subplot(2,3,2);
    imshow(mat2gray(logB));
    title('Log-FFT Spektrum');

    subplot(2,3,3);
    imagesc(theta, xp, R);
    colormap(hot);
    xlabel('Winkel (°)');
    ylabel('r');
    title('Radon\{logB\}');
    colorbar;

    subplot(2,3,4);
    plot(theta, WR, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('Winkel (°)');
    ylabel('WR');
    title(sprintf('WR-Profil (θ* = %.1f°, GEBA = %.2f)', v_angle_est, GEBA));

    subplot(2,3,5);
    text(0.1, 0.6, sprintf('Geschätzter Winkel: %.1f°', v_angle_est), 'FontSize', 12);
    text(0.1, 0.3, sprintf('GEBA = %.2f', GEBA),               'FontSize', 12);
    axis off;

    subplot(2,3,6);
    axis off;
end
