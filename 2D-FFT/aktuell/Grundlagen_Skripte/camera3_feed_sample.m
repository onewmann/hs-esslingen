function camera3_feed_sample_10ms()
    % Kamera 3 auswählen
    camIdx = 3;
    cam    = webcam(camIdx);

    % Zell-Array zum Speichern gesampelter Frames
    samples = {};

    % Flag zum Beenden der Schleife
    stopLoop = false;

    % Figure mit KeyPressFcn und CloseRequestFcn
    hFig = figure( ...
        'Name', 'Live Feed – s = sample, q = quit', ...
        'NumberTitle', 'off', ...
        'KeyPressFcn', @keyHandler, ...
        'CloseRequestFcn', @closeFigure ...
    );

    % Zeitsteuerung: alle 10 ms ein neues Frame holen
    interval   = 0.01;     % 10 ms
    lastUpdate = tic;      % Startzeitpunkt

    % Haupt-Loop
    while ~stopLoop && ishandle(hFig)
        if toc(lastUpdate) >= interval
            frame = snapshot(cam);
            imshow(frame);
            title('Drücken Sie s zum Samplen, q zum Beenden');
            lastUpdate = tic;    % Zeit zurücksetzen
        end
        drawnow;                 % erlaubt KeyPressFcn-Abfrage
    end

    % Aufräumen und Samples in den Workspace geben
    clear cam;
    assignin('base', 'cameraSamples', samples);
    fprintf('\nBeendet. %d Frames gespeichert in ''cameraSamples''.\n', numel(samples));

    % --- Nested Functions ---

    function keyHandler(~, event)
        switch event.Key
            case 's'
                % Aktuelles Frame ins Zell-Array packen
                samples{end+1} = frame;  %#ok<AGROW>
                fprintf('Frame %d gesampelt\n', numel(samples));
            case 'q'
                % Beenden
                stopLoop = true;
                if ishandle(hFig), close(hFig); end
        end
    end

    function closeFigure(src, ~)
        stopLoop = true;
        delete(src);
    end
end

