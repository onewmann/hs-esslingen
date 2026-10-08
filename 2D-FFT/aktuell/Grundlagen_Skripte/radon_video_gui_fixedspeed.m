function radon_video_gui_fixedspeed()
% RADON_VIDEO_GUI_FIXEDSPEED
% Radon-Analyse aus einem Video mit fester Abspielgeschwindigkeit.
% - Period des Timers wird hart auf >= 1 ms begrenzt
% - Warnung zur Sub-ms-Period wird global deaktiviert bis zum Cleanup
% - Automatisches Stoppen und Auswertung am Videoende

  %% 1. Parameter
  L       = 4;              % Halbfenster für Ringpuffer-Mittelung
  bufLen  = 2*L + 1;        % Pufferlänge
  theta   = 0:179;          % Radon-Winkel
  N       = 256;            % Arbeitsauflösung
  padDim  = 2 * N;          % Zero-Padding für FFT
  clipLim = [0.01 0.99];    % Kontrast-Stretch-Limits

  % Feste Abspielgeschwindigkeit
  speedFactor = 0.5;        % 0.5 = halb so schnell, 1.0 = normal, 2.0 = doppelt so schnell

  %% 2. Video initialisieren
  [file, path] = uigetfile({'*.mp4;*.avi','Video-Dateien'}, 'Video wählen');
  if isequal(file,0), return; end
  vidObj = VideoReader(fullfile(path,file));

  %% 3. Warnungen gezielt deaktivieren (bis zum Cleanup)
  % Sichere alten Warnzustand und deaktiviere relevante IDs
  warnState_badPeriod   = warning('query','MATLAB:timer:badPeriod');
  warnState_altPrecision = warning('query','MATLAB:timer:periodPrecision'); % falls andere ID
  warning('off','MATLAB:timer:badPeriod');
  warning('off','MATLAB:timer:periodPrecision');

  %% 4. Ring-Puffer & Mess-Arrays
  frameBuf   = zeros(N, N, bufLen);
  bufIdx     = 1;
  count      = 0;
  timeSec    = [];
  beta_vals  = [];
  velocities = [];
  t0 = tic;

  %% 5. GUI aufsetzen
  hFig = figure('Name','Radon Video','NumberTitle','off', 'CloseRequestFcn', @onClose);

  subplot(2,2,1);
    hLive = imshow(zeros(N,N));
    title('Einzelbild');

  subplot(2,2,2);
    hSum = imshow(zeros(N,N));
    title('Summenbild');

  hAx = subplot(2,2,[3 4]);
    hPlot = plot(hAx, theta, zeros(size(theta)), 'LineWidth', 1.5);
    xlabel(hAx,'\theta (°)'); ylabel(hAx,'WR');
    title(hAx,'WR-Profil');
    grid(hAx,'on');
    set(hAx,'Visible','off'); set(hPlot,'Visible','off');

  hToggle = uicontrol('Style','togglebutton','String','Start', ...
                      'Position',[10 10 60 30],'Parent',hFig, ...
                      'Callback',@onToggle);

  %% 6. Timer konfigurieren (Period hart begrenzen)
  desiredFPS = vidObj.FrameRate * speedFactor;
  targetFPS  = min(desiredFPS, 1000);     % Begrenzung auf 1000 FPS
  period     = max(1/targetFPS, 0.001);   % niemals kleiner als 1 ms

  t = timer('ExecutionMode','fixedRate', 'Period', period, ...
            'TimerFcn', @onTimer, ...
            'ErrorFcn', @(~,e) disp(e.Data.message));

  %% Callback Start/Stop
  function onToggle(src,~)
    if get(src,'Value')
      set(src,'String','Stop');
      set(hAx,'Visible','on'); set(hPlot,'Visible','on');
      start(t);
    else
      set(src,'String','Start');
      if isvalid(t), stop(t); end
      plotResults();
    end
  end

  %% TimerFcn: Frame holen und verarbeiten
  function onTimer(~,~)
    if hasFrame(vidObj)
      frameRaw = readFrame(vidObj);
      % RGB -> Graustufen
      if size(frameRaw,3) == 3
        frameRaw = rgb2gray(frameRaw);
      end
    else
      % Video zu Ende -> automatisch stoppen und Ergebnisse anzeigen
      if isvalid(t), stop(t); end
      set(hToggle,'Value',0); set(hToggle,'String','Start');
      plotResults();
      return;
    end

    % Preprocessing
    I = im2double(frameRaw);
    I = imresize(I,[N N]);
    I = imadjust(I, stretchlim(I,clipLim), []);
    set(hLive,'CData',I);

    % Ringpuffer
    frameBuf(:,:,bufIdx) = I;
    bufIdx = bufIdx + 1; if bufIdx > bufLen, bufIdx = 1; end
    count  = min(count + 1, bufLen);

    if count == bufLen
      bp = mean(frameBuf,3);
      set(hSum,'CData',bp);

      bg    = imopen(bp, strel('disk', round(N/20)));
      bp_dc = max(bp - bg, 0);
      bp_dn = wiener2(bp_dc, [5 5]);

      Hwin = hann(N); H = Hwin * Hwin.';
      Ipad = padarray(bp_dn .* H, [(padDim-N)/2 (padDim-N)/2], 0, 'both');

      B      = fftshift(fft2(Ipad));
      logB   = log(1 + abs(B));
      [R, xp] = radon(logB, theta);

      WR      = sum(abs(diff(R,1,1)),1);
      [~, idx] = max(WR);
      beta    = theta(idx);

      proj   = R(:,idx);
      mask   = xp >= 0;
      xp_p   = xp(mask);
      proj_p = proj(mask);

      [~, locs] = findpeaks(proj_p, xp_p, 'MinPeakProminence', 0.05, 'MinPeakDistance', 3);

      locs_sub = locs;
      for k = 1:numel(locs)
        [~, ix] = min(abs(xp_p - locs(k)));
        if ix > 1 && ix < length(proj_p)
          y1 = proj_p(ix-1); y0 = proj_p(ix); y2 = proj_p(ix+1);
          locs_sub(k) = locs(k) + (y1 - y2) / (2 * (y1 - 2*y0 + y2));
        end
      end

      if numel(locs_sub) >= 2
        v_mag = padDim / mean(diff(sort(locs_sub)));
      else
        v_mag = NaN;
      end

      set(hPlot,'YData',WR);
      title(hAx, sprintf('β=%.1f°  v=%.1f px/frame', beta, v_mag));
      drawnow;

      elapsed = toc(t0);
      timeSec(end+1)    = elapsed;
      beta_vals(end+1)  = beta;
      velocities(end+1) = v_mag;
    end
  end

  %% Ergebnisse plotten
  function plotResults()
    if isempty(timeSec), return; end
    tmax = max(timeSec);
    if tmax <= 10
      xt = 0:1:ceil(tmax);
    elseif tmax <= 60
      xt = 0:5:ceil(tmax/5)*5;
    else
      xt = 0:10:ceil(tmax/10)*10;
    end

    figure('Name','Winkel über Zeit');
    plot(timeSec, beta_vals, '-o', 'LineWidth', 1.5);
    xlabel('Zeit (s)'); ylabel('β (°)');
    xlim([0 tmax]); ylim([0 179]);
    xticks(xt); grid on;

    figure('Name','Geschwindigkeit über Zeit');
    plot(timeSec, velocities, '-o', 'LineWidth', 1.5);
    xlabel('Zeit (s)'); ylabel('v (px/frame)');
    xlim([0 tmax]); xticks(xt); grid on;
  end

  %% Cleanup
  function onClose(~,~)
    % Timer sauber schließen
    if exist('t','var') && isvalid(t), stop(t); delete(t); end
    % Ursprünglichen Warnzustand wiederherstellen
    if ~isempty(warnState_badPeriod)
      warning(warnState_badPeriod.state, 'MATLAB:timer:badPeriod');
    else
      warning('on','MATLAB:timer:badPeriod');
    end
    if ~isempty(warnState_altPrecision)
      warning(warnState_altPrecision.state, 'MATLAB:timer:periodPrecision');
    else
      warning('on','MATLAB:timer:periodPrecision');
    end
    delete(hFig);
  end

end

