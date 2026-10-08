function radon_video_gui_withnoise()
% RADON_VIDEO_GUI_WITHNOISE
% Radon-Analyse aus einem Video mit fester Abspielgeschwindigkeit.
% Erweiterung: optionales Hinzufügen von Rauschen ins Bild oder in die Radon-Projektionen.
%
% Nach Stop oder Video-Ende werden zwei Diagramme geplottet:
%   • Winkel β über Zeit
%   • Geschwindigkeit v über Zeit

  %% 1. Parameter
  L       = 4;              % Halbfenster für Ringpuffer-Mittelung
  bufLen  = 2*L + 1;        % Pufferlänge
  theta   = 0:179;          % Radon-Winkel
  N       = 256;            % Arbeitsauflösung
  padDim  = 2 * N;          % Zero-Padding für FFT
  clipLim = [0.01 0.99];    % Kontrast-Stretch-Limits

  % Feste Abspielgeschwindigkeit
  speedFactor = 0.5;        % 0.5 = halb so schnell, 1.0 = normal, 2.0 = doppelt so schnell

  % --- NEU: Rausch-Optionen ---
  applyImageNoise = false;   % true = Rauschen ins Bild
  applyRadonNoise = false;  % true = Rauschen in Radon-Projektionen
  noiseLevelImage = 10.0;   % Standardabweichung für Bildrauschen
  noiseLevelRadon = 90.0;   % Standardabweichung für Radonrauschen

  %% 2. Video initialisieren
  [file, path] = uigetfile({'*.mp4;*.avi','Video-Dateien'}, 'Video wählen');
  if isequal(file,0), return; end
  vidObj = VideoReader(fullfile(path,file));

  %% 3. Ring-Puffer & Mess-Arrays
  frameBuf   = zeros(N, N, bufLen);
  bufIdx     = 1;
  count      = 0;
  timeSec    = [];
  beta_vals  = [];
  velocities = [];
  t0 = tic;

  %% 4. GUI aufsetzen
  hFig = figure('Name','Radon Video mit Rauschen','NumberTitle','off', ...
                'CloseRequestFcn', @onClose);

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

  %% 5. Timer konfigurieren
  desiredFPS = vidObj.FrameRate * speedFactor;
  targetFPS  = min(desiredFPS, 1000);
  period     = max(1/targetFPS, 0.001);

  warning('off','MATLAB:timer:badPeriod');
  t = timer('ExecutionMode','fixedRate', 'Period', period, ...
            'TimerFcn', @onTimer, ...
            'ErrorFcn', @(~,e) disp(e.Data.message));
  warning('on','MATLAB:timer:badPeriod');

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
      if size(frameRaw,3) == 3
        frameRaw = rgb2gray(frameRaw);
      end
    else
      if isvalid(t), stop(t); end
      set(hToggle,'Value',0); set(hToggle,'String','Start');
      plotResults();
      return;
    end

    % Preprocessing
    I = im2double(frameRaw);
    I = imresize(I,[N N]);
    I = imadjust(I, stretchlim(I,clipLim), []);

    % --- NEU: Bildrauschen ---
    if applyImageNoise
        I = imnoise(I, 'gaussian', 0, noiseLevelImage^2);
    end

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

      % --- NEU: Radonrauschen ---
      if applyRadonNoise
          R = R + noiseLevelRadon * randn(size(R));
      end

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
    if exist('t','var') && isvalid(t), stop(t); delete(t); end
    delete(hFig);
  end

end
