function live_radon_basler_gui_save_seconds()
% LIVE_RADON_BASLER_GUI_SAVE_SECONDS
% Echtzeit-Radon-Velocity mit Basler-Kamera (gentl, Mono8) in einer GUI.
% Misst Winkel und Geschwindigkeit, speichert Zeit (s seit Start).
% Nach Stop zwei Diagramme mit:
%   • Winkel β über Zeit [0…179°]
%   • Geschwindigkeit v über Zeit
% Die Live-WR-Achse bleibt bis Start ausgeblendet.
  %% 1. Parameter
  L       = 4;              % Halbfenster für Ringpuffer-Mittelung (Summe aus 2L+1 Frames)
  bufLen  = 2*L + 1;        % Pufferlänge, ungerade für zentrierte Mittelung
  theta   = 0:179;          % Radon-Winkel in Grad (1°-Raster, volle Halbkreisabdeckung)
  N       = 256;            % Arbeitsauflösung (Bild wird auf NxN rescaled)
  padDim  = 2 * N;          % Zero-Padding-Dimension für FFT (Hann-Fenster auf NxN, danach Pad auf padDim)
  clipLim = [0.01 0.99];    % Kontrast-Stretch-Limits (egalisiert globale Helligkeit)
  %% 2. Basler-Kamera initialisieren
  hw = imaqhwinfo;                                                      % Systemweite IMAQ-Info
  assert(ismember('gentl',hw.InstalledAdaptors), 'gentl-Adapter fehlt.'); % Sicherstellen, dass GenTL verfügbar ist
  camInfo = imaqhwinfo('gentl');                                        % Adapter-spezifische Info
  devID   = camInfo.DeviceIDs{1};                                       % Erstes Gerät (vereinfachend)
  fmts    = imaqhwinfo('gentl',devID).SupportedFormats;                 % Liste unterstützter Formate
  idxMono = find(contains(fmts,'Mono8'),1,'first');                      % Erstes Mono8-Format wählen
  fmt     = fmts{idxMono};                                              % Format-String
  vid = videoinput('gentl',devID,fmt);                                  % Video-Input-Objekt erstellen
  
  % ROI sofort auf volle Auflösung setzen (vermeidet implizite Crops)
  res = vid.VideoResolution;             
  vid.ROIPosition = [0 0 res(1) res(2)]; 
  % Trigger auf manuell, 1 Frame pro Trigger, unendliche Wiederholung
  triggerconfig(vid,'manual');
  vid.FramesPerTrigger   = 1;
  vid.TriggerRepeat      = Inf;
  vid.ReturnedColorspace = 'grayscale';  % Farbraum als Graustufen
  start(vid); flushdata(vid);            % Stream starten und Puffer leeren
  %% 3. Ring-Puffer & Mess-Arrays
  frameBuf   = zeros(N, N, bufLen);  % 3D-Puffer für gleitende Mittelung
  bufIdx     = 1;                    % Schreibindex (rotierend)
  count      = 0;                    % Anzahl aktuell gefüllter Frames
  timeSec    = [];                   % Zeitstempel (s seit Start)
  beta_vals  = [];                   % Winkel-Messungen β (Grad)
  velocities = [];                   % Geschwindigkeits-Messungen v (px/frame)
  t0 = tic;                          % Start-Timer für Zeitmessung
  %% 4. GUI aufsetzen
  hFig = figure('Name','Live Radon (Basler)','NumberTitle','off',...
                'CloseRequestFcn',@onClose); % Cleanup beim Schließen
  % Live-Einzelbild-Ansicht (nach Preprocessing)
  subplot(2,2,1);
    hLive = imshow(zeros(N,N));
    title('Einzelbild');
  % Gepuffertes Summen-/Mittelbild (Basis für Radon)
  subplot(2,2,2);
    hSum = imshow(zeros(N,N));
    title('Summenbild');
  % WR-Profil (Gewichtete Radon-Differenzen über Winkel)
  hAx = subplot(2,2,[3 4]);
    hPlot = plot(hAx, theta, zeros(size(theta)), 'LineWidth',1.5);
    xlabel(hAx,'\theta (°)'); ylabel(hAx,'WR');
    title(hAx,'WR-Profil');
    grid(hAx,'on');
    % Ausblenden bis Start (UI bleibt clean bis Aufnahme beginnt)
    set(hAx,'Visible','off');
    set(hPlot,'Visible','off');
  % Start/Stop Toggle-Button
  uicontrol('Style','togglebutton','String','Start',...
            'Position',[10 10 60 30],'Parent',hFig,...
            'Callback',@onToggle);
  %% 5. Timer konfigurieren
  % Periodische Abtastung (alle 0.1 s) für Trigger + Verarbeitung
  t = timer('ExecutionMode','fixedRate','Period',0.1,...
            'TimerFcn',@onTimer,'ErrorFcn',@(~,e) disp(e.Data.message));
  %% Callback Start/Stop
  function onToggle(src,~)
    if get(src,'Value')
      set(src,'String','Stop');
      % Achse anzeigen, sobald Messung startet
      set(hAx,'Visible','on');
      set(hPlot,'Visible','on');
      start(t);                   % Timer starten
    else
      set(src,'String','Start');
      stop(t);                    % Timer stoppen
      plotResults();              % Nachlauf: Zeitreihen darstellen
    end
  end
  %% TimerFcn: Live-Loop
  % Holt jeweils 1 Frame, verarbeitet, puffert und führt Radon-Analyse aus, wenn Puffer voll.
  function onTimer(~,~)
    trigger(vid);                 % Manuell einen Frame anfordern
    frameRaw = getdata(vid,1);    % 1 Frame lesen
    % Resize & Kontrastnormalisierung
    I = im2double(frameRaw);                          % In double [0..1] konvertieren
    I = imresize(I,[N N]);                            % Auf Arbeitsauflösung bringen
    I = imadjust(I, stretchlim(I,clipLim), []);       % Kontrast spreizen anhand Percentiles
    set(hLive,'CData',I);                             % Live-Ansicht aktualisieren
    % Ringpuffer füllen (rotierender Index)
    frameBuf(:,:,bufIdx) = I;
    bufIdx = bufIdx+1;
    if bufIdx>bufLen, bufIdx=1; end
    count = min(count+1, bufLen);                     % bis Puffer voll
    % Wenn Puffer voll, Radon-Analyse auf dem gemittelten Bild
    if count==bufLen
      bp = mean(frameBuf,3);                          % Gleitende Mittelung über Puffer
      set(hSum,'CData',bp);                           % Summen-/Mittelbild zeigen
      % Hintergrund entfernen (Morph. Öffnung) und leicht denoisen
      bg    = imopen(bp,strel('disk',round(N/20)));   % Großskalige Strukturen als Hintergrund
      bp_dc = max(bp-bg,0);                           % DC-Anteil entfernen, negativ clampen
      bp_dn = wiener2(bp_dc,[5 5]);                   % Rauschminderung (adaptiv)
      % Hann-Fenster zur Randglättung + Zero-Padding für FFT
      Hwin = hann(N); H = Hwin*Hwin.';
      Ipad = padarray(bp_dn.*H,[(padDim-N)/2 (padDim-N)/2],0,'both');
      % FFT und Log-Magnitude (betont Frequenzstrukturen)
      B      = fftshift(fft2(Ipad));
      logB   = log(1+abs(B));
      % Radon-Transformation über Winkelraster
      [R,xp] = radon(logB, theta);
      % WR-Metrik: Summe der Absolutdifferenzen über Projektionen => Richtungsdominanz
      WR       = sum(abs(diff(R,1,1)),1);
      [~,idx]  = max(WR);           % stärkste Richtung indexieren
      beta     = theta(idx);        % Winkel β in Grad
      % Projektion in der dominanten Richtung
      proj     = R(:,idx);
      mask     = xp>=0;             % nur positive Radien (eine Halbachse) nutzen
      xp_p     = xp(mask);
      proj_p   = proj(mask);
      % Peak-Suche (Abstände korrelieren mit Frequenz/Bewegung)
      [~,locs] = findpeaks(proj_p,xp_p,'MinPeakProminence',0.05,'MinPeakDistance',3);
      % Subpixel-Verfeinerung der Peak-Position durch quadratische Interpolation
      locs_sub = locs;
      for k=1:numel(locs)
        [~,ix] = min(abs(xp_p-locs(k)));
        if ix>1 && ix<length(proj_p)
          y1=proj_p(ix-1); y0=proj_p(ix); y2=proj_p(ix+1);
          % Parabolische Scheitelpunkt-Korrektur: Δx = (y1 - y2) / (2*(y1 - 2*y0 + y2))
          locs_sub(k) = locs(k)+(y1-y2)/(2*(y1-2*y0+y2));
        end
      end
      % Geschwindigkeits-Schätzer: padDim / mittlerer Peakabstand (px/frame)
      % Annahme: regelmäßige Struktur -> Geschwindigkeit ∝ inverse Abstand
      v_mag = padDim / mean(diff(sort(locs_sub)));
      % Live-Profil updaten (WR-Kurve + Titel mit aktuellen Schätzwerten)
      set(hPlot,'YData',WR);
      title(hAx,sprintf('β=%.1f°  v=%.1f px/frame',beta,v_mag));
      drawnow;
      % Zeit & Messwerte speichern (fortlaufende Zeitbasis ab Start)
      elapsed = toc(t0);
      timeSec(end+1)    = elapsed;
      beta_vals(end+1)  = beta;
      velocities(end+1) = v_mag;
    end
  end
  %% Nach Stop: Auswertung plotten
  % Erzeugt zwei Figuren: β(t) und v(t) mit adaptiver X-Tick-Skalierung.
  function plotResults()
    tmax = max(timeSec);                           % Maximale Zeit für Achsenlimit
    % Adaptive X-Ticks (fein für kurz, gröber für lang)
    if tmax<=10
      xt = 0:1:ceil(tmax);
    elseif tmax<=60
      xt = 0:5:ceil(tmax/5)*5;
    else
      xt = 0:10:ceil(tmax/10)*10;
    end
    % Winkel über Zeit
    figure('Name','Winkel über Zeit (s)');
    plot(timeSec,beta_vals,'-o','LineWidth',1.5);
    xlabel('Zeit seit Start (s)');
    ylabel('Winkel β (°)');
    xlim([0 tmax]); ylim([0 179]);
    xticks(xt); set(gca,'YTick',0:30:180);
    grid on;
    % Punkt-Index zur schnellen Zuordnung (optional)
    for i=1:numel(timeSec)
      text(timeSec(i),beta_vals(i),sprintf('%d',i),...
           'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
    % Geschwindigkeit über Zeit
    figure('Name','Geschwindigkeit über Zeit (s)');
    plot(timeSec,velocities,'-o','LineWidth',1.5);
    xlabel('Zeit seit Start (s)');
    ylabel('v (px/frame)');
    xlim([0 tmax]); xticks(xt);
    grid on;
    % Punkt-Index markieren (optional)
    for i=1:numel(timeSec)
      text(timeSec(i),velocities(i),sprintf('%d',i),...
           'VerticalAlignment','bottom','HorizontalAlignment','right');
    end
  end
  %% Cleanup beim Schließen
  % Sorgt für geordnetes Stoppen/Löschen von Timer und Kamera sowie der GUI.
  function onClose(~,~)
    if isvalid(t), stop(t); delete(t); end        % Timer stoppen/entsorgen
    if isvalid(vid), stop(vid); delete(vid); end  % Kamera stoppen/entsorgen
    delete(hFig);                                 % GUI schließen
  end
end
