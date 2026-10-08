% single_frame_velocity_estimation.m
% Vollständiges, robustes Skript zur Einzelbild-Analyse (Kap. 3)
% - lädt data/vid4.mat (oder anderes .mat mit 3D-Aufnahme)
% - berechnet Summenbild (Summation oder Mittelung)
% - berechnet Spektrum (B = FFT), verwendet |B| für Radon
% - berechnet WR (summierter Gradientbetrag) zur Winkelschätzung
% - lokalisiert Peaks in der Radon-Projektion (xp-Einheiten) zur Betragsschätzung
% - bietet Fallback über 1D-FFT der Projektion
% - enthält einfache Kalibrierungsformel (Platzhalter): Benutzer muss fL bzw. Messparameter setzen
% - robustere GEBA-Normalisierung (Median)
%
% Aufrufen: Im selben Verzeichnis wie data/vid4.mat ausführen oder image_file anpassen.
% Hinweis: Kalibrierungs-Konstante C_v (siehe Abschnitt "Kalibrierung") ist erforderlich
%         um das Ergebnis in px/frame oder m/s zu bekommen.

clear; close all; clc;
rng(0); % reproduzierbar

%% --------- Benutzereinstellungen ---------
image_file = 'data/vid4.mat';   % Pfad zur .mat-Datei; leer '' für synthetisch
requestedVarName = '';         % optional: Name der Variable in .mat, leer = automatisch
N = 256;                       % Zielgröße (N x N) für Analyse
useMeanSummation = true;       % true: mean(sequence,3); false: sum(sequence,3)
theta = 0:1:179;               % Winkel für Radon (deg)
peakProminenceRel = 0.05;      % findpeaks MinPeakProminence (relativ)
wrSmoothKernel = 5;            % glättungsfenster für WR (samples)
GEBA_normalizer = 'median';    % 'median' oder 'mean'
% Kalibrierung: Benötigte Messtechnik-Parameter (Benutzer muss prüfen)
fL = [];   % Pulsfrequenz in Hz (falls gepulste Belichtung) oder Bildrate bei Summation (Hz)
C_v = [];  % optionale Kalibrierungskonstante: v [px/frame] = C_v * spacing_pixels
% Wenn C_v leer bleibt, meldet Skript "uncalibrated" und liefert spacing in px.

%% --------- Lade oder generiere Sequenz ---------
if isempty(image_file)
    fprintf('Kein image_file angegeben — Erzeuge synthetische Sequenz.\n');
    % einfache synthetische Textursequenz (fallback)
    base = imresize(randn(512), [N N]);
    % Erzeuge kleine Translationen
    v_test = [10, 5]; % px per frame (synthetisch)
    numFrames = 7;
    sequence = zeros(N, N, numFrames);
    center = ceil(numFrames/2);
    for k = 1:numFrames
        shift = (k - center) * v_test;
        sequence(:,:,k) = imtranslate(base, shift, 'linear', 'FillValues', 0);
    end
else
    % Lade .mat und finde geeignete Variable
    if ~isfile(image_file)
        error('Datei nicht gefunden: %s', image_file);
    end
    S = load(image_file);
    vars = fieldnames(S);
    recording = [];
    if ~isempty(requestedVarName) && isfield(S, requestedVarName)
        recording = double(S.(requestedVarName));
    else
        % Suche erstes plausibles 3D-Array
        for i = 1:numel(vars)
            vtmp = S.(vars{i});
            if isnumeric(vtmp) && ndims(vtmp) >= 3
                recording = double(vtmp);
                fprintf('Verwende Variable "%s" aus %s\n', vars{i}, image_file);
                break;
            end
        end
    end
    if isempty(recording)
        error('Keine geeignete 3D-Variable in %s gefunden.', image_file);
    end

    % Vorbereitung: center-crop, wähle frames (gleichmäßig verteilt), resize auf N x N
    sz = size(recording);
    minSide = min(sz(1), sz(2));
    r0 = floor((sz(1)-minSide)/2)+1; c0 = floor((sz(2)-minSide)/2)+1;
    crop = recording(r0:r0+minSide-1, c0:c0+minSide-1, :);

    % Wähle gleichmäßig numFrames aus verfügbaren Frames (max 2L+1 aus Thesis-Beispiel)
    % Ziel: möglichst viele Frames aber begrenze auf 2*L+1 ~ 7 wenn vorhanden.
    maxDesired = 7;
    avail = size(crop,3);
    numFrames = min(avail, maxDesired);
    idx = round(linspace(1, avail, numFrames));
    seq = crop(:,:,idx);

    % Resize auf N x N und Normierung je Bild (DC entfernen, Normierung)
    sequence = zeros(N, N, numFrames);
    for k = 1:numFrames
        img = im2double(seq(:,:,k));
        img = imresize(img, [N N]);
        img = img - mean(img(:));               % DC entfernen
        imgStd = std(img(:));
        sequence(:,:,k) = img ./ (imgStd + eps); % Normierung
    end
end

% Validierung
assert(ndims(sequence)==3 && size(sequence,1)==N && size(sequence,2)==N, 'sequence dims invalid');

%% --------- Summenbild (by) ---------
if useMeanSummation
    by = mean(sequence, 3);   % skalenstabil
else
    by = sum(sequence, 3);
end
figure('Name','Summed image'); imshow(mat2gray(by)); title('Summed image by(x;v)');

%% --------- Spektrum berechnen (magnitude) ---------
B = fftshift(fft2(by));
absB = abs(B);
logAbsB = log(1 + absB); % nur für Visualisierung optional
figure('Name','Spectrum (magnitude)');
subplot(1,2,1); imshow(mat2gray(absB)); title('|B| (magnitude)');
subplot(1,2,2); imshow(mat2gray(logAbsB)); title('log(1+|B|)');

%% --------- Radon-Transformation auf |B| (nicht auf Log) ---------
% Für radon() benötigen wir ein 2D-Image; radon nimmt x-y Bild -> Projektionsraum
[R, xp] = radon(absB, theta);   % R: rows correspond to xp, cols to theta
figure('Name','Radon of |B|'); imagesc(theta, xp, R); axis xy; colorbar;
xlabel('Angle (deg)'); ylabel('r (freq-index)'); title('Radon(|B|)');

%% --------- WR: summierter Gradientbetrag über Radon-Spalten ---------
numAngles = numel(theta);
WR = zeros(1, numAngles);
for ia = 1:numAngles
    proj = R(:, ia);
    gradProj = abs(diff(proj));
    WR(ia) = sum(gradProj);
end
WRs = smoothdata(WR, 'gaussian', wrSmoothKernel);

% Robust Normalizer
switch GEBA_normalizer
    case 'median'
        normVal = median(WRs) + eps;
    otherwise
        normVal = mean(WRs) + eps;
end
[WR_max, bestIdx] = max(WRs);
beta_est = theta(bestIdx);      % Geschwindigkeitswinkel (bis 180° Ambig.)
GEBA = WR_max / normVal;

fprintf('Estimated angle (beta) : %.2f deg\n', beta_est);
fprintf('GEBA (quality measure) : %.4f\n', GEBA);

% Plot WR
figure('Name','WR vs angle');
plot(theta, WRs, '-b', 'LineWidth', 1.6); hold on;
plot(beta_est, WRs(bestIdx), 'ro', 'MarkerSize',8, 'LineWidth',1.5);
xlabel('Angle (deg)'); ylabel('WR'); title('WR(\phi) with estimated \beta'); grid on;

%% --------- Projektion bei bestem Winkel und Peak-Analyse (xp-Einheiten) ---------
selectedProj = R(:, bestIdx);  % vector with same length as xp
figure('Name','Selected Radon projection');
plot(xp, selectedProj, 'k-', 'LineWidth', 1.4); grid on;
xlabel('r (projection coordinate xp)'); ylabel('Projection amplitude');
title(sprintf('Radon projection at %.2f°', beta_est));

% Peak-Find mit xp-Koordinaten (so sind Abstandseinheiten interpretable)
minProm = peakProminenceRel * range(selectedProj);
[pks, locs_xp] = findpeaks(selectedProj, xp, 'MinPeakProminence', max(minProm, eps));

if numel(locs_xp) >= 2
    spacings_xp = diff(locs_xp);           % Abstand in xp-Einheiten
    spacing_avg_xp = mean(spacings_xp);
    fprintf('Detected %d peaks in projection. Avg spacing (xp units): %.4f\n', numel(locs_xp), spacing_avg_xp);
else
    spacing_avg_xp = NaN;
    fprintf('Weniger als 2 Peaks detektiert in Projektionssignal -> fallback zur FFT-Methode\n');
end

%% --------- Fallback: 1D-FFT der Projektion -> Peak in Frequenzdomäne ---------
% Diese Methode liefert eine periodische Grundfrequenz in Einheiten cycles/xpUnit
proj_fft = abs(fftshift(fft(selectedProj - mean(selectedProj))));
nP = numel(proj_fft);
% Frequenzachse in cycles per xpUnit: xp axis is symmetric; construct relative freq
freqs = linspace(-0.5, 0.5, nP);         % normalized cycles per sample (per xp index)
[~, kmax] = max(proj_fft);
peakFreq = freqs(kmax);                  % cycles per xp-sample
if abs(peakFreq) > 1e-12
    spacing_pixels_est = 1 / abs(peakFreq);  % samples (xp-samples) per cycle
    fprintf('1D-FFT peakFreq (cycles/xp-sample): %.6f -> spacing (xp-samples per cycle): %.4f\n', peakFreq, spacing_pixels_est);
else
    spacing_pixels_est = NaN;
    fprintf('Kein deutlicher Frequenzpeak in 1D-FFT gefunden.\n');
end

% Decide which spacing to use: prefer xp-based peak spacing if available
if ~isnan(spacing_avg_xp)
    spacing_used = spacing_avg_xp;   % xp-units (which correspond to projection coordinate scale)
else
    spacing_used = spacing_pixels_est; % fallback (relative to xp-sample spacing)
end

%% --------- Kalibrierung: spacing -> Geschwindigkeit (px/frame)
% Hinweis: Ohne exakte Messparameter ist das Ergebnis nicht in physikalischen Einheiten kalibriert.
% Benutzer muss C_v definieren, oder die Beziehung aus Kap. 3 (fL, Sampling, Summation) einsetzen.
if isempty(C_v)
    % Versuche heuristische Umrechnung falls fL gegeben:
    if ~isempty(fL) && ~isnan(spacing_used)
        % Beispielannahme (benutzer prüfen!):
        % Für Summation: impulses spacing in frequency <-> spatial period ~ 1/peakFreq
        % Eine mögliche relation (Platzhalter): v [px/frame] = spacing_used * fL
        % Diese Annahme ist projekt-spezifisch und meist falsch ohne Modell.
        v_mag_est = spacing_used * fL; % PLATZHALTER
        calibrated = true;
        fprintf('Kalibrierung mit fL angewendet (PLATZHALTER-Formel): v_est = spacing * fL\n');
    else
        v_mag_est = spacing_used;
        calibrated = false;
        fprintf('Keine Kalibrierung angewendet: gebe spacing in xp-units (oder px) zurück. Setze C_v/fL für kalibrierte Werte.\n');
    end
else
    if ~isnan(spacing_used)
        v_mag_est = C_v * spacing_used;
        calibrated = true;
        fprintf('Kalibrierung über C_v angewendet: v_est = C_v * spacing\n');
    else
        v_mag_est = NaN;
        calibrated = false;
    end
end

%% --------- Ausgabe der Resultate ---------
fprintf('\n--- Result summary ---\n');
fprintf('Estimated angle (beta)           : %.2f deg\n', beta_est);
if calibrated
    fprintf('Estimated speed (calibrated)    : %.4f (units depend on C_v/fL)\n', v_mag_est);
else
    fprintf('Estimated spacing (un-calibrated): %.4f (xp-units or px)\n', v_mag_est);
end
fprintf('Quality measure GEBA            : %.4f\n', GEBA);

%% --------- Plots zur Diagnose ---------
figure('Name','Diagnostics: projection FFT & peaks');
subplot(2,1,1);
plot(freqs, proj_fft, 'LineWidth',1.2); grid on;
xlabel('Normalized frequency (cycles per xp-sample)');
ylabel('|FFT|');
title('1D-FFT magnitude of selected projection');

subplot(2,1,2);
plot(xp, selectedProj, '-k'); hold on;
if ~isempty(locs_xp)
    plot(locs_xp, pks, 'ro', 'MarkerFaceColor','r');
end
xlabel('xp (projection coord)'); ylabel('Projection amplitude');
title('Selected projection with detected peaks');

%% --------- Optional: Visualisierung der Richtung im Bildzentrum ---------
figure('Name','Direction overlay');
imagesc(sequence(:,:,1)); colormap gray; axis image off; hold on;
cx = N/2; cy = N/2;
plot(cx, cy, 'ro', 'MarkerSize',6);
if ~isnan(v_mag_est)
    % angle convention: beta_est from Radon (deg). We draw vector for visualization.
    a = deg2rad(beta_est);
    dx = v_mag_est * cos(a);
    dy = v_mag_est * sin(a);
    quiver(cx, cy, dx, -dy, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2); % -dy weil Bild-y nach unten wächst
end
title('Estimated motion direction (overlay)');

%% --------- Ende Skript ---------
fprintf('Fertig. Prüfen Sie die Diagnostik-Plots und Kalibrierungsparameter (C_v / fL).\n');
