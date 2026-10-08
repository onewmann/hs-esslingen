%% Chapter3_4_commented.m
% MATLAB-Skript zu Kapitel 3.4: Selbstbeurteilung in der Einzelbild-Analyse
clear; close all; clc;

%% 1. Parameter
N = 512;                 % N        : Bildgröße (N x N)             
A = 0.3;                 % A        : Texturparameter, steuert Spektralbreite der generierten Textur
v = [153, 53];          % v        : wahrer Geschwindigkeitsvektor [vx, vy] in Pixel pro Frame 
L = 80;                  % L        : Anzahl halb-breiten der Sequenz; numFrames = 2*L+1
dt = 1;                  % dt       : Zeit zwischen Frames (hier 1 frame -> dt = 1)
numFrames = 2 * L + 1;   % Gesamtzahl der Frames (hier 161)

%% 2. Statisches Texturbild & Bildfolge
% simulateTexture: erzeugt ein zufallsphasen-basiertes Texturbild im
% Frequenzraum mit radialer Amplitudenhülle. Ergebnis normiert auf [-1,1].
b = simulateTexture(N, A);                         % statisches Texturbild

% generateImageSequence: erzeugt eine Folge von Frames durch zyklische
% Verschiebung (circshift) des Basisbildes b entsprechend v und numFrames.
sequence = generateImageSequence(b, v, numFrames); % N x N x numFrames

% bp: mittleres Bild (arithmetisches Mittel über die Frames) -> Bewegungsunschärfe
bp = mean(sequence, 3);

% Anzeige des gemittelten Bildes
figure('Name','Summiertes Bild');
imshow(mat2gray(bp));  % mat2gray skaliert auf [0,1] für Anzeige
title('Summiertes Bild bp(x;v)');

%% 3. 2D-FFT & logarithmierte Darstellung des Spektrums
% Ziel: gerichtete Strukturen im Frequenzbereich sichtbar machen
B    = fftshift(fft2(bp));         % 2D-FFT, fftshift zentriert die Nullfrequenz
logB = log(1 + abs(B));            % logarithmierte Betragsdarstellung

figure('Name','Log-Spektrum');
imshow(mat2gray(logB));
title('Logarithmiertes Betragsspektrum von bp(x;v)');

%% 4. Radon-Transformation des log-Spektrums
% Radon transformiert das 2D-Bild entlang Geraden verschiedener Winkel theta.
% Hohe Werte in R bei einem Winkel zeigen dominante Linienorientierung.
theta = 0:1:359;
[R, xp] = radon(logB, theta);      % R: Projektionen, xp: Radialkoordinaten

figure('Name','Radon des Log-Spektrums');
imagesc(theta, xp, R);             % Winkel entlang x-Achse, r entlang y-Achse
xlabel('Winkel °'); ylabel('r');
title('Radon-Transformation des log-Spektrums');
colorbar;

%% 5. Berechnung des summierten Gradientenbetrags WR
% Für jeden Winkel wird die Projektion entlang r genommen, ihr Gradientenbetrag
% summiert und als "Linienstärke" WR(θ) gespeichert.
numAngles = length(theta);
WR = zeros(1, numAngles);
for i = 1:numAngles
    proj     = R(:, i);            % Projektion bei Winkel theta(i)
    gradProj = abs(diff(proj));    % absoluter Gradient entlang r
    WR(i)    = sum(gradProj);      % summierter Betrag -> Maß für Linienstärke
end

%% 6. Bestimme Winkel & Gütemaß (GEBA)
% WR_max: maximaler Summierter Gradient; bestIdx: Winkelindex; GEBA: Qualitätsmaß
[WR_max, bestIdx] = max(WR);
v_angle_radon       = theta(bestIdx);       % Winkel mit maximaler Linienstärke
GEBA                = WR_max / mean(WR);    % GEBA = max(WR) / mean(WR)

% Darstellung der Radon-Projektion am besten Winkel
selectedProj = R(:, bestIdx);
figure('Name','Beste Radon-Projektion');
plot(selectedProj, 'LineWidth', 2);
xlabel('r index'); ylabel('Projektion');
title(sprintf('Radon-Projektion bei %.1f°', v_angle_radon));
grid on;

% Peaks in der Projektion finden, um periodische Struktur (Spacing) zu messen
% Zuerst glätten wir leicht, dann findpeaks mit sinnvollen Parameterwerten
sel = movmedian(selectedProj, 5);  % Median-Glättung, reduziert Ausreißer
[pks, locs] = findpeaks(sel, 'MinPeakProminence', 0.02, 'MinPeakDistance', 3);
if numel(locs) > 1
    spacings    = diff(locs);             % Abstand zwischen Peak-Indizes
    spacing_avg = mean(spacings);         % mittlerer Abstand -> heuristischer Maßstab
else
    spacing_avg = NaN;                    % keine ausreichenden Peaks gefunden
end

% Heuristische Radon-basierte Magnitude (nur indikativ)
v_mag_est_radon = spacing_avg;

%% 7. Robuste Verschiebungsschätzung per Phase-Correlation
% Phase-Correlation liefert eine robuste (auch subpixelfähige) Verschiebung
% zwischen zwei Bilder I1 und I2. Wir verwenden das mittlere Frame und das
% darauf folgende Frame (mittleres Frame = Mitte der Sequenz).
mid = ceil((numFrames+1)/2);
I1 = sequence(:,:,mid);
I2 = sequence(:,:,mid+1);   % ein Frame später

% estimateShiftPhaseCorr gibt [dx, dy] zurück: Verschiebung von I1 zu I2
shift = estimateShiftPhaseCorr(I1, I2);   % dx positiv -> rechts, dy positiv -> down
v_est_vec = shift / dt;                   % Pixel pro Frame
v_mag_est_pc = norm(v_est_vec);           % Betrag (Pixel/frame)
% Winkelbestimmung: atan2d verwendet y nach oben, Bildkoordinaten haben y nach unten,
% daher invertieren wir das y-Zeichen für die Anzeige
v_angle_est_pc = atan2d(-v_est_vec(2), v_est_vec(1));

%% 8. Vergleich, Fallback-Logik und Ausgabe
% Bevorzugte Methode: Phase-Correlation (robust). Falls ungültig, Radon-Fallback.
if all(isfinite(v_est_vec)) && any(abs(v_est_vec) > 0)
    v_mag_est = v_mag_est_pc;
    v_angle_est = v_angle_est_pc;
    method_used = 'Phase Correlation';
else
    v_mag_est = v_mag_est_radon;
    v_angle_est = v_angle_radon;
    method_used = 'Radon spacing fallback';
end

v_true_mag = norm(v);    % wahrer Betrag des vorgegebenen Vektors

% Ausgabe in Konsole
fprintf('Methode: %s\n', method_used);
fprintf('Wahrer Geschwindigkeitsbetrag (Pixel/frame): %.2f\n', v_true_mag);
fprintf('Geschätzter Geschwindigkeitsbetrag (Pixel/frame): %.2f\n', v_mag_est);
fprintf('Geschätzter Geschwindigkeitswinkel: %.2f°\n', v_angle_est);
fprintf('Selbstbeurteilungsmaß GEBA: %.3f\n', GEBA);

% Plot WR-Profil (Linienstärke vs Winkel)
figure('Name','WR Profil');
plot(theta, WR, 'b-', 'LineWidth', 2);
xlabel('Winkel °');
ylabel('Summierter Gradient, WR');
title(sprintf('Welligkeitsprofil WR (GEBA = %.3f)', GEBA));
grid on;

%% 9. Zusätzliche Visualisierung: Frames mit Richtungspfeil
% Zeige mittleres Frame und nächstes Frame; zeichne geschätzten Richtungsvektor ein
figure('Name','Beispiel-Frames');
subplot(1,2,1)
imagesc(sequence(:,:,mid))
colormap('gray')
axis image off
hold on
sz = size(sequence(:,:,1));
center = [sz(2)/2, sz(1)/2];  % [x, y] Mittelpunkt
plot(center(1), center(2), 'ro', 'MarkerSize', 8, 'LineWidth', 1.2)
len = 40;                      % Längenfaktor für Visualisierung
dx = len * cosd(v_angle_est);
dy = -len * sind(v_angle_est);  % invertiert, weil y-Bildkoordinaten nach unten zeigen
quiver(center(1), center(2), dx, dy, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2)
title('Mittleres Frame mit geschätzter Bewegungsrichtung')
hold off

subplot(1,2,2)
imagesc(sequence(:,:,mid+1))
colormap('gray')
axis image off
title('Darauffolgendes Frame')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hilfsfunktionen mit Kommentaren
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = simulateTexture(N, A)
% simulateTexture Erzeugt ein statisches, isotropes Texturbild b(x).
% - Erzeugt ein Frequenzgitter (u,v) und berechnet radialen Abstand rho.
% - Definiert eine radiale Amplitudenhülle amp (z.B. Gauß), moduliert mit
%   einer Zufallsphase für natürliche Texturen.
% - Führt inverse FFT durch und normiert das Ergebnis auf [-1,1].
    if nargin < 2, A = 0.3; end

    % Frequenzgitter (gleichmäßig, zentriert)
    [u, v] = meshgrid( (-floor(N/2)):(ceil(N/2)-1) );
    rho = sqrt(u.^2 + v.^2);
    rho = rho / max(rho(:));        % Normierung auf [0,1]

    % Radiale Spektralform: Gaußförmige Hüllkurve mit Breite sigma = A
    sigma = A;
    amp = exp(- (rho.^2) / (2 * sigma^2));
    amp(1,1) = 0;                    % optional: DC-Komponente reduzieren

    % Zufallsphase für natürliche Textur (reproduzierbar durch rng)
    rng(0);                          % konstante Seed für Reproduzierbarkeit
    phase = exp(1i * 2 * pi * rand(N));

    % Kombiniertes Spektrum und Rücktransformation in den Ortsraum
    S = fftshift(amp) .* phase;      % amp zentriert; multiplikation mit Phase
    b_complex = ifft2(ifftshift(S)) * N^2; % Normierungsfaktor für ifft2/fftshift
    b = real(b_complex);

    % Normalisieren: Mittelwert entfernen, auf maximale Amplitude skalieren
    b = b - mean(b(:));
    b = b / max(abs(b(:)));
end

function seq = generateImageSequence(b, v, numFrames)
% generateImageSequence Erzeugt eine Bildfolge durch Translation von b
% - b: N x N Basisbild
% - v: [vx, vy] Verschiebung in Pixel pro Frame
% - numFrames: Anzahl Frames (ungerade empfohlen, Mitte = Referenz)
%
% Randbehandlung: circshift -> zyklische (wrap-around) Verschiebung.
    if nargin < 3, numFrames = 7; end
    [N, M] = size(b);
    seq = zeros(N, M, numFrames);

    % Mitte der Sequenz so wählen, dass das mittlere Frame Verschiebung 0 hat
    mid = ceil((numFrames+1)/2);

    for k = 1:numFrames
        shift_frames = (k - mid);                    % relative Frame-Index zur Mitte
        shift_x = round( shift_frames * v(1) );      % Verschiebung in x-Richtung (Spalten)
        shift_y = round( shift_frames * v(2) );      % Verschiebung in y-Richtung (Zeilen)

        % circshift interpretiert [shift_y, shift_x] als [rows, cols]
        seq(:,:,k) = circshift(b, [shift_y, shift_x]);
    end
end

function shift = estimateShiftPhaseCorr(I1, I2)
% estimateShiftPhaseCorr Subpixelverschiebung zwischen I1 und I2 per Phase Correlation
% Rückgabe: shift = [dx, dy], Verschiebung von I1 nach I2 (Pixel).
% - Multiplikation mit einem Hanning-Fenster reduziert Kantenartefakte.
% - Normierte Kreuzspektralterm sorgt für Peak-sharpness.
% - Nach Peak-Lokalisierung optional 3x3-Quadratic-Fit für Subpixelkorrektur.
    I1 = double(I1); I2 = double(I2);

    % Hanning-Fenster zum Unterdrücken von Kantenartefakten
    win = hann(size(I1,1)) * hann(size(I1,2))';
    F1 = fft2(I1 .* win);
    F2 = fft2(I2 .* win);

    % Kreuzspektralprodukt, normalisiert auf Einheitsamplitude
    R = F1 .* conj(F2);
    R = R ./ (abs(R) + eps);    % Normalisierung, eps verhindert Division durch 0

    % Inverse FFT ergibt Korrelationsoberfläche; fftshift zentriert Peak
    r = ifft2(R);
    r = fftshift(real(r));

    % Grobe Peak-Lokalisierung (Integer-Pixel)
    [~, idx] = max(r(:));
    [py, px] = ind2sub(size(r), idx);

    % Abstand zum Zentrum -> Verschiebung in Pixel (positive dy -> down)
    center = floor(size(r)/2) + 1;
    dy = py - center(1);
    dx = px - center(2);

    % Subpixelkorrektur: 3x3-Quadratic-Fit an Log-Amplitude des Peak-Patches
    patchRadius = 1;
    yRange = max(py-patchRadius,1):min(py+patchRadius,size(r,1));
    xRange = max(px-patchRadius,1):min(px+patchRadius,size(r,2));
    P = r(yRange, xRange);
    if all(size(P) == [3,3])
        Z = log(abs(P(:)) + eps);   % logarithmieren für stabilere Quadratik-Anpassung
        [Xg, Yg] = meshgrid(-1:1, -1:1);
        % Designmatrix für Quadrik: ax^2 + by^2 + dxy + ex + fy + c
        A = [Xg(:).^2, Yg(:).^2, Xg(:).*Yg(:), Xg(:), Yg(:), ones(9,1)];
        coeff = A \ Z;
        a = coeff(1); b = coeff(2); d = coeff(3);
        e = coeff(4); f = coeff(5);
        % Subpixelmaximum approximativ bestimmen (vereinfachte Herleitung)
        % Avoid division by zero checks
        if abs(2*a) > eps
            dx_sub = -e / (2*a);
        else
            dx_sub = 0;
        end
        if abs(2*b) > eps
            dy_sub = -f / (2*b);
        else
            dy_sub = 0;
        end
        if isfinite(dx_sub), dx = dx + dx_sub; end
        if isfinite(dy_sub), dy = dy + dy_sub; end
    end

    % Rückgabe als [dx, dy] (dx -> Spalten/Richtung rechts, dy -> Zeilen/Richtung down)
    shift = [dx, dy];
end
