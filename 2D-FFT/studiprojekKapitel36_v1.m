%% Chapter3_4_complete_fixed.m
% MATLAB-Skript zu Kapitel 3.4: Selbstbeurteilung in der Einzelbild-Analyse
% Vollständige, robuste Version mit Fehlerkorrekturen und Fallback-Logik
clear; close all; clc;

%% 1. Parameter
N           = 512;             % Bildgröße: N x N Pixel
A           = 4;             % Texturparameter
v           = [16, 14];          % Geschwindigkeitsvektor [vx, vy] in Pixel/frame
v_true_mag  = norm(v);         % Wahrer Betrag (Pixel/frame)
L           = 50;               % Summierung von 2L+1 = 7 Frames
dt          = 1;               % Zeitabstand (Referenz: 1 frame)
numFrames   = 2*L + 1;         % Anzahl der Einzelbilder
theta       = 0:1:179;         % Winkel für Radon-Transformation (°)
image_file  = '';              % leer = synthetische Bilder; sonst Pfad zur .mat Datei

%% 2. Statisches Texturbild & Bildfolge
if isempty(image_file)
    fprintf('Verwende synthetische Bilder.\n');
    b        = simulateTexture(N, A);
    sequence = generateImageSequence(b, v, numFrames);
else
    fprintf('Lade %s.\n', image_file);
    s = load(image_file);
    fn = fieldnames(s);
    recording = double(s.(fn{1}));
    sz = size(recording);
    new_img_sz = min(sz(1), sz(2));

    % sichere Auswahl der Frame-Indizes (alle 2. Frames beginnend bei 1)
    totalFrames = size(recording, 4);
    desiredIndices = 1:2:(2*numFrames-1);
    desiredIndices = desiredIndices(desiredIndices <= totalFrames); % beschränken

    if isempty(desiredIndices)
        error('Keine geeigneten Frames in der geladenen Datei vorhanden.');
    end

    recording = recording(1:new_img_sz, 1:new_img_sz, :, desiredIndices);
    sequence = imresize(recording, [N, N]);

    for i=1:size(sequence,3)
        img_tmp = double(sequence(:,:,i));
        img_tmp = img_tmp - mean(img_tmp(:));
        sigma = std(img_tmp(:));
        sequence(:,:,i) = img_tmp / (sigma + eps);
        sequence(:,:,i) = exponentialFilter(A, sequence(:,:,i));
    end
    sequence = squeeze(sequence);
end

assert((ndims(sequence)==3) && ...
    (size(sequence, 1)==N) && (size(sequence, 2)==N) && ...
    (size(sequence, 3)==numFrames), ...
    'falsche Input Dimension');

bp = mean(sequence, 3);  % Mittelung erzeugt Bewegungsunschärfe

figure;
imshow(mat2gray(bp));
title('Summiertes Bild bp(x;v)');

%% 3. 2D-FFT & logarithmierte Darstellung des Spektrums
B    = fftshift(fft2(bp));
logB = log(1 + abs(B));

figure;
imshow(mat2gray(logB));
title('Logarithmiertes Betragsspektrum von bp(x;v)');

%% 4. Radon-Transformation des log-Spektrums
[R, xp] = radon(logB, theta);

figure;
imagesc(theta, xp, R);
xlabel('Winkel (°)');
ylabel('r');
title('Radon-Transformation des log-Spektrums');
colorbar;

%% 5. Berechnung des summierten Gradientenbetrags WR(theta)
numAngles = numel(theta);
WR        = zeros(1, numAngles);
for i = 1:numAngles
    proj     = R(:, i);
    gradProj = abs(diff(proj));
    WR(i)    = sum(gradProj);
end

%% 6. Bestimme dominanten Winkel & Gütemaß GEBA
[WR_max, bestIdx] = max(WR);
v_angle_est       = theta(bestIdx);
GEBA              = WR_max / mean(WR + eps);

%% 7. Schätzung des Geschwindigkeitsbetrags
selectedProj = R(:, bestIdx);

% Plot der Projektion am besten Winkel
figure;
plot(xp, selectedProj, 'LineWidth', 2);
xlabel('r (Spektrum-Einheit)');
ylabel('Projektion');
title(sprintf('Radon-Projektion bei %.1f°', v_angle_est));

% Peaks suchen und mittleren Abstand bestimmen (robust)
minProm = 0.05 * (max(selectedProj)-min(selectedProj));
minProm = max(minProm, 1e-6); % Sicherheitsuntergrenze
[pks, locs] = findpeaks(selectedProj, 'MinPeakProminence', minProm, 'MinPeakDistance', 2);

if numel(locs) > 1
    % Wenn xp existiert, nutze reale r-Positionen
    if exist('xp','var') && numel(xp) >= max(locs)
        r_positions = xp(locs);
        spacings = diff(r_positions);
    else
        spacings = diff(locs); % fallback: Index-Abstände
    end
    spacing_avg = mean(spacings);
    % Skalierungsfaktor k: hier heuristisch 1; anpassbar durch Kalibrierung
    k = 1;
    v_mag_est = k * spacing_avg;
else
    spacing_avg = NaN;
    v_mag_est = NaN;
    warning('Nicht genügend Peaks gefunden: v_mag_est wird auf NaN gesetzt.');
end

%% 8. Ausgabe der Ergebnisse
if ~exist('v_true_mag','var')
    v_true_mag = norm(v);
end

fprintf('Wahrer Geschwindigkeitsbetrag   : %.2f Pixel/frame\n', v_true_mag);
if ~isnan(v_mag_est)
    fprintf('Geschätzter Geschwindigkeitsbetrag: %.2f Pixel/frame\n', v_mag_est);
else
    fprintf('Geschätzter Geschwindigkeitsbetrag: NaN (nicht genügend Peaks)\n');
end
fprintf('Geschätzter Geschwindigkeitswinkel: %.2f°\n',         v_angle_est);
fprintf('Selbstbeurteilungsmaß GEBA        : %.3f\n',         GEBA);

figure;
plot(theta, WR, 'b-', 'LineWidth', 2);
xlabel('Winkel (°)');
ylabel('Summierter Gradient, WR');
title(sprintf('Welligkeitsprofil WR (GEBA = %.3f)', GEBA));
grid on;

%% 9. Visualisierung von zwei Frames & Bewegungsrichtung
figure;
sz = size(sequence(:,:,1));
center_x = sz(2)/2;
center_y = sz(1)/2;

subplot(1,2,1);
imagesc(sequence(:,:,1));
colormap('gray');
axis image off;
hold on;
plot(center_x, center_y, 'ro', 'MarkerSize', 8);

if ~isnan(v_mag_est)
    angle_rad = deg2rad(v_angle_est);
    dx = v_mag_est * cos(angle_rad);
    dy = v_mag_est * sin(angle_rad);
    quiver(center_x, center_y, dx, -dy, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', 2);
end
title('Frame k0');
hold off;

subplot(1,2,2);
imagesc(sequence(:,:,2));
colormap('gray');
axis image off;
title('Frame k0+1');

figure;
imagesc(abs(sequence(:,:,2)-sequence(:,:,1)));
colormap('hsv');
axis image;
title('Differenzbild');
colorbar;

%% --- Hilfsfunktionen -----------------------------------------------------
function b = simulateTexture(N, A)
% Erzeugt ein normalisiertes, texturiertes Bild mit spektraler Charakteristik
    rng(0);
    white = randn(N);
    % Gauß-Filterbreite proportional zu A (min 1)
    sigma = max(1, round(A * 10));
    H = fspecial('gaussian', 6*sigma+1, sigma);
    b = conv2(white, H, 'same');
    b = b - mean(b(:));
    b = b / (std(b(:)) + eps);
end

function seq = generateImageSequence(b, v, numFrames)
% Verschiebt das statische Bild b subpixelgenau über numFrames (zentriert bei Mittelframe)
    [N, ~] = size(b);
    mid = ceil(numFrames/2);
    seq = zeros(N, N, numFrames);
    % Erzeuge Gitter für Fouriershift fallback
    [X, Y] = meshgrid(0:N-1, 0:N-1);
    for k = 1:numFrames
        dt = (k - mid);
        shiftx = v(1) * dt;
        shifty = v(2) * dt;
        % Nutze imtranslate falls vorhanden
        try
            seq(:,:,k) = imtranslate(b, [shiftx, shifty], 'OutputView', 'same', 'FillValues', 0);
        catch
            % Fourier-basiertes subpixel shift (Fallback)
            F = fft2(b);
            % Phasenfaktor: exp(-i*2*pi*(u*dx/N + v*dy/N))
            [u, vgrid] = meshgrid(0:N-1, 0:N-1);
            phase = exp(-1i*2*pi*( (u * shiftx)/N + (vgrid * shifty)/N ));
            seq(:,:,k) = real(ifft2(F .* phase));
        end
    end
    % normalize each frame
    for k = 1:numFrames
        f = seq(:,:,k);
        f = f - mean(f(:));
        seq(:,:,k) = f / (std(f(:)) + eps);
    end
end

function out = exponentialFilter(A, img)
% Einfacher Exponentialfilter zur Betonung niederfrequenter Anteile
    sz = size(img,1);
    [u,v] = meshgrid( (-sz/2):(sz/2-1) );
    r = sqrt(u.^2 + v.^2);
    H = exp(-A * (r / max(r(:))));
    IMG = fftshift(fft2(img));
    out = real(ifft2(ifftshift(IMG .* H)));
end
