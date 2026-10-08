%% Chapter3_4_complete.m
% MATLAB-Skript zu Kapitel 3.4: Selbstbeurteilung in der Einzelbild-Analyse
%
% Ziel: Bestimmung eines Gütemaßes (GEBA) und Schätzung des
%       Geschwindigkeitsbetrags und -winkels bei der Einzelbild-Analyse.

clear; close all; clc;

%% 1. Parameter
N           = 256;             % Bildgröße: N x N Pixel
A           = 0.3;             % Texturparameter
v           = [150, -150];       % Geschwindigkeitsvektor [vx, vy] in Pixel/frame
v_true_mag  = norm(v);         % Wahrer Betrag (Pixel/frame)
L           = 2;               % Summierung von 2L+1 = 7 Frames
dt          = 1;               % Zeitabstand (Referenz: 1 frame)
numFrames   = 2*L + 1;         % Anzahl der Einzelbilder
theta       = 0:1:179;         % Winkel für Radon-Transformation
image_file  = 'data/vid4.mat'; % leer lassen für synthethische Bilder

%% 2. Statisches Texturbild & Bildfolge
if isempty(image_file)
    fprintf('Verwende synthetische Bilder.\n')
    b        = simulateTexture(N, A);
    sequence = generateImageSequence(b, v, numFrames);
else
    fprintf('Lade %s.\n', image_file)
    load(image_file)
    recording = double(recording);
    sz = size(recording);
    new_img_sz = min(sz(1), sz(2));
    recording = recording(1:new_img_sz, 1:new_img_sz, :, 1:2:2*numFrames-1);

    sequence = imresize(recording, [N, N]);

    for i=1:size(sequence,3)
        sequence(:,:,i) = sequence(:,:,i) - mean2(sequence(:,:,i));
        img_tmp = sequence(:,:,i);
        std = sqrt(var(img_tmp(:)));
        sequence(:,:,i) = sequence(:,:,i)/(std+eps);
        sequence(:,:,i) = exponentialFilter(A, sequence(:,:,i));
    end
    sequence = squeeze(sequence);
end
assert((ndims(sequence)==3) && ...
    (size(sequence, 1)==N) && (size(sequence, 2)==N) && ...
    (size(sequence, 3)==numFrames), ...
    'falsche Input Dimension')
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

%% 5. Berechnung des summierten Gradientenbetrags WR(θ)
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
GEBA              = WR_max / mean(WR);

%% 7. Schätzung des Geschwindigkeitsbetrags
selectedProj = R(:, bestIdx);

% Plot der Projektion am besten Winkel
figure;
plot(selectedProj, 'LineWidth', 2);
xlabel('r index');
ylabel('Projektion');
title(sprintf('Radon-Projektion bei %.1f°', v_angle_est));

% Peaks suchen und mittleren Abstand bestimmen
[pks, locs] = findpeaks(selectedProj, 'MinPeakProminence', 0.05);
if numel(locs) > 1
    spacings    = diff(locs);
    spacing_avg = mean(spacings);
else
    spacing_avg = NaN;
end

% Direkte Proportionalität: Abstand = Geschwindigkeit (in Einheiten/frame)
v_mag_est = spacing_avg;

%% 8. Ausgabe der Ergebnisse
fprintf('Wahrer Geschwindigkeitsbetrag   : %.2f Pixel/frame\n', v_true_mag);
fprintf('Geschätzter Geschwindigkeitsbetrag: %.2f Pixel/frame\n', v_mag_est);
fprintf('Geschätzter Geschwindigkeitswinkel: %.2f°\n',         v_angle_est);
fprintf('Selbstbeurteilungsmaß GEBA        : %.3f\n',         GEBA);

%% 9. Visualisierung von zwei Frames & Bewegungsrichtung
figure;
sz = size(sequence(:,:,1));
dx = v_mag_est * sin(v_angle_est*pi/180);
dy = v_mag_est * cos(v_angle_est*pi/180);

subplot(1,2,1);
imagesc(sequence(:,:,1));
colormap('gray');
axis image off;
hold on;
plot(sz(2)/2, sz(1)/2, 'ro', 'MarkerSize', 8);
plot([sz(2)/2, sz(2)/2+dx], [sz(1)/2, sz(1)/2-dy], 'r-', 'LineWidth', 2);
title('Frame k0');

subplot(1,2,2);
imagesc(sequence(:,:,2));
colormap('gray');
axis image off;
title('Frame k0+1');

figure
imagesc(abs(sequence(:,:,2)-sequence(:,:,1)))
colormap('hsv')
grid on

