%% masterAnalysis.m
% Großes MATLAB-Programm zur vollumfänglichen Einzelbild-Analyse
% Kapitel 3.2–3.4 aus der Dissertation:
%   • Fluss-Analyse: Extraktion von Richtung und Betrag
%   • Auflösungsreduktion: Pyramid-Down vs. Crop
%   • Selbstbeurteilung: Gütemaß GEBA
%
% Wissenschaftlicher Kontext:
%   Die vorliegende Routine demonstriert die Methode, ein gepulst
%   belichtetes Bild (bp) mittels Fourier- und Radon-Transformation
%   auszuwerten. Der Winkel des Geschwindigkeitsvektors ergibt sich
%   als Maximum eines Welligkeitsprofils WR(θ), der Betrag über den
%   mittleren Abstand der sich wiederholenden Strukturen in der Radon-Domäne.
%   Abschließend wird ein dimensionsloses Messgüte-Kriterium GEBA
%   (Güte der Einzelbild-Analyse) definiert als Verhältnis aus
%   Maximum und Mittelwert von WR(θ).

clear; close all; clc;

%% 1. Parameterdefinition
% Bildgröße, Texturparameter und Bewegungsvektor
params.N           = 256;              % Auflösung: N×N Pixel
params.A           = 0.3;              % Bandbreite des Exponentialfilters
params.v            = [8, 8];          % True velocity [vx, vy] in Pixel/Frame
params.L           = 3;                % Halbe Impulsanzahl → Gesamt 2L+1 Frames
params.dt          = 1;                % Zeitdiskretisierung (Frameabstand)
params.numFrames   = 2*params.L + 1;   % Gesamtzahl der Einzelbilder
params.downFactor  = 0.5;              % Herunterskalierung (Pyramid-Down)
params.cropSize    = 128;              % Seitenlänge des zentralen Ausschnitts

%% 2. Statisches Texturbild & Bildfolge
% 2.1 Simulation eines stochastisch geprägten Texturbilds gemäß 
%     Abschnitt 3.1.1 (WhiteNoise * ExponentialFilter)
b = simulateTexture(params.N, params.A);

function b = simulateTexture(N, A)
    % simulateTexture: Generiert synthetisches Texturbild
    %   Eingabe:
    %     N : Bildkantenlänge in Pixeln
    %     A : Filterbandbreite (Exponentialkernel)
    %   Verfahren:
    %     1. Weißes Rauschen generieren
    %     2. Filterung mit isotropem Exponentialkernel:
    %        h(r) = exp(−r/A) normalisiert auf ∑h = 1
    %
    kernel_radius = ceil(3 * A);
    [X, Y]       = meshgrid(-kernel_radius:kernel_radius);
    h            = exp(-sqrt(X.^2 + Y.^2) / A);
    h            = h / sum(h(:));
    noise        = randn(N, N);
    b            = imfilter(noise, h, 'replicate');
end

figure;
imshow(mat2gray(b));
title('Statisches Texturbild b(x)');

% 2.2 Generierung einer Bildfolge {b_k(x)} mit linearer Verschiebung
%     b_k(x) = b(x − v·(k − k0)), k0 = ceil(numFrames/2)
sequence = generateImageSequence(b, params.v, params.numFrames);

function sequence = generateImageSequence(b, v, numFrames)
    % generateImageSequence: Erstellt verschobene Bildfolge
    %   Eingabe:
    %     b         : Statisches Texturbild (N×N)
    %     v         : [vx, vy] Bewegung pro Frame (Pixel)
    %     numFrames : Anzahl Frames (2L+1)
    %   Ausgabe:
    %     sequence  : 3D-Array (N×N×numFrames)
    %
    [N, ~]      = size(b);
    sequence    = zeros(N, N, numFrames);
    k0          = ceil(numFrames / 2);
    for k = 1:numFrames
        shift       = v * (k - k0);            
        sequence(:,:,k) = imtranslate(b, ...
                               [shift(1), shift(2)], ...
                               'linear', 'FillValues', 0);
    end
end

% 2.3 Arithmetische Mittelung (Software-Pulsing) nach Gleichung (3.3)
bp = mean(sequence, 3);

figure;
imshow(mat2gray(bp));
title('Gepulst belichtetes Bild bp(x;v)');

%% 3. Fourier-Analyse & logarithmisches Spektrum
% 3.1 Berechnung der 2D-FFT und Zentrierung via fftshift
B    = fftshift(fft2(bp));                 
absB = abs(B);                             
logB = log(1 + absB);  % Dynamikanpassung für Visualisierung

figure;
imshow(mat2gray(logB));
title('Logarithmiertes Betragsspektrum von bp');

%% 4. Radon-Transformation des Log-Spektrums
% Wendet die Radon-Transformation R(r, θ) auf logB an und erzeugt ein
% sinogramm-artiges Bild, in dem impulsartige Frequenzstrukturen als
% Linien zu erkennen sind.
theta = 0:1:179;                   % Diskrete Winkelauflösung
[R, xp] = radon(logB, theta);

figure;
imagesc(theta, xp, R);
xlabel('θ [°]');
ylabel('r [Pixel]');
title('Radon-Transformation des log-Spektrums');
colorbar;

%% 5. Fluss-Analyse: Bestimmung von Winkel und Betrag
% 5.1 Welligkeitsprofil WR(θ):
%     WR(θ) = ∑_r |∂/∂r R(r, θ)|  – Maß für die "Geradlinigkeit" der Energie
numAngles = numel(theta);
WR = zeros(1, numAngles);
for i = 1:numAngles
    proj      = R(:, i);
    gradProj  = abs(diff(proj));       % Finite Differenzen entlang r
    WR(i)     = sum(gradProj);         % Summierte Kantendichte
end

% 5.2 Geschwindigkeitsrichtung als Argmax von WR(θ)
[~, bestIdx]    = max(WR);
v_angle_est     = theta(bestIdx);

% 5.3 Abschätzung des Geschwindigkeitsbetrags:
%     Abstand der Spitzen in R(r, θ_best) ∼ 1/|v|
selectedProj    = R(:, bestIdx);
[pks, locs]     = findpeaks(selectedProj, 'MinPeakProminence', 0.05);

if numel(locs) > 1
    spacings    = diff(locs);
    peakSpacing = mean(spacings);     % Mittlerer Peak-Abstand
else
    peakSpacing = NaN;
end

% Darstellung und Ausgabe
fprintf('--- Fluss-Analyse (Kap. 3.2) ---\n');
fprintf('Geschätzter Winkel: %.2f°\n', v_angle_est);
fprintf('Mittlerer Peak-Abstand: %.2f Pixel-Indizes\n\n', peakSpacing);

figure;
subplot(2,1,1);
plot(theta, WR, 'LineWidth', 1.5);
xlabel('θ [°]');
ylabel('WR(θ)');
title('Welligkeitsprofil WR vs. Winkel');
grid on;

subplot(2,1,2);
plot(selectedProj, 'r-', 'LineWidth', 1.5);
hold on;
plot(locs, pks, 'ko', 'MarkerFaceColor','y');
xlabel('r [Pixel]');
ylabel('R(r, θ_{best})');
title(sprintf('Radon-Projektion bei %.0f° mit Peaks', v_angle_est));
grid on;

%% … nach deiner Flow-Analyse (Kap. 3.2) …

% 5.4 Wahrer Betrag des Geschwindigkeitsvektors
v_true_mag = norm(params.v);

% 5.5 Geschätzter Betrag (relativ, invers proportional zum mittleren Peak-Abstand)
if ~isnan(peakSpacing)
    v_est_mag = 1 / peakSpacing;
else
    v_est_mag = NaN;
end

% 5.6 Ausgabe von Winkel und Beträgen
fprintf('--- Geschwindigkeitsergebnisse ---\n');
fprintf('Wahrer Betrag       : %.2f Pixel/Frame\n', v_true_mag);
fprintf('Geschätzter Winkel  : %.2f°\n', v_angle_est);
fprintf('Geschätzter Betrag  : %.2f (relative Einheit)\n\n', v_est_mag);

%% 6. Auflösungsreduktion & FFT-Vergleich (Kap. 3.3)
% 6.1 Pyramid-Down: Anti-Aliasing geführtes Herunterskalieren auf 128×128
bp_pyr = imresize(bp, params.downFactor, 'bilinear');

% 6.2 Zentraler Ausschnitt (Crop) mit identischer Zielgröße
startIdx = floor((params.N - params.cropSize)/2) + 1;
bp_crop  = bp(startIdx:startIdx+params.cropSize-1, ...
              startIdx:startIdx+params.cropSize-1);

% 6.3 Erneute Fourier-Analyse für alle drei Varianten
logB_full = log(1 + abs(fftshift(fft2(bp))));
logB_pyr  = log(1 + abs(fftshift(fft2(bp_pyr))));
logB_crop = log(1 + abs(fftshift(fft2(bp_crop))));

figure;
subplot(1,3,1); imshow(mat2gray(logB_full));  title('FFT (256×256)');
subplot(1,3,2); imshow(mat2gray(logB_pyr));   title('FFT (Pyramid-Down)');
subplot(1,3,3); imshow(mat2gray(logB_crop));  title('FFT (Crop)');

%% 7. Selbstbeurteilung: Gütemaß GEBA (Kap. 3.4)
% Definiert als Verhältnis aus dominante Linie vs. durchschnittliche
% Welligkeit im Sinogramm:
WR_max = max(WR);
GEBA   = WR_max / mean(WR);

fprintf('--- Selbstbeurteilungsmaß GEBA (Kap. 3.4) ---\n');
fprintf('GEBA = %.3f (max WR / mean WR)\n', GEBA);

figure;
plot(theta, WR, 'm-', 'LineWidth', 2);
xlabel('θ [°]');
ylabel('WR(θ)');
title(sprintf('Welligkeitsprofil mit GEBA = %.3f', GEBA));
grid on;



