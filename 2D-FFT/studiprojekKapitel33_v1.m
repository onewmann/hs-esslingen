%% Chapter3_1_3.m
% MATLAB-Skript zu Kapitel 3.1.3: Summation aufeinander folgender Bilder
% Dieses Skript simuliert:
%  - ein statisches Texturbild b(x)
%  - eine Bildfolge {b_k(x)}, bei der das Bild jeweils um v*(k-k0) verschoben wird,
%  - die Summation (Mittelwert) der Bilder als "Software-Pulsing".

clear; close all; clc;

%% 1. Parameter
N = 256;                   % Bildgröße: N x N Pixel
A = 0.3;                   % Texturparameter (steuert die Breite der Autokorrelationsfunktion)
v = [8, 4];                % Geschwindigkeitsvektor [vx, vy] in Pixel pro Frame
L = 2;                     % Es werden 2L+1 Bilder summiert (hier L=2: 5 Bilder)
dt = 1;                    % Zeitabstand zwischen den Aufnahmen (hier als 1 angenommen)

numFrames = 2 * L + 1;     % Anzahl der Einzelbilder

%% 2. Statisches Textur Bild
%b = simulateTexture(N, A);
%% 2. Statisches Texturbild & Bildfolge
b = simulateTexture(N, A);
sequence = generateImageSequence(b, v, numFrames);
bp = mean(sequence, 3);  % Summiertes (bewegungsunschärfes) Bild

figure;
imshow(mat2gray(bp));
title('Summiertes Bild bp(x;v)');
%% 3. Erzeugung der Bildfolge {b_k(x)}
% Mittlerer Frame (k0) soll unverschoben sein.
sequence = generateImageSequence(b, v, numFrames);
figure;
nShow = min(5, numFrames);
for k = 1:nShow
    subplot(1, nShow, k);
    imshow(mat2gray(sequence(:, :, k)));
    title(sprintf('Frame %d', k));
end
sgtitle('Bildfolge {b_k(x)}');

%% 4. Summation der Einzelbilder (Software-Pulsing)
bp_sum = mean(sequence, 3);
figure;
imshow(mat2gray(bp_sum));
title('Summiertes Bild bp(x;v) (Software-Pulsing)');

%% 5. Spektralanalyse des summierten Bildes
B = fftshift(fft2(bp_sum));
absB = abs(B);
logB = log(1+absB);
figure;
imshow(mat2gray(logB));
title('Logarithmiertes Betragsspektrum von bp(x;v)');

