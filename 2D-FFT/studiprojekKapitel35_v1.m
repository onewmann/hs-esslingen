%% Chapter3_3.m
% MATLAB-Skript zu Kapitel 3.3: Reduktion der Pixelanzahl in der Einzelbild-Analyse
%
% Ziel: Demonstration, wie eine Auflösungsreduktion (zum Beispiel zur Kostenreduktion)
% Einfluss auf das summierte Bild bp(x;v) und dessen spektrale Darstellung hat.
%
% Vorgehensweise:
%   1. Erzeuge ein statisches Texturbild b(x) (mittels simulateTexture).
%   2. Simuliere eine Bildfolge {b_k(x)} bei Übersetzung um einen Geschwindigkeitsvektor v.
%   3. Berechne das summierte (gepulst belichtete) Bild bp(x;v) durch Mittelung.
%   4. Reduziere die Pixelanzahl mit zwei Ansätzen: 
%        a) Pyramid-Down: mittels imresize.
%        b) Bildausschnitt: zentraler Crop.
%   5. Vergleiche die Bilder und deren 2D-FFT-Spektren.

clear; close all; clc;

%% 1. Parameter
N = 256;           % Originalbildgröße: N x N Pixel
A = 0.3;           % Texturparameter (Steuerung der Breite der AKF)
v = [8, 4];        % Geschwindigkeitsvektor [vx,vy] (Pixel pro Frame)
L = 3;             % Es werden 2L+1 Bilder summiert (hier L=3 -> 7 Bilder)
dt = 1;            % dt=1 (als Referenz, da wir reine relative Verschiebungen simulieren)
numFrames = 2 * L + 1;

%% 2. Erzeuge statisches Texturbild
b = simulateTexture(N, A);
figure;
imshow(mat2gray(b));
title('Statisches Texturbild b(x)');

%% 3. Bildfolge und Summation (Software-Pulsing)
sequence = generateImageSequence(b, v, numFrames);
bp = mean(sequence, 3);

figure;
imshow(mat2gray(bp));
title('Summiertes Bild bp(x;v) (256 x 256)');

%% 4. Reduktion der Pixelanzahl

% Option A: Pyramid-Down
downFactor = 0.5;  % Faktor 0.5: Reduktion auf 128x128 Pixel
bp_pyr = imresize(bp, downFactor, 'bilinear'); % Verwendet integriertes Anti-Aliasing
figure;
imshow(mat2gray(bp_pyr));
title('Gepulst belichtetes Bild bp(x;v) - Pyramid Down (128 x 128)');

% Option B: Bildausschnitt (Crop)
cropSize = 128;    % Zielgröße: 128 x 128 Pixel
startIdx = floor((N - cropSize)/2) + 1;
bp_crop = bp(startIdx:startIdx+cropSize-1, startIdx:startIdx+cropSize-1);
figure;
imshow(mat2gray(bp_crop));
title('Gepulst belichtetes Bild bp(x;v) - Bildausschnitt (128 x 128)');

%% 5. Spektralanalyse: Vergleich der FFT-Spektren
% Full resolution FFT:
B_full = fftshift(fft2(bp));
logB_full = log(1 + abs(B_full));

% Pyramid-Down FFT:
B_pyr = fftshift(fft2(bp_pyr));
logB_pyr = log(1 + abs(B_pyr));

% Crop FFT:
B_crop = fftshift(fft2(bp_crop));
logB_crop = log(1 + abs(B_crop));

figure;
subplot(1,3,1);
imshow(mat2gray(logB_full));
title('FFT-Spektrum (256x256)');

subplot(1,3,2);
imshow(mat2gray(logB_pyr));
title('FFT-Spektrum (Pyramid Down, 128x128)');

subplot(1,3,3);
imshow(mat2gray(logB_crop));
title('FFT-Spektrum (Crop, 128x128)');

