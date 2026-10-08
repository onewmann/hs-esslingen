%% Chapter3_2.m
% MATLAB-Skript zu Kapitel 3.2: Extraktion von Betrag und Richtung
%
% Dieses Skript simuliert den Ablauf, mit dem aus einem
% gepulst belichteten Bild (b_p) der Geschwindigkeitsvektor in Betrag und
% Richtung extrahiert wird.
%
% Vorgehensweise:
%   1. Generiere ein statisches Texturbild b(x)
%   2. Erzeuge ein Bewegungsbild bp(x;v) via gepulster Belichtung:
%         bp(x;v)= (1/(2L+1)) * sum_{k=-L}^{L} b(x - v*(k*dt))
%   3. Berechne mittels 2D‑FFT das logarithmierte Betragsspektrum.
%   4. Wende die radon-Transformation auf das Spektrum an.
%   5. Für jeden Winkel berechne einen "Welligkeitswert" (WR), indem
%      die Differenzen entlang der \(r\)-Richtung summiert werden.
%   6. Der Winkel mit maximalem WR gibt den geschätzten Geschwindigkeitswinkel.
%   7. In der entsprechenden Radon‑Spalte werden die Peak-Abstände (z.B. per findpeaks)
%      gemessen, was proportional zum Geschwindigkeitsbetrag ist.
%
% Hinweis: Die exakte Skalierung (Umrechnung in physikalische Einheiten)
% ist hier nicht voll implementiert – das Skript zeigt den prinzipiellen Ablauf.
% Autor: Oliver Neumann
% Datum: 22.09


clear; close all; clc;

%% 1. Parameter und statisches Texturbild
N = 256;                % Bildgröße: N x N Pixel
A = 0.3;                % Texturparameter (steuert die Breite der Textur-AKF)
v_true = [-9, 9];        % Wahrer Geschwindigkeitsvektor in Pixel pro Sekunde
% (bei der Simulation hier angenommen als "Pixel pro Frame", da dt = 1 verwendet wird)
L = 6;                  % Anzahl der Impulse: 2L+1 (hier 7 Bilder)
dt = 1;                 % Zeitabstand zwischen den Impulsen (wir setzen dt=1, da wir
                        % in diesem Teil ausschließlich den relativen Effekt simulieren)

% Erzeuge statisches Texturbild
b = simulateTexture(N, A);
figure;
imshow(mat2gray(b));
title('Statisches Texturbild b(x)');

%% 2. Erzeugung eines Bildes via Software-Summation
% In Kapitel 3.1.3 (Software-Pulsing) wird eine Bildfolge {b_k(x)} erzeugt,
% wobei der mittlere Frame unverschoben ist und die restlichen Frames
% um v*(k - k0) verschoben werden.
numFrames = 2 * L + 1;  % Zahl der Einzelbilder
sequence = generateImageSequence(b, v_true, numFrames);

% Summiere die Einzelbilder (arithmetisches Mittel)
bp = mean(sequence, 3);
figure;
imshow(mat2gray(bp));
title('Summiertes Bild bp(x;v) (Software-Pulsing)');

%% 3. Fourier Transformation und Log Spektrum
B = fftshift(fft2(bp));
absB = abs(B);
logB = log(1 + absB);  % logarithmierte Darstellung zur besseren Sichtbarkeit

figure;
imshow(mat2gray(logB));
title('Logarithmiertes Betragsspektrum von bp(x;v)');

%% 4. Radon-Transformation auf das log-Spektrum
% Wähle den Winkelbereich (0 bis 179° in 1°-Schritten)
theta = 0:1:179;
[R, xp] = radon(logB, theta);

% Visualisierung der Radon-Transformation (optional)
figure;
imagesc(theta, xp, R);
xlabel('Winkel (°)'); ylabel('r'); colorbar;
title('Radon-Transformation des log-Spektrums');

%% 5. Extraktion des Geschwindigkeitswinkels
% Für jede Radon-Projektion berechnen wir den summierten absoluten
% Gradienten (WR). Dies liefert ein Maß dafür, wie stark
% die Impulslinien in der jeweiligen Projektion ausgeprägt sind.
numAngles = length(theta);
WR = zeros(1, numAngles);
for i = 1:numAngles
    proj = R(:, i);
    gradProj = abs(diff(proj));   % absolute Ableitung entlang r
    WR(i) = sum(gradProj);        % Summe als Maß der "Kantenstärke"
end

% Finde den Winkel, bei dem WR maximal ist – dieser entspricht dem geschätzten Geschwindigkeitswinkel
[~, bestIdx] = max(WR);
v_angle_est = theta(bestIdx);
fprintf('Geschätzter Geschwindigkeitswinkel: %.2f°\n', v_angle_est);

% Plot WR versus Winkel
figure;
plot(theta, WR, 'b-', 'LineWidth', 2);
xlabel('Winkel (°)');
ylabel('Summierter Gradient WR');
title('Auswertung der Radon-Projektionen');

%% 6. Extraktion des Geschwindigkeitsbetrags
% In der Radon-Projektion bei bestIdx messen wir nun die Abstände zwischen
% Peaks. Diese Abstände sind (in diesem Modell) proportional zum
% Betragsanteil v_mag.
selectedProj = R(:, bestIdx);
% Optional: Plot der Projektion
figure;
plot(selectedProj, 'LineWidth', 2);
xlabel('r index'); ylabel('Projektion');
title(sprintf('Radon-Projektion bei %.0f°', v_angle_est));

% Finde Peaks in der Projektion (Parameter können je nach Textur angepasst werden)
[pks, locs] = findpeaks(selectedProj, 'MinPeakProminence', 0.05);
if length(locs) > 1
    spacings = diff(locs);
    spacing_avg = mean(spacings);
else
    spacing_avg = NaN;
end
fprintf('Durchschnittlicher Abstand der Peaks in der Radon-Projektion: %.2f (willkürliche Einheit)\n', spacing_avg);

% Hinweis: Die exakte Beziehung zwischen Peak-Abstand und Geschwindigkeitsbetrag
% ist in der Dissertation theoretisch hergeleitet; hier vergleichen wir nur
% den relativen Wert zum wahren Geschwindigkeitsbetrag
v_true_mag = norm(v_true);
fprintf('Wahrer Geschwindigkeitsbetrag (Pixel/s): %.2f\n', v_true_mag);

%% 7. Interpretation und Zusammenfassung
% In diesem Beispiel stellst du fest:
%  - v_angle_est entspricht dem geschätzten Winkel des Geschwindigkeitsvektors.
%  - spacing_avg liefert einen Hinweis auf den Geschwindigkeitsbetrag.
% Die Skalierung darf hier ohne exakte Korrekturfaktoren als relativer Vergleich
% interpretiert werden.
