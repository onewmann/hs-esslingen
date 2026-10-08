%% Kapitel 2.1: Bilderfassung – Vergleich von konstanter und gepulster Belichtung
% Dieses Skript simuliert die Bildaufnahme einer texturierten Oberfläche in zwei Modi:
%
%   1. Konstante Belichtung: Es erfolgt eine numerische Integration über die 
%      Gesamtbelichtungszeit T_total, sodass das aufgenommene Bild einer kontinuierlichen
%      Belichtung entspricht.
%
%   2. Gepulste Belichtung: Die Beleuchtung wird nur in diskreten Impulsen aktiviert.
%      Für jeden Lichtimpuls wird ein eigenes Bild aufgenommen, das in einem separaten
%      Fenster ausgegeben wird.
%
% Das vorliegende Modell nimmt an, dass sich die Textur durch Translation (definiert durch den
% Geschwindigkeitsvektor v) bewegt. Die Modellierung erfolgt im Rahmen eines einfachen
% pinhole-Kameramodells, wobei ein normiertes Koordinatengitter den kontinuierlichen Raum repräsentiert.
%
% Autor: [Dein Name]
% Datum: [Datum]

clear; close all; clc;

%% Parameterdefinition
% Bild- und Raumparameter
N = 256;  % Bildgröße: N x N Pixel
% Erzeuge ein normiertes kontinuierliches Koordinatengitter im Intervall [-1, 1] für beide Dim.
[x, y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));

%% Modellierung der Textur
% Die Textur wird mit einer Gaußschen Verteilungsfunktion modelliert, welche eine
% isotrope Oberflächentextur approximiert. Additives Rauschen simuliert natürliche Störeinflüsse.
sigma = 0.3;  % Standardabweichung der Gaußfunktion (räumliche Ausdehnung)
texture = exp(-((x.^2 + y.^2) / (2 * sigma^2)));

% Additiver, normalverteilter Rauschterm
noiseLevel = 0.1;
textureNoisy = texture + noiseLevel * randn(size(texture));

% Normierung der Intensitätswerte auf den Intervall [0, 1]
textureNoisy = mat2gray(textureNoisy);

%% Bewegung und Belichtungsparameter
T_total = 1.0;           % Gesamtbelichtungszeit in Sekunden
v = [20, 10];            % Geschwindigkeitsvektor in Pixel pro Sekunde [dx, dy]
                         % (d.h., die Textur wird während T_total um v*T_total verschoben)

% Parameter für konstant belichteten Fall: 
numSteps_const = 50;     % Anzahl der kleinen Zeitschritte für die numerische Integration
t_const = linspace(0, T_total, numSteps_const);

% Parameter für gepulste Belichtung:
numPulses = 5;           % Anzahl der Lichtpulse innerhalb von T_total
t_pulses = linspace(0, T_total, numPulses);

%% A) Aufnahme bei konstanter Belichtung (kontinuierliche Integration)
% Modell: l(t) = 1/T_total ist konstant während T_total.
% Das resultierende Bild entspricht der gewichteten Summe der zur jeweiligen Zeit t verschobenen Bilder.
constantExposureImage = zeros(N, N);
dt_const = T_total / numSteps_const;
for idx = 1:numSteps_const
    t = t_const(idx);
    shiftedImage = imtranslate(textureNoisy, v * t, 'linear', 'FillValues', 0);
    constantExposureImage = constantExposureImage + (1/T_total)*shiftedImage * dt_const;
end

%% B) Aufnahme bei gepulster Belichtung: Ein Bild pro Lichtpuls
% Es wird angenommen, dass die Beleuchtung nur in diskreten Impulsen erfolgt.
% Für jeden Impuls (Zeitpunkt t_puls) wird ein eigenes Bild aufgenommen.
pulsedExposureImages = zeros(N, N, numPulses);  % 3D-Array: [Zeile, Spalte, Impuls]
for idx = 1:numPulses
    t = t_pulses(idx);
    pulsedExposureImages(:,:,idx) = imtranslate(textureNoisy, v * t, 'linear', 'FillValues', 0);
end

%% Ausgabe in separaten Fenstern (Extrafenster)
% Statisches Texturbild in einem eigenen Fenster:
figure;
imshow(textureNoisy);
title('Statisches Texturbild (Gaußsche Textur mit Rauschen)');

% Konstante Belichtung: Ergebnisbild in einem eigenen Fenster:
figure;
imshow(constantExposureImage);
title(['Bild bei konstanter Belichtung, T_{total} = ' num2str(T_total) ' s']);

% Gepulste Belichtung: Für jeden Puls ein separates Fenster
for idx = 1:numPulses
    figure;
    imshow(pulsedExposureImages(:,:,idx));
    title(sprintf('Bild bei gepulster Belichtung: Puls %d (t = %.2f s)', idx, t_pulses(idx)));
end

%% Optionale Animationen
% Animation zur Visualisierung der kontinuierlichen Translation (konstante Belichtung)
figure;
for idx = 1:numSteps_const
    t = t_const(idx);
    shiftedImage = imtranslate(textureNoisy, v * t, 'linear', 'FillValues', 0);
    imshow(shiftedImage);
    title(sprintf('Konstante Belichtung (Simulation): t = %.2f s', t));
    pause(0.1);
end

% Animation der gepulsten Bilder
figure;
for idx = 1:numPulses
    imshow(pulsedExposureImages(:,:,idx));
    title(sprintf('Gepulster Puls: %d, t = %.2f s', idx, t_pulses(idx)));
    pause(0.5);
end
