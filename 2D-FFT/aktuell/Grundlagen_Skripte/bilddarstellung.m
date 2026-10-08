% system_theory_image_acquisition_demo_commented.m
% Ausführlich kommentiertes MATLAB-Demo-Skript zur Illustration der
% systemtheoretischen Formulierung der Bildaufnahme:
%    b2(x;v) = b(x) * h(x;v)
% und im Frequenzraum:
%    |BΣ(f;v)| ≈ |B(f)| · |H(f;v)|
%
% Autor: für dich angepasst (Oliver)
% Hinweise: MATLAB R2018b oder neuer empfohlen. Benötigt keine zusätzlichen Toolboxes.

clear; close all; clc;

%% -------------------- Parameter --------------------
imgSize      = 512;          % Größe des quadratischen Bildes (Pixel)
textureType  = 'white'; % 'white' für weißes Rauschen, 'structured' für Muster + Rauschen
velocity     = [5, 5];       % Bewegungs-/Aufnahmevektor nu = [vx, vy] (einfaches Modell)
kernelLength = 61;           % Länge des linearen Kernels (ungerade Zahl empfohlen)
sigma_spread = 0.8;          % "Splat"-Breite beim Rasterisieren der Linie (Nähe zu Delta)

% Numerische Stabilität für Log/Anzeige
epss = 1e-12;

%% -------------------- 1) Erzeuge Textur b(x) --------------------
% Zwei Varianten: 'white' ergibt annähernd flaches Spektrum (gut zum Sichtbarmachen von H)
% 'structured' erzeugt sichtbare Texturfrequenzen (realistischere Situation)
switch textureType
    case 'white'
        % Gaußsches weißes Rauschen (Nullmittel, Varianz 1)
        B = randn(imgSize);
    case 'structured'
        % Kombination aus zwei Sinusmustern (unterschiedliche Frequenzen) plus Rauschen
        [X, Y] = meshgrid(linspace(0,8,imgSize));         % räumliche Koordinaten
        B = 0.6*sin(2*pi*0.5*X) + 0.4*sin(2*pi*1.2*Y);   % zwei dominante Textur-Komponenten
        B = B + 0.25*randn(imgSize);                     % etwas Rauschen dazu
    otherwise
        error('Unbekannter textureType: benutze ''white'' oder ''structured''.');
end

% Normierung: Mittelwert 0, Standardabweichung 1 (macht Spektren vergleichbar)
B = B - mean(B(:));
B = B / max(1e-12, std(B(:)));

%% -------------------- 2) Konstruiere Kernel h(x;v) --------------------
% Ziel: approximativer "Delta"-Kern entlang der Bewegungsrichtung.
% Im Buch wird h(x;v) oft mit Dirac-Impulsen modelliert; hier nutzen wir eine
% diskrete Näherung: eine dünne Linie (subpixel) "gesplatet" mit Gauß-Weights.
vx = velocity(1);
vy = velocity(2);
theta = atan2(vy, vx);         % Winkel phi der Bewegungsrichtung (Bogenmaß)

% Kernel-Größe: len x len (kleiner als das Bild). Später normalisieren wir Energie.
len = kernelLength;
cx = (len+1)/2;                % Kernzentrum (Mittelpunkt)

% Leeres Kernel anlegen
h = zeros(len, len);

% Parameter t läuft entlang der Linie; wir platzieren len Stützpunkte
t = linspace(-floor(len/2), floor(len/2), len);

% Subpixel-Positionen der Linienpunkte im Kernel-Koordinatensystem
x_line = cx + (t * cos(theta));
y_line = cx + (t * sin(theta));

% Gausssche "Splat"-Näherung: jedes Linien-Stützpunkt trägt eine Gauß-Glocke im Kernel bei
% Dadurch entsteht eine glatte, diskrete Näherung eines Linienimpulses.
[Xk, Yk] = meshgrid(1:len, 1:len);
for i = 1:len
    dx = Xk - x_line(i);
    dy = Yk - y_line(i);
    h = h + exp(-(dx.^2 + dy.^2) / (2 * sigma_spread^2));
end

% Normiere Kernel-Energie auf 1 (Summe = 1), damit Kernel eine Impulsantwort ist
h = h / sum(h(:));

%% -------------------- 3) Ortsraum-Faltung: b2 = b * h --------------------
% Wir verwenden imfilter mit 'conv' für Konvolution; 'replicate' vermeidet harte Ränder.
% Resultat b2 ist das beobachtete Bild b_2(x;v)
b2 = imfilter(B, h, 'conv', 'same', 'replicate');

%% -------------------- 4) Frequenzbereich: FFTs und Darstellungen --------------------
% Wir definieren eine kleine Helferfunktion für zentrierte (shifted) FFT mit Padding
fftCentered = @(I, sz) fftshift( fft2(I, sz(1), sz(2)) );

% FFT des Originalbildes (auf Bildgröße)
FB  = fftCentered(B,  [imgSize imgSize]);   % Fourier von b(x), zentriert
% FFT des Kernels: Kernel vor der FFT auf Bildgröße gepadded (entspricht Zero-Padding)
FH  = fftCentered(h,  [imgSize imgSize]);   % Fourier von h(x;v) (padded)
% FFT des gefalteten Bildes
FB2 = fftCentered(b2, [imgSize imgSize]);   % Fourier von b2(x;v)

% Betragsspektren (Amplitude)
magFB  = abs(FB);
magFH  = abs(FH);
magFB2 = abs(FB2);

% Logarithmische Darstellung für bessere Visualisierung (mit kleinen Offset)
S_B  = log(magFB  + epss);
S_H  = log(magFH  + epss);
S_B2 = log(magFB2 + epss);

% Numerischer Test der Beziehung |BΣ| ≈ |B| · |H|
prodSpec = magFB .* magFH;
% Relative Fehlernorm (L2-Fehler / Norm von magFB2)
relError = norm(magFB2(:) - prodSpec(:)) / max(1e-12, norm(magFB2(:)));

%% -------------------- 5) Visualisierung --------------------
figure('Name','Systemtheorie Bildaufnahme (ausführlich kommentiert)','Units','normalized','Position',[0.05 0.05 0.9 0.85]);

% Links: Bilder; Mitte/ rechts: Spektren / Prüfungen
subplot(3,4,1);
imagesc(B); axis image off; colormap gray; title('b(x) Textur');
subplot(3,4,2);
imagesc(h); axis image off; colormap gray; title('h(x;v) Kernel (linear approx.)');
subplot(3,4,3);
imagesc(b2); axis image off; colormap gray; title('b_2(x;v) Gefaltetes Bild');
subplot(3,4,4);
text(0.02,0.5, sprintf(['velocity = [%.1f, %.1f]\n|nu| = %.2f\nangle = %.1f°\nrelErr FFT-prod = %.3e'], ...
    vx, vy, norm(velocity), rad2deg(theta), relError), 'FontSize', 10);
axis off;

% Log-Betragsspektren
subplot(3,4,5);
imagesc(S_B); axis image off; colormap jet; title('|B(f)| (log)');
subplot(3,4,6);
imagesc(S_H); axis image off; colormap jet; title('|H(f;v)| (log)');
subplot(3,4,7);
imagesc(S_B2); axis image off; colormap jet; title('|B_{\\Sigma}(f;v)| (log)');
subplot(3,4,8);
imagesc(log(prodSpec + epss)); axis image off; colormap jet; title('log(|B|·|H|)');

% 1D-Spektralschnitt entlang Richtung phi (zeigt Ausrichtung/Skalierungseffekt)
cxI = imgSize/2 + 1;
cyI = imgSize/2 + 1;
L = 300;                                 % Anzahl Samples entlang Schnitt
t = linspace(-120, 120, L);              % Parameter entlang Linie im Frequenzraum
fx = round(cxI + t * cos(theta));        % x-Indices entlang Linie
fy = round(cyI + t * sin(theta));        % y-Indices entlang Linie
valid = fx>=1 & fx<=imgSize & fy>=1 & fy<=imgSize;
fx = fx(valid); fy = fy(valid);
specB_line  = magFB (sub2ind(size(magFB),  fy, fx));
specH_line  = magFH (sub2ind(size(magFH),  fy, fx));
specB2_line = magFB2(sub2ind(size(magFB2), fy, fx));

subplot(3,4,9);
plot(specB_line,'b-','LineWidth',1.0); hold on;
plot(specH_line,'r-','LineWidth',1.0);
plot(specB2_line,'k-','LineWidth',1.0);
legend('|B|','|H|','|BΣ|','Location','best'); title('Spektralschnitte entlang \phi'); xlabel('Samples'); grid on;

% Zoom auf Kernel-Spektrum (zentraler Bereich) — macht Richtungscharakter sichtbar
subplot(3,4,10);
imagesc(S_H); axis image off; colormap jet; title('Zoom |H(f;v)| (log)');

% Differenzkarte: |BΣ| - |B|·|H| (zeigt numerische Diskrepanzen)
diffMap = magFB2 - prodSpec;
subplot(3,4,11);
imagesc(diffMap); axis image off; colormap parula; colorbar; title('|BΣ| - |B|·|H|');

% Histogramm der Differenzen (Statistik)
subplot(3,4,12);
histogram(diffMap(:), 80); title('Histogramm(|BΣ|-|B|·|H|)'); xlabel('Differenz');

sgtitle('Demo: b_2(x;v)=b(x)*h(x;v) und |BΣ| ≈ |B|·|H| (ausführlich kommentiert)', 'FontSize', 14);

%% -------------------- 6) Kurze numerische Ausgabe / Hinweise --------------------
fprintf('Relative Fehler zwischen |FFT(b2)| und |FFT(b)|·|FFT(h)|: %.3e\n', relError);
fprintf('Erläuterung: Numerische Fehler kommen von Diskretisierung, Padding, Ränder, Rundung.\n');
fprintf('Tipp: Setze textureType=''white'' um |H| deutlicher zu sehen (nahezu flaches |B|).\n');

% Ende des Skripts
