%% Chapter3_4_preproc_radon_contrast_fullscript.m
% Einzelbild-Analyse mit Background-Subtraktion, Rauschminderung,
% Kontrastverstärkung, Windowing, Zero-Padding, Radon-Transform und
% subpixelgenauer Peak-Fitting

clear; close all; clc;

%% 0. IMAGE-PATH VALIDIERUNG UND AUSWAHL
% Platzhalter-Pfad, wird abgefragt falls nicht existent
imgPath = "C:\Users\onzen\OneDrive\Dokumente\2D-FFT_Studienprojekt\data\image001.png";

if ~isfile(imgPath)
    [fname, fpath] = uigetfile( ...
        {'*.png;*.jpg;*.tif','Bilddateien (*.png,*.jpg,*.tif)'; ...
         '*.*','Alle Dateien'}, ...
        'Wählen Sie eine Input-Bilddatei aus');
    if isequal(fname,0)
        error('Kein Bild ausgewählt. Abbruch.');
    end
    imgPath = fullfile(fpath, fname);
end

%% 1. PARAMETER
v_true    = [50,  53];    % Wahrer Geschw.-Vektor [vx, vy] in px/Frame
L         = 80;            % Halbbreite der Sequenz
numFrames = 2*L + 1;       % Gesamtzahl der Frames
theta     = 0:1:179;       % Radon-Winkelbereich [0°,180)
N         = 512;           % Quadratbildgröße

%% 2. BILD LADEN & SEQUENZ ERZEUGEN
I = im2double(imread(imgPath));
if size(I,3)==3
    I = rgb2gray(I);
end
I = imresize(I, [N, N]);

sequence = generateImageSequence(I, v_true, numFrames);

% Summenbild erzeugen
bp = mean(sequence, 3);

%% 3. VORVERARBEITUNG
% 3.1 Background-Subtraktion (DC-Entfernen)
bg      = imopen(bp, strel('disk', round(N/20)));
bp_dc   = bp - bg;
bp_dc(bp_dc<0) = 0;

% 3.2 Rauschminderung (Wiener-Filter)
bp_denoised = wiener2(bp_dc, [5 5]);

% 3.3 Kontrastverstärkung (Stretchlim + imadjust)
lims        = stretchlim(bp_denoised, [0.01 0.99]);
bp_contrast = imadjust(bp_denoised, lims, []);

figure('Name','Vorverarbeitung');
subplot(1,3,1), imshow(bp),            title('Original Summenbild');
subplot(1,3,2), imshow(bp_denoised),   title('DC-entfernt + Wiener');
subplot(1,3,3), imshow(bp_contrast),   title('Kontrastverstärkt');

%% 4. WINDOWING & ZERO-PADDING
% 4.1 Hanning-Fenster
h1     = hann(N);
h2d    = h1 * h1.';
bp_win = bp_contrast .* h2d;

% 4.2 Zero-Padding auf doppelte Größe
M      = 2 * N;
pad    = (M - N) / 2;
bp_pad = padarray(bp_win, [pad pad], 0, 'both');

%% 5. 2D-FFT & LOG-BETRAGSSPEKTRUM
B    = fftshift(fft2(bp_pad));
logB = log(1 + abs(B));

figure('Name','Log-Spektrum');
imshow(mat2gray(logB));
title('log(1 + |FFT\{bp\}|)');

%% 6. RADON-TRANSFORMATION
[R, xp] = radon(logB, theta);

figure('Name','Radon(logB)');
imagesc(theta, xp, R);
xlabel('Winkel (°)'); ylabel('r (Pixel)');
title('Radon-Transformierte des Spektrums');
colorbar;

%% 7. WELLIGKEITSPROFIL WR(θ)
WR = sum(abs(diff(R,1,1)), 1);

% Winkel schätzen
[WR_max, idxR] = max(WR);
ridge_angle    = theta(idxR);
v_angle_est    = mod(ridge_angle, 180);   % kein +90°
GEBA           = WR_max / mean(WR);

figure('Name','WR Profil');
plot(theta, WR, 'LineWidth',2);
xlabel('θ (°)'); ylabel('WR');
title(sprintf('Welligkeitsprofil WR (GEBA = %.3f)', GEBA));
grid on;

%% 8. BETRAGSSCHÄTZUNG MIT SUBPIXEL-PEAK-FITTING
proj     = R(:, idxR);
mask     = xp >= 0;
xp_pos   = xp(mask);
proj_pos = proj(mask);

% Peaks detektieren
[pks, locs] = findpeaks(proj_pos, xp_pos, ...
    'MinPeakProminence', 0.05, ...
    'MinPeakDistance',    3);

% Subpixel-Parabel-Fit pro Peak
locs_sub = nan(size(locs));
for i = 1:numel(locs)
    [~, ix] = min(abs(xp_pos - locs(i)));
    if ix>1 && ix<length(proj_pos)
        y_m1 = proj_pos(ix-1);
        y_0  = proj_pos(ix);
        y_p1 = proj_pos(ix+1);
        denom   = 2*(y_m1 - 2*y_0 + y_p1);
        offset  = (y_m1 - y_p1) / denom;
        locs_sub(i) = locs(i) + offset;
    else
        locs_sub(i) = locs(i);
    end
end

% Mittlerer Abstand und Geschwindigkeit (M statt N!)
delta_r   = mean(diff(sort(locs_sub)));
v_mag_est = (M / delta_r);  % px/frame

figure('Name','Radon-Projektion + Subpixel');
plot(xp_pos, proj_pos, 'b-','LineWidth',1.5); hold on;
plot(locs,   pks,        'ro','MarkerSize',8);
plot(locs_sub, interp1(xp_pos, proj_pos, locs_sub), ...
     'kx','MarkerSize',8);
xlabel('r (Pixel)'); ylabel('Projektion');
title(sprintf('Peaks & Subpixel-Fitting bei %.1f°', ridge_angle));
legend('Projektion','Orig Peaks','Subpixel');

%% 9. AUSGABE
fprintf('Wahrer Betrag      : %.2f px/frame\n', norm(v_true));
fprintf('Geschätzter Betrag : %.2f px/frame\n', v_mag_est);
fprintf('Geschätzter Winkel : %.2f°\n',          v_angle_est);
fprintf('Selbstbeurteilungsmaß GEBA: %.3f\n',    GEBA);

%% Hilfsfunktion: Sequenz per circshift
function seq = generateImageSequence(img, v, M)
    [N, ~] = size(img);
    seq    = zeros(N, N, M);
    mid    = ceil((M+1)/2);
    for k = 1:M
        dx = round((k-mid) * v(1));
        dy = round((k-mid) * v(2));
        seq(:,:,k) = circshift(img, [dy, dx]);
    end
end
