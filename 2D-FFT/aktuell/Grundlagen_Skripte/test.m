%% Chapter3_4_radon_video_subsample5.m
% Einzelbild-Analyse aus Video, nutzt nur jeden 5. Frame
% • Video-Auswahl ohne harten Abbruch  
% • Subsampling: frameStep = 5  
% • DC-Entfernung, Wiener-Denoising, Kontrast  
% • Hanning-Window, Zero-Padding  
% • 2D-FFT + Radon  
% • subpixelgenaue Peak-Fitting

clear; close all; clc;

%% 0. VIDEO-PFAD VALIDIERUNG UND AUSWAHL
videoPath = '...\data\video.mp4';
if ~isfile(videoPath)
    while true
        [fname,fpath] = uigetfile( ...
          {'*.mp4;*.avi;*.mov','Video-Dateien';'*.*','Alle Dateien'}, ...
          'Wählen Sie eine Video-Datei aus');
        if isequal(fname,0)
            if strcmp( questdlg( ...
                   'Kein Video ausgewählt. Skript beenden?', ...
                   'Abbruch bestätigen','Ja','Nein','Ja'), 'Ja')
                return
            else
                continue
            end
        end
        videoPath = fullfile(fpath,fname);
        break
    end
end

%% 1. PARAMETER
frameStep = 5;           % nur jeder 5. Frame verwenden
theta     = 0:179;       % Radon-Winkel [0°…179°]
N         = 512;         % Quadratgröße vor FFT
padDim    = 2 * N;       % Bildgröße nach Zero-Padding

%% 2. VIDEO EINLESEN + EVERY 5th FRAME SAMPLEN
vidObj      = VideoReader(videoPath);
totalFrames = floor(vidObj.Duration * vidObj.FrameRate);

% Indizes der zu nutzenden Frames
sampleIdx = 1:frameStep:totalFrames;
numFrames = numel(sampleIdx);

% Speicher für subsampled Frames
frames = zeros(N, N, numFrames);

for k = 1:numFrames
    % Springe zu Frame sampleIdx(k)
    vidObj.CurrentTime = (sampleIdx(k)-1) / vidObj.FrameRate;
    fr = readFrame(vidObj);
    
    % Graustufen + Resize
    gray = im2double(fr);
    if size(gray,3)==3
        gray = rgb2gray(gray);
    end
    gray = imresize(gray, [N, N]);
    
    % Kontraststreckung
    lims = stretchlim(gray, [0.01 0.99]);
    gray = imadjust(gray, lims, []);
    
    frames(:,:,k) = gray;
end

fprintf('Total Frames im Video: %d\n', totalFrames);
fprintf('Verwendete Subsampled Frames: %d (jeder %d.)\n', numFrames, frameStep);

% Summenbild
bp = mean(frames, 3);

%% 3. VORVERARBEITUNG
% 3.1 DC-Entfernung
bg    = imopen(bp, strel('disk', round(N/20)));
bp_dc = bp - bg;
bp_dc(bp_dc<0) = 0;

% 3.2 Wiener-Denoising
bp_dn = wiener2(bp_dc, [5 5]);

figure('Name','Summenbild Vorverarbeitet');
subplot(1,2,1), imshow(bp),    title('Roh Summenbild');
subplot(1,2,2), imshow(bp_dn), title('DC entfernt + Wiener');

%% 4. WINDOWING & ZERO-PADDING
h1      = hann(N);
h2d     = h1 * h1.';
bp_win  = bp_dn .* h2d;

pad     = (padDim - N)/2;
bp_pad  = padarray(bp_win, [pad pad], 0, 'both');

%% 5. 2D-FFT & LOG-SPEKTRUM
B    = fftshift(fft2(bp_pad));
logB = log(1 + abs(B));

figure('Name','Log-Spektrum');
imshow(mat2gray(logB));
title('log(1 + |FFT\{Summenbild\}|)');

%% 6. RADON-TRANSFORMATION
[R, xp] = radon(logB, theta);

figure('Name','Radon-Transformierte');
imagesc(theta, xp, R); axis xy;
xlabel('θ (°)'); ylabel('r (Pixel)');
title('Radon(log-Spektrum)');
colorbar;

%% 7. WELLIGKEITSPROFIL WR(θ)
WR = sum(abs(diff(R,1,1)), 1);
[WR_max, idxR] = max(WR);

ridge_angle = theta(idxR);
v_angle_est = mod(ridge_angle, 180);
GEBA        = WR_max / mean(WR);

figure('Name','WR-Profil');
plot(theta, WR, 'LineWidth',2);
xlabel('θ (°)'); ylabel('WR');
title(sprintf('WR-Profil (GEBA = %.3f)', GEBA));
grid on;

%% 8. BETRAGSSCHÄTZUNG MIT SUBPIXEL-FITTING
proj     = R(:, idxR);
mask     = xp >= 0;
xp_pos   = xp(mask);
proj_pos = proj(mask);

[pks, locs] = findpeaks(proj_pos, xp_pos, ...
    'MinPeakProminence',0.05, 'MinPeakDistance',3);

locs_sub = nan(size(locs));
for i = 1:numel(locs)
    [~, ix] = min(abs(xp_pos - locs(i)));
    if ix>1 && ix<length(proj_pos)
        y_m1 = proj_pos(ix-1);
        y_0  = proj_pos(ix);
        y_p1 = proj_pos(ix+1);
        denom  = 2*(y_m1 - 2*y_0 + y_p1);
        offset = (y_m1 - y_p1) / denom;
        locs_sub(i) = locs(i) + offset;
    else
        locs_sub(i) = locs(i);
    end
end

delta_r   = mean(diff(sort(locs_sub)));
v_mag_est = padDim / delta_r;

figure('Name','Radon-Projektion + Subpixel');
plot(xp_pos, proj_pos,'b-','LineWidth',1.5); hold on;
plot(locs,     pks,     'ro','MarkerSize',8);
plot(locs_sub, interp1(xp_pos,proj_pos,locs_sub),...
     'kx','MarkerSize',8);
xlabel('r (Pixel)'); ylabel('Projektion');
title(sprintf('Subpixel-Fit bei %.1f°',ridge_angle));
legend('Projektion','Peaks','Subpixel');

%% 9. ERGEBNISAUSGABE
fprintf('Geschätzter Betrag     : %.2f px/frame\n', v_mag_est);
fprintf('Geschätzter Winkel     : %.2f°\n',          v_angle_est);
fprintf('Selbstbeurteilungsmaß GEBA: %.3f\n',        GEBA);
