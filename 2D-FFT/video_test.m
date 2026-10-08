%% Chapter3_4_GEBA.m
% Selbstbeurteilung der Messqualität (GEBA) in Einzelbild-Analyse
% - Erzeugt eine Textur, simuliert eine Bildfolge bei Translation v
% - Bildfolge wird gemittelt, Spektrum (log) berechnet, Radon angewendet
% - WR(θ) = sum(|d/d r Radon(logB,θ)|) berechnet, GEBA = max(WR)/mean(WR)
% - Zusätzlich: robuste Verschiebungsschätzung per Phase-Correlation als Vergleich
% Hinweise: Testen Sie zuerst mit moderatem N (z.B. N=256) wegen RAM/CPU.

clear; close all; clc;

%% ------------------------
%% 1. Parameter (anpassen)
%% ------------------------
N = 512;                 % Bildgröße N x N
A = 0.3;                 % Texturparameter (Spektralbreite)
v_true = [-13, -13];     % Wahrer Geschwindigkeitsvektor [vx, vy] in Pixel/frame
L = 80;                  % Halbanzahl der Frames: numFrames = 2*L + 1
dt = 1;                  % Zeit zwischen Frames (frames)
numFrames = 2*L + 1;

%% ------------------------
%% 2. Erzeuge Textur & Sequenz
%% ------------------------
b = simulateTexture(N, A);                       % Basis-Texturbild (normiert)
sequence = generateImageSequence(b, v_true, numFrames); % N x N x numFrames
bp = mean(sequence, 3);                           % Summiertes (Bewegungsunschärfe) Bild

figure('Name','Summiertes Bild'); imshow(mat2gray(bp)); title('bp (gemittelt)');

%% ------------------------
%% 3. Spektrum (FFT) und Log-Darstellung
%% ------------------------
B = fftshift(fft2(bp));
logB = log(1 + abs(B));

figure('Name','Log-Spektrum'); imshow(mat2gray(logB)); title('Log |B|');

%% ------------------------
%% 4. Radon-Transformation & WR-Berechnung
%% ------------------------
theta = 0:1:179;
[R, xp] = radon(logB, theta);   % Radon-Projektionen

% Visualisierung Radon
figure('Name','Radon map'); imagesc(theta, xp, R); xlabel('θ (deg)'); ylabel('r'); title('Radon(logB)'); colorbar;

% WR(θ): summierter absoluter Gradient entlang r (Linienstärke)
numAngles = numel(theta);
WR = zeros(1, numAngles);
for i = 1:numAngles
    proj = R(:, i);
    gradProj = abs(diff(proj));    % |d/d r proj|
    WR(i) = sum(gradProj);
end

% Bestimmung maximaler Richtung & GEBA
[WR_max, bestIdx] = max(WR);
v_angle_radon = theta(bestIdx);      % Winkel aus Radon (deg)
GEBA = WR_max / mean(WR);

% Visualisierung ausgewählte Projektion
selectedProj = R(:, bestIdx);
figure('Name','Best Radon Projection'); plot(selectedProj, 'LineWidth', 2); xlabel('r index'); ylabel('Projection'); title(sprintf('Radon at %.1f°', v_angle_radon)); grid on;

% Peaks & spacing (heuristisch)
sel_sm = movmedian(selectedProj, 5);
[pks, locs] = findpeaks(sel_sm, 'MinPeakProminence', 0.02, 'MinPeakDistance', 3);
if numel(locs) > 1
    spacings = diff(locs);
    spacing_avg = mean(spacings);
else
    spacing_avg = NaN;
end
v_mag_est_radon = spacing_avg;  % nur heuristisch (Index-Abstand)

%% ------------------------
%% 5. Robuste Schätzung: Phase Correlation (empfohlen)
%% ------------------------
mid = ceil((numFrames + 1)/2);
I1 = sequence(:,:,mid);
I2 = sequence(:,:,mid+1);   % ein Frame später
shift = estimateShiftPhaseCorr(I1, I2);  % [dx, dy] I2 relativ zu I1
v_est_vec = shift / dt;
v_mag_est_pc = norm(v_est_vec);
% Winkel: Bild y-Achse nach unten -> invertieren für Anzeige
v_angle_est_pc = atan2d(-v_est_vec(2), v_est_vec(1));

% Wähle Methode: Phase-Corr falls gültig, sonst Radon-Fallback
if all(isfinite(v_est_vec)) && any(abs(v_est_vec) > 0)
    v_mag_est = v_mag_est_pc;
    v_angle_est = v_angle_est_pc;
    method = 'Phase Correlation';
else
    v_mag_est = v_mag_est_radon;
    v_angle_est = v_angle_radon;
    method = 'Radon spacing fallback';
end

%% ------------------------
%% 6. Ausgabe & Visualisierung
%% ------------------------
v_true_mag = norm(v_true);
fprintf('Methode: %s\n', method);
fprintf('Wahrer Betrag (px/frame): %.2f\n', v_true_mag);
fprintf('Geschätzter Betrag (px/frame): %.2f\n', v_mag_est);
fprintf('Geschätzter Winkel (deg): %.2f\n', v_angle_est);
fprintf('GEBA = max(WR)/mean(WR): %.3f\n', GEBA);

% WR Profil
figure('Name','WR profile'); plot(theta, WR, 'LineWidth', 2); xlabel('θ (deg)'); ylabel('WR'); title(sprintf('WR (GEBA=%.3f)', GEBA)); grid on;

% Beispiel-Frames mit Richtungs-Vektor
figure('Name','Frames and motion');
subplot(1,2,1); imagesc(I1); colormap gray; axis image off; title('Frame mid');
hold on; c = size(I1); ctr = [c(2)/2, c(1)/2]; plot(ctr(1), ctr(2), 'ro');
len = 40; dx = len * cosd(v_angle_est); dy = -len * sind(v_angle_est);
quiver(ctr(1), ctr(2), dx, dy, 0, 'r', 'LineWidth', 1.5, 'MaxHeadSize', 2); hold off;
subplot(1,2,2); imagesc(I2); colormap gray; axis image off; title('Frame mid+1');

%% ------------------------
%% 7. Hilfsfunktionen (unten einfügen)
%% ------------------------
% Füge die folgenden Funktionen ans Ende derselben Datei oder als separate
% Dateien mit denselben Funktionsnamen ein: simulateTexture, generateImageSequence,
% estimateShiftPhaseCorr. Die Implementationen siehe unten.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IMPLEMENTATION HILFSFUNKTIONEN (kopiere in dieselbe Datei unterhalb
% oder speichere als separate .m-Dateien: simulateTexture.m, generateImageSequence.m,
% estimateShiftPhaseCorr.m)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function b = simulateTexture(N, A)
% Erzeugt isotrope Textur über radiale Amplitude + Zufallsphase.
% Ausgabe normiert auf [-1,1].
    if nargin < 2, A = 0.3; end
    [u, v] = meshgrid((-floor(N/2)):(ceil(N/2)-1));
    rho = sqrt(u.^2 + v.^2);
    rho = rho / max(rho(:));
    sigma = A;
    amp = exp(- (rho.^2) / (2 * sigma^2));
    amp(1,1) = 0;
    rng(0);
    phase = exp(1i * 2 * pi * rand(N));
    S = fftshift(amp) .* phase;
    b_complex = ifft2(ifftshift(S)) * N^2;
    b = real(b_complex);
    b = b - mean(b(:));
    b = b / max(abs(b(:)));
end

function seq = generateImageSequence(b, v, numFrames)
% Erzeugt Sequenz durch Translation: b_k(x) = b(x - v*(k - k0))
% Verwende imtranslate wenn verfügbar, sonst circshift-Fallback.
    k0 = ceil((numFrames + 1)/2);
    [N, M] = size(b);
    seq = zeros(N, M, numFrames, 'like', b);
    hasImtranslate = exist('imtranslate', 'file') == 2;
    for k = 1:numFrames
        shift = v * (k - k0);
        if hasImtranslate
            frame = imtranslate(b, [shift(1), shift(2)], 'linear', 'FillValues', 0, 'OutputView', 'same');
        else
            dx = round(shift(1)); dy = round(shift(2));
            frame = circshift(b, [dy, dx]);
            if dx > 0, frame(:,1:dx)=0; elseif dx<0, frame(:,end+dx+1:end)=0; end
            if dy > 0, frame(1:dy,:)=0; elseif dy<0, frame(end+dy+1:end,:)=0; end
        end
        seq(:,:,k) = frame;
    end
end

function shift = estimateShiftPhaseCorr(I1, I2)
% Subpixelverschiebung per Phase Correlation
% Rückgabe shift = [dx, dy] : Verschiebung I2 relativ zu I1 (px)
    I1 = double(I1); I2 = double(I2);
    win = hann(size(I1,1)) * hann(size(I1,2))';
    F1 = fft2(I1 .* win);
    F2 = fft2(I2 .* win);
    R = F1 .* conj(F2);
    R = R ./ (abs(R) + eps);
    r = ifft2(R);
    r = fftshift(real(r));
    [~, idx] = max(r(:));
    [py, px] = ind2sub(size(r), idx);
    center = floor(size(r)/2) + 1;
    dy = py - center(1);
    dx = px - center(2);
    % Subpixelkorrektur 3x3 Quadratic Fit falls möglich
    patchRadius = 1;
    yRange = max(py-patchRadius,1):min(py+patchRadius,size(r,1));
    xRange = max(px-patchRadius,1):min(px+patchRadius,size(r,2));
    P = r(yRange, xRange);
    if all(size(P) == [3,3])
        Z = log(abs(P(:)) + eps);
        [Xg, Yg] = meshgrid(-1:1, -1:1);
        A = [Xg(:).^2, Yg(:).^2, Xg(:).*Yg(:), Xg(:), Yg(:), ones(9,1)];
        coeff = A \ Z;
        a = coeff(1); b = coeff(2); e = coeff(4); f = coeff(5);
        if abs(2*a) > eps, dx_sub = -e/(2*a); else dx_sub = 0; end
        if abs(2*b) > eps, dy_sub = -f/(2*b); else dy_sub = 0; end
        if isfinite(dx_sub), dx = dx + dx_sub; end
        if isfinite(dy_sub), dy = dy + dy_sub; end
    end
    shift = [dx, dy];
end
