% main_single_image_analysis.m
% Orchestriert Simulation, Summation, Spektralanalyse, Radon-Auswertung,
% Winkel- und Betragsschätzung, Auflösungsreduktion und Gütemaß.
clearvars; close all; clc;

% Parameter (anpassbar)
imgSize = 256;          % quadratische Bildgröße (Pixel)
nFrames = 5;            % Anzahl zu summierender Einzelbilder (2L+1)
frameRate = 1000;       % Hz (nur relativ wichtig für Simulation)
motion_px_per_frame = [31, 10]; % [vx, vy] in Pixel pro Frame (Beispiel)
theta_deg = atan2d(motion_px_per_frame(2), motion_px_per_frame(1));
SNR_db = -6;            % Signal-to-noise ratio für Simulation

% Erzeuge Modellsequenz (textur + Bewegung)
[b_frames, true_v] = sim_texture_and_motion(imgSize, nFrames, motion_px_per_frame, SNR_db);

% Summation (Kapitel 3.1.3)
b_sum = compute_sum_image(b_frames);

% Spektralanalyse (Hanning + FFT + Betrag)
[Bmag, Blog] = compute_magnitude_spectrum(b_sum);

% Radon via Zentralschnitt-Theorem und WR-Berechnung
angles = 0:179; % diskretisierte Winkel Pd
[RadonMat, WR] = radon_via_central_slice(Bmag, angles);

% Winkelabschätzung (Kapitel 3.2.1 / 3.2.3)
beta_est = estimate_angle_from_WR(WR, angles);

% Betragabschätzung (Kapitel 3.2.1 / 3.2.2)
[v_px_per_frame, peakInfo] = estimate_speed_from_projection(RadonMat, beta_est, frameRate);

% Gütemaß (Kapitel 3.4)
GEBA = compute_GEBA(WR);

% Ausgabe / Visualisierung
fprintf('True angle [deg]: %.2f, Estimated beta [deg]: %.2f\n', theta_deg, beta_est);
fprintf('Estimated speed [px/frame]: %.3f\n', v_px_per_frame);
fprintf('GEBA = %.3f\n', GEBA);

% Plots
figure('Name','Inputs and Spectra','NumberTitle','off');
subplot(2,3,1); imagesc(b_frames(:,:,1)); axis image off; title('Frame 1 (example)');
subplot(2,3,2); imagesc(b_sum); axis image off; title('Summenbild b_\Sigma');
subplot(2,3,3); imagesc(Blog); axis image off; title('Log |B_\Sigma|');
subplot(2,3,4); imagesc(RadonMat); colormap hot; axis xy; xlabel('angle index'); ylabel('r index'); title('Radon via ZST');
subplot(2,3,5); plot(angles, WR); xlabel('\phi (deg)'); title('WR(\phi)'); grid on;
subplot(2,3,6); plot(peakInfo.freqs, peakInfo.amps, 'o-'); xlabel('frequency index'); title('Peaks in projection spectrum');

% Experiment: Auflösungsreduktion (Kapitel 3.3)
resolutions = [256 128 64 32];
sigma_beta = zeros(size(resolutions));
for i=1:numel(resolutions)
    ds = imgSize / resolutions(i);
    b_frames_ds = downsample_pyramid(b_frames, ds);
    b_sum_ds = compute_sum_image(b_frames_ds);
    [Bmag_ds, ~] = compute_magnitude_spectrum(b_sum_ds);
    [Rad_ds, WR_ds] = radon_via_central_slice(Bmag_ds, angles);
    beta_ds = estimate_angle_from_WR(WR_ds, angles);
    sigma_beta(i) = wrapTo180(beta_ds - theta_deg); % single-sample diff here for demo
    fprintf('Res %d -> beta_est diff = %.3f deg\n', resolutions(i), sigma_beta(i));
end
