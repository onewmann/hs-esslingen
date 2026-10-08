%% Polyphasen-Filterbank (Analyse) – Visualisierung
clear; close all; clc;

%% Parameter
M   = 8;          % Anzahl der Teilbänder (Polyphasen-Zweige)
Lpb = 12;         % Taps pro Band (Gesamtordnung ~ M*Lpb)
Fs  = 1;          % Normierte Abtastrate

%% Prototyp-Tiefpass entwerfen (z.B. mit Kaiser-Fenster)
% Übergangsbereich und Dämpfung grob wählen
Fp  = 1/(2*M);    % Passbandgrenze
Fst = 1/M;        % Stopbandgrenze
Ap  = 0.1;        % Passband-Ripple (dB)
Ast = 80;         % Stopband-Dämpfung (dB)

d = designfilt('lowpassfir', ...
    'PassbandFrequency',Fp, ...
    'StopbandFrequency',Fst, ...
    'PassbandRipple',Ap, ...
    'StopbandAttenuation',Ast, ...
    'DesignMethod','kaiserwin', ...
    'SampleRate',Fs);

h = d.Coefficients(:).';          % Prototyp-Impulsantwort
N = length(h);

% Länge an M anpassen (für saubere Polyphasen-Zerlegung)
Lpad = M*ceil(N/M);
h = [h zeros(1,Lpad-N)];
N  = length(h);

%% Polyphasen-Zerlegung
% E_k(z) = Sum_n h[nM + k] z^{-n}
E = zeros(M, N/M);
for k = 1:M
    E(k,:) = h(k:M:end);
end

%% Modulationsfilterbank (M bandbegrenzte Filter)
% H_k(z) = E(z) * e^{-j*2*pi*k*n/M} (hier über Frequenzgang realisiert)
nfft = 4096;
f    = linspace(0, Fs, nfft);

Hk = zeros(M, nfft);
for k = 0:M-1
    % Frequenzgang des Prototyps
    H = freqz(h, 1, 2*pi*f/Fs);
    % Modulation (Frequenzverschiebung)
    Hk(k+1,:) = H .* exp(-1j*2*pi*k*(0:length(H)-1)/M);
end

%% Visualisierung der Teilband-Frequenzgänge
figure;
plot(f, 20*log10(abs(Hk.')));
grid on;
xlabel('Frequenz (normiert)');
ylabel('Betrag [dB]');
title(sprintf('Frequenzgänge der %d-Teilband-Polyphasen-Filterbank', M));
ylim([-120 10]);

%% Beispielsignal: Mehrton-Signal
Nsig = 4096;
n    = 0:Nsig-1;
x = cos(2*pi*0.05*n) + ...   % Ton 1
    cos(2*pi*0.18*n) + ...   % Ton 2
    cos(2*pi*0.32*n);        % Ton 3

%% Analyse-Filterbank mit Polyphasenstruktur
% 1) Polyphasen-Eingang: in M Unterfolgen zerlegen
xpad = [x zeros(1, M*ceil(length(x)/M) - length(x))];
Xmat = reshape(xpad, M, []).';   % Zeilen: Zeit, Spalten: Polyphasen-Zweige

% 2) Faltung mit Polyphasenfiltern (FIR) und anschließende M‑Punkt‑FFT
Lblk = size(Xmat,1);
y_sub = zeros(Lblk, M);          % Subband-Signale (niedrige Rate)

for nblk = 1:Lblk
    % Polyphasen-Eingangsvektor (M x 1)
    x_vec = Xmat(nblk,:).';
    % Faltung mit Polyphasenfiltern: y_k[n] = sum_l E_k[l] * x[n-l]
    % Hier einfache FIR-Filterung über einen Zustandsvektor
    % (für Demo: direkte Faltung über alle Lpb-Taps)
    if nblk >= size(E,2)
        x_hist = Xmat(nblk-size(E,2)+1:nblk,:).';  % (M x Lpb)
    else
        x_hist = [zeros(M, size(E,2)-nblk), Xmat(1:nblk,:).'].'; %#ok<NASGU>
    end
end

% Für eine kompakte Demo nutzen wir stattdessen direkt die Standard-Implementierung:
% (Polyphasenstruktur ist oben über E sichtbar, hier nur funktionale Demo)
chan = dsp.Channelizer(M);   % Polyphase FFT Analyse-Filterbank
y = chan(x.');               % y: M x Nout

%% Subband-Spektren visualisieren
Y = fftshift(fft(y, nfft, 2), 2);
f_sub = linspace(-0.5, 0.5, nfft);

figure;
imagesc(f_sub, 1:M, 20*log10(abs(Y)));
axis xy;
xlabel('Frequenz (normiert, Subband)');
ylabel('Subband-Index');
title('Spektren der Subband-Signale (Polyphasen-Analysefilterbank)');
colorbar;
caxis([-120 10]);
