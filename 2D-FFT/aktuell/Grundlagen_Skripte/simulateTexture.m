function b = simulateTexture(N, A)
% simulateTexture: Erzeugt ein synthetisches Texturbild b(x)
% Eingaben:
%   N : Bildgröße (N x N)
%   A : Texturparameter, der die Breite der
%       Autokorrelationsfunktion steuert
% Vorgehen: Weißes Rauschen wird mit einem isotropen (radialsymmetrischen)
% Tiefpassfilter (Exponentialfilter) gefiltert.

noise = randn(N,N);

b = exponentialFilter(A,noise);