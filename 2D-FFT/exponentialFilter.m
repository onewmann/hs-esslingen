function b = exponentialFilter(A, img)
% exponentialFilter: Filterung mit einem isotropen (radialsymmetrischen)
% Tiefpassfilter (Exponentialfilter)
% Eingaben:
%   img : Bild
%   A : Texturparameter, der die Breite der
%       Autokorrelationsfunktion steuert

kernel_radius = ceil(3 * A);
[X, Y] = meshgrid(-kernel_radius:kernel_radius, -kernel_radius:kernel_radius);
h = exp(-sqrt(X.^2+Y.^2)/A);
h = h / sum(h(:));
b = imfilter(img, h, 'replicate');
