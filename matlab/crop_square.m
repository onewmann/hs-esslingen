function [x, r0, c0] = crop_square(img, n)
%CROP_SQUARE  Central N-by-N crop of a 2-D image.
%   X = CROP_SQUARE(IMG) returns the largest central square of IMG.
%   X = CROP_SQUARE(IMG, N) returns the central N-by-N square.
%   [X, R0, C0] = CROP_SQUARE(...) also returns the zero-based offsets of
%   the crop. Odd margins put the extra pixel after the crop.

    if ndims(img) ~= 2
        error('crop_square:dims', 'Expected a 2-D grey-scale image.');
    end
    [h, w] = size(img);
    if nargin < 2
        n = min(h, w);
    end
    if n > min(h, w)
        error('crop_square:size', 'Crop of %d px does not fit a %dx%d image.', n, h, w);
    end
    r0 = floor((h - n)/2);
    c0 = floor((w - n)/2);
    x = img(r0 + (1:n), c0 + (1:n));
end
