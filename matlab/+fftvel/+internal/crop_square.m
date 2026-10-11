function [x, r0, c0] = crop_square(img, n)
%CROP_SQUARE  Central N-by-N crop of a 2-D image.
%   X = fftvel.internal.crop_square(IMG) returns the largest central square
%   of IMG.
%   X = fftvel.internal.crop_square(IMG, N) returns the central N-by-N
%   square (N empty: the largest).
%   [X, R0, C0] = fftvel.internal.crop_square(...) also returns the
%   zero-based offsets of the crop. Odd margins put the extra pixel after
%   the crop.

    if ndims(img) ~= 2
        error('crop_square:dims', ['Expected a 2-D grey-scale image, got size %s. ' ...
            'Convert a colour image to grey; for a frame stack use ' ...
            'fftvel.pulsed_from_stack(frames), which averages exactly 2L+1 frames ' ...
            '(getdata returns H x W x 1 x F).'], ...
            mat2str(size(img)));
    end
    [h, w] = size(img);
    if nargin < 2 || isempty(n)
        n = min(h, w);
    end
    if n > min(h, w)
        error('crop_square:size', 'Crop of %d px does not fit a %dx%d image.', n, h, w);
    end
    r0 = floor((h - n)/2);
    c0 = floor((w - n)/2);
    x = img(r0 + (1:n), c0 + (1:n));
end
