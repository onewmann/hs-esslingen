function s = saturated_share(img)
%SATURATED_SHARE  Share of pixels at the maximum, if that looks like saturation.
%   S = fftvel.internal.saturated_share(IMG) is, for an integer image
%   (int8 ... uint64) whose maximum is INTMAX of its class (255 for uint8,
%   65535 for uint16), the share of pixels at that value, however small.
%   Otherwise (double, single or logical images, and integer images whose
%   maximum lies below INTMAX, such as 12-bit data in uint16) it is the
%   share of pixels at the maximum when that share is at least SAT_SHARE
%   (0.5 %) and the image is not constant, else 0. IMG is used as it is,
%   before any conversion; fftvel.estimate_velocity passes the analysed
%   square. Same as saturated_share in python/src/fftvel/core.py.

    c = fftvel.internal.constants();
    top = max(img(:));
    s = 0;
    if isinteger(img) && top == intmax(class(img))
        s = nnz(img == top)/numel(img);
    elseif top > min(img(:))
        share = nnz(img == top)/numel(img);
        if share >= c.SAT_SHARE
            s = share;
        end
    end
end
