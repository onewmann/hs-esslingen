function out = dilate(mask, radius)
%DILATE  Logical mask grown by a disc of the given radius.
%   OUT = fftvel.internal.dilate(MASK, RADIUS) is true where a true pixel
%   of MASK lies within RADIUS pixels (dy^2 + dx^2 <= RADIUS^2). Nothing is
%   wrapped around the borders. No toolbox is needed (imdilate would need
%   the Image Processing Toolbox). ClipOutliers uses it so that the faint
%   rim of a glint goes with its clipped core. Same as dilate in
%   python/src/fftvel/core.py.

    mask = logical(mask);
    out = mask;
    [h, w] = size(mask);
    for dy = -radius:radius
        for dx = -radius:radius
            if dy*dy + dx*dx > 0 && dy*dy + dx*dx <= radius*radius
                rows = max(dy, 0) + 1:h + min(dy, 0);
                cols = max(dx, 0) + 1:w + min(dx, 0);
                out(rows, cols) = out(rows, cols) | ...
                    mask(max(-dy, 0) + 1:h + min(-dy, 0), max(-dx, 0) + 1:w + min(-dx, 0));
            end
        end
    end
end
