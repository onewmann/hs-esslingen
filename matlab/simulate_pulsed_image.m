function [bp, frames] = simulate_pulsed_image(n, v, L, varargin)
%SIMULATE_PULSED_IMAGE  Pulsed-exposure image of a texture moving under a camera.
%   BP = SIMULATE_PULSED_IMAGE(N, V, L) returns the N-by-N mean of 2L+1
%   exposures k = -L..L of a random texture shifted by k*V, V = [vx vy] in
%   pixels per step (x along columns, y along rows pointing down). The
%   camera window stays fixed while the texture moves underneath, so
%   content enters and leaves the image as it would in a real recording.
%   There are no zero-filled borders.
%
%   Shifts are exact band-limited (Fourier) shifts of a texture that is
%   large enough that no wrapped content reaches the window, so V may be
%   fractional.
%
%   Options (name-value pairs):
%     'CorrLength'  correlation length of the random texture (default 1.5)
%     'Texture'     use this image instead of a random texture; it must be
%                   at least N + 2*L*max(|V|) + 32 pixels in both directions
%     'Noise'       standard deviation of white noise added to BP, relative
%                   to the standard deviation of the noise-free BP
%                   (default 0)
%
%   [BP, FRAMES] = SIMULATE_PULSED_IMAGE(...) also returns the individual
%   exposures as an N-by-N-by-(2L+1) array.

    opts = parse_options(varargin, struct('CorrLength', 1.5, ...
        'Texture', [], 'Noise', 0));
    vx = v(1);
    vy = v(2);
    margin = 2*L*max(abs([vx vy])) + 32;
    if isempty(opts.Texture)
        m = 2^ceil(log2(n + margin));
        tex = make_texture(m, opts.CorrLength);
    else
        tex = double(opts.Texture);
        if size(tex, 1) < n + margin || size(tex, 2) < n + margin
            error('simulate_pulsed_image:texture', ...
                'Texture too small for N=%d and |v|=%g with L=%d.', n, max(abs([vx vy])), L);
        end
        side = min(size(tex));
        tex = crop_square(tex, side - mod(side, 2));
        tex = tex - mean(tex(:));
    end
    m = size(tex, 1);
    f = [0:ceil(m/2)-1, -floor(m/2):-1]/m;
    [FX, FY] = meshgrid(f, f);
    T = fft2(tex);
    r0 = floor((m - n)/2);
    win = r0 + (1:n);

    want_frames = nargout > 1;
    if want_frames
        frames = zeros(n, n, 2*L + 1);
    end
    D = zeros(m);
    for k = -L:L
        e = exp(-2i*pi*(FX*vx + FY*vy)*k);
        D = D + e;
        if want_frames
            fr = real(ifft2(T .* e));
            frames(:, :, k + L + 1) = fr(win, win);
        end
    end
    big = real(ifft2(T .* D/(2*L + 1)));
    bp = big(win, win);
    if opts.Noise > 0
        bp = bp + opts.Noise*std(bp(:))*randn(n);
    end
end
