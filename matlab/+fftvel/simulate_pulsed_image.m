function [bp, frames] = simulate_pulsed_image(n, v, L, varargin)
%SIMULATE_PULSED_IMAGE  Pulsed-exposure image of a texture moving under a camera.
%   BP = fftvel.simulate_pulsed_image(N, V, L) returns the N-by-N mean of
%   2L+1 exposures k = -L..L of a random texture shifted by k*V,
%   V = [vx vy] in pixels per step (x along columns, y along rows pointing
%   down). L is the half-count of Horn's thesis: L = 4 gives K = 2L+1 = 9
%   flashes. The camera window stays fixed while the texture moves
%   underneath, so content enters and leaves the image as it would in a
%   real recording. There are no zero-filled borders.
%
%   Shifts are exact band-limited (Fourier) shifts of a texture that is
%   large enough that no wrapped content reaches the window, so V may be
%   fractional.
%
%   Options (name-value pairs):
%     'CorrLength'  correlation length of the random texture (default 1.5)
%     'Texture'     use this image instead of a random texture; it must be
%                   at least N + (K-1+Duty)*max(|V|) + 32 pixels in both
%                   directions. An RGB(A) image is converted to grey with
%                   the weights 0.2989, 0.5870, 0.1140.
%     'Noise'       standard deviation of white noise added to BP, relative
%                   to the standard deviation of the noise-free BP
%                   (default 0)
%     'Exposures'   number of exposures K (default [] = 2L+1), overrides
%                   L; an even K puts the exposures at half-integer steps
%                   around the centre
%     'Duty'        duration of each exposure divided by the time step, in
%                   [0, 1): the texture also moves during each exposure
%                   (box blur of Duty*|V| pixels along V). Default 0:
%                   instantaneous exposures.
%
%   [BP, FRAMES] = fftvel.simulate_pulsed_image(...) also returns the
%   individual exposures as an N-by-N-by-K array.
%
%   Example:
%     % L = 4, so K = 2L+1 = 9 flashes
%     [bp, frames] = fftvel.simulate_pulsed_image(256, [5 -3], 4);
%     est = fftvel.estimate_velocity(bp, 'Exposures', 9, 'Diagnostics', true);
%     fftvel.plot_estimate(bp, est, frames(:, :, 5));   % the middle exposure, k = 0

    c = fftvel.internal.constants();
    opts = fftvel.internal.parse_options(varargin, struct('CorrLength', c.DEFAULT_CORR_LEN, ...
        'Texture', [], 'Noise', 0, 'Exposures', [], 'Duty', 0));
    L = fftvel.internal.check_option(L, 'L', 'positive', 'simulate_pulsed_image:L');
    if isempty(opts.Exposures)
        K = 2*L + 1;
    else
        K = fftvel.internal.check_option(opts.Exposures, 'Exposures', 'exposures', ...
            'simulate_pulsed_image:Exposures');
    end
    duty = fftvel.internal.check_option(opts.Duty, 'Duty', 'duty', 'simulate_pulsed_image:Duty');
    % same rule as _velocity in python/src/fftvel/simulate.py
    if ~isnumeric(v) || ~isreal(v) || numel(v) ~= 2 || ~all(isfinite(v(:)))
        error('simulate_pulsed_image:v', 'v must hold two finite real numbers [vx vy].');
    end
    vx = double(v(1));
    vy = double(v(2));
    margin = ((K - 1) + duty)*max(abs([vx vy])) + c.SIM_MARGIN;
    if isempty(opts.Texture)
        m = 2^ceil(log2(n + margin));
        tex = fftvel.make_texture(m, opts.CorrLength);
    else
        tex = double(opts.Texture);
        if ndims(tex) == 3 && (size(tex, 3) == 3 || size(tex, 3) == 4)
            % RGB(A) -> grey, rgb2gray weights
            gw = c.GREY_WEIGHTS;
            tex = gw(1)*tex(:, :, 1) + gw(2)*tex(:, :, 2) + gw(3)*tex(:, :, 3);
        end
        if ndims(tex) ~= 2
            error('simulate_pulsed_image:texture', ...
                'Texture must be a 2-D grey-scale or an RGB(A) image, got size %s.', mat2str(size(tex)));
        end
        if size(tex, 1) < n + margin || size(tex, 2) < n + margin
            error('simulate_pulsed_image:texture', ...
                'Texture too small for N=%d and |v|=%g with %d exposures.', n, max(abs([vx vy])), K);
        end
        side = min(size(tex));
        tex = fftvel.internal.crop_square(tex, side - mod(side, 2));
        tex = tex - mean(tex(:));
    end
    m = size(tex, 1);
    f = [0:ceil(m/2)-1, -floor(m/2):-1]/m;
    [FX, FY] = meshgrid(f, f);
    T = fft2(tex);
    if duty > 0
        % motion during each exposure: box blur along v (sinc without toolbox)
        x = duty*(FX*vx + FY*vy);
        B = ones(size(x));
        nz = x ~= 0;
        B(nz) = sin(pi*x(nz))./(pi*x(nz));
        T = T .* B;
    end
    r0 = floor((m - n)/2);
    win = r0 + (1:n);

    want_frames = nargout > 1;
    if want_frames
        frames = zeros(n, n, K);
    end
    D = zeros(m);
    for q = 1:K
        k = (q - 1) - (K - 1)/2;       % -L..L for K = 2L+1
        e = exp(-2i*pi*(FX*vx + FY*vy)*k);
        D = D + e;
        if want_frames
            fr = real(ifft2(T .* e));
            frames(:, :, q) = fr(win, win);
        end
    end
    big = real(ifft2(T .* D/K));
    bp = big(win, win);
    if opts.Noise > 0
        bp = bp + opts.Noise*std(bp(:))*randn(n);
    end
end
