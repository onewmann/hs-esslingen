function t = make_texture(m, corr_len)
%MAKE_TEXTURE  Random isotropic texture with exponential autocorrelation.
%   T = MAKE_TEXTURE(M, CORR_LEN) returns an M-by-M texture with zero mean,
%   unit standard deviation and autocorrelation close to exp(-r/CORR_LEN)
%   (default CORR_LEN = 1.5 px). White noise is shaped in the Fourier
%   domain with the matching power spectrum
%   (1 + (2*pi*CORR_LEN*f)^2)^(-3/2). The texture is periodic.

    if nargin < 2
        corr_len = 1.5;
    end
    f = [0:ceil(m/2)-1, -floor(m/2):-1]/m;      % same order as fft
    [FX, FY] = meshgrid(f, f);
    psd = (1 + (2*pi*corr_len)^2*(FX.^2 + FY.^2)).^(-1.5);
    t = real(ifft2(fft2(randn(m)) .* sqrt(psd)));
    t = (t - mean(t(:)))/std(t(:));
end
