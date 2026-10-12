% FFTVEL  Surface velocity from a single multi-flash image (2-D FFT and Radon transform).
% Version 0.1.0
%
% A multi-flash image holds K flashes of a moving surface in one camera
% frame. The method comes from J. Horn's doctoral thesis (2006, see
% README.md), and like the thesis the function names call such an image a
% pulsed-exposure image. The option 'Exposures' gives K: the number of
% flashes, or of the frames that pulsed_from_stack averages. 'L' is the
% half-count of the thesis, K = 2L+1.
%
% Estimation
%   estimate_velocity     - Velocity of a moving texture from one pulsed-exposure image.
%   pulsed_from_stack     - Pulsed-exposure image from a stack of single exposures.
%
% Units and angles
%   to_physical           - Velocity in length units per second from an estimate in pixels.
%   angle_difference      - Smallest difference of two directions modulo 180 degrees.
%
% Simulation
%   simulate_pulsed_image - Pulsed-exposure image of a texture moving under a camera.
%   make_texture          - Random isotropic texture with exponential autocorrelation.
%
% Phase correlation
%   phase_correlation     - Translation between two images.
%
% Plotting
%   plot_estimate         - Show every stage of the estimate in one figure.
%
% Example (from the repository root):
%   addpath('matlab')
%   bp = fftvel.simulate_pulsed_image(256, [5 -3], 4);   % L = 4, so K = 2L+1 = 9 flashes
%   est = fftvel.estimate_velocity(bp, 'Exposures', 9);
%   fprintf('%.1f deg, %.2f px per flash\n', est.angle_deg, est.speed)
%   % about 31.0 deg, 5.83 px per flash
%
% The live camera script live_basler and the demos are in matlab/examples.
% The functions in fftvel.internal are helpers of the functions above.
% They are not part of the API and can change without notice.
%
% In MATLAB, help fftvel shows this text; in GNU Octave, help fftvel.Contents.
