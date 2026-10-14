% FFTVEL  Surface velocity from averaged camera frames (2-D FFT and Radon transform).
% Version 0.1.0
%
% A camera takes K frames at a fixed frame rate while a textured surface
% moves; their mean (pulsed_from_stack) shows the surface K times, shifted
% by equal steps. A strobe flashing K times during one exposure gives the
% same image in a single frame, a multi-flash image, and the same code
% measures it. The method comes from J. Horn's doctoral thesis (2006, see
% README.md), and like the thesis the function names call such an image a
% pulsed-exposure image. The option 'Exposures' gives K: the number of
% averaged frames, or of flashes. 'L' is the half-count of the thesis,
% K = 2L+1.
%
% Estimation
%   pulsed_from_stack     - Pulsed-exposure image from a burst of camera frames: their mean.
%   estimate_velocity     - Velocity of a moving texture from one pulsed-exposure image.
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
%   [~, frames] = fftvel.simulate_pulsed_image(256, [5 -3], 4);   % 9 frames (L = 4, K = 2L+1)
%   [bp, K] = fftvel.pulsed_from_stack(frames);                   % their mean, K = 9
%   est = fftvel.estimate_velocity(bp, 'Exposures', K);
%   fprintf('%.1f deg, %.2f px per frame\n', est.angle_deg, est.speed)
%   % about 31.0 deg, 5.83 px per frame
%
% The live camera script live_basler and the demos are in matlab/examples.
% The functions in fftvel.internal are helpers of the functions above.
% They are not part of the API and can change without notice.
%
% In MATLAB, help fftvel shows this text; in GNU Octave, help fftvel.Contents.
