function p = to_physical(est, gsd, period)
%TO_PHYSICAL  Velocity in length units per second from an estimate in pixels.
%   P = fftvel.to_physical(EST, GSD, PERIOD) converts EST from
%   fftvel.estimate_velocity, whose speed is in pixels per step, into a
%   velocity per second:
%
%     P.speed = EST.speed*GSD/PERIOD
%
%   and P.vx, P.vy likewise; P.angle_deg is EST.angle_deg. GSD is the
%   ground sample distance, the length on the moving surface that one
%   pixel covers (for example mm per px, measured from an image of a ruler
%   or a checkerboard in the plane of the surface). PERIOD is the time
%   step in seconds: the time between two frames for averaged frames, or
%   between two flashes for a multi-flash image. The result is in the length
%   unit of GSD per second (mm/s for GSD in mm per px). Both must be
%   finite and greater than 0.
%
%   EST is converted whether or not it is valid, so check EST.valid first.
%   Any struct with the fields speed, vx, vy and angle_deg works. Same as
%   fftvel.to_physical in Python.
%
%   Example:
%     est = fftvel.estimate_velocity(bp, 'Exposures', 9);
%     p = fftvel.to_physical(est, 0.04, 1e-3);    % 0.04 mm/px, 1 ms between frames
%     fprintf('%.1f mm/s\n', p.speed);
%
%   See also fftvel.estimate_velocity.

    gsd = fftvel.internal.check_option(gsd, 'gsd', 'positive_real', 'to_physical:gsd');
    period = fftvel.internal.check_option(period, 'period', 'positive_real', 'to_physical:period');
    need = {'speed', 'vx', 'vy', 'angle_deg'};
    if ~isstruct(est) || ~all(isfield(est, need))
        error('to_physical:est', ['EST must be a result of fftvel.estimate_velocity ' ...
            '(a struct with the fields speed, vx, vy and angle_deg).']);
    end
    p = struct();
    p.speed = est.speed*gsd/period;
    p.vx = est.vx*gsd/period;
    p.vy = est.vy*gsd/period;
    p.angle_deg = est.angle_deg;
end
