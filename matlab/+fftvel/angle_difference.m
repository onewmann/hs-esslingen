function d = angle_difference(a, b)
%ANGLE_DIFFERENCE  Smallest difference of two directions modulo 180 degrees.
%   D = fftvel.angle_difference(A, B) returns A - B in degrees, wrapped to
%   [-90, 90), element-wise:
%
%     D = mod(A - B + 90, 180) - 90
%
%   A direction and its opposite count as the same, as for EST.angle_deg
%   of fftvel.estimate_velocity, so abs(D) is the angle between the two
%   directions, at most 90 degrees; D > 0 when A lies counter-clockwise
%   of B. A and B are arrays of the same size, or one is a scalar. Both
%   are converted to double first, so integer inputs neither saturate nor
%   round, and D is double. Same as fftvel.angle_difference in Python.
%
%   Example:
%     fftvel.angle_difference(179, 1)      % -2
%     fftvel.angle_difference(1, 179)      %  2
%     abs(fftvel.angle_difference(est.angle_deg, 135)) < 0.5
%
%   See also fftvel.estimate_velocity.

    d = mod(double(a) - double(b) + 90, 180) - 90;
end
