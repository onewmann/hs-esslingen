function d = angle_error(a, b)
%ANGLE_ERROR  Difference of two directions modulo 180 degrees, in [-90, 90).
    d = mod(a - b + 90, 180) - 90;
end
