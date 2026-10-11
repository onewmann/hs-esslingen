function x = check_option(x, name, kind, id)
%CHECK_OPTION  Check the value of an option and return it as double.
%   X = fftvel.internal.check_option(X, NAME, KIND, ID) returns X as a
%   double if it is a real, finite numeric scalar (not text, not logical)
%   of the given KIND, and raises the error with identifier ID otherwise.
%   NAME is the option name used in the message. KIND is one of
%     'positive'   integer >= 1
%     'exposures'  integer >= 2
%     'roi'        integer >= 16 (pixels)
%     'unit'       real number in [0, 1]
%     'duty'       real number in [0, 1)
%     'flag'       true or false (logical, or numeric 0 or 1); X is
%                  returned as logical
%   Text such as '4' is rejected, because double('4') would be 52.
%   Same rules and messages as _check_option in python/src/fftvel/core.py.

    switch kind
        case 'positive'
            what = 'a positive integer';
        case 'exposures'
            what = 'an integer of at least 2';
        case 'roi'
            what = 'an integer of at least 16 (px)';
        case 'unit'
            what = 'a finite real number in [0, 1]';
        case 'duty'
            what = 'a finite real number in [0, 1)';
        case 'flag'
            what = 'true or false';
        otherwise
            error('check_option:kind', 'Unknown kind ''%s''.', kind);
    end

    if ischar(x) || isstring(x)
        c = char(x);
        if size(c, 1) == 1
            got = sprintf('text ''%s''', c);
        else
            got = 'text';
        end
    elseif ~isscalar(x)
        got = sprintf('size %s', mat2str(size(x)));
    elseif islogical(x)
        if strcmp(kind, 'flag')
            return
        end
        got = 'a logical value';
    elseif isnumeric(x) && ~isreal(x)
        got = 'a complex value';
    elseif ~isnumeric(x)
        got = sprintf('a value of class %s', class(x));
    else
        v = double(x);
        switch kind
            case 'flag'
                ok = v == 0 || v == 1;
            case 'unit'
                ok = v >= 0 && v <= 1;
            case 'duty'
                ok = v >= 0 && v < 1;
            otherwise
                lowest = struct('positive', 1, 'exposures', 2, 'roi', 16);
                ok = isfinite(v) && v == floor(v) && v >= lowest.(kind);
        end
        if ok
            if strcmp(kind, 'flag')
                x = v == 1;
            else
                x = v;
            end
            return
        end
        got = sprintf('%g', v);
    end
    error(id, '%s must be %s, got %s.', name, what, got);
end
