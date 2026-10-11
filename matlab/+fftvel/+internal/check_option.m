function x = check_option(x, name, kind, id)
%CHECK_OPTION  Check the value of an option and return it as double.
%   X = fftvel.internal.check_option(X, NAME, KIND, ID) returns X as a
%   double if it is a real, finite numeric scalar (not text, not logical)
%   of the given KIND, and raises the error with identifier ID otherwise.
%   NAME is the option name used in the message. KIND is one of
%     'positive'       integer >= 1
%     'exposures'      integer >= MIN_EXPOSURES (2)
%     'roi'            integer >= MIN_SIZE (16 pixels)
%     'unit'           real number in [0, 1]
%     'duty'           real number in [0, 1)
%     'positive_real'  finite real number > 0
%     'flag'           true or false (logical, or numeric 0 or 1); X is
%                      returned as logical
%   MIN_EXPOSURES and MIN_SIZE come from fftvel.internal.constants.
%   Text such as '4' is rejected, because double('4') would be 52.
%   Same rules and messages as _check_option in python/src/fftvel/core.py.

    lowest = [];
    switch kind
        case 'positive'
            what = 'a positive integer';
            lowest = 1;
        case 'exposures'
            cst = fftvel.internal.constants();
            lowest = cst.MIN_EXPOSURES;
            what = sprintf('an integer of at least %d', lowest);
        case 'roi'
            cst = fftvel.internal.constants();
            lowest = cst.MIN_SIZE;
            what = sprintf('an integer of at least %d (px)', lowest);
        case 'unit'
            what = 'a finite real number in [0, 1]';
        case 'duty'
            what = 'a finite real number in [0, 1)';
        case 'positive_real'
            what = 'a finite real number > 0';
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
            case 'positive_real'
                ok = isfinite(v) && v > 0;
            otherwise
                ok = isfinite(v) && v == floor(v) && v >= lowest;
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
