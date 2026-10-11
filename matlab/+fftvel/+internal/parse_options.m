function opts = parse_options(args, defaults)
%PARSE_OPTIONS  Merge name-value pairs into a struct of defaults.
%   OPTS = fftvel.internal.parse_options(ARGS, DEFAULTS) with ARGS a cell
%   array of 'Name', Value pairs. Names are matched case-insensitively
%   against the fields of DEFAULTS; unknown names raise an error. The
%   values are not checked here (see fftvel.internal.check_option).

    opts = defaults;
    if mod(numel(args), 2) ~= 0
        error('parse_options:pairs', 'Options must be given as name-value pairs.');
    end
    names = fieldnames(defaults);
    for k = 1:2:numel(args)
        name = args{k};
        if ~((ischar(name) && size(name, 1) == 1) || (isstring(name) && isscalar(name)))
            error('parse_options:name', 'Option name %d is not text (expected one of: %s).', ...
                (k + 1)/2, strjoin(names.', ', '));
        end
        hit = find(strcmpi(names, char(name)), 1);
        if isempty(hit)
            error('parse_options:unknown', 'Unknown option ''%s'' (expected one of: %s).', ...
                char(name), strjoin(names.', ', '));
        end
        opts.(names{hit}) = args{k + 1};
    end
end
