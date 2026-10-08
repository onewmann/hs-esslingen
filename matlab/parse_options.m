function opts = parse_options(args, defaults)
%PARSE_OPTIONS  Merge name-value pairs into a struct of defaults.
%   OPTS = PARSE_OPTIONS(ARGS, DEFAULTS) with ARGS a cell array of
%   'Name', Value pairs. Names are matched case-insensitively against the
%   fields of DEFAULTS; unknown names raise an error.

    opts = defaults;
    if mod(numel(args), 2) ~= 0
        error('parse_options:pairs', 'Options must be given as name-value pairs.');
    end
    names = fieldnames(defaults);
    for k = 1:2:numel(args)
        hit = find(strcmpi(names, args{k}), 1);
        if isempty(hit)
            error('parse_options:unknown', 'Unknown option ''%s''.', args{k});
        end
        opts.(names{hit}) = args{k + 1};
    end
end
