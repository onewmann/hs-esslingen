function test_constants()
%TEST_CONSTANTS  fftvel.internal.constants equals tests/constants.json.
%   tests/python/test_constants.py checks the module constants of
%   python/src/fftvel/core.py against the same file, so both
%   implementations use the same numbers.
    here = fileparts(mfilename('fullpath'));
    fid = fopen(fullfile(here, '..', 'constants.json'), 'r');
    assert(fid > 0, 'cannot open tests/constants.json');
    txt = fread(fid, Inf, '*char').';
    fclose(fid);
    spec = jsondecode(txt);
    spec = spec.constants;
    c = fftvel.internal.constants();

    % a new constant must be added to the JSON file and to core.py as well
    names = sort(fieldnames(c));
    expected = sort(fieldnames(spec));
    assert(isequal(names, expected), sprintf('names differ: %s', ...
        strjoin(setxor(names, expected).', ', ')));
    for k = 1:numel(expected)
        name = expected{k};
        value = spec.(name).value;
        assert(isa(c.(name), 'double') && isequal(c.(name)(:), value(:)), ...
            sprintf('%s is %s, tests/constants.json says %s', name, mat2str(c.(name)), ...
            mat2str(value(:).')));
        assert(~isempty(spec.(name).meaning) && ~isempty(spec.(name).used_in), ...
            sprintf('%s: meaning or used_in missing', name));
    end

    % the null-quality offset counts the points taken by the polynomial and by the fit
    assert(c.NULL_OFFSET == c.TREND_DEGREE + 2, 'NULL_OFFSET and TREND_DEGREE');
    assert(mod(180, c.ANGLE_STEP) == 0, 'ANGLE_STEP must divide 180');

    % the struct is built once and then returned as it is
    assert(isequal(fftvel.internal.constants(), c), 'constants changed between calls');
end
