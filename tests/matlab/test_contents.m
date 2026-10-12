function test_contents()
%TEST_CONTENTS  matlab/+fftvel/Contents.m lists every public function and the version.
%   Every function file in matlab/+fftvel has a line 'name - summary' in
%   Contents.m whose summary is the H1 line of the function, and the
%   version in Contents.m is the one in python/pyproject.toml and
%   CITATION.cff.
    here = fileparts(mfilename('fullpath'));
    root = fullfile(here, '..', '..');
    pkg = fullfile(root, 'matlab', '+fftvel');
    txt = read_text(fullfile(pkg, 'Contents.m'));

    tok = regexp(txt, '^%[ \t]+(\w+)[ \t]+-[ \t]+([^\n]*?)[ \t]*$', 'tokens', 'lineanchors');
    listed = cellfun(@(t) t{1}, tok, 'UniformOutput', false);
    summary = cellfun(@(t) t{2}, tok, 'UniformOutput', false);
    files = dir(fullfile(pkg, '*.m'));
    names = regexprep({files.name}, '\.m$', '');
    names = names(~strcmp(names, 'Contents'));
    assert(numel(unique(listed)) == numel(listed), 'Contents.m lists a function twice');
    assert(isequal(sort(listed), sort(names)), sprintf(['Contents.m and the files in ' ...
        'matlab/+fftvel differ: %s'], strjoin(setxor(listed, names), ', ')));
    for k = 1:numel(listed)
        src = read_text(fullfile(pkg, [listed{k} '.m']));
        h1 = regexp(src, ['^%' upper(listed{k}) '[ \t]+([^\n]*?)[ \t]*$'], 'tokens', 'once', ...
            'lineanchors');
        assert(~isempty(h1), sprintf('%s.m has no H1 line', listed{k}));
        assert(strcmp(h1{1}, summary{k}), sprintf(['Contents.m: the summary of %s is ''%s'', ' ...
            'its H1 line says ''%s'''], listed{k}, summary{k}, h1{1}));
    end

    % the same version everywhere
    v = regexp(txt, '^% Version (\S+)[ \t]*$', 'tokens', 'once', 'lineanchors');
    assert(~isempty(v), 'Contents.m has no line ''% Version x.y.z''');
    py = regexp(read_text(fullfile(root, 'python', 'pyproject.toml')), ...
        '^version[ \t]*=[ \t]*"([^"]+)"', 'tokens', 'once', 'lineanchors');
    cff = regexp(read_text(fullfile(root, 'CITATION.cff')), '^version:[ \t]*(\S+)', ...
        'tokens', 'once', 'lineanchors');
    assert(~isempty(py) && ~isempty(cff), 'no version in python/pyproject.toml or CITATION.cff');
    cff = strrep(strrep(cff{1}, '"', ''), '''', '');
    assert(strcmp(v{1}, py{1}), sprintf('Contents.m says version %s, pyproject.toml %s', v{1}, py{1}));
    assert(strcmp(v{1}, cff), sprintf('Contents.m says version %s, CITATION.cff %s', v{1}, cff));
end

function txt = read_text(file)
    fid = fopen(file, 'r');
    assert(fid > 0, sprintf('cannot open %s', file));
    txt = fread(fid, Inf, '*char').';
    fclose(fid);
    txt = strrep(txt, char(13), '');
end
