function run_tests()
%RUN_TESTS  Run the MATLAB/Octave test suite.
%   From the repository root:
%     octave --no-gui --eval "cd tests/matlab; run_tests"
%     matlab -batch "cd tests/matlab; run_tests"
%   A failing test is reported with its message and the file and line of
%   every call that led to the error, innermost first. If any test fails,
%   run_tests throws an error at the end, so octave --eval and matlab -batch
%   exit with a non-zero status.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'matlab'));
    addpath(here);
    tests = {'test_conventions', 'test_options', 'test_known_velocity', ...
             'test_model_options', 'test_phase_correlation', 'test_testsequences', ...
             'test_parity'};
    failed = {};
    for k = 1:numel(tests)
        t0 = tic;
        try
            feval(tests{k});
            fprintf('PASS  %-24s %6.1f s\n', tests{k}, toc(t0));
        catch err
            fprintf('FAIL  %-24s %s\n', tests{k}, err.message);
            print_stack(err);
            failed{end + 1} = tests{k}; %#ok<AGROW>
        end
    end
    if isempty(failed)
        fprintf('All %d tests passed.\n', numel(tests));
    else
        error('run_tests:failed', '%d of %d tests failed: %s', numel(failed), ...
            numel(tests), strjoin(failed, ', '));
    end
end

function print_stack(err)
%PRINT_STACK  Function, file and line of each call below run_tests.
    for s = err.stack(:).'
        if strcmp(s.name, 'run_tests')
            break;
        end
        [~, file, ext] = fileparts(s.file);
        fprintf('      at %s (%s%s, line %d)\n', s.name, file, ext, s.line);
    end
end
