function run_tests()
%RUN_TESTS  Run the MATLAB/Octave test suite.
%   From the repository root:
%     octave --no-gui --eval "cd tests/matlab; run_tests"
%     matlab -batch "cd tests/matlab; run_tests"
%   Throws an error (non-zero exit status) if any test fails.

    here = fileparts(mfilename('fullpath'));
    addpath(fullfile(here, '..', '..', 'matlab'));
    addpath(here);
    tests = {'test_conventions', 'test_options', 'test_known_velocity', ...
             'test_phase_correlation', 'test_testsequences', 'test_parity'};
    failed = {};
    for k = 1:numel(tests)
        t0 = tic;
        try
            feval(tests{k});
            fprintf('PASS  %-24s %6.1f s\n', tests{k}, toc(t0));
        catch err
            fprintf('FAIL  %-24s %s\n', tests{k}, err.message);
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
