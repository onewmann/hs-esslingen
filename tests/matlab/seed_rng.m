function seed_rng(s)
%SEED_RNG  Seed the random generator in MATLAB and in old and new Octave.
    if exist('rng', 'file') || exist('rng', 'builtin')
        rng(s);
    else
        randn('state', s); %#ok<RAND>
        rand('state', s); %#ok<RAND>
    end
end
