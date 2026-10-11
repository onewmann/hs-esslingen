function test_physical()
%TEST_PHYSICAL  fftvel.to_physical, with the numbers of tests/python/test_physical.py.
    est = struct('speed', 3.6053444133199184, 'vx', 2.997450239305243, ...
        'vy', -2.0034471297107115, 'angle_deg', 33.758125216267615);
    gsd = 0.0123;             % mm per px
    period = 1/750;           % s
    p = fftvel.to_physical(est, gsd, period);
    assert(abs(p.speed/33.25930221287625 - 1) < 1e-12, 'speed');
    assert(abs(p.vx/27.65147845759087 - 1) < 1e-12, 'vx');
    assert(abs(p.vy/-18.481799771581315 - 1) < 1e-12, 'vy');
    assert(p.angle_deg == est.angle_deg, 'angle unchanged');

    % 2.5 px per flash step, 0.04 mm per px, 1 ms between flashes: 100 mm/s
    e = struct('speed', 2.5, 'vx', 0, 'vy', -2.5, 'angle_deg', 90, 'valid', false);
    p = fftvel.to_physical(e, 0.04, 1e-3);
    assert(abs(p.speed - 100) < 1e-12 && abs(p.vy + 100) < 1e-12 && p.vx == 0, '100 mm/s');
    assert(isequal(sort(fieldnames(p)), sort({'speed'; 'vx'; 'vy'; 'angle_deg'})), 'fields');

    % parity cases: estimate as in test_parity.m, convert, compare with the expected speed and angle
    here = fileparts(mfilename('fullpath'));
    pdir = fullfile(here, '..', 'parity');
    d = load(fullfile(pdir, 'cases.mat'));
    ex = dlmread(fullfile(pdir, 'expected.csv'), ',', 1, 0);
    for i = [1 6 10]          % odd n with PadFactor 3, and Duty
        e = fftvel.estimate_velocity(double(d.(sprintf('bp%02d', i))), 'L', ex(i, 2), ...
            'PadFactor', ex(i, 3), 'Duty', ex(i, 10));
        p = fftvel.to_physical(e, gsd, period);
        speed = ex(i, 5)*gsd/period;
        assert(abs(p.speed/speed - 1) < 1e-9, sprintf('case %d: speed', i));
        assert(abs(p.vx - speed*cosd(ex(i, 4))) < 1e-9*speed, sprintf('case %d: vx', i));
        assert(abs(p.vy + speed*sind(ex(i, 4))) < 1e-9*speed, sprintf('case %d: vy', i));
        assert(abs(p.angle_deg - ex(i, 4)) < 1e-9, sprintf('case %d: angle', i));
    end

    % the same checks as in Python
    bad = {0, -0.1, NaN, Inf, '0.1', [], true, [0.1 0.2]};
    for k = 1:numel(bad)
        expect_error(@() fftvel.to_physical(est, bad{k}, period), 'to_physical:gsd');
        expect_error(@() fftvel.to_physical(est, gsd, bad{k}), 'to_physical:period');
    end
    expect_error(@() fftvel.to_physical(struct('speed', 1, 'angle_deg', 0), gsd, period), ...
        'to_physical:est');
    expect_error(@() fftvel.to_physical(2.5, gsd, period), 'to_physical:est');
end

function expect_error(f, id)
    try
        f();
    catch err
        if ~strcmp(err.identifier, id)
            error('test_physical:id', 'expected %s, got %s (%s)', id, err.identifier, err.message);
        end
        return
    end
    error('test_physical:noerror', 'expected the error %s', id);
end
