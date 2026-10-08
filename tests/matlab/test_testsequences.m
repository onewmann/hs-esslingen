function test_testsequences()
%TEST_TESTSEQUENCES  Ground truth from data/testsequences.mat.
%   'model' moves by (5, 5) px/frame and 'tyre' by (4, 4) px/frame (x right,
%   y down), i.e. 135 degrees. Every frame is the previous one shifted with
%   zero fill, so a window of 9 frames is cropped to the region that holds
%   image content in all of them.
    here = fileparts(mfilename('fullpath'));
    d = load(fullfile(here, '..', '..', 'data', 'testsequences.mat'));
    names = {'model', 'tyre'};
    step = [5 4];
    for k = 1:2
        s = double(d.(names{k}))/255;
        for start = [1 6 11]
            last = start + 8;
            b = step(k)*(last - 1);
            bp = mean(s(b+1:end, b+1:end, start:last), 3);
            e = estimate_velocity(bp, 'L', 4);
            assert(abs(angle_error(e.angle_deg, 135)) < 0.1, ...
                sprintf('%s frames %d-%d: angle %.3f', names{k}, start, last, e.angle_deg));
            assert(abs(e.speed/(step(k)*sqrt(2)) - 1) < 0.005, ...
                sprintf('%s frames %d-%d: speed %.4f', names{k}, start, last, e.speed));
            assert(e.valid, sprintf('%s frames %d-%d not valid', names{k}, start, last));
        end
    end
end
