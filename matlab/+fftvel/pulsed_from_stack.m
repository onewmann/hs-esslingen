function [bp, K] = pulsed_from_stack(frames, varargin)
%PULSED_FROM_STACK  Pulsed-exposure image from a stack of single exposures.
%   [BP, K] = fftvel.pulsed_from_stack(FRAMES) averages the frames of
%   FRAMES, taken along its last dimension, into one pulsed-exposure image
%   BP (double). All other dimensions of length 1 are dropped, so the
%   H-by-W-by-1-by-F arrays that GETDATA returns work as they are. The
%   number of frames F must then be odd, and K = F.
%
%   [BP, K] = fftvel.pulsed_from_stack(FRAMES, 'Name', Value, ...) sets
%   options:
%     'Exposures'  average the first K frames (K >= 2, also even); wins
%                  over 'L'
%     'L'          average the first 2L+1 frames (L is the half-count of
%                  Horn's thesis)
%
%   Pass K on: EST = fftvel.estimate_velocity(BP, 'Exposures', K).
%   Same as fftvel.pulsed_from_stack in Python.
%
%   Example (Image Acquisition Toolbox):
%     frames = getdata(vid, 9);              % H x W x 1 x 9
%     [bp, K] = fftvel.pulsed_from_stack(frames, 'Exposures', 9);
%     est = fftvel.estimate_velocity(bp, 'Exposures', K);
%
%   See also fftvel.estimate_velocity.

    opts = fftvel.internal.parse_options(varargin, struct('L', [], 'Exposures', []));
    sz = size(frames);
    nf = sz(end);
    keep = sz(1:end-1);
    keep = keep(keep ~= 1);
    if numel(keep) ~= 2
        error('pulsed_from_stack:dims', ...
            'Expected a stack of 2-D frames (H x W x F or H x W x 1 x F), got size %s.', mat2str(sz));
    end
    frames = reshape(frames, [keep nf]);
    if nf < 2
        error('pulsed_from_stack:frames', 'A frame stack needs at least 2 frames, got size %s.', mat2str(sz));
    end
    if ~isempty(opts.L)
        L = fftvel.internal.check_option(opts.L, 'L', 'positive', 'pulsed_from_stack:L');
    end
    if ~isempty(opts.Exposures)
        K = fftvel.internal.check_option(opts.Exposures, 'Exposures', 'exposures', ...
            'pulsed_from_stack:Exposures');
    elseif ~isempty(opts.L)
        K = 2*L + 1;
    elseif mod(nf, 2) == 0
        error('pulsed_from_stack:even', ['The stack has %d frames, an even number; pass ''L'' ' ...
            'or ''Exposures'' to choose how many to average.'], nf);
    else
        K = nf;
    end
    if K > nf
        error('pulsed_from_stack:frames', '%d exposures need %d frames, the stack has %d.', K, K, nf);
    end
    bp = mean(double(frames(:, :, 1:K)), 3);
end
