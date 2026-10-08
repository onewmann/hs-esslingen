function wr = ripple_profile(M, j, rfrac)
%RIPPLE_PROFILE  Total variation of each projection, WR(theta).
%   WR = RIPPLE_PROFILE(M, J, RFRAC) returns sum(|diff(M)|) over the bins
%   with |J| <= RFRAC*max(J), one value per column of M. The motion lines
%   in the spectrum turn into a sharp comb only in the projection along v,
%   so WR peaks at the direction of motion.

    if nargin < 3
        rfrac = 0.9;
    end
    sel = abs(j) <= rfrac*max(j);
    wr = sum(abs(diff(M(sel, :), 1, 1)), 1);
end
