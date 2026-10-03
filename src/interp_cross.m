function x = interp_cross(xv, yv, level)
%INTERP_CROSS  Where a monotonic curve first crosses a level, interpolated.
%
%   X = INTERP_CROSS(XV, YV, LEVEL) finds the first index where YV reaches
%   LEVEL and linearly interpolates between that point and the one before
%   it.  Returns NaN if the curve never reaches the level.
%
%   Used to read "the SNR at which Pd = 0.9" off a measured ROC curve
%   without pretending the sampled grid landed exactly on it.

i = find(yv(:).' >= level, 1);
if isempty(i)
    x = NaN;
elseif i == 1
    x = xv(1);
else
    dy = yv(i) - yv(i - 1);
    if dy == 0
        x = xv(i);
    else
        x = xv(i - 1) + (level - yv(i - 1)) * (xv(i) - xv(i - 1)) / dy;
    end
end
end
