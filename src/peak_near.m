function [lvl, pos, idx] = peak_near(vals, axis_vals, target, halfwidth)
%PEAK_NEAR  Largest value within HALFWIDTH of TARGET on a given axis.
%
%   [LVL, POS, IDX] = PEAK_NEAR(VALS, AXIS_VALS, TARGET, HALFWIDTH)
%
%   Used to read an ambiguity cut at a predicted false peak position: it
%   returns the level found, where exactly it was found, and its index, so
%   that a test can report the offset from the prediction instead of merely
%   asserting that something is there.
%
%   Returns -Inf and NaN if the window contains no samples.

sel = abs(axis_vals(:) - target) <= halfwidth;
v   = vals(:);
v(~sel) = -Inf;
[lvl, idx] = max(v);
if isinf(lvl)
    pos = NaN;
    idx = NaN;
else
    pos = axis_vals(idx);
end
end
