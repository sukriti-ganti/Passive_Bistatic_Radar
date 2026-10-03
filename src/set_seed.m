function set_seed(s)
%SET_SEED  Seed the random generators so every experiment is repeatable.
%
%   SET_SEED(S) uses rng where it exists and falls back to the legacy
%   generator syntax otherwise, so the same scripts run unchanged on MATLAB
%   and on GNU Octave.

try
    rng(s, 'twister');
catch
    rand('twister', s);      %#ok<RAND>
    randn('state', s);       %#ok<RAND>
end
end
