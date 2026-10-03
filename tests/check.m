function ok = check(name, cond, fmt, varargin)
%CHECK  Print one acceptance test line and return its outcome.
%
%   OK = CHECK(NAME, COND)
%   OK = CHECK(NAME, COND, FMT, ...)
%
%   Prints "  [PASS] name : detail" or "  [FAIL] name : detail" and returns
%   the logical outcome, so a test file can accumulate with
%
%       ok = check(...) && ok;
%
%   The detail string is where the measured number goes.  A test that only
%   prints PASS is much less useful than one that prints PASS and the margin
%   it passed by, especially at a viva.

ok = logical(cond);
if ok
    tag = 'PASS';
else
    tag = 'FAIL';
end
if nargin < 3
    fprintf('  [%s] %s\n', tag, name);
else
    fprintf('  [%s] %s : %s\n', tag, name, sprintf(fmt, varargin{:}));
end
end
