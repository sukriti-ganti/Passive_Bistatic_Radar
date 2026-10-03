function ok = run_all_tests()
%RUN_ALL_TESTS  Run the eight acceptance tests and summarise.
%
%   OK = RUN_ALL_TESTS() runs test_at1 through test_at8 in order and returns
%   true only if all of them pass.  A failure in one test does not stop the
%   others; the summary at the end says which ones failed.
%
%   AT1 first is not an accident.  A subtle indexing error in the waveform
%   generator would invalidate every result downstream without ever throwing
%   an error, so nothing else is worth reading until AT1 is green.
%
%   Runtime: three to four minutes, most of it AT5 and AT7.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

names = {'test_at1', 'test_at2', 'test_at3', 'test_at4', ...
         'test_at5', 'test_at6', 'test_at7', 'test_at8'};
what  = {'waveform generator', 'ambiguity function', 'bistatic geometry', ...
         'received signal synthesis', 'direct signal cancellation', ...
         'range-Doppler map', 'CFAR detection', 'link budget'};

results = false(1, numel(names));
times   = zeros(1, numel(names));

for i = 1:numel(names)
    f = str2func(names{i});
    t0 = tic;
    try
        results(i) = f();
    catch err
        results(i) = false;
        fprintf('  [ERROR] %s threw: %s\n', names{i}, err.message);
    end
    times(i) = toc(t0);
end

fprintf('\n================= summary =================\n');
for i = 1:numel(names)
    if results(i), tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('  %-9s %-28s %s   %6.1f s\n', names{i}, what{i}, tag, times(i));
end
ok = all(results);
fprintf('  ----------------------------------------\n');
if ok
    fprintf('  all %d acceptance tests pass (%.0f s total)\n', numel(names), sum(times));
else
    fprintf('  %d of %d acceptance tests FAILED\n', sum(~results), numel(names));
end
fprintf('===========================================\n');
end
