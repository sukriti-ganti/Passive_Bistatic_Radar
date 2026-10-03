function ok = test_at7()
%TEST_AT7  Acceptance test for M7, the two dimensional CFAR detector.
%
%   OK = TEST_AT7() returns true when every check passes.
%
%   Checks, from the PRD
%     noise only input over more than 1e6 cells gives a measured false alarm
%       rate within a factor of two of the design Pfa
%     a strong target is detected
%     the ROC curve rises monotonically with SNR
%
%   The false alarm test uses a synthetic complex Gaussian map rather than a
%   radar map on purpose.  CA-CFAR's threshold factor assumes independent
%   exponentially distributed training cells, and that assumption is what is
%   being tested here.  A real range-Doppler map has correlated neighbouring
%   cells, so its measured rate would be testing the map, not the detector.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT7  CFAR detection ===\n');
ok = true;

% ------------------------------------------------- false alarm rate
fprintf('-- false alarm rate on noise only\n');
set_seed(51);
n = 1024;                                   % 1048576 cells
x = (randn(n) + 1j * randn(n)) / sqrt(2);
p = abs(x).^2;

for pfa = [1e-3 1e-4 1e-5]
    [det_mask, ~, s] = cfar2d(p, struct('pfa', pfa));
    meas = s.n_det / s.n_tested;
    ratio = meas / pfa;
    ok = check(sprintf('measured Pfa within a factor of 2 of %.0e', pfa), ...
        ratio > 0.5 && ratio < 2, ...
        '%d alarms in %d tested cells, measured %.3e, ratio %.2f', ...
        s.n_det, s.n_tested, meas, ratio) && ok;
end

ok = check('more than 1e6 cells were used', numel(p) >= 1e6, ...
    '%d cells, %d of them tested', numel(p), s.n_tested) && ok;

% The threshold factor itself.
np = 416;
alpha_expected = np * (1e-4^(-1 / np) - 1);
[~, ~, s] = cfar2d(p, struct('pfa', 1e-4));
ok = check('alpha = N*(Pfa^(-1/N) - 1) with the default window', ...
    abs(s.alpha - alpha_expected) < 1e-9 && s.n_train == np, ...
    'alpha = %.4f (%.2f dB) with N = %d training cells', ...
    s.alpha, 10 * log10(s.alpha), s.n_train) && ok;

% Guard cells must not be counted as training cells.
[~, ~, s2] = cfar2d(p, struct('guard_range', 3, 'guard_dopp', 1, ...
                              'train_range', 5, 'train_dopp', 4));
expected = (2 * 8 + 1) * (2 * 5 + 1) - (2 * 3 + 1) * (2 * 1 + 1);
ok = check('the guard ring is excluded from the training cells', ...
    s2.n_train == expected, 'N = %d, expected %d', s2.n_train, expected) && ok;

% ------------------------------------------------- strong target
fprintf('-- a strong target\n');
p2 = p;
p2(500, 700) = 200;                          % 23 dB above the mean noise
det_mask = cfar2d(p2, struct('pfa', 1e-4));
ok = check('a strong target is detected', det_mask(500, 700), ...
    'cell under test at %.1f dB above the noise mean', ...
    10 * log10(200 / mean(p(:)))) && ok;

% A target sitting in the guard ring must not raise its own threshold.
p3 = p;
p3(500, 700) = 200;
p3(501, 700) = 60;
det_mask = cfar2d(p3, struct('pfa', 1e-4));
ok = check('a nearby second return does not mask the target', ...
    det_mask(500, 700) && det_mask(501, 700), ...
    'both cells declared') && ok;

% ------------------------------------------------- ROC
fprintf('-- ROC over the full chain\n');
out = roc_curve(struct('nrb', 25, 'nsubframes', 20, ...
                       'snr_db', -55:5:-30, 'ntrials', 60, 'seed', 52));
for i = 1:numel(out.snr_db)
    fprintf('   in-band SNR %5.1f dB, after %.1f dB of processing gain %5.1f dB : Pd = %.3f\n', ...
        out.snr_db(i), out.gp_db, out.snr_post_db(i), out.pd(i));
end

ok = check('Pd rises monotonically with SNR', ...
    all(diff(out.pd) >= -0.05), 'largest backward step %.3f', ...
    max(-diff(out.pd))) && ok;
ok = check('Pd spans the full range from no detection to certain detection', ...
    min(out.pd) < 0.05 && max(out.pd) > 0.95, ...
    'from %.3f to %.3f', min(out.pd), max(out.pd)) && ok;
ok = check('the false alarm rate on the radar map stays near the design value', ...
    out.pfa_meas < 5e-4, 'measured %.2e against a design %.0e', ...
    out.pfa_meas, out.cfg.cfar.pfa) && ok;

% The PRD's link budget assumes 13 dB of post-processing SNR is enough.
% Check that against the measured curve rather than taking it on trust.
i13 = find(out.snr_post_db >= 13, 1);
if ~isempty(i13)
    ok = check('13 dB after processing gives a high probability of detection', ...
        out.pd(i13) > 0.8, 'Pd = %.3f at %.1f dB post-processing SNR', ...
        out.pd(i13), out.snr_post_db(i13)) && ok;
end

fprintf('AT7 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
