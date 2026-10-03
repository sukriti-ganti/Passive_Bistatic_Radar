function ok = test_at5()
%TEST_AT5  Acceptance test for M5, direct signal interference cancellation.
%
%   OK = TEST_AT5() returns true when every check passes.
%
%   Checks, from the PRD
%     cancellation depth of at least 40 dB on the direct path
%     the target peak loses less than 3 dB to the cancellation
%     depth against filter order K improves monotonically then plateaus
%
%   The K sweep is also written to results/cancellation.csv by
%   run_experiment3_cancellation, which plots it.  Two things set the
%   plateau and both are visible in the numbers this test prints:
%
%     the filter has to be long enough to reach the furthest clutter
%       scatterer, K >= clutter delay in samples, or the residual clutter
%       sets the floor
%     the residual still contains the receiver noise, so the depth can
%       never exceed the direct to noise ratio of the input
%
%   The scenario therefore fixes the direct to noise ratio at 60 dB, which
%   is a strong but ordinary direct path, and puts the clutter inside 3 km.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT5  direct signal cancellation ===\n');
ok = true;
C = 299792458;

set_seed(31);
[wave, ~, info] = lte_generate(25, 1.0, 17, 20);

cfg = struct( ...
    'fc', 1800e6, 'tx_pos', [0 0], 'rx_pos', [30e3 0], ...
    'der_db', 60, 'dnr_db', 60, ...
    'n_clutter', 5, 'clutter_max_range', 3000, 'clutter_db', -25, ...
    'targets', struct('pos', {[12e3 9e3]}, 'vel', {[80 -40]}));
[s_surv, s_ref, truth] = synth_rx(wave, info, cfg);

t = truth.targets(1);
clutter_max = max([truth.clutter.delay_samples]);
k_needed = ceil(cfg.clutter_max_range / C * info.fs);

fprintf('-- scenario: DER %.0f dB, DNR %.1f dB, target at %.0f m / %.1f Hz\n', ...
    cfg.der_db, truth.dnr_db, t.Rb, t.fd);
fprintf('   clutter out to %.0f samples, so K must reach %d\n', clutter_max, k_needed);

% ------------------------------------------------- depth against K
k_list = [0 8 16 24 32 48 64 80 96 112 128 160 192];
depth  = zeros(size(k_list));
loss   = zeros(size(k_list));
for i = 1:numel(k_list)
    [s_clean, ~, ci] = eca(s_surv, s_ref, k_list(i));
    a = match_amp(s_clean, s_ref, info.fs, t.delay_samples, t.fd);
    depth(i) = ci.depth_db;
    loss(i)  = 20 * log10(t.amp / abs(a));
    fprintf('   K = %3d : depth %6.2f dB, target loss %5.2f dB\n', ...
        k_list(i), depth(i), loss(i));
end

k_full = find(k_list >= k_needed, 1);

ok = check('cancellation depth reaches at least 40 dB', ...
    max(depth) >= 40, 'best %.2f dB at K = %d, ceiling set by DNR is %.1f dB', ...
    max(depth), k_list(find(depth == max(depth), 1)), truth.dnr_db) && ok;

ok = check('target loss stays under 3 dB once K covers the clutter', ...
    all(loss(k_full:end) < 3), 'worst %.2f dB for K >= %d', ...
    max(loss(k_full:end)), k_list(k_full)) && ok;

% Monotonic improvement, allowing a small numerical wobble.
ok = check('depth improves monotonically with K', ...
    all(diff(depth) > -0.5), 'largest backward step %.3f dB', ...
    -min(diff(depth))) && ok;

ok = check('and then plateaus', abs(depth(end) - depth(end - 1)) < 1, ...
    'last two points differ by %.3f dB', abs(depth(end) - depth(end - 1))) && ok;

ok = check('a filter too short to reach the clutter is much worse', ...
    depth(1) < max(depth) - 20, ...
    'K = 0 gives %.2f dB against %.2f dB at K = %d', ...
    depth(1), max(depth), k_list(find(depth == max(depth), 1))) && ok;

% ------------------------------------------------- the target becomes visible
kbest = k_list(find(depth == max(depth), 1));
[s_clean, ~, ci] = eca(s_surv, s_ref, kbest);
nb   = info.samples_per_subframe;
mlag = 400;
[chi_b, tau_b, fd_b] = caf(s_surv,  s_ref, info.fs, nb, mlag);
[chi_a, tau_a, fd_a] = caf(s_clean, s_ref, info.fs, nb, mlag);

[~, rb] = min(abs(tau_b * C - t.Rb));
[~, cb] = min(abs(fd_b - t.fd));
snr_before = 10 * log10(abs(chi_b(rb, cb))^2 / median(abs(chi_b(:)).^2));
snr_after  = 10 * log10(abs(chi_a(rb, cb))^2 / median(abs(chi_a(:)).^2));
ok = check('the target cell rises above the map floor after cancellation', ...
    snr_after > snr_before + 20, ...
    '%.1f dB above the map median before, %.1f dB after', snr_before, snr_after) && ok;

% The strongest cell on the cancelled map should now be the target, not the
% direct path.
p = abs(chi_a).^2;
[~, k] = max(p(:));
[r, c] = ind2sub(size(p), k);
ok = check('the strongest remaining peak is the target', ...
    abs(tau_a(r) * C - t.Rb) < 2 * C / info.fs && abs(fd_a(c) - t.fd) < 2 * (fd_a(2) - fd_a(1)), ...
    'peak at %.0f m / %.0f Hz, target at %.0f m / %.1f Hz', ...
    tau_a(r) * C, fd_a(c), t.Rb, t.fd) && ok;

% ------------------------------------------------- NLMS comparison
fprintf('-- NLMS, the simpler alternative\n');
for mu = [0.003 0.01 0.03]
    [s_n, ~, ni] = nlms(s_surv, s_ref, kbest, mu);
    a = match_amp(s_n, s_ref, info.fs, t.delay_samples, t.fd);
    fprintf('   mu = %5.3f : depth %6.2f dB, target loss %6.2f dB, notch about %.0f Hz wide\n', ...
        mu, ni.depth_db, 20 * log10(t.amp / abs(a)), mu * info.fs / (2 * pi * (kbest + 1)));
end
[s_n, ~, ni] = nlms(s_surv, s_ref, kbest, 0.01);
ok = check('ECA beats NLMS at the same filter order', ...
    ci.depth_db > ni.depth_db + 10, ...
    'ECA %.1f dB against NLMS %.1f dB at mu = 0.01', ci.depth_db, ni.depth_db) && ok;
ok = check('ECA drives the residual correlation with the reference to zero', ...
    ci.dsi_residual_db < -100 && ni.dsi_residual_db > ci.dsi_residual_db, ...
    'ECA %.0f dB, NLMS %.0f dB', ci.dsi_residual_db, ni.dsi_residual_db) && ok;

fprintf('AT5 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
