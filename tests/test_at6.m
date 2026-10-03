function ok = test_at6()
%TEST_AT6  Acceptance test for M6, the range-Doppler map.
%
%   OK = TEST_AT6() returns true when every check passes.
%
%   M6 is CAF applied to the cancelled surveillance signal instead of to the
%   reference, so this test exercises the whole chain M1 -> M4 -> M5 -> M6.
%
%   Checks, from the PRD
%     a target injected in M4 appears within one resolution cell of its
%       predicted position
%     two targets separated by more than dR give two distinct peaks
%     two targets separated by less than dR give one merged peak
%
%   The separation is applied along the bistatic bisector, because that is
%   the direction the resolution formula dR = c/(2*B*cos(beta/2)) refers to:
%   a ground displacement s along the bisector moves the bistatic range by
%   2*s*cos(beta/2).  Displacing a target along the iso-range ellipse
%   instead would not change its bistatic range at all, however far it went.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT6  range-Doppler map ===\n');
ok = true;
C = 299792458;

set_seed(41);
[wave, ~, info] = lte_generate(50, 1.0, 17, 20);
nb   = info.samples_per_subframe;
mlag = 600;
b    = info.occupied_bw;

tx = [0 0];
rx = [30e3 0];
fc = 1800e6;

% ------------------------------------------------- full chain, one target
fprintf('-- full chain: generate, synthesise, cancel, correlate\n');
cfg = struct('fc', fc, 'tx_pos', tx, 'rx_pos', rx, ...
             'der_db', 60, 'dnr_db', 60, ...
             'n_clutter', 5, 'clutter_max_range', 3000, ...
             'targets', struct('pos', {[12e3 9e3]}, 'vel', {[80 -40]}));
[s_surv, s_ref, truth] = synth_rx(wave, info, cfg);
t = truth.targets(1);

k_order = ceil(cfg.clutter_max_range / C * info.fs) + 16;
s_clean = eca(s_surv, s_ref, k_order);
[chi, tau, fd] = caf(s_clean, s_ref, info.fs, nb, mlag);

p = abs(chi).^2;
[~, k] = max(p(:));
[r, c] = ind2sub(size(p), k);
range_cell = C / info.fs;
dopp_cell  = fd(2) - fd(1);
ok = check('the injected target appears within one range cell', ...
    abs(tau(r) * C - t.Rb) <= range_cell, ...
    'measured %.1f m, predicted %.1f m, cell %.1f m', tau(r) * C, t.Rb, range_cell) && ok;
ok = check('and within one Doppler cell', abs(fd(c) - t.fd) <= dopp_cell, ...
    'measured %.1f Hz, predicted %.1f Hz, cell %.1f Hz', fd(c), t.fd, dopp_cell) && ok;

% ------------------------------------------------- resolution limit
% Both targets get the same velocity so that they share a Doppler bin and
% the test is purely about range resolution.  No direct path, no clutter and
% no noise, so that nothing but the resolution limit can merge the peaks.
p1  = [12e3 9e3];
vel = [80 -40];
g1  = geometry(tx, rx, p1(:), vel(:), fc, b);
dR  = g1.dR;
fprintf('-- target 1 at (%.0f, %.0f) m, beta = %.1f deg, dR = %.1f m along the bisector\n', ...
    p1(1), p1(2), g1.beta_deg, dR);

for factor = [3.0 0.3]
    p2 = p1(:) + factor * dR * g1.bisector;
    g2 = geometry(tx, rx, p2, vel(:), fc, b);
    drb = abs(g2.Rb - g1.Rb);

    cfg2 = struct('fc', fc, 'tx_pos', tx, 'rx_pos', rx, ...
                  'der_db', 60, 'snr_db', 20, ...
                  'direct_on', false, 'clutter_on', false, ...
                  'targets', struct('pos', {p1, p2(:).'}, 'vel', {vel, vel}));
    [s2, r2, tr2] = synth_rx(wave, info, cfg2);
    [chi2, tau2, fd2] = caf(s2, r2, info.fs, nb, mlag);

    % Range cut through the Doppler bin the targets share.
    [~, cbin] = min(abs(fd2 - tr2.targets(1).fd));
    cut = abs(chi2(:, cbin));
    sel = abs(tau2 * C - g1.Rb) < 6 * dR * cos(g1.beta / 2) * 2 + 200;
    [npk, pk_ranges] = count_peaks(cut(sel), tau2(sel) * C, 6);

    fprintf('   separation %.1f * dR = %.1f m on the ground, %.1f m in bistatic range\n', ...
        factor, factor * dR, drb);
    if factor > 1
        ok = check('two targets separated by more than dR give two peaks', ...
            npk == 2, '%d peaks at %s m, delay resolution c/B = %.1f m', ...
            npk, mat2str(round(pk_ranges)), C / b) && ok;
    else
        ok = check('two targets separated by less than dR give one merged peak', ...
            npk == 1, '%d peak at %s m, delay resolution c/B = %.1f m', ...
            npk, mat2str(round(pk_ranges)), C / b) && ok;
    end
end

% ------------------------------------------------- two targets, two Dopplers
fprintf('-- two targets at different Doppler\n');
cfg3 = struct('fc', fc, 'tx_pos', tx, 'rx_pos', rx, ...
              'der_db', 60, 'dnr_db', 60, 'n_clutter', 5, ...
              'targets', struct('pos', {[12e3 9e3], [20e3 6e3]}, ...
                                'vel', {[80 -40], [-70 25]}));
[s3, r3, tr3] = synth_rx(wave, info, cfg3);
s3c = eca(s3, r3, k_order);
[chi3, tau3, fd3] = caf(s3c, r3, info.fs, nb, mlag);
p3 = abs(chi3).^2;

found = true;
for i = 1:2
    ti = tr3.targets(i);
    [~, ri] = min(abs(tau3 * C - ti.Rb));
    [~, ci] = min(abs(fd3 - ti.fd));
    win = p3(max(1, ri - 2):min(end, ri + 2), max(1, ci - 1):min(end, ci + 1));
    found = found && (max(win(:)) > 20 * median(p3(:)));
    fprintf('   target %d at %.0f m / %.1f Hz reads %.1f dB above the map median\n', ...
        i, ti.Rb, ti.fd, 10 * log10(max(win(:)) / median(p3(:))));
end
ok = check('both targets are visible on the cancelled map', found, ...
    'two peaks, %.0f m apart in range and %.0f Hz apart in Doppler', ...
    abs(tr3.targets(1).Rb - tr3.targets(2).Rb), ...
    abs(tr3.targets(1).fd - tr3.targets(2).fd)) && ok;

fprintf('AT6 result: %s\n', pass_string(ok));
end

% ======================================================================

function [n, positions] = count_peaks(v, axis_vals, db_below)
%COUNT_PEAKS  Local maxima within DB_BELOW dB of the largest value.
v = abs(v(:));
thr = max(v) * 10^(-db_below / 20);
is_pk = [false; v(2:end - 1) > v(1:end - 2) & v(2:end - 1) >= v(3:end); false] & v > thr;
positions = axis_vals(is_pk);
n = numel(positions);
end

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
