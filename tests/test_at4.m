function ok = test_at4()
%TEST_AT4  Acceptance test for M4, the received signal synthesis.
%
%   OK = TEST_AT4() returns true when every check passes.
%
%   Checks, from the PRD
%     direct path only, no targets, no noise: the range-Doppler map shows a
%       single peak at zero delay and zero Doppler
%     one target, no direct path, no noise: the peak lands at the delay and
%       Doppler that M3 predicts
%     the measured direct to echo ratio matches the configured value to
%       within 1 dB
%
%   The direct to echo ratio is measured with a matched replica rather than
%   from the map, because a peak that falls between grid points reads one to
%   two dB low from straddle loss alone.  MATCH_AMP has no such quantisation.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT4  received signal synthesis ===\n');
ok = true;
C = 299792458;

set_seed(21);
[wave, ~, info] = lte_generate(50, 1.0, 17, 20);
nb   = info.samples_per_subframe;
mlag = 400;

base = struct('fc', 1800e6, 'tx_pos', [0 0], 'rx_pos', [30e3 0], ...
              'der_db', 60, 'clutter_on', false, 'noise_on', false);

% ------------------------------------------------- direct path only
cfg = base;
cfg.targets   = [];
cfg.direct_on = true;
[s_surv, s_ref] = synth_rx(wave, info, cfg);

ok = check('direct path only reproduces the reference exactly', ...
    max(abs(s_surv - s_ref)) == 0, 'max difference %.1e', max(abs(s_surv - s_ref))) && ok;

[chi, tau, fd] = caf(s_surv, s_ref, info.fs, nb, mlag);
p = abs(chi).^2;
[~, k] = max(p(:));
[r, c] = ind2sub(size(p), k);
ok = check('its range-Doppler map peaks at zero delay and zero Doppler', ...
    tau(r) == 0 && fd(c) == 0, 'peak at %.0f m, %.0f Hz', tau(r) * C, fd(c)) && ok;

% A single peak, in the sense that matters: no second peak that could be
% mistaken for a target.  The floor around it is the LTE signal's own
% ambiguity sidelobe structure, which D2 characterises; the point here is
% that nothing else comes close to the direct path itself.
m = af_metrics(chi, tau, fd, 8, 2);
ok = check('and nothing else on the map comes within 15 dB of it', ...
    m.psl_db < -15, 'worst sidelobe %.1f dB at %.0f m, %.0f Hz', ...
    m.psl_db, m.psl_tau * C, m.psl_fd) && ok;

% ------------------------------------------------- one target, no direct
cfg = base;
cfg.direct_on = false;
cfg.targets = struct('pos', [12e3 9e3], 'vel', [80 -40]);
[s_surv, s_ref, truth] = synth_rx(wave, info, cfg);
t = truth.targets(1);

g = geometry(cfg.tx_pos, cfg.rx_pos, t.pos(:), t.vel(:), cfg.fc, info.occupied_bw);
fprintf('-- target at (%.0f, %.0f) m: M3 predicts Rb = %.1f m, fd = %.1f Hz\n', ...
    t.pos(1), t.pos(2), g.Rb, g.fd);

[chi, tau, fd] = caf(s_surv, s_ref, info.fs, nb, mlag);
p = abs(chi).^2;
[~, k] = max(p(:));
[r, c] = ind2sub(size(p), k);

range_cell_m = C / info.fs;                 % delay axis sample spacing
dopp_cell_hz = fd(2) - fd(1);
ok = check('the peak is within one range cell of the M3 prediction', ...
    abs(tau(r) * C - g.Rb) <= range_cell_m, ...
    'measured %.1f m, predicted %.1f m, cell %.1f m', tau(r) * C, g.Rb, range_cell_m) && ok;
ok = check('the peak is within one Doppler cell of the M3 prediction', ...
    abs(fd(c) - g.fd) <= dopp_cell_hz, ...
    'measured %.1f Hz, predicted %.1f Hz, cell %.1f Hz', fd(c), g.fd, dopp_cell_hz) && ok;

% ------------------------------------------------- direct to echo ratio
for der = [40 60 80]
    cfg = base;
    cfg.der_db = der;
    cfg.direct_on = false;
    cfg.targets = struct('pos', [12e3 9e3], 'vel', [80 -40]);
    [s_t, s_r, tr] = synth_rx(wave, info, cfg);
    a_echo = match_amp(s_t, s_r, info.fs, tr.targets(1).delay_samples, tr.targets(1).fd);

    cfg.direct_on = true;
    cfg.targets = [];
    [s_d, s_r2] = synth_rx(wave, info, cfg);
    a_dir = match_amp(s_d, s_r2, info.fs, 0, 0);

    meas = 20 * log10(abs(a_dir) / abs(a_echo));
    ok = check('measured direct to echo ratio within 1 dB', abs(meas - der) < 1, ...
        'configured %d dB, measured %.3f dB', der, meas) && ok;
end

% ------------------------------------------------- fractional delay
fprintf('-- fractional delay interpolation\n');
d_int = 37;
y = frac_delay(wave, d_int);
ok = check('an integer delay is exact', ...
    max(abs(y(200:end - 200) - wave(200 - d_int:end - 200 - d_int))) < 1e-9, ...
    'max error %.1e', max(abs(y(200:end - 200) - wave(200 - d_int:end - 200 - d_int)))) && ok;

% Against an analytic signal, where the delayed waveform is known in closed
% form: a tone at f must come back phase shifted by exactly -2*pi*f*d/fs.
n  = 8192;
f0 = 0.17 * info.fs;
nn = (0:n - 1).';
tone = exp(2j * pi * f0 * nn / info.fs);
d    = 3.27;
yt   = frac_delay(tone, d);
want = exp(2j * pi * f0 * (nn - d) / info.fs);
mid  = round(n / 2);
err  = max(abs(yt(mid - 100:mid + 100) - want(mid - 100:mid + 100)));
ok = check('a fractional delay of a tone matches the closed form', err < 1e-5, ...
    'delay %.2f samples, worst error %.1e in the middle of the record', d, err) && ok;

% The ideal interpolator has infinite support, so the ends of a finite
% record are reconstructed from data that is not there.  That error must
% decay as 1/distance from the edge, and it must stay small.
e = abs(yt - want);
ok = check('the edge error decays as 1/distance and is negligible by 200 samples in', ...
    e(200) < 2e-3 && e(200) > 10 * e(2000), ...
    '%.1e at 200 samples, %.1e at 2000, %.1e in the middle', e(200), e(2000), err) && ok;

% Delay resolution finer than a sample: sub sample shifts must move the peak.
cfgA = base; cfgA.direct_on = false;
cfgA.targets = struct('pos', [12e3 9e3], 'vel', [0 0]);
[sa, ra, ta] = synth_rx(wave, info, cfgA);
aa = match_amp(sa, ra, info.fs, ta.targets(1).delay_samples, 0);
ab = match_amp(sa, ra, info.fs, round(ta.targets(1).delay_samples), 0);
ok = check('a rounded delay loses matched filter gain, a fractional one does not', ...
    abs(aa) > abs(ab), ...
    'fractional %.4e, rounded %.4e, loss %.2f dB', ...
    abs(aa), abs(ab), 20 * log10(abs(aa) / abs(ab))) && ok;

% ------------------------------------------------- clutter and noise
fprintf('-- clutter and noise\n');
cfg = base;
cfg.clutter_on = true;
cfg.noise_on   = true;
cfg.n_clutter  = 5;
cfg.dnr_db     = 60;
cfg.targets    = struct('pos', {[12e3 9e3], [20e3 5e3]}, 'vel', {[80 -40], [-60 20]});
[s_surv, s_ref, truth] = synth_rx(wave, info, cfg);
ok = check('the requested number of clutter scatterers is present', ...
    numel(truth.clutter) == 5, '%d scatterers, delays %s samples', ...
    numel(truth.clutter), mat2str(round([truth.clutter.delay_samples]))) && ok;
ok = check('the direct to noise ratio is as configured', ...
    abs(truth.dnr_db - 60) < 0.5, 'measured %.2f dB', truth.dnr_db) && ok;
ok = check('two moving targets are present with distinct Doppler', ...
    numel(truth.targets) == 2 && abs(truth.targets(1).fd - truth.targets(2).fd) > 1, ...
    'fd = %.1f Hz and %.1f Hz', truth.targets(1).fd, truth.targets(2).fd) && ok;

fprintf('AT4 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
