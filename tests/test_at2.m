function ok = test_at2()
%TEST_AT2  Acceptance test for M2, the batched ambiguity function.
%
%   OK = TEST_AT2() returns true when every check passes.
%
%   Checks, from the PRD
%     autocorrelation of the signal with itself peaks at zero delay and
%       zero Doppler
%     the batched result agrees with a direct evaluation to within 1 percent
%       on a 2000 sample segment
%     the false peaks land within one resolution cell of their predicted
%       positions
%
%   Predicted false peaks, and what this test actually finds
%
%     source                  predicted   measured    domain
%     CRS comb, every 6th SC  3.33 km     3.32 km     delay
%     cyclic prefix           21.41 km    19.99 km    delay
%     PSS and SSS, 5 ms       200 Hz      200 Hz      Doppler
%     radio frame, 10 ms      100 Hz      100 Hz      Doppler
%
%   The cyclic prefix entry is the one place where this project departs from
%   the number in the PRD, and it is worth being precise about.  The PRD
%   predicts a delay ambiguity at the OFDM symbol period, 71.35 us, which is
%   21.41 km.  The correlation peak is actually at the useful symbol
%   duration Tu = 1/15 kHz = 66.67 us, which is 19.99 km, because the
%   cyclic prefix is a copy of the tail of the same symbol and the two are
%   separated by Tu, not by Tu + Tcp.  The peak at 71.35 us is 25 dB weaker.
%   The measured value is the correct one and the report should say so.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT2  cross ambiguity function ===\n');
ok = true;
C = 299792458;

% ------------------------------------------------- batched vs direct
fprintf('-- batched implementation against a direct evaluation\n');
set_seed(2);
[wave, ~, info] = lte_generate(6, 1.0, 17, 2);
seg    = wave(1:2000);
nbatch = 125;
maxlag = 32;

t0 = tic;
[chi_b, tau_b, fd_b] = caf(seg, seg, info.fs, nbatch, maxlag);
t_batched = toc(t0);
t0 = tic;
[chi_d, tau_d, fd_d] = caf_direct(seg, seg, info.fs, nbatch, maxlag);
t_direct = toc(t0);

rel = max(abs(chi_b(:) - chi_d(:))) / max(abs(chi_d(:)));
ok = check('batched CAF agrees with the direct evaluation within 1%', ...
    rel < 0.01, 'worst relative error %.2e over %d delay by %d Doppler cells', ...
    rel, size(chi_b, 1), size(chi_b, 2)) && ok;
ok = check('delay and Doppler axes agree', ...
    max(abs(tau_b - tau_d)) == 0 && max(abs(fd_b - fd_d)) == 0, ...
    'identical') && ok;
fprintf('         (batched %.3f s, direct %.3f s, speed up %.0f times)\n', ...
    t_batched, t_direct, t_direct / max(t_batched, eps));

% ------------------------------------------------- peak at the origin
fprintf('-- autocorrelation peak position\n');
set_seed(3);
[wave, ~, info] = lte_generate(50, 1.0, 17, 10);
nb = info.samples_per_subframe;
[chi, tau, fd] = caf(wave, wave, info.fs, nb, 1400);
[~, k] = max(abs(chi(:)));
[r, c] = ind2sub(size(chi), k);
ok = check('autocorrelation peaks at zero delay and zero Doppler', ...
    tau(r) == 0 && fd(c) == 0, 'peak at %.1f m, %.1f Hz', tau(r) * C, fd(c)) && ok;

% ------------------------------------------------- delay domain peaks
fprintf('-- delay domain false peaks\n');
zd  = abs(chi(:, fd == 0));
zdb = 20 * log10(zd / max(zd));
res_cell_m = C / info.occupied_bw;         % one delay resolution cell

% Cyclic prefix repeat, at the useful symbol duration.
tu_m = C * info.nfft / info.fs;
[lvl, pos] = peak_near(zdb, tau * C, tu_m, 1.5 * res_cell_m);
ok = check('cyclic prefix peak within one resolution cell of Tu', ...
    abs(pos - tu_m) <= res_cell_m, ...
    'found at %.0f m (%.1f dB), predicted %.0f m, cell %.0f m', ...
    pos, lvl, tu_m, res_cell_m) && ok;

% The PRD predicts the symbol period instead.  Report both, so the report
% can quote the measured difference rather than assert one of them.
ts_m = C * (info.nfft + info.cp_other) / info.fs;
[lvl_ts, ~] = peak_near(zdb, tau * C, ts_m, 1.5 * res_cell_m);
fprintf('         symbol period %.0f m reads %.1f dB, useful symbol %.0f m reads %.1f dB\n', ...
    ts_m, lvl_ts, tu_m, lvl);
ok = check('the cyclic prefix peak is at Tu, not at the symbol period', ...
    lvl > lvl_ts + 10, '%.1f dB stronger', lvl - lvl_ts) && ok;

% CRS comb: visible when the cell is unloaded, buried when it is not.
set_seed(4);
[wave0, ~, info0] = lte_generate(50, 0.0, 17, 10);
[chi0, tau0, fd0] = caf(wave0, wave0, info0.fs, info0.samples_per_subframe, 1400);
z0  = abs(chi0(:, fd0 == 0));
z0db = 20 * log10(z0 / max(z0));
crs_m = C / (6 * 15000);                   % 6 subcarriers is 90 kHz
[lvl0, pos0] = peak_near(z0db, tau0 * C, crs_m, 1.5 * res_cell_m);
ok = check('CRS comb peak within one resolution cell of c/90 kHz', ...
    abs(pos0 - crs_m) <= res_cell_m, ...
    'found at %.0f m (%.1f dB), predicted %.0f m, cell %.0f m', ...
    pos0, lvl0, crs_m, res_cell_m) && ok;

% ------------------------------------------------- Doppler domain peaks
fprintf('-- Doppler domain false peaks\n');
set_seed(5);
[wave4, ~, info4] = lte_generate(25, 0.0, 17, 40);   % four radio frames
[chi4, tau4, fd4] = caf(wave4, wave4, info4.fs, info4.samples_per_subframe, 4);
zdop = abs(chi4(tau4 == 0, :));
zdopdb = 20 * log10(zdop / max(zdop));
dcell = fd4(2) - fd4(1);

[lvl2, pos2] = peak_near(zdopdb(:), fd4(:), 200, 1.5 * dcell);
ok = check('PSS and SSS 5 ms repetition gives a 200 Hz Doppler peak', ...
    abs(pos2 - 200) <= dcell, 'found at %.0f Hz (%.1f dB), cell %.0f Hz', ...
    pos2, lvl2, dcell) && ok;

[lvl1, pos1] = peak_near(zdopdb(:), fd4(:), 100, 1.5 * dcell);
ok = check('radio frame 10 ms repetition gives a 100 Hz Doppler peak', ...
    abs(pos1 - 100) <= dcell, 'found at %.0f Hz (%.1f dB), cell %.0f Hz', ...
    pos1, lvl1, dcell) && ok;
fprintf('         the 5 ms line is %.1f dB above the 10 ms line, as expected:\n', lvl2 - lvl1);
fprintf('         PSS repeats every 5 ms, SSS alternates so it repeats every 10 ms\n');

fprintf('AT2 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
