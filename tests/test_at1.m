function ok = test_at1()
%TEST_AT1  Acceptance test for M1, the LTE downlink waveform generator.
%
%   OK = TEST_AT1() returns true when every check passes.
%
%   Checks, from the PRD
%     peak to average power ratio between 8 and 12 dB at rho = 1
%     occupied bandwidth from the spectrum, within 5 percent of the table
%     DC subcarrier power zero
%     number of samples exactly nsubframes * 0.001 * fs
%     CRS positions form a regular lattice at the expected spacing
%
%   AT1 has to pass before any other module is trusted: an indexing error in
%   the generator would quietly invalidate every result downstream.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT1  LTE downlink waveform generator ===\n');
ok = true;

table_nrb  = [6 25 50 100];
table_nfft = [128 512 1024 2048];
table_nsc  = [72 300 600 1200];
table_bw   = [1.08 4.50 9.00 18.00] * 1e6;

nsubframes = 2;
ncellid    = 17;

for i = 1:numel(table_nrb)
    nrb = table_nrb(i);
    set_seed(100 + nrb);
    [wave, grid, info] = lte_generate(nrb, 1.0, ncellid, nsubframes);

    fprintf('-- nrb = %d (%d point FFT, %.2f MHz sample rate)\n', ...
        nrb, info.nfft, info.fs / 1e6);

    % --- parameter table ------------------------------------------------
    ok = check('nfft matches the table', info.nfft == table_nfft(i), ...
        '%d', info.nfft) && ok;
    ok = check('fs = nfft * 15000', info.fs == info.nfft * 15000, ...
        '%.2f MHz', info.fs / 1e6) && ok;
    ok = check('occupied subcarriers', info.nsc == table_nsc(i), ...
        '%d', info.nsc) && ok;
    ok = check('nominal occupied bandwidth', ...
        abs(info.occupied_bw - table_bw(i)) < 1, '%.2f MHz', info.occupied_bw / 1e6) && ok;

    % --- sample count ---------------------------------------------------
    expected_n = nsubframes * 0.001 * info.fs;
    ok = check('sample count is exact', numel(wave) == expected_n, ...
        '%d samples, expected %d', numel(wave), expected_n) && ok;

    % --- grid size ------------------------------------------------------
    ok = check('grid is 12*nrb by 14*nsubframes', ...
        all(size(grid) == [12 * nrb, 14 * nsubframes]), ...
        '%d by %d', size(grid, 1), size(grid, 2)) && ok;

    % --- PAPR -----------------------------------------------------------
    papr = info.papr_db;
    ok = check('PAPR between 8 and 12 dB at rho = 1', ...
        papr >= 8 && papr <= 12, '%.2f dB', papr) && ok;

    % --- DC subcarrier --------------------------------------------------
    % Strip the cyclic prefix from every symbol and look at the DC bin of
    % the transform.  This is the bin the modulator is responsible for; the
    % spectrum of the whole waveform has CP leakage from the neighbours in
    % it and would not be a fair test of the mapping.
    dc = 0; pk = 0; p = 1;
    for s = 1:(14 * nsubframes)
        cp = info.cplen(mod(s - 1, 14) + 1);
        x  = wave(p + cp + (0:info.nfft - 1));
        xf = fft(x);
        dc = max(dc, abs(xf(1)));
        pk = max(pk, max(abs(xf)));
        p  = p + cp + info.nfft;
    end
    ok = check('DC subcarrier power is zero', dc < 1e-9 * pk, ...
        '%.1f dB below the strongest subcarrier', 20 * log10(pk / max(dc, realmin))) && ok;

    % --- occupied bandwidth measured from the spectrum ------------------
    meas_bw = measured_occupied_bw(wave, info.fs, 0.99);
    err_pct = 100 * abs(meas_bw - info.occupied_bw) / info.occupied_bw;
    ok = check('measured 99% occupied bandwidth within 5%', err_pct < 5, ...
        '%.3f MHz measured vs %.3f MHz nominal, %.2f%% error', ...
        meas_bw / 1e6, info.occupied_bw / 1e6, err_pct) && ok;

    % --- CRS lattice ----------------------------------------------------
    % Regenerate with an empty PDSCH so that only the physical signals are
    % present, then confirm the comb.
    [~, g0] = lte_generate(nrb, 0.0, ncellid, nsubframes);
    crs_cols = [1 5 8 12];                      % symbols 0, 4, 7, 11, 1-based
    lattice_ok = true;
    for cc = crs_cols
        rows = find(abs(g0(:, cc)) > 0);
        lattice_ok = lattice_ok && numel(rows) == 2 * nrb ...
                  && rows(1) == mod(ncellid, 6) + 1 ...
                  && all(diff(rows) == 6);
    end
    rows = find(abs(g0(:, 1)) > 0);
    ok = check('CRS lattice: 2*nrb per symbol, spacing 6, shift mod(ncellid,6)', ...
        lattice_ok, 'first row %d, spacing %d, %d per symbol, symbols %s', ...
        rows(1), rows(2) - rows(1), numel(rows), mat2str(crs_cols - 1)) && ok;

    % --- CRS amplitude is unit power QPSK -------------------------------
    ok = check('CRS symbols are unit power QPSK', ...
        max(abs(abs(g0(rows, 1)) - 1)) < 1e-12, ...
        'max magnitude error %.1e', max(abs(abs(g0(rows, 1)) - 1))) && ok;

    % --- PSS and SSS placement ------------------------------------------
    pss_rows = find(abs(g0(:, 7)) > 0);
    sss_rows = find(abs(g0(:, 6)) > 0);
    centre   = (info.nsc / 2 - 30):(info.nsc / 2 + 31);
    ok = check('PSS on the central 62 subcarriers of the last symbol of slot 0', ...
        numel(pss_rows) == 62 && all(pss_rows(:).' == centre), ...
        'rows %d to %d', pss_rows(1), pss_rows(end)) && ok;
    ok = check('SSS on the central 62 subcarriers of the previous symbol', ...
        numel(sss_rows) == 62 && all(sss_rows(:).' == centre), ...
        'rows %d to %d', sss_rows(1), sss_rows(end)) && ok;
end

% --- load factor definition --------------------------------------------
fprintf('-- load factor\n');
set_seed(7);
[~, ~, i0] = lte_generate(50, 0.0, 17, 2);
[~, ~, i5] = lte_generate(50, 0.5, 17, 2);
[~, ~, i1] = lte_generate(50, 1.0, 17, 2);
ok = check('rho = 0 leaves the PDSCH empty', i0.n_pdsch_used == 0, ...
    '%d of %d PDSCH resource elements filled', i0.n_pdsch_used, i0.n_pdsch_avail) && ok;
ok = check('rho = 0.5 fills half the PDSCH', abs(i5.rho_actual - 0.5) < 1e-3, ...
    '%.4f', i5.rho_actual) && ok;
ok = check('rho = 1 fills the PDSCH', abs(i1.rho_actual - 1) < 1e-12, ...
    '%.4f', i1.rho_actual) && ok;

fprintf('AT1 result: %s\n', pass_string(ok));
end

% ======================================================================

function bw = measured_occupied_bw(wave, fs, frac)
%MEASURED_OCCUPIED_BW  Two sided bandwidth holding FRAC of the total power.
n  = numel(wave);
p  = abs(fft(wave)).^2;
f  = (0:n - 1) / n * fs;
f(f >= fs / 2) = f(f >= fs / 2) - fs;
[fa, k] = sort(abs(f));
cp = cumsum(p(k));
cp = cp / cp(end);
i  = find(cp >= frac, 1);
bw = 2 * fa(i);
end

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
