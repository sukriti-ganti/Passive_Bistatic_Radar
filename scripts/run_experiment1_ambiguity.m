% RUN_EXPERIMENT1_AMBIGUITY  Deliverable D2: the LTE ambiguity function.
%
% Produces the full ambiguity surface for each of the four LTE bandwidths,
% the zero-Doppler and zero-delay cuts, and a labelled version of the cuts
% in which every false peak is annotated with the feature of the LTE frame
% that causes it.
%
% Predicted false peaks and their causes
%
%   feature                     period        domain    ambiguity
%   cyclic prefix               Tu = 66.67 us delay     19.99 km
%   CRS comb, every 6th SC      1/90 kHz      delay      3.33 km
%   PSS and SSS                 5 ms          Doppler     200 Hz
%   radio frame                 10 ms         Doppler     100 Hz
%
% The cyclic prefix figure needs care.  The cyclic prefix is a copy of the
% tail of its own symbol, and the two copies are separated by the useful
% symbol duration Tu = 1/15 kHz, not by the whole symbol period Tu + Tcp.
% The correlation peak is therefore at 66.67 us, 19.99 km, and the peak at
% the 71.35 us symbol period is about 25 dB weaker.  The PRD predicted the
% latter; the measurement in test_at2 and in this script says the former.
%
% Set MAKE_FIGURES = false before running to compute without plotting.
%
% Runtime: about a minute.  Waveforms are cached in data/.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end

C          = 299792458;
NCELLID    = 17;
NSUBFRAMES = 10;
NRB_LIST   = [6 25 50 100];
MAX_RANGE  = 25e3;            % delay span of the surface, metres
DYN_DB     = 60;

fprintf('\n===== Experiment 1: ambiguity function (D2) =====\n');

% ---------------------------------------------------------------- surfaces
for i = 1:numel(NRB_LIST)
    nrb = NRB_LIST(i);
    set_seed(1000 + nrb);
    [wave, ~, info] = wave_cache(nrb, 1.0, NCELLID, NSUBFRAMES);

    maxlag = ceil(MAX_RANGE / C * info.fs);
    nbatch = info.samples_per_subframe;
    [chi, tau, fd] = caf(wave, wave, info.fs, nbatch, maxlag);

    excl = ceil(2 * info.fs / info.occupied_bw);
    m = af_metrics(chi, tau, fd, excl, 2);
    fprintf('nrb = %3d, B = %5.2f MHz : PSL %6.2f dB, ISL %6.2f dB (main lobe +-%d delay cells)\n', ...
        nrb, info.occupied_bw / 1e6, m.psl_db, m.isl_db, excl);
    fprintf('              worst sidelobe at %8.0f m, %6.0f Hz\n', m.psl_tau * C, m.psl_fd);

    if MAKE_FIGURES
        h = figure('Name', sprintf('ambiguity nrb=%d', nrb));
        plot_rdmap(chi, tau, fd, sprintf(...
            'LTE ambiguity surface, %d RB (%.2f MHz), \\rho = 1, %d ms', ...
            nrb, info.occupied_bw / 1e6, NSUBFRAMES), DYN_DB);
        save_fig(h, sprintf('d2_ambiguity_surface_nrb%d', nrb));
        close(h);
    end
end

% ------------------------------------------------------- zero Doppler cuts
% Two loads, so that the CRS comb can be seen emerging as the traffic falls.
nrb = 50;
cuts = struct('rho', {0.0, 1.0}, 'tau', {[], []}, 'db', {[], []});
for i = 1:numel(cuts)
    set_seed(2000 + i);
    [wave, ~, info] = wave_cache(nrb, cuts(i).rho, NCELLID, NSUBFRAMES);
    maxlag = ceil(MAX_RANGE / C * info.fs);
    [chi, tau, fd] = caf(wave, wave, info.fs, info.samples_per_subframe, maxlag);
    z = abs(chi(:, fd == 0));
    cuts(i).tau = tau;
    cuts(i).db  = 20 * log10(z / max(z));
    cuts(i).info = info;
end

tu_km  = C * cuts(1).info.nfft / cuts(1).info.fs / 1e3;
ts_km  = C * (cuts(1).info.nfft + cuts(1).info.cp_other) / cuts(1).info.fs / 1e3;
crs_km = C / (6 * 15000) / 1e3;

fprintf('\nzero-Doppler cut, %d RB, level at the predicted false peaks\n', nrb);
fprintf('  %-32s %10s %10s %12s\n', 'feature', 'rho = 0', 'rho = 1', 'position');
delay_features = { ...
    'cyclic prefix, Tu = 66.67 us', tu_km; ...
    'OFDM symbol period, 71.35 us', ts_km; ...
    'CRS comb, 90 kHz spacing',     crs_km; ...
    'CRS comb, second harmonic',    2 * crs_km};
for i = 1:size(delay_features, 1)
    lvl = zeros(1, 2); pos = zeros(1, 2);
    for j = 1:2
        [lvl(j), pos(j)] = peak_near(cuts(j).db, cuts(j).tau * C / 1e3, ...
                                     delay_features{i, 2}, 0.06);
    end
    fprintf('  %-32s %7.1f dB %7.1f dB %9.2f km\n', ...
        delay_features{i, 1}, lvl(1), lvl(2), pos(2));
end

if MAKE_FIGURES
    h = figure('Name', 'zero Doppler cut');
    plot(cuts(1).tau * C / 1e3, cuts(1).db, 'LineWidth', 1); hold on;
    plot(cuts(2).tau * C / 1e3, cuts(2).db, 'LineWidth', 1);
    yl = [-70 2];
    for x = [tu_km, crs_km, 2 * crs_km, -tu_km, -crs_km, -2 * crs_km]
        plot([x x], yl, 'k:');
    end
    text(tu_km,  -4, ' cyclic prefix, T_u');
    text(crs_km, -12, ' CRS comb');
    axis([-25 25 yl]);
    grid on;
    xlabel('bistatic range (km)');
    ylabel('normalised |\chi| (dB)');
    title(sprintf('Zero Doppler cut, %d RB (%.2f MHz), %d ms integration', ...
        nrb, cuts(1).info.occupied_bw / 1e6, NSUBFRAMES));
    legend('\rho = 0, reference and sync only', '\rho = 1, fully loaded', ...
        'Location', 'southeast');
    save_fig(h, 'd2_zero_doppler_cut');
    close(h);
end

% --------------------------------------------------------- zero delay cuts
% Four radio frames, so that the Doppler axis resolves 25 Hz and the 5 ms
% and 10 ms repetition lines are separated.
nrb_d = 25;
dcuts = struct('rho', {0.0, 1.0}, 'fd', {[], []}, 'db', {[], []});
for i = 1:numel(dcuts)
    set_seed(3000 + i);
    [wave, ~, info] = wave_cache(nrb_d, dcuts(i).rho, NCELLID, 40);
    [chi, tau, fd] = caf(wave, wave, info.fs, info.samples_per_subframe, 4);
    z = abs(chi(tau == 0, :));
    dcuts(i).fd = fd;
    dcuts(i).db = 20 * log10(z / max(z));
end

fprintf('\nzero-delay cut, %d RB, four radio frames, Doppler cell %.0f Hz\n', ...
    nrb_d, dcuts(1).fd(2) - dcuts(1).fd(1));
fprintf('  %-32s %10s %10s %12s\n', 'feature', 'rho = 0', 'rho = 1', 'position');
dopp_features = { ...
    'PSS and SSS, 5 ms period',     200; ...
    'radio frame, 10 ms period',    100; ...
    'PSS and SSS, second harmonic', 400};
for i = 1:size(dopp_features, 1)
    lvl = zeros(1, 2); pos = zeros(1, 2);
    for j = 1:2
        [lvl(j), pos(j)] = peak_near(dcuts(j).db, dcuts(j).fd, dopp_features{i, 2}, 13);
    end
    fprintf('  %-32s %7.1f dB %7.1f dB %9.0f Hz\n', ...
        dopp_features{i, 1}, lvl(1), lvl(2), pos(2));
end

if MAKE_FIGURES
    h = figure('Name', 'zero delay cut');
    plot(dcuts(1).fd, dcuts(1).db, 'LineWidth', 1); hold on;
    plot(dcuts(2).fd, dcuts(2).db, 'LineWidth', 1);
    yl = [-80 2];
    for x = [-400 -200 200 400]
        plot([x x], yl, 'k:');
    end
    for x = [-300 -100 100 300]
        plot([x x], yl, 'k--');
    end
    axis([-500 500 yl]);
    grid on;
    xlabel('bistatic Doppler (Hz)');
    ylabel('normalised |\chi| (dB)');
    title(sprintf('Zero delay cut, %d RB, 40 ms integration', nrb_d));
    legend('\rho = 0, reference and sync only', '\rho = 1, fully loaded', ...
        'Location', 'southeast');
    save_fig(h, 'd2_zero_delay_cut');
    close(h);
end

fprintf('\nExperiment 1 complete.\n');
