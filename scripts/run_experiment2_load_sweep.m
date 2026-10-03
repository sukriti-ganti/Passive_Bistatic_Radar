% RUN_EXPERIMENT2_LOAD_SWEEP  Deliverables D1 and D3.
%
% D1  the waveform library: four bandwidths by five load factors, generated
%     and written to data/wave_nrb{N}_rho{RRR}.mat
% D3  the load curve: peak sidelobe level and integrated sidelobe level
%     against the resource occupancy rho, plus results/af_metrics.csv
%
% The question this answers is the second of the three the project asks:
% how does detection performance depend on the traffic the cell is carrying?
% An LTE downlink is only noise-like where it carries data.  With the cell
% idle, all that is left is the reference and synchronisation signals, and
% those are strictly periodic: the CRS comb repeats every sixth subcarrier
% and puts a false peak at c/90 kHz = 3.33 km that is within a few dB of the
% true peak.  A passive radar looking at an idle cell would report a ghost
% target at 3.33 km, and another at 6.66 km, for every real one.
%
% Set MAKE_FIGURES = false before running to compute without plotting.
%
% Runtime: a couple of minutes on the first run, seconds afterwards because
% wave_cache reloads the library from data/.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end

C          = 299792458;
NCELLID    = 17;
NSUBFRAMES = 10;
NRB_LIST   = [6 25 50 100];
RHO_LIST   = [0 0.25 0.5 0.75 1.0];
MAX_RANGE  = 25e3;

fprintf('\n===== Experiment 2: load sweep (D1, D3) =====\n');

n_nrb = numel(NRB_LIST);
n_rho = numel(RHO_LIST);
psl   = zeros(n_nrb, n_rho);
isl   = zeros(n_nrb, n_rho);
crs   = zeros(n_nrb, n_rho);      % CRS comb false peak level, dB
rows  = cell(n_nrb * n_rho, 4);
k     = 0;

crs_km = C / (6 * 15000) / 1e3;

for i = 1:n_nrb
    nrb = NRB_LIST(i);
    for j = 1:n_rho
        rho = RHO_LIST(j);
        set_seed(4000 + 10 * i + j);
        [wave, ~, info] = wave_cache(nrb, rho, NCELLID, NSUBFRAMES);

        maxlag = ceil(MAX_RANGE / C * info.fs);
        [chi, tau, fd] = caf(wave, wave, info.fs, info.samples_per_subframe, maxlag);

        % Main lobe is two delay resolution cells either side of the peak,
        % and two Doppler bins.  Tying the exclusion to fs/B rather than to
        % a fixed number of samples is what makes the four bandwidths
        % comparable: they are oversampled by different factors.
        excl = ceil(2 * info.fs / info.occupied_bw);
        m = af_metrics(chi, tau, fd, excl, 2);

        psl(i, j) = m.psl_db;
        isl(i, j) = m.isl_db;

        z = abs(chi(:, fd == 0));
        zdb = 20 * log10(z / max(z));
        crs(i, j) = peak_near(zdb, tau * C / 1e3, crs_km, 0.06);

        k = k + 1;
        rows(k, :) = {nrb, rho, m.psl_db, m.isl_db};

        fprintf('nrb %3d, rho %.2f : PSL %6.2f dB, ISL %6.2f dB, CRS ghost at 3.33 km %6.1f dB\n', ...
            nrb, rho, m.psl_db, m.isl_db, crs(i, j));
    end
end

write_csv(fullfile(root, 'results', 'af_metrics.csv'), ...
    {'nrb', 'rho', 'psl_db', 'isl_db'}, rows);

fprintf('\nSummary: going from an idle cell to a fully loaded one\n');
for i = 1:n_nrb
    fprintf('  nrb %3d : PSL %6.2f -> %6.2f dB, ISL %6.2f -> %6.2f dB, CRS ghost %6.1f -> %6.1f dB\n', ...
        NRB_LIST(i), psl(i, 1), psl(i, end), isl(i, 1), isl(i, end), crs(i, 1), crs(i, end));
end

if MAKE_FIGURES
    h = figure('Name', 'load curve');

    subplot(1, 3, 1);
    plot(RHO_LIST, psl.', '-o', 'LineWidth', 1);
    grid on;
    xlabel('resource occupancy \rho');
    ylabel('peak sidelobe level (dB)');
    title('PSL against load');
    legend(arrayfun(@(n) sprintf('%d RB', n), NRB_LIST, 'UniformOutput', false), ...
        'Location', 'best');

    subplot(1, 3, 2);
    plot(RHO_LIST, isl.', '-o', 'LineWidth', 1);
    grid on;
    xlabel('resource occupancy \rho');
    ylabel('integrated sidelobe level (dB)');
    title('ISL against load');

    subplot(1, 3, 3);
    plot(RHO_LIST, crs.', '-o', 'LineWidth', 1);
    grid on;
    xlabel('resource occupancy \rho');
    ylabel('CRS ghost level at 3.33 km (dB)');
    title('The false peak the CRS comb makes');

    save_fig(h, 'd3_load_curve');
    close(h);
end

fprintf('\nExperiment 2 complete.\n');
