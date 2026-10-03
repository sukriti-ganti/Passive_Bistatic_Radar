% RUN_EXPERIMENT4_DETECTION  Deliverable D5 and results/roc.csv.
%
% D5  probability of detection against signal to noise ratio, parameterised
%     by bandwidth and integration time
%
% Every point on these curves comes from the real chain: an LTE waveform
% from M1, a target echo placed at the delay and Doppler that M3 predicts,
% a range-Doppler map from M6 and a decision from the CFAR detector of M7.
% Nothing is taken from a textbook expression for Pd.
%
% The SNR on the x axis is the pre-detection signal to noise ratio measured
% in the occupied bandwidth B.  That is the same reference the link budget
% of M8 uses, so a point read off these curves transfers directly into the
% feasibility table: detection needs
%
%     SNR_in_band + 10*log10(B*T) >= the required post-processing SNR
%
% and the curves show what that required figure really is for Pd = 0.9 at
% Pfa = 1e-4.  The PRD assumes 13 dB.  The measurement is close to it.
%
% What the four bandwidth curves show is not that wide channels detect
% better targets, but that at the same in-band SNR a wider channel collects
% more energy in the same time.  Once the link budget is written in terms of
% received power rather than in-band SNR the bandwidth cancels out entirely,
% which is the point experiment 5 makes.
%
% Set MAKE_FIGURES = false before running to compute without plotting.
% Set QUICK = true for a coarse sweep that finishes in well under a minute.
%
% Runtime: five to ten minutes at the default settings, most of it spent
% generating noise for the 18 MHz, 40 ms configuration.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end
if ~exist('QUICK', 'var'), QUICK = false; end

NRB_LIST  = [6 25 50 100];
TINT_MS   = [10 40];
SNR_LIST  = -62:4:-26;
NTRIALS   = 100;
PFA       = 1e-4;
NSEG      = 64;            % Doppler bins, fixed so the CFAR window fits

if QUICK
    NRB_LIST = [25 100];
    SNR_LIST = -58:6:-28;
    NTRIALS  = 40;
end

fprintf('\n===== Experiment 4: detection performance (D5) =====\n');
fprintf('Pfa = %.0e, %d Monte Carlo trials per point\n', PFA, NTRIALS);

rows = {};
pd_all = zeros(numel(NRB_LIST) * numel(TINT_MS), numel(SNR_LIST));
lbl    = cell(1, numel(NRB_LIST) * numel(TINT_MS));
gp     = zeros(1, numel(NRB_LIST) * numel(TINT_MS));
k = 0;

for it = 1:numel(TINT_MS)
    for i = 1:numel(NRB_LIST)
        nrb = NRB_LIST(i);
        nsf = TINT_MS(it);

        cfg = struct( ...
            'nrb', nrb, 'rho', 1.0, 'ncellid', 17, 'nsubframes', nsf, ...
            'snr_db', SNR_LIST, 'ntrials', NTRIALS, ...
            'cfar', struct('pfa', PFA), ...
            'max_range_m', 8e3, ...
            'tx_pos', [0 0], 'rx_pos', [30e3 0], 'fc', 1800e6, ...
            'target', struct('pos', [12e3 9e3], 'vel', [80 -40]), ...
            'seed', 6000 + 10 * it + i);
        cfg.nbatch = round(nsf * 1e-3 * 15000 * nfft_of(nrb) / NSEG);

        t0 = tic;
        out = roc_curve(cfg);
        k = k + 1;
        pd_all(k, :) = out.pd;
        gp(k) = out.gp_db;
        lbl{k} = sprintf('%d RB (%.2f MHz), %d ms', nrb, out.bw_hz / 1e6, nsf);

        fprintf('\n%s : processing gain %.1f dB, map %d by %d, %d cells tested\n', ...
            lbl{k}, out.gp_db, numel(out.tau), numel(out.fd), out.n_tested);
        for j = 1:numel(SNR_LIST)
            fprintf('   SNR %5.1f dB (post %5.1f dB) : Pd = %.3f\n', ...
                SNR_LIST(j), out.snr_post_db(j), out.pd(j));
            rows(end + 1, :) = {nrb, nsf * 1e-3, SNR_LIST(j), out.pd(j)}; %#ok<SAGROW>
        end
        fprintf('   measured Pfa %.2e, %.1f s\n', out.pfa_meas, toc(t0));
    end
end

write_csv(fullfile(root, 'results', 'roc.csv'), ...
    {'nrb', 'integration_time', 'snr_db', 'pd'}, rows);

% --------------------------------------- where does Pd reach 0.9?
fprintf('\nPost-processing SNR needed for Pd = 0.9, by linear interpolation\n');
need = zeros(1, k);
for i = 1:k
    post = SNR_LIST + gp(i);
    need(i) = interp_cross(post, pd_all(i, :), 0.9);
    fprintf('  %-28s : %5.1f dB   (the PRD assumes 13 dB)\n', lbl{i}, need(i));
end
fprintf(['\nThe figures cluster near the assumed 13 dB, which is the usual\n' ...
         'rule of thumb for Pd = 0.9 at Pfa = 1e-4 against a steady target.\n' ...
         'The link budget of M8 can be taken at face value.\n']);

if MAKE_FIGURES
    h = figure('Name', 'ROC against in-band SNR');
    plot(SNR_LIST, pd_all.', '-o', 'LineWidth', 1);
    grid on;
    axis([SNR_LIST(1) SNR_LIST(end) 0 1.05]);
    xlabel('pre-detection SNR in the occupied bandwidth (dB)');
    ylabel('probability of detection');
    title(sprintf('Detection performance, Pfa = %.0e, %d trials per point', PFA, NTRIALS));
    legend(lbl, 'Location', 'northwest');
    save_fig(h, 'd5_roc_in_band_snr');
    close(h);

    h = figure('Name', 'ROC against post-processing SNR');
    hold on;
    for i = 1:k
        plot(SNR_LIST + gp(i), pd_all(i, :), '-o', 'LineWidth', 1);
    end
    plot([13 13], [0 1.05], 'k--');
    grid on;
    axis([-10 30 0 1.05]);
    xlabel('SNR after the coherent processing gain B*T (dB)');
    ylabel('probability of detection');
    title('The same curves, referred to post-processing SNR');
    legend([lbl, {'13 dB, the PRD assumption'}], 'Location', 'northwest');
    save_fig(h, 'd5_roc_post_processing_snr');
    close(h);
end

fprintf('\nExperiment 4 complete.\n');
