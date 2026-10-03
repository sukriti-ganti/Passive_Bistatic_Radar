% RUN_EXPERIMENT3_CANCELLATION  Deliverable D4 and results/cancellation.csv.
%
% D4  range-Doppler map before and after direct signal cancellation
%     cancellation.csv: filter order against cancellation depth and target loss
%
% The direct path arrives 60 dB above the target echoes here, and the LTE
% ambiguity function has a peak sidelobe level of about -21 dB at full load.
% Untreated, the direct path alone therefore paints the whole map at roughly
% 39 dB above the target, and no amount of integration helps because the
% sidelobes integrate coherently along with the peak.  The two maps this
% script draws are the single most convincing figure in the report.
%
% The script also compares ECA with NLMS, which is the report section the
% PRD asks for.  The comparison is not close.  NLMS has a single knob, the
% step size, and it sets both how fast the filter converges and how large a
% misadjustment floor it settles at.  Over a 20 ms record no step size gets
% NLMS within 25 dB of ECA: small steps have not converged by the end of the
% record, large steps converge quickly to a noisy solution.  ECA has no such
% coupling because it solves for the weights in closed form.
%
% Read the NLMS target loss column together with its depth column.  The loss
% is measured by matching a replica against the canceller output, and while
% a strong direct path is still present that estimate is biased, so a loss
% figure taken from a run with only 6 dB of cancellation means nothing.
%
% The last section measures the Doppler notch that cancellation digs around
% zero Doppler.  This is the mechanism behind the target loss requirement in
% AT5: because the weights are fitted to data that contains the target, the
% fit absorbs part of any echo that looks stationary over a batch.  Note the
% precondition - ECA spans delays 0 to K samples and nothing outside it, so
% a target beyond K*c/fs of bistatic range is untouched whatever its
% Doppler.  The notch probe is deliberately placed at 1 km, well inside the
% span, or the measurement would report a comfortable zero that means
% nothing.
%
% Set MAKE_FIGURES = false before running to compute without plotting.
%
% Runtime: two to three minutes, most of it the NLMS sample loop.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end

C          = 299792458;
NRB        = 50;
RHO        = 1.0;
NCELLID    = 17;
NSUBFRAMES = 20;
MAX_RANGE  = 10e3;
K_LIST     = [0 16 32 48 64 80 96 112 128 144 160 192 224];
MU_LIST    = [0.001 0.003 0.01 0.03 0.1];
CLUTTER_RANGE = 3000;

fprintf('\n===== Experiment 3: direct signal cancellation (D4) =====\n');

set_seed(5001);
[wave, ~, info] = wave_cache(NRB, RHO, NCELLID, NSUBFRAMES);

cfg = struct( ...
    'fc', 1800e6, 'tx_pos', [0 0], 'rx_pos', [30e3 0], ...
    'der_db', 60, 'dnr_db', 60, ...
    'n_clutter', 5, 'clutter_max_range', CLUTTER_RANGE, 'clutter_db', -25, ...
    'targets', struct('pos', {[12e3 9e3], [20e3 6e3]}, ...
                      'vel', {[80 -40],   [-70 25]}), ...
    'seed', 5002);
[s_surv, s_ref, truth] = synth_rx(wave, info, cfg);

maxlag = ceil(MAX_RANGE / C * info.fs);
nbatch = info.samples_per_subframe;
k_need = ceil(CLUTTER_RANGE / C * info.fs);

fprintf('scenario: %d RB, %.0f ms, DER %.0f dB, DNR %.1f dB, %d clutter scatterers\n', ...
    NRB, NSUBFRAMES, cfg.der_db, truth.dnr_db, cfg.n_clutter);
for i = 1:numel(truth.targets)
    fprintf('  target %d : %8.1f m bistatic range, %7.2f Hz Doppler, SNR %.1f dB per sample\n', ...
        i, truth.targets(i).Rb, truth.targets(i).fd, truth.targets(i).snr_db);
end
fprintf('  clutter reaches %.0f samples, so K must reach at least %d\n', ...
    max([truth.clutter.delay_samples]), k_need);

% -------------------------------------------------- cancellation depth vs K
fprintf('\nECA, filter order sweep\n');
depth = zeros(size(K_LIST));
loss  = zeros(numel(K_LIST), numel(truth.targets));
rows  = cell(numel(K_LIST), 3);

for i = 1:numel(K_LIST)
    [s_clean, ~, ci] = eca(s_surv, s_ref, K_LIST(i));
    depth(i) = ci.depth_db;
    for t = 1:numel(truth.targets)
        a = match_amp(s_clean, s_ref, info.fs, ...
                      truth.targets(t).delay_samples, truth.targets(t).fd);
        loss(i, t) = 20 * log10(truth.targets(t).amp / abs(a));
    end
    rows(i, :) = {K_LIST(i), depth(i), max(loss(i, :))};
    fprintf('  K = %3d : depth %6.2f dB, target loss %5.2f / %5.2f dB\n', ...
        K_LIST(i), depth(i), loss(i, 1), loss(i, 2));
end

write_csv(fullfile(root, 'results', 'cancellation.csv'), ...
    {'filter_order', 'depth_db', 'target_loss_db'}, rows);

fprintf(['\nThe depth climbs while K is still short of the furthest clutter\n' ...
         'scatterer, then flattens at %.1f dB.  The ceiling is the direct to\n' ...
         'noise ratio of %.1f dB: the residual still contains the receiver\n' ...
         'noise and no projection can remove that.\n'], max(depth), truth.dnr_db);
fprintf(['A negative target loss at very small K is not the canceller\n' ...
         'amplifying anything.  It is the matched amplitude estimate being\n' ...
         'biased by the direct path that is still there, and it disappears\n' ...
         'as soon as the cancellation works.\n']);

% -------------------------------------------------- before and after maps
k_best = K_LIST(find(K_LIST >= k_need + 16, 1));
if isempty(k_best), k_best = K_LIST(end); end
[s_clean, ~, ci] = eca(s_surv, s_ref, k_best);

[chi_b, tau_b, fd_b] = caf(s_surv,  s_ref, info.fs, nbatch, maxlag);
[chi_a, tau_a, fd_a] = caf(s_clean, s_ref, info.fs, nbatch, maxlag);
pk = max(abs(chi_b(:)));

fprintf('\nrange-Doppler maps at K = %d (depth %.1f dB)\n', k_best, ci.depth_db);
for t = 1:numel(truth.targets)
    [~, r] = min(abs(tau_b - truth.targets(t).tau));
    [~, c] = min(abs(fd_b - truth.targets(t).fd));
    fprintf(['  target %d cell stands %5.1f dB above the map median before, ' ...
             '%5.1f dB after\n'], t, ...
        10 * log10(abs(chi_b(r, c))^2 / median(abs(chi_b(:)).^2)), ...
        10 * log10(abs(chi_a(r, c))^2 / median(abs(chi_a(:)).^2)));
end

if MAKE_FIGURES
    h = figure('Name', 'cancellation before and after');
    subplot(1, 2, 1);
    plot_rdmap(chi_b, tau_b, fd_b, 'Before cancellation', 80, pk);
    subplot(1, 2, 2);
    plot_rdmap(chi_a, tau_a, fd_a, sprintf('After ECA, K = %d', k_best), 80, pk);
    save_fig(h, 'd4_rdmap_before_after');
    close(h);

    h = figure('Name', 'cancellation depth against filter order');
    subplot(1, 2, 1);
    plot(K_LIST, depth, '-o', 'LineWidth', 1); hold on;
    plot([k_need k_need], [0 max(depth) + 5], 'k--');
    plot([K_LIST(1) K_LIST(end)], [truth.dnr_db truth.dnr_db], 'r:');
    grid on;
    xlabel('filter order K (samples)');
    ylabel('cancellation depth (dB)');
    title('Depth against filter order');
    legend('measured', 'furthest clutter', 'direct to noise ratio', 'Location', 'southeast');

    subplot(1, 2, 2);
    plot(K_LIST, loss, '-o', 'LineWidth', 1); hold on;
    plot([K_LIST(1) K_LIST(end)], [3 3], 'r:');
    grid on;
    xlabel('filter order K (samples)');
    ylabel('target loss (dB)');
    title('Target loss against filter order');
    legend('target 1', 'target 2', '3 dB limit', 'Location', 'northeast');
    save_fig(h, 'd4_cancellation_depth');
    close(h);
end

% -------------------------------------------------- ECA against NLMS
fprintf(['\nNLMS at the same filter order, K = %d.  Read the loss column\n' ...
         'together with the depth column: while the direct path is still\n' ...
         'largely present the matched amplitude estimate is biased.\n'], k_best);
n_depth = zeros(size(MU_LIST));
n_loss  = zeros(size(MU_LIST));
for i = 1:numel(MU_LIST)
    [s_n, ~, ni] = nlms(s_surv, s_ref, k_best, MU_LIST(i));
    a = match_amp(s_n, s_ref, info.fs, ...
                  truth.targets(1).delay_samples, truth.targets(1).fd);
    n_depth(i) = ni.depth_db;
    n_loss(i)  = 20 * log10(truth.targets(1).amp / abs(a));
    fprintf('  mu = %6.3f : depth %6.2f dB, target 1 loss %6.2f dB, converged after %5.2f ms\n', ...
        MU_LIST(i), n_depth(i), n_loss(i), ni.converge_samples / info.fs * 1e3);
end
fprintf('  ECA        : depth %6.2f dB, target 1 loss %6.2f dB\n', ...
    ci.depth_db, loss(K_LIST == k_best, 1));

if MAKE_FIGURES
    h = figure('Name', 'ECA against NLMS');
    subplot(1, 2, 1);
    semilogx(MU_LIST, n_depth, '-o', 'LineWidth', 1); hold on;
    semilogx([MU_LIST(1) MU_LIST(end)], [ci.depth_db ci.depth_db], 'k--');
    grid on;
    xlabel('NLMS step size \mu');
    ylabel('cancellation depth (dB)');
    title(sprintf('Cancellation depth, K = %d', k_best));
    legend('NLMS', 'ECA', 'Location', 'northwest');

    subplot(1, 2, 2);
    semilogx(MU_LIST, n_loss, '-s', 'LineWidth', 1); hold on;
    semilogx([MU_LIST(1) MU_LIST(end)], [3 3], 'r:');
    grid on;
    xlabel('NLMS step size \mu');
    ylabel('measured target loss (dB)');
    title('Target loss, unreliable where the depth is small');
    legend('NLMS', '3 dB limit', 'Location', 'northeast');
    save_fig(h, 'd4_eca_vs_nlms');
    close(h);
end

% -------------------------------------------------- the Doppler notch
% Measured on the 5 MHz channel to keep the run time down.  One target at a
% time, its Doppler swept, with the batch length as the parameter.  The
% notch has to be measured rather than assumed, because it is the reason
% AT5 asks for the target loss explicitly.
fprintf('\nDoppler notch of ECA, %d RB, %d ms record\n', 25, NSUBFRAMES);
set_seed(5003);
[wave2, ~, info2] = wave_cache(25, RHO, NCELLID, NSUBFRAMES);
k2 = ceil(CLUTTER_RANGE / C * info2.fs) + 16;
n2 = numel(wave2);
FD_LIST    = [0 10 20 40 60 100 200 400];
BATCH_MS   = [5 20];
notch = zeros(numel(BATCH_MS), numel(FD_LIST));

for b = 1:numel(BATCH_MS)
    nb2 = round(BATCH_MS(b) * 1e-3 * info2.fs);
    nb2 = min(nb2, n2);
    for i = 1:numel(FD_LIST)
        % Place a target at a fixed delay inside the cancelled delay span,
        % and force the Doppler by choosing the velocity along the bisector.
        pos = [15e3 3905];
        g   = geometry(cfg.tx_pos, cfg.rx_pos, pos(:), [0; 0], cfg.fc, info2.occupied_bw);
        v   = FD_LIST(i) * g.lambda / 2 / max(cos(g.beta / 2), eps);
        vel = (v * g.bisector).';

        c2 = cfg;
        c2.targets = struct('pos', pos, 'vel', vel);
        c2.seed = 5004;
        [ss2, sr2, tr2] = synth_rx(wave2, info2, c2);
        probe_rb = tr2.targets(1).Rb;
        sc2 = eca(ss2, sr2, k2, nb2);
        a = match_amp(sc2, sr2, info2.fs, tr2.targets(1).delay_samples, tr2.targets(1).fd);
        notch(b, i) = 20 * log10(tr2.targets(1).amp / abs(a));
    end
    fprintf('  batch %2d ms (notch about %4.0f Hz wide) : loss ', ...
        BATCH_MS(b), info2.fs / nb2);
    fprintf('%5.1f ', notch(b, :));
    fprintf(' dB\n');
end
fprintf('  target Doppler                              :       ');
fprintf('%5.0f ', FD_LIST);
fprintf(' Hz\n');
fprintf(['  The probe sits at %.0f m of bistatic range, %.0f samples, inside the\n' ...
         '  K = %d span the canceller covers.  A target beyond that span is not\n' ...
         '  touched at all, whatever its Doppler - which is why the two targets\n' ...
         '  in the sweep above lose almost nothing.\n'], ...
        probe_rb, probe_rb / C * info2.fs, k2);

if MAKE_FIGURES
    h = figure('Name', 'ECA Doppler notch');
    plot(FD_LIST, notch.', '-o', 'LineWidth', 1); hold on;
    plot([FD_LIST(1) FD_LIST(end)], [3 3], 'r:');
    grid on;
    xlabel('target bistatic Doppler (Hz)');
    ylabel('target loss caused by cancellation (dB)');
    title('Cancellation digs a notch around zero Doppler');
    lbl = arrayfun(@(m) sprintf('%d ms batch', m), BATCH_MS, 'UniformOutput', false);
    lbl{end + 1} = '3 dB limit';
    legend(lbl, 'Location', 'northeast');
    save_fig(h, 'd4_eca_doppler_notch');
    close(h);
end

fprintf('\nExperiment 3 complete.\n');
