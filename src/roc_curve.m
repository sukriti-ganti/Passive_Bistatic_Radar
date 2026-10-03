function out = roc_curve(cfg)
%ROC_CURVE  Monte Carlo probability of detection against SNR (module M7).
%
%   OUT = ROC_CURVE(CFG)
%
%   Measures the detection performance of the CFAR detector of M7 operating
%   on the range-Doppler map of M6, over the LTE waveform of M1.  Nothing is
%   modelled analytically: the target echo is synthesised by M4 at the true
%   bistatic delay and Doppler, the map is formed by CAF, and CFAR2D decides.
%
%   CFG fields, all optional
%     nrb, rho, ncellid       waveform, default 50, 1.0, 17
%     nsubframes              record length in ms, default 20
%     nbatch                  CAF batch length in samples, default fs*0.5e-3
%     max_range_m             bistatic range span of the map, default 8 km.
%                             Converted to a lag count at the sample rate;
%                             it must cover the target or nothing is found
%     snr_db                  vector of in-band pre-detection SNRs to sweep,
%                             default -60:5:-25
%     ntrials                 Monte Carlo trials per SNR, default 200
%     cfar                    params struct passed to CFAR2D
%     tx_pos, rx_pos, fc      geometry, defaults [0 0], [30e3 0], 1800e6
%     target                  struct with pos and vel, default one crossing
%                             target at (12 km, 9 km)
%     seed                    default 1
%
%   OUT fields
%     snr_db        the SNR axis, as requested
%     snr_post_db   the same SNRs after the coherent processing gain B*T
%     pd            measured probability of detection at the target cell
%     pfa_meas      measured false alarm rate over every tested cell of every
%                   trial, excluding a guard box around the target
%     pfa_by_snr    the same, resolved per SNR point
%     gp_db         processing gain 10*log10(B*T)
%     n_tested, alpha, n_train, target, info, cfg
%
%   Definition of SNR.  The x axis is the pre-detection signal to noise
%   ratio measured in the occupied bandwidth B, not in the sample rate.
%   That is the convention the link budget of M8 uses, so a point on this
%   curve can be read straight into the feasibility table: a target is
%   detected when its in-band SNR plus the processing gain 10*log10(B*T)
%   exceeds the required post-processing SNR.  The receiver oversamples,
%   fs > B, so the per sample SNR used inside SYNTH_RX is lower than this by
%   10*log10(fs/B), and the code converts between them explicitly.
%
%   Speed.  The CAF is linear in the surveillance signal, so the noise free
%   target map is computed once and only the noise map is recomputed for
%   each trial.  That is exact, not an approximation, and it makes a 200
%   trial sweep take seconds instead of minutes.
%
%   Example
%     out = roc_curve(struct('nrb', 50, 'ntrials', 100));
%     plot(out.snr_db, out.pd);
%
%   See also CFAR2D, CAF, SYNTH_RX, LINK_BUDGET.

if nargin < 1, cfg = struct(); end
cfg = default_field(cfg, 'nrb',        50);
cfg = default_field(cfg, 'rho',        1.0);
cfg = default_field(cfg, 'ncellid',    17);
cfg = default_field(cfg, 'nsubframes', 20);
cfg = default_field(cfg, 'max_range_m', 8e3);
cfg = default_field(cfg, 'snr_db',     -60:5:-25);
cfg = default_field(cfg, 'ntrials',    200);
cfg = default_field(cfg, 'cfar',       struct('pfa', 1e-4));
cfg = default_field(cfg, 'tx_pos',     [0 0]);
cfg = default_field(cfg, 'rx_pos',     [30e3 0]);
cfg = default_field(cfg, 'fc',         1800e6);
cfg = default_field(cfg, 'target',     struct('pos', [12e3 9e3], 'vel', [80 -40]));
cfg = default_field(cfg, 'seed',       1);

set_seed(cfg.seed);
[wave, ~, info] = lte_generate(cfg.nrb, cfg.rho, cfg.ncellid, cfg.nsubframes);
cfg = default_field(cfg, 'nbatch', round(info.fs * 0.5e-3));
if ~isfield(cfg, 'maxlag') || isempty(cfg.maxlag)
    cfg.maxlag = ceil(cfg.max_range_m / 299792458 * info.fs);
end

bw   = info.occupied_bw;
tint = cfg.nsubframes * 1e-3;
gp   = bw * tint;

% Noise free target echo at unit amplitude.  der_db = 0 with the direct path
% switched off gives an echo of amplitude 1.
scfg = struct('fc', cfg.fc, 'tx_pos', cfg.tx_pos, 'rx_pos', cfg.rx_pos, ...
              'der_db', 0, 'direct_on', false, 'clutter_on', false, ...
              'noise_on', false, 'targets', cfg.target, 'seed', cfg.seed);
[s_t, s_ref, truth] = synth_rx(wave, info, scfg);
t = truth.targets(1);

if t.tau > cfg.maxlag / info.fs
    error('roc_curve:range', ...
        ['The target is at %.0f m bistatic range but the map only reaches %.0f m. ' ...
         'Raise max_range_m or move the target.'], ...
        t.Rb, cfg.maxlag / info.fs * 299792458);
end

[chi_sig, tau, fd] = caf(s_t, s_ref, info.fs, cfg.nbatch, cfg.maxlag);
[~, trow] = min(abs(tau - t.tau));
[~, tcol] = min(abs(fd - t.fd));

n     = numel(s_ref);
pref  = truth.pref;
nsnr  = numel(cfg.snr_db);
pd    = zeros(1, nsnr);
pfa   = zeros(1, nsnr);
n_det_fa = 0;
n_tested_total = 0;
alpha = NaN; n_train = NaN; n_tested = NaN;

for i = 1:nsnr
    % In-band SNR gamma means noise power per sample sigma2 such that
    % pref / (sigma2 * B/fs) = gamma.
    gamma  = 10^(cfg.snr_db(i) / 10);
    sigma2 = pref * info.fs / (bw * gamma);

    hits = 0;
    fa_i = 0;
    tested_i = 0;
    for trial = 1:cfg.ntrials
        noise = sqrt(sigma2 / 2) * (randn(n, 1) + 1j * randn(n, 1));
        chi_n = caf(noise, s_ref, info.fs, cfg.nbatch, cfg.maxlag);
        p = abs(chi_sig + chi_n).^2;

        [det_mask, ~, ci] = cfar2d(p, cfg.cfar);
        alpha = ci.alpha; n_train = ci.n_train; n_tested = ci.n_tested;

        hits = hits + double(det_mask(trow, tcol));

        % False alarms: every tested cell except a guard box round the target.
        fa = det_mask & ci.tested;
        r0 = max(1, trow - 6); r1 = min(size(p, 1), trow + 6);
        c0 = max(1, tcol - 3); c1 = min(size(p, 2), tcol + 3);
        fa(r0:r1, c0:c1) = false;
        fa_i     = fa_i + sum(fa(:));
        tested_i = tested_i + ci.n_tested - (r1 - r0 + 1) * (c1 - c0 + 1);
    end
    pd(i)  = hits / cfg.ntrials;
    pfa(i) = fa_i / max(tested_i, 1);
    n_det_fa = n_det_fa + fa_i;
    n_tested_total = n_tested_total + tested_i;
end

out = struct();
out.snr_db      = cfg.snr_db;
out.gp_db       = 10 * log10(gp);
out.snr_post_db = cfg.snr_db + out.gp_db;
out.pd          = pd;
out.pfa_meas    = n_det_fa / max(n_tested_total, 1);
out.pfa_by_snr  = pfa;
out.alpha       = alpha;
out.n_train     = n_train;
out.n_tested    = n_tested;
out.target      = t;
out.tau         = tau;
out.fd          = fd;
out.target_cell = [trow, tcol];
out.info        = info;
out.tint_s      = tint;
out.bw_hz       = bw;
out.cfg         = cfg;
end

% ======================================================================

function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
