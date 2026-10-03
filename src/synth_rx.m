function [s_surv, s_ref, truth] = synth_rx(wave, info, cfg)
%SYNTH_RX  Synthesise the surveillance and reference channels (module M4).
%
%   [S_SURV, S_REF, TRUTH] = SYNTH_RX(WAVE, INFO, CFG)
%
%   Builds the signal actually seen by a passive bistatic receiver:
%
%     s_surv(t) =  A0 * s_ref(t)                                direct path
%                + sum_i c_i * s_ref(t - tau_i)                 static clutter
%                + sum_k a_k * s_ref(t - tau_k) * e^{j2*pi*fd_k*t}  targets
%                + n(t)                                         noise
%
%   The reference channel is assumed ideal, S_REF = WAVE.  Reconstructing a
%   real reference channel is explicitly out of scope for this project and
%   that assumption has to be repeated in the report.
%
%   Amplitudes are all referred to the direct path, which is given unit
%   amplitude.  Target 1 sits DER_DB below it, and the noise is set from the
%   requested per sample signal to noise ratio of target 1.  Nothing here is
%   scaled by the processing gain; that gain appears later, in the CAF.
%
%   Inputs
%     WAVE  reference waveform from LTE_GENERATE (column vector)
%     INFO  the info struct from LTE_GENERATE (fs is the only field used)
%     CFG   configuration struct, every field optional:
%
%       fc          carrier frequency, Hz                      1800e6
%       tx_pos      illuminator position [x y], m              [0 0]
%       rx_pos      receiver position [x y], m                 [10e3 0]
%       targets     struct array, fields pos [x y], vel [vx vy]
%                   and optional rcs_db (relative, dB)         two defaults
%                   pass [] for a direct path only scenario
%       der_db      direct to echo ratio for target 1, dB      60
%       snr_db      per sample SNR of target 1 echo, dB        -20
%       n_clutter   number of static clutter scatterers        5
%       clutter_db  mean clutter to direct ratio per scatterer -25
%       clutter_max_range  largest clutter bistatic range, m   3000
%       direct_on   include the direct path                    true
%       clutter_on  include the clutter                        true
%       noise_on    include the noise                          true
%       seed        random seed, for repeatability             1
%
%   Outputs
%     S_SURV  surveillance channel, complex column vector
%     S_REF   reference channel, complex column vector
%     TRUTH   ground truth struct:
%       targets   struct array with pos, vel, Rb, tau, fd, delay_samples,
%                 amp, snr_db, and the full geometry output
%       clutter   struct array with tau, delay_samples, amp
%       a0        direct path amplitude
%       sigma2    noise power per sample
%       fs, fc, cfg, der_db_actual
%
%   Because the direct path defines the zero of the delay axis, target
%   delays in TRUTH are bistatic delays Rb/c, not absolute propagation
%   delays.  That is exactly what the range-Doppler map measures.
%
%   Example
%     [w, ~, info] = lte_generate(50, 1.0, 17, 20);
%     [ss, sr, truth] = synth_rx(w, info, struct('der_db', 60, 'snr_db', -25));
%
%   See also LTE_GENERATE, GEOMETRY, FRAC_DELAY, ECA, CAF.

if nargin < 3, cfg = struct(); end
if nargin < 2 || ~isfield(info, 'fs')
    error('synth_rx:info', 'info must be the struct returned by lte_generate.');
end

C  = 299792458;
fs = info.fs;

cfg = default_field(cfg, 'fc',          1800e6);
cfg = default_field(cfg, 'tx_pos',      [0 0]);
cfg = default_field(cfg, 'rx_pos',      [10e3 0]);
cfg = default_field(cfg, 'der_db',      60);
cfg = default_field(cfg, 'snr_db',      -20);
cfg = default_field(cfg, 'n_clutter',   5);
cfg = default_field(cfg, 'clutter_db',  -25);
cfg = default_field(cfg, 'clutter_max_range', 3000);
cfg = default_field(cfg, 'direct_on',   true);
cfg = default_field(cfg, 'clutter_on',  true);
cfg = default_field(cfg, 'noise_on',    true);
cfg = default_field(cfg, 'seed',        1);
if ~isfield(cfg, 'targets')
    cfg.targets = struct( ...
        'pos', {[5000 3000],   [7000 1500]}, ...
        'vel', {[  40   30],   [ -60   10]});
end
if cfg.direct_on && (cfg.der_db < 40 || cfg.der_db > 80)
    warning('synth_rx:der', ...
        'der_db = %g dB is outside the 40 to 80 dB range the PRD calls realistic.', cfg.der_db);
end
if cfg.n_clutter < 0 || cfg.n_clutter > 20
    error('synth_rx:nclutter', 'n_clutter must be between 0 and 20.');
end

set_seed(cfg.seed);

s_ref  = wave(:);
n      = numel(s_ref);
t      = (0:n - 1).' / fs;
s_surv = zeros(n, 1);

% ------------------------------------------------------- amplitude budget
a0 = 1;                                       % direct path is the reference
a1 = 10^(-cfg.der_db / 20);                   % target 1 echo amplitude
pref = mean(abs(s_ref).^2);
if isfield(cfg, 'dnr_db') && ~isempty(cfg.dnr_db)
    % Anchor on the direct path instead.  The two anchors are linked by
    % snr_db = dnr_db - der_db; the direct to noise ratio is what caps the
    % cancellation depth that M5 can reach, so it is the natural handle for
    % cancellation experiments.
    sigma2 = a0^2 * pref / 10^(cfg.dnr_db / 10);
    cfg.snr_db = cfg.dnr_db - cfg.der_db;
else
    sigma2 = a1^2 * pref / 10^(cfg.snr_db / 10);
end

% ------------------------------------------------------------ direct path
if cfg.direct_on
    s_surv = s_surv + a0 * s_ref;
end

% --------------------------------------------------------- static clutter
clutter = struct('tau', {}, 'delay_samples', {}, 'amp', {}, 'range_m', {});
if cfg.clutter_on && cfg.n_clutter > 0
    for i = 1:cfg.n_clutter
        rng_m = cfg.clutter_max_range * rand();
        tau   = rng_m / C;
        d     = tau * fs;
        % Log-uniform spread of +-5 dB about the nominal clutter level, so
        % the scatterers are not all identical.
        amp   = 10^((cfg.clutter_db + 10 * (rand() - 0.5)) / 20) * exp(2j * pi * rand());
        s_surv = s_surv + amp * frac_delay(s_ref, d);
        clutter(i).tau = tau;
        clutter(i).delay_samples = d;
        clutter(i).amp = amp;
        clutter(i).range_m = rng_m;
    end
end

% ----------------------------------------------------------- moving targets
nt = numel(cfg.targets);
if nt > 4
    error('synth_rx:ntargets', 'At most 4 moving targets are supported.');
end
targets = struct('pos', {}, 'vel', {}, 'Rb', {}, 'tau', {}, 'fd', {}, ...
                 'delay_samples', {}, 'amp', {}, 'snr_db', {}, 'geom', {});
for k = 1:nt
    tg = cfg.targets(k);
    g  = geometry(cfg.tx_pos, cfg.rx_pos, tg.pos(:), tg.vel(:), cfg.fc, info.occupied_bw);

    if isfield(tg, 'rcs_db') && ~isempty(tg.rcs_db)
        rel = tg.rcs_db;
    else
        rel = 0;
    end
    amp = a1 * 10^(rel / 20);
    d   = g.tau * fs;

    if d >= n
        error('synth_rx:toolong', ...
            'Target %d has a bistatic delay of %.0f samples, longer than the waveform.', k, d);
    end

    s_surv = s_surv + amp * frac_delay(s_ref, d) .* exp(2j * pi * g.fd * t);

    targets(k).pos = tg.pos;
    targets(k).vel = tg.vel;
    targets(k).Rb  = g.Rb;
    targets(k).tau = g.tau;
    targets(k).fd  = g.fd;
    targets(k).delay_samples = d;
    targets(k).amp = amp;
    targets(k).snr_db = 10 * log10(amp^2 * pref / sigma2);
    targets(k).geom = g;
end

% ------------------------------------------------------------------ noise
if cfg.noise_on && sigma2 > 0
    s_surv = s_surv + sqrt(sigma2 / 2) * (randn(n, 1) + 1j * randn(n, 1));
end

truth = struct();
truth.targets = targets;
truth.clutter = clutter;
truth.a0      = a0;
truth.a1      = a1;
truth.sigma2  = sigma2;
truth.fs      = fs;
truth.fc      = cfg.fc;
truth.pref    = pref;
truth.dnr_db  = 10 * log10(a0^2 * pref / sigma2);
if nt > 0
    truth.der_db_actual = 20 * log10(a0 / targets(1).amp);
else
    truth.der_db_actual = NaN;
end
truth.cfg     = cfg;
end

% ======================================================================

function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
