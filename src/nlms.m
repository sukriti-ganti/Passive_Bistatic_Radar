function [s_clean, w, info] = nlms(s_surv, s_ref, k_order, mu, delta)
%NLMS  Normalised least mean squares direct signal canceller (module M5).
%
%   [S_CLEAN, W, INFO] = NLMS(S_SURV, S_REF, K_ORDER)
%   [S_CLEAN, W, INFO] = NLMS(S_SURV, S_REF, K_ORDER, MU, DELTA)
%
%   The simple alternative to ECA that the PRD asks for as a comparison.
%   Where ECA solves the normal equations in one shot per batch, NLMS walks
%   the weights downhill one sample at a time:
%
%       u(n) = [s_ref(n) ... s_ref(n-K)]'
%       e(n) = s_surv(n) - w'*u(n)
%       w    = w + mu * u(n) * conj(e(n)) / (delta + u(n)'*u(n))
%
%   S_CLEAN is the error sequence e(n).
%
%   Inputs
%     K_ORDER  filter order, K+1 taps
%     MU       step size, 0 < mu < 2 for stability.  Default 0.5.  Larger mu
%              converges faster and leaves a larger misadjustment floor,
%              which is the whole trade-off
%     DELTA    regularisation in the denominator.  Default 1e-6
%
%   INFO carries depth_db, dsi_residual_db, converge_samples (the index at
%   which the running error power first falls within 3 dB of its final
%   value) and the settings used.  The convergence transient at the start of
%   the record is the practical price of NLMS: the first few thousand
%   samples are barely cancelled at all, and they still contain a full
%   strength direct path that will smear across the range-Doppler map.
%
%   NLMS is O(n*K) and runs sample by sample, so on a 20 ms record it is
%   noticeably slower in MATLAB than the vectorised ECA despite the lower
%   operation count.
%
%   See also ECA.

if nargin < 3
    error('nlms:nargin', 'At least three input arguments are required.');
end
if nargin < 4 || isempty(mu),    mu = 0.5;    end
if nargin < 5 || isempty(delta), delta = 1e-6; end
if mu <= 0 || mu >= 2
    error('nlms:mu', 'mu must lie in (0, 2) for stability, got %g.', mu);
end

s_surv = s_surv(:);
s_ref  = s_ref(:);
n = min(numel(s_surv), numel(s_ref));
s_surv = s_surv(1:n);
s_ref  = s_ref(1:n);

kk = k_order + 1;
w  = zeros(kk, 1);
s_clean = zeros(n, 1);
refpad  = [zeros(k_order, 1); s_ref];

for i = 1:n
    u = refpad(i + k_order:-1:i);         % [ref(i) ref(i-1) ... ref(i-K)]
    e = s_surv(i) - w' * u;
    w = w + (mu / (delta + real(u' * u))) * u * conj(e);
    s_clean(i) = e;
end

pin  = mean(abs(s_surv).^2);
pout = mean(abs(s_clean).^2);

% Convergence point: the first sample at which a full length moving average
% of the error power has come within 3 dB of the steady state level.  The
% search starts at win, because before that the moving average is still
% filling and reads artificially low.
win  = max(1, floor(n / 100));
pr   = filter(ones(win, 1) / win, 1, abs(s_clean).^2);
tail = mean(abs(s_clean(max(1, end - 10 * win):end)).^2);
conv_idx = find(pr(win:end) <= 2 * tail, 1);
if isempty(conv_idx)
    conv_idx = n;
else
    conv_idx = conv_idx + win - 1;
end

info = struct();
info.depth_db        = 10 * log10(pin / pout);
info.dsi_residual_db = 20 * log10( ...
    max(abs(s_ref' * s_clean), eps * abs(s_ref' * s_surv)) / abs(s_ref' * s_surv));
info.converge_samples = conv_idx;
info.k_order = k_order;
info.mu      = mu;
info.delta   = delta;
info.input_power  = pin;
info.output_power = pout;
end
