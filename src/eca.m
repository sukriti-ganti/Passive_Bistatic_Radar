function [s_clean, w, info] = eca(s_surv, s_ref, k_order, nbatch, loading)
%ECA  Extensive Cancellation Algorithm, batched (module M5).
%
%   [S_CLEAN, W, INFO] = ECA(S_SURV, S_REF, K_ORDER)
%   [S_CLEAN, W, INFO] = ECA(S_SURV, S_REF, K_ORDER, NBATCH, LOADING)
%
%   The direct path from the illuminator arrives 60 to 90 dB above any
%   target echo.  Its ambiguity sidelobes sit around -24 dB for a fully
%   loaded LTE signal, so untreated they bury every target on the map.  ECA
%   removes it, together with any static clutter, by projecting the
%   surveillance signal onto the orthogonal complement of the subspace
%   spanned by delayed copies of the reference:
%
%       X = [s_ref(t), s_ref(t-1), ... , s_ref(t-K)]
%       w = (X'*X) \ (X'*s_surv)
%       s_clean = s_surv - X*w
%
%   Everything at zero Doppler and at a delay of 0..K samples is cancelled.
%   A moving target survives because its Doppler shift makes it (nearly)
%   orthogonal to the zero Doppler subspace, which is also why the target
%   loss has to be measured rather than assumed - see test_at5.
%
%   Inputs
%     S_SURV   surveillance channel
%     S_REF    reference channel
%     K_ORDER  filter order K.  K+1 delayed copies, delays 0..K samples.
%              K must cover the longest clutter delay of interest:
%              K = ceil(range_m * fs / c)
%     NBATCH   samples per batch.  Default: the whole signal in one batch.
%              Batching keeps the (K+1)-by-(K+1) normal equations small and
%              lets the canceller track slow changes, at the cost of some
%              Doppler resolution inside each batch
%     LOADING  diagonal loading, as a fraction of the mean eigenvalue.
%              Default 1e-10, enough to keep an ill conditioned batch
%              solvable without measurably changing the answer
%
%   Outputs
%     S_CLEAN  cancelled surveillance signal, same length as the input
%     W        (K+1) by NBATCHES matrix of cancellation weights
%     INFO     struct with
%                depth_db        total power reduction, dB.  This is the
%                                headline cancellation depth
%                dsi_residual_db residual correlation with the reference,
%                                relative to the input.  ECA drives this to
%                                the numerical floor because the projection
%                                is exact; NLMS does not, which is the point
%                                of comparing them
%                nbatches, k_order, nbatch, loading
%                input_power, output_power
%
%   How deep can it go.  The residual still contains the receiver noise, so
%   the total power reduction can never exceed the direct to noise ratio of
%   the input.  A K sweep therefore rises, then flattens at that ceiling -
%   which is exactly the shape test_at5 looks for.  It also cannot exceed
%   what the filter length covers: K must reach past the longest clutter
%   delay or the uncancelled clutter sets the floor instead.
%
%   Doppler notch.  The weights are fitted to data that contains the target
%   as well as the direct path, so the fit absorbs part of any echo that
%   looks stationary over a batch.  The notch that results is about
%   +-1/(NBATCH/fs) wide, and run_experiment3_cancellation measures it: with
%   a 20 ms batch a target above roughly 50 Hz of Doppler is untouched,
%   while one below that loses several dB.  Keep the batch long compared
%   with 1/fd of the slowest target of interest.  The default, one batch
%   over the whole record, is usually right for a simulation of a few tens
%   of milliseconds.
%
%   Complexity is O(n*K^2 + nbatches*K^3).  Doubling K roughly quadruples
%   the run time, which is why test_at5 sweeps K and looks for the knee
%   rather than simply using the largest K affordable.  The normal equations
%   are accumulated over row chunks so that the delayed copy matrix never
%   has to exist in memory all at once.
%
%   Example
%     [s_clean, ~, ci] = eca(s_surv, s_ref, 64, info.samples_per_subframe);
%     fprintf('cancellation depth %.1f dB\n', ci.depth_db);
%
%   See also NLMS, SYNTH_RX, CAF.

if nargin < 3
    error('eca:nargin', 'At least three input arguments are required.');
end
s_surv = s_surv(:);
s_ref  = s_ref(:);
n = min(numel(s_surv), numel(s_ref));
s_surv = s_surv(1:n);
s_ref  = s_ref(1:n);

if nargin < 4 || isempty(nbatch), nbatch = n; end
if nargin < 5 || isempty(loading), loading = 1e-10; end
if k_order < 0 || k_order ~= fix(k_order)
    error('eca:order', 'k_order must be a non-negative integer.');
end
if nbatch <= k_order
    error('eca:nbatch', ...
        'nbatch (%d) must exceed the filter order (%d), otherwise the normal equations are rank deficient.', ...
        nbatch, k_order);
end

kk   = k_order + 1;
nbat = max(1, floor(n / nbatch));

% Zero padded reference so that delayed copies can reach back before the
% start of a batch.  refpad(i + k_order) == s_ref(i).
refpad = [zeros(k_order, 1); s_ref];

% Row chunk size that keeps the temporary delayed copy matrix near 32 MB.
chunk = max(kk, floor(2^21 / kk));

s_clean = s_surv;
w = zeros(kk, nbat);

for b = 1:nbat
    i0 = (b - 1) * nbatch + 1;
    if b == nbat
        i1 = n;                      % last batch absorbs the remainder
    else
        i1 = b * nbatch;
    end

    % --- normal equations, accumulated chunk by chunk -------------------
    r = zeros(kk, kk);
    p = zeros(kk, 1);
    for c0 = i0:chunk:i1
        c1  = min(c0 + chunk - 1, i1);
        idx = (c0:c1).';
        x   = zeros(numel(idx), kk);
        for k = 0:k_order
            x(:, k + 1) = refpad(idx + k_order - k);
        end
        r = r + x' * x;
        p = p + x' * s_surv(idx);
    end
    r = r + (loading * real(trace(r)) / kk) * eye(kk);

    wb = r \ p;
    w(:, b) = wb;

    % --- residual, as an FIR filtering of the reference -----------------
    % filter() with the batch's leading k_order samples of history is
    % identical to forming X*wb, without ever building X.
    ext = refpad(i0:i1 + k_order);
    y   = filter(wb, 1, ext);
    s_clean(i0:i1) = s_surv(i0:i1) - y(k_order + 1:end);
end

pin  = mean(abs(s_surv).^2);
pout = mean(abs(s_clean).^2);

info = struct();
info.depth_db     = 10 * log10(pin / pout);
info.dsi_residual_db = 20 * log10( ...
    max(abs(s_ref' * s_clean), eps * abs(s_ref' * s_surv)) / abs(s_ref' * s_surv));
info.nbatches     = nbat;
info.k_order      = k_order;
info.nbatch       = nbatch;
info.loading      = loading;
info.input_power  = pin;
info.output_power = pout;
end
