function [det_mask, thr, out] = cfar2d(pow, params)
%CFAR2D  Two dimensional cell averaging CFAR detector (module M7).
%
%   [DET_MASK, THR, OUT] = CFAR2D(POW)
%   [DET_MASK, THR, OUT] = CFAR2D(POW, PARAMS)
%
%   POW is a power map, normally ABS(CHI).^2 from CAF, with bistatic range
%   down the rows and Doppler across the columns.  For every cell under test
%   the detector
%
%     1. excludes a ring of guard cells immediately around it, so that a
%        strong target cannot raise its own threshold
%     2. averages the power over the ring of training cells outside the
%        guard ring
%     3. multiplies that average by alpha = N*(Pfa^(-1/N) - 1), the factor
%        that holds the false alarm rate at Pfa for exponentially
%        distributed noise power and N training cells
%     4. declares a detection where the cell under test exceeds it
%
%   PARAMS fields, all optional
%     guard_range   guard cells either side in range     2
%     guard_dopp    guard cells either side in Doppler   2
%     train_range   training cells either side in range  8
%     train_dopp    training cells either side in Doppler 8
%     pfa           design probability of false alarm    1e-4
%
%   Outputs
%     DET_MASK  logical, same size as POW, true where a detection is declared
%     THR       the threshold applied, NaN in the border where the full
%               window does not fit
%     OUT       struct with
%                 alpha       threshold multiplier used
%                 n_train     number of training cells N
%                 tested      logical mask of cells that were tested
%                 n_tested    number of tested cells, the denominator for a
%                             measured false alarm rate
%                 n_det       number of detections
%                 noise_est   the cell averaged noise estimate
%                 params      the settings actually used
%
%   Border cells cannot be tested because their training ring falls off the
%   map; they are reported as not detected and excluded from N_TESTED, so a
%   measured false alarm rate is N_DET/N_TESTED and not N_DET/NUMEL(POW).
%
%   The window sums are computed from a summed area table, so the cost is
%   independent of the window size: two cumsums and eight lookups per cell.
%
%   Example
%     [d, thr, s] = cfar2d(abs(chi).^2, struct('pfa', 1e-4));
%     fprintf('%d detections in %d tested cells\n', s.n_det, s.n_tested);
%
%   See also CAF, ECA.

if nargin < 2, params = struct(); end
params = default_field(params, 'guard_range', 2);
params = default_field(params, 'guard_dopp',  2);
params = default_field(params, 'train_range', 8);
params = default_field(params, 'train_dopp',  8);
params = default_field(params, 'pfa',         1e-4);

gr = params.guard_range;
gd = params.guard_dopp;
tr = params.train_range;
td = params.train_dopp;
pfa = params.pfa;

if pfa <= 0 || pfa >= 1
    error('cfar2d:pfa', 'pfa must lie strictly between 0 and 1.');
end
if tr < 1 && td < 1
    error('cfar2d:train', 'At least one training dimension must be non-zero.');
end

hr = gr + tr;
hd = gd + td;
[nr, nd] = size(pow);
if nr < 2 * hr + 1 || nd < 2 * hd + 1
    error('cfar2d:size', ...
        'Map is %d by %d, too small for a %d by %d CFAR window.', ...
        nr, nd, 2 * hr + 1, 2 * hd + 1);
end

n_train = (2 * hr + 1) * (2 * hd + 1) - (2 * gr + 1) * (2 * gd + 1);
alpha   = n_train * (pfa^(-1 / n_train) - 1);

% Summed area table.  ii(i+1,j+1) is the sum of pow(1:i, 1:j).
ii = zeros(nr + 1, nd + 1);
ii(2:end, 2:end) = cumsum(cumsum(pow, 1), 2);

r = (hr + 1):(nr - hr);
c = (hd + 1):(nd - hd);

outer = ii(r + hr + 1, c + hd + 1) - ii(r - hr, c + hd + 1) ...
      - ii(r + hr + 1, c - hd)     + ii(r - hr, c - hd);
inner = ii(r + gr + 1, c + gd + 1) - ii(r - gr, c + gd + 1) ...
      - ii(r + gr + 1, c - gd)     + ii(r - gr, c - gd);

noise_est = (outer - inner) / n_train;
thr_in    = alpha * noise_est;

thr = nan(nr, nd);
thr(r, c) = thr_in;

tested = false(nr, nd);
tested(r, c) = true;

det_mask = false(nr, nd);
det_mask(r, c) = pow(r, c) > thr_in;

out = struct();
out.alpha     = alpha;
out.n_train   = n_train;
out.tested    = tested;
out.n_tested  = numel(thr_in);
out.n_det     = sum(det_mask(:));
out.noise_est = nan(nr, nd);
out.noise_est(r, c) = noise_est;
out.params    = params;
end

% ======================================================================

function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
