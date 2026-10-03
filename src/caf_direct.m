function [chi, tau_axis, fd_axis] = caf_direct(surv, ref, fs, nbatch, maxlag)
%CAF_DIRECT  Brute force cross ambiguity function, validation reference only.
%
%   [CHI, TAU_AXIS, FD_AXIS] = CAF_DIRECT(SURV, REF, FS, NBATCH, MAXLAG)
%
%   Same arguments, same output convention and same mathematical definition
%   as CAF, but every correlation lag is accumulated with an explicit loop
%   and the Doppler transform is an explicit discrete Fourier sum.  Neither
%   fft nor fftshift is used anywhere in this file, so agreement between
%   CAF and CAF_DIRECT is genuine evidence that the fast implementation,
%   its lag extraction and its Doppler axis are all correct.
%
%   This routine is O(nseg * nbatch * maxlag * nseg) and must only ever be
%   run on a short segment.  A few thousand samples takes a second, a full
%   10 ms waveform would take hours.  Production code calls CAF.
%
%   See also CAF, TEST_AT2.

if nargin < 5
    error('caf_direct:nargin', 'Five input arguments are required.');
end
surv = surv(:);
ref  = ref(:);
n    = min(numel(surv), numel(ref));
nseg = floor(n / nbatch);
if nseg < 1
    error('caf_direct:short', 'At least one full batch is required.');
end
if maxlag >= nbatch
    error('caf_direct:maxlag', 'maxlag must be smaller than nbatch.');
end

lags = -maxlag:maxlag;
c    = zeros(numel(lags), nseg);          % lag by batch, before Doppler

for m = 1:nseg
    off = (m - 1) * nbatch;
    for li = 1:numel(lags)
        tau = lags(li);
        % Sample range for which both surv(n) and ref(n-tau) fall inside
        % this batch.  Identical support to the FFT correlation in CAF.
        p0 = max(1, 1 + tau);
        p1 = min(nbatch, nbatch + tau);
        acc = 0;
        for p = p0:p1
            acc = acc + surv(off + p) * conj(ref(off + p - tau));
        end
        c(li, m) = acc;
    end
end

% Explicit DFT across batches, written straight into fftshifted order.
chi = zeros(numel(lags), nseg);
m   = 0:(nseg - 1);
for i = 1:nseg
    k = i - 1 - floor(nseg / 2);          % signed Doppler bin index
    chi(:, i) = c * exp(-2j * pi * k * m(:) / nseg);
end

tau_axis = lags / fs;
prf      = fs / nbatch;
fd_axis  = ((0:nseg - 1) - floor(nseg / 2)) * (prf / nseg);
end
