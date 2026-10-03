function [chi, tau_axis, fd_axis] = caf(surv, ref, fs, nbatch, maxlag)
%CAF  Batched cross ambiguity function / range-Doppler map (modules M2, M6).
%
%   [CHI, TAU_AXIS, FD_AXIS] = CAF(SURV, REF, FS, NBATCH, MAXLAG)
%
%   Computes
%
%       chi(tau, f) = sum_n  surv(n) * conj(ref(n - tau)) * exp(-j2*pi*f*n/fs)
%
%   using the batched (segmented) algorithm:
%
%     1. both signals are split into NSEG batches of NBATCH samples
%     2. inside each batch the surveillance signal is correlated against the
%        reference over lags -MAXLAG .. +MAXLAG, by FFT multiplication
%        (xcorr lives in the Signal Processing Toolbox and is not used; the
%        FFT route is also considerably faster)
%     3. the resulting lag-by-batch matrix is Fourier transformed across the
%        batch index to resolve Doppler
%     4. the Doppler axis is fftshifted so that zero Doppler sits in the
%        middle of the output
%
%   Inputs
%     SURV    surveillance channel samples (vector)
%     REF     reference channel samples (vector)
%     FS      sample rate in Hz
%     NBATCH  samples per batch.  Sets the pulse repetition frequency of the
%             equivalent pulsed radar, PRF = FS/NBATCH, and therefore the
%             unambiguous Doppler span +-PRF/2
%     MAXLAG  largest lag in samples.  Bistatic range covered is
%             MAXLAG*c/FS metres
%
%   Outputs
%     CHI       (2*MAXLAG+1) by NSEG complex matrix, delay down, Doppler across
%     TAU_AXIS  1 by (2*MAXLAG+1) delay axis in seconds.  Bistatic range is
%               TAU_AXIS * 299792458
%     FD_AXIS   1 by NSEG Doppler axis in Hz, already fftshifted
%
%   Accuracy note.  The batched algorithm treats the Doppler phase as
%   constant inside a batch.  The resulting amplitude error is the factor
%   sinc(fd/PRF), which is under 1 percent for |fd| < PRF/12 and reaches
%   -3.9 dB at the edge of the unambiguous band.  Choose NBATCH so that the
%   Doppler of interest sits well inside the band.  CAF_DIRECT computes the
%   same quantity by brute force and is used by test_at2 to verify this
%   implementation.
%
%   Example
%     [chi, tau, fd] = caf(s_surv, s_ref, info.fs, 15360, 512);
%     imagesc(fd, tau*3e8, 20*log10(abs(chi)/max(abs(chi(:)))));
%
%   See also CAF_DIRECT, AF_METRICS, LTE_GENERATE.

if nargin < 5
    error('caf:nargin', 'Five input arguments are required.');
end
surv = surv(:);
ref  = ref(:);
if ~isscalar(nbatch) || nbatch < 1 || nbatch ~= fix(nbatch)
    error('caf:nbatch', 'nbatch must be a positive integer.');
end
if ~isscalar(maxlag) || maxlag < 0 || maxlag ~= fix(maxlag)
    error('caf:maxlag', 'maxlag must be a non-negative integer.');
end
if maxlag >= nbatch
    error('caf:maxlag', 'maxlag (%d) must be smaller than nbatch (%d).', maxlag, nbatch);
end

n    = min(numel(surv), numel(ref));
nseg = floor(n / nbatch);
if nseg < 1
    error('caf:short', ...
        'Signals hold %d samples, at least one batch of %d is needed.', n, nbatch);
end

% Correlation FFT length.  Nc >= nbatch + maxlag + 1 keeps the wrapped
% negative lags clear of the positive lags we extract.
nc = 2^nextpow2(nbatch + maxlag + 1);

% Batch the signals.  Columns are batches.
sm = reshape(surv(1:nseg * nbatch), nbatch, nseg);
rm = reshape(ref(1:nseg * nbatch),  nbatch, nseg);

% Cross correlate every batch, in blocks so that the temporary nc-by-block
% matrices stay small enough for a 4 GB laptop.
chi = zeros(2 * maxlag + 1, nseg);
block = max(1, floor(4e6 / nc));
for a = 1:block:nseg
    b = min(a + block - 1, nseg);
    s = fft(sm(:, a:b), nc, 1);
    r = fft(rm(:, a:b), nc, 1);
    c = ifft(s .* conj(r), [], 1);
    chi(:, a:b) = [c(nc - maxlag + 1:nc, :); c(1:maxlag + 1, :)];
end

% Doppler transform across the batch index.
chi = fftshift(fft(chi, nseg, 2), 2);

% Axes.
tau_axis = (-maxlag:maxlag) / fs;
prf      = fs / nbatch;
fd_axis  = ((0:nseg - 1) - floor(nseg / 2)) * (prf / nseg);
end
