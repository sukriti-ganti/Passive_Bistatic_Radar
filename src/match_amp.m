function [a, replica] = match_amp(s, ref, fs, delay_samples, fd)
%MATCH_AMP  Least squares estimate of one echo amplitude from ground truth.
%
%   [A, REPLICA] = MATCH_AMP(S, REF, FS, DELAY_SAMPLES, FD)
%
%   Builds the replica REF delayed by DELAY_SAMPLES samples and shifted by
%   FD hertz, then returns the complex amplitude A that best explains S in a
%   least squares sense, A = (u'*s)/(u'*u).
%
%   Used by the acceptance tests to measure an echo amplitude without the
%   straddle loss that a gridded range-Doppler peak suffers from: the CAF
%   samples delay at 1/fs and Doppler at PRF/nseg, so a peak that falls
%   between bins reads one to two dB low.  A matched replica has no such
%   quantisation, which is why the direct to echo ratio is verified this way.
%
%   The estimate is only unbiased when the other components of S are
%   uncorrelated with the replica, so the tests apply it to scenarios that
%   contain a single echo at a time.
%
%   See also SYNTH_RX, FRAC_DELAY, TEST_AT4.

s = s(:);
ref = ref(:);
n = numel(s);
t = (0:n - 1).' / fs;
replica = frac_delay(ref, delay_samples) .* exp(2j * pi * fd * t);
a = (replica' * s) / (replica' * replica);
end
