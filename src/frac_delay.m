function y = frac_delay(x, d)
%FRAC_DELAY  Delay a signal by a fractional number of samples.
%
%   Y = FRAC_DELAY(X, D) delays X by D samples, where D need not be an
%   integer.  The delay is applied as a linear phase ramp in the frequency
%   domain, which is the ideal band limited (sinc) interpolator:
%
%       Y(f) = X(f) * exp(-j*2*pi*f*D)
%
%   Rounding a target delay to the nearest sample would quantise bistatic
%   range to c/(2*fs), about 10 m at 15.36 MHz, and would put an artificial
%   staircase into every range-Doppler measurement.  Hence this routine.
%
%   X is zero padded by the delay plus a 64 sample margin before the ramp is
%   applied, so the interpolator tails run into the padding instead of
%   wrapping onto the start of the signal.  The output is truncated back to
%   the length of X.  For even transform lengths the Nyquist bin is forced
%   real, otherwise the output of a real input would pick up a spurious
%   imaginary part.
%
%   Edge accuracy.  The ideal interpolator has infinite support, so a sample
%   near either end of a finite record is reconstructed from data that does
%   not exist.  The resulting error decays as 1/(pi*distance from the edge):
%   on an 8192 sample tone it is about 2e-2 ten samples in, 1e-3 two hundred
%   samples in and 1e-6 in the middle.  For this project that is far below
%   anything that matters, because the error scales with the echo being
%   delayed and the echoes are already 60 dB below the direct path, but it
%   is the reason the acceptance test checks the interior of the record.
%
%   Example
%     y = frac_delay(x, 12.37);
%
%   See also SYNTH_RX.

if nargin < 2
    error('frac_delay:nargin', 'Two input arguments are required.');
end
x = x(:);
n = numel(x);
if d == 0
    y = x;
    return
end

pad = ceil(abs(d)) + 64;
np  = n + pad;
np  = np + mod(-np, 64);                 % round up to a friendly FFT length

xf = fft(x, np);
k  = [0:ceil(np / 2) - 1, -floor(np / 2):-1].';
h  = exp(-2j * pi * k * d / np);
if mod(np, 2) == 0
    h(np / 2 + 1) = cos(pi * d);         % keep the Nyquist bin real
end

y = ifft(xf .* h);
y = y(1:n);
end
