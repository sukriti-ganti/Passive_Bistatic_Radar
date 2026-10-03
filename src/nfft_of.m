function n = nfft_of(nrb)
%NFFT_OF  FFT length for an LTE downlink bandwidth.
%
%   N = NFFT_OF(NRB) returns 128, 512, 1024 or 2048 for NRB of 6, 25, 50 or
%   100.  The sample rate is always N*15000.  Handy where a script needs the
%   sample rate before it has generated a waveform.

switch nrb
    case 6,   n = 128;
    case 25,  n = 512;
    case 50,  n = 1024;
    case 100, n = 2048;
    otherwise
        error('nfft_of:nrb', 'nrb must be one of 6, 25, 50, 100 (got %g).', nrb);
end
end
