function plot_rdmap(chi, tau_axis, fd_axis, ttl, dyn_db, ref_peak)
%PLOT_RDMAP  Draw a range-Doppler map or an ambiguity surface.
%
%   PLOT_RDMAP(CHI, TAU_AXIS, FD_AXIS, TTL)
%   PLOT_RDMAP(CHI, TAU_AXIS, FD_AXIS, TTL, DYN_DB, REF_PEAK)
%
%   Plots 20*log10(|CHI|) normalised to its own peak, with bistatic range in
%   kilometres on the vertical axis and Doppler in hertz on the horizontal.
%   DYN_DB sets the displayed dynamic range below the peak, default 60 dB.
%   REF_PEAK normalises against another map's peak instead of this one's,
%   which is how a before and after cancellation pair is drawn on one scale.

if nargin < 4, ttl = ''; end
if nargin < 5 || isempty(dyn_db), dyn_db = 60; end

a = abs(chi);
if nargin < 6 || isempty(ref_peak)
    ref_peak = max(a(:));
end
z = 20 * log10(a / ref_peak + eps);

imagesc(fd_axis, tau_axis * 299792458 / 1e3, z);
set(gca, 'YDir', 'normal');
caxis([-dyn_db 0]);
colorbar;
xlabel('bistatic Doppler (Hz)');
ylabel('bistatic range (km)');
title(ttl);
end
