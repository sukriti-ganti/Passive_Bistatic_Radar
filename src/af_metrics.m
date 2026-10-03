function m = af_metrics(chi, tau_axis, fd_axis, excl_delay, excl_dopp)
%AF_METRICS  Peak sidelobe and integrated sidelobe level of an ambiguity surface.
%
%   M = AF_METRICS(CHI, TAU_AXIS, FD_AXIS)
%   M = AF_METRICS(CHI, TAU_AXIS, FD_AXIS, EXCL_DELAY, EXCL_DOPP)
%
%   CHI is the complex ambiguity surface produced by CAF, delay down and
%   Doppler across.  The main lobe is taken to be the rectangle of
%   +-EXCL_DELAY delay cells by +-EXCL_DOPP Doppler cells centred on the
%   global peak; everything else is sidelobe.  Defaults are 2 and 2, which
%   is one resolution cell either side for an LTE signal sampled at the
%   Nyquist rate of its own occupied bandwidth.
%
%   Returned fields
%     psl_db      peak sidelobe level, 10*log10(max sidelobe power / peak power)
%     isl_db      integrated sidelobe level,
%                 10*log10(total sidelobe power / total main lobe power)
%     peak_power  |chi|^2 at the peak
%     peak_tau    delay of the peak, seconds
%     peak_fd     Doppler of the peak, Hz
%     peak_row, peak_col   indices of the peak
%     psl_tau, psl_fd      location of the worst sidelobe
%     mainlobe    logical mask of the excluded region
%
%   Both metrics are ratios, so they are invariant to any scaling of the
%   input signals.
%
%   See also CAF.

if nargin < 4 || isempty(excl_delay), excl_delay = 2; end
if nargin < 5 || isempty(excl_dopp),  excl_dopp  = 2; end

p = abs(chi).^2;
[nd, nf] = size(p);

[peak_power, pk] = max(p(:));
[pr, pc] = ind2sub([nd, nf], pk);

rows = max(1, pr - excl_delay):min(nd, pr + excl_delay);
cols = max(1, pc - excl_dopp):min(nf, pc + excl_dopp);

mainlobe = false(nd, nf);
mainlobe(rows, cols) = true;

side = p;
side(mainlobe) = 0;

[psl_power, sk] = max(side(:));
[sr, sc] = ind2sub([nd, nf], sk);

m = struct();
m.peak_power = peak_power;
m.peak_row   = pr;
m.peak_col   = pc;
m.peak_tau   = tau_axis(pr);
m.peak_fd    = fd_axis(pc);
m.psl_db     = 10 * log10(psl_power / peak_power);
m.isl_db     = 10 * log10(sum(side(:)) / sum(p(mainlobe)));
m.psl_tau    = tau_axis(sr);
m.psl_fd     = fd_axis(sc);
m.mainlobe   = mainlobe;
end
