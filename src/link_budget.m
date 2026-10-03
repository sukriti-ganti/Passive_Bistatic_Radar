function out = link_budget(params)
%LINK_BUDGET  Bistatic link budget and maximum detection range (module M8).
%
%   OUT = LINK_BUDGET(PARAMS)
%
%   Bistatic radar equation
%
%       Pr = Pt*Gt*Gr*lambda^2*sigma_b / ((4*pi)^3 * Rt^2 * Rr^2 * L_sys)
%
%   with coherent processing gain Gp = B*T applied afterwards.  At B = 9 MHz
%   and T = 0.2 s, Gp = 1.8e6, or 62.5 dB.  That gain is the only reason a
%   20 W illuminator is usable at all: the raw echo sits far below the noise
%   and only the correlation lifts it out.
%
%   PARAMS fields, all optional, scalars or equally sized arrays
%     eirp_dbm     eNodeB effective isotropic radiated power   60
%     fc           carrier frequency, Hz                       1800e6
%     gr_dbi       receive antenna gain                        10
%     nf_db        receiver noise figure                       3
%     lsys_db      system losses                               3
%     snr_req_db   SNR required after processing               13
%     bw_hz        signal bandwidth                            9e6
%     t_int_s      coherent integration time                   0.2
%     rcs_m2       bistatic radar cross section                10
%     temp_k       reference temperature                       290
%     load_factor  fraction of resource elements transmitted   1
%     beta_deg     bistatic angle used for the range resolution 0
%
%   OUT fields (same shape as the broadcast inputs)
%     max_range_m   maximum detection range for the symmetric case Rt = Rr.
%                   This is the usual way a bistatic budget is quoted: the
%                   real constraint is on the product Rt*Rr, and every
%                   target on the Cassini oval Rt*Rr = max_range_m^2 is
%                   equally detectable
%     rt_rr_max     that product, in square metres
%     pr_dbm        received echo power at max_range_m
%     noise_dbm     receiver noise power in the signal bandwidth
%     snr_pre_db    SNR before processing at max_range_m
%     gp_db         coherent processing gain, 10*log10(B*T)
%     range_res_m   bistatic range resolution c/(2*B*cos(beta/2))
%     lambda_m, and the parameters actually used
%
%   Why bandwidth does not change the range.  Widening B raises the noise
%   power by B and raises the processing gain by B, and the two cancel
%   exactly: SNR after processing is Pr*T/(k*T0*F), independent of B.
%   Bandwidth buys range resolution, not detection range.  The feasibility
%   table shows the same max_range_m in all four bandwidth columns and a
%   range resolution that improves by 16 times across them; that is a
%   result, not a bug, and it should be stated as such in the report.
%
%   Free space propagation only.  No terrain, no buildings, no atmospheric
%   or rain loss, and no multipath.  Real ranges over a city will be well
%   short of these numbers, which makes every figure here an upper bound.
%
%   Example
%     out = link_budget(struct('rcs_m2', [0.1 1 10 100], 't_int_s', 0.2));
%     disp(out.max_range_m);
%
%   See also GEOMETRY, CFAR2D.

if nargin < 1, params = struct(); end
params = default_field(params, 'eirp_dbm',    60);
params = default_field(params, 'fc',          1800e6);
params = default_field(params, 'gr_dbi',      10);
params = default_field(params, 'nf_db',       3);
params = default_field(params, 'lsys_db',     3);
params = default_field(params, 'snr_req_db',  13);
params = default_field(params, 'bw_hz',       9e6);
params = default_field(params, 't_int_s',     0.2);
params = default_field(params, 'rcs_m2',      10);
params = default_field(params, 'temp_k',      290);
params = default_field(params, 'load_factor', 1);
params = default_field(params, 'beta_deg',    0);

C  = 299792458;
KB = 1.380649e-23;

lambda  = C ./ params.fc;
eirp_w  = 10.^((params.eirp_dbm - 30) / 10) .* params.load_factor;
gr      = 10.^(params.gr_dbi / 10);
f_lin   = 10.^(params.nf_db / 10);
lsys    = 10.^(params.lsys_db / 10);
snr_req = 10.^(params.snr_req_db / 10);

noise_w = KB * params.temp_k .* params.bw_hz .* f_lin;
gp      = params.bw_hz .* params.t_int_s;

% Required echo power: SNR after processing = (Pr/N) * Gp = Pr*T/(k*T0*F).
pr_req = snr_req .* KB .* params.temp_k .* f_lin ./ params.t_int_s;

% Rt^2 * Rr^2 = numerator / pr_req, and Rt = Rr = R gives R^4.
num  = eirp_w .* gr .* lambda.^2 .* params.rcs_m2 ./ ((4 * pi)^3 * lsys);
r4   = num ./ pr_req;
rmax = r4.^0.25;

cosb = cos(params.beta_deg * pi / 360);      % cos(beta/2)

out = struct();
out.max_range_m = rmax;
out.rt_rr_max   = sqrt(r4);
out.pr_dbm      = 10 * log10(pr_req) + 30;
out.noise_dbm   = 10 * log10(noise_w) + 30;
out.snr_pre_db  = 10 * log10(pr_req ./ noise_w);
out.gp_db       = 10 * log10(gp);
out.range_res_m = C ./ (2 * params.bw_hz .* max(cosb, realmin));
out.lambda_m    = lambda;
out.params      = params;
end

% ======================================================================

function s = default_field(s, name, value)
if ~isfield(s, name) || isempty(s.(name))
    s.(name) = value;
end
end
