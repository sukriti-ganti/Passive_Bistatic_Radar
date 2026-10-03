function out = geometry(tx_pos, rx_pos, tgt_pos, tgt_vel, fc, bw)
%GEOMETRY  Two dimensional bistatic radar geometry (module M3).
%
%   OUT = GEOMETRY(TX_POS, RX_POS, TGT_POS, TGT_VEL, FC)
%   OUT = GEOMETRY(TX_POS, RX_POS, TGT_POS, TGT_VEL, FC, BW)
%
%   Inputs
%     TX_POS   [x y] position of the illuminator (the eNodeB), metres
%     RX_POS   [x y] position of the passive receiver, metres
%     TGT_POS  [x y] target position, metres.  May be a 2 by M matrix of
%              column vectors, in which case every output is 1 by M
%     TGT_VEL  [vx vy] target velocity, metres per second.  Either one
%              column or one column per target
%     FC       carrier frequency, Hz
%     BW       signal bandwidth used for the range resolution, Hz.
%              Optional, defaults to 9e6 (the 10 MHz LTE channel)
%
%   Output struct OUT
%     L        baseline, transmitter to receiver, metres
%     Rt       transmitter to target range, metres
%     Rr       target to receiver range, metres
%     Rb       bistatic range Rt + Rr - L, metres.  This is the quantity a
%              passive radar actually measures, because the direct path
%              arrival defines the zero of the delay axis
%     tau      bistatic delay Rb/c, seconds
%     beta     bistatic angle at the target, radians
%     beta_deg the same in degrees
%     delta    angle between the velocity vector and the bistatic bisector
%     fd       bistatic Doppler, Hz, positive when the bistatic range closes
%     dR       bistatic range resolution c/(2*B*cos(beta/2)), metres
%     lambda   wavelength, metres
%     v        target speed, metres per second
%     u_tx     2 by M unit vectors pointing from the target to the transmitter
%     u_rx     2 by M unit vectors pointing from the target to the receiver
%     bisector 2 by M unit vectors along the bistatic bisector, the direction
%              in which a displacement changes bistatic range fastest.  A
%              ground displacement s along it changes Rb by 2*s*cos(beta/2),
%              which is why the ground range resolution is c/(2*B*cos(beta/2))
%     c        speed of light used, metres per second
%     fc, bw   echoed back for traceability
%
%   Doppler convention.  fd = -(1/lambda) * d(Rb)/dt, which expands to the
%   textbook form fd = (2*v/lambda)*cos(delta)*cos(beta/2).  It is computed
%   here from the projection of the velocity onto the two line of sight unit
%   vectors, which is the same number and stays well conditioned when the
%   bisector direction is undefined.
%
%   Geometric blindness.  As the target approaches the baseline, beta tends
%   to 180 degrees, cos(beta/2) tends to zero, and dR grows without bound.
%   Close to the baseline the radar cannot resolve range at all, no matter
%   how much bandwidth it has.  That region is plotted by
%   run_experiment5_feasibility.
%
%   Example
%     out = geometry([0 0], [10e3 0], [5e3 4e3], [50 0], 1800e6, 9e6);
%
%   See also SYNTH_RX, LINK_BUDGET.

if nargin < 5
    error('geometry:nargin', 'At least five input arguments are required.');
end
if nargin < 6 || isempty(bw)
    bw = 9e6;
end

C = 299792458;

tx = tx_pos(:);
rx = rx_pos(:);
if numel(tx) ~= 2 || numel(rx) ~= 2
    error('geometry:pos2d', 'Positions must be two element vectors, this is a 2-D model.');
end

tgt = tgt_pos;
if isvector(tgt), tgt = tgt(:); end
if size(tgt, 1) ~= 2
    error('geometry:tgt', 'tgt_pos must be 2 by M.');
end
m = size(tgt, 2);

vel = tgt_vel;
if isvector(vel), vel = vel(:); end
if size(vel, 1) ~= 2
    error('geometry:vel', 'tgt_vel must be 2 by M.');
end
if size(vel, 2) == 1 && m > 1
    vel = repmat(vel, 1, m);
end
if size(vel, 2) ~= m
    error('geometry:velcount', 'tgt_vel must have one column, or one per target.');
end

lambda = C / fc;

dt = tgt - repmat(tx, 1, m);        % transmitter to target
dr = repmat(rx, 1, m) - tgt;        % target to receiver

Rt = sqrt(sum(dt.^2, 1));
Rr = sqrt(sum(dr.^2, 1));
L  = sqrt(sum((rx - tx).^2));

Rb  = Rt + Rr - L;
tau = Rb / C;

% Unit vectors pointing from the target towards each site.
ut = -dt ./ repmat(max(Rt, realmin), 2, 1);
ur =  dr ./ repmat(max(Rr, realmin), 2, 1);

cosbeta = sum(ut .* ur, 1);
cosbeta = min(1, max(-1, cosbeta));            % guard against rounding
beta    = acos(cosbeta);

% d(Rb)/dt = -v . (ut + ur), and fd = -(1/lambda) d(Rb)/dt.
bisum = ut + ur;                               % magnitude 2*cos(beta/2)
fd    = sum(vel .* bisum, 1) / lambda;

v      = sqrt(sum(vel.^2, 1));
cosh2  = cos(beta / 2);
denom  = v .* max(cosh2, realmin) * 2;
cosdel = (fd * lambda) ./ max(denom, realmin);
cosdel(v == 0) = 0;
cosdel = min(1, max(-1, cosdel));
delta  = acos(cosdel);

dR = C ./ (2 * bw * max(cosh2, realmin));
dR(cosh2 <= 0) = Inf;

out = struct();
out.L        = L;
out.Rt       = Rt;
out.Rr       = Rr;
out.Rb       = Rb;
out.tau      = tau;
out.beta     = beta;
out.beta_deg = beta * 180 / pi;
out.delta    = delta;
out.fd       = fd;
out.dR       = dR;
out.lambda   = lambda;
out.v        = v;
out.u_tx     = ut;
out.u_rx     = ur;
bnorm        = sqrt(sum(bisum.^2, 1));
out.bisector = bisum ./ repmat(max(bnorm, realmin), 2, 1);
out.c        = C;
out.fc       = fc;
out.bw       = bw;
end
