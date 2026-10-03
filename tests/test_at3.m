function ok = test_at3()
%TEST_AT3  Acceptance test for M3, the bistatic geometry.
%
%   OK = TEST_AT3() returns true when every check passes.
%
%   Checks, from the PRD
%     a target on the baseline midpoint has zero bistatic range
%     a very distant target on the perpendicular bisector has beta near zero
%       and range resolution tending to c/2B
%     points computed to have equal bistatic range lie on an ellipse with
%       the transmitter and receiver at its foci
%
%   Plus the Doppler formula, checked against its textbook closed form, and
%   the blindness of the baseline region, which is the result the report has
%   to show as a figure.

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT3  bistatic geometry ===\n');
ok = true;

C   = 299792458;
tx  = [0 0];
rx  = [30e3 0];
fc  = 1800e6;
bw  = 9e6;
L   = 30e3;

% ------------------------------------------------- baseline midpoint
g = geometry(tx, rx, [L / 2, 0], [0 0], fc, bw);
ok = check('target on the baseline midpoint has Rb = 0', ...
    abs(g.Rb) < 1e-6, 'Rb = %.2e m', g.Rb) && ok;
ok = check('and a bistatic angle of 180 degrees', ...
    abs(g.beta_deg - 180) < 1e-6, 'beta = %.6f deg', g.beta_deg) && ok;
ok = check('so its range resolution is unusable', ...
    g.dR > 1e6, 'dR = %.2e m against c/2B = %.1f m', g.dR, C / (2 * bw)) && ok;

% Every point on the baseline between the sites has Rb = 0.
xs = linspace(0.05, 0.95, 9) * L;
gb = geometry(tx, rx, [xs; zeros(1, numel(xs))], [0 0], fc, bw);
ok = check('the whole baseline is one zero range cell', ...
    max(abs(gb.Rb)) < 1e-6, 'max |Rb| = %.2e m over 9 points', max(abs(gb.Rb))) && ok;

% ------------------------------------------------- distant target
g = geometry(tx, rx, [L / 2, 1e9], [0 0], fc, bw);
ok = check('a very distant target has beta tending to zero', ...
    g.beta_deg < 1e-2, 'beta = %.2e deg', g.beta_deg) && ok;
ok = check('and dR tending to c/2B', ...
    abs(g.dR - C / (2 * bw)) / (C / (2 * bw)) < 1e-6, ...
    'dR = %.6f m against c/2B = %.6f m', g.dR, C / (2 * bw)) && ok;

% ------------------------------------------------- iso-range ellipse
% An ellipse with foci at tx and rx has Rt + Rr constant by definition, so
% Rb = Rt + Rr - L is constant too.  Build one and confirm.
rb_target = 8e3;
a  = (L + rb_target) / 2;                     % semi major axis
b  = sqrt(a^2 - (L / 2)^2);                   % semi minor axis
th = linspace(0, 2 * pi, 37);
pts = [a * cos(th) + L / 2; b * sin(th)];
ge  = geometry(tx, rx, pts, [0 0], fc, bw);
spread = max(ge.Rb) - min(ge.Rb);
ok = check('points on an ellipse with tx and rx at the foci share one Rb', ...
    spread < 1e-6 && abs(mean(ge.Rb) - rb_target) < 1e-6, ...
    'Rb = %.4f m, spread %.2e m over 37 points', mean(ge.Rb), spread) && ok;

% And the converse: a point off the ellipse has a different Rb.
goff = geometry(tx, rx, [L / 2; b + 500], [0 0], fc, bw);
ok = check('a point off that ellipse has a different Rb', ...
    goff.Rb > rb_target + 100, 'Rb = %.1f m against %.1f m', goff.Rb, rb_target) && ok;

% ------------------------------------------------- Doppler
% The implementation projects the velocity onto the two lines of sight.
% The PRD quotes fd = (2v/lambda)*cos(delta)*cos(beta/2).  They must agree.
set_seed(11);
pts = [rand(1, 200) * 60e3 - 15e3; rand(1, 200) * 40e3 + 1e3];
vel = [rand(1, 200) * 200 - 100; rand(1, 200) * 200 - 100];
gd  = geometry(tx, rx, pts, vel, fc, bw);
closed = (2 * gd.v / gd.lambda) .* cos(gd.delta) .* cos(gd.beta / 2);
ok = check('Doppler matches (2v/lambda)*cos(delta)*cos(beta/2)', ...
    max(abs(gd.fd - closed)) < 1e-6 * max(abs(gd.fd)), ...
    'worst difference %.2e Hz over 200 random targets', max(abs(gd.fd - closed))) && ok;

% A target closing along the bisector has the largest possible Doppler.
g = geometry(tx, rx, [L / 2, 20e3], [0 -100], fc, bw);
fd_max = 2 * 100 / g.lambda * cos(g.beta / 2);
ok = check('closing along the bisector gives the maximum Doppler', ...
    abs(g.fd - fd_max) < 1e-6 && g.fd > 0, ...
    'fd = %.1f Hz, delta = %.1f deg', g.fd, g.delta * 180 / pi) && ok;

g = geometry(tx, rx, [L / 2, 20e3], [100 0], fc, bw);
ok = check('crossing the bisector gives zero Doppler', ...
    abs(g.fd) < 1e-9, 'fd = %.2e Hz, delta = %.1f deg', g.fd, g.delta * 180 / pi) && ok;

% ------------------------------------------------- blind region
% Range resolution degrades as the target approaches the baseline.  Find how
% close a target can get before dR exceeds 10 times its far field value.
ys = logspace(1, 5, 400);
gg = geometry(tx, rx, [repmat(L / 2, 1, numel(ys)); ys], [0 0], fc, bw);
i  = find(gg.dR < 10 * C / (2 * bw), 1);
ok = check('the baseline blind region is bounded and measurable', ...
    ~isempty(i) && all(diff(gg.dR) < 0), ...
    'dR falls below 10x c/2B only beyond %.0f m from the baseline', ys(i)) && ok;

fprintf('AT3 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
