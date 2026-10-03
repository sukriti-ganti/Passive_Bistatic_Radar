function ok = test_at8()
%TEST_AT8  Acceptance test for M8, the link budget and feasibility table.
%
%   OK = TEST_AT8() returns true when every check passes.
%
%   The PRD asks for a table of maximum detection range for each target
%   class, for each of the four bandwidths, at two integration times, with
%   every cell traceable to the bistatic radar equation under the stated
%   assumptions.  This test builds that table and then re-derives every cell
%   independently, in decibels, from the equation written out longhand.  If
%   the two disagree the table is not traceable and the test fails.
%
%   Stated assumptions
%     eNodeB EIRP                60 dBm
%     carrier                    1800 MHz, LTE band 3
%     receive antenna gain       10 dBi
%     receiver noise figure      3 dB
%     system losses              3 dB
%     required SNR after processing  13 dB
%     free space propagation, no terrain and no atmospheric loss

here = fileparts(mfilename('fullpath'));
addpath(here);
addpath(fullfile(fileparts(here), 'src'));

fprintf('\n=== AT8  link budget and feasibility ===\n');
ok = true;

C  = 299792458;
KB = 1.380649e-23;
T0 = 290;

targets   = {'pedestrian', 'small drone', 'car', 'truck'};
rcs       = [0.1 1 10 100];
nrb_list  = [6 25 50 100];
bw_list   = [1.08 4.50 9.00 18.00] * 1e6;
tint_list = [0.05 0.2];

p = struct('eirp_dbm', 60, 'fc', 1800e6, 'gr_dbi', 10, ...
           'nf_db', 3, 'lsys_db', 3, 'snr_req_db', 13);

fprintf('-- maximum detection range, metres\n');
fprintf('   %-12s %8s %10s %10s\n', 'target', 'T (s)', 'B (MHz)', 'Rmax (m)');

worst = 0;
for it = 1:numel(tint_list)
    for ib = 1:numel(bw_list)
        for ir = 1:numel(rcs)
            q = p;
            q.bw_hz  = bw_list(ib);
            q.t_int_s = tint_list(it);
            q.rcs_m2 = rcs(ir);
            out = link_budget(q);

            % Independent re-derivation, entirely in dB.
            lambda   = C / q.fc;
            eirp_dbw = q.eirp_dbm - 30;
            num_db   = eirp_dbw + q.gr_dbi + 20 * log10(lambda) + 10 * log10(q.rcs_m2) ...
                     - 30 * log10(4 * pi) - q.lsys_db;
            pr_req_db = q.snr_req_db + 10 * log10(KB * T0) + q.nf_db ...
                      - 10 * log10(q.t_int_s);
            r4_db  = num_db - pr_req_db;
            rmax_h = 10^(r4_db / 40);

            err = abs(rmax_h - out.max_range_m) / out.max_range_m;
            worst = max(worst, err);

            if ib == 3
                fprintf('   %-12s %8.2f %10.2f %10.0f\n', ...
                    targets{ir}, q.t_int_s, q.bw_hz / 1e6, out.max_range_m);
            end
        end
    end
end

ok = check('every table cell reproduces a longhand dB calculation', worst < 1e-9, ...
    'worst relative difference %.2e over %d cells', ...
    worst, numel(tint_list) * numel(bw_list) * numel(rcs)) && ok;

% ------------------------------------------------- scaling laws
fprintf('-- scaling laws\n');
a = link_budget(setfields(p, 'rcs_m2', 1,   'bw_hz', 9e6, 't_int_s', 0.2));
b = link_budget(setfields(p, 'rcs_m2', 100, 'bw_hz', 9e6, 't_int_s', 0.2));
ok = check('range scales as the fourth root of RCS', ...
    abs(b.max_range_m / a.max_range_m - 100^0.25) < 1e-9, ...
    '100x RCS gives %.4f times the range, expected %.4f', ...
    b.max_range_m / a.max_range_m, 100^0.25) && ok;

c = link_budget(setfields(p, 'rcs_m2', 1, 'bw_hz', 9e6, 't_int_s', 0.05));
ok = check('range scales as the fourth root of integration time', ...
    abs(a.max_range_m / c.max_range_m - 4^0.25) < 1e-9, ...
    '4x integration gives %.4f times the range, expected %.4f', ...
    a.max_range_m / c.max_range_m, 4^0.25) && ok;

% Bandwidth: noise power and processing gain both scale with B, so they
% cancel and the range does not move.  Resolution is what improves.
r = zeros(1, 4); dres = zeros(1, 4);
for ib = 1:4
    o = link_budget(setfields(p, 'rcs_m2', 10, 'bw_hz', bw_list(ib), 't_int_s', 0.2));
    r(ib) = o.max_range_m;
    dres(ib) = o.range_res_m;
end
ok = check('bandwidth does not change the detection range', ...
    (max(r) - min(r)) / mean(r) < 1e-12, ...
    'range %.0f m at every bandwidth from %.2f to %.2f MHz', ...
    mean(r), bw_list(1) / 1e6, bw_list(end) / 1e6) && ok;
ok = check('but it does change the range resolution, proportionally', ...
    abs(dres(1) / dres(4) - bw_list(4) / bw_list(1)) < 1e-9, ...
    'resolution improves from %.1f m to %.1f m, a factor of %.1f', ...
    dres(1), dres(4), dres(1) / dres(4)) && ok;

% ------------------------------------------------- processing gain
o = link_budget(setfields(p, 'bw_hz', 9e6, 't_int_s', 0.2, 'rcs_m2', 10));
ok = check('processing gain B*T at 9 MHz and 0.2 s is 62.5 dB', ...
    abs(o.gp_db - 62.55) < 0.1, '%.2f dB, a gain of %.2e', o.gp_db, 9e6 * 0.2) && ok;

% ------------------------------------------------- sanity of the numbers
ok = check('noise power in 9 MHz with a 3 dB noise figure is about -101 dBm', ...
    abs(o.noise_dbm + 101.4) < 0.5, '%.2f dBm', o.noise_dbm) && ok;
ok = check('the raw echo is far below the noise, as expected', ...
    o.snr_pre_db < -40, 'pre-detection SNR %.1f dB, recovered by %.1f dB of gain', ...
    o.snr_pre_db, o.gp_db) && ok;
ok = check('a car is detectable at tens of kilometres, a pedestrian at under ten', ...
    o.max_range_m > 25e3 && o.max_range_m < 40e3, ...
    'car at %.1f km, pedestrian at %.1f km', o.max_range_m / 1e3, ...
    link_budget(setfields(p, 'bw_hz', 9e6, 't_int_s', 0.2, 'rcs_m2', 0.1)).max_range_m / 1e3) && ok;

% ------------------------------------------------- geometry interaction
% The equation constrains the product Rt*Rr, so the equal range figure is
% one point on a Cassini oval, not a circle.
ok = check('rt_rr_max is the square of the equal range figure', ...
    abs(o.rt_rr_max - o.max_range_m^2) / o.rt_rr_max < 1e-12, ...
    'Rt*Rr <= %.3e m^2', o.rt_rr_max) && ok;

% Cell load reduces the radiated power and so the range.
half = link_budget(setfields(p, 'bw_hz', 9e6, 't_int_s', 0.2, 'rcs_m2', 10, 'load_factor', 0.5));
ok = check('halving the cell load costs 16% of the range', ...
    abs(half.max_range_m / o.max_range_m - 0.5^0.25) < 1e-9, ...
    '%.1f km against %.1f km, a factor of %.4f', ...
    half.max_range_m / 1e3, o.max_range_m / 1e3, half.max_range_m / o.max_range_m) && ok;

fprintf('AT8 result: %s\n', pass_string(ok));
end

% ======================================================================

function s = setfields(s, varargin)
for i = 1:2:numel(varargin)
    s.(varargin{i}) = varargin{i + 1};
end
end

function s = pass_string(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end
