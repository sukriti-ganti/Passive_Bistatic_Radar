% RUN_EXPERIMENT5_FEASIBILITY  Deliverable D6 and results/feasibility.csv.
%
% D6  maximum detection range for four target radar cross sections, at four
%     bandwidths and two integration times
%
% Assumptions, all of which belong in the report next to the table
%
%   eNodeB EIRP                     60 dBm
%   carrier                         1800 MHz, LTE band 3
%   receive antenna gain            10 dBi
%   receiver noise figure           3 dB
%   system losses                   3 dB
%   required SNR after processing   13 dB
%   propagation                     free space, no terrain, no buildings,
%                                   no atmospheric or rain loss
%
% Two results in this table need explaining rather than just reporting.
%
% First, the maximum range does not depend on the bandwidth.  Widening the
% channel raises the noise power by B and raises the coherent processing
% gain by B, and the two cancel exactly.  Bandwidth buys range resolution,
% which improves by a factor of sixteen and a half across the four channels,
% and it buys nothing at all in range.  The table shows both columns so the
% trade is visible.
%
% Second, the quoted range is the symmetric case, Rt = Rr = R.  The radar
% equation constrains the product Rt*Rr, so the true detection boundary is
% the Cassini oval Rt*Rr = R^2, not a circle.  Near the baseline that oval
% bulges outwards, which looks like good news until the range resolution
% figure is laid over it: the same region where the oval is widest is the
% region where cos(beta/2) collapses and the radar cannot resolve range at
% all.  Both figures are drawn here, on the same axes, for that reason.
%
% Set MAKE_FIGURES = false before running to compute without plotting.
%
% Runtime: seconds.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end

C = 299792458;

TARGETS   = {'pedestrian', 'small drone', 'car', 'truck'};
RCS       = [0.1 1 10 100];
NRB_LIST  = [6 25 50 100];
BW_LIST   = [1.08 4.50 9.00 18.00] * 1e6;
TINT_LIST = [0.05 0.2];

base = struct('eirp_dbm', 60, 'fc', 1800e6, 'gr_dbi', 10, ...
              'nf_db', 3, 'lsys_db', 3, 'snr_req_db', 13);

fprintf('\n===== Experiment 5: link budget and feasibility (D6) =====\n');

o = link_budget(setfield_local(base, {'bw_hz', 9e6, 't_int_s', 0.2, 'rcs_m2', 10}));
fprintf('reference case, 10 MHz channel, 0.2 s integration, 10 m^2 target\n');
fprintf('  wavelength            %.4f m\n', o.lambda_m);
fprintf('  noise power           %.2f dBm in %.2f MHz\n', o.noise_dbm, 9);
fprintf('  required echo power   %.2f dBm\n', o.pr_dbm);
fprintf('  pre-detection SNR     %.2f dB\n', o.snr_pre_db);
fprintf('  processing gain B*T   %.2f dB (%.2e)\n', o.gp_db, 9e6 * 0.2);
fprintf('  post-processing SNR   %.2f dB, the requirement\n', o.snr_pre_db + o.gp_db);
fprintf('  maximum range         %.2f km, and Rt*Rr <= %.3e m^2\n', ...
    o.max_range_m / 1e3, o.rt_rr_max);

% ------------------------------------------------------------ the table
rows = cell(numel(TARGETS) * numel(NRB_LIST) * numel(TINT_LIST), 5);
k = 0;
rmax = zeros(numel(RCS), numel(TINT_LIST));

fprintf('\n%-12s %8s %6s %10s %12s %12s\n', ...
    'target', 'RCS m^2', 'nrb', 'T (s)', 'Rmax (km)', 'dR (m)');
for it = 1:numel(TINT_LIST)
    for ir = 1:numel(RCS)
        for ib = 1:numel(NRB_LIST)
            q = base;
            q.bw_hz   = BW_LIST(ib);
            q.t_int_s = TINT_LIST(it);
            q.rcs_m2  = RCS(ir);
            out = link_budget(q);

            k = k + 1;
            rows(k, :) = {TARGETS{ir}, RCS(ir), NRB_LIST(ib), ...
                          TINT_LIST(it), out.max_range_m};
            rmax(ir, it) = out.max_range_m;

            fprintf('%-12s %8.1f %6d %10.2f %12.2f %12.2f\n', ...
                TARGETS{ir}, RCS(ir), NRB_LIST(ib), TINT_LIST(it), ...
                out.max_range_m / 1e3, out.range_res_m);
        end
    end
end

write_csv(fullfile(root, 'results', 'feasibility.csv'), ...
    {'target', 'rcs', 'nrb', 'integration_time', 'max_range_m'}, rows);

fprintf(['\nThe range column repeats across the four bandwidths because the\n' ...
         'noise power and the processing gain both scale with B and cancel.\n' ...
         'The resolution column is what the bandwidth actually buys.\n']);

% ---------------------------------------------------- geometry over the plane
tx = [0 0];
rx = [30e3 0];
bw = 9e6;
L  = 30e3;

x = linspace(-20e3, 50e3, 351);
y = linspace(-25e3, 25e3, 251);
[xx, yy] = meshgrid(x, y);
pts = [xx(:).'; yy(:).'];
g = geometry(tx, rx, pts, [0; 0], base.fc, bw);

dr_map = reshape(g.dR, size(xx));
rt_rr  = reshape(g.Rt .* g.Rr, size(xx));
rb_map = reshape(g.Rb, size(xx));

car_rmax = rmax(3, 2);                         % car, 0.2 s
frac_blind = mean(g.dR(:) > 10 * C / (2 * bw));
fprintf(['\n%.1f%% of the mapped plane has a range resolution more than ten\n' ...
         'times worse than c/2B, purely because of where it sits relative to\n' ...
         'the baseline.  That is the geometric blind region.\n'], 100 * frac_blind);

if MAKE_FIGURES
    h = figure('Name', 'maximum range against RCS');
    loglog(RCS, rmax / 1e3, '-o', 'LineWidth', 1);
    grid on;
    xlabel('bistatic radar cross section (m^2)');
    ylabel('maximum detection range, R_t = R_r (km)');
    title('Maximum detection range, 60 dBm EIRP at 1800 MHz, free space');
    legend(arrayfun(@(t) sprintf('T = %.2f s', t), TINT_LIST, ...
        'UniformOutput', false), 'Location', 'northwest');
    save_fig(h, 'd6_max_range_vs_rcs');
    close(h);

    h = figure('Name', 'geometry over the plane');

    subplot(1, 2, 1);
    z = 10 * log10(dr_map / (C / (2 * bw)));
    imagesc(x / 1e3, y / 1e3, z);
    set(gca, 'YDir', 'normal');
    caxis([0 20]);
    colorbar;
    hold on;
    plot(tx(1) / 1e3, tx(2) / 1e3, 'kv', 'MarkerFaceColor', 'k');
    plot(rx(1) / 1e3, rx(2) / 1e3, 'k^', 'MarkerFaceColor', 'k');
    plot([tx(1) rx(1)] / 1e3, [tx(2) rx(2)] / 1e3, 'k-');
    xlabel('x (km)'); ylabel('y (km)');
    title({'Range resolution loss, dB relative to c/2B', ...
           'the radar is blind along the baseline'});

    subplot(1, 2, 2);
    contour(x / 1e3, y / 1e3, sqrt(rt_rr) / 1e3, ...
        [5 10 car_rmax / 1e3 40 60 80], 'LineWidth', 1);
    hold on;
    contour(x / 1e3, y / 1e3, rb_map / 1e3, [2 5 10 20 40], 'k:');
    plot(tx(1) / 1e3, tx(2) / 1e3, 'kv', 'MarkerFaceColor', 'k');
    plot(rx(1) / 1e3, rx(2) / 1e3, 'k^', 'MarkerFaceColor', 'k');
    axis equal;
    axis([min(x) max(x) min(y) max(y)] / 1e3);
    grid on;
    colorbar;
    xlabel('x (km)'); ylabel('y (km)');
    title({sprintf('Coverage: contours of sqrt(R_t R_r) in km,'), ...
           sprintf('a car is detected inside %.0f km; dotted are iso-range ellipses', ...
           car_rmax / 1e3)});
    save_fig(h, 'd6_geometry_coverage');
    close(h);
end

fprintf('\nExperiment 5 complete.\n');
