% RUN_ALL  Regenerate every deliverable from scratch.
%
% Runs the five experiments in order, which between them produce
%
%   data/wave_nrb{N}_rho{RRR}.mat   D1  the waveform library
%   figures/d2_*.pdf                D2  ambiguity surfaces and cuts
%   figures/d3_load_curve.pdf       D3  sidelobe level against cell load
%   results/af_metrics.csv
%   figures/d4_*.pdf                D4  cancellation, before and after
%   results/cancellation.csv
%   figures/d5_*.pdf                D5  ROC curves
%   results/roc.csv
%   figures/d6_*.pdf                D6  feasibility and coverage
%   results/feasibility.csv
%
% Set MAKE_FIGURES = false to compute the numbers without plotting, which is
% what a headless run needs.  Set QUICK = true to shorten experiment 4.
%
% Runtime: ten to twenty minutes on the free MATLAB Online tier, most of it
% the Monte Carlo in experiment 4.  Everything is written to disk as it goes,
% so an interrupted run can be resumed by running the remaining scripts
% individually.

here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(fullfile(root, 'src'));

if ~exist('MAKE_FIGURES', 'var'), MAKE_FIGURES = true; end
if ~exist('QUICK', 'var'), QUICK = false; end

t_all = tic;
run(fullfile(here, 'run_experiment1_ambiguity.m'));
run(fullfile(here, 'run_experiment2_load_sweep.m'));
run(fullfile(here, 'run_experiment3_cancellation.m'));
run(fullfile(here, 'run_experiment4_detection.m'));
run(fullfile(here, 'run_experiment5_feasibility.m'));

fprintf('\nAll experiments complete in %.1f minutes.\n', toc(t_all) / 60);
fprintf('Now run tests/run_all_tests.m to confirm the acceptance tests still pass.\n');
