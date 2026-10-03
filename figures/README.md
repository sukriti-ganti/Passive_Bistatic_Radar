# figures/

Committed as PDF, one file per deliverable, produced by the scripts in
`scripts/`. Nothing in here is written by hand.

| File | Deliverable | Produced by |
|---|---|---|
| `d2_ambiguity_surface_nrb{6,25,50,100}.pdf` | D2 | `run_experiment1_ambiguity.m` |
| `d2_zero_doppler_cut.pdf` | D2 | `run_experiment1_ambiguity.m` |
| `d2_zero_delay_cut.pdf` | D2 | `run_experiment1_ambiguity.m` |
| `d3_load_curve.pdf` | D3 | `run_experiment2_load_sweep.m` |
| `d4_rdmap_before_after.pdf` | D4 | `run_experiment3_cancellation.m` |
| `d4_cancellation_depth.pdf` | D4 | `run_experiment3_cancellation.m` |
| `d4_eca_vs_nlms.pdf` | D4 | `run_experiment3_cancellation.m` |
| `d4_eca_doppler_notch.pdf` | D4 | `run_experiment3_cancellation.m` |
| `d5_roc_in_band_snr.pdf` | D5 | `run_experiment4_detection.m` |
| `d5_roc_post_processing_snr.pdf` | D5 | `run_experiment4_detection.m` |
| `d6_max_range_vs_rcs.pdf` | D6 | `run_experiment5_feasibility.m` |
| `d6_geometry_coverage.pdf` | D6 | `run_experiment5_feasibility.m` |

To regenerate all of them, run `scripts/run_all.m` in MATLAB.

Before the report is submitted, every figure needs axis labels with units and
a caption. The scripts set the labels; the captions live in the Overleaf
document.
