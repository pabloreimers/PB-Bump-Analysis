# ionto_celltypes

Effect of iontophoresis on each imaged cell type, pooled across every dataset (built 2026-10-29):

- **Dopamine ejection** (all `.data/dopamine_ionto*.mat`): LPsP, EPG (7f / syt8m / syt8s), EPG>GRAB-DA2m, plus ACh controls.
- **LPsP stimulation** (all `.data/lpsp_p2x2*.mat` + the never-saved redo/syt8s trials): lpsp>P2X2 minus a stim-matched
  empty>P2X2 control, ATP ejection, EPG imaging.
- Per-fly **ATP vs dopamine** call with a confidence (dye SNR, EPG green response, file-name label).
- **Stim order**: do later ejections (same fly, same side) have a bigger dye-adjusted effect?

## Files

| | |
|---|---|
| `fig_ionto_celltypes.m` | Main script: loads everything, writes the 7 PDFs and `tables/`. Header documents every rule. |
| `lpsp_p2x2_unsaved_build.m` | Builds `.data/lpsp_p2x2_redo_syt8s_20261029.mat` from `Z:\pablo\lpsp_p2x2_redo` and `_syt8s` (never saved before). |
| `fig01_timecourses.pdf` | Peri-stim green + red traces per condition. |
| `fig02_summary_per_fly.pdf` | Per-fly effect, 0-3 s / 3-15 s / ejection side / opposite side. |
| `fig03_ejection_qc.pdf` | Per-stim ejection z and success rate. |
| `fig04_artifact_checks.pdf` | Fast (0.5 s) response, kinetics vs dye, dye vs response. |
| `fig05_glomerulus_heatmaps.pdf` | Glomerulus heatmaps, ejection side on top. |
| `fig06_atp_vs_dopamine_calls.pdf` | Per-fly substance evidence and calls. |
| `fig07_stim_order.pdf` | Dye and dye-adjusted response vs ejection number (count restarts at each side switch). |
| `fig_figure3_p2x2_traces.m` | Makes fig08 from the two Figure 3 lpsp_p2x2 datasets (reredo, walking post-7/8). |
| `fig09_figure3_p2x2_lme.m` / `fig09_figure3_p2x2_lme.pdf` | Mixed models on the fig08 ejection-trials data (one row per stim): `resp ~ dye + sat * lpsp + (1 \| fly)`, `sat = 1 - exp(-x/tau)`, tau profiled; rounds pooled vs continuous, compared with linear / log2 shapes by AIC (`tables/fig09_model_comparison.csv`). Best: pooled, tau ~ 1 stim; lpsp - empty after stim 1 = +0.0093 dF/F, p = 0.0005. |
| `fig08_figure3_p2x2_lpsp_minus_empty_ejtrials.pdf` | As fig08, but only trials with >= 1 confirmed ejection (all their commanded stims), split into round 1 / round 2 = first / second such trial of a same-side run (dotted line between). `ej_trials = true;` before running. |
| `fig08_figure3_p2x2_lpsp_minus_empty_ejtrials_pooled.pdf` | Same trial selection, all such trials pooled by stim # (no round split). `ej_trials = true; split_rounds = false;` |
| `fig08_figure3_p2x2_lpsp_minus_empty.pdf` | lpsp>P2X2, empty>P2X2 and lpsp - empty, ejection-side GCaMP + dye: reredo 2 s pulses (left), walking 0.5 s pulses in 5-stim bins (middle), per-fly mean over the 0-10 s inter-pulse interval per bin (right). 600 s / 59-pulse walking trials excluded (no dye response). Rightmost panel: build-up slope from `resp ~ dye + stim*lpsp + (1 | fly)` on the per-fly 5-stim-bin values (stim # = bin middle stim within the trial); lpsp and empty stim slopes at matched dye, 95% CIs. |
| `tables/p2x2_fly_ledger.csv`, `tables/p2x2_trial_ledger.csv` | Every p2x2 fly / trial on Z:, .mat coverage, duplicates, used or why not, substance call. |
| `tables/p2x2_z_trial_status.csv` | Current Z: folder status + genotype per trial; **input** to the main script. |
| `tables/per_fly_ejection_features.csv` | Per-fly dye / green features, votes, call, confidence. |
| `tables/p2x2_trials_used.csv`, `tables/stim_order_effects.csv` | Trials used per stratum; mixed-model order effects. |
| `z_audit/` | How the Z: audit was made (scan -> match -> status -> ledger), steps 1-6. |

## Running

From the repo root in MATLAB (R2023a, batch is fine):

```matlab
addpath('ugly_figures/ionto_celltypes'); fig_ionto_celltypes          % ~4 min with the cache
rebuild_p2x2 = true; fig_ionto_celltypes                              % ~20 min, reloads all lpsp_p2x2 .mat files
```

The p2x2 per-stim cache lives at `.data/fig_dopamine_ionto_p2x2_cache.mat` (gitignored). Rerun `z_audit/` only if trials
are added or moved on Z:; step 5 rewrites `tables/p2x2_z_trial_status.csv`, step 6 the ledgers.

## Decisions baked in (see the script header for detail)

- Trials deduplicated by dataset root + trial folder name; fly = root + date + fly number.
- Excluded: `open loop`, `exclusions`, `pilot (dont use)` folders (current Z: location) and 2024 `_low` trials.
- Blinded flies with no tag on the fly folder = lpsp>P2X2. th>P2X2 (2024) = LPsP stimulation.
- Walking p2x2 0.5 s flies are ATP (confirmed; low ejection strength), regardless of their dim dye.
- p2x2 flies called dopamine with high confidence get their own column, not the ATP comparison.
