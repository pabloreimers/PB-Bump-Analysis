# PB-Bump-Analysis

MATLAB analysis code for imaging experiments on the fly protocerebral bridge (PB) — EPG
"bump" dynamics, GRAB-DA dopamine sensor imaging, dopamine iontophoresis, LPsP silencing,
gain-change/closed-loop behavior, and related in-vivo two-photon imaging paradigms in
*Drosophila*.

## Layout

- **Top-level `*.m` / `*.asv` scripts** — legacy, mostly interactive analysis scripts
  written as MATLAB cell-mode sections (`%%`), meant to be run one section at a time,
  not executed top-to-bottom as a function. Filenames describe the experiment
  (`dopamine_ionto_rnai.m`, `epg_grab_script.m`, `native_gain_analysis.m`, etc.). `.asv`
  files are MATLAB autosave backups of the matching `.m` file — ignore unless the `.m`
  file looks broken.
- **`claude/`** — newer, modular helper functions with documented headers, `pb_`-prefixed
  (`pb_load_raw_tif.m`, `pb_register.m`, `pb_remove_scan_noise.m`, `pb_process_trial.m`,
  `pb_make_mask.m`, ...). This is the current preprocessing pipeline: raw ScanImage tif →
  scan-noise removal → optional shot-noise reduction → motion correction. Prefer building
  on these over the top-level legacy scripts. `pb_process_trial.m` is the entry point and
  documents the scan-noise vs. shot-noise distinction in its header.
- **`data scripts/`** — one analysis script per experiment/dataset, similar cell-mode
  style to the top-level scripts.
- **`data/`** — processed `.mat` datasets (gitignored via `*.mat`; these are multi-GB and
  never committed). Scripts `load(uigetfile(...))` these interactively.
- **`ugly_figures/`** — figure-generation workspace (new). See below.
- **External/vendored toolboxes** — `circ_stats/` (circular statistics), `28790/`
  (FileExchange `colorspace`), `carl files/` (tif I/O helpers from a collaborator). Treat
  as third-party, don't modify.
- **`example_movies/`** — sample `.avi`s used by some scripts for illustration.

## Conventions

- Scripts are written for interactive use in the MATLAB Editor (`%% section`, run with
  Ctrl+Enter), not as batch/CLI scripts. When adding new analysis, follow this style
  unless writing a reusable function.
- New reusable pipeline code goes in `claude/` with a `pb_` prefix and a header comment
  block explaining what it does and why (see `pb_process_trial.m` / `pb_remove_scan_noise.m`
  for the expected level of detail on non-obvious tradeoffs).
- Data files (`*.mat`, `*.avi` in `example_movies` excepted) are gitignored — don't fight
  this by force-adding data.

## Figures workspace (`ugly_figures/`)

New figure-making work for this analysis lives in `ugly_figures/`:

- `ugly_figures/scripts/` — one `%%`-sectioned MATLAB script per figure or figure set,
  prefixed `fig_` (e.g. `fig_bump_amplitude_vs_velocity.m`).
- `ugly_figures/exports/` — rendered output (`.png`/`.pdf`/`.eps`); gitignored by default
  since these are regenerable from the scripts.
