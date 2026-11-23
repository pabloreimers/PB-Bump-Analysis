# V2 figures — running notes

Decisions and to-dos for the V2 manuscript figures. Scripts: `ugly_figures/scripts/Figure_N_V2.m`; outputs: `ugly_figures/all_figs/V2/`, copied to Dropbox `Pablo Reimers/LPsP/V2/`.

## To include

- **Last figure (long-term effects of LPsP activation): change in cue following after LPsP>P2X2 stims.** In the iontophoresis dataset (`lpsp_p2x2_reredo_20260518`), the bright bar sweeps in open loop between stims. After stims, LPsP>P2X2 bumps follow those sweeps more. Tracking gain goes from 0.13 to 0.30, p = 0.016 within flies, while Empty>P2X2 goes from 0.06 to 0.10, n.s. The change vs Empty is not yet significant: p = 0.21.
  - The increase doesn't depend on the stimulated side, or on whether the bar sweeps toward or away from the push.
  - Scripts: `fig_p2x2_sweep_tracking.m`, `fig_p2x2_side_effects.m` (bottom row).
  - To tighten before using it: the "pre" baseline pools the start-of-trial sweep of every bright-bar trial, but only a fly's first trial is truly ATP-naive. Try using only that one.

## Come back to

- **Figure 1, Empty>P2X2 flies with large "displacement" (e.g. 20260505 fly 1, about +3 rad per stim).** This isn't a visible bump rotating. It's the PVA phase estimate traveling:
  - During the 2 s ATP pulse it's pulled into the pipette-site bright spot.
  - As the spot fades (2–5 s), it sweeps diagonally across the PB to static bright glomeruli (the medial edge, bins 16/17).
  - In weak-bump stims it wanders slowly through static stripes.
  - rho stays around 0.5, so rho masking doesn't catch it.
  - Example: `exports/V2/scratch/fig_p2x2_empty_example_stims.png`; that plot was made inline, not yet scripted.
  - Resolved in part: Figure 1 now starts the window at 4 s, which removes the bright-spot sweep. What's left is the slow wandering through static stripes in weak-bump stims.
  - It adds noise to both genotypes but doesn't create the LPsP effect. Starting the window later keeps it: 3–15 s p = 0.007, 4–15 s p = 0.008, 5–15 s p = 0.008; LPsP +1.0 to +1.4 rad, Empty about 0.
  - Options when we return:
    - a bump-quality filter, e.g. a single-peak von Mises fit, or bump contrast over the static pattern, instead of rho
    - leaving it as noise
  - Don't use a circular start-to-end measure. It can't count more than half a turn and wipes out both groups.

## Conventions (all V2 figures)

- Time is shown as a scale bar, never a time axis.
- Colorbars label only their min and max.
- Panels to be drawn later get a dashed placeholder box; final placement is done in Illustrator.

## Decided

- **Figure 1, P2X2 net displacement:** scored 4–15 s after each confirmed stim (changed from 2 s on 2026-11-23, so the window starts after the ATP-evoked bright spot has faded), while the bar is still (it holds from −0.5 to +18.6 s around each stim). Bump position is masked where rho < 0.1 and interpolated, and signs are flipped onto the left-stim convention. All 4 stims per trial are included. Stims 3–4 alone give a larger effect, which fits the build-up across stims.
- **Figure 1 example mean images:** made from the raw movies in Dropbox `Pablo Reimers/LPsP/<trial>/`. The PNG-crop fallback must not ship.
