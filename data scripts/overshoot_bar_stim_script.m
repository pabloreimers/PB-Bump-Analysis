%% overshoot_bar_stim_script
% Single-trial analysis for the "sweeping bar, then held bar + LED stim"
% protocol in Z:\noah_np123\Data\flyg\overshoot\<fly>\trial_NNN\ (first
% instance: 20260930-1_epg_syt8s_lpsp_cschrimson\trial_002). Same raw-data
% format as overshoot_drug_stim_script.m / led_stim_berg4_script.m (raw
% ScanImage fastZ tif + one WaveSurfer h5), so the tif -> z-summed video ->
% pb_register -> mask front end is carried over from those, and the
% glomerulus clustering / PVA bump estimate follows this repo's process_im
% convention (lpsp_cschrimson_reredo_script.m: 16 glomeruli per hemisphere,
% alpha repeated over the two halves, PVA on the z-scored traces).
%
% ---- trial structure (checked against the protocol csv + h5) ----------------
% trial_002, dangling_burn_in_sweep36_hold270_913s: bar sweeps back and forth
%   for VR t = 0-612.5 s, then is HELD at 270 deg for 612.5-912.5 s. VR t=0 is
%   DAQ 39.4 s; imaging = 9704 volumes, DAQ 30.0-953.5 s (13 frames/volume =
%   10 valid + 3 flyback, 10.50 Hz); LED = 30 pulses of exactly 5 s on / 5 s
%   off, DAQ 640.0-935.0 s -- only during the hold (apart from the first
%   pulse or two, see note 4). Matches the "~10 min sweep, then 5 min hold
%   with 5 s on / 5 s off" description.
% trial_003, dangling_burn_in_sweep36_190s: the same sweep for 190 s only.
%   2111 volumes (DAQ 30.0-230.8 s), NO held bar, NO LED pulses -> sections
%   5-6 run (heatmap, bump vs bar), section 7 (stim diff image) is skipped.
% trial_004, dangling_sweep_15_36_two_dirs_holds_613s: 30 s sweeps of 450 deg
%   alternating direction, each followed by a ~10 s hold (at 0/90/180/270
%   deg). 6551 volumes (DAQ 30.0-653.4 s); LED = 50 pulses of 2 s every 10 s
%   from DAQ 60-552 s, running through sweeps and holds alike. Section 7
%   therefore falls back to using ALL pulses (LED-off time between pulses as
%   the control) when fewer than 3 pulses sit fully inside a hold.
% 20260930-2 / trial_001 (second fly): sweep -> ~7 min hold at 270 deg (DAQ
%   ~349-769 s) -> sweep again. 10955 volumes (142415 frames, i.e. well past
%   the 65535-frame libtiff limit, see note 1). LED = 30 pulses of 5 s on / 5 s
%   off, DAQ 400-695 s, ALL inside the hold. This trial folder has NO
%   scan_params.json -> obs_fallback_scan_params is used; the tif header text
%   (SI.hStackManager.numFramesPerVolume etc.) was checked and matches fly 1.
% 20260930-2 / trial_002, dangling_burn_in_sweep36_hold270_sweep36_1032s: same
%   protocol and LED timing as trial_001 of this fly (hold DAQ 349.9-768.2 s,
%   30 x 5 s pulses 400-695 s, all inside the hold), 10954 volumes, again no
%   scan_params.json. 187 volumes (1.7 %) excluded for |shift| > 10 px.
% 20261001-1_epg_syt8s_cyoTM2 / trial_001 (NO-OPSIN CONTROL, same 1032 s
%   protocol): hold DAQ 352.1-770.4 s, 30 x 5 s pulses 400-695 s inside it,
%   10956 volumes, zoom 5 (not 4) -> obs_scan_params now reads zoom/planes/
%   pixels from the tif header and rescales um_per_px (0.2531 um/px). 405
%   volumes (3.7 %) excluded. NB: the LED still modulates PB fluorescence here
%   (rho peaks at every LED onset, bilaterally symmetric diff image) despite
%   no CsChrimson -> a light artifact / visual response confound to keep in
%   mind when interpreting the CsChrimson flies. (Checked with
%   overshoot_bar_stim_led_bleed_check.m: the LED adds a spatially FLAT
%   offset to the background, which the per-volume subtraction removes; the
%   PB response itself is real, i.e. the fly responds to seeing the LED.)
% 20261001-1_epg_syt8s_cyoTM2 / trial_002 (control, hold at 90 deg,
%   dangling_burn_in_sweep36_hold90_sweep36_1027s): has scan_params.json
%   (zoom 5). 10903 volumes, hold DAQ 345.8-764.1 s, 30 x 5 s pulses 400-695
%   s, 14 volumes excluded. Every LED pulse flips the bump ~165 deg and it
%   returns between pulses -- the same signature as the CsChrimson fly 1, in
%   a fly with no opsin.
% 20261001-1_epg_syt8s_cyoTM2 / trial_003 (control, hold at 270 deg again,
%   no scan_params.json): 10955 volumes, hold DAQ 350.8-769.1 s, pulses
%   400-695 s, 0 excluded. Tracking poor throughout (offset drifts 173 -> 90
%   deg over the pre-stim sweeps, circ std 45-70 deg) and the LED does almost
%   nothing: the bump sits at ~0 deg through the first 7 pulses, jumps once
%   to ~-130 deg at ~465 s and stays; diff image is ~10x weaker than trials
%   1-2. The lag estimate (1.33 s) is unreliable here.
% 20261002-1_epg_syt8s_lpsp_cschjrimson / trial_001 (CsChrimson fly 3, hold
%   at 270 deg, no scan_params.json, zoom 5): 10955 volumes, hold DAQ
%   345.0-763.3 s, pulses 394.4-689.4 s, 0 excluded. The medial PB (cols
%   ~100-150) is out of the imaging volume, so the two hemispheres come out
%   as separate blobs and the imclose bridge between them carries no signal:
%   glomeruli 13-18 are NaN'd (maskCore logic in sections 3/5). During the
%   LED block the bump ROTATES one full turn per 10-s on/off cycle (cf.
%   20260930-2 trial_001); post-stim tracking is poor (circ std 59-95 deg)
%   and the offset changes from ~-90 to ~+20 deg.
% 20261002-1 / trial_002 (same fly, hold at 270 deg): 10954 volumes, hold DAQ
%   349.1-767.5 s, pulses 400-695 s, 0 excluded; medial PB again missing ->
%   glomeruli 13-19 NaN'd. Tracking bimodal throughout (main diagonal at
%   offset ~0 plus an anti-phase branch; circ std 55-89 deg pre AND post),
%   bump parks at -160 deg before the LED, then rotates/drifts through the
%   LED block (less regular than trial_001); rho does not collapse.
% 20261005-1_epg_syt8s_lpsp_cschrimson / trial_001 (CsChrimson fly 4,
%   dangling_burn_in_sweep36_190s, zoom 4): 2111 volumes, no hold / no LED
%   -> section 7 skipped. Clean mask (both hemispheres contiguous). The bump
%   is coherent (L/R hemisphere PVAs agree to 19 deg) but does NOT track the
%   bar: it rotates ~3x faster than the sweep, loosely in the same direction
%   (velocity r = 0.51), circ std 126 deg even after the lag fit (the 2.84 s
%   "lag" is meaningless here). The fly was essentially stationary on the
%   ball (iter.ft_heading flat, forward velocity ~0), so this is not
%   self-motion-driven either: the bump rotates on its own at 2-3 rad/s in
%   ~20-30 s bursts, loosely in the sweep direction. Burn-in trial.
% 20261005-1 / trial_002: ABORTED after 30 s (319 volumes, 222 MB tif) with
%   NO fluorescence (mean image ~0 raw units everywhere vs ~7000 at the PB in
%   trials 001/003 -- shutter/laser/PMT off). Registration blows up (+-45 px)
%   and the mask becomes the whole frame; outputs were deleted. Not analysable.
% 20261005-1 / trial_003, dangling_sweep_12_24_36_prehold270_347s: 12
%   constant-speed sweeps at 12/24/36 deg/s (alternating direction), each
%   followed by a 5-s hold, then a hold at 270 deg. 3757 volumes (DAQ 30-387
%   s). Run with obs_min_hold_sec = 2.5 (5-s holds erode to ~3.5 s in the
%   detector) -> 12 holds found. The one LED pulse (DAQ 400-405 s) is AFTER
%   imaging ended -> diff image skipped (guard added in section 7). Bump
%   again free-running: ~110-135 deg/s regardless of bar speed (ratio 8.8 /
%   4.9 / 3.6 at 12 / 24 / 36 deg/s), ~75 % in the bar's direction, and it
%   keeps rotating at ~110 deg/s through the holds. Fly stationary on the
%   ball. Protocol -> DAQ offset 39.6 s (scratchpad probe_speeds.m).
%   (NB those per-volume velocity ratios are inflated by PVA jitter; the
%   fit-based numbers for trial_004 below are the reliable kind.)
% 20261005-1 / trial_004, dangling_sweep_6_two_dirs_holds_300s: four 6 deg/s
%   sweeps alternating direction (65/60 s each), each followed by a 10-s
%   hold; no LED. 3258 volumes, 0 excluded, 4 holds detected. Cleanest
%   demonstration of this fly's behaviour: during sweeps the bump rotates at
%   39-60 deg/s = 6.5-10x the bar speed, same sign in all four sweeps (6-10
%   turns per sweep), and during the holds it (almost) stops (|vel| 3-30
%   deg/s, median 13). Velocity-coupled to the bar, not position-coupled:
%   hold positions do not map onto bar angle. Fit on a 0.5-s-smoothed,
%   unwrapped PVA (scratchpad probe_sweeps_slow.m). Fly still mostly
%   stationary (0.6 turns of heading over the trial).
% 20261005-1 / trial_005, dangling_sweep_12_24_36_hold270_sweep_12_24_36_1097s:
%   speed series (0-342 s) -> hold at 270 deg (DAQ 382-801 s) with 30 x 5 s
%   LED pulses 400-695 s -> same speed series again. 11636 volumes, 51
%   excluded. (The 5-s holds erode below min_hold_sec = 4, so only the long
%   hold is detected -- intended, it gives the pre/post split.) PRE: bump
%   free-running as in trials 001/003/004 (circ std 88-126 deg). LED block:
%   bump parks at ~-160 deg for pulses 1-5, then slowly rotates, then rho
%   collapses to ~0.15 from ~560 s to the end of the hold. POST: rho jumps
%   to 0.6-0.9 and the bump TRACKS the bar (1-min circ std 68/64/74/28/24/98
%   deg, offset ~0). First time this fly tracked the cue.
% 20261005-1 / trial_006 (same protocol, ScanImage epoch 12:28 vs 12:08 for
%   trial_005; both copied to the server at 13:15): 11636 volumes, 2
%   excluded. Raw PB fluorescence is higher than in trial_005 (1586 vs 1073
%   pre-hold) but rho is much lower throughout. Both trials show the raw PB
%   mean stepping with the LED (005: +9 %, 006: -5 %) while the background
%   barely moves -- real, but spatially uniform, so no bump information.
%   PRE: bump now tracks loosely from the start (circ std
%   74-99 deg, offset ~0, slope 1 -- no lapping), rho 0.3-0.9 declining.
%   LED block: rho ~0.15-0.3 with a small rise on each LED onset; bump
%   wanders; diff image is a uniform DIMMING of the whole PB with LED on
%   (-200..-400 raw units; the opposite sign to trial_005's uniform
%   brightening). POST: rho stays low (0.1-0.3), tracking loose (65-122 deg).
%   So the trial_005 post-stim tracking did NOT persist into this trial.
%   ** TTX was in the bath for trial_006 ** (user). The LED-locked dimming
%   is therefore the direct (graded) LPsP->EPG effect. Quantified in
%   overshoot_ttx_inhibition_map.m: -8.5 % dF/F on average, non-uniform
%   (-3.5 .. -13.9 %, ANOVA p = 3e-6, L/R r = 0.56; deepest in glomeruli
%   7-11 and 26-29, weakest at the lateral tips), rebound after offset, and
%   NOT correlated with where the trial_005 bump sat during stimulation
%   (r ~ 0.1, perm p > 0.7). The response is BIPHASIC
%   (overshoot_ttx_early_vs_late.m): early 0-2 s -0.7 % vs late 2-5 s
%   -7.7 % pooled (paired p = 5e-4; 17/32 glomeruli early > late, Bonf.);
%   glomeruli 2-4 and 19-23 show a transient +6-7 % EXCITATION peaking ~1 s
%   after onset (early > 0 in 22/30 pulses), while 7-12 / 26-31 are
%   inhibited from onset. In saline (trial_005) the same glomeruli show no
%   such split (+8.5 % early, +5.7 % late, uniform). Neither the early nor
%   the late TTX map correlates with trial_005 bump occupancy during the
%   pulses or in the LED-off gaps (|r| <= 0.25, perm p > 0.6;
%   overshoot_ttx_occupancy_vs_phase.m) -- the bump occupancy there was
%   near-uniform (0.55-2.2x) because the bump rotated / collapsed. The only
%   structured relation is with the pre-LED parked position (glomeruli
%   13-18 = "neutral" group, KW p = 0.01), i.e. the bump parked between the
%   excited and inhibited zones -- n = 1 hold, so anecdotal.
%   Split by halves of the LED block (evp_block_part): in the FIRST half the
%   trial_005 bump sat preferentially in the TTX early-EXCITED glomeruli
%   (occupancy 1.29 / 1.18 / 0.51 x uniform for excited / neutral /
%   inhibited, KW p < 0.01); in the SECOND half it had moved to the
%   early-INHIBITED glomeruli (0.27 / 0.89 / 1.73, KW p < 0.01; r(early,
%   occupancy) = -0.64, r(early, mean z) = -0.81, perm p 0.06-0.19) and
%   stayed there after the LED (0.08 / 1.13 / 1.55). One slow relocation,
%   so the per-glomerulus correlations are not independent evidence.
% 20261005-1 / trial_007 (same protocol, ScanImage epoch 13:15, ~47 min
%   after the TTX trial; bath condition to be confirmed): 11636 volumes, 0
%   excluded, has scan_params.json. Bump amplitude low throughout (rho
%   0.18-0.28), no bar tracking at any point: the PVA sits in a fixed band
%   (-60..-90 deg pre and during the hold, then jumps to ~+110 deg at ~960 s
%   and stays there with rho rising to 0.38 while the bar sweeps). The LED
%   does essentially NOTHING here: per-glomerulus dF/F early +0.2 %, late
%   +1.5 %, no glomerulus significant, raw PB step +1.9 % (vs -7.4 % under
%   TTX in trial_006 and +4.1 % in saline trial_005); the trial_006 late-
%   inhibition profile is not reproduced (r = -0.66, i.e. if anything
%   inverted). Consistent with spiking still blocked and the direct LED
%   effect having run down (desensitisation / depletion / prep decline).
% 20261005-1_epg_syt8_tmp2 (NO-OPSIN CONTROL #2, EPG>syt-GCaMP8 / TM2),
%   trials 001 and 002 (epochs 15:10 / 15:30), both the 1097-s
%   sweep->hold270+LED->sweep protocol, zoom 5 (scan_params.json present).
%   11635 / 11634 volumes, 0 excluded; medial glomeruli 16 (t1) and 15-16
%   (t2) on the imclose bridge -> NaN. Good tracking both trials (pre circ
%   std 21-78 / 28-98 deg, post 28-62 / 20-45 deg; offset -40..-100 deg,
%   stable) -- post is as good or better than pre, so no LED after-effect.
%   Bump dims in the hold (rho 0.45 -> 0.18) irrespective of LED. The LED
%   drives a LARGE, PB-wide, purely excitatory response in both trials
%   (raw PB +19 % / +27 %, background ~0): per-glomerulus dF/F early +15 /
%   +33 %, late +11 / +26 %, early > 0 in 22 / 30 of 32 glomeruli (Bonf.),
%   NO glomerulus inhibited in either window; mild early > late in the
%   lateral tips (1-3, 18-19, 31-32). Square-wave time course, fast on/off,
%   ~uniform (range 3-23 % t1, 13-56 % t2 -- t2 larger medially). This is
%   the visual response to the stim LED without any CsChrimson: sets the
%   control against which the CsChrimson flies' inhibition / biphasic
%   responses must be read (overshoot_ttx_early_vs_late.m, evl_trials).
% 20261005-3_epg_syt8_lpsp_cschrimson / trial_001 (CsChrimson fly 5, 190-s
%   burn-in, zoom 5, no scan_params.json): 2110 volumes, 0 excluded, glom 16
%   on the bridge -> NaN. Textbook tracking from ~45 s on: offset +46 deg,
%   circ std 19 deg, lag 0.28 s, rho 0.4-0.7. Best-tracking CsChrimson fly
%   so far at the burn-in stage.
% 20261005-3 / trial_003 (ScanImage epoch 16:45; the trial_002 folder is an
%   aborted start with only the csv): 1097-s LED protocol, 11633 volumes, 0
%   excluded. LED pulses 430-725 s (30 s later than usual), all in the hold.
%   The 5-s protocol hold at the very end (1133-1137 s) survives the
%   detector -> fig_bump_vs_cue_scatter_blocks now splits on the LONGEST
%   hold. THE CLEAN CASE: tracking tight pre AND post (circ std 19-27 deg,
%   offset +50, identical after stim); bump parks at -90 deg (bar) in the
%   pre-LED hold; EVERY pulse flips it ~165 deg to ~+80 deg and it returns
%   within ~1 s of offset, 30/30, with rho dipping only at the transitions;
%   after the LED it sits at ~+45 deg (offset +130) until the sweeps resume
%   and re-anchor it. Per-glomerulus dF/F: lateral tips + medial (1-5,
%   15-20, 30-32) +30..+50 %, shoulders (7-13, 23-28) -20..-45 %, square-
%   wave, no early/late asymmetry (early-late r = 0.99) -> pure relocation,
%   nothing like the controls' uniform excitation. Same flip-and-return
%   signature as CsChrimson fly 1 t002 and control 20261001-1 t002.
% 20261005-3 / trial_004 (epoch 17:18, 33 min after t3; bath condition to
%   be confirmed -- looks like TTX): 11633 volumes, 0 excluded. Tracking
%   GONE pre and post (circ std ~117 deg; PVA drifts in a slow band, slope
%   0, rho 0.2-0.4 with slow waves). LED does essentially nothing: dF/F
%   +1.2 % early / +1.7 % late, no glomerulus significant, no flips, no
%   relocation map (diff image = noise). Compare fly 20261005-1 t006 (TTX)
%   which still had -8 % inhibition -- here even that is absent.
%   The small residual excitation (glomeruli 9, 11, 13, 28, +2..3.5 %)
%   sits where the trial_003 bump rested BETWEEN pulses (gap occupancy
%   1.9x vs 0.9x uniform, KW p = 0.03) and NOT where it jumped to during
%   pulses (0.04x vs 1.1x, p = 0.02) -- i.e. on the glomeruli the LED
%   suppressed in trial_003 (overshoot_ttx_occupancy_vs_phase.m with
%   evp_ttxDir/evp_salDir; perm p 0.06-0.26 on the correlations).
%   CROSS-TRIAL INDEX CAVEAT (overshoot_crosstrial_glomerulus_alignment.m):
%   masks/skeletons are per trial. For t3 vs t4 the FOV shifted 5,5 px
%   (1.8 um) and 31/32 glomeruli keep the same index after alignment
%   (median Dice 0.75) -> the result above holds with one segmentation
%   (aligned: excited 11,13,28; LED-on occupancy 0.03 vs 1.10, p = 0.03).
%   For fly 20261005-1 t5 vs t6 the shift was 11,3 px and 18/32 glomeruli
%   were OFF BY ONE in the index-matched comparison; with one segmentation
%   the on/gap correlations stay null and the pre/post hold grouping stays
%   (pre 1.75 vs 0.75, p < 0.01; post 0.11 vs 1.30, p = 0.01).
%   Mean-z version (crossTrialAlignment_meanZ): fly -3 t4 early response
%   vs t3 mean z: r = -0.50 (LED on), +0.49 (gaps, perm p 0.06) -- the
%   residual LED excitation sits on glomeruli that were ACTIVE between
%   pulses and SILENT during pulses in t3. Fly -1 TTX t6 vs t5: no
%   relation to on/gap z (|r| <= 0.4), only to the pre/post hold position.
%   STATISTICS NOW via claude\pb_hemi_corr.m: profiles averaged across
%   hemispheres (16 angles), exact 16-rotation null (p floor 1/16 =
%   0.0625), left/right numbering offset removed first (fly -1 needed -1
%   glomerulus; fly -3 none). Fly -3: early vs z(gap) r = +0.67, p = 0.062
%   (best of 16); early vs z(on) r = -0.66, p = 0.125. Fly -1: early vs
%   z(on/gap) r ~ 0, p = 1.0; late vs z(post hold) r = -0.88, p = 0.062.
% The script is written per trial; set obs_trialDir (section 0) to switch.
%
% ---- tricky things worth flagging ------------------------------------------
% 1. scan_params.json's nframes_total (132026) is NOT the real frame count:
%    it equals filesize / (100*256*2 bytes) = 132026.4, i.e. it was estimated
%    from the tif's byte size and overcounts TIFF header/IFD overhead. The
%    DAQ saw 126152 frame-clock edges = 9704 volumes x 13 exactly. On top of
%    that the tif (a BigTIFF, 6.8 GB) ends in a dangling IFD: Tiff's
%    nextDirectory() errors before lastDirectory() ever turns true, so the
%    usual count-then-read approach fails. The reader below reads in one
%    pass until the chain ends and drops any trailing partial volume.
% 2. The raw frame buffer for this trial would be ~6.5 GB as int16 if read
%    the overshoot_drug_stim_script.m way (all frames at once), so the reader
%    here sums each volume as it goes and stores the result as single
%    (100 x 256 x 9704 x 4 bytes ~ 1 GB). Summing is done in single, not in
%    int16, so 10 planes can't saturate the way an int16 sum would.
% 3. Bar position on the LED display is the ANALOG WaveSurfer channel
%    vr_ao0_copy (the LED is the digital led_stim_feedback line). Its scale
%    is 1 V per radian, 0..2pi: fitting it against the VR's own per-iteration
%    record (sync_info.mat: iter.output_vec row 5 / iter.dot_deg_ego sampled
%    at iter_times_on_daq) gives slope 0.998, offset 0.011 V, circular
%    residual 0.4 deg, and the hold plateau (4.71 V = 270 deg) matches the
%    protocol's dot_deg_ego = 630 = 270 mod 360. vr_ao1_copy and
%    vr_mfc_feedback are flat (~0 V) for the whole trial.
% 4. Timing offsets between clocks: VR time 0 = DAQ 39.4 s, imaging starts at
%    DAQ 30.0 s, and the LED alternation starts at DAQ 640.0 s -- i.e. ~12 s
%    BEFORE the bar stops sweeping (hold starts at DAQ ~651.9 s). So LED pulse
%    1 (640-645 s) happens during the sweep and pulse 2 (650-655 s) straddles
%    the sweep->hold transition; section 7's on-vs-off diff image only uses
%    pulses that fall entirely inside the hold (pulses 3-30).
% 5. Glomerulus numbering / alpha = -pi is set by which end of the mask
%    skeleton graph_sort starts from, so the bump can come out tracking the
%    bar directly or mirrored. Section 6 picks the sign that gives the
%    tighter bump-bar offset during the sweep and states which in the figure.
% 6. pb_register's rigid registration ran away on one ~10 s block near the
%    start of the sweep (104 volumes, DAQ 51.6-61.5 s, column shifts pinned
%    at the +/-45 px maxShift bound; everything else is < 5 px). Those
%    volumes show up as an all-glomerulus dip if left in, so section 4 flags
%    any volume with |shift| > 10 px and sections 5-7 NaN it out / exclude it.
%    None of them fall in the hold, so the LED results are unaffected.
%
% Run one %% section at a time.

%% 0. parameters
claudeDir = fullfile(fileparts(mfilename('fullpath')), '..', 'claude');
if isempty(which('pb_register')) && isfolder(claudeDir)
    addpath(claudeDir);
end
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort')) && isfolder(repoRootDir)
    addpath(repoRootDir);
end
if isempty(which('circ_mean')) && isfolder(fullfile(repoRootDir, 'circ_stats'))
    addpath(fullfile(repoRootDir, 'circ_stats'));
end

trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20260930-1_epg_syt8s_lpsp_cschrimson\trial_002\';
% batch override: matlab -batch "obs_trialDir='Z:\...\trial_003\'; run('overshoot_bar_stim_script.m')"
if exist('obs_trialDir', 'var') && ~isempty(obs_trialDir)
    trialDir = obs_trialDir;
end

% held-bar detection (section 4): a hold is any stretch >= min_hold_sec where
% the (median-filtered) bar channel stays within hold_flat_tol_V (0.1 V =
% 5.7 deg). Trials differ: trial_002 has one 5-min hold, trial_004 has ~18
% ten-second holds between sweeps, trial_003 has none.
min_hold_sec    = 4;   % NB: the 2-s flatness window erodes ~1 s off each end, so 5-s holds need <= 3 (use obs_min_hold_sec)
if exist('obs_min_hold_sec', 'var') && ~isempty(obs_min_hold_sec), min_hold_sec = obs_min_hold_sec; end
hold_flat_tol_V = 0.1;
lag_range_sec   = [-2 5]; % cue-bump lag scan range (s) for section 6; positive = bump follows the cue

% only used if scan_params.json is missing (it is present for this trial)
obs_fallback_scan_params = struct( ...
    'nplanes', 13, 'nplanes_valid', 10, 'flyback_planes', [10,11,12], ...
    'nchannels', 1, 'channels_saved', 1, 'px_height', 100, 'px_width', 256, ...
    'um_per_px', [0.31640625, 0.31640625, 4]);

do_motion_correction = true;
overwrite_video = false;
overwrite_reg   = false;
overwrite_mask  = false;

n_per_hemisphere = 16; % process_im convention: 16 per hemisphere -> 32 glomeruli, alpha repeated over the two halves

% mask smoothing/morphology in physical units. overshoot_drug_stim_script.m
% used sigma 2.5 um + a 60th-percentile threshold; at this FOV (0.32 um/px,
% PB arch ~25 px thick) that gave a mask ballooning well past the tissue, so
% both are tightened here.
mask_smooth_sigma_um = 1.5;
mask_open_radius_um  = 1   * 1.0125;
mask_join_line_um    = 15  * 1.0125;
mask_thresh_pct      = 80; % percentile of the smoothed, clipped mean image kept as "PB"

pathParts = strsplit(trialDir, filesep);
pathParts = pathParts(~cellfun(@isempty, pathParts));
trialName = [pathParts{end-1} filesep pathParts{end}];

%% 1. raw tif -> z-summed video (imgData_sum.mat, cached)
sumFile = fullfile(trialDir, 'imgData_sum.mat');
if ~isfile(sumFile) || overwrite_video
    fprintf('[%s] reading raw tif...\n', trialName);
    [imgData_sum, sp] = obs_read_raw_tif_summed(trialDir, obs_fallback_scan_params);
    save(sumFile, 'imgData_sum', 'sp', '-v7.3');
    fprintf('  saved %s  (%d x %d x %d volumes)\n', sumFile, size(imgData_sum,1), size(imgData_sum,2), size(imgData_sum,3));
else
    fprintf('skip (exists): %s\n', sumFile);
end

%% 2. motion correction (imgData_sum_reg.mat, cached)
regFile = fullfile(trialDir, 'imgData_sum_reg.mat');
if do_motion_correction && (~isfile(regFile) || overwrite_reg)
    S = load(sumFile, 'imgData_sum');
    fprintf('[%s] motion-correcting (%d volumes)...\n', trialName, size(S.imgData_sum,3));
    [imgData_sum_reg, regShifts, regDiag] = pb_register(S.imgData_sum, ...
        'smoothWindow', 5, 'maxShift', [30 45], 'verbose', false, 'figureVisible', false); %#ok<ASGLU>
    close(regDiag.figHandles);
    save(regFile, 'imgData_sum_reg', 'regShifts', '-v7.3');
    fprintf('  saved %s\n', regFile);
else
    fprintf('skip (exists): %s\n', regFile);
end

%% 3. PB mask (mask.mat + mask_qc.png, cached)
maskFile = fullfile(trialDir, 'mask.mat');
S = load(sumFile, 'sp');
umPerPx = S.sp.um_per_px(1);
if ~isfile(maskFile) || overwrite_mask
    S = load(regFile, 'imgData_sum_reg');
    projMean = mean(S.imgData_sum_reg, 3);

    sigma_px = max(mask_smooth_sigma_um / umPerPx, 0.5);
    open_radius_px = max(round(mask_open_radius_um / umPerPx), 1);
    join_line_px = max(round(mask_join_line_um / umPerPx), 3);
    fprintf('mask params in px (um_per_px=%.4f): sigma=%.2f, open_radius=%d, join_line=%d\n', ...
        umPerPx, sigma_px, open_radius_px, join_line_px);

    projSmooth = imgaussfilt(projMean, sigma_px);
    top_pct = prctile(projSmooth, 98, 'all');
    bot_pct = prctile(projSmooth, 5, 'all');
    projClipped = min(max(projSmooth, bot_pct), top_pct);

    mask = projClipped > prctile(projClipped(:), mask_thresh_pct);
    mask = imopen(mask, strel('disk', open_radius_px));
    mask = imfill(mask, 'holes');
    areas = sort([regionprops(mask).Area], 'descend');
    if numel(areas) < 2 || (areas(1) / areas(2)) > 1.5
        mask = bwareafilt(mask, 1);
        maskCore = mask;
    else
        mask = bwareafilt(mask, 2);
        maskCore = mask; % the two real blobs, before bridging
        mask = imclose(mask, strel('line', join_line_px, 0));
    end
    mask = imfill(mask, 'holes');
    maskCore = maskCore & mask;
    % maskCore = pixels with real PB signal; mask may additionally contain the
    % bridge drawn by imclose between two separate hemisphere blobs (needed so
    % the skeleton is one path). Glomeruli that fall only on the bridge have
    % no data and are NaN'd in section 5.
    save(maskFile, 'mask', 'maskCore');
    fprintf('saved %s\n', maskFile);
else
    fprintf('skip (exists): %s -- redrawing the QC figure only\n', maskFile);
    M = load(maskFile); mask = M.mask;
    if isfield(M, 'maskCore'), maskCore = M.maskCore; else, maskCore = mask; end % older mask.mat
    clear M
    S = load(regFile, 'imgData_sum_reg');
    projMean = mean(S.imgData_sum_reg, 3);
    clear S
end

% QC: mask contour + the glomerulus numbering it implies (same overlay as
% led_stim_berg4_script.m), so a bad skeleton is visible here rather than
% as a scrambled heatmap later. Always redrawn (cheap) so the .fig/.png stay
% in sync with the mask on disk.
    figure(1); clf
    imagesc(projMean); axis equal tight; colormap(bone); hold on
    contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
    if any(mask(:) & ~maskCore(:))
        contour(maskCore, [0.5 0.5], 'c', 'LineWidth', 1); % real-signal region; outside it = imclose bridge
    end
    try
        qcClusters = obs_pb_glomeruli(mask, n_per_hemisphere);
        nC = 2 * n_per_hemisphere;
        cen = zeros(nC, 2);
        for c = 1:nC
            [cy, cx] = find(qcClusters == c);
            cen(c,:) = [mean(cx), mean(cy)];
        end
        plot(cen(:,1), cen(:,2), 'y-', 'LineWidth', 1)
        plot(cen(:,1), cen(:,2), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 4)
        for c = 1:nC
            text(cen(c,1), cen(c,2), num2str(c), 'Color', 'g', 'FontSize', 6, ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
    catch ME
        fprintf('  note: glomerulus QC overlay failed (%s) -- mask itself still saved.\n', ME.message);
    end
    title([trialName ' mask QC'], 'Interpreter', 'none');
    obs_export(figure(1), fullfile(trialDir, 'mask_qc.png'));
M = load(maskFile); mask = M.mask;
if isfield(M, 'maskCore'), maskCore = M.maskCore; else, maskCore = mask; end
clear M

%% 4. sync from the WaveSurfer h5: volume times, LED, bar position, hold window
h5list = dir(fullfile(trialDir, '*.h5'));
if numel(h5list) ~= 1
    error('overshoot_bar_stim_script:h5', 'Expected exactly one .h5 in %s, found %d.', trialDir, numel(h5list));
end
sync = obs_h5_sync(fullfile(h5list(1).folder, h5list(1).name));

S = load(regFile, 'imgData_sum_reg', 'regShifts');
img = single(S.imgData_sum_reg);
% Volumes where rigid registration ran away (|shift| > max_ok_shift_px --
% in this trial 104 volumes in one ~10 s block at 51.6-61.5 s DAQ pinned to
% the +/-45 px maxShift bound) are junk for glomerulus traces: they get
% NaN'd out of f_cluster / the PVA below and excluded from the diff image.
regShiftXY = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), S.regShifts(:), 'UniformOutput', false));
max_ok_shift_px = 10;
badVol = any(abs(regShiftXY) > max_ok_shift_px, 2)';
clear S
nVol = size(img, 3);
if numel(sync.t_volume) ~= nVol
    fprintf('note: video has %d volumes but the h5 has %d volume-clock edges -- truncating to the shorter.\n', nVol, numel(sync.t_volume));
    n = min(nVol, numel(sync.t_volume));
    img = img(:,:,1:n); nVol = n; badVol = badVol(1:n);
    sync.t_volume = sync.t_volume(1:n);
end
t_vol = sync.t_volume(:)'; % DAQ seconds, one per imaging volume
fprintf('%d/%d volumes (%.1f%%) have a registration shift > %d px and are excluded (DAQ %.1f-%.1f s)\n', ...
    sum(badVol), nVol, 100*mean(badVol), max_ok_shift_px, min([t_vol(badVol) NaN]), max([t_vol(badVol) NaN]));

led_vol     = interp1(sync.t_fine, double(sync.led), t_vol, 'nearest', 'extrap') > 0.5;
bar_rad_vol = mod(interp1(sync.t_fine, sync.bar_rad, t_vol, 'nearest', 'extrap'), 2*pi);
bar_deg_vol = rad2deg(bar_rad_vol);

% Held-bar segments, found from the analog trace itself (not from protocol
% timestamps, since VR time 0 sits ~40-50 s into the DAQ record): every
% stretch >= min_hold_sec where the median-filtered bar stays within
% hold_flat_tol_V. holdSegs is n x 2 (DAQ s), possibly empty.
holdSegs = obs_flat_segments(sync.t_fine, sync.bar_rad, sync.vr_start, sync.fs, hold_flat_tol_V, min_hold_sec);
% drop the dark display after the protocol ends (bar channel parks at 0 V,
% which otherwise looks like a long "hold" at 0 deg) and anything after the
% last imaging volume
holdSegs = holdSegs(holdSegs(:,1) < sync.vr_end - 0.5 & holdSegs(:,1) < t_vol(end), :);
inHold_vol = false(1, nVol);
pulseInHold = false(size(sync.led_on));
for k = 1:size(holdSegs, 1)
    inHold_vol  = inHold_vol  | (t_vol >= holdSegs(k,1) & t_vol <= holdSegs(k,2));
    pulseInHold = pulseInHold | (sync.led_on >= holdSegs(k,1) & sync.led_off <= holdSegs(k,2));
end
fprintf('%d held-bar segment(s) >= %g s:\n', size(holdSegs, 1), min_hold_sec);
for k = 1:size(holdSegs, 1)
    segAng = rad2deg(median(sync.bar_rad(sync.t_fine >= holdSegs(k,1) & sync.t_fine <= holdSegs(k,2))));
    fprintf('  %7.1f - %7.1f s (%6.1f s) at %6.1f deg\n', holdSegs(k,1), holdSegs(k,2), diff(holdSegs(k,:)), segAng);
end
if isempty(sync.led_on)
    fprintf('LED: no pulses in this trial.\n');
else
    fprintf('LED: %d pulses of %.2f s (%.1f-%.1f s DAQ), %d fully inside a held-bar segment.\n', ...
        numel(sync.led_on), median(sync.led_off - sync.led_on), sync.led_on(1), sync.led_off(end), sum(pulseInHold));
end
% zoom window for the "hold" figures: spans the holds (with margin) unless
% that is essentially the whole trial, in which case the zoom is skipped
zoomWin = [];
if ~isempty(holdSegs)
    zoomWin = [max(holdSegs(1,1) - 30, t_vol(1)), min(holdSegs(end,2) + 10, t_vol(end))];
    if diff(zoomWin) > 0.9 * (t_vol(end) - t_vol(1)); zoomWin = []; end
end

% per-frame background subtraction: each frame minus its own mean OUTSIDE the
% mask (same convention as overshoot_drug_stim_script.m). Matters here
% because any LED bleed-through into the frame would otherwise show up as a
% stim-locked change in every glomerulus at once.
img_2d = reshape(img, [], nVol);
bg = mean(img_2d(~mask(:), :), 1);
img = img - reshape(single(bg), 1, 1, nVol);
clear img_2d

%% 5. glomeruli -> z-scored activity heatmap, with bar position and LED
clusterIdx = obs_pb_glomeruli(mask, n_per_hemisphere);
nClusters  = 2 * n_per_hemisphere;
img_2d = reshape(img, [], nVol);
f_cluster = nan(nClusters, nVol);
nCorePx = zeros(nClusters, 1);
for c = 1:nClusters
    pix = clusterIdx(:) == c & maskCore(:); % only pixels with real PB signal (not the imclose bridge)
    nCorePx(c) = nnz(pix);
    if nCorePx(c) > 0, f_cluster(c,:) = mean(img_2d(pix, :), 1); end
end
clear img_2d
% glomeruli that are (almost) entirely on the bridge between two separate
% hemisphere blobs carry no signal -> drop them from z-scoring and the PVA
emptyGlom = nCorePx < max(10, 0.2 * median(nCorePx(nCorePx > 0)));
if any(emptyGlom)
    f_cluster(emptyGlom, :) = NaN;
    fprintf('%d/%d glomeruli have (almost) no real-signal pixels (%s) -> NaN, excluded from the PVA\n', ...
        nnz(emptyGlom), nClusters, mat2str(find(emptyGlom)'));
end
f_cluster(:, badVol) = NaN; % registration-failure volumes (section 4)
% z-score per glomerulus over the whole trial. This equals zscore(dF/F)
% exactly (dF/F is an affine transform of F within a glomerulus), so no F0
% choice is needed here -- same quantity process_im calls im.z. Done by hand
% rather than zscore() so the NaN'd volumes are skipped instead of poisoning
% every glomerulus's mean/std.
z_cluster = (f_cluster - mean(f_cluster, 2, 'omitnan')) ./ std(f_cluster, 0, 2, 'omitnan');
alpha = repmat(linspace(-pi, pi - 2*pi/n_per_hemisphere, n_per_hemisphere), 1, 2); % process_im convention

rwb = obs_redwhiteblue(256);
bar_deg_plot = bar_deg_vol; bar_deg_plot(abs(diff([bar_deg_plot bar_deg_plot(end)])) > 180) = NaN; % break the line at wraps

for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1
        xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial';
    else
        xl = zoomWin; tag = 'hold';
    end
    figure(2); clf
    set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(6,1,1:3);
    imagesc(t_vol, 1:nClusters, z_cluster, 'AlphaData', double(~isnan(z_cluster))); colormap(ax1, rwb); clim([-3 3]); hold on
    set(ax1, 'Color', [0.7 0.7 0.7]) % NaN'd (registration-failure) volumes render gray, not as the colormap's first color
    obs_mark_holds(holdSegs)
    ylabel('glomerulus')
    title(sprintf('%s: z-scored glomerulus activity (background-subtracted F, %d glomeruli); dashed = held-bar segments (%d)', trialName, nClusters, size(holdSegs,1)), 'Interpreter', 'none')
    cb = colorbar(ax1); cb.Label.String = 'z-score';
    ax2 = subplot(6,1,4:5);
    obs_shade_pulses(sync.led_on, sync.led_off, [0 360]); hold on
    plot(t_vol, bar_deg_plot, 'k', 'LineWidth', 1)
    ylim([0 360]); yticks(0:90:360); ylabel('bar position (deg)')
    title('bar position on the LED display (red shading = stim LED on)')
    ax3 = subplot(6,1,6);
    plot(t_vol, double(led_vol), 'r', 'LineWidth', 1); ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'off','on'})
    ylabel('stim LED'); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3], 'x'); xlim(xl)
    ax1.Position([1 3]) = ax2.Position([1 3]); % undo the colorbar's shrink so the time axes stay aligned
    outFile = fullfile(trialDir, sprintf('glomHeatmap_zscore_%s.png', tag));
    obs_export(figure(2), outFile);
end

%% 6. bump position (PVA of the z-scored glomerulus activity) vs. bar position
% process_im's estimate: each glomerulus is a unit vector at angle alpha
% scaled by its z-score; the population vector's angle is the bump position
% (mu), its length the bump strength (rho). No extra smoothing (im_win = 1,
% as in lpsp_cschrimson_reredo_script.m).
[x_tmp, y_tmp] = pol2cart(alpha, z_cluster');
[mu, rho] = cart2pol(mean(x_tmp, 2, 'omitnan'), mean(y_tmp, 2, 'omitnan')); % omitnan: skips empty (bridge) glomeruli
mu = mu(:)'; rho = rho(:)';
rho_thresh = 0.1; % this repo's usual low-confidence cutoff (rho_thresh in the reredo/debug scripts)

% Which glomerulus is alpha = -pi (and whether the numbering runs with or
% against the display) is set by the mask skeleton, not anatomy, so the bump
% can track the bar directly or mirrored. Pick whichever sign gives the
% tighter bump-bar offset during the sweep and say so in the figure.
bar_pm = angle(exp(1i * bar_rad_vol)); % bar wrapped to (-pi, pi], like mu
sweepOK = ~inHold_vol & rho >= rho_thresh & t_vol > sync.vr_start;
cs_same = obs_circ_std(mu(sweepOK) - bar_pm(sweepOK));
cs_flip = obs_circ_std(mu(sweepOK) + bar_pm(sweepOK));
if cs_flip < cs_same
    bar_sign = -1; signNote = 'glomerulus order is MIRRORED relative to the display (bump ~ -bar)';
else
    bar_sign = 1;  signNote = 'glomerulus order runs WITH the display (bump ~ +bar)';
end
fprintf('bump-vs-bar: circ std of offset during sweep = %.2f rad (same sign) vs %.2f rad (flipped) -> %s\n', cs_same, cs_flip, signNote);

% Cue-bump lag. The indicator / circuit delay the bump by well under a
% second, but at 15 deg/s sweeps that is a 10-15 deg bias whose sign flips
% with sweep direction, so the unlagged offset trace zig-zags. Scan lags
% (cue at t - lag vs bump at t, interpolated on the unit circle so wrap-
% around is handled) and keep the one that tightens the sweep offset most.
% Everything downstream (offset trace, bar line in the figure, saved
% results) uses the lagged cue; lag 0 is used if the scan finds nothing.
bar_c = exp(1i * bar_rad_vol);
dt_vol = median(diff(t_vol));
lagScan = lag_range_sec(1):dt_vol:lag_range_sec(2);
lagScore = nan(size(lagScan));
for L = 1:numel(lagScan)
    barL = angle(interp1(t_vol, bar_c, t_vol - lagScan(L), 'linear', 'extrap'));
    lagScore(L) = obs_circ_std(mu(sweepOK) - bar_sign * barL(sweepOK));
end
[~, iBest] = min(lagScore);
lag_sec = lagScan(iBest);
if isnan(lag_sec) || lag_sec < 0
    % a bump that LEADS the cue is not physical; it only happens when the
    % tracking is too poor for the scan to mean anything -> no lag
    lag_sec = 0;
end
bar_rad_lag = angle(interp1(t_vol, bar_c, t_vol - lag_sec, 'linear', 'extrap'));
bar_pm = angle(exp(1i * bar_rad_lag)); % lagged cue, wrapped to (-pi, pi] like mu
iZero = find(abs(lagScan) == min(abs(lagScan)), 1);
iUsed = find(abs(lagScan - lag_sec) == min(abs(lagScan - lag_sec)), 1);
fprintf('  cue-bump lag scan (%.1f..%.1f s): using lag = %.2f s (sweep offset circ std %.1f deg, vs %.1f deg at lag 0)%s\n', ...
    lag_range_sec, lag_sec, rad2deg(lagScore(iUsed)), rad2deg(lagScore(iZero)), ...
    repmat(sprintf(' -- scan minimum was at %.2f s, rejected as non-physical', lagScan(iBest)), 1, lagScan(iBest) < 0));

offset = angle(exp(1i * (mu - bar_sign * bar_pm))); % bump - (signed, lagged) bar, in (-pi, pi]
fprintf('  mean offset during sweep = %.1f deg (circ std %.1f deg); hold, LED off = %.1f deg (%.1f); hold, LED on = %.1f deg (%.1f)\n', ...
    rad2deg(obs_circ_mean(offset(sweepOK))), rad2deg(obs_circ_std(offset(sweepOK))), ...
    rad2deg(obs_circ_mean(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), rad2deg(obs_circ_std(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), ...
    rad2deg(obs_circ_mean(offset(inHold_vol &  led_vol & rho >= rho_thresh))), rad2deg(obs_circ_std(offset(inHold_vol &  led_vol & rho >= rho_thresh))));

% Per-volume results for downstream figure scripts (ugly_figures/scripts/
% fig_bump_vs_cue_scatter_blocks.m etc.), so they need not redo the PVA.
bump = struct('t_vol', t_vol, 'mu', mu, 'rho', rho, 'bar_rad', bar_rad_vol, ...
    'bar_sign', bar_sign, 'bar_rad_lag', bar_rad_lag, 'lag_sec', lag_sec, 'offset', offset, 'led', led_vol, 'inHold', inHold_vol, ...
    'badVol', badVol, 'holdSegs', holdSegs, 'led_on', sync.led_on, 'led_off', sync.led_off, ...
    'vr_start', sync.vr_start, 'vr_end', sync.vr_end, 'rho_thresh', rho_thresh, ...
    'z_cluster', z_cluster, 'f_cluster', f_cluster, 'clusterIdx', clusterIdx, 'nCorePx', nCorePx, ...
    'alpha', alpha, 'trialName', trialName, 'signNote', signNote); %#ok<NASGU>
% f_cluster = per-glomerulus mean of the background-subtracted F (raw units),
% so downstream scripts can compute dF/F; z_cluster is its z-score.
save(fullfile(trialDir, 'bump_results.mat'), 'bump');
fprintf('saved %s\n', fullfile(trialDir, 'bump_results.mat'));

mu_plot = rad2deg(mu); mu_plot(rho < rho_thresh) = NaN;
mu_plot(abs(diff([mu_plot mu_plot(end)])) > 180) = NaN;
barS_plot = rad2deg(angle(exp(1i * bar_sign * bar_pm)));
barS_plot(abs(diff([barS_plot barS_plot(end)])) > 180) = NaN;
off_plot = rad2deg(offset); off_plot(rho < rho_thresh) = NaN;

for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1
        xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial';
    else
        xl = zoomWin; tag = 'hold';
    end
    figure(3); clf
    set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(3,1,1);
    obs_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, barS_plot, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5)
    plot(t_vol, mu_plot, 'b.', 'MarkerSize', 4)
    obs_mark_holds(holdSegs)
    ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg)')
    legend({sprintf('bar (signed, lagged %.2f s)', lag_sec), 'bump (PVA), rho >= 0.1'}, 'Location', 'northeastoutside')
    title(sprintf('%s: bump position vs bar position -- %s; bar lagged by %.2f s (best fit during sweeps); red = stim LED on', ...
        trialName, signNote, lag_sec), 'Interpreter', 'none')
    ax2 = subplot(3,1,2);
    obs_shade_pulses(sync.led_on, sync.led_off, [0 max(rho)*1.05]); hold on
    plot(t_vol, rho, 'k'); yline(rho_thresh, 'r:');
    ylabel('bump strength (rho)')
    ax3 = subplot(3,1,3);
    obs_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, off_plot, 'b.', 'MarkerSize', 4); yline(0, 'k:');
    obs_mark_holds(holdSegs)
    ylim([-180 180]); yticks(-180:90:180); ylabel(sprintf('bump - bar(t - %.2f s) (deg)', lag_sec)); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3], 'x'); xlim(xl)
    for ax = [ax2 ax3]; ax.Position([1 3]) = ax1.Position([1 3]); end % match the legend-shrunk top axis
    outFile = fullfile(trialDir, sprintf('bumpVsBar_%s.png', tag));
    obs_export(figure(3), outFile);
end

%% 7. hold-period mean diff image: stim LED on - stim LED off
% Only pulses that fall entirely inside the hold count as "on"; "off" is the
% interleaved LED-off time between those pulses (so both sides come from the
% same held-bar block, not from the sweep or from after the VR ended).
if isempty(sync.led_on)
    fprintf('no LED pulses in this trial -- skipping the stim on/off diff image.\n');
    return
end
if sum(pulseInHold) >= 3
    pulseIdx = find(pulseInHold); offMask = inHold_vol;
    diffNote = sprintf('bar held: %d pulses inside held-bar segments', numel(pulseIdx));
else
    % not enough pulses inside holds (e.g. LED running through the sweeps):
    % fall back to every pulse, with the LED-off time between them as control
    pulseIdx = 1:numel(sync.led_on); offMask = true(1, nVol);
    diffNote = sprintf('ALL %d pulses (only %d inside held-bar segments)', numel(pulseIdx), sum(pulseInHold));
end
onIdx = false(1, nVol);
for k = pulseIdx
    onIdx = onIdx | (t_vol >= sync.led_on(k) & t_vol < sync.led_off(k));
end
offIdx = offMask & ~led_vol & t_vol >= sync.led_on(pulseIdx(1)) & t_vol <= sync.led_off(pulseIdx(end));
onIdx = onIdx & ~badVol; offIdx = offIdx & ~badVol;
if sum(onIdx) < 5 || sum(offIdx) < 5
    % e.g. the LED fired after imaging stopped (20261005-1 trial_003: one
    % pulse at DAQ 400 s, imaging ended at 388 s) -> nothing to compare
    fprintf('diff image skipped: only %d LED-on / %d LED-off imaging volumes (%s; LED %.0f-%.0f s, imaging %.0f-%.0f s).\n', ...
        sum(onIdx), sum(offIdx), diffNote, sync.led_on(1), sync.led_off(end), t_vol(1), t_vol(end));
    return
end
meanOn  = mean(img(:,:,onIdx), 3);
meanOff = mean(img(:,:,offIdx), 3);
diffImg = meanOn - meanOff;
projMean = mean(img, 3);
fprintf('diff image (%s): %d LED-on volumes vs %d LED-off volumes\n', diffNote, sum(onIdx), sum(offIdx));

figure(4); clf
set(gcf, 'Position', [50 50 1400 450], 'Color', 'w')
subplot(1,2,1)
imagesc(projMean); axis equal tight; colormap(gca, bone); hold on
contour(mask, [0.5 0.5], 'r', 'LineWidth', 1); title('mean image (background-subtracted) + mask'); colorbar
subplot(1,2,2)
dlim = prctile(abs(diffImg(mask)), 99);
imagesc(diffImg); axis equal tight; colormap(gca, rwb); clim([-dlim dlim]); hold on
contour(mask, [0.5 0.5], 'k', 'LineWidth', 0.75); colorbar
title(sprintf('mean(LED on) - mean(LED off): %s; %d on / %d off volumes', diffNote, sum(onIdx), sum(offIdx)), 'FontSize', 9)
sgtitle([trialName ': stim on - off diff image (red = brighter with LED on, blue = dimmer)'], 'Interpreter', 'none')
outFile = fullfile(trialDir, 'holdDiff_stimOn_minus_off.png');
obs_export(figure(4), outFile);

%% local functions

function s = obs_h5_sync(h5Path)
% Everything this script needs from the WaveSurfer file, on two clocks:
% the fine DAQ clock (t_fine, fs) and the imaging volume clock (t_volume =
% si_volumeclock rising edges). LED = led_stim_feedback digital line;
% bar = vr_ao0_copy analog line, which is the bar angle at 1 V per radian
% (checked against the VR's own per-iteration record in sync_info.mat:
% slope 0.998, offset 0.011 V, circular residual 0.4 deg).
info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
if numel(sweepNames) ~= 1
    error('overshoot_bar_stim_script:h5sweep', 'Expected exactly one sweep group in %s, found %d.', h5Path, numel(sweepNames));
end
s.fs = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
aiNames = strtrim(string(h5read(h5Path, '/header/AIChannelNames')));
coef = double(h5read(h5Path, '/header/AIScalingCoefficients')); % polynomial (c0..c3) per AI channel, raw counts -> volts

digi = int32(h5read(h5Path, [sweepNames{1} '/digitalScans']));
volCh = find(diNames == "si_volumeclock", 1);
ledCh = find(diNames == "led_stim_feedback", 1);
barCh = find(aiNames == "vr_ao0_copy", 1);
if isempty(volCh) || isempty(ledCh) || isempty(barCh)
    error('overshoot_bar_stim_script:h5chan', 'Missing si_volumeclock / led_stim_feedback / vr_ao0_copy in %s (DI: %s; AI: %s).', ...
        h5Path, strjoin(diNames, ', '), strjoin(aiNames, ', '));
end
N = numel(digi);
s.t_fine = (0:N-1)' / s.fs;
volBit = bitget(digi, volCh);
s.t_volume = s.t_fine(find(diff(volBit) > 0) + 1);
s.led = bitget(digi, ledCh) > 0;
onI  = find(diff([0; double(s.led)]) == 1);
offI = find(diff([double(s.led); 0]) == -1);
m = min(numel(onI), numel(offI));
s.led_on  = s.t_fine(onI(1:m))';
s.led_off = s.t_fine(offI(1:m))';

ana = h5read(h5Path, [sweepNames{1} '/analogScans']);
if size(ana, 1) < size(ana, 2); ana = ana'; end
s.bar_rad = polyval(flipud(coef(:, barCh)), double(ana(:, barCh))); % volts == radians for this channel
% VR start: first time the bar channel leaves ~0 V (the display is dark
% before the protocol starts). Used only to exclude the pre-VR dead time.
s.vr_start = s.t_fine(find(s.bar_rad > 0.2, 1));
s.vr_end   = s.t_fine(find(s.bar_rad > 0.2, 1, 'last')); % display goes dark (0 V) after the protocol ends
end

function segs = obs_flat_segments(t, bar, tStart, fs, tolV, minSec)
% Held-bar segments: stretches (after tStart) where the bar channel, after a
% 0.25 s median filter to kill ADC quantization noise, moves less than tolV
% within any 2 s window, lasting at least minSec. Returns n x 2 [start end]
% in the units of t. (Differentiating the raw 2 kHz trace does NOT work:
% one ADC count of jitter is ~0.6 V/s.)
barF = movmedian(bar(:), round(0.25 * fs));
win = round(2 * fs);
flat = (movmax(barF, win) - movmin(barF, win)) < tolV & t(:) > tStart;
d = diff([0; flat; 0]);
s = find(d == 1); e = find(d == -1) - 1;
keep = (e - s + 1) / fs >= minSec;
segs = [t(s(keep)), t(e(keep))];
if isempty(segs); segs = zeros(0, 2); end
end

function obs_export(fig, outFile)
% exportgraphics otherwise bakes the interactive axes toolbar into the PNG
% when run headless with linked axes ("Exported image displays axes toolbar").
set(findall(fig, 'Type', 'axestoolbar'), 'Visible', 'off');
drawnow
exportgraphics(fig, outFile, 'Resolution', 150);
fprintf('saved %s\n', outFile);
% Also save an interactive .fig next to the PNG so it can be reopened and
% zoomed in MATLAB. Headless figures are created invisible, so force Visible
% on in the saved file or openfig() would bring up an invisible window.
[p, n] = fileparts(outFile);
figFile = fullfile(p, [n '.fig']);
vis = get(fig, 'Visible');
set(fig, 'Visible', 'on');
savefig(fig, figFile, 'compact');
set(fig, 'Visible', vis);
fprintf('saved %s\n', figFile);
end

function obs_mark_holds(segs)
% Dashed vertical lines at the start/end of every held-bar segment.
for k = 1:size(segs, 1)
    xline(segs(k,1), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
    xline(segs(k,2), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off');
end
end

function obs_shade_pulses(onsets, offsets, yl)
% Light red patches over each [onset offset] interval on the current axes.
for k = 1:numel(onsets)
    patch([onsets(k) offsets(k) offsets(k) onsets(k)], [yl(1) yl(1) yl(2) yl(2)], [1 0.55 0.55], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off');
end
ylim(yl)
end

function m = obs_circ_mean(x)
m = angle(mean(exp(1i * x(:)), 'omitnan'));
end

function s = obs_circ_std(x)
x = x(~isnan(x));
if isempty(x); s = NaN; return; end
R = abs(mean(exp(1i * x(:))));
s = sqrt(-2 * log(max(R, eps)));
end

function sp = obs_scan_params(trialDir, fallback)
jsonFile = fullfile(trialDir, 'scan_params.json');
if isfile(jsonFile)
    sp = jsondecode(fileread(jsonFile));
    sp.provenance = 'scan_params.json';
else
    sp = fallback;
    sp.provenance = 'FALLBACK (scan_params.json missing)';
    fprintf('  note: %s has no scan_params.json; using fallback scan params', trialDir);
    % Cross-check / correct the fallback against the ScanImage header text at
    % the start of the tif (plain "SI.x.y = value" lines). The zoom in
    % particular differs between flies (4 on 20260930, 5 on 20261001) and sets
    % um_per_px, which the mask parameters (in um) depend on.
    tifList = dir(fullfile(trialDir, '*.tif'));
    if ~isempty(tifList)
        fid = fopen(fullfile(tifList(1).folder, tifList(1).name), 'r');
        hdr = fread(fid, 60000, 'uint8=>char')';
        fclose(fid);
        getNum = @(key) str2double(regexp(hdr, ['SI\.' regexptranslate('escape', key) '\s*=\s*([-0-9.eE]+)'], 'tokens', 'once'));
        nWith   = getNum('hStackManager.numFramesPerVolumeWithFlyback');
        nValid  = getNum('hStackManager.numFramesPerVolume');
        nFly    = getNum('hFastZ.numDiscardFlybackFrames');
        lines   = getNum('hRoiManager.linesPerFrame');
        pixels  = getNum('hRoiManager.pixelsPerLine');
        zoom    = getNum('hRoiManager.scanZoomFactor');
        if ~isnan(nWith) && ~isnan(nValid) && ~isnan(nFly) && nWith == nValid + nFly
            sp.nplanes = nWith; sp.nplanes_valid = nValid;
            sp.flyback_planes = nValid:(nWith - 1); % 0-indexed, flyback frames come last
        end
        if ~isnan(lines),  sp.px_height = lines;  end
        if ~isnan(pixels), sp.px_width  = pixels; end
        if ~isnan(zoom)
            % fallback um_per_px is calibrated at zoom 4 (0.3164 um/px); FOV
            % scales inversely with zoom.
            um_xy = fallback.um_per_px(1) * 4 / zoom;
            sp.um_per_px = [um_xy, um_xy, fallback.um_per_px(3)];
            sp.zoom = zoom;
        end
        fprintf(' (tif header: %dx%d px, %d planes = %d valid + %d flyback, zoom %g -> %.4f um/px)', ...
            sp.px_height, sp.px_width, sp.nplanes, sp.nplanes_valid, sp.nplanes - sp.nplanes_valid, zoom, sp.um_per_px(1));
    end
    fprintf('.\n');
end
end

function [imgSum, sp] = obs_read_raw_tif_summed(trialDir, fallbackSp)
% Like overshoot_drug_stim_script.m's ovs_read_raw_tif_summed, but it does
% NOT go through MATLAB's Tiff class: libtiff indexes directories with a
% 16-bit counter, so Tiff/nextDirectory dies ("Unable to read the next
% directory") after exactly 65535 frames -- on this 126k-frame trial that
% silently drops everything after ~510 s, i.e. the whole bar-hold + LED
% block. Instead this walks the BigTIFF IFD chain itself with fread (each
% IFD: uint64 entry count, 20-byte entries, uint64 next-IFD offset; tags 273
% / 279 give each frame's data offset / byte count) and reads every frame's
% pixels by offset, summing each volume as it goes into a single-precision
% stack (never holding all raw frames at once). Any trailing partial volume
% is dropped. Validated frame-for-frame against Tiff.read() on frames inside
% the 65535 range before use. (ScanImageTiffReader, on this machine via the
% flyg repo, would also work but loads the entire 6.5 GB stack at once.)
sp = obs_scan_params(trialDir, fallbackSp);
if sp.nchannels ~= 1
    error('overshoot_bar_stim_script:multiChannel', '%s has %d saved channels; this script assumes 1.', trialDir, sp.nchannels);
end
validPlanes = setdiff(1:sp.nplanes, sp.flyback_planes(:)' + 1); % json flyback is 0-indexed
if numel(validPlanes) ~= sp.nplanes_valid
    error('overshoot_bar_stim_script:flybackMismatch', ...
        '%s: nplanes_valid=%d but nplanes/flyback_planes leaves %d.', trialDir, sp.nplanes_valid, numel(validPlanes));
end

tifList = dir(fullfile(trialDir, '*.tif'));
if numel(tifList) ~= 1
    error('overshoot_bar_stim_script:tif', 'Expected exactly one .tif in %s, found %d.', trialDir, numel(tifList));
end
tifPath = fullfile(tifList(1).folder, tifList(1).name);

framesPerVol = sp.nplanes * sp.nchannels;
H = sp.px_height; W = sp.px_width;
d = dir(tifPath); fileSize = d.bytes;
% json's nframes_total is a filesize-based OVERestimate (header note 1), so
% it is a safe upper bound for preallocation; trimmed to what was read below.
% The fallback params (no scan_params.json) have no nframes_total -> use the
% same filesize estimate the acquisition code uses.
if ~isfield(sp, 'nframes_total') || isempty(sp.nframes_total)
    sp.nframes_total = ceil(fileSize / (H * W * 2));
end
nVolAlloc = ceil(sp.nframes_total / framesPerVol) + 10;
fid = fopen(tifPath, 'r', 'ieee-le');
cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid, 2, 'uint8=>char')';
version = fread(fid, 1, 'uint16');
if ~strcmp(magic, 'II') || version ~= 43
    error('overshoot_bar_stim_script:notBigTiff', ...
        '%s is not a little-endian BigTIFF (magic ''%s'', version %d) -- this reader only handles ScanImage BigTIFFs.', tifPath, magic, version);
end
offsetSize = fread(fid, 1, 'uint16'); % 8 for BigTIFF
fread(fid, 1, 'uint16');              % reserved, always 0
if offsetSize ~= 8
    error('overshoot_bar_stim_script:notBigTiff', '%s: BigTIFF offset size is %d, expected 8.', tifPath, offsetSize);
end
p = fread(fid, 1, 'uint64'); % first IFD offset (bytes 8-15 of the header)

imgSum = zeros(H, W, nVolAlloc, 'single');
vol = zeros(H, W, framesPerVol, 'single');
tic
k = 0; v = 0; f = 0; chainOk = true; checkedFormat = false;
while p > 0
    if p + 8 > fileSize
        chainOk = false;
        fprintf('  note: next-IFD pointer (%d) points past end of file (%d) after %d frames -- treating that as end of file.\n', p, fileSize, k);
        break
    end
    fseek(fid, p, 'bof');
    nEntries = fread(fid, 1, 'uint64');
    if isempty(nEntries) || nEntries < 1 || nEntries > 500
        chainOk = false;
        fprintf('  note: IFD at %d has an implausible entry count after %d frames -- treating that as end of file.\n', p, k);
        break
    end
    ent = fread(fid, [10 nEntries], 'uint16=>double'); % 20-byte entries: tag, type, count(uint64), value/offset(uint64)
    tags   = ent(1,:);
    counts = ent(3,:) + ent(4,:)*2^16 + ent(5,:)*2^32 + ent(6,:)*2^48;
    vals   = ent(7,:) + ent(8,:)*2^16 + ent(9,:)*2^32 + ent(10,:)*2^48;
    nextIFD = fread(fid, 1, 'uint64');
    dataOff = vals(tags == 273); nStrips = counts(tags == 273);
    if isempty(dataOff) || nStrips ~= 1
        error('overshoot_bar_stim_script:strips', 'Frame %d: expected a single strip per frame (StripOffsets count 1), got %d.', k+1, nStrips);
    end
    if ~checkedFormat
        w = vals(tags == 256); h = vals(tags == 257); bps = vals(tags == 258); fmt = vals(tags == 339);
        if w ~= W || h ~= H || bps ~= 16 || (~isempty(fmt) && fmt ~= 2)
            error('overshoot_bar_stim_script:format', ...
                'tif frame is %dx%d, %d bits, SampleFormat %s; expected %dx%d int16.', w, h, bps, mat2str(fmt), W, H);
        end
        checkedFormat = true;
    end
    fseek(fid, dataOff, 'bof');
    frame = fread(fid, [W H], 'int16=>single'); % TIFF rows are contiguous -> read as W x H, transpose below
    if numel(frame) ~= W*H
        chainOk = false;
        fprintf('  note: frame %d is truncated (%d of %d pixels) -- treating that as end of file.\n', k+1, numel(frame), W*H);
        break
    end
    k = k + 1; f = f + 1;
    vol(:,:,f) = frame';
    if f == framesPerVol
        v = v + 1;
        imgSum(:,:,v) = sum(vol(:,:,validPlanes), 3);
        f = 0;
        if mod(v, 1000) == 0
            fprintf('  read %d volumes so far (%.0f frames/s)\n', v, k/toc);
        end
    end
    p = nextIFD;
end
if v == 0
    error('overshoot_bar_stim_script:noFrames', 'No complete volumes could be read from %s -- not caching an empty stack.', tifPath);
end
imgSum = imgSum(:,:,1:v);
leftover = k - v * framesPerVol;
chainStr = {'ABNORMALLY', 'cleanly (next-IFD pointer = 0)'};
fprintf('  %d frames read -> %d complete volumes (%d leftover frame(s) dropped; json nframes_total=%d; IFD chain ended %s)\n', ...
    k, v, leftover, sp.nframes_total, chainStr{chainOk + 1});
sp.nVolumes = v;
sp.nFramesRead = k;
sp.tifChainOk = chainOk;
end

function clusterIdx = obs_pb_glomeruli(mask, nPerHemisphere)
% Skeletonize -> order (graph_sort) -> resample into 2*nPerHemisphere evenly
% spaced centroids -> nearest-centroid assignment. Same as process_im's
% clustering step (lpsp_cschrimson_reredo_script.m) and the *_pb_glomeruli
% copies in led_stim_berg4_script.m / overshoot_drug_stim_script.m. Assumes
% the mask is a single open arch. clusterIdx: 0 outside the mask, 1..2N
% inside, numbered along the skeleton from one end to the other.
[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('obs_pb_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
end
D    = bwdistgeodesic(skel, epIdx);
[~, idxFar] = max(D(:));
D2   = bwdistgeodesic(skel, idxFar);
mid  = (D + D2) == mode(D + D2, 'all');

[y_mid, x_mid] = find(mid);
[x_mid, y_mid] = graph_sort(x_mid, y_mid);

xq    = -min_axis:(length(x_mid) + min_axis);
x_mid = round(interp1(1:length(x_mid), x_mid, xq, 'linear', 'extrap'));
y_mid = round(interp1(1:length(y_mid), y_mid, xq, 'linear', 'extrap'));

keep  = ismember([x_mid', y_mid'], [x_mask, y_mask], 'rows');
x_mid = x_mid(keep);
y_mid = y_mid(keep);

nClusters = 2 * nPerHemisphere;
xq2 = linspace(1, length(y_mid), 2*nClusters + 1)';
centroids = [interp1(1:length(y_mid), y_mid, xq2), interp1(1:length(x_mid), x_mid, xq2)];
centroids = centroids(2:2:end-1, :);

[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);

clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end

function cmap = obs_redwhiteblue(n)
half = n / 2;
blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
whiteToRed  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
cmap = [blueToWhite; whiteToRed];
end
