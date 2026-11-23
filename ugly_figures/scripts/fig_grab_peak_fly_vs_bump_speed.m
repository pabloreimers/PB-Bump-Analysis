%% fig_grab_peak_fly_vs_bump_speed
% Q: does the GRAB(DA) peak scale with the FLY's rotational speed or with
% the BUMP's rotational speed? (They are highly correlated.)
%
% Data: EPG>GRAB(DA2m) trials from lpsp_cl_reproc_nobg_20261117.mat (21
% closed loop, 21 dark, 13 flies), PMT-offset-corrected per-bin dF/F
% (same preprocessing as fig_summary_lpsp_grab_bump_model reproc_nobg).
% Bump position/speed is taken from GRAB itself (PVA phase of each PB half
% and of both halves combined). EPG>GCaMP7f from the same file is only a
% benchmark for how well a strong EPG bump's velocity can be measured.
% (The jRGECO1a/GRAB co-imaging set was dropped: too noisy for this.)
%
% Using the same signal for bump and peak needs care: when DA is low the
% bump is harder to localize, which distorts measured bump speed exactly
% when the peak is small. So the within-trial test is also run CROSS-HALF
% (bump speed from one PB half, GRAB peak from the other, averaged both
% ways), and the across-trial test normalizes each trial's GRAB slope by
% its own GRAB range.
%
% Tests:
%  (1) within trials: non-overlapping 2/4/8 s windows; |fly rotation| and
%      |bump displacement| per window vs mean GRAB peak (PB max, lagged by
%      GRAB's measured 0.53 s); within-trial z; joint regression.
%  (2) across dark trials: without a cue the bump's gain relative to the
%      fly varies between trials. If GRAB follows bump speed, trials with
%      a higher bump gain should show a steeper GRAB-vs-fly-speed slope.
%      (Noise in bump displacement doesn't bias the gain slope.)
%
% Result when written:
%  - GRAB-derived bump velocity vs cue velocity |r| 0.23-0.50 (0.5-4 s
%    smoothing; EPG GCaMP7f benchmark 0.47-0.76); the two PB halves agree
%    on bump displacement only weakly (r 0.06-0.19). So per-window bump
%    speed is still noisy.
%  - Within trials: GRAB peak vs fly speed r 0.59-0.75; with both in the
%    model, fly 0.58-0.76 and bump -0.07..0.00 (same cross-half), in closed
%    loop and in the dark.
%  - Across dark trials: bump gain 0.01-0.75; vs GRAB-vs-fly slope
%    normalized to GRAB range: Spearman 0.24 (p ~0.3); the raw slope
%    correlates (0.53) but that is expected from weak-GRAB trials having
%    both a poorly measured bump (gain biased to 0) and small raw slopes.
%  - Lag scan (section 3): letting fly speed and GRAB bump speed each take
%    their own lag (-3..+3 s) vs the GRAB peak: fly speed peaks at +0.4 s
%    (CL, r 0.81) / +0.5 s (dark, 0.74); frame-wise bump speed never
%    exceeds r 0.04 at any lag; 4 s window displacement <= 0.09. In frames
%    where the predictor half's bump is clear (cross-half), bump speed alone
%    reaches 0.26 in closed loop at the SAME +0.4 s lag as fly speed, but
%    jointly fly 0.70 vs bump <= 0.07 at the best bump lag; dark 0.01.
%  - Folded, well-observed bumps (section 4): halves averaged after
%    aligning by each trial's left-right offset, bump phase/rho (as im.rho;
%    r 0.95 with it) from the folded profile, only frames with rho > 0.4
%    throughout the 1 s velocity window (~10% CL, ~5% dark). There, bump
%    speed tracks cue speed better (r 0.49 vs 0.29) and DOES relate to the
%    GRAB peak: CL fly r 0.63 (+0.4 s) vs bump 0.46 (+0.2 s), jointly fly
%    0.47 / bump 0.20; dark fly 0.59 vs bump 0.38, jointly 0.46 / 0.19.
%    Bump term positive in 12/13 CL and 18/21 dark trials (signrank p 0.01
%    / 0.001); fly term larger in 11/13 and 19/21; time-shifted control ~0.
%    No extra dark-specific dissociation: same split as closed loop.
%  - Matched pairs (section 5): within trials, non-overlapping 1.5 s
%    stretches with fly speed within 15% (> 1 rad/s), clear bump (rho > 0.4
%    in >= 80% of frames), bump speed >= 1.5x different: the faster-bump
%    stretch has the brighter GRAB peak in a per-trial median 63% (CL, 13
%    trials, p 0.002) / 66% (dark, 20 trials, p 0.0003) of pairs; median
%    GRAB difference +8% / +6%. Examples: CL fly 1.68 vs 1.69 rad/s, bump
%    0.11 vs 0.63, GRAB peak +34%; dark fly 1.45 vs 1.46, bump 0.18 vs 0.63,
%    GRAB +56% (but GRAB was already high before that window -> slower
%    state co-variation; direction of influence not determined).
%  - Wiggles / forward walking (section 6): could the slow-bump, high-|vr|
%    moments just be back-and-forth corrective turns or forward walking?
%    Over the 1.5 s windows turning is consistent in both (net / total ~0.94-
%    1.0) and net rotation similar, but slow-bump stretches have more forward
%    walking (~0.5 vs 0.3). Re-matching on NET rotation (consistent turning
%    >= 0.8 in both): faster-bump brighter in 65% (CL, p 0.03) / 68% (dark,
%    p 0.001); also matching forward speed: 66% (p 0.007) / 69% (p 0.01),
%    GRAB +6-8%. 4-predictor model (|vr|, net rotation, forward, bump): bump
%    0.21 (CL) / 0.17 (dark), > 0 in 12/13 and 18/21 trials (p 0.005 each);
%    forward speed slightly negative (-0.08 / -0.19).
%  - LPsP syt7f (section 7; 24 CL / 20 dark trials, ~1 fly per trial):
%    folded rho ~half of GRAB's, so threshold rho > 0.15. Joint model on the
%    PB max: fly 0.52 / bump 0.13 (CL, bump > 0 in 13/14 trials, p 0.01),
%    fly 0.62 / bump 0.08 (dark, 17/19, p 0.0007); same with the 4-predictor
%    controls and at rho > 0.10 / 0.20 (bump 0.07-0.13). Matched on |vr|:
%    faster bump higher in 60% of pairs (CL p 0.005, dark p 0.08); on net
%    rotation + forward speed few pairs (CL 66%, p 0.1; dark 53%).
%    The bump-speed term sits in the PB MEAN (global; LPsP 0.10-0.14, GRAB
%    0.26-0.34), not in max - mean (local bump height: ~0, -0.06..+0.06) --
%    for LPsP and for GRAB alike. BUT for GRAB that additive split is the
%    wrong frame: its trough is ~0 (trough/peak 0.12), so pure scaling moves
%    mean and max - mean together. Multiplicative split (7b), 4-predictor
%    model: fly speed scales GRAB without changing its shape (flatness
%    -0.14 / +0.03), while bump speed raises the gain (log peak +0.23 CL /
%    +0.16 dark, 12/13 and 18/21 trials) AND flattens/widens the profile
%    (mean/peak +0.34 / +0.22, width +0.32 / +0.15, trough +0.40 / +0.20;
%    13/13 and 20/21 trials), the flattening peaking later (+0.6-1.0 s).
%    A broader, raised profile after fast bump motion is what slow sensor
%    kinetics predict (the old position is still decaying while the new
%    one rises); smearing alone would LOWER the peak, so the gain increase
%    is, if anything, underestimated. LPsP: gain +0.12 / +0.08, shape
%    +0.08..+0.16 (but |vr| also flattens LPsP via its global component).
%  - Summary figure (section 8; fly + bump joint model, per-trial medians):
%    peak: GRAB fly 0.52 / 0.45, bump 0.16 / 0.19 (CL / dark); LPsP fly
%    0.52 / 0.61, bump 0.11 / 0.08. Trough: GRAB fly 0.15 / 0.42, bump
%    0.36 / 0.17; LPsP (global readout) fly 0.39 / 0.56, bump 0.12 / 0.07
%    (13/14 and 18/19 trials > 0). Unchanged with net rotation + forward
%    speed added. -> fig_grab_peak_fly_vs_bump_speed_peak_trough.png
%    EPG GCaMP7f (23 trials): peak fly 0.28 / 0.43, bump 0.23 / 0.13;
%    trough fly 0.14 / 0.20, bump 0.33 / 0.21 (10/11-11/11 trials > 0). Its
%    bump term is slow: ~0 at lag <= 0, rising to a broad plateau at +1 to
%    +2 s (so the +1.0 s window edge isn't hiding a larger effect).
%  - Widening (section 9): far field = folded bins > 90 deg (or 135) from
%    the bump throughout the last 2 s (or 4 s), untouched by bump or trail.
%    LPsP far field still scales with bump speed beyond the ball (4-pred
%    model: bump +0.14 / +0.12 CL / dark; 135 deg +0.15 / +0.09; 4 s +0.15
%    / +0.13), also cross-half (+0.09 / +0.12), with bump concentration as
%    covariate (+0.11 / +0.13) and as net 1.5 s bump displacement (+0.12 /
%    +0.08) -> not trail, not phase jitter. Given that readout, bump speed
%    adds NOTHING to the LPsP peak (-0.02 / -0.01; cross-half the same) and
%    is slightly negative for peak - far (-0.06 / -0.02). So in LPsP the
%    bump-speed effect on brightness is carried entirely by the global
%    (speed) signal, which tracks bump rotation beyond ball rotation.
%    Caveat: EPG7f far field also rises with bump speed (+0.2-0.35 z) but
%    in raw units it is small (0.04-0.05 dF/F per rad/s vs LPsP 0.11-0.18),
%    so inherited EPG input explains at most a minor part. GRAB: far field
%    rises too (+0.28-0.48, partly its multiplicative gain); with it in the
%    model the GRAB bump term drops to -0.02..+0.12 (inconsistent). EPG7f
%    peak keeps its bump term (+0.16-0.21) with far in the model.
%  - Strict far-field trough (section 10; min over bins > 135 deg from the
%    bump for the last 4 s; same frames in all rows): LPsP trough ~ fly
%    0.44 / 0.53, bump +0.15 / +0.14 (13/14, 16/16 trials); peak ~ fly 0.50
%    / 0.56, bump +0.11 / +0.12; peak with trough in the model: trough 0.89
%    / 0.93, fly 0.07 / 0.04, bump -0.03 / -0.01 (1/14, 5/15 > 0). EPG7f
%    keeps a bump term (+0.09 / +0.10) given its trough; GRAB keeps only
%    ~9% of clear frames (4-7 trials): +0.11 / -0.01, inconclusive.
%    -> fig_grab_peak_fly_vs_bump_speed_strict_trough.png
%  - Per-frame anti-phase trough (section 11): trough > 135 deg from the
%    CURRENT bump, dropping only bins the trail came within 90 deg of over
%    the last 2 s (atr2): kept in 54% (GRAB), 95-100% (LPsP), 89% (EPG7f)
%    of clear frames vs 6% / 49-74% / 37-45% for section 10's version.
%    LPsP unchanged: trough bump +0.14 / +0.11; peak given trough, bump
%    -0.03 / -0.00. GRAB trough bump +0.32 / +0.19 (12/12, 13/17 trials);
%    peak given its trough, bump +0.06 / +0.07 (8/12, 10/17; n.s.). EPG7f
%    keeps a bump term (+0.15 / +0.21) given its trough.
%    Trail window from the measured decay (11b) instead of 2/4 s: GRAB
%    tau ~1.1 s, EPG7f ~2.2 s, LPsP ~2.6 s (poor fit; upper bounds). With a
%    3-tau trail GRAB keeps 30-35% of clear frames (section 10: 6%); trough
%    bump +0.35 / +0.17; peak given trough, bump +0.07 / +0.08 (6/11, 9/13
%    trials > 0). LPsP and EPG7f unchanged across 1-3 tau.
%  - Bump position from "commanded" speed (section 12): ball velocity
%    scaled by the LPsP trough's deviation from ball speed gets a small,
%    consistent positive weight in closed loop (+0.07..+0.14, 10-12 of
%    11-14 trials; shuffled ~0) but does NOT improve out-of-sample
%    prediction of LPsP bump displacement (CV R2 0.22-0.27 ball vs 0.18-
%    0.25); trough-set speed alone is worse (0.14-0.22). Dark: few clear
%    windows, no gain. Limits: LPsP bump phase is noisy (ball explains only
%    ~25%), and the closed-loop bump is held to a ball-driven cue.
%    Dark, more frames (12b: clarity only at window endpoints, rho > 0.10,
%    overlapping windows, block CV): 19 trials, ~100-220 independent
%    windows, but ball R2 only 0.05-0.09 and no gain from the trough (dev
%    coef ~0, 10/19). Reason: the dark LPsP bump trajectory is mostly
%    noise -- the two PB halves agree on 2 s displacement at r = 0.21
%    (CL 0.39; GRAB 0.84-0.95, EPG7f 0.98). Needs a strong bump in the
%    same fly (co-imaging) to test. [Overstated -- see 12c: much of the poor
%    agreement was 2 pi slips in the UNWRAPPED phase, and the GRAB numbers
%    were clear frames only.] With WRAPPED endpoint displacement (12c):
%    dark still no gain (CV R2 ball 0.01-0.11, +dev the same; dev coef
%    -0.00..+0.10, not above its shuffled control +0.00..+0.09); closed
%    loop dev coef +0.06..+0.10 (12-13/14 trials, shuffled ~0), CV R2 gain
%    <= +0.02.
%  - Cross-centred profiles by fly speed (section 13; all frames, each half
%    centred on the other half's bump): LPsP has a clear bump in the dark,
%    0.12-0.22 dF/F peak - trough (CL 0.18-0.29), as large in dF/F as
%    GRAB's (0.04-0.26 dark), but on a high, speed-driven baseline and with
%    ~2x GRAB's frame noise (0.056 vs 0.028 dF/F per bin): bump height /
%    frame noise 2.2-3.9 in the dark (GRAB 1.5-9.5, EPG7f 16-35). LPsP's
%    noisy dark bump position = low per-frame SNR, not a missing bump.
%    LPsP bump height barely grows with speed (its baseline does); GRAB's
%    grows ~6x, EPG7f's ~2x (dark).
%    -> fig_grab_peak_fly_vs_bump_speed_profiles_by_speed.png
%  - GRAB as "motion energy" of LPsP calcium (section 14)? Oriented energy
%    can't see LPsP's uniform speed signal and is centred like the product
%    (fly speed dominating argues against it). Temporal-change energy
%    ([dL/dt]+) predicts release at the bump's leading edge and turn-onset
%    transients. Closed loop, profiles aligned on the cue-predicted bump
%    (lag from EPG7f, -0.44 s): GRAB centroid -1, -6, -10, -18, -29 deg (0.12
%    -3 rad/s; LPsP similar, EPG7f ~0) = a centred bump delayed by sensor
%    kinetics (product model, tau 0.3-0.5 s: -2..-30 deg); [dL/dt]+ predicts
%    +8..+80 deg ahead at any tau <= 1.1 s. Turn onsets: GRAB residual vs a
%    low-passed |vr| model +0.08 SD (CL, 19/21) / +0.04 (dark, n.s.); no
%    fast offset drop; LPsP trough no overshoot (dark -0.15 SD).
%    -> fig_grab_peak_fly_vs_bump_speed_motion_energy.png
%  - Depletion test (section 15, closed loop, bump position from the cue):
%    recent occupancy of the bump's current location lowers the GRAB peak
%    at matched speed: basic model -0.22 (tau 1 s) .. -0.14 (8 s), 13/13
%    trials; with cubic speed terms, recent-speed history, reversal
%    covariates and the first 30 s excluded -0.16 / -0.08 / -0.05 / -0.05
%    (10-12/13, p 0.0005-0.02). Location-specific: occupancy of the
%    OPPOSITE location ~0 (6-8/13). Recovery curve (moving bump): revisit
%    within 1-2 s -0.11 SD, 2-30 s ~0, > 30 s +0.13 SD; fresh (> 8 s) minus
%    revisit (< 2 s) +0.12 SD (9/10, p 0.014). Effect ~0.1 SD; consistent
%    with fast-recovering depletion (~1-2 s) on top of the speed gain.
%    CORRECTION (section 17): the strict controls lacked fly-speed history;
%    a product model with no depletion already shows -0.09 (tau 1 s) from
%    its ~1 s integration of fly speed. With mean |vr| over the past 1/2/4 s
%    (and squares) added, the data effect is -0.07 .. -0.09 (tau 1 s, 12/13),
%    -0.03 .. -0.04 (2 s), ~0 beyond; data minus product residual -0.11 /
%    -0.06 / -0.05 (tau 1/2/4 s, p 0.01-0.04). Recovery curve: revisit < 2 s
%    -0.18 SD, fresh minus revisit +0.24 SD (9/10, p 0.01).
%    -> fig_grab_peak_fly_vs_bump_speed_depletion.png
%  - Release models from LPsP (section 16; synthetic LPsP built per GRAB
%    fly from an LPsP input model fitted in LPsP flies: F/F0 = 0.99 + 0.73
%    |vr| + 0.20 |bump v| + 0.17 bump(theta), tau 1 s, R2 0.57). CV R2 of
%    the full GRAB(theta,t), closed loop: product (ball gain) 0.53; product
%    with LPsP gain incl. bump speed 0.53 (-0.006); per-glomerulus power
%    Lhat^n 0.39 and depletion Lhat^n R 0.39 (-0.11, worse in 13/13; depl.
%    vs power +0.003, 11/13). Reason: GRAB's bump / trough contrast RISES
%    with speed (0.74 -> 1.7-2.1), LPsP's FALLS (0.44 -> 0.16), and any per-
%    glomerulus transform of the summed signal inherits the fall (power /
%    depletion 0.29 -> 0.09); only the global x local product rises (0.83 ->
%    1.84). Inside the product, a free bump-speed weight (w ~ -0.2, CV
%    -0.02) or depletion of the product's drive (u -> 0 in one half, CV
%    -0.004) don't help. In closed loop bump and fly speed are nearly
%    collinear, and the section-15 effect (~0.1 SD) is too small to show in
%    a full-field R2.
%    -> fig_grab_peak_fly_vs_bump_speed_release_models.png
%  - EPG gate x adapting LPsP (section 17): Rel = [dL - beta LP_tauA(dL)]+
%    (1 + k2 EPG). Free fit: beta ~0.5-0.75, tauA ~0.3-0.4 s, CV R2 0.542 vs
%    product 0.532 (+0.010, 10/13, p 0.13); beats PB-wide-only adaptation
%    (+0.009, 12/13, p 0.0005), which is no better than the product. Pure
%    temporal change (beta = 1) 0.40. But the fitted model does not capture
%    the occupancy effect (prediction -0.04 vs data -0.09 at tau 1 s; the
%    residual keeps -0.12), nor does adaptation held at 1-2 s (worse fit).
%    -> fig_grab_peak_fly_vs_bump_speed_gated_adapting.png
%  => With well-observed bumps, bump rotation adds a modest independent
%     contribution, but fly rotation remains the stronger predictor.
%  (Earlier, all-frame summary:) No sign that GRAB follows bump speed beyond fly speed; GRAB tracks
%     the fly's rotation. Limit: residual bump motion (beyond fly turning)
%     is mostly measurement noise here, so only a large bump-speed effect
%     is ruled out. Decisive test: closed-loop gain manipulations.

%% load + per-trial signals
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_reproc_nobg_20261117.mat')), data_dir = '.data'; else, data_dir = 'data'; end
B = load(fullfile(data_dir,'lpsp_cl_reproc_nobg_20261117.mat'),'all_data'); B = B.all_data(:);
ind = arrayfun(@(s) s.info.indicator, B, 'UniformOutput', false);
gi = find(strcmp(ind,'EPG>GRAB(DA2m)')); gi = gi(:)';
snr_frac = 0.5; a16 = linspace(-pi, pi-2*pi/16, 16)';
S = arrayfun(@(k) grab_trial_signals(B(k), snr_frac, a16), gi);
dk = [S.dark];
lagG = 0.53;                     % GRAB peak lag behind |fly speed| (median of per-trial xcorr peaks)

%% (0) how well is bump velocity measured? GRAB (closed loop) vs EPG GCaMP7f benchmark
wins = [0.5 1 2 4];
lagbest = @(a,b,L) max(arrayfun(@(t) abs(corr(a(max(1,1+t):min(end,end+t))', b(max(1,1-t):min(end,end-t))')), -L:L));
Rk = nan(1,numel(S)); Vc = nan(numel(S),numel(wins)); Rel = nan(numel(S),numel(wins));
for j = 1:numel(S)
    t = S(j); fr = t.fr; n = numel(t.pA);
    if ~t.dark
        ph = angle(exp(1i*t.pA)); cue = angle(exp(1i*t.cue)); best = 0;
        for tt = -round(2.5*fr):round(1*fr)
            idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n);
            for sg = [-1 1], best = max(best, abs(mean(exp(1i*(ph(idx) + sg*cue(idx+tt)))))); end
        end
        Rk(j) = best;
    end
    for w = 1:numel(wins)
        L = max(2,round(wins(w)*fr)); e = 1:L:n-L; e = e(arrayfun(@(a) ~any(t.bad(a:a+L-1)), e));
        Rel(j,w) = corr((t.pL(e+L-1)-t.pL(e))', (t.pR(e+L-1)-t.pR(e))');
        if ~t.dark
            sm = @(x) smoothdata(x,'gaussian',max(3,round(wins(w)*fr)));
            Vc(j,w) = lagbest(gradient(sm(t.pA))*fr, gradient(sm(t.cue))*fr, round(0.8*fr));
        end
    end
end
clT = ~dk & Rk >= 0.4;           % closed-loop trials whose GRAB bump tracks the cue
selE = find(strcmp(ind,'EPG>GCaMP7f') & ~arrayfun(@(s) s.info.dark, B)); selE = selE(:)';
Vep = nan(numel(selE),numel(wins));
for j = 1:numel(selE)
    s = B(selE(j)); xb = s.ft.xb(:)'; xf = s.ft.xf(:); n = numel(xb); fr = s.info.fps;
    cue = unwrap(angle(interp1(xf, exp(1i*s.ft.cue(:)), xb, 'linear', 'extrap')));
    ph = angle(sum(exp(1i*a16)*ones(1,n) .* zscore(s.im.f(1:16,:),0,2), 1));
    for w = 1:numel(wins)
        sm = @(x) smoothdata(x,'gaussian',max(3,round(wins(w)*fr)));
        Vep(j,w) = lagbest(gradient(sm(unwrap(ph)))*fr, gradient(sm(cue))*fr, round(0.6*fr));
    end
end
fprintf('bump velocity vs cue velocity |r| by smoothing %s s: GRAB (n=%d tracking CL) %s | EPG GCaMP7f %s\n', mat2str(wins), sum(clT), mat2str(median(Vc(clT,:),1),2), mat2str(median(Vep,1),2));
fprintf('half-1 vs half-2 bump displacement agreement r: CL %s | dark %s\n', mat2str(median(Rel(clT,:),1),2), mat2str(median(Rel(dk,:),1),2));

%% (1) within trials: windows, joint regression (and cross-half)
winS = [2 4 8]; JW = nan(2,numel(winS),6);   % cond x window x [r(fly,bump) rFly rBump bFly bBump bBumpXhalf]
for d_ = 0:1
    for wi = 1:numel(winS)
        ws = winS(wi); W = [];
        for j = find(dk == d_ & (d_ == 1 | Rk >= 0.4))
            t = S(j); fr = t.fr; n = numel(t.pA); L = round(ws*fr); lag = round(lagG*fr);
            e = 1:L:n-L-lag; e = e(arrayfun(@(a) ~any(t.bad(a:a+L-1+lag)), e)); if numel(e) < 10, continue, end
            fly = abs(t.head(e+L-1)-t.head(e))'/ws; bmp = abs(t.pA(e+L-1)-t.pA(e))'/ws;
            bL = abs(t.pL(e+L-1)-t.pL(e))'/ws; bR = abs(t.pR(e+L-1)-t.pR(e))'/ws;
            gp = arrayfun(@(a) mean(t.pk(a+lag:a+L-1+lag),'omitnan'), e)';
            gL = arrayfun(@(a) mean(t.pkL(a+lag:a+L-1+lag),'omitnan'), e)'; gR = arrayfun(@(a) mean(t.pkR(a+lag:a+L-1+lag),'omitnan'), e)';
            R_ = [fly bmp bL bR gp gL gR]; W = [W; (R_ - mean(R_,1,'omitnan'))./std(R_,0,1,'omitnan')]; %#ok<AGROW>
        end
        W = W(all(~isnan(W),2),:); X1 = ones(size(W,1),1);
        b  = [X1 W(:,1) W(:,2)] \ W(:,5);
        bx = ([X1 W(:,1) W(:,3)] \ W(:,7) + [X1 W(:,1) W(:,4)] \ W(:,6))/2;   % bump from one half, peak from the other
        JW(d_+1,wi,:) = [corr(W(:,1),W(:,2)) corr(W(:,1),W(:,5)) corr(W(:,2),W(:,5)) b(2) b(3) bx(3)];
        fprintf('%-11s %d s (n=%4d): r(fly,bump) %.2f | GRAB vs fly %.2f, vs bump %.2f | joint fly %.2f, bump %.2f | cross-half bump %.2f\n', ...
            ifelse_c(d_), ws, size(W,1), squeeze(JW(d_+1,wi,:)));
    end
end

%% (2) across dark trials: does GRAB's speed tuning scale with the bump's gain?
TR = nan(numel(S),3);
for j = find(dk)
    t = S(j); fr = t.fr; n = numel(t.pA); L = round(4*fr);
    e = 1:L:n-L; e = e(arrayfun(@(a) ~any(t.bad(a:a+L-1)), e));
    dB = (t.pA(e+L-1)-t.pA(e))'; dH = (t.head(e+L-1)-t.head(e))';
    p = [dH ones(size(dH))] \ dB;                         % bump gain = slope of bump displacement on fly rotation
    lag = round(lagG*fr); a = abs(t.vr(1:n-lag)); y = t.pk(1+lag:n); m = ~t.bad(1+lag:n) & ~isnan(y);
    q = polyfit(a(m), y(m), 1); rngG = prctile(y(m),95) - prctile(y(m),5);
    TR(j,:) = [abs(p(1)) q(1)/rngG q(1)];
end
[rN,pN] = corr(TR(dk,1), TR(dk,2), 'type','Spearman'); [rR,pR] = corr(TR(dk,1), TR(dk,3), 'type','Spearman');
fprintf('dark trials (n=%d): bump gain %.2f-%.2f; Spearman(bump gain, GRAB-vs-fly slope / GRAB range) = %.2f (p=%.2g); raw slope %.2f (p=%.2g, confounded by GRAB quality)\n', ...
    sum(dk), min(TR(dk,1)), max(TR(dk,1)), rN, pN, rR, pR);

%% (3) lag scan: each predictor gets its own lag relative to the GRAB peak
% Frame-wise: fly |vr| (0.5 s smoothing) and GRAB bump speed (|d phase/dt|,
% 1 s smoothing) shifted by -3..+3 s (positive = predictor leads). Also the
% cross-half version restricted to frames where the PREDICTOR half's own
% bump is clear (top 50% PVA strength), with a joint model (fly at its best
% lag, bump at each lag) to see if any bump lag adds to fly speed.
lagsS = -3:0.1:3; LS = struct();
for d_ = 0:1
    trials = find(dk == d_ & (d_ == 1 | Rk >= 0.4));
    CF = nan(numel(trials),numel(lagsS)); CB = CF; CQ = CF; JB = CF;
    for jj = 1:numel(trials)
        t = S(trials(jj)); fr = t.fr; n = numel(t.pA); sm = @(x) smoothdata(x,'gaussian',max(3,round(1*fr)));
        fly = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr)));
        bA = abs(gradient(sm(t.pA))*fr); bL = abs(gradient(sm(t.pL))*fr); bR = abs(gradient(sm(t.pR))*fr);
        qL = t.strL > median(t.strL); qR = t.strR > median(t.strR); Lf = round(0.45*fr);
        for li = 1:numel(lagsS)
            L = round(lagsS(li)*fr); i1 = max(1+max(L,Lf),1):min(n,n+min(L,Lf)); i0 = i1 - L; iF = i1 - Lf;
            base = ~t.bad(i1) & ~t.bad(i0) & ~isnan(t.pk(i1)) & ~isnan(t.pkL(i1)) & ~isnan(t.pkR(i1));
            CF(jj,li) = corr(fly(i0(base))', t.pk(i1(base))'); CB(jj,li) = corr(bA(i0(base))', t.pk(i1(base))');
            o = nan(2,2);
            for hh = 1:2
                if hh == 1, m = base & qL(i0); x = bL(i0(m)); y = t.pkR(i1(m)); else, m = base & qR(i0); x = bR(i0(m)); y = t.pkL(i1(m)); end
                Zq = zscore([reshape(fly(iF(m)),[],1) x(:) y(:)]); b = [ones(size(Zq,1),1) Zq(:,1:2)] \ Zq(:,3);
                o(hh,:) = [corr(Zq(:,2),Zq(:,3)) b(3)];
            end
            CQ(jj,li) = mean(o(:,1)); JB(jj,li) = mean(o(:,2));
        end
    end
    LS(d_+1).CF = CF; LS(d_+1).CB = CB; LS(d_+1).CQ = CQ; LS(d_+1).JB = JB;
    [mf,jf] = max(median(CF,1)); [mb,jb] = max(median(CB,1)); [mq,jq] = max(median(CQ,1)); [mj,jm] = max(median(JB,1));
    fprintf('%-11s lag scan: fly best %+.1f s (r %.2f) | bump (all frames) %+.1f s (r %.2f) | bump, clear-half cross-half %+.1f s (r %.2f) | joint bump coefficient max %.2f at %+.1f s\n', ...
        ifelse_c(d_), lagsS(jf), mf, lagsS(jb), mb, lagsS(jq), mq, mj, lagsS(jm));
end

%% (4) folded halves, well-observed bumps only (rho > 0.4)
% Halves averaged after aligning by the trial's left-right phase offset
% (fractional circular shift) -> 16-bin profile with more signal; bump
% phase and rho (= |mean(z e^{i alpha})|, as im.rho) from it. Bump speed
% (1 s smoothing) is used only where rho > rthr throughout that window.
% NB: rho comes from GRAB too, so these frames are also higher-DA moments;
% that narrows the GRAB range and lowers BOTH correlations, but the fly-vs-
% bump comparison within the same frames stays fair.
rthr = 0.4; Fd = arrayfun(@(k) grab_folded(B(k), snr_frac, a16), gi);
RkF = nan(1,numel(Fd));
for j = find(~[Fd.dark])
    t = Fd(j); n = numel(t.ph); ph = angle(exp(1i*t.ph)); cue = angle(exp(1i*t.cue)); best = 0;
    for tt = -round(2.5*t.fr):round(1*t.fr)
        idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n);
        for sg = [-1 1], best = max(best, abs(mean(exp(1i*(ph(idx) + sg*cue(idx+tt)))))); end
    end
    RkF(j) = best;
end
FO = struct();
for d_ = 0:1
    trials = find([Fd.dark] == d_ & (d_ == 1 | RkF >= 0.4));
    CF = nan(numel(trials),numel(lagsS)); CB = CF; JF = CF; JB = CF; NF = CF;
    for jj = 1:numel(trials)
        t = Fd(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(1*fr)); Lf = round(0.45*fr);
        fly = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr)));
        bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
        good = movmin(t.rho_raw,w) > rthr & ~movmax(t.bad,w);
        for li = 1:numel(lagsS)
            L = round(lagsS(li)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
            m = good(i0) & ~t.bad(i1) & ~isnan(t.pk(i1)); if sum(m) < 50, continue, end
            CF(jj,li) = corr(fly(i0(m))', t.pk(i1(m))'); CB(jj,li) = corr(bs(i0(m))', t.pk(i1(m))'); NF(jj,li) = sum(m);
            Zq = zscore([fly(iF(m))' bs(i0(m))' t.pk(i1(m))']); b = [ones(sum(m),1) Zq(:,1:2)] \ Zq(:,3);
            JF(jj,li) = b(2); JB(jj,li) = b(3);
        end
    end
    wm = @(X) sum(X.*NF,1,'omitnan')./sum(NF.*~isnan(X),1,'omitnan');    % frame-weighted mean across trials
    FO(d_+1).CF = wm(CF); FO(d_+1).CB = wm(CB); FO(d_+1).JF = wm(JF); FO(d_+1).JB = wm(JB);
    [mf,jf] = max(FO(d_+1).CF); [mb,jb] = max(FO(d_+1).CB); [mj,jm] = max(FO(d_+1).JB);
    jbT = JB(:,jm); jfT = JF(:,jm); ok = ~isnan(jbT);
    fprintf('%-11s folded, rho > %.1f (%d trials): fly best %+.1f s (r %.2f) | bump best %+.1f s (r %.2f) | joint at %+.1f s: fly %.2f, bump %.2f | bump > 0 in %d/%d trials (signrank p=%.2g), fly > bump in %d/%d\n', ...
        ifelse_c(d_), rthr, numel(trials), lagsS(jf), mf, lagsS(jb), mb, lagsS(jm), FO(d_+1).JF(jm), mj, sum(jbT(ok)>0), sum(ok), signrank(jbT(ok)), sum(jfT(ok)>jbT(ok)), sum(ok));
end

%% (5) matched pairs: same fly speed, different bump speed
% Within each trial: non-overlapping 1.5 s GRAB windows; fly speed (lag
% 0.45 s) and folded-bump speed (lag 0.2 s) averaged over the matching
% earlier stretches; clear bump (rho > 0.4) in >= 80% of the bump frames.
% Pairs: fly speed within 15% (both > 1 rad/s), bump speed >= 1.5x apart.
% Population: per trial, fraction of pairs where the faster-bump stretch
% has the brighter mean GRAB peak. Examples: typical trials (60-80%), fly
% within 5%, bump >= 2.5x, GRAB difference near that condition's 75th pct.
wsP = 1.5; Pp = [];
for j = 1:numel(Fd)
    t = Fd(j); if ~t.dark && ~(RkF(j) >= 0.4), continue, end
    fr = t.fr; n = numel(t.ph); L = round(wsP*fr); st = round(0.5*fr); Lf = round(0.45*fr); Lb = round(0.2*fr);
    fly = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',max(3,round(fr))))*fr);
    good = t.rho_raw > rthr & ~t.bad; rows = [];
    for a = 1+max(Lf,Lb):st:n-L
        g_ = a:a+L-1; fi = g_ - Lf; bi = g_ - Lb;
        if mean(good(bi)) < 0.8 || any(t.bad(g_)) || any(isnan(t.pk(g_))), continue, end
        rows(end+1,:) = [a mean(fly(fi)) mean(bs(bi(good(bi)))) mean(t.pk(g_))]; %#ok<AGROW>
    end
    for p1 = 1:size(rows,1)-1
        for p2 = p1+1:size(rows,1)
            if abs(rows(p1,1)-rows(p2,1)) < L, continue, end
            fA = rows(p1,2); fB = rows(p2,2); if min(fA,fB) < 1 || abs(fA-fB)/mean([fA fB]) > 0.15, continue, end
            if max(rows(p1,3),rows(p2,3))/max(min(rows(p1,3),rows(p2,3)),eps) < 1.5, continue, end
            if rows(p1,3) > rows(p2,3), q = [p2 p1]; else, q = [p1 p2]; end
            Pp(end+1,:) = [j t.dark rows(q(1),1) rows(q(2),1) rows(q(1),2) rows(q(2),2) rows(q(1),3) rows(q(2),3) rows(q(1),4) rows(q(2),4)]; %#ok<AGROW>
        end
    end
end
EXP = []; PTs = cell(1,2);
for d_ = 0:1
    sel = Pp(:,2) == d_; tr = unique(Pp(sel,1));
    pt = arrayfun(@(k) mean(Pp(sel & Pp(:,1)==k,10) > Pp(sel & Pp(:,1)==k,9)), tr); npt = arrayfun(@(k) sum(sel & Pp(:,1)==k), tr);
    ptr = arrayfun(@(k) median((Pp(sel & Pp(:,1)==k,10) - Pp(sel & Pp(:,1)==k,9))./Pp(sel & Pp(:,1)==k,9)), tr);
    keep = npt >= 5; PTs{d_+1} = [tr pt npt];
    fprintf('%-11s matched pairs: %d in %d trials; faster-bump stretch brighter in per-trial median %.0f%% (signrank vs 50%% p=%.2g, %d trials); per-trial median GRAB difference %+.0f%%\n', ...
        ifelse_c(d_), sum(sel), numel(tr), 100*median(pt(keep)), signrank(pt(keep)-0.5), sum(keep), 100*median(ptr(keep)));
    okT = npt >= 20 & pt >= 0.6 & pt <= 0.8; cs = find(sel & ismember(Pp(:,1), tr(okT)));
    gd = (Pp(cs,10)-Pp(cs,9))./Pp(cs,9);
    cand = cs(abs(Pp(cs,5)-Pp(cs,6))./mean(Pp(cs,5:6),2) <= 0.05 & Pp(cs,8)./Pp(cs,7) >= 2.5 & gd > 0);
    if ~isempty(cand)
        g = (Pp(cand,10)-Pp(cand,9))./Pp(cand,9); [~,o] = sort(abs(g - prctile(gd(gd>0),75))); EXP(end+1,:) = Pp(cand(o(1)),:); %#ok<AGROW>
    end
end

%% (6) are slow-bump moments just wiggles or forward walking?
% (a) re-match pairs on NET rotation |mean(vr)| over the fly-lag-shifted
%     window (both stretches turning consistently: net/total >= 0.8), then
%     additionally on forward speed (within 0.15); (b) frame-wise joint model
%     with |vr|, net rotation (|vr smoothed 1 s, signed|), forward speed and
%     bump speed.
for j = 1:numel(Fd)
    s = B(gi(j)); Fd(j).fw = interp1(s.ft.xf(:), smoothdata(s.ft.f_speed(:),'gaussian',30), s.ft.xb(:)', 'linear', 'extrap');
end
mlab = {'net rotation','net rotation + forward speed'};
for mm = 1:2
    P2 = [];
    for j = 1:numel(Fd)
        t = Fd(j); if ~t.dark && ~(RkF(j) >= 0.4), continue, end
        fr = t.fr; n = numel(t.ph); L = round(wsP*fr); st = round(0.5*fr); Lf = round(0.45*fr); Lb = round(0.2*fr);
        bs = abs(gradient(smoothdata(t.ph,'gaussian',max(3,round(fr))))*fr); good = t.rho_raw > rthr & ~t.bad; rows = [];
        for a = 1+max(Lf,Lb):st:n-L
            g_ = a:a+L-1; fi = g_ - Lf; bi = g_ - Lb;
            if mean(good(bi)) < 0.8 || any(t.bad(g_)) || any(isnan(t.pk(g_))), continue, end
            rows(end+1,:) = [a abs(mean(t.vr(fi))) mean(t.fw(fi)) mean(bs(bi(good(bi)))) mean(t.pk(g_)) abs(mean(t.vr(fi)))/mean(abs(t.vr(fi)))]; %#ok<AGROW>
        end
        for p1 = 1:size(rows,1)-1
            for p2 = p1+1:size(rows,1)
                if abs(rows(p1,1)-rows(p2,1)) < L, continue, end
                nA = rows(p1,2); nB = rows(p2,2); if min(nA,nB) < 1 || abs(nA-nB)/mean([nA nB]) > 0.15, continue, end
                if min(rows([p1 p2],6)) < 0.8, continue, end
                if mm == 2 && abs(rows(p1,3)-rows(p2,3)) > 0.15, continue, end
                if max(rows(p1,4),rows(p2,4))/max(min(rows(p1,4),rows(p2,4)),eps) < 1.5, continue, end
                if rows(p1,4) > rows(p2,4), q = [p2 p1]; else, q = [p1 p2]; end
                P2(end+1,:) = [j t.dark rows(q(1),5) rows(q(2),5)]; %#ok<AGROW>
            end
        end
    end
    for d_ = 0:1
        sel = P2(:,2) == d_; tr = unique(P2(sel,1)); npt = arrayfun(@(k) sum(sel & P2(:,1)==k), tr); keep = npt >= 5;
        pt = arrayfun(@(k) mean(P2(sel & P2(:,1)==k,4) > P2(sel & P2(:,1)==k,3)), tr);
        ptr = arrayfun(@(k) median((P2(sel & P2(:,1)==k,4) - P2(sel & P2(:,1)==k,3))./P2(sel & P2(:,1)==k,3)), tr);
        if sum(keep) < 3, continue, end
        fprintf('matched on %-29s %-11s %5d pairs, %2d trials: faster-bump brighter in per-trial median %.0f%% (p=%.2g), GRAB diff %+.0f%%\n', ...
            mlab{mm}, ifelse_c(d_), sum(sel), sum(keep), 100*median(pt(keep)), signrank(pt(keep)-0.5), 100*median(ptr(keep)));
    end
end
for d_ = 0:1
    trials = find([Fd.dark] == d_ & (d_ == 1 | RkF >= 0.4)); Bc = nan(numel(trials),numel(lagsS),4); NF = nan(numel(trials),numel(lagsS));
    for jj = 1:numel(trials)
        t = Fd(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(0.45*fr);
        sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); net = abs(smoothdata(t.vr,'gaussian',max(3,round(fr))));
        fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
        good = movmin(t.rho_raw,w) > rthr & ~movmax(t.bad,w);
        for li = 1:numel(lagsS)
            L = round(lagsS(li)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
            m = good(i0) & ~t.bad(i1) & ~isnan(t.pk(i1)); if sum(m) < 50, continue, end
            Zq = zscore([sp(iF(m))' net(iF(m))' fw(iF(m))' bs(i0(m))' t.pk(i1(m))']);
            b = [ones(sum(m),1) Zq(:,1:4)] \ Zq(:,5); Bc(jj,li,:) = b(2:5); NF(jj,li) = sum(m);
        end
    end
    wm = @(X) sum(X.*NF,1,'omitnan')./sum(NF.*~isnan(X),1,'omitnan');
    [mb,jb] = max(wm(Bc(:,:,4))); co = arrayfun(@(k) subsref(wm(Bc(:,:,k)), struct('type','()','subs',{{jb}})), 1:4);
    pt = Bc(:,jb,4); ok = ~isnan(pt);
    fprintf('%-11s 4-predictor model at bump lag %+.1f s: |vr| %.2f, net rotation %.2f, forward %.2f, bump %.2f | bump > 0 in %d/%d trials (signrank p=%.2g)\n', ...
        ifelse_c(d_), lagsS(jb), co, sum(pt(ok) > 0), sum(ok), signrank(pt(ok)));
end

%% (7) same question for LPsP syt7f (presynaptic Ca): does bump speed add to fly speed?
% LPsP>syt7f trials from the same file (24 CL, 20 dark). Folded bump as in
% (4), but the LPsP bump is weak on top of the speed-driven global signal:
% folded rho is ~half of GRAB's (median 0.11 vs 0.22, almost never > 0.4).
% So the clarity threshold is rho > 0.15 (bump-vs-cue velocity |r| 0.47 in
% those frames, ~30% CL / ~13% dark frames); 0.10 and 0.20 for robustness.
% Fly-speed lag from the LPsP lag scan (+0.4 s CL, +0.3 s dark). Outcomes:
% PB max (as for GRAB), and its split into PB mean (global) and
% max - mean (local bump height); GRAB split the same way for comparison.
li = find(strcmp(ind,'LPsP>syt7f')); li = li(:)'; lfL = [0.4 0.3];
Fl = arrayfun(@(k) grab_folded(B(k), snr_frac, a16), li);
for j = 1:numel(Fl)
    s = B(li(j)); Tg = grab_trial_signals(s, snr_frac, a16);
    Fl(j).fw = interp1(s.ft.xf(:), smoothdata(s.ft.f_speed(:),'gaussian',30), s.ft.xb(:)', 'linear', 'extrap');
    Fl(j).mn = (Tg.mnL + Tg.mnR)/2; Fl(j).loc = Fl(j).pk - Fl(j).mn;
end
for j = 1:numel(Fd), Tg = grab_trial_signals(B(gi(j)), snr_frac, a16); Fd(j).mn = (Tg.mnL + Tg.mnR)/2; Fd(j).loc = Fd(j).pk - Fd(j).mn; end
RkL = nan(1,numel(Fl));
for j = find(~[Fl.dark])
    t = Fl(j); n = numel(t.ph); ph = angle(exp(1i*t.ph)); cue = angle(exp(1i*t.cue)); best = 0;
    for tt = -round(2.5*t.fr):round(1*t.fr)
        idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n);
        for sg = [-1 1], best = max(best, abs(mean(exp(1i*(ph(idx) + sg*cue(idx+tt)))))); end
    end
    RkL(j) = best;
end
onm = {'pk','mn','loc'}; olab = {'PB max','PB mean (global)','max - mean (local)'};
% frame-wise joint models: [fly + bump] and [|vr|, net rotation, forward, bump]
for sig = 1:2
    if sig == 1, FF = Fl; RR = RkL; ths = [0.10 0.15 0.20]; nmS = 'LPsP'; else, FF = Fd; RR = RkF; ths = rthr; nmS = 'GRAB'; end
    for th = ths, for oo = 1:3, for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | RR >= 0.4));
        if sig == 1, lf = lfL(d_+1); else, lf = 0.45; end
        B2 = nan(numel(trials),numel(lagsS),2); B4 = nan(numel(trials),numel(lagsS),4); NF = nan(numel(trials),numel(lagsS));
        for jj = 1:numel(trials)
            t = FF(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(lf*fr); y = t.(onm{oo});
            sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); net = abs(smoothdata(t.vr,'gaussian',w));
            fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
            good = movmin(t.rho_raw,w) > th & ~movmax(t.bad,w);
            for li_ = 1:numel(lagsS)
                L = round(lagsS(li_)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
                m = good(i0) & ~t.bad(i1) & ~isnan(y(i1)); if sum(m) < 50, continue, end
                Zq = zscore([sp(iF(m))' net(iF(m))' fw(iF(m))' bs(i0(m))' y(i1(m))']); X1 = ones(sum(m),1);
                b = [X1 Zq(:,[1 4])] \ Zq(:,5); B2(jj,li_,:) = b(2:3);
                b = [X1 Zq(:,1:4)] \ Zq(:,5);   B4(jj,li_,:) = b(2:5); NF(jj,li_) = sum(m);
            end
        end
        wm = @(X) sum(X.*NF,1,'omitnan')./sum(NF.*~isnan(X),1,'omitnan');
        W2 = wm(B2(:,:,2)); W1 = wm(B2(:,:,1)); win = find(lagsS >= -0.5 & lagsS <= 1); [~,k] = max(abs(W2(win))); jb = win(k);
        p2 = B2(:,jb,2); p4 = B4(:,jb,4); ok = ~isnan(p2); c4 = arrayfun(@(q) subsref(wm(B4(:,:,q)), struct('type','()','subs',{{jb}})), 1:4);
        fprintf('%s rho>%.2f %-19s %-11s bump lag %+.1f | fly %.2f bump %.2f (>0 %d/%d, p=%.2g) | 4-pred |vr| %.2f net %.2f fwd %.2f bump %.2f (>0 %d/%d, p=%.2g)\n', ...
            nmS, th, olab{oo}, ifelse_c(d_), lagsS(jb), W1(jb), W2(jb), sum(p2(ok)>0), sum(ok), signrank(p2(ok)), c4, sum(p4(ok)>0), sum(ok), signrank(p4(ok)));
    end, end, end
end
% matched pairs (as (5)/(6)), LPsP, rho > 0.15 in >= 80% of the bump frames:
% matched on |vr|, then on net rotation (consistency >= 0.8) + forward speed
thL = 0.15; mlabL = {'|vr|','net rot + fwd'};
for mm = 1:2
    P3 = [];
    for j = 1:numel(Fl)
        t = Fl(j); if ~t.dark && ~(RkL(j) >= 0.4), continue, end
        fr = t.fr; n = numel(t.ph); L = round(wsP*fr); st = round(0.5*fr); Lf = round(lfL(t.dark+1)*fr); Lb = round(0.2*fr);
        fly = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',max(3,round(fr))))*fr);
        good = t.rho_raw > thL & ~t.bad; rows = [];
        for a = 1+max(Lf,Lb):st:n-L
            g_ = a:a+L-1; fi = g_ - Lf; bi = g_ - Lb;
            if mean(good(bi)) < 0.8 || any(t.bad(g_)) || any(isnan(t.pk(g_))), continue, end
            rows(end+1,:) = [a mean(fly(fi)) abs(mean(t.vr(fi))) mean(t.fw(fi)) mean(bs(bi(good(bi)))) mean(t.pk(g_)) ...
                abs(mean(t.vr(fi)))/mean(abs(t.vr(fi))) mean(t.mn(g_)) mean(t.loc(g_))]; %#ok<AGROW>
        end
        for p1 = 1:size(rows,1)-1
            for p2 = p1+1:size(rows,1)
                if abs(rows(p1,1)-rows(p2,1)) < L, continue, end
                c = 2 + (mm == 2); A = rows(p1,c); Bq = rows(p2,c); if min(A,Bq) < 1 || abs(A-Bq)/mean([A Bq]) > 0.15, continue, end
                if mm == 2 && (min(rows([p1 p2],7)) < 0.8 || abs(rows(p1,4)-rows(p2,4)) > 0.15), continue, end
                if max(rows(p1,5),rows(p2,5))/max(min(rows(p1,5),rows(p2,5)),eps) < 1.5, continue, end
                if rows(p1,5) > rows(p2,5), q = [p2 p1]; else, q = [p1 p2]; end
                P3(end+1,:) = [j t.dark rows(q,6)' rows(q,8)' rows(q,9)' rows(q,4)' rows(q,7)']; %#ok<AGROW>
            end
        end
    end
    for d_ = 0:1
        sel = P3(:,2) == d_; tr = unique(P3(sel,1)); npt = arrayfun(@(k) sum(sel & P3(:,1)==k), tr); keep = npt >= 5;
        if sum(keep) < 3, continue, end
        f = @(c) arrayfun(@(k) mean(P3(sel & P3(:,1)==k,c+1) > P3(sel & P3(:,1)==k,c)), tr);
        rd = arrayfun(@(k) median((P3(sel & P3(:,1)==k,4) - P3(sel & P3(:,1)==k,3))./abs(P3(sel & P3(:,1)==k,3))), tr);
        pk_ = f(3); mn_ = f(5); lc_ = f(7);
        fprintf('LPsP pairs on %-13s %-11s %5d pairs %2d trials: faster bump higher PB max %.0f%% (p=%.2g, %+.0f%%) | PB mean %.0f%% (p=%.2g) | local %.0f%% (p=%.2g) | fwd %.2f/%.2f, consistency %.2f/%.2f (slow/fast)\n', ...
            mlabL{mm}, ifelse_c(d_), sum(sel), sum(keep), 100*median(pk_(keep)), signrank(pk_(keep)-0.5), 100*median(rd(keep)), ...
            100*median(mn_(keep)), signrank(mn_(keep)-0.5), 100*median(lc_(keep)), signrank(lc_(keep)-0.5), median(P3(sel,9)), median(P3(sel,10)), median(P3(sel,11)), median(P3(sel,12)));
    end
end

% (7b) multiplicative split. GRAB's trough is ~0 (trough/peak 0.12 in high-
% peak frames; LPsP 0.37), so a GRAB profile that just scales (gain x fixed
% shape) changes peak, mean AND max - mean together, with mean/peak fixed;
% the additive global/local split above doesn't fit GRAB. Here: gain = log
% PB max; shape = mean/peak (flatness), width (fraction of bins > half
% max), trough (min of the 3-bin-smoothed profile).
for sig = 1:2
    if sig == 1, FF = Fd; idxs = gi; else, FF = Fl; idxs = li; end
    for j = 1:numel(FF)
        s = B(idxs(j)); bad = logical(s.info.badVol(:)'); off = prctile(double(s.flash.meanImg(:)),0.1);
        bgs = smoothdata(double(s.flash.bg_pb(:)'),'gaussian',s.flash.params.smoothWindow);
        V = s.im.f - off - (bgs - median(bgs(~bad))); F0 = prctile(V(:,~bad),7,2); Ds = smoothdata((V-F0)./F0,1,'movmean',3,'omitnan');
        FF(j).tr = min(Ds,[],1); FF(j).wid = mean(Ds > 0.5*max(Ds,[],1),1,'omitnan');
        FF(j).ratio = FF(j).mn./FF(j).pk; FF(j).lpk = log(max(FF(j).pk,1e-3));
    end
    if sig == 1, Fd = FF; else, Fl = FF; end
end
onm2 = {'lpk','ratio','wid','tr'}; olab2 = {'log peak (gain)','mean/peak (flatness)','width >half max','trough'};
for sig = 1:2
    if sig == 1, FF = Fd; RR = RkF; th = rthr; nmS = 'GRAB'; else, FF = Fl; RR = RkL; th = thL; nmS = 'LPsP'; end
    for oo = 1:4, for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | RR >= 0.4)); if sig == 1, lf = 0.45; else, lf = lfL(d_+1); end
        B4 = nan(numel(trials),numel(lagsS),4); NF = nan(numel(trials),numel(lagsS));
        for jj = 1:numel(trials)
            t = FF(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(lf*fr); y = t.(onm2{oo});
            sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); net = abs(smoothdata(t.vr,'gaussian',w));
            fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
            good = movmin(t.rho_raw,w) > th & ~movmax(t.bad,w);
            for li_ = 1:numel(lagsS)
                L = round(lagsS(li_)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
                m = good(i0) & ~t.bad(i1) & isfinite(y(i1)) & t.pk(i1) > 0; if sum(m) < 50, continue, end
                Zq = zscore([sp(iF(m))' net(iF(m))' fw(iF(m))' bs(i0(m))' y(i1(m))']);
                b = [ones(sum(m),1) Zq(:,1:4)] \ Zq(:,5); B4(jj,li_,:) = b(2:5); NF(jj,li_) = sum(m);
            end
        end
        wm = @(X) sum(X.*NF,1,'omitnan')./sum(NF.*~isnan(X),1,'omitnan');
        W = wm(B4(:,:,4)); c1 = wm(B4(:,:,1)); win = find(lagsS >= -0.5 & lagsS <= 1); [~,k] = max(abs(W(win))); jb = win(k); p4 = B4(:,jb,4); ok = ~isnan(p4);
        fprintf('%s %-21s %-11s bump lag %+.1f | |vr| %+.2f | bump %+.2f (>0 %d/%d, p=%.2g)\n', nmS, olab2{oo}, ifelse_c(d_), lagsS(jb), c1(jb), W(jb), sum(p4(ok)>0), sum(ok), signrank(p4(ok)));
    end, end
end

%% (8) summary figure: fly speed vs bump speed -> peak and trough, per indicator
% Per trial, frame-wise joint model (standardized): outcome ~ |fly speed|
% (at the indicator's fly lag) + |bump speed| (folded bump, 1 s smoothing,
% clear-bump frames only: rho > 0.4 GRAB, EPG / > 0.15 LPsP). Outcomes: PB max
% (peak) and PB min of the 3-bin-smoothed profile (trough; for LPsP a
% readout of the global signal, for GRAB ~0 unless the profile broadens).
% Bump lag: the one maximizing the frame-weighted mean |coefficient| within
% -0.5..+1 s, per indicator x outcome x condition. Needs sections 4, 7, 7b.
% EPG GCaMP7f added for completeness: 23 trials (12 CL, 11 dark, 12 flies)
% after dropping the 9 with meanImg/bg_pb on a different scale from im.f
% (as in fig_summary_lpsp_grab_bump_model); rho > 0.4 (66% of frames);
% fly lag from its lag scan (+0.7 s CL, +0.6 s dark).
ei = find(strcmp(ind,'EPG>GCaMP7f')); ei = ei(:)'; keepE = true(size(ei));
for k = 1:numel(ei), s = B(ei(k)); bad = logical(s.info.badVol(:)'); keepE(k) = prctile(double(s.flash.meanImg(:)),0.1) < median(prctile(s.im.f(:,~bad),7,2)); end
ei = ei(keepE); lfE = [0.7 0.6]; Fe = arrayfun(@(k) grab_folded(B(k), snr_frac, a16), ei);
for j = 1:numel(Fe)
    s = B(ei(j)); bad = logical(s.info.badVol(:)'); off = prctile(double(s.flash.meanImg(:)),0.1);
    Fe(j).fw = interp1(s.ft.xf(:), smoothdata(s.ft.f_speed(:),'gaussian',30), s.ft.xb(:)', 'linear', 'extrap');
    bgs = smoothdata(double(s.flash.bg_pb(:)'),'gaussian',s.flash.params.smoothWindow);
    V = s.im.f - off - (bgs - median(bgs(~bad))); F0 = prctile(V(:,~bad),7,2); Fe(j).tr = min(smoothdata((V-F0)./F0,1,'movmean',3,'omitnan'),[],1);
end
RkE = nan(1,numel(Fe));
for j = find(~[Fe.dark])
    t = Fe(j); n = numel(t.ph); ph = angle(exp(1i*t.ph)); cue = angle(exp(1i*t.cue)); best = 0;
    for tt = -round(2.5*t.fr):round(1*t.fr)
        idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n);
        for sg = [-1 1], best = max(best, abs(mean(exp(1i*(ph(idx) + sg*cue(idx+tt)))))); end
    end
    RkE(j) = best;
end
S8 = struct(); onm8 = {'pk','tr'}; olab8 = {'peak (PB max)','trough (PB min)'}; ind8 = {'EPG GRAB(DA)','LPsP syt7f','EPG GCaMP7f'};
for sig = 1:3
    if sig == 1, FF = Fd; RR = RkF; th = rthr; elseif sig == 2, FF = Fl; RR = RkL; th = thL; else, FF = Fe; RR = RkE; th = rthr; end
    for oo = 1:2, for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | RR >= 0.4)); lfs = {[0.45 0.45], lfL, lfE}; lf = lfs{sig}(d_+1);
        B2 = nan(numel(trials),numel(lagsS),2); B4 = nan(numel(trials),numel(lagsS)); NF = nan(numel(trials),numel(lagsS));
        for jj = 1:numel(trials)
            t = FF(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(lf*fr); y = t.(onm8{oo});
            sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); net = abs(smoothdata(t.vr,'gaussian',w));
            fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
            good = movmin(t.rho_raw,w) > th & ~movmax(t.bad,w);
            for li_ = 1:numel(lagsS)
                L = round(lagsS(li_)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
                m = good(i0) & ~t.bad(i1) & isfinite(y(i1)); if sum(m) < 50, continue, end
                Zq = zscore([sp(iF(m))' bs(i0(m))' net(iF(m))' fw(iF(m))' y(i1(m))']); X1 = ones(sum(m),1);
                b = [X1 Zq(:,1:2)] \ Zq(:,5); B2(jj,li_,:) = b(2:3);
                b = [X1 Zq(:,1:4)] \ Zq(:,5); B4(jj,li_) = b(3); NF(jj,li_) = sum(m);
            end
        end
        wm = @(X) sum(X.*NF,1,'omitnan')./sum(NF.*~isnan(X),1,'omitnan');
        W = wm(B2(:,:,2)); win = find(lagsS >= -0.5 & lagsS <= 1); [~,k] = max(abs(W(win))); jb = win(k);
        S8(sig,oo,d_+1).fly = B2(:,jb,1); S8(sig,oo,d_+1).bump = B2(:,jb,2); S8(sig,oo,d_+1).bump4 = B4(:,jb); S8(sig,oo,d_+1).lag = lagsS(jb);
        ok = ~isnan(B2(:,jb,2)); S8(sig,oo,d_+1).p = signrank(B2(ok,jb,2)); S8(sig,oo,d_+1).pf = signrank(B2(ok,jb,1));
        fprintf('%-12s %-15s %-11s bump lag %+.1f | fly %.2f bump %.2f (>0 %d/%d, p=%.2g) | bump, 4-pred %.2f\n', ind8{sig}, olab8{oo}, ifelse_c(d_), lagsS(jb), ...
            median(B2(ok,jb,1)), median(B2(ok,jb,2)), sum(B2(ok,jb,2)>0), sum(ok), S8(sig,oo,d_+1).p, median(B4(ok,jb),'omitnan'));
    end, end
end
figure('Color','w','Position',[80 80 1400 680]); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cF = [0.15 0.15 0.15]; cB = [0 0.45 0.74]; pstr = @(p) ifelse_p(p);
for oo = 1:2, for sig = 1:3
    ax = nexttile(tl, (oo-1)*3 + sig); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]);
    for d_ = 0:1
        Q = S8(sig,oo,d_+1); xs = d_*3 + [1 2]; vals = {Q.fly, Q.bump}; cc = {cF, cB};
        for e = 1:2
            v = vals{e}; v = v(~isnan(v)); jx = xs(e) + 0.12*randn(size(v)); jx = min(max(jx, xs(e)-0.3), xs(e)+0.3);
            scatter(ax, jx, v, 18, cc{e}, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
            plot(ax, xs(e) + [-0.32 0.32], median(v)*[1 1], '-', 'Color', cc{e}, 'LineWidth', 3);
        end
        text(ax, xs(2), 0.97, sprintf('%s\nlag %+.1f s', pstr(Q.p), Q.lag), 'Units','data','HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color',cB);
    end
    xlim(ax, [0.3 5.7]); ylim(ax, [-0.5 1]); set(ax, 'XTick', [1 2 4 5], 'XTickLabel', {'fly','bump','fly','bump'}, 'Box','off', 'TickDir','out');
    text(ax, 1.5, -0.47, 'closed loop', 'HorizontalAlignment','center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]); text(ax, 4.5, -0.47, 'dark', 'HorizontalAlignment','center', 'FontSize', 9, 'Color', [0.4 0.4 0.4]);
    if sig == 1, ylabel(ax, sprintf('%s\njoint coefficient (z)', olab8{oo})); end
    if oo == 1, title(ax, ind8{sig}, 'FontWeight','normal'); end
end, end
title(tl, 'Fly speed vs bump speed as joint predictors of PB fluorescence (dots = trials, bar = median)', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_peak_trough.png'),'Resolution',150);

%% (9) is the trough effect just bump widening? far-field readout of the speed signal
% Trough/flatness can rise with bump speed simply because a moving bump
% (plus slow indicator decay) leaves a trail. Far field = mean folded dF/F
% over bins that stayed > 90 deg (or 135) from the bump phase throughout
% the last 2 s (or 4 s), i.e. never touched by the bump or its trail. For
% LPsP this is the candidate readout of the neural (commanded) speed signal.
% Model as section 6/8: outcome ~ |vr| (fly lag) + net rotation + forward
% speed + bump speed (+ extras), per trial, clear-bump frames only.
% (a) does the far field still scale with bump speed beyond the ball?
% (b) controls: cross-half (bump speed from one PB half, far field of the
%     other), bump concentration (rho) as a covariate, and net bump
%     displacement over 1.5 s instead of speed (phase jitter mostly cancels).
% (c) given that neural readout, does bump speed still predict the peak?
for j = 1:numel(Fd), Fd(j).net = abs(smoothdata(Fd(j).vr,'gaussian',max(3,round(Fd(j).fr)))); end
for j = 1:numel(Fl), Fl(j).net = abs(smoothdata(Fl(j).vr,'gaussian',max(3,round(Fl(j).fr)))); end
for j = 1:numel(Fe), Fe(j).net = abs(smoothdata(Fe(j).vr,'gaussian',max(3,round(Fe(j).fr)))); end
SETS = {Fd, gi, RkF, rthr, [0.45 0.45], 'GRAB'; Fl, li, RkL, thL, lfL, 'LPsP'; Fe, ei, RkE, rthr, lfE, 'EPG7f'};
ctl = {'fly:net','fly:fw'};
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2}; th = SETS{sgi,4};
    for j = 1:numel(FF)
        [FF(j).far, ~, FF(j).anti] = far_field(B(idx(j)), FF(j), snr_frac, 2, 90);
        FF(j).far135 = far_field(B(idx(j)), FF(j), snr_frac, 2, 135); FF(j).far4 = far_field(B(idx(j)), FF(j), snr_frac, 4, 90);
        H = far_half(B(idx(j)), FF(j).fr, snr_frac, 2, 90); FF(j).H = H;
        FF(j).pk1 = H(1).pk; FF(j).pk2 = H(2).pk; FF(j).far1 = H(1).far; FF(j).far2 = H(2).far; FF(j).loc = FF(j).pk - FF(j).far;
    end
    SETS{sgi,1} = FF;
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); lf = SETS{sgi,5}(d_+1); nmS = sprintf('%-5s %-11s', SETS{sgi,6}, ifelse_c(d_));
        pr = @(lab, C, k) fprintf('%s %-34s fly %+.2f bump %+.2f (>0 %d/%d, p=%.2g)%s\n', nmS, lab, median(C(~isnan(C(:,k)),1)), median(C(~isnan(C(:,k)),k)), ...
            sum(C(~isnan(C(:,k)),k) > 0), sum(~isnan(C(:,k))), signrank(C(~isnan(C(:,k)),k)), '');
        % (a) far field vs bump speed
        for f = {'tr','far','far135','far4'}
            C = joint_coef(FF, trials, f{1}, lf, th, ctl, lagsS, [-0.5 2]); pr(['(a) ' f{1}], C, 2);
        end
        % (b) cross-half, + rho control, displacement
        R = nan(numel(trials),3,2);
        for dirn = 1:2
            G = FF(trials);
            for jj = 1:numel(G)
                a = dirn; b = 3 - dirn; fr = G(jj).fr; W = round(1.5*fr); H = G(jj).H;
                G(jj).ph = H(a).ph; G(jj).rho_raw = H(a).rho; G(jj).farx = H(b).far;
                G(jj).rhoP = smoothdata(H(a).rho,'gaussian',max(3,round(fr)));
                G(jj).disp = [nan(1,W) abs(H(a).ph(1+W:end) - H(a).ph(1:end-W))]/1.5;
            end
            C = joint_coef(G, 1:numel(G), 'farx', lf, th, ctl, lagsS, [-0.5 2]);             R(:,1,dirn) = C(:,2);
            C = joint_coef(G, 1:numel(G), 'farx', lf, th, [ctl {'rhoP'}], lagsS, [-0.5 2]);  R(:,2,dirn) = C(:,2);
            C = joint_coef(G, 1:numel(G), 'farx', lf, th, [ctl {'disp'}], lagsS, [-0.5 2]);  R(:,3,dirn) = C(:,5);
        end
        Rm = mean(R,3,'omitnan'); lab = {'(b) cross-half far','(b) cross-half far + rho ctrl','(b) cross-half far, net displacement'};
        for k = 1:3, ok = ~isnan(Rm(:,k)); fprintf('%s %-34s bump %+.2f (>0 %d/%d, p=%.2g)\n', nmS, lab{k}, median(Rm(ok,k)), sum(Rm(ok,k)>0), sum(ok), signrank(Rm(ok,k))); end
        % (c) peak given the far-field readout (same-signal and cross-half)
        C1 = joint_coef(FF, trials, 'pk', lf, th, [ctl {'far'}], lagsS);
        Cx = (joint_coef(FF, trials, 'pk1', lf, th, [ctl {'far2'}], lagsS) + joint_coef(FF, trials, 'pk2', lf, th, [ctl {'far1'}], lagsS))/2;
        ok = ~isnan(C1(:,2)) & ~isnan(Cx(:,2));
        fprintf('%s (c) peak ~ ... + far: far %+.2f, bump %+.2f (>0 %d/%d, p=%.2g) | cross-half far %+.2f, bump %+.2f (>0 %d/%d, p=%.2g)\n', nmS, ...
            median(C1(ok,5)), median(C1(ok,2)), sum(C1(ok,2)>0), sum(ok), signrank(C1(ok,2)), median(Cx(ok,5)), median(Cx(ok,2)), sum(Cx(ok,2)>0), sum(ok), signrank(Cx(ok,2)));
        if sgi == 2
            C2 = joint_coef(FF, trials, 'loc', lf, th, [ctl {'far'}], lagsS); ok = ~isnan(C2(:,2));
            fprintf('%s (c) local (peak - far) ~ ... + far:  far %+.2f, bump %+.2f (>0 %d/%d, p=%.2g)\n', nmS, median(C2(ok,5)), median(C2(ok,2)), sum(C2(ok,2)>0), sum(ok), signrank(C2(ok,2)));
        end
    end
end
Fd = SETS{1,1}; Fl = SETS{2,1}; Fe = SETS{3,1};

%% (10) strict far-field trough: peak vs trough, fly vs bump speed (figure)
% Strict trough = min of the 3-bin (circular) smoothed folded profile over
% bins that stayed > 135 deg from the bump phase throughout the last 4 s
% (far_field, 5th output). Only frames where it is defined are used, in all
% three rows (so rows are comparable): (1) strict trough ~ |fly speed| +
% bump speed; (2) peak ~ same; (3) peak ~ same + strict trough (the speed-
% signal readout). Figure: fly + bump only (with net rotation and forward
% speed added, |vr| and net split one shared effect); the bump term with
% those controls is printed alongside ('ctrl'). Coefficients per trial; bump
% lag chosen as in joint_coef (-0.5..+1 s; trough rows -0.5..+2 s).
% Where the bump term is ~0 (row 3, LPsP/GRAB dark) the selected lag is
% arbitrary (window edge). NB: requiring untouched bins over 4 s drops fast, wide sweeps, so the
% frames are biased to slower bump motion; GRAB keeps only ~9% of its
% clear frames (4-7 trials) -> shown, but underpowered.
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2};
    for j = 1:numel(FF), [~,~,~,~,FF(j).ftr4] = far_field(B(idx(j)), FF(j), snr_frac, 4, 135); end
    SETS{sgi,1} = FF;
end
S10 = trough_fig(SETS, 'ftr4', lagsS, 'Strict far-field trough (bins > 135{\circ} from the bump for the last 4 s): does bump speed add beyond the speed signal?', ...
    fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_strict_trough.png'));
Fd = SETS{1,1}; Fl = SETS{2,1}; Fe = SETS{3,1};

%% (11) per-frame anti-phase trough (moves with the bump) + trail-only exclusion
% Section 10's strict trough required bins to stay > 135 deg from EVERY bump
% position in the last 4 s: allowed arc = 90 deg minus the 4 s sweep, so any
% sweep > ~60 deg leaves no bins. Clear GRAB frames sweep a median 138 deg
% (they are big turns: fly speed 4.9x the other frames), so only ~6% kept.
% Here the trough is taken opposite the CURRENT bump (> 135 deg from it,
% recomputed every frame) and only bins the recent trail came within 90 deg
% of are dropped: atr0 = no trail exclusion, atr2 / atr4 = trail over the
% last 2 / 4 s. Rows as in section 10 (trough ~ fly + bump; peak ~ fly +
% trough + bump, same frames).
vars11 = {'atr0', 0, 0; 'atr2', 2, 90; 'atr4', 4, 90};
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2}; th = SETS{sgi,4};
    for j = 1:numel(FF)
        [~,~,~,Df] = far_field(B(idx(j)), FF(j), snr_frac, 2, 90);
        for v = 1:3, FF(j).(vars11{v,1}) = anti_trough(Df, FF(j).ph, FF(j).fr, 135, vars11{v,2}, vars11{v,3}); end
    end
    SETS{sgi,1} = FF;
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); lf = SETS{sgi,5}(d_+1);
        for v = {'atr0','atr2','atr4','ftr4'}
            cv = arrayfun(@(k) mean(~isnan(FF(k).(v{1})(movmin(FF(k).rho_raw,max(3,round(FF(k).fr))) > th & ~FF(k).bad))), trials);
            C = joint_coef(FF, trials, v{1}, lf, th, {}, lagsS, [-0.5 2]); ok = ~isnan(C(:,2));
            for j = 1:numel(FF), FF(j).pkV = FF(j).pk; FF(j).pkV(isnan(FF(j).(v{1}))) = nan; end
            C3 = joint_coef(FF, trials, 'pkV', lf, th, v, lagsS); ok3 = ~isnan(C3(:,2));
            fprintf('%-5s %-11s %-4s defined %3.0f%% of clear | trough: fly %+.2f bump %+.2f (>0 %d/%d, p=%.2g) | peak+trough: trough %+.2f bump %+.2f (>0 %d/%d, p=%.2g)\n', ...
                SETS{sgi,6}, ifelse_c(d_), v{1}, 100*median(cv), median(C(ok,1)), median(C(ok,2)), sum(C(ok,2)>0), sum(ok), signrank(C(ok,2)), ...
                median(C3(ok3,3)), median(C3(ok3,2)), sum(C3(ok3,2)>0), sum(ok3), signrank(C3(ok3,2)));
        end
    end
end
Fd = SETS{1,1}; Fl = SETS{2,1}; Fe = SETS{3,1};

% (11b) trail window from each indicator's measured decay instead of a
% fixed 2/4 s. Decay: for every bin currently > 90 deg from the bump,
% dF/F (normalized to the trial's 90th-pct peak) vs time since the bump
% (clear frames) was last within 22.5 deg of it; median across trials;
% exponential + floor fit. Bins > 90 deg away can't have been visited
% < ~0.9 s ago, so the fit starts there. Upper bound on sensor decay: a
% far bin visited recently also means the fly was just turning (and LPsP
% far bins carry the speed-driven global signal) -> LPsP fit is poor.
a16 = linspace(-pi, pi-2*pi/16, 16)'; tb = 0:0.25:8; TAU = nan(1,3);
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2}; th = SETS{sgi,4}; prof = nan(numel(FF), numel(tb)-1);
    for j = 1:numel(FF)
        t = FF(j); fr = t.fr; n = numel(t.ph); [~,~,~,Df] = far_field(B(idx(j)), t, snr_frac, 2, 90);
        Dn = Df / prctile(max(Df,[],1), 90); ph = angle(exp(1i*t.ph)); good = movmin(t.rho_raw,max(3,round(fr))) > th & ~t.bad;
        d = abs(angle(exp(1i*(a16 - ph)))); last = nan(16,1); since = nan(16,n);
        for k = 1:n, if good(k), last(d(:,k) < pi/8) = k; end; since(:,k) = (k - last)/fr; end
        ok = d > pi/2 & repmat(~t.bad,16,1) & ~isnan(since); [~,~,bi] = histcounts(since(ok), tb); v = Dn(ok); m = bi > 0 & ~isnan(v);
        cnt = accumarray(bi(m), 1, [numel(tb)-1 1])'; pj = accumarray(bi(m), v(m), [numel(tb)-1 1], @mean, nan)'; pj(cnt < 20) = nan; prof(j,:) = pj;
    end
    y = median(prof,1,'omitnan'); y(sum(~isnan(prof),1) < 0.5*numel(FF)) = nan; x = tb(1:end-1) + 0.125; okx = ~isnan(y); xs = x(okx); ys = y(okx);
    p = lsqcurvefit(@(p,xx) p(1)*exp(-(xx - xs(1))/p(2)) + p(3), [ys(1)-ys(end) 1 ys(end)], xs, ys, [0 0.05 -1], [5 20 2], optimoptions('lsqcurvefit','Display','off'));
    TAU(sgi) = p(2);
    fprintf('%-5s decay after the bump leaves a bin: tau = %.2f s (fit from %.1f s; floor %.2f of peak)\n', SETS{sgi,6}, p(2), xs(1), p(3));
end
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2}; th = SETS{sgi,4}; tau = TAU(sgi);
    for j = 1:numel(FF)
        [~,~,~,Df] = far_field(B(idx(j)), FF(j), snr_frac, 2, 90);
        for m_ = 1:3, FF(j).(sprintf('atau%d',m_)) = anti_trough(Df, FF(j).ph, FF(j).fr, 135, m_*tau, 90); end
    end
    SETS{sgi,1} = FF;
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); lf = SETS{sgi,5}(d_+1);
        for m_ = 1:3
            v = sprintf('atau%d',m_);
            cv = arrayfun(@(k) mean(~isnan(FF(k).(v)(movmin(FF(k).rho_raw,max(3,round(FF(k).fr))) > th & ~FF(k).bad))), trials);
            C = joint_coef(FF, trials, v, lf, th, {}, lagsS, [-0.5 2]); ok = ~isnan(C(:,2));
            for j = 1:numel(FF), FF(j).pkV = FF(j).pk; FF(j).pkV(isnan(FF(j).(v))) = nan; end
            C3 = joint_coef(FF, trials, 'pkV', lf, th, {v}, lagsS); ok3 = ~isnan(C3(:,2));
            fprintf('%-5s %-11s trail %d tau (%.1f s): defined %3.0f%% | trough: fly %+.2f bump %+.2f (>0 %d/%d, p=%.2g) | peak+trough: trough %+.2f bump %+.2f (>0 %d/%d, p=%.2g)\n', ...
                SETS{sgi,6}, ifelse_c(d_), m_, m_*tau, 100*median(cv), median(C(ok,1)), median(C(ok,2)), sum(C(ok,2)>0), sum(ok), signrank(C(ok,2)), ...
                median(C3(ok3,3)), median(C3(ok3,2)), sum(C3(ok3,2)>0), sum(ok3), signrank(C3(ok3,2)));
        end
    end
end
S11 = trough_fig(SETS, 'atau3', lagsS, sprintf('Trough opposite the current bump (> 135{\\circ}), excluding bins the trail came within 90{\\circ} of in the last 3\\tau (GRAB %.1f s, LPsP %.1f s, EPG7f %.1f s)', 3*TAU), ...
    fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_anti_trough_3tau.png'));
Fd = SETS{1,1}; Fl = SETS{2,1}; Fe = SETS{3,1};

%% (12) does the LPsP trough ("commanded speed") predict bump position better than the ball?
% LPsP bump displacement over non-overlapping 1 / 2 s windows (clear bump,
% rho > 0.15, in >= 80% of the window) predicted from: ball = integrated
% ball velocity; + dev = also integrated ball velocity x trough deviation
% (trough above / below what ball speed predicts, z; trough advanced by its
% fly lag); trough-set = integrated sign(vr) x trough-calibrated speed.
% Trough = far field (90 deg / 2 s) or the 3-tau anti-phase trough. 2-fold
% CV over windows (bump lag chosen on training windows); in-sample dev
% coefficient (standardized, with ball) vs a time-shifted (1/3 trial) dev.
FF = SETS{2,1};
for trF = {'far','atau3'}, for W = [1 2], for d_ = 0:1
    trials = find([FF.dark] == d_ & (d_ == 1 | RkL >= 0.4)); lf = lfL(d_+1); out = []; outS = [];
    for j = trials
        R = bump_from_speed(FF(j), trF{1}, lf, thL, W, false); S_ = bump_from_speed(FF(j), trF{1}, lf, thL, W, true);
        if isempty(R) || isempty(S_), continue, end
        out(end+1,:) = [R.r2 R.b R.bball R.n]; outS(end+1,:) = [S_.r2 S_.b]; %#ok<AGROW>
    end
    dR = out(:,2) - out(:,1); okR = ~isnan(dR);
    fprintf('%-5s W=%d s %-11s %2d trials | CV R2 ball %.2f -> +dev %.2f (better %d/%d) | trough-set speed %.2f | dev coef %+.2f (>0 %d/%d, p=%.2g) vs shuffled %+.2f\n', ...
        trF{1}, W, ifelse_c(d_), size(out,1), median(out(okR,1)), median(out(okR,2)), sum(dR(okR)>0), sum(okR), median(out(:,3),'omitnan'), ...
        median(out(:,4)), sum(out(:,4)>0), size(out,1), signrank(out(:,4)), median(outS(:,4)));
end, end, end

% (12b) dark: recover frames. Clarity only at the two window endpoints
% (rho over +-0.25 s > 0.10 / 0.15), overlapping windows (0.5 s step),
% 1-8 s windows, leave-one-block-out CV (4 contiguous blocks).
FF = SETS{2,1}; trials = find([FF.dark] == 1); lf = lfL(2);
for th_ = [0.10 0.15], for trF = {'far','atau3'}, for W = [1 2 4 8]
    out = []; outS = [];
    for j = trials
        R = bump_from_speed2(FF(j), trF{1}, lf, th_, W, false, 4); S_ = bump_from_speed2(FF(j), trF{1}, lf, th_, W, true, 4);
        if isempty(R) || isempty(S_), continue, end
        out(end+1,:) = [R.r2 R.b R.bball R.nind]; outS(end+1,:) = [S_.r2 S_.b]; %#ok<AGROW>
    end
    if size(out,1) < 3, continue, end
    dR = out(:,2) - out(:,1);
    fprintf('dark rho>%.2f %-5s W=%d s: %2d trials (~%3.0f indep. windows) | CV R2 ball %.2f -> +dev %.2f (better %d/%d) | trough-set %.2f | dev coef %+.2f (>0 %d/%d, p=%.2g) vs shuffled %+.2f\n', ...
        th_, trF{1}, W, size(out,1), median(out(:,6)), median(out(:,1)), median(out(:,2)), sum(dR>0), numel(dR), median(out(:,3)), ...
        median(out(:,4)), sum(out(:,4)>0), size(out,1), signrank(out(:,4)), median(outS(:,4)));
end, end, end
% is there a bump trajectory to predict? 2 s displacement: PB half 1 vs half 2
for sgi = 1:3
    FF = SETS{sgi,1}; th = SETS{sgi,4}; if sgi == 2, th = 0.10; end
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); rr = []; rb = [];
        for j = trials
            H = FF(j).H; t = FF(j); fr = t.fr; Lw = round(2*fr); h = round(0.25*fr);
            c = movmean(H(1).rho,2*h+1) > th & movmean(H(2).rho,2*h+1) > th & ~movmax(t.bad,2*h+1);
            st = 1:round(0.5*fr):numel(t.ph)-Lw; e = st + Lw - 1; k = c(st) & c(e); if sum(k) < 15, continue, end
            d1 = H(1).ph(e(k)) - H(1).ph(st(k)); d2 = H(2).ph(e(k)) - H(2).ph(st(k)); hd = cumsum(t.vr)/fr; db = hd(e(k)) - hd(st(k));
            rr(end+1) = corr(d1', d2'); rb(end+1) = corr(((d1+d2)/2)', db'); %#ok<AGROW>
        end
        fprintf('%-5s %-11s 2 s bump displacement: half 1 vs half 2 r = %.2f | mean of halves vs ball r = %.2f (%d trials)\n', SETS{sgi,6}, ifelse_c(d_), median(rr), median(rb), numel(rr));
    end
end

% (12c) as 12b but displacement = WRAPPED endpoint phase difference
% (unwrapped phase slips by 2 pi wherever the bump is faint mid-window,
% adding fake displacement even with clear endpoints); 1 / 2 s windows,
% windows with |ball heading change| > 120 deg dropped. Endpoint clarity
% rho > 0 (none) / 0.10 / 0.15. Trough = far field or the anti-phase
% trough without trail exclusion (atr0; defined in all frames).
% Half 1 vs half 2 agreement on wrapped displacement, all frames: LPsP
% 0.28 / 0.44 (CL 1 / 2 s), 0.15 / 0.28 (dark); GRAB 0.34 / 0.40, 0.23 /
% 0.28; EPG7f 0.96 / 0.97, 0.84 / 0.88.
FF = SETS{2,1}; for j = 1:numel(FF), FF(j).atau0 = FF(j).atr0; end
for d_ = [1 0], for th_ = [0 0.10 0.15], for trF = {'far','atau0'}, for W = [1 2]
    trials = find([FF.dark] == d_ & (d_ == 1 | RkL >= 0.4)); lf = lfL(d_+1); out = []; outS = [];
    for j = trials
        R = bump_from_speed2(FF(j), trF{1}, lf, th_, W, false, 4, true); S_ = bump_from_speed2(FF(j), trF{1}, lf, th_, W, true, 4, true);
        if isempty(R) || isempty(S_), continue, end
        out(end+1,:) = [R.r2 R.b R.bball R.nind]; outS(end+1,:) = [S_.r2 S_.b]; %#ok<AGROW>
    end
    if size(out,1) < 3, continue, end
    dR = out(:,2) - out(:,1);
    fprintf('%-11s rho>%.2f %-5s W=%d s: %2d trials (~%3.0f windows) | CV R2 ball %.2f -> +dev %.2f (better %d/%d) | trough-set %.2f | dev coef %+.2f (>0 %d/%d, p=%.2g) vs shuffled %+.2f | ball coef %.2f\n', ...
        ifelse_c(d_), th_, trF{1}, W, size(out,1), median(out(:,6)), median(out(:,1)), median(out(:,2)), sum(dR>0), numel(dR), median(out(:,3)), ...
        median(out(:,4)), sum(out(:,4)>0), size(out,1), signrank(out(:,4)), median(outS(:,4)), median(out(:,5)));
end, end, end, end

%% (13) cross-centered PB profiles by fly speed, per indicator, CL vs dark (figure)
% Is the LPsP bump less bump-shaped in the dark? All non-bad frames (no
% clarity selection); each PB half centred on the OTHER half's bump phase
% (+/- the trial's left-right offset) so noise can't create a bump; binned
% by |fly speed| at each indicator's fly lag; mean over frames per trial,
% then mean +- SEM over trials. Phase errors blur the cross-centred average,
% so bump height here is a lower bound. Also: per-frame noise per bin
% (std of frame-to-frame differences / sqrt 2) and bump height / noise.
edges = [0 0.25 0.5 1 2 inf]; nb = numel(edges) - 1; PR = cell(3,2); NZ = nan(3,2);
for sgi = 1:3
    FF = SETS{sgi,1}; idx = SETS{sgi,2};
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); lf = SETS{sgi,5}(d_+1);
        Pall = nan(16,nb,numel(trials)); nz = nan(1,numel(trials));
        for jj = 1:numel(trials)
            j = trials(jj); s = B(idx(j)); Pall(:,:,jj) = xcentered_profiles(s, FF(j).fr, FF(j).vr, snr_frac, lf, edges);
            bad = logical(s.info.badVol(:)'); off = prctile(double(s.flash.meanImg(:)),0.1);
            bgs = smoothdata(double(s.flash.bg_pb(:)'),'gaussian',s.flash.params.smoothWindow);
            V = s.im.f - off - (bgs - median(bgs(~bad))); F0 = prctile(V(:,~bad),7,2); D = (V-F0)./F0;
            nz(jj) = median(std(diff(D(:,~bad),1,2),0,2,'omitnan'))/sqrt(2);
        end
        PR{sgi,d_+1} = Pall; NZ(sgi,d_+1) = median(nz);
        M = mean(Pall,3,'omitnan'); con = max(M,[],1) - min(M,[],1);
        fprintf('%-5s %-11s %2d trials | bump height (peak - trough, dF/F) by speed bin %s | / mean %s | / frame noise (%.3f) %s\n', SETS{sgi,6}, ifelse_c(d_), numel(trials), ...
            mat2str(round(con,2)), mat2str(round(con./mean(M,1),2)), NZ(sgi,d_+1), mat2str(round(con/NZ(sgi,d_+1),1)));
    end
end
figure('Color','w','Position',[60 60 1300 760]); tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cmap = [0.78 0.85 0.95; 0.55 0.70 0.90; 0.30 0.52 0.82; 0.12 0.33 0.68; 0.03 0.15 0.45];
xd = [-180:22.5:157.5 180]; ttl13 = {'EPG GRAB(DA)','LPsP syt7f','EPG GCaMP7f'};
blab = arrayfun(@(k) sprintf('%.2g-%.2g rad/s', edges(k), edges(k+1)), 1:nb, 'UniformOutput', false); blab{end} = sprintf('> %g rad/s', edges(end-1));
for sgi = 1:3
    yl = [inf -inf];
    for d_ = 0:1
        ax(d_+1) = nexttile(tl, d_*3 + sgi); hold(ax(d_+1),'on'); %#ok<SAGROW>
        Pall = PR{sgi,d_+1}; M = mean(Pall,3,'omitnan'); E = std(Pall,0,3,'omitnan')./sqrt(sum(~isnan(Pall),3));
        for k = 1:nb
            m = [M(:,k); M(1,k)]'; e = [E(:,k); E(1,k)]';
            fill(ax(d_+1), [xd fliplr(xd)], [m-e fliplr(m+e)], cmap(k,:), 'EdgeColor','none','FaceAlpha',0.25,'HandleVisibility','off');
            plot(ax(d_+1), xd, m, '-', 'Color', cmap(k,:), 'LineWidth', 2);
        end
        yl = [min(yl(1), min(M-E,[],'all')) max(yl(2), max(M+E,[],'all'))];
        con = max(M,[],1) - min(M,[],1);
        text(ax(d_+1), 0.02, 0.97, sprintf('bump height / frame noise:\n%s', strjoin(arrayfun(@(v) sprintf('%.1f',v), con/NZ(sgi,d_+1), 'UniformOutput', false), ', ')), ...
            'Units','normalized','VerticalAlignment','top','FontSize',8,'Color',[0.3 0.3 0.3]);
        xlim(ax(d_+1), [-180 180]); set(ax(d_+1),'XTick',-180:90:180,'Box','off','TickDir','out');
        if d_ == 1, xlabel(ax(d_+1), 'position relative to the other half''s bump (deg)'); end
        if sgi == 1, ylabel(ax(d_+1), sprintf('%s\ndF/F', ifelse_c(d_))); end
        if d_ == 0, title(ax(d_+1), ttl13{sgi}, 'FontWeight','normal'); end
    end
    pad = 0.05*diff(yl); set(ax, 'YLim', [yl(1)-pad yl(2)+3*pad]);
    if sgi == 3, legend(ax(1), blab, 'Box','off','Location','northeast','FontSize',8); end
end
title(tl, 'Cross-centred PB profile by |fly speed| (all frames; mean \pm SEM over trials)', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_profiles_by_speed.png'),'Resolution',150);

%% (14) is GRAB a "motion energy" of LPsP calcium? leading edge and turn-onset tests (figure)
% Alternative to GRAB = speed x bump: release driven by CHANGE in LPsP
% calcium. Two kinds of energy on L(theta,t) ~ G(t) + k E(theta - phi(t)):
%  A oriented (Adelson-Bergen / Reichardt): blind to the uniform G(t);
%    centred on the bump, scales with bump speed -> same position as the
%    product model, separable only by fly- vs bump-speed dependence
%    (sections 4-8: fly speed dominates).
%  B temporal change ([dL/dt]+, adapting synapse): a PB-wide transient when
%    G rises (turn onset overshoot, fast offset) plus release at the bump's
%    LEADING edge.
% Test 1 (closed loop): GRAB profile aligned on the CUE-predicted bump
%   position (time lag from EPG7f, -0.44 s; only sign and a constant offset
%   fitted per trial -- a lag fitted on GRAB would absorb a speed-dependent
%   shift), mirrored so travel is +, binned by bump (cue) speed; centroid
%   of profile - min. Compared with models built on the same trajectories
%   (template = GRAB's slow-speed profile): product (template x |speed|)
%   and rectified d/dt, each low-passed with tau 0-1.1 s.
% Test 2 (CL + dark): turn onsets / offsets (|vr| crossing 1 rad/s, < 0.3
%   for 1 s before / > 1 for 1 s after, or the mirror); residual of y vs
%   a fitted low-passed |vr| model (lag 0-1 s, tau 0.05-1.2 s), in SD of y.
edges = [0 0.25 0.5 1 2 inf]; vc = [0.12 0.37 0.75 1.5 3];
QE = []; for j = find(~[Fe.dark]), [r_, l_] = cue_track(Fe(j)); QE(end+1,:) = [r_ l_]; end %#ok<AGROW>
lagE = median(QE(QE(:,1) >= 0.4, 2)); fprintf('cue -> bump lag (EPG7f, closed loop): %.2f s\n', lagE);
SETS14 = {Fd, gi, RkF, 'GRAB'; Fl, li, RkL, 'LPsP'; Fe, ei, RkE, 'EPG7f'}; CR = cell(3,1);
for sgi = 1:3
    FF = SETS14{sgi,1}; idx = SETS14{sgi,2}; trials = find(~[FF.dark] & SETS14{sgi,3} >= 0.4); C = nan(numel(trials), 5); PP = nan(16,5,numel(trials));
    for jj = 1:numel(trials)
        j = trials(jj); [PP(:,:,jj), ~, C(jj,:)] = cue_ref_profiles(B(idx(j)), FF(j), snr_frac, lagE, edges);
    end
    CR{sgi} = struct('C', C, 'P', PP);
    fprintf('%-5s CL (%2d trials) centroid offset vs cue-predicted bump (deg, + = ahead): %s\n', SETS14{sgi,4}, numel(trials), mat2str(round(median(C,1,'omitnan'))));
end
M = mean(CR{1}.P,3,'omitnan'); tp = M(:,2) - min(M(:,2)); tp = (tp + tp(mod(17-(1:16),16)+1))/2; tp = tp/max(tp); [~,i0] = max(tp); tp = circshift(tp, 9 - i0);
tmpl = @(th) interp1([a16; pi], [tp; tp(1)], mod(th + pi, 2*pi) - pi, 'linear');
taus = [0.1 0.3 0.5 1.1]; trials = find(~[Fd.dark] & RkF >= 0.4); SIM = nan(numel(taus), 5, 2);
for it = 1:numel(taus)
    Cm = nan(numel(trials),5,2);
    for jj = 1:numel(trials)
        j = trials(jj); t = Fd(j); s = B(gi(j)); n = numel(t.ph); fr = t.fr;
        [~,~,~,ref,v] = cue_ref_profiles(s, t, snr_frac, lagE, edges, zeros(16,n));
        E = tmpl(a16 - fillmissing(ref,'nearest')); spd = fillmissing(smoothdata(abs(v),'gaussian',max(3,round(0.5*fr))),'nearest');
        mods = {E .* spd, max(gradient(E)*fr, 0)};
        k = exp(-(0:round(5*taus(it)*fr))/(taus(it)*fr)); k = k/sum(k);
        for mm = 1:2, [~,~,Cm(jj,:,mm)] = cue_ref_profiles(s, t, snr_frac, lagE, edges, filter(k, 1, mods{mm}, [], 2)); end
    end
    SIM(it,:,:) = median(Cm,1,'omitnan');
    fprintf('model, tau %.1f s | product %s | rectified d/dt %s\n', taus(it), mat2str(round(SIM(it,:,1))), mat2str(round(SIM(it,:,2))));
end
win = [-2 4]; tg = win(1):0.1:win(2); onoff = {'onsets','offsets'};
for j = 1:numel(Fl), [~,~,~,Df] = far_field(B(li(j)), Fl(j), snr_frac, 2, 90); Fl(j).atr0 = anti_trough(Df, Fl(j).ph, Fl(j).fr, 135, 0, 0); end
TT = {Fd, 'pk', 'GRAB peak'; Fl, 'atr0', 'LPsP trough'; Fl, 'pk', 'LPsP peak'; Fe, 'pk', 'EPG7f peak'}; TR = cell(size(TT,1),2);
for q = 1:size(TT,1)
    FF = TT{q,1};
    for d_ = 0:1
        R_ = struct('e',{},'r',{},'y',{},'m',{},'sp',{}); ne = [0 0]; fq = [];
        for j = find([FF.dark] == d_)
            t = FF(j); [E_, nev, f_] = turn_transients(t.(TT{q,2}), t.vr, t.fr, t.bad, win); if sum(nev) == 0, continue, end
            ne = ne + nev; fq(end+1,:) = f_; ip = @(z) interp1(E_.t, z, tg, 'linear'); %#ok<AGROW>
            for e = 1:2, if nev(e) > 0, R_(end+1) = struct('e',e,'r',ip(E_.r{e}),'y',ip(E_.y{e}),'m',ip(E_.m{e}),'sp',ip(E_.sp{e})); end, end %#ok<AGROW>
        end
        TR{q,d_+1} = struct('R', R_, 'ne', ne, 'fq', fq);
        for e = 1:2
            rr = vertcat(R_([R_.e]==e).r); ea = mean(rr(:, tg > 0 & tg <= 1.5),2,'omitnan');
            fprintf('%-11s %-11s %-7s (%3d events, %2d trials; model R2 %.2f): residual 0-1.5 s %+.2f SD (>0 %d/%d, p=%.2g)\n', TT{q,3}, ifelse_c(d_), onoff{e}, ne(e), size(rr,1), median(fq(:,3)), median(ea), sum(ea>0), numel(ea), signrank(ea));
        end
    end
end
% figure
figure('Color','w','Position',[60 60 1400 820]); tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
cmap = [0.78 0.85 0.95; 0.55 0.70 0.90; 0.30 0.52 0.82; 0.12 0.33 0.68; 0.03 0.15 0.45]; xd = [-180:22.5:157.5 180];
blab = {'0-0.25','0.25-0.5','0.5-1','1-2','> 2 rad/s'}; ttl14 = {'GRAB(DA)','LPsP syt7f','EPG GCaMP7f'};
for sgi = 1:3
    ax = nexttile(tl, sgi); hold(ax,'on'); M = mean(CR{sgi}.P,3,'omitnan');
    for k = 1:5, plot(ax, xd, [M(:,k); M(1,k)], '-', 'Color', cmap(k,:), 'LineWidth', 1.8); end
    xline(ax, 0, ':'); xlim(ax, [-180 180]); set(ax,'XTick',-180:90:180,'Box','off','TickDir','out');
    xlabel(ax, 'relative to cue-predicted bump (deg; + = ahead)'); if sgi == 1, ylabel(ax, 'dF/F (closed loop)'); end
    title(ax, ttl14{sgi}, 'FontWeight','normal'); if sgi == 3, legend(ax, blab, 'Box','off','FontSize',8,'Location','northeast'); end
end
ax = nexttile(tl, 4); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]);
cS = [0.85 0.4 0.1; 0.2 0.6 0.3; 0.5 0.5 0.5];
for sgi = 1:3, plot(ax, vc, median(CR{sgi}.C,1,'omitnan'), 'o-', 'Color', [0 0 0]*(sgi==1) + cS(sgi,:)*(sgi>1), 'LineWidth', 2, 'MarkerFaceColor', 'w'); end
ls_ = {':','--','-.','-'};
for it = [2 3], plot(ax, vc, SIM(it,:,1), ls_{it}, 'Color', [0 0.45 0.74], 'LineWidth', 1.5); end
for it = [2 4], plot(ax, vc, SIM(it,:,2), ls_{it}, 'Color', [0.75 0.2 0.2], 'LineWidth', 1.5); end
set(ax,'XScale','log','XTick',vc,'XTickLabel',{'0.12','0.37','0.75','1.5','3'},'Box','off','TickDir','out'); xlim(ax,[0.1 3.5]);
xlabel(ax,'bump (cue) speed (rad/s)'); ylabel(ax,'profile centroid (deg; + = ahead)'); title(ax,'where the signal sits','FontWeight','normal');
legend(ax, {'','GRAB','LPsP','EPG7f',sprintf('product, \\tau %.1f s',taus(2)),sprintf('product, \\tau %.1f s',taus(3)),sprintf('[dL/dt]_+, \\tau %.1f s',taus(2)),sprintf('[dL/dt]_+, \\tau %.1f s',taus(4))}, 'Box','off','FontSize',7,'Location','south','NumColumns',2); ylim(ax,[-75 75]);
qs = [1 2]; lab = {'GRAB peak','LPsP trough'};
for qq = 1:2, for d_ = 0:1
    ax = nexttile(tl, 4 + (qq-1)*2 + d_ + 1); hold(ax,'on'); R_ = TR{qs(qq),d_+1}.R;
    for e = 1:2
        cc_ = [0 0 0; 0.5 0.5 0.5]; yy = vertcat(R_([R_.e]==e).y); mm = vertcat(R_([R_.e]==e).m);
        plot(ax, tg, mean(yy,1,'omitnan'), '-', 'Color', cc_(e,:), 'LineWidth', 2); plot(ax, tg, mean(mm,1,'omitnan'), '--', 'Color', cc_(e,:), 'LineWidth', 1.2);
    end
    xline(ax, 0, ':'); xlim(ax, win); set(ax,'Box','off','TickDir','out'); xlabel(ax,'time from turn onset / offset (s)');
    if qq == 1 && d_ == 0, ylabel(ax,'signal (SD units)'); end
    title(ax, sprintf('%s, %s', lab{qq}, ifelse_c(d_)), 'FontWeight','normal');
    if qq == 1 && d_ == 0, legend(ax, {'onset: data','onset: low-passed |v_r| model','offset: data','offset: model'}, 'Box','off','FontSize',7,'Location','northwest'); end
end, end
title(tl, 'GRAB as motion energy of LPsP calcium? No leading edge, no onset transient: a centred bump delayed by sensor kinetics fits', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_motion_energy.png'),'Resolution',150);

%% (15) depletion test: is GRAB lower where the bump has recently been? (figure)
% Synaptic depression ("new glomeruli") predicts that, at matched speed,
% release is lower when the bump sits where it has recently been and higher
% when it arrives somewhere it hasn't been for a while. Closed loop only:
% the bump position comes from the cue (ref = cue-predicted position, lag
% from EPG7f), so the bump history is known even when GRAB is dark and is
% independent of GRAB's own phase noise. (No equivalent reference in the
% dark.) Per trial, GRAB peak (t) regressed on fly variables at t - 0.45 s
% and bump variables at t - 0.3 s:
%  basic:    |vr|, net rotation, forward speed, |bump speed|, occupancy
%  strict:   cubic terms of each speed, mean bump speed over the past 1/2/4
%            s, reversal covariates (fraction of the past 2 s turning the
%            other way, # reversals), first 30 s of each trial excluded
%  placebo:  strict + occupancy of the location OPPOSITE the bump
% Occupancy = exponentially weighted (tau 1-8 s) closeness of past bump
% positions to the current one (c(45 deg) = 0.5). Recovery curve: GRAB
% residual after the strict controls (no occupancy), bump moving > 0.5
% rad/s, by time since the bump last overlapped its current location
% (after being > 45 deg away).
taus = [1 2 4 8]; trials = find(~[Fd.dark] & RkF >= 0.4); HG = cell(1,numel(Fd)); Lfly = 0.45; cap = 30; tb = [0 1 2 4 8 16 30];
for j = 1:numel(Fd), s = B(gi(j)); Fd(j).fw = interp1(s.ft.xf(:), smoothdata(s.ft.f_speed(:),'gaussian',30), s.ft.xb(:)', 'linear', 'extrap'); end
for j = trials
    t = Fd(j); [~,~,~,ref,v] = cue_ref_profiles(B(gi(j)), t, snr_frac, lagE, edges, zeros(16,numel(t.ph)));
    H = bump_history(ref, t.fr, taus); H.ref = ref; H.v = v; HG{j} = H;
end
pb = @(z) [z(:) z(:).^2 z(:).^3]; kk = log(2)/(1 - cos(pi/4));
CB = nan(numel(trials), numel(taus)); CS = CB; CPa = CB; CPo = CB; RC = nan(numel(trials), numel(tb));
for jj = 1:numel(trials)
    j = trials(jj); t = Fd(j); H = HG{j}; fr = t.fr; n = numel(t.ph); y = t.pk; Lf = round(Lfly*fr); Lb = round(0.3*fr); sh = @(z,L) [nan(1,L) z(1:n-L)];
    sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); net = abs(smoothdata(t.vr,'gaussian',max(3,round(fr)))); fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bsp = abs(H.v);
    vs = smoothdata(H.v,'gaussian',max(3,round(0.3*fr))); cur = sign(vs); opp = nan(1,n); nrev = nan(1,n); W = round(2*fr);
    for i = W+1:n, seg = vs(i-W:i-1); opp(i) = mean(sign(seg) == -cur(i) & abs(seg) > 0.3); nrev(i) = sum(abs(diff(sign(seg(abs(seg) > 0.3)))) > 0); end
    Cb = [sh(sp,Lf); sh(net,Lf); sh(fw,Lf); sh(bsp,Lb)]';
    Cs = [pb(sh(sp,Lf)) pb(sh(net,Lf)) pb(sh(fw,Lf)) pb(sh(bsp,Lb)) sh(opp,Lb)' sh(nrev,Lb)'];
    for wsec = [1 2 4], Cs = [Cs movmean(sh(bsp,Lb), [round(wsec*fr) 0])']; end %#ok<AGROW>
    for wsec = [1 2 4], h_ = movmean(sh(sp,Lf), [round(wsec*fr) 0])'; Cs = [Cs h_ h_.^2]; end %#ok<AGROW>   % fly-speed history (added later: see header)
    early = (1:n)' <= round(30*fr); refi = fillmissing(H.ref,'nearest');
    for it = 1:numel(taus)
        oc = sh(H.occ(it,:),Lb)';
        X = [Cb oc]; m = all(isfinite([X y']),2) & ~t.bad'; Z = zscore([X(m,:) y(m)']); b = [ones(sum(m),1) Z(:,1:end-1)] \ Z(:,end); CB(jj,it) = b(end);
        X = [Cs oc]; m = all(isfinite([X y']),2) & ~t.bad' & ~early; Z = zscore([X(m,:) y(m)']); b = [ones(sum(m),1) Z(:,1:end-1)] \ Z(:,end); CS(jj,it) = b(end);
        K = round(5*taus(it)*fr); w = exp(-(1:K)/(taus(it)*fr)); w = w/sum(w); oO = zeros(1,n);
        for q = 1:K, past = [nan(1,q) refi(1:n-q)]; cO = exp(kk*(cos(past - refi - pi) - 1)); cO(isnan(cO)) = 0; oO = oO + w(q)*cO; end
        oO(1:K) = nan; X = [Cs oc sh(oO,Lb)']; m = all(isfinite([X y']),2) & ~t.bad' & ~early;
        Z = zscore([X(m,:) y(m)']); b = [ones(sum(m),1) Z(:,1:end-1)] \ Z(:,end); CPa(jj,it) = b(end-1); CPo(jj,it) = b(end);
    end
    m = all(isfinite([Cs y']),2) & ~t.bad' & ~early; b = [ones(sum(m),1) Cs(m,:)] \ y(m)'; r = nan(1,n); r(m) = y(m) - ([ones(sum(m),1) Cs(m,:)]*b)'; r = r/std(r(m));
    tsl = sh(H.tsl,Lb); tsl(isnan(tsl) & ~isnan(r)) = cap + 1; g = ~isnan(r) & sh(bsp,Lb) > 0.5;
    for k = 1:numel(tb)
        if k < numel(tb), mm = g & tsl >= tb(k) & tsl < tb(k+1); else, mm = g & tsl > cap; end
        if sum(mm) >= 20, RC(jj,k) = mean(r(mm)); end
    end
end
lab15 = {'basic','strict','strict: own location','strict: OPPOSITE location'}; CC15 = {CB, CS, CPa, CPo};
for c = 1:4, fprintf('%-26s occupancy coef (tau 1/2/4/8 s) %s | >0 in %s of %d | p %s\n', lab15{c}, mat2str(round(median(CC15{c},1),3)), mat2str(sum(CC15{c}>0,1)), numel(trials), mat2str(arrayfun(@(k) signrank(CC15{c}(:,k)), 1:4),2)); end
fprintf('recovery curve (moving bump), GRAB residual by time since last visit [0-1 1-2 2-4 4-8 8-16 16-30 >30 s]: %s (trials %s)\n', mat2str(round(median(RC,1,'omitnan'),3)), mat2str(sum(~isnan(RC),1)));
d_ = mean(RC(:,5:7),2,'omitnan') - mean(RC(:,1:2),2,'omitnan'); fprintf('fresh (> 8 s) minus revisit (< 2 s): %+.3f SD (>0 %d/%d, p=%.2g)\n', median(d_,'omitnan'), sum(d_>0), sum(~isnan(d_)), signrank(d_));
% figure
figure('Color','w','Position',[80 80 1250 430]); tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl, [1 2]); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off'); cc15 = [0.5 0.5 0.5; 0 0 0; 0 0.45 0.74; 0.85 0.4 0.1];
for c = 1:4
    for it = 1:numel(taus)
        x0 = it + (c - 2.5)*0.18; v_ = CC15{c}(:,it);
        scatter(ax, x0 + 0.03*randn(size(v_)), v_, 14, cc15(c,:), 'filled', 'MarkerFaceAlpha', 0.35, 'HandleVisibility','off');
        plot(ax, x0 + [-0.07 0.07], median(v_)*[1 1], '-', 'Color', cc15(c,:), 'LineWidth', 3, 'HandleVisibility', ifelse_on(it == 1));
    end
end
set(ax,'XTick',1:numel(taus),'XTickLabel',arrayfun(@(x) sprintf('%g s',x), taus, 'UniformOutput', false),'Box','off','TickDir','out'); xlim(ax,[0.4 numel(taus)+0.6]);
xlabel(ax,'occupancy time window \tau'); ylabel(ax,'partial coefficient on GRAB peak (z)');
legend(ax, lab15, 'Box','off','Location','southeast','FontSize',8); title(ax,'recent occupancy of the bump''s location lowers GRAB (closed loop, dots = trials)','FontWeight','normal');
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]); xc = [0.5 1.5 3 6 12 23 40];
for jj = 1:size(RC,1), plot(ax, xc, RC(jj,:), '-', 'Color', [0.8 0.8 0.8]); end
errorbar(ax, xc, median(RC,1,'omitnan'), std(RC,0,1,'omitnan')./sqrt(sum(~isnan(RC),1)), 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w');
set(ax,'XScale','log','XTick',xc(2:end),'XTickLabel',{'1-2','2-4','4-8','8-16','16-30','>30'},'Box','off','TickDir','out'); xlim(ax,[1 50]); ylim(ax,[-0.6 0.7]);
xlabel(ax,'time since the bump last visited its current location (s)'); ylabel(ax,'GRAB peak residual (SD)'); title(ax,'recovery after a visit (moving bump)','FontWeight','normal');
title(tl, 'Depletion test: GRAB vs the bump''s recent history at its current location (speed, reversals, recent speed controlled)', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_depletion.png'),'Resolution',150);

%% (16) modelling GRAB from LPsP: product vs per-glomerulus release (power, depletion) (figure)
% LPsP and GRAB are in different flies, so (1) an input model of LPsP
% calcium is fitted in LPsP closed-loop flies:
%   F/F0(theta,t) = c0 + gx x(t) + gb xb(t) + k T_L(theta - ref(t)),
% x / xb = |fly| / |cue-predicted bump| speed, low-passed (tau, lag on a
% grid); ref = cue-predicted bump position (EPG7f lag, section 14);
% (2) a synthetic LPsP Lhat is built for each GRAB closed-loop fly from its
% own behaviour and cue trajectory; (3) release models are fitted to the
% full GRAB(theta,t), each through a GRAB sensor kernel (exp, tauS):
%   product: G(t) (1 + k2 T(theta - ref)), G = gx x (ball) or gx x + gb xb
%   power: Lhat^n        depletion: Lhat^n R, dR/dt = (1-R)/tauRec - u Lhat^n R
%   product variants: free bump-speed weight w, depletion acting on the
%   product's drive, both.
% Shared parameters by fminsearch; per-trial offset and scale by least
% squares; 2-fold CV over trial halves (fit on one half, test the other).
% Needs sections 14 (lagE, CR) and 15 not required.
M = mean(CR{2}.P,3,'omitnan'); tpL = mean(M,2) - min(mean(M,2)); tpL = (tpL + tpL(mod(17-(1:16),16)+1))/2; [~,i0] = max(tpL); tpL = circshift(tpL, 9-i0); tpL = tpL/max(tpL);
trL = find(~[Fl.dark] & RkL >= 0.4); tgrid = [0.3 0; 0.6 0; 1.0 0; 1.5 0; 2.5 0; 0.6 0.2; 1.0 0.2; 1.5 0.2];
[PL_, ~, flog] = lpsp_input_fit(B(li), Fl, trL, snr_frac, lagE, edges, tpL, tgrid);
fprintf('LPsP input model (F/F0): %.2f + %.3f |vr| + %.3f |bump v| + %.3f bump(theta) | tau %.1f s lag %.1f s | median R2 %.2f (%d trials)\n', PL_.c0, PL_.gx, PL_.gb, PL_.k, PL_.tauG, PL_.lagG, median(PL_.r2), numel(trL));
tmplL = @(th) interp1([a16; pi], [tpL; tpL(1)], mod(th + pi, 2*pi) - pi, 'linear');
trG = find(~[Fd.dark] & RkF >= 0.4); DAT = struct('Y',{},'L',{},'G',{},'E',{},'fr',{},'ok',{},'x',{},'xb',{}); DATb = DAT;
for jj = 1:numel(trG)
    j = trG(jj); t = Fd(j); s = B(gi(j)); n = numel(t.ph);
    [~,~,~,Df] = far_field(s, t, snr_frac, 2, 90); [~,~,~,ref,v] = cue_ref_profiles(s, t, snr_frac, lagE, edges, zeros(16,n));
    [x, xb] = input_drives(t.vr, v, t.fr, PL_.tauG, PL_.lagG); E = tmplL(a16 - fillmissing(ref,'nearest')); G = PL_.gx*x + PL_.gb*xb;
    ok = ~t.bad & ~isnan(ref) & ~isnan(x);
    DAT(jj) = struct('Y', Df, 'L', PL_.c0 + repmat(G,16,1) + PL_.k*E, 'G', repmat(G,16,1), 'E', E, 'fr', t.fr, 'ok', ok, 'x', x, 'xb', xb);
    DATb(jj) = DAT(jj); DATb(jj).G = repmat(PL_.gx*x,16,1);
end
opts = optimset('Display','off','MaxFunEvals',600,'MaxIter',600,'TolX',1e-3,'TolFun',1e-5);
MOD = {'product (ball gain)', 'product', DATb, log([0.4 3]); 'product (LPsP gain, incl. bump speed)', 'product', DAT, log([0.4 3]); ...
       'power L^n', 'power', DAT, log([0.4 4]); 'depletion L^n R', 'depletion', DAT, log([0.4 4 1.5 0.05]); ...
       'product, free bump-speed w', 'product_w', DAT, [log(0.3) log(3) 0]; 'product + depletion', 'product_dep', DATb, [log(0.3) log(3) log(1.5) log(0.5)]; ...
       'product, free w + depletion', 'product_wdep', DAT, [log(0.3) log(3) 0 log(1.5) log(0.5)]};
FIT = struct('cv',{},'thA',{});
for mi = 1:size(MOD,1)
    th1 = fminsearch(@(th) grab_model_eval(th, MOD{mi,2}, MOD{mi,3}, 1), MOD{mi,4}, opts);
    th2 = fminsearch(@(th) grab_model_eval(th, MOD{mi,2}, MOD{mi,3}, 2), MOD{mi,4}, opts);
    [~,o12] = grab_model_eval(th1, MOD{mi,2}, MOD{mi,3}, 2); [~,o21] = grab_model_eval(th2, MOD{mi,2}, MOD{mi,3}, 1);
    thA = fminsearch(@(th) grab_model_eval(th, MOD{mi,2}, MOD{mi,3}, 0), (th1+th2)/2, opts);
    FIT(mi).cv = (o12.r2 + o21.r2)/2; FIT(mi).thA = thA; d = FIT(mi).cv - FIT(1).cv*(mi > 1);
    fprintf('%-38s CV R2 %.3f | vs ball-gain product %+.4f (better %d/%d) | params (all data, log except w) %s\n', MOD{mi,1}, median(FIT(mi).cv), median(d)*(mi > 1), sum(d>0)*(mi > 1), numel(d), mat2str(round(thA,2)));
end
RELC = nan(numel(trG), 5, 4); mi_ = [1 3 4];
for q = 1:3
    [~, o] = grab_model_eval(FIT(mi_(q)).thA, MOD{mi_(q),2}, MOD{mi_(q),3}, 0);
    for jj = 1:numel(trG), P = cue_ref_profiles(B(gi(trG(jj))), Fd(trG(jj)), snr_frac, lagE, edges, o.pred{jj}); RELC(jj,:,q) = (max(P,[],1) - min(P,[],1))./min(P,[],1); end
end
for jj = 1:numel(trG), P = CR{1}.P(:,:,jj); RELC(jj,:,4) = (max(P,[],1) - min(P,[],1))./min(P,[],1); end
PL2 = CR{2}.P; relL = squeeze((max(PL2,[],1) - min(PL2,[],1))./min(PL2,[],1))';
nm_ = {'product model','power model','depletion model','GRAB data'};
for q = [4 1 2 3], fprintf('%-16s bump height / trough by speed bin: %s\n', nm_{q}, mat2str(round(median(RELC(:,:,q),1,'omitnan'),2))); end
fprintf('%-16s bump height / trough by speed bin: %s\n', 'LPsP data', mat2str(round(median(relL,1,'omitnan'),2)));
figure('Color','w','Position',[80 80 1250 460]); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); mlab = {'product (ball gain)','product (LPsP gain, incl. bump speed)','power L^n','depletion L^n R','product, free bump-speed w','product + depletion','product, free w + depletion'};
ord = [1 2 5 6 7 3 4]; CVm = cell2mat(arrayfun(@(k) FIT(k).cv(:), ord, 'UniformOutput', false));
plot(ax, 1:numel(ord), CVm', '-', 'Color', [0.85 0.85 0.85]);
for q = 1:numel(ord)
    c_ = [0 0 0]; if ord(q) >= 3 && ord(q) <= 4, c_ = [0.75 0.2 0.2]; end
    scatter(ax, q*ones(size(CVm,1),1), CVm(:,q), 18, c_, 'filled', 'MarkerFaceAlpha', 0.5);
    plot(ax, q + [-0.25 0.25], median(CVm(:,q))*[1 1], '-', 'Color', c_, 'LineWidth', 3);
end
set(ax,'XTick',1:numel(ord),'XTickLabel',mlab(ord),'XTickLabelRotation',30,'Box','off','TickDir','out'); xlim(ax,[0.5 numel(ord)+0.5]);
ylabel(ax,'cross-validated R^2 of GRAB(\theta,t)'); title(ax,'GRAB from synthetic LPsP (closed loop, lines = trials)','FontWeight','normal');
ax = nexttile(tl); hold(ax,'on'); vc = [0.12 0.37 0.75 1.5 3];
cc16 = {[0 0 0], [0.2 0.6 0.3], [0 0.45 0.74], [0.75 0.2 0.2], [0.95 0.55 0.55]};
plot(ax, vc, median(RELC(:,:,4),1,'omitnan'), 'o-', 'Color', cc16{1}, 'LineWidth', 2.5, 'MarkerFaceColor', 'w');
plot(ax, vc, median(relL,1,'omitnan'), 'o-', 'Color', cc16{2}, 'LineWidth', 2.5, 'MarkerFaceColor', 'w');
plot(ax, vc, median(RELC(:,:,1),1,'omitnan'), '--', 'Color', cc16{3}, 'LineWidth', 2);
plot(ax, vc, median(RELC(:,:,2),1,'omitnan'), '--', 'Color', cc16{4}, 'LineWidth', 2);
plot(ax, vc, median(RELC(:,:,3),1,'omitnan'), ':', 'Color', cc16{5}, 'LineWidth', 2);
set(ax,'XScale','log','XTick',vc,'XTickLabel',{'0.12','0.37','0.75','1.5','3'},'Box','off','TickDir','out'); xlim(ax,[0.1 3.5]);
xlabel(ax,'bump (cue) speed (rad/s)'); ylabel(ax,'bump height / trough'); title(ax,'relative bump contrast vs speed','FontWeight','normal');
legend(ax, {'GRAB data','LPsP data','product model','power model','depletion model'}, 'Box','off','Location','northwest','FontSize',8);
title(tl, 'Modelling GRAB from LPsP: a per-glomerulus transform (power, depletion) can''t make contrast rise with speed; the global x local product can', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_release_models.png'),'Resolution',150);

%% (17) EPG-gated "motion energy" of LPsP calcium: an adapting synapse per glomerulus (figure)
% Instead of gain = global LPsP calcium, release = a per-glomerulus temporal-
% change signal of (synthetic) LPsP calcium, still gated by the EPG bump:
%   Rel = [dL - beta * LP_tauA(dL)]+ (1 + k2 T(theta - ref)),  GRAB = sensor * Rel
% dL = LPsP F/F0 minus rest (section 16 input model), LP = causal exponential
% low-pass per glomerulus. beta = 0: ~ LPsP x EPG; beta = 1: pure temporal
% change ("motion energy"), the EPG gate confining any global transient to
% the bump. Adaptation lowers release where LPsP has been high (the bump's
% recent location), so it is also a depletion-like mechanism. Control: only
% the PB-wide part adapts. 2-fold CV over trial halves as in section 16.
% Then the section-15 occupancy test (now with fly-speed history controls,
% see header) on the data, on each model's predicted peak, and on the
% data - model residual: a model that captures the depletion effect should
% reproduce it in its prediction and remove it from the residual.
for jj = 1:numel(DAT)
    j = trG(jj); [~,~,~,ref,v] = cue_ref_profiles(B(gi(j)), Fd(j), snr_frac, lagE, edges, zeros(16,numel(Fd(j).ph)));
    DAT(jj).dL = DAT(jj).L - PL_.c0; DAT(jj).ref = ref; DAT(jj).v = v; DATb(jj).dL = DAT(jj).dL;
end
MODg = {'product (ball gain)', 'product', DATb, log([0.3 3]), []; 'gated LPsP (beta = 0)', 'gated_L', DAT, log([0.3 3]), []; ...
        'gated adapting LPsP', 'gated_adapt', DAT, [log(0.3) log(3) 0 log(2)], []; 'only PB-wide part adapts', 'gated_adapt_glob', DAT, [log(0.3) log(3) 0 log(2)], []; ...
        'gated adapting, tauA 1 s', 'gated_adapt', DAT, [log(0.3) log(3) 0], log(1); 'gated adapting, tauA 2 s', 'gated_adapt', DAT, [log(0.3) log(3) 0], log(2); ...
        'gated motion energy (beta = 1)', 'gated_me', DAT, [log(0.3) log(3) log(2)], []};
FG = struct('cv',{},'thA',{},'th',{});
for mi = 1:size(MODg,1)
    fx = MODg{mi,5}; full_ = @(th) [th fx];
    cst = @(th, h) grab_model_eval(full_(th), MODg{mi,2}, MODg{mi,3}, h);
    th1 = fminsearch(@(th) cst(th,1), MODg{mi,4}, opts); th2 = fminsearch(@(th) cst(th,2), MODg{mi,4}, opts);
    [~,o12] = grab_model_eval(full_(th1), MODg{mi,2}, MODg{mi,3}, 2); [~,o21] = grab_model_eval(full_(th2), MODg{mi,2}, MODg{mi,3}, 1);
    thA = fminsearch(@(th) cst(th,0), (th1+th2)/2, opts);
    FG(mi).cv = (o12.r2 + o21.r2)/2; FG(mi).th = full_(thA); d = FG(mi).cv - FG(1).cv;
    p_ = exp(FG(mi).th); if any(strcmp(MODg{mi,2}, {'gated_adapt','gated_adapt_glob'})), p_(3) = 1/(1+exp(-FG(mi).th(3))); end
    fprintf('%-32s CV R2 %.3f | vs product %+.4f (better %d/%d, p=%.2g) | params [tauS k2 (beta) tauA] %s\n', MODg{mi,1}, median(FG(mi).cv), median(d), sum(d>0), numel(d), signrank(d+eps*(mi==1)), mat2str(round(p_,2)));
end
d = FG(3).cv - FG(4).cv; fprintf('per-glomerulus minus PB-wide-only adaptation: %+.4f (better %d/%d, p=%.2g)\n', median(d), sum(d>0), numel(d), signrank(d));
taus17 = [1 2 4]; mq = [1 3 5]; [~, oP] = grab_model_eval(FG(1).th, 'product', DATb, 0); PRD = cell(1,3); PRD{1} = oP.pred;
for q = 2:3, [~, o_] = grab_model_eval(FG(mq(q)).th, MODg{mq(q),2}, DAT, 0); PRD{q} = o_.pred; end
OC17 = nan(numel(trG), numel(taus17), 6);
for jj = 1:numel(trG)
    j = trG(jj); y = max(DAT(jj).Y,[],1); pp = cellfun(@(P) max(P{jj},[],1), PRD, 'UniformOutput', false);
    ys = {y, pp{1}, pp{2}, pp{3}, y - pp{1}, y - pp{2}};
    for q = 1:6, OC17(jj,:,q) = occ_coef(ys{q}, Fd(j), DAT(jj).ref, DAT(jj).v, taus17, true, false); end
end
lab17 = {'GRAB data','product prediction','gated adapting prediction','gated adapting (tauA 1 s) prediction','residual: data - product','residual: data - gated adapting'};
for q = 1:6, fprintf('%-38s occupancy coef (tau 1/2/4 s) %s | <0 in %s of %d\n', lab17{q}, mat2str(round(median(OC17(:,:,q),1,'omitnan'),3)), mat2str(sum(OC17(:,:,q)<0,1)), numel(trG)); end
% figure
figure('Color','w','Position',[80 80 1250 470]); tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]); ord = 2:size(MODg,1);
D17 = cell2mat(arrayfun(@(k) FG(k).cv(:) - FG(1).cv(:), ord, 'UniformOutput', false));
for q = 1:numel(ord)
    c_ = [0 0.45 0.74]; if ord(q) == 4, c_ = [0.5 0.5 0.5]; elseif ord(q) == 7, c_ = [0.75 0.2 0.2]; end
    scatter(ax, q + 0.06*randn(size(D17,1),1), D17(:,q), 18, c_, 'filled', 'MarkerFaceAlpha', 0.5);
    plot(ax, q + [-0.25 0.25], median(D17(:,q))*[1 1], '-', 'Color', c_, 'LineWidth', 3);
end
set(ax,'XTick',1:numel(ord),'XTickLabel',MODg(ord,1),'XTickLabelRotation',30,'Box','off','TickDir','out'); xlim(ax,[0.5 numel(ord)+0.5]);
ylabel(ax,'\Delta cross-validated R^2 vs product (per trial)'); title(ax,'GRAB(\theta,t) from EPG-gated LPsP change signals (closed loop)','FontWeight','normal');
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]); qs = [1 2 3 5 6]; cc17 = {[0 0 0],[0 0.45 0.74],[0.85 0.4 0.1],[0 0.45 0.74],[0.85 0.4 0.1]};
for q = 1:numel(qs)
    v_ = OC17(:,1,qs(q)); scatter(ax, q + 0.06*randn(size(v_)), v_, 18, cc17{q}, 'filled', 'MarkerFaceAlpha', 0.5);
    plot(ax, q + [-0.25 0.25], median(v_)*[1 1], '-', 'Color', cc17{q}, 'LineWidth', 3);
end
xline(ax, 3.5, ':'); set(ax,'XTick',1:numel(qs),'XTickLabel',{'data','product pred.','adapting pred.','data - product','data - adapting'},'XTickLabelRotation',30,'Box','off','TickDir','out'); xlim(ax,[0.5 numel(qs)+0.5]);
ylabel(ax,'occupancy coefficient on GRAB peak (\tau 1 s)'); title(ax,'does the model capture the depletion-like occupancy effect?','FontWeight','normal');
title(tl, 'EPG gate x adapting LPsP: a small, glomerulus-specific gain in fit, but the occupancy effect stays in the residual', 'FontSize', 11);
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_gated_adapting.png'),'Resolution',150);

%% figure
c = [0.30 0.30 0.30; 0.00 0.45 0.74; 0.85 0.33 0.10];
figure(1); clf; set(gcf,'Color','w','Position',[100 100 1400 400]); tl = tiledlayout(1,3,'TileSpacing','compact');
nexttile; hold on
errorbar(wins, median(Vep,1), iqr(Vep)/2, 'o-', 'Color', c(1,:), 'MarkerFaceColor', c(1,:), 'LineWidth', 1.5)
errorbar(wins, median(Vc(clT,:),1), iqr(Vc(clT,:))/2, 's-', 'Color', c(2,:), 'MarkerFaceColor', c(2,:), 'LineWidth', 1.5)
plot(wins, median(Rel(clT,:),1), 's:', 'Color', c(2,:), 'LineWidth', 1.2)
set(gca,'XScale','log','XTick',wins); ylim([0 1]); xlabel('velocity smoothing / window (s)'); ylabel('|r|')
legend({'EPG GCaMP7f bump vel vs cue vel', sprintf('GRAB bump vel vs cue vel (%d CL trials)', sum(clT)), 'GRAB: half 1 vs half 2 displacement'}, 'Box','off','Location','northwest')
title('A  how well GRAB resolves bump motion','FontWeight','normal')
nexttile; hold on
for d_ = 0:1
    plot(winS, squeeze(JW(d_+1,:,4)), 'o-', 'Color', c(1,:)*(1-0.4*d_), 'MarkerFaceColor', c(1,:)*(1-0.4*d_), 'LineWidth', 1.5, 'LineStyle', ifelse_s(d_))
    plot(winS, squeeze(JW(d_+1,:,5)), 's-', 'Color', c(2,:), 'MarkerFaceColor', c(2,:), 'LineWidth', 1.5, 'LineStyle', ifelse_s(d_))
    plot(winS, squeeze(JW(d_+1,:,6)), 'd', 'Color', c(3,:), 'LineWidth', 1.2, 'LineStyle', ifelse_s(d_))
end
yline(0,':'); set(gca,'XTick',winS); xlabel('window (s)'); ylabel('standardized coefficient, joint model'); ylim([-0.2 0.9])
legend({'fly speed (CL)','bump speed (CL)','bump speed, cross-half (CL)','fly speed (dark)','bump speed (dark)','bump speed, cross-half (dark)'},'Box','off','Location','east')
title('B  GRAB peak ~ fly speed + bump speed (within trials)','FontWeight','normal')
nexttile
scatter(TR(dk,1), TR(dk,2), 40, c(2,:), 'filled'); lsline
xlabel('bump gain (bump displacement / fly rotation), dark trial'); ylabel('GRAB-vs-|fly speed| slope / GRAB range')
title(sprintf('C  dark trials: Spearman %.2f (p = %.2g)', rN, pN),'FontWeight','normal')
title(tl, {'Does the GRAB peak follow the fly''s rotational speed or the bump''s? (EPG GRAB, bump from GRAB itself)', ...
    'GRAB tracks fly speed; bump motion beyond the fly''s turning adds nothing -- but per-window bump speed is noisy (A), so only a large bump effect is excluded'})
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed.png'),'Resolution',150);

figure(2); clf; set(gcf,'Color','w','Position',[100 550 1100 380]); tl2 = tiledlayout(1,2,'TileSpacing','compact');
for d_ = 0:1
    nexttile; hold on
    plot(lagsS, median(LS(d_+1).CF,1), '-', 'Color', c(1,:), 'LineWidth', 2)
    plot(lagsS, median(LS(d_+1).CB,1), '-', 'Color', c(2,:), 'LineWidth', 1.5)
    plot(lagsS, median(LS(d_+1).CQ,1), '--', 'Color', c(2,:), 'LineWidth', 1.5)
    plot(lagsS, median(LS(d_+1).JB,1), ':', 'Color', c(3,:), 'LineWidth', 2)
    xline(0,':'); yline(0,':'); ylim([-0.2 0.9]); xlabel('lag of predictor before GRAB peak (s)'); ylabel('correlation / coefficient with GRAB peak')
    title(ifelse_c(d_),'FontWeight','normal')
    if d_ == 0
        legend({'fly speed','GRAB bump speed (all frames)','bump speed, clear-bump frames (cross-half)','bump coefficient with fly speed in the model'},'Box','off','Location','northwest')
    end
end
title(tl2, 'Lag scan: no lag makes GRAB bump speed explain the GRAB peak beyond fly speed')
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_lagscan.png'),'Resolution',150);

figure(3); clf; set(gcf,'Color','w','Position',[100 950 1100 380]); tl3 = tiledlayout(1,2,'TileSpacing','compact');
for d_ = 0:1
    nexttile; hold on
    plot(lagsS, FO(d_+1).CF, '-', 'Color', c(1,:), 'LineWidth', 2)
    plot(lagsS, FO(d_+1).CB, '-', 'Color', c(2,:), 'LineWidth', 2)
    plot(lagsS, FO(d_+1).JF, '--', 'Color', c(1,:), 'LineWidth', 1.3)
    plot(lagsS, FO(d_+1).JB, '--', 'Color', c(2,:), 'LineWidth', 1.3)
    xline(0,':'); yline(0,':'); ylim([-0.2 0.8]); xlabel('lag of predictor before GRAB peak (s)'); ylabel('correlation / coefficient with GRAB peak')
    title(ifelse_c(d_),'FontWeight','normal')
    if d_ == 0, legend({'fly speed (alone)','bump speed (alone)','fly speed (joint)','bump speed (joint)'},'Box','off','Location','northwest'); end
end
title(tl3, sprintf('Folded halves, well-observed bumps only (rho > %.1f): bump speed adds a modest share, fly speed stays stronger', rthr))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_folded_rho.png'),'Resolution',150);

if ~isempty(EXP)
    figure(4); clf; set(gcf,'Color','w','Position',[50 50 1300 1250]); tl4 = tiledlayout(6,2,'TileSpacing','compact');
    for r = 1:size(EXP,1), plot_pair(tl4, r, B(gi(EXP(r,1))), Fd(EXP(r,1)), EXP(r,:), ifelse_c(EXP(r,2)), [0.45 0.45 0.45; 0.85 0.33 0.10], wsP); end
    title(tl4, {'Matched fly rotation, different bump rotation (EPG GRAB; bump from folded GRAB, rho > 0.4)', ...
        'each row: two stretches from the SAME trial (same fly); shared axes and colour scale within a row', ...
        'white line: folded bump phase (half-1 coordinates); dashed: 1.5 s window; speeds shifted by their best lags to align with GRAB'})
    fontsize(gcf,11,'pixels')
    exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_peak_fly_vs_bump_speed_matched_pairs.png'),'Resolution',150);
end

%% functions
function plot_pair(tl, rowIdx, s, t, pr, ttl, cols, wsP)
% two excerpts (slower / faster bump) from one GRAB trial at matched fly speed
fr = t.fr; L = round(wsP*fr); pad = round(1.5*fr); Lf = round(0.45*fr); Lb = round(0.2*fr);
bad = logical(s.info.badVol(:)'); off = prctile(double(s.flash.meanImg(:)),0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'),'gaussian',s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad))); F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
fly = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',max(3,round(fr))))*fr);
lab = {'slower bump','faster bump'}; n = numel(t.ph); LK = gobjects(1,2); LH = LK; LS = LK;
IDX = arrayfun(@(e) max(1,pr(2+e)-pad):min(n,pr(2+e)+L-1+pad), 1:2, 'UniformOutput', false);
cl = prctile(D(:,[IDX{:}]),[2 99.5],'all');                    % one colour scale for both excerpts of the pair
toBin = @(p) mod(p+pi,2*pi)/(2*pi)*16 + 0.5;
for e = 1:2
    a = pr(2+e); idx = IDX{e}; tt = (idx - a)/fr;
    ax = nexttile(tl, (rowIdx-1)*6 + e); hold(ax,'on')
    imagesc(ax, tt, 1:32, D(:,idx)); clim(ax, cl); axis(ax,'tight'); set(ax,'YDir','reverse')
    pb = toBin(t.ph(idx)); pb(abs(diff([pb pb(end)])) > 8) = nan; plot(ax, tt, pb, 'w-', 'LineWidth', 1.3)
    xline(ax, 0, 'w--'); xline(ax, (L-1)/fr, 'w--'); set(ax,'XTick',[],'YTick',[8 24],'YTickLabel',{'half 1','half 2'})
    title(ax, sprintf('%s: %s', ttl, lab{e}), 'FontWeight','normal','Color',cols(e,:))
    ax2 = nexttile(tl, (rowIdx-1)*6 + 2 + e); hold(ax2,'on')
    ok = idx - Lf >= 1; plot(ax2, tt(ok), fly(idx(ok)-Lf), 'k-', 'LineWidth', 1.6)
    ok = idx - Lb >= 1; plot(ax2, tt(ok), bs(idx(ok)-Lb), '-', 'Color', [0 0.45 0.74], 'LineWidth', 1.6)
    xline(ax2, 0, ':'); xline(ax2, (L-1)/fr, ':'); xlim(ax2, [tt(1) tt(end)]); set(ax2,'XTick',[]); ylabel(ax2,'rad/s'); LS(e) = ax2;
    text(ax2, 0.02, 0.95, sprintf('window means: fly %.2f | bump %.2f rad/s', pr(4+e), pr(6+e)), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    if e == 1, legend(ax2, {'fly speed (shifted +0.45 s)','GRAB bump speed (shifted +0.2 s)'}, 'Box','off','Location','east','FontSize',8); end
    ax3 = nexttile(tl, (rowIdx-1)*6 + 4 + e); hold(ax3,'on')
    plot(ax3, tt, t.pk(idx), '-', 'Color', cols(e,:), 'LineWidth', 1.8)
    xline(ax3, 0, ':'); xline(ax3, (L-1)/fr, ':'); xlim(ax3, [tt(1) tt(end)]); xlabel(ax3,'time from window start (s)'); ylabel(ax3,'GRAB peak dF/F')
    text(ax3, 0.02, 0.95, sprintf('window mean %.3f', pr(8+e)), 'Units','normalized','VerticalAlignment','top','FontSize',9)
    LK(e) = ax3; LH(e) = ax;
end
linkaxes(LK,'y'); linkaxes(LS,'y'); ylim(LS(1), [0 max(cellfun(@(ax_) max(ylim(ax_)), num2cell(LS)))]);
colorbar(LH(2),'eastoutside');
end

function T = grab_folded(s, snr_frac, a16)
% folded GRAB bump: halves aligned by the trial's mean left-right phase
% offset and averaged -> 16-bin profile; phase + rho (|mean(z e^{i alpha})|)
xb = s.ft.xb(:)'; xf = s.ft.xf(:); fr = s.info.fps; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z0 = Z; Z0(isnan(Z0)) = 0;
vL = exp(1i*a16).'*Z0(1:16,:); vR = exp(1i*a16).'*Z0(17:32,:);
dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
Zf = mean(cat(3, Z(1:16,:), ring_shift(Z(17:32,:), dLR)), 3, 'omitnan'); Zf(isnan(Zf)) = 0;
v = exp(1i*a16).'*Zf / 16;
T.fr = fr; T.dark = s.info.dark; T.bad = bad;
T.ph = smoothdata(unwrap(angle(v)), 'gaussian', max(3, round(0.5*fr))); T.rho_raw = abs(v);
T.vr = interp1(xf, smoothdata(s.ft.r_speed(:),'gaussian',30), xb, 'linear', 'extrap');
T.pk = max(smoothdata(D,1,'movmean',3,'omitnan'),[],1);
T.cue = unwrap(angle(interp1(xf, exp(1i*s.ft.cue(:)), xb, 'linear', 'extrap')));
end

function Y = ring_shift(X, d)
% Y(alpha) = X(alpha - d) on the 16-bin ring, linear interpolation
nb = size(X,1); step = 2*pi/nb; a = linspace(-pi, pi-step, nb)';
pos = mod((a - d + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
Y = (1-w).*X(lo+1,:) + w.*X(hi+1,:);
end

function T = grab_trial_signals(s, snr_frac, a16)
% offset-corrected per-bin dF/F (PMT offset = darkest-pixel level; background
% fluctuations removed, level kept), low-SNR bins dropped; GRAB bump phase
% per PB half and combined; GRAB PB-max per half and overall; fly heading
xb = s.ft.xb(:)'; xf = s.ft.xf(:); fr = s.info.fps; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z(isnan(Z)) = 0;
vL = exp(1i*a16).'*Z(1:16,:); vR = exp(1i*a16).'*Z(17:32,:);
dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
sm = @(y) smoothdata(y, 'gaussian', max(3, round(0.5*fr)));
T.fr = fr; T.dark = s.info.dark; T.bad = bad;
T.vr = interp1(xf, smoothdata(s.ft.r_speed(:),'gaussian',30), xb, 'linear', 'extrap');
T.head = cumsum(T.vr)/fr;
T.pL = sm(unwrap(angle(vL))); T.pR = sm(unwrap(angle(vR))); T.pA = sm(unwrap(angle(vL + vR.*exp(1i*dLR))));
T.strL = abs(vL)./sum(abs(Z(1:16,:)),1); T.strR = abs(vR)./sum(abs(Z(17:32,:)),1);
Ds = smoothdata(D,1,'movmean',3,'omitnan');
T.pkL = max(Ds(1:16,:),[],1); T.pkR = max(Ds(17:32,:),[],1); T.pk = max(Ds,[],1);
T.cue = unwrap(angle(interp1(xf, exp(1i*s.ft.cue(:)), xb, 'linear', 'extrap')));
end

function [far, nfar, anti, Df, ftr] = far_field(s, t, snr_frac, Tw, minDeg)
% far-field PB signal: mean folded dF/F over bins that stayed > minDeg from
% the bump phase throughout the last Tw s (so not touched by the bump or
% its decaying trail). nfar = # such bins; anti = mean of bins > minDeg from
% the CURRENT phase only (no trail exclusion). ftr = far-field TROUGH: min
% over those untouched bins of the 3-bin (circular) smoothed profile.
% t = grab_folded output.
a16 = linspace(-pi, pi-2*pi/16, 16)'; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z(isnan(Z)) = 0;
vL = exp(1i*a16).'*Z(1:16,:); vR = exp(1i*a16).'*Z(17:32,:); dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
Df = mean(cat(3, D(1:16,:), ring_shift(D(17:32,:), dLR)), 3, 'omitnan');
ph = angle(exp(1i*t.ph)); n = numel(ph); K = round(Tw*t.fr); th = minDeg*pi/180;
dmin = inf(16, n);
for k = 0:K
    p = [nan(1,k) ph(1:n-k)];                       % phase k frames ago
    d = abs(angle(exp(1i*(a16 - p))));            % 16 x n distance
    d(isnan(d)) = inf; dmin = min(dmin, d);
end
dmin(:, 1:K) = nan;                                % trail not defined at start
M = dmin > th; Dm = Df; Dm(~M) = nan;
far = mean(Dm, 1, 'omitnan'); nfar = sum(M & ~isnan(Df), 1); far(nfar < 2) = nan;
Ds = (circshift(Df,1,1) + Df + circshift(Df,-1,1))/3; Ds(~M) = nan; ftr = min(Ds, [], 1, 'omitnan'); ftr(nfar < 2) = nan;
dc = abs(angle(exp(1i*(a16 - ph)))); Da = Df; Da(dc <= th) = nan; anti = mean(Da, 1, 'omitnan');
end


function H = far_half(s, fr, snr_frac, Tw, minDeg)
% per PB half: unwrapped phase (own PVA), rho, and far-field dF/F (bins
% > minDeg from that half's own phase throughout the last Tw s).
a16 = linspace(-pi, pi-2*pi/16, 16)'; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z(isnan(Z)) = 0;
K = round(Tw*fr); th = minDeg*pi/180; n = size(V,2);
for h = 1:2
    r = (h-1)*16 + (1:16); v = exp(1i*a16).'*Z(r,:)/16; ph = angle(v);
    H(h).ph = smoothdata(unwrap(ph),'gaussian',max(3,round(0.5*fr))); H(h).rho = abs(v); %#ok<AGROW>
    dmin = inf(16,n);
    for k = 0:K, p = [nan(1,k) ph(1:n-k)]; d = abs(angle(exp(1i*(a16 - p)))); d(isnan(d)) = inf; dmin = min(dmin,d); end
    dmin(:,1:K) = nan; Dh = D(r,:); Dh(~(dmin > th)) = nan; far = mean(Dh,1,'omitnan'); far(sum(~isnan(Dh),1) < 2) = nan;
    H(h).far = far; H(h).pk = max(smoothdata(D(r,:),1,'movmean',3,'omitnan'),[],1);
end
end

function [C, lag, N] = joint_coef(FF, trials, yname, lf, th, xnames, lagsS, bwin)
% per-trial standardized joint regression of FF(j).(yname) on |fly speed|
% (lag lf s) + |folded bump speed| (scanned lag) + optional extra regressors
% xnames (fields of FF, same-time as outcome unless named 'fly:<field>' =
% taken at the fly lag). Only clear-bump frames (rho_raw > th over the 1 s
% velocity window). Bump lag = max |frame-weighted mean coefficient| in bwin.
% C: trials x (2 + numel(xnames)) at that lag, order [fly bump extras].
if nargin < 8, bwin = [-0.5 1]; end
nx = numel(xnames); P = 2 + nx;
Ca = nan(numel(trials), numel(lagsS), P); NF = nan(numel(trials), numel(lagsS));
for jj = 1:numel(trials)
    t = FF(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(lf*fr); y = t.(yname);
    sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
    good = movmin(t.rho_raw,w) > th & ~movmax(t.bad,w);
    for li = 1:numel(lagsS)
        L = round(lagsS(li)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
        X = [sp(iF)' bs(i0)'];
        for q = 1:nx
            nm = xnames{q};
            if startsWith(nm,'fly:'), v = t.(nm(5:end)); X(:,end+1) = v(iF)'; else, v = t.(nm); X(:,end+1) = v(i1)'; end %#ok<AGROW>
        end
        m = good(i0) & ~t.bad(i1) & isfinite(y(i1)) & all(isfinite(X),2)'; if sum(m) < 50, continue, end
        Zq = zscore([X(m,:) y(i1(m))']); b = [ones(sum(m),1) Zq(:,1:P)] \ Zq(:,P+1);
        Ca(jj,li,:) = b(2:end); NF(jj,li) = sum(m);
    end
end
W = sum(Ca(:,:,2).*NF,2-1,'omitnan')./sum(NF.*~isnan(Ca(:,:,2)),1,'omitnan');
win = find(lagsS >= bwin(1) & lagsS <= bwin(2)); [~,k] = max(abs(W(win))); jb = win(k);
C = squeeze(Ca(:,jb,:)); if numel(trials) == 1, C = C(:)'; end
lag = lagsS(jb); N = NF(:,jb);
end

function [atr, cov] = anti_trough(Df, ph, fr, curDeg, Tw, trailDeg)
% per-frame trough opposite the CURRENT bump: min of the 3-bin (circular)
% smoothed folded profile over bins > curDeg from the current phase, also
% dropping bins that came within trailDeg of the bump at any time in the
% last Tw s (trail exclusion; Tw = 0 -> none). cov = # bins used.
a16 = linspace(-pi, pi-2*pi/16, 16)'; ph = angle(exp(1i*ph(:)')); n = numel(ph);
Ds = (circshift(Df,1,1) + Df + circshift(Df,-1,1))/3;
M = abs(angle(exp(1i*(a16 - ph)))) > curDeg*pi/180;
K = round(Tw*fr);
for k = 1:K
    p = [nan(1,k) ph(1:n-k)]; d = abs(angle(exp(1i*(a16 - p)))); d(isnan(d)) = inf;
    M = M & d > trailDeg*pi/180;
end
if K > 0, M(:,1:K) = false; end
Ds(~M) = nan; atr = min(Ds, [], 1, 'omitnan'); cov = sum(M & ~isnan(Df), 1); atr(cov < 2) = nan;
end

function S10 = trough_fig(SETS, trF, lagsS, ttlStr, fname)
% 3 x 3 figure (rows: trough ~ fly + bump; peak ~ fly + bump; peak ~ fly +
% trough + bump; columns: GRAB, LPsP, EPG7f; per-trial coefficients), all
% rows on the frames where trough field trF is defined. Prints the bump
% term also with net rotation + forward speed as controls ('ctrl').
ctl = {'fly:net','fly:fw'}; S10 = struct(); rng(1);
for sgi = 1:3
    FF = SETS{sgi,1};
    for j = 1:numel(FF), FF(j).pkS = FF(j).pk; FF(j).pkS(isnan(FF(j).(trF))) = nan; end
    for d_ = 0:1
        trials = find([FF.dark] == d_ & (d_ == 1 | SETS{sgi,3} >= 0.4)); lf = SETS{sgi,5}(d_+1); th = SETS{sgi,4};
        [C1, l1] = joint_coef(FF, trials, trF, lf, th, {}, lagsS, [-0.5 2]);
        [C2, l2] = joint_coef(FF, trials, 'pkS', lf, th, {}, lagsS);
        [C3, l3] = joint_coef(FF, trials, 'pkS', lf, th, {trF}, lagsS);
        K1 = joint_coef(FF, trials, trF, lf, th, ctl, lagsS, [-0.5 2]); K2 = joint_coef(FF, trials, 'pkS', lf, th, ctl, lagsS);
        K3 = joint_coef(FF, trials, 'pkS', lf, th, [ctl {trF}], lagsS); KK = {K1(:,2), K2(:,2), K3(:,2)};
        S10(sgi,d_+1).C = {C1(:,[1 2]), C2(:,[1 2]), C3(:,[1 3 2])}; S10(sgi,d_+1).lag = [l1 l2 l3];
        for r = 1:3
            Cr = S10(sgi,d_+1).C{r}; ok = ~isnan(Cr(:,end));
            k_ = KK{r}(~isnan(KK{r}));
            fprintf('%-5s %-11s row %d: %s | bump %+.2f (>0 %d/%d, p=%.2g) | ctrl bump %+.2f (>0 %d/%d)\n', SETS{sgi,6}, ifelse_c(d_), r, mat2str(round(median(Cr(ok,:),1),2)), ...
                median(Cr(ok,end)), sum(Cr(ok,end)>0), sum(ok), signrank(Cr(ok,end)), median(k_), sum(k_>0), numel(k_));
        end
    end
end
figure('Color','w','Position',[60 60 1400 950]); tl = tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
cF = [0.15 0.15 0.15]; cT = [0.85 0.4 0.1]; cB = [0 0.45 0.74]; ttl10 = {'EPG GRAB(DA)','LPsP syt7f','EPG GCaMP7f'};
rlab = {sprintf('trough\n(fly + bump)'), sprintf('peak\n(fly + bump)'), sprintf('peak\n(fly + trough + bump)')};
for r = 1:3, for sgi = 1:3
    ax = nexttile(tl, (r-1)*3 + sgi); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]);
    for d_ = 0:1
        Cr = S10(sgi,d_+1).C{r}; ok = ~isnan(Cr(:,end)); Cr = Cr(ok,:); np = size(Cr,2);
        if np == 2, cc = {cF, cB}; xl = {'fly','bump'}; else, cc = {cF, cT, cB}; xl = {'fly','trough','bump'}; end
        x0 = d_*4; xs = x0 + (1:np);
        for e = 1:np
            v = Cr(:,e); jx = xs(e) + min(max(0.1*randn(size(v)), -0.25), 0.25);
            scatter(ax, jx, v, 16, cc{e}, 'filled', 'MarkerFaceAlpha', 0.35, 'MarkerEdgeColor', 'none');
            plot(ax, xs(e) + [-0.28 0.28], median(v)*[1 1], '-', 'Color', cc{e}, 'LineWidth', 3);
        end
        if ~isempty(Cr), pv = signrank(Cr(:,end)); else, pv = nan; end
        if ~isempty(Cr), mb = median(Cr(:,end)); else, mb = nan; end
        text(ax, x0 + (np+1)/2, 1.42, sprintf('bump %+.2f, %s\nn=%d, lag %+.1f s', mb, ifelse_p(pv), sum(ok), S10(sgi,d_+1).lag(r)), 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color',cB);
        text(ax, x0 + (np+1)/2, -0.55, ifelse_c(d_), 'HorizontalAlignment','center','FontSize',9,'Color',[0.4 0.4 0.4]);
        XT{d_+1} = xs; XL{d_+1} = xl; %#ok<SAGROW>
    end
    xlim(ax, [0.3 7.7]); ylim(ax, [-0.6 1.45]); set(ax, 'XTick', [XT{:}], 'XTickLabel', [XL{:}], 'Box','off', 'TickDir','out');
    if sgi == 1, ylabel(ax, sprintf('%s\njoint coefficient (z)', rlab{r})); end
    if r == 1, title(ax, ttl10{sgi}, 'FontWeight','normal'); end
end, end
title(tl, ttlStr, 'FontSize', 11);
exportgraphics(gcf, fname, 'Resolution', 150);
end

function R = bump_from_speed(t, trF, lf, th, W, shuf)
% Predict LPsP bump displacement over non-overlapping W-s windows from
% (1) integrated ball velocity, (2) + integrated ball velocity x trough
% deviation (trough above/below what the ball speed predicts = more/less
% commanded rotation), (3) integrated sign(vr) x trough-calibrated speed
% (trough alone sets the speed). Trough is advanced by its fly lag lf to
% align with commanded timing; bump lag chosen on the training folds.
% 2-fold CV (odd/even windows). shuf = true circularly shifts the trough
% deviation by 1/3 trial (control). Returns R.r2 = [ball, ball+dev, trough-
% only] out-of-sample R^2, R.b = coefficient on the dev term (all windows).
fr = t.fr; n = numel(t.ph); vr = t.vr(:)'; ph = t.ph(:)';
tr = t.(trF)(:)'; ok = ~isnan(tr); if sum(ok) < 0.5*n, R = []; return, end
tr = interp1(find(ok), tr(ok), 1:n, 'linear', 'extrap');
Lf = round(lf*fr); tra = [tr(1+Lf:end) nan(1,Lf)];                % trough(t + lf) ~ commanded speed at t
sp = smoothdata(abs(vr),'gaussian',max(3,round(0.5*fr)));
m = ~isnan(tra) & ~t.bad; c = [ones(sum(m),1) sp(m)'] \ tra(m)';      % trough predicted from ball speed
dev = tra - (c(1) + c(2)*sp); dev = dev / std(dev(m),'omitnan');      % neural-minus-ball speed (z)
if shuf, dev = circshift(dev, round(n/3)); end
cal = [ones(sum(m),1) tra(m)'] \ sp(m)'; shat = max(cal(1) + cal(2)*tra, 0);   % trough in ball-speed units
Lw = round(W*fr); lagsB = round((0:0.1:1.5)*fr); good = t.rho_raw > th & ~t.bad;
starts = 1:Lw:n-Lw-max(lagsB)-Lf;
X = nan(numel(starts), 3, numel(lagsB)); Y = X(:,1,:); keep = false(numel(starts),1);
for k = 1:numel(starts)
    g = starts(k):starts(k)+Lw-1;
    if any(isnan(dev(g))) || any(t.bad(g)), continue, end
    for q = 1:numel(lagsB)
        gb = g + lagsB(q); if mean(good(gb)) < 0.8, continue, end
        Y(k,1,q) = ph(gb(end)) - ph(gb(1));
        X(k,:,q) = [sum(vr(g)) sum(vr(g).*dev(g)) sum(sign(vr(g)).*shat(g))]/fr;
        keep(k) = true;
    end
end
fold = mod(1:numel(starts),2)' + 1; r2 = nan(2,3); R.b = nan;
if sum(keep) < 20, R = []; return, end
for f = 1:2
    tr_ = fold ~= f; te = fold == f;
    % choose bump lag on training windows (ball-only fit)
    sc = arrayfun(@(q) corr_ok(X(tr_,1,q), Y(tr_,1,q)), 1:numel(lagsB)); [~,qb] = max(sc);
    x = X(:,:,qb); y = Y(:,1,qb); a = tr_ & all(~isnan([x y]),2); b = te & all(~isnan([x y]),2);
    if sum(a) < 10 || sum(b) < 10, continue, end
    mods = {1, [1 2], 3};
    for mm = 1:3
        A = [ones(sum(a),1) x(a,mods{mm})]; w = A \ y(a);
        yh = [ones(sum(b),1) x(b,mods{mm})] * w; r2(f,mm) = 1 - sum((y(b)-yh).^2)/sum((y(b)-mean(y(b))).^2);
    end
end
R.r2 = mean(r2,1,'omitnan');
sc = arrayfun(@(q) corr_ok(X(:,1,q), Y(:,1,q)), 1:numel(lagsB)); [~,qb] = max(sc); x = X(:,:,qb); y = Y(:,1,qb); a = all(~isnan([x y]),2);
z = zscore([x(a,1:2) y(a)]); w = [ones(sum(a),1) z(:,1:2)] \ z(:,3); R.b = w(3); R.bball = w(2); R.n = sum(a); R.lagB = lagsB(qb)/fr;
end

function r = corr_ok(x, y)
k = ~isnan(x) & ~isnan(y); if sum(k) < 10, r = -inf; else, r = corr(x(k), y(k)); end
end

function R = bump_from_speed2(t, trF, lf, th, W, shuf, nblk, wrapY)
% As bump_from_speed, but: (i) clarity required only at the two window
% endpoints (mean rho_raw over +-0.25 s > th, no bad frames there), since
% displacement only needs the phase at start and end; (ii) overlapping
% windows (step 0.5 s); (iii) blocked CV: trial split into nblk contiguous
% blocks, windows assigned by start time, leave-one-block-out, windows
% straddling a block edge dropped. Bump lag chosen on training blocks.
% wrapY = true: displacement = WRAPPED endpoint phase difference (immune to
% phase slips in faint mid-window frames); windows with |ball heading
% change| > 120 deg dropped so the true displacement stays < 180 deg.
% R.r2 = [ball, ball+dev, trough-set]; R.b = in-sample dev coef (z, with ball).
if nargin < 8, wrapY = false; end
fr = t.fr; n = numel(t.ph); vr = t.vr(:)'; ph = t.ph(:)';
tr = t.(trF)(:)'; ok = ~isnan(tr); if sum(ok) < 0.5*n, R = []; return, end
tr = interp1(find(ok), tr(ok), 1:n, 'linear', 'extrap');
Lf = round(lf*fr); tra = [tr(1+Lf:end) nan(1,Lf)];
sp = smoothdata(abs(vr),'gaussian',max(3,round(0.5*fr)));
m = ~isnan(tra) & ~t.bad; c = [ones(sum(m),1) sp(m)'] \ tra(m)';
dev = tra - (c(1) + c(2)*sp); dev = dev / std(dev(m),'omitnan');
if shuf, dev = circshift(dev, round(n/3)); end
cal = [ones(sum(m),1) tra(m)'] \ sp(m)'; shat = max(cal(1) + cal(2)*tra, 0);
h = round(0.25*fr); clr = movmean(t.rho_raw, 2*h+1) > th & ~movmax(t.bad, 2*h+1);
Lw = round(W*fr); lagsB = round((0:0.1:1.5)*fr); st = round(0.5*fr);
starts = 1:st:n-Lw-max(lagsB)-Lf; nq = numel(lagsB);
X = nan(numel(starts),3,nq); Y = nan(numel(starts),nq);
for k = 1:numel(starts)
    g = starts(k):starts(k)+Lw-1; if any(isnan(dev(g))), continue, end
    xv = [sum(vr(g)) sum(vr(g).*dev(g)) sum(sign(vr(g)).*shat(g))]/fr;
    for q = 1:nq
        e0 = g(1) + lagsB(q); e1 = g(end) + lagsB(q);
        if wrapY && abs(xv(1)) > 2*pi/3, continue, end
        if clr(e0) && clr(e1)
            if wrapY, Y(k,q) = angle(exp(1i*(ph(e1) - ph(e0)))); else, Y(k,q) = ph(e1) - ph(e0); end
            X(k,:,q) = xv;
        end
    end
end
edges = round(linspace(1, n+1, nblk+1)); blk = discretize(starts, edges)';
straddle = discretize(starts + Lw - 1 + max(lagsB), edges)' ~= blk; blk(straddle) = nan;
r2 = nan(nblk,3); yb = cell(1,3); yy = [];
for f = 1:nblk
    a0 = blk ~= f & ~isnan(blk); b0 = blk == f;
    sc = arrayfun(@(q) corr_cc(X(a0,1,q), Y(a0,q)), 1:nq); [~,qb] = max(sc);
    x = X(:,:,qb); y = Y(:,qb); a = a0 & all(~isnan([x y]),2); b = b0 & all(~isnan([x y]),2);
    if sum(a) < 10 || sum(b) < 3, continue, end
    mods = {1, [1 2], 3};
    for mm = 1:3
        w = [ones(sum(a),1) x(a,mods{mm})] \ y(a); yb{mm} = [yb{mm}; [ones(sum(b),1) x(b,mods{mm})]*w];
    end
    yy = [yy; y(b)]; %#ok<AGROW>
end
if numel(yy) < 15, R = []; return, end
R.r2 = cellfun(@(yh) 1 - sum((yy-yh).^2)/sum((yy-mean(yy)).^2), yb);     % pooled over held-out blocks
sc = arrayfun(@(q) corr_cc(X(:,1,q), Y(:,q)), 1:nq); [~,qb] = max(sc); x = X(:,:,qb); y = Y(:,qb); a = all(~isnan([x y]),2);
z = zscore([x(a,1:2) y(a)]); w = [ones(sum(a),1) z(:,1:2)] \ z(:,3);
R.b = w(3); R.bball = w(2); R.n = sum(a); R.nind = sum(a)*st/Lw; R.lagB = lagsB(qb)/fr;
end

function r = corr_cc(x, y)
k = ~isnan(x) & ~isnan(y); if sum(k) < 10, r = -inf; else, r = corr(x(k), y(k)); end
end

function [P, N] = xcentered_profiles(s, fr, vr, snr_frac, lf, edges)
% Mean cross-centered PB profile (dF/F) per |fly speed| bin. Each half is
% shifted so the OTHER half's bump phase (+/- the trial's left-right
% offset) sits at 0, so noise can't create a bump. All non-bad frames.
% |fly speed| (0.5 s smoothed) taken lf s before the imaging frame.
% P: 16 x nbins (bins at a16 relative to the bump), N: frames per bin.
a16 = linspace(-pi, pi-2*pi/16, 16)'; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z(isnan(Z)) = 0;
vL = exp(1i*a16).'*Z(1:16,:); vR = exp(1i*a16).'*Z(17:32,:); dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
w = max(3, round(0.5*fr));
p1 = angle(exp(1i*smoothdata(unwrap(angle(vL)),'gaussian',w)));   % half-1 phase (half-1 coords)
p2 = angle(exp(1i*smoothdata(unwrap(angle(vR)),'gaussian',w)));   % half-2 phase (half-2 coords)
C1 = ring_shift_c(D(1:16,:), -(p2 + dLR));      % half 1 centred on half 2's bump
C2 = ring_shift_c(D(17:32,:), -(p1 - dLR));     % half 2 centred on half 1's bump
n = size(D,2); Lf = round(lf*fr); sp = smoothdata(abs(vr(:)'),'gaussian',w);
spl = [nan(1,Lf) sp(1:n-Lf)];
[~,~,b] = histcounts(spl, edges); nb = numel(edges) - 1; P = nan(16,nb); N = zeros(1,nb);
for k = 1:nb
    m = b == k & ~bad; N(k) = sum(m); if N(k) < 20, continue, end
    P(:,k) = mean([mean(C1(:,m),2,'omitnan') mean(C2(:,m),2,'omitnan')],2,'omitnan');
end
end

function Y = ring_shift_c(X, d)
% per-column fractional circular shift: Y(alpha,t) = X(alpha - d(t), t)
nb = size(X,1); step = 2*pi/nb; a = linspace(-pi, pi-step, nb)';
pos = mod((a - d(:)' + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
cols = repmat(1:size(X,2), nb, 1);
Y = (1-w).*X(sub2ind(size(X), lo+1, cols)) + w.*X(sub2ind(size(X), hi+1, cols));
end

function [r, lag, sg, off] = cue_track(t, lagFix)
% How well the bump phase follows the visual cue: bump(t) ~ off - sg*cue(t+lag).
% Scans lag in -2.5..+1 s and sign (or uses lagFix s if given); returns the
% tracking strength r = |mean exp(i(bump + sg*cue(t+lag)))|, the lag, the
% sign and the constant offset off = angle of that mean.
n = numel(t.ph); ph = angle(exp(1i*t.ph)); cue = angle(exp(1i*t.cue)); r = 0; lag = nan; sg = nan; off = nan;
if nargin > 1 && ~isempty(lagFix), tts = round(lagFix*t.fr); else, tts = -round(2.5*t.fr):round(1*t.fr); end
for tt = tts
    idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n & ~t.bad(idx));
    for s_ = [-1 1]
        m = mean(exp(1i*(ph(idx) + s_*cue(idx+tt))));
        if abs(m) > r, r = abs(m); lag = tt/t.fr; sg = s_; off = angle(m); end
    end
end
end

function [P, N, cen, ref, v] = cue_ref_profiles(s, t, snr_frac, lagRef, edges, Df)
% Closed loop: folded dF/F profile aligned on the CUE-predicted bump
% position ref(t) = off - sg*cue(t + lagRef) (lagRef fixed, e.g. from EPG
% 7f; sign and constant offset fitted per trial), mirrored so the bump's
% direction of travel is +, binned by |ref speed| (rad/s). P: 16 x nbins
% (a16 = bins relative to ref, + = ahead), N: frames per bin, cen: circular
% centroid (deg) of P - min(P) per bin. Optional Df (16 x n) replaces the
% measured profile (e.g. a model), aligned with the same reference.
a16 = linspace(-pi, pi-2*pi/16, 16)';
if nargin < 6 || isempty(Df), [~,~,~,Df] = far_field(s, t, snr_frac, 2, 90); end
[~,~,sg,off] = cue_track(t, lagRef);
n = numel(t.ph); tt = round(lagRef*t.fr); cu = unwrap(t.cue(:)');
idx = (1:n) + tt; ok = idx >= 1 & idx <= n; ref = nan(1,n); ref(ok) = off - sg*cu(idx(ok));
v = gradient(smoothdata(ref,'gaussian',max(3,round(0.5*t.fr)),'omitnan'))*t.fr;
Y = shift_cols(Df, -ref);                                   % ref at bin 0
mir = mod(17 - (1:16), 16) + 1; Ym = Y(mir,:);              % theta -> -theta
Y(:, v < 0) = Ym(:, v < 0);                                 % direction of travel = +
[~,~,b] = histcounts(abs(v), edges); nb = numel(edges) - 1; P = nan(16,nb); N = zeros(1,nb); cen = nan(1,nb);
for k = 1:nb
    m = b == k & ~t.bad & ok & ~isnan(v); N(k) = sum(m); if N(k) < 30, continue, end
    P(:,k) = mean(Y(:,m), 2, 'omitnan');
    q = P(:,k) - min(P(:,k)); cen(k) = angle(sum(q.*exp(1i*a16)))*180/pi;
end
end

function Y = shift_cols(X, d)
% per-column fractional circular shift: Y(alpha,t) = X(alpha - d(t), t); NaN d -> NaN column
nb = size(X,1); step = 2*pi/nb; a = linspace(-pi, pi-step, nb)'; dd = d(:)'; bad = isnan(dd); dd(bad) = 0;
pos = mod((a - dd + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
cols = repmat(1:size(X,2), nb, 1);
Y = (1-w).*X(sub2ind(size(X), lo+1, cols)) + w.*X(sub2ind(size(X), hi+1, cols)); Y(:,bad) = nan;
end

function [ETA, nev, fitq] = turn_transients(y, vr, fr, bad, win)
% Turn onset / offset transients relative to a fitted speed model.
% Model: y ~ a + b * (|vr| delayed by L and low-passed by exp(-t/tau)),
% (L, tau) on a grid, chosen by R^2 over all non-bad frames. Events:
% onsets = smoothed |vr| crosses 1 rad/s upward with mean < 0.3 over the
% previous 1 s and > 1 over the next 1 s; offsets = the mirror image.
% ETA: struct with time axis and event-triggered means of y, model and
% residual (y - model) in units of the trial's y SD; nev = [#on #off].
sp = smoothdata(abs(vr(:)'),'gaussian',max(3,round(0.3*fr))); n = numel(sp); y = y(:)';
Ls = round((0:0.1:1)*fr); taus = [0.05 0.1 0.2 0.3 0.5 0.8 1.2]; best = -inf; ok = ~bad & isfinite(y);
for L = Ls
    for tau = taus
        k = exp(-(0:round(5*tau*fr))/(tau*fr)); k = k/sum(k);
        x = filter(k, 1, [nan(1,L) sp(1:n-L)]); x(isnan(x)) = 0; x(1:L) = nan;
        m = ok & ~isnan(x); c = [ones(sum(m),1) x(m)'] \ y(m)'; r2 = 1 - var(y(m) - (c(1)+c(2)*x(m)))/var(y(m));
        if r2 > best, best = r2; mdl = c(1) + c(2)*x; fitq = [L/fr tau r2]; end
    end
end
sd = std(y(ok)); W = round(1*fr); up = find(sp(1:end-1) < 1 & sp(2:end) >= 1) + 1; dn = find(sp(1:end-1) >= 1 & sp(2:end) < 1) + 1;
isOn  = arrayfun(@(i) i > W && i + W <= n && mean(sp(i-W:i-1)) < 0.3 && mean(sp(i:i+W-1)) > 1, up);
isOff = arrayfun(@(i) i > W && i + W <= n && mean(sp(i-W:i-1)) > 1 && mean(sp(i:i+W-1)) < 0.3, dn);
ev = {up(isOn), dn(isOff)}; nev = [numel(ev{1}) numel(ev{2})];
lagW = round(win*fr); ETA.t = (lagW(1):lagW(2))/fr;
for e = 1:2
    Yy = nan(numel(ev{e}), numel(ETA.t)); Mm = Yy; Ss = Yy;
    for q = 1:numel(ev{e})
        ii = ev{e}(q) + (lagW(1):lagW(2)); g = ii >= 1 & ii <= n; ii2 = ii(g);
        bb = bad(ii2); yy = y(ii2); yy(bb) = nan; mm = mdl(ii2); mm(bb) = nan;
        Yy(q,g) = yy/sd; Mm(q,g) = mm/sd; Ss(q,g) = sp(ii2);
    end
    ETA.y{e} = mean(Yy,1,'omitnan'); ETA.m{e} = mean(Mm,1,'omitnan'); ETA.r{e} = ETA.y{e} - ETA.m{e}; ETA.sp{e} = mean(Ss,1,'omitnan');
end
end

function H = bump_history(ref, fr, taus)
% Bump-history ("depletion") variables from a bump-position trace ref (rad,
% e.g. the cue-predicted position in closed loop; may contain NaN).
%  H.occ(i,:)  recent occupancy of the CURRENT location over the past
%              taus(i) s: exponentially weighted mean over past frames of
%              closeness c(dphi) = exp(k (cos dphi - 1)), c(45 deg) = 0.5
%              (excludes the current frame); high = bump has been here.
%  H.tsl       time since the bump last overlapped (|dphi| < 22.5 deg) its
%              current location BEFORE the current stay (the bump must have
%              been > 45 deg away in between); NaN if not within 30 s.
%  H.dwell     time the bump has stayed within 22.5 deg of its current spot.
n = numel(ref); k = log(2)/(1 - cos(pi/4)); refi = fillmissing(ref(:)','nearest');
H.occ = nan(numel(taus), n);
for it = 1:numel(taus)
    K = round(5*taus(it)*fr); w = exp(-(1:K)/(taus(it)*fr)); w = w/sum(w); acc = zeros(1,n); ws = zeros(1,n);
    for q = 1:K
        past = [nan(1,q) refi(1:n-q)]; c = exp(k*(cos(past - refi) - 1)); g = ~isnan(c);
        acc(g) = acc(g) + w(q)*c(g); ws(g) = ws(g) + w(q);
    end
    o = acc./ws; o(ws < 0.9) = nan; H.occ(it,:) = o;
end
H.tsl = nan(1,n); H.dwell = nan(1,n); Lmax = round(30*fr);
for i = 2:n
    if isnan(ref(i)), continue, end
    lo = max(1, i - Lmax); d = abs(angle(exp(1i*(refi(i-1:-1:lo) - refi(i)))));   % back in time
    away = find(d > pi/4, 1);                                                     % end of current stay
    inn = find(d >= pi/8, 1); if isempty(inn), H.dwell(i) = (i-lo)/fr; else, H.dwell(i) = (inn-1)/fr; end
    if isempty(away), continue, end
    prev = find(d(away:end) < pi/8, 1); if ~isempty(prev), H.tsl(i) = (away + prev - 1)/fr; end
end
end

function y = ifelse_on(b), if b, y = 'on'; else, y = 'off'; end, end

function [P, T, fitlog] = lpsp_input_fit(S, FF, trials, snr_frac, lagRef, edges, T, grid)
% Fit an input model of LPsP calcium in LPsP closed-loop trials:
%   F/F0(theta,t) = c0 + gx * x(t) + gb * xb(t) + k * T(theta - ref(t))
% x = |fly speed| and xb = |cue-predicted bump speed|, each delayed by lag
% and low-passed by exp(-t/tauG); ref = cue-predicted bump position (lag
% lagRef); T = LPsP bump template (16 x 1, peak 1 at bin 9, min 0).
% (tauG, lag) from grid (rows [tauG lag]) by pooled SSE; coefficients per
% trial, then median across trials. F/F0 = 1 + folded dF/F.
a16 = linspace(-pi, pi-2*pi/16, 16)'; tmpl = @(th) interp1([a16; pi], [T; T(1)], mod(th + pi, 2*pi) - pi, 'linear');
best = inf; fitlog = [];
for gi_ = 1:size(grid,1)
    tauG = grid(gi_,1); lagG = grid(gi_,2); sse = 0; coefs = nan(numel(trials),4); r2 = nan(numel(trials),1);
    for jj = 1:numel(trials)
        j = trials(jj); t = FF(j); s = S(j); n = numel(t.ph); fr = t.fr;
        [~,~,~,Df] = far_field(s, t, snr_frac, 2, 90);
        [~,~,~,ref,v] = cue_ref_profiles(s, t, snr_frac, lagRef, edges, zeros(16,n));
        [x, xb] = input_drives(t.vr, v, fr, tauG, lagG);
        E = tmpl(a16 - fillmissing(ref,'nearest'));
        Y = 1 + Df; m = ~t.bad & ~isnan(ref) & ~isnan(x) & ~isnan(xb);
        Ym = Y(:,m); A = [ones(numel(Ym),1) reshape(repmat(x(m),16,1),[],1) reshape(repmat(xb(m),16,1),[],1) reshape(E(:,m),[],1)];
        g = ~isnan(Ym(:)); c = A(g,:) \ Ym(g); res = Ym(g) - A(g,:)*c;
        sse = sse + sum(res.^2)/var(Ym(g)); coefs(jj,:) = c'; r2(jj) = 1 - var(res)/var(Ym(g));
    end
    fitlog(end+1,:) = [tauG lagG sse median(r2)]; %#ok<AGROW>
    if sse < best, best = sse; P = struct('c0', median(coefs(:,1)), 'gx', median(coefs(:,2)), 'gb', median(coefs(:,3)), 'k', median(coefs(:,4)), ...
            'tauG', tauG, 'lagG', lagG, 'coefs', coefs, 'r2', r2); end
end
end

function [x, xb] = input_drives(vr, v, fr, tauG, lagG)
% |fly speed| and |bump (cue) speed| (0.3 s smoothed), delayed by lagG s
% and low-passed by a causal exponential with time constant tauG s.
n = numel(vr); L = round(lagG*fr); k = exp(-(0:round(5*tauG*fr))/(tauG*fr)); k = k/sum(k);
a = smoothdata(abs(vr(:)'),'gaussian',max(3,round(0.3*fr))); b = smoothdata(abs(v(:)'),'gaussian',max(3,round(0.3*fr)),'omitnan');
b = fillmissing(b,'nearest');
a = filter(k, 1, [a(1)*ones(1,L) a(1:n-L)]); b = filter(k, 1, [b(1)*ones(1,L) b(1:n-L)]);
x = a; xb = b; x(1:numel(k)) = nan; xb(1:numel(k)) = nan;
end

function [Rel, R] = release_sim(L, fr, n_, tauRec, u)
% Depleting synapse per glomerulus (rows of L, 16 x T; L = presynaptic
% calcium proxy, e.g. F/F0). Drive D = L.^n_; vesicle availability R:
%   dR/dt = (1 - R)/tauRec - u*D*R,  release = D .* R.
% u = 0 -> no depletion (pure power law). Euler steps at the frame rate,
% R starts at its steady state for the first frame's drive.
D = max(L,0).^n_; [nb, T] = size(D); dt = 1/fr; R = nan(nb,T);
if u == 0, R(:) = 1; Rel = D; return, end
r = 1 ./ (1 + u*tauRec*D(:,1)); r(isnan(r)) = 1;
for k = 1:T
    d = D(:,k); d(isnan(d)) = 0; R(:,k) = r;
    r = r + dt*((1 - r)/tauRec - u*d.*r); r = min(max(r,0),1);
end
Rel = D .* R;
end

function [cost, out] = grab_model_eval(th, model, DAT, half)
% Cost of a release model on GRAB closed-loop data (sum over trials of
% normalized SSE on the frames of 'half': 1 = first half of each trial, 2 =
% second, 0 = all), GRAB(theta,t) ~ a + b * (sensor kernel * Rel), a and b
% least-squares per trial on the same frames. th (log-parameterized):
%  'product'  : [log tauS, log k2]           Rel = G(t) * (1 + k2 T(theta-ref))
%  'power'    : [log tauS, log n]            Rel = Lhat^n
%  'depletion': [log tauS, log n, log tauRec, log u]   Rel = Lhat^n R
%  'product_w': [log tauS, log k2, w]       Rel = (x + w xb) (1 + k2 T)  (free bump-speed weight)
%  'product_dep': [log tauS, log k2, log tauRec, log u]  drive = G (1 + k2 T), Rel = drive R (depleting)
%  'product_wdep': [log tauS, log k2, w, log tauRec, log u]
%  'gated_L'     : [log tauS, log k2]                      Rel = dL (1 + k2 T)
%  'gated_adapt' : [log tauS, log k2, logit beta, log tauA] Rel = [dL - beta LP_tauA(dL)]+ (1 + k2 T)
%  'gated_me'    : [log tauS, log k2, log tauA]           Rel = [dL - LP_tauA(dL)]+ (1 + k2 T)  (beta = 1)
%  'gated_adapt_glob': as gated_adapt but only the PB-wide part adapts:
%                Rel = [G - beta LP(G) + (dL - G)]+ (1 + k2 T)   (control)
%  dL = synthetic LPsP F/F0 minus its resting level (.dL), LP = causal
%  exponential low-pass per glomerulus (an adapting synapse; beta = the
%  adapting fraction, 1 = pure temporal-change / motion energy).
% DAT(j): .Y (16 x T GRAB dF/F), .L (16 x T synthetic LPsP F/F0), .G (1 x T
% LPsP-derived gain, F/F0 - 1 without the bump), .E (16 x T template at
% ref), .x/.xb (1 x T fly / bump speed drives), .dL (16 x T, L minus rest),
% .fr, .ok (1 x T usable frames).
tauS = exp(th(1)); cost = 0; out = struct('r2',[],'pred',[]);
for j = 1:numel(DAT)
    d = DAT(j); fr = d.fr;
    switch model
        case 'product',   Rel = d.G .* (1 + exp(th(2))*d.E);
        case 'power',     Rel = release_sim(d.L, fr, exp(th(2)), 1, 0);
        case 'depletion', Rel = release_sim(d.L, fr, exp(th(2)), exp(th(3)), exp(th(4)));
        case 'product_w', Rel = max(repmat(d.x + th(3)*d.xb,16,1),0) .* (1 + exp(th(2))*d.E);
        case 'product_dep', Rel = release_sim(d.G .* (1 + exp(th(2))*d.E), fr, 1, exp(th(3)), exp(th(4)));
        case 'product_wdep', Rel = release_sim(max(repmat(d.x + th(3)*d.xb,16,1),0) .* (1 + exp(th(2))*d.E), fr, 1, exp(th(4)), exp(th(5)));
        case 'gated_L',     Rel = d.dL .* (1 + exp(th(2))*d.E);
        case 'gated_adapt', Rel = max(d.dL - 1/(1+exp(-th(3)))*lp_exp(d.dL, exp(th(4))*fr), 0) .* (1 + exp(th(2))*d.E);
        case 'gated_adapt_glob', Rel = max(d.G - 1/(1+exp(-th(3)))*lp_exp(d.G, exp(th(4))*fr) + (d.dL - d.G), 0) .* (1 + exp(th(2))*d.E);
        case 'gated_me',    Rel = max(d.dL - lp_exp(d.dL, exp(th(3))*fr), 0) .* (1 + exp(th(2))*d.E);
    end
    k = exp(-(0:round(5*tauS*fr))/(tauS*fr)); k = k/sum(k); X = filter(k, 1, Rel, [], 2);
    T = size(X,2); sel = d.ok; if half == 1, sel = sel & (1:T) <= T/2; elseif half == 2, sel = sel & (1:T) > T/2; end
    sel(1:numel(k)) = false;
    Ys = d.Y(:,sel); Xs = X(:,sel); g = ~isnan(Ys) & ~isnan(Xs);
    A = [ones(sum(g(:)),1) Xs(g)]; c = A \ Ys(g); res = Ys(g) - A*c; v = var(Ys(g));
    cost = cost + sum(res.^2)/v/numel(res);
    out.r2(j) = 1 - var(res)/v; out.c(j,:) = c';
    if nargout > 1 && nargin > 3, P = nan(size(d.Y)); P(:,sel) = c(1) + c(2)*Xs; out.pred{j} = P; end
end
cost = cost/numel(DAT);
end

function Y = lp_exp(X, tauFr)
% causal exponential low-pass along time (columns), time constant tauFr frames;
% starts from the first column's value; NaN filled with the nearest value
X = fillmissing(X, 'nearest', 2); a = exp(-1/max(tauFr, 1e-3));
Y = filter(1 - a, [1 -a], X, (X(:,1)*a)', 2);
end

function c = occ_coef(y, t, ref, v, taus, flyHist, placebo)
% Section-15 "strict" occupancy test for one closed-loop trial: partial
% (standardized) coefficient of the recent occupancy of the bump's current
% location on y (e.g. GRAB peak or a model's predicted peak), with cubic
% terms of |fly speed|, net rotation, forward speed and bump speed, mean
% bump speed over the past 1/2/4 s, reversal covariates; first 30 s out.
% t needs .fr .vr .fw .net .bad; ref / v = cue-predicted bump position and
% velocity. flyHist = true adds mean |fly speed| over the past 1/2/4 s (and
% its square) as controls; placebo = true uses occupancy of the location
% OPPOSITE the bump instead. Returns one coefficient per tau.
if nargin < 6, flyHist = false; end, if nargin < 7, placebo = false; end
fr = t.fr; n = numel(y); y = y(:)'; Lf = round(0.45*fr); Lb = round(0.3*fr); sh = @(z,L) [nan(1,L) z(1:n-L)];
H = bump_history(ref, fr, taus); pb = @(z) [z(:) z(:).^2 z(:).^3];
if placebo   % occupancy of the location opposite the bump: closeness of past positions to (current + pi)
    kk = log(2)/(1 - cos(pi/4)); refi = fillmissing(ref(:)','nearest');
    for it = 1:numel(taus)
        K = round(5*taus(it)*fr); w = exp(-(1:K)/(taus(it)*fr)); w = w/sum(w); o = zeros(1,n);
        for q = 1:K, past = [nan(1,q) refi(1:n-q)]; cO = exp(kk*(cos(past - refi - pi) - 1)); cO(isnan(cO)) = 0; o = o + w(q)*cO; end
        o(1:K) = nan; H.occ(it,:) = o;
    end
end
sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bsp = abs(v);
vs = smoothdata(v,'gaussian',max(3,round(0.3*fr))); cur = sign(vs); opp = nan(1,n); nrev = nan(1,n); W = round(2*fr);
for i = W+1:n, seg = vs(i-W:i-1); opp(i) = mean(sign(seg) == -cur(i) & abs(seg) > 0.3); nrev(i) = sum(abs(diff(sign(seg(abs(seg) > 0.3)))) > 0); end
Cs = [pb(sh(sp,Lf)) pb(sh(t.net,Lf)) pb(sh(fw,Lf)) pb(sh(bsp,Lb)) sh(opp,Lb)' sh(nrev,Lb)'];
for wsec = [1 2 4], Cs = [Cs movmean(sh(bsp,Lb), [round(wsec*fr) 0])']; end %#ok<AGROW>
if flyHist, for wsec = [1 2 4], h_ = movmean(sh(sp,Lf), [round(wsec*fr) 0])'; Cs = [Cs h_ h_.^2]; end, end %#ok<AGROW>
early = (1:n)' <= round(30*fr); c = nan(1, numel(taus));
for it = 1:numel(taus)
    X = [Cs sh(H.occ(it,:),Lb)']; m = all(isfinite([X y']),2) & ~t.bad(:) & ~early;
    Z = zscore([X(m,:) y(m)']); b = [ones(sum(m),1) Z(:,1:end-1)] \ Z(:,end); c(it) = b(end);
end
end

function y = ifelse_p(p), if p < 0.001, y = 'p<0.001'; else, y = sprintf('p=%.2g', p); end, end
function y = ifelse_c(d), if d, y = 'dark'; else, y = 'closed loop'; end, end
function y = ifelse_s(d), if d, y = '--'; else, y = '-'; end, end
