%% Figure_4_claude
% Example EMPTY>KIR (control) fly/trial/snippet showing a well-tracking
% EPG bump: fly heading velocity vs. bump angular velocity, correlated as
% strongly as this dataset gets, using the FULLY-SMOOTHED pipeline
% "data scripts/lpsp_kir_claude.m" itself settles on in its own
% "RNAi-style velocity gain" section -- not the simpler raw-gradient
% method that script's earlier gain-scatter figures (its figures 4-9) use.
%
% SOURCE DATASET: lpsp_kir_redo_data_20240206.mat, empty>kir trials only
% (the LPsP>Kir dataset's control genotype -- see lpsp_kir_claude.m's own
% header for how genotype/fly ID/light-condition are recovered from
% all_data.meta and each trial's own trialSettings.csv).
%
% ==THIS DATASET IS ELLIPSOID BODY (EB), NOT PROTOCEREBRAL BRIDGE (PB)==
% Per explicit note: despite variable names inherited from the PB pipeline
% throughout this codebase (im.alpha, "PB angle", etc.), this specific
% dataset's im.z/im.mu/im.rho describe EB wedges. The EB carries exactly
% ONE copy of the EPG bump around its ring of wedges (unlike the PB, which
% carries the SAME bump duplicated across its two mirror-image
% hemispheres) -- so, unlike every PB figure in this folder, there is no
% left/right symmetry or hemisphere-pairing consideration here: im.mu is
% already a single, unambiguous bump-position estimate.
%
% ============================================================
% ENVIRONMENT NOTES (this script was authored, then actually run, on a Mac
% -- NOT the Windows machine lpsp_kir_claude.m itself assumes)
% ============================================================
% lpsp_kir_redo_data_20240206.mat was not present locally when this script
% was first written (data/ only had the OTHER datasets Figure_1-3_claude.m
% use) and has since been added directly to data/. Two Mac-specific
% consequences of running lpsp_kir_claude.m's own logic somewhere other
% than its native Windows environment, both fixed here rather than worked
% around:
%   1) meta stores literal Windows paths ("Z:\pablo\..."), but lpsp_kir_
%      claude.m's own trial_folder_name/trial_fly_id split on filesep,
%      which is '/' on a Mac -- so on a Mac that split silently no-ops and
%      every trial falls back to being treated as its own unique "fly"
%      (confirmed directly: this happened the first time this script ran,
%      collapsing 34 empty>kir trials into 34 "flies" instead of the real
%      16, breaking the "one point per fly" pooling this whole pipeline
%      depends on). Fixed with split_path (splits on '\' OR '/'), the same
%      fix Figure_1_claude.m / Figure_2_claude.m already apply for this
%      exact reason.
%   2) base_dir ("Z:\pablo\lpsp_kir_redo\", used to read each trial's own
%      trialSettings.csv for light condition) is a Windows-only network
%      share this Mac cannot reach, and this dataset's ft never has a
%      pattern field at all (0/74 trials, confirmed directly). See the
%      light-condition section below for the fallback this forced.
%
% ================================================================
% THE SMOOTHING PIPELINE, IN FULL (ported from lpsp_kir_claude.m's
% "RNAi-style velocity gain" section, itself adapted from
% lpsp_rnai_claude_v2.m / gain_scratch_claude.m)
% ================================================================
% Two raw signals feed everything below:
%   heading  = -ft.cue        (the fly's own closed-loop cue position,
%                               sign-flipped -- lpsp_kir_claude.m confirmed
%                               directly, corr(mu,-cue)=+0.64 on a sample
%                               trial, that this is a valid heading proxy
%                               in BOTH closed loop and dark for this
%                               dataset, since ft.cue keeps tracking the
%                               fly's own rotation even with no visible
%                               pattern)
%   bump pos = im.mu          (the fitted EPG bump position, one value per
%                               imaging frame, already on the EB's own
%                               wedge-angle circular scale)
%
% Step 1 -- UNWRAP before anything else. Both mu and heading are circular
%   (wrap at +/-pi); unwrap() is applied BEFORE any smoothing or
%   differentiation so a legitimate crossing of the wrap boundary doesn't
%   register as a spurious +/-2*pi-in-one-frame velocity spike. (The
%   dataset's ft.cue has zero NaN samples -- confirmed directly in
%   lpsp_kir_claude.m -- so, unlike lpsp_rnai_claude_v2.m's own dataset,
%   no NaN-gap-filling step is needed before this unwrap.)
%
% Step 2 -- GAUSSIAN SMOOTH each unwrapped trace, separately, each on its
%   own native timebase:
%     mu_smoothed      = smoothdata(unwrap(im.mu), 'gaussian', win_im)
%       where win_im = round(mu_smooth_s / dt_im), dt_im = the trial's own
%       imaging frame period.
%     heading_smoothed = smoothdata(unwrap(-ft.cue), 'gaussian', win_hd)
%       where win_hd = round(heading_smooth_s / dt_ft), dt_ft = the
%       trial's own fictrac sample period.
%   NON-OBVIOUS GOTCHA (flagged explicitly because it's easy to misread
%   mu_smooth_s/heading_smooth_s as "the Gaussian's standard deviation in
%   seconds" -- they are NOT): MATLAB's smoothdata('gaussian', window)
%   treats its second argument as a WINDOW LENGTH, and defines the
%   Gaussian's actual standard deviation internally as window/5. So the
%   EFFECTIVE smoothing std in real time is mu_smooth_s/5 seconds (resp.
%   heading_smooth_s/5) -- e.g. a "winning" mu_smooth_s of 0.5s from the
%   sweep below is a ~0.1s-std Gaussian smooth, not a ~0.5s-std one.
%
% Step 3 -- RESAMPLE mu onto the fictrac timebase, THEN DIFFERENTIATE:
%     bump_vel = gradient( interp1(xb, mu_smoothed, xf) ) / dt_ft
%     heading_vel = gradient( heading_smoothed ) / dt_ft   (already on xf)
%   Interpolating the ALREADY-SMOOTHED mu (not raw mu) onto xf keeps the
%   smoothing genuinely anti-aliasing the imaging-rate signal rather than
%   being partially undone by a subsequent upsample.
%
% Step 4 -- LAG: bump_vel is shifted `lag_frames` (fictrac frames) LATER
%   than heading_vel (i.e. the bump is expected to follow the fly's own
%   turning after a short sensorimotor delay), by trimming both vectors
%   to overlap (heading_vel = heading_vel_full(1:end-lag),
%   bump_vel = bump_vel_full(lag+1:end) for lag>0; mirrored for lag<0;
%   untouched for lag==0) -- the same shift-and-trim convention used
%   throughout this codebase (e.g. lpsp_compartments_claude_script.m's
%   lagged_corr). rho (im.rho, interpolated the same way as mu) is trimmed
%   identically so its samples still line up with the lag-shifted pair.
%
% Step 5 -- VALIDITY GATE (thresholds ported unchanged from
%   lpsp_kir_claude.m's own ORIGINAL, no-smoothing gain-scatter section --
%   not re-tuned for the smoothed pipeline):
%     valid = |heading_vel| > vel_thresh   (0.2 rad/s -- "the fly is
%                                            actually turning")
%           & |heading_vel| < vel_max      (5 rad/s -- exclude fictrac
%                                            tracking glitches)
%           & |bump_vel|    < bump_thresh  (10 rad/s -- exclude
%                                            gradient/unwrap artifacts at
%                                            wrap-around points)
%           &  rho          > rho_thresh   (0.2 -- "the bump is actually
%                                            discernable")
%
% Step 6 -- FIT: on valid samples only, bump_vel ~ 1 + heading_vel (a
%   type-2/y-on-x-style linear fit WITH intercept: b = [1,heading_vel] \
%   bump_vel; gain = b(2)). A fly/trial only gets a gain/r value at all if
%   n_valid >= 50 AND mean(|heading_vel|>vel_thresh) >= 0.01 (this last
%   "some real turning happened" activity floor, and the n_valid=50 floor,
%   are both ported from lpsp_rnai_claude_v2.m's own fly_gain_cached --
%   lpsp_kir_claude.m's ORIGINAL gain fit had no activity floor at all).
%
% Step 7 -- CHOOSE (mu_smooth_s, heading_smooth_s, lag_frames) THEMSELVES
%   by a 3-stage sequential/greedy sweep, exactly as lpsp_kir_claude.m
%   does: fix the two not-yet-swept stages at neutral defaults
%   (heading_smooth_s=0.1s, lag=0), sweep mu_smooth_s over
%   [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3] s and keep the value minimizing
%   mean((fly_gain-1)^2) pooled over empty>kir CLOSED-LOOP flies (one gain
%   value per fly, that fly's own trials pooled first) -- i.e. the
%   smoothing/lag choice that makes the control genotype's own gain sit as
%   close to the biologically-expected value of 1 as possible, NOT the
%   choice that maximizes correlation directly (a subtly different
%   objective, ported as-is for consistency with lpsp_kir_claude.m). Then
%   fix mu_smooth_s at its winner and sweep heading_smooth_s over the same
%   grid; then fix both and sweep lag_frames over -10:1:40 (fictrac
%   frames, this rig's own frame grid, per lpsp_kir_claude.m).
%
% ================================================================
% WHAT THIS SCRIPT ADDS ON TOP OF THAT PIPELINE
% ================================================================
% lpsp_kir_claude.m applies the winning pipeline to every genotype x light
% condition GROUP and plots one summary gain value per fly (its figure
% 35). This script instead asks a narrower, illustrative question: among
% CLOSED-LOOP TRIALS of a given genotype specifically (one trial at a
% time, not a whole fly's pooled trials, so there is a single well-defined
% trace to show), which is REPRESENTATIVE of that genotype's own typical
% heading-velocity/bump-velocity gain and correlation, using that exact
% winning pipeline (itself optimized on empty>kir flies only, then applied
% UNCHANGED to lpsp>kir trials too, same as lpsp_kir_claude.m's own
% convention) -- per explicit request, NOT the single highest-r/highest-
% gain trial (a best-case pick), but the trial whose own gain and r sit
% jointly closest to that genotype's own fly-level mean gain and mean r
% (see the "pick TWO example flies" section below for the exact distance
% metric). Within that representative trial, which snippet_s-second window
% is shown is chosen separately, by widest CIRCULAR COVERAGE of heading
% position among its discernible-bump samples (pick_wide_heading_snippet,
% same objective as Figure_1_claude.m's own pick_snippet) -- so the
% snippet shows the bump tracking through many different headings, per
% explicit request, rather than whichever window happens to correlate
% best locally. Representativeness (gain/r vs. genotype mean) applies to
% WHICH TRIAL is picked; heading coverage applies to WHICH WINDOW within
% it is shown -- two separate criteria. TWO such trials are shown, side by
% side as two columns, one per genotype: representative empty>kir (column
% 1) and representative lpsp>kir (column 2).
%
% The "fly velocity" signal used throughout (for ranking, the position/
% velocity traces, and the scatter) is the cue-derived heading velocity --
% i.e. -ft.cue's own velocity, which reflects the fly's speed relative to
% the closed-loop VR environment, not ft.r_speed (the fictrac ball's own
% raw rotational speed) directly. An earlier version of this script also
% compared against ft.r_speed directly and found the two traces track each
% other closely in shape but differ in gain (bump-vs-cue-heading gain
% =0.70 vs. bump-vs-r_speed gain=0.51 for one example trial) -- since the
% VR-relative signal is the one that matters here, that r_speed comparison
% was removed rather than kept as a second, unused code path.

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
data_dir    = fullfile(repo_root,'data');
export_dir  = fullfile(repo_root,'ugly_figures','exports');
all_figs_dir = fullfile(repo_root,'ugly_figures','all_figs');
base_dir    = 'Z:\pablo\lpsp_kir_redo\'; % raw-session tree; needed for trialSettings.csv (light condition) -- see below

%% load data
source_file = 'lpsp_kir_redo_data_20240206.mat';
tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype per trial (ported from lpsp_kir_claude.m)
genotype = cell(1,n_trials);
for i = 1:n_trials
    tname = trial_folder_name(all_data(i).meta);
    if contains(tname,'_LPsP_kir','IgnoreCase',true)
        genotype{i} = 'lpsp>kir';
    elseif contains(tname,'_empty_kir','IgnoreCase',true)
        genotype{i} = 'empty>kir';
    else
        genotype{i} = '';
    end
end
empty_rows = find(strcmp(genotype,'empty>kir'));
lpsp_rows  = find(strcmp(genotype,'lpsp>kir'));
fprintf('%d/%d trials are empty>kir\n', numel(empty_rows), n_trials);

%% light condition per trial: ft.pattern when stored, else this trial's own
% trialSettings.csv on the raw-session tree (base_dir), else (this Mac's
% actual situation: ft.pattern is never populated in this dataset -- 0/74
% trials, confirmed again directly here -- and base_dir, "Z:\pablo\...",
% is a Windows-only network path unreachable from here) a THIRD fallback:
% the trial-number-PARITY heuristic lpsp_kir_script_minimal.m itself uses
% (its own line ~81: regex out every digit run in the trial folder name,
% e.g. "20231109-3_EPG_7f_empty_kir" -> {"20231109","3","7"}, take the
% SECOND run -- the number right after the date-dash, i.e. the trial
% number -- and call it EVEN=dark, ODD=closed loop). lpsp_kir_claude.m's
% own header calls this "fragile if a session ever skips or repeats a
% trial" and prefers the CSV instead -- but checked directly against every
% one of this dataset's 74 trials (grouped by fly): EVERY fly's own trial
% numbers form clean, unbroken odd/even CL/dark pairs (e.g. fly with
% trials [9,10,11,12] -> CL,dark,CL,dark) with no skips or repeats
% anywhere, so the fragility this heuristic warns about does not actually
% occur in this dataset -- confirmed exhaustively, not assumed.
has_pattern = arrayfun(@(s) isfield(s.ft,'pattern') && ~isempty(s.ft.pattern), all_data);
is_dark = nan(1,n_trials);
is_dark_source = repmat({'none'},1,n_trials);

if ~all(has_pattern) && isfolder(base_dir)
    d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
    csv_by_name = containers.Map('KeyType','char','ValueType','char');
    for k = 1:numel(d)
        [~,tname] = fileparts(fileparts(d(k).folder)); % strip trailing \csv
        csv_by_name(lower(tname)) = fullfile(d(k).folder,d(k).name);
    end
    fprintf('indexed %d trialSettings.csv under %s\n', numel(d), base_dir);
end

for i = 1:n_trials
    if has_pattern(i)
        pattern_str = char(all_data(i).ft.pattern);
        is_dark(i) = contains(pattern_str,'background','IgnoreCase',true);
        is_dark_source{i} = 'ft.pattern';
    elseif isfolder(base_dir) && isKey(csv_by_name,lower(trial_folder_name(all_data(i).meta)))
        T = readtable(csv_by_name(lower(trial_folder_name(all_data(i).meta))),'TextType','string');
        pattern_str = char(T.patternPath(1));
        is_dark(i) = contains(pattern_str,'background','IgnoreCase',true);
        is_dark_source{i} = 'trialSettings.csv';
    else
        trial_num = trial_number(all_data(i).meta);
        if ~isnan(trial_num)
            is_dark(i) = mod(trial_num,2) == 0; % even -> dark, odd -> closed loop
            is_dark_source{i} = 'trial-number parity (fallback)';
        end
    end
end

if any(strcmp(is_dark_source,'trial-number parity (fallback)'))
    warning(['light condition for %d/%d trials came from the trial-number-parity FALLBACK heuristic ', ...
        '(ported from lpsp_kir_script_minimal.m), not ft.pattern or trialSettings.csv, since neither is reachable ', ...
        'from this Mac -- verified directly that every fly''s own trial numbers form clean, unbroken odd/even pairs ', ...
        'in this dataset, so the heuristic''s own known failure mode (a skipped/repeated trial number) does not occur here.'], ...
        sum(strcmp(is_dark_source,'trial-number parity (fallback)')), n_trials)
end

empty_cl_rows   = empty_rows(is_dark(empty_rows)==0);
lpsp_cl_rows    = lpsp_rows(is_dark(lpsp_rows)==0);
empty_dark_rows = empty_rows(is_dark(empty_rows)==1);
lpsp_dark_rows  = lpsp_rows(is_dark(lpsp_rows)==1);
fprintf('%d empty>kir closed-loop candidate trials for the example (column 1)\n', numel(empty_cl_rows));
fprintf('%d lpsp>kir closed-loop candidate trials for the example (column 2)\n', numel(lpsp_cl_rows));
fprintf('%d empty>kir dark trials, %d lpsp>kir dark trials (row 5 group comparison only)\n', numel(empty_dark_rows), numel(lpsp_dark_rows));

%% fly ID per trial (ported from lpsp_kir_claude.m)
fly_id = cell(1,n_trials);
for i = 1:n_trials
    fly_id{i} = trial_fly_id(all_data(i).meta);
end

%% shared thresholds -- ported UNCHANGED from lpsp_kir_claude.m's original
% (non-smoothed) gain-scatter section; the RNAi-style pipeline inherits
% these as-is, only adding the smoothing/lag stages on top.
vel_thresh  = 0.2;  % rad/s, "the fly is turning"
bump_thresh = 10;   % rad/s, exclude gradient/unwrap artifacts at wrap points
rho_thresh  = 0.2;  % im.rho, "the bump is discernable"
vel_max     = 5;    % rad/s, exclude fictrac tracking glitches

min_vel_n = 100; % minimum valid samples to trust a single TRIAL's own r (matches Figure_1_claude.m's min_vel_n)

%% ================= sweep 1/3: bump (mu) Gaussian smoothing (s) =================
% empty>kir, closed loop only; heading smoothing + lag held at neutral
% defaults while this stage is swept (ported from lpsp_kir_claude.m).
opt_flies = unique(fly_id(empty_cl_rows));
opt_fly_trials = cellfun(@(f) empty_cl_rows(strcmp(fly_id(empty_cl_rows),f)), opt_flies, 'UniformOutput',false);
fprintf('\noptimizing RNAi-style smoothing/lag pipeline on %d empty>kir closed-loop flies (%d trials)\n', numel(opt_flies), numel(empty_cl_rows));

default_heading_smooth_s = 0.1;
default_lag_frames       = 0;

mu_smooth_candidates_s = [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];
fly_gc_mu = nan(numel(opt_flies),numel(mu_smooth_candidates_s));
fprintf('\n=== sweep 1/3: bump (mu) Gaussian smoothing (s) ===\n');
for c = 1:numel(mu_smooth_candidates_s)
    msm = mu_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        fly_gc_mu(k,c) = fly_gain_v2(all_data,opt_fly_trials{k},msm,default_heading_smooth_s,default_lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
    end
    crit_c = mean((fly_gc_mu(:,c)-1).^2,'omitnan');
    fprintf('  mu_smooth=%.2fs: mean gain=%.3f, MSE-from-1=%.4f (n=%d flies)\n', ...
        msm, mean(fly_gc_mu(:,c),'omitnan'), crit_c, sum(~isnan(fly_gc_mu(:,c))));
end
crit_mu = mean((fly_gc_mu-1).^2,1,'omitnan');
[~,best_muc] = min(crit_mu);
mu_smooth_opt_s = mu_smooth_candidates_s(best_muc);
fprintf('winner: bump (mu) smoothing = %.2fs\n', mu_smooth_opt_s);

%% ================= sweep 2/3: heading (cue) Gaussian smoothing (s) =================
heading_smooth_candidates_s = [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];
fly_gc_hd = nan(numel(opt_flies),numel(heading_smooth_candidates_s));
fprintf('\n=== sweep 2/3: heading (cue) Gaussian smoothing (s), bump=%.2fs fixed ===\n',mu_smooth_opt_s);
for c = 1:numel(heading_smooth_candidates_s)
    hsm = heading_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        fly_gc_hd(k,c) = fly_gain_v2(all_data,opt_fly_trials{k},mu_smooth_opt_s,hsm,default_lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
    end
    crit_c = mean((fly_gc_hd(:,c)-1).^2,'omitnan');
    fprintf('  heading_smooth=%.2fs: mean gain=%.3f, MSE-from-1=%.4f (n=%d flies)\n', ...
        hsm, mean(fly_gc_hd(:,c),'omitnan'), crit_c, sum(~isnan(fly_gc_hd(:,c))));
end
crit_hd = mean((fly_gc_hd-1).^2,1,'omitnan');
[~,best_hdc] = min(crit_hd);
heading_smooth_opt_s = heading_smooth_candidates_s(best_hdc);
fprintf('winner: heading (cue) smoothing = %.2fs\n', heading_smooth_opt_s);

%% ================= sweep 3/3: lag (fictrac frames) =================
lag_candidates = -10:1:40;
fly_gc_lag = nan(numel(opt_flies),numel(lag_candidates));
fprintf('\n=== sweep 3/3: lag (frames), bump=%.2fs, heading=%.2fs fixed ===\n',mu_smooth_opt_s,heading_smooth_opt_s);
for c = 1:numel(lag_candidates)
    lg = lag_candidates(c);
    for k = 1:numel(opt_flies)
        fly_gc_lag(k,c) = fly_gain_v2(all_data,opt_fly_trials{k},mu_smooth_opt_s,heading_smooth_opt_s,lg, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
    end
end
crit_lag = mean((fly_gc_lag-1).^2,1,'omitnan');
[~,best_lagc] = min(crit_lag);
lag_frames_opt = lag_candidates(best_lagc);
fprintf('winner: lag = %d frames\n', lag_frames_opt);

fprintf('\n=== winning pipeline ===\n  bump (mu) smoothing: %.2fs\n  heading smoothing:   %.2fs\n  lag: %d frames\n', ...
    mu_smooth_opt_s, heading_smooth_opt_s, lag_frames_opt);

%% ================= rank CLOSED-LOOP TRIALS by r, using the winning pipeline =================
% one r (and gain) per TRIAL (not per fly-pooled-across-trials, unlike
% lpsp_kir_claude.m's own summary figures) so each eventual example is a
% single, nameable trial with a single trace to show. The winning pipeline
% itself was optimized on empty>kir flies only (see above); it is APPLIED
% UNCHANGED to lpsp>kir trials here too, same as lpsp_kir_claude.m's own
% approach of applying one dataset-wide winning pipeline to every genotype.
trial_r    = nan(1,n_trials);
trial_gain = nan(1,n_trials);
trial_n_valid = zeros(1,n_trials);
for i = [empty_cl_rows(:); lpsp_cl_rows(:)]'
    [hv,bv,vd] = trial_gain_vectors_v2(all_data(i),all_data(i).im.mu,all_data(i).im.rho, ...
        mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt,vel_thresh,bump_thresh,rho_thresh,vel_max);
    trial_n_valid(i) = sum(vd);
    if trial_n_valid(i) >= min_vel_n
        trial_r(i) = corr(hv(vd),bv(vd));
        b = [ones(trial_n_valid(i),1),hv(vd)] \ bv(vd);
        trial_gain(i) = b(2);
    end
end

%% ================= pick TWO example flies, one per genotype (both closed loop) =================
% Per explicit request: NOT the single highest-r/highest-gain trial (which
% would cherry-pick each genotype's best case) -- instead, the trial whose
% OWN r and gain sit jointly closest to that genotype's own mean r and
% mean gain, i.e. a REPRESENTATIVE example rather than a best-case one.
%
% The reference mean/std is computed per FLY (pooling each fly's own
% closed-loop trials first via fly_gain_corr_v2, same "one point per fly"
% convention used throughout this codebase/script -- e.g. the smoothing
% sweep above already pools trials per fly before averaging), so a fly
% with 2 closed-loop trials doesn't pull the genotype mean toward itself
% twice. "Closest" is then evaluated per TRIAL (not per fly) against that
% mean/std, since a single continuous trial is what gets displayed;
% r and gain are z-scored by the fly-level population's own std before
% combining into one Euclidean distance, so neither metric dominates
% just because it happens to have a larger raw scale/spread.
geno_cl_rows = {empty_cl_rows, lpsp_cl_rows};
example_geno_label = {'empty>kir','lpsp>kir'};
example_i = nan(1,2);
gallery   = cell(1,2); % per genotype: ALL its candidate trial indices, ranked z-dist-ascending (closest to genotype mean first) -- feeds the browsing gallery further below
gallery_z = cell(1,2); % per genotype: those trials' own z-dist, same order as gallery{g}
group_fly_gain = cell(1,4); % per-fly gain/r for all 4 genotype x light-condition groups (1=empty CL, 2=lpsp CL, 3=empty dark, 4=lpsp dark) -- CL groups (1-2) filled here, dark groups (3-4) filled further below; feeds row 5's group comparison
group_fly_r    = cell(1,4);
for g = 1:2
    rows_g = geno_cl_rows{g};
    flies_g = unique(fly_id(rows_g));
    fly_gain_g = nan(numel(flies_g),1);
    fly_r_g    = nan(numel(flies_g),1);
    for f = 1:numel(flies_g)
        trial_list = rows_g(strcmp(fly_id(rows_g),flies_g{f}));
        [fly_gain_g(f),fly_r_g(f)] = fly_gain_corr_v2(all_data,trial_list, ...
            mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt,vel_thresh,bump_thresh,rho_thresh,vel_max);
    end
    group_fly_gain{g} = fly_gain_g;
    group_fly_r{g}    = fly_r_g;

    mean_gain_g = mean(fly_gain_g,'omitnan'); std_gain_g = std(fly_gain_g,'omitnan');
    mean_r_g    = mean(fly_r_g,'omitnan');    std_r_g    = std(fly_r_g,'omitnan');
    fprintf('\n%s: fly-level mean gain=%.3f (std=%.3f), mean r=%.3f (std=%.3f), n=%d flies\n', ...
        example_geno_label{g},mean_gain_g,std_gain_g,mean_r_g,std_r_g,sum(~isnan(fly_gain_g)));

    z_dist = sqrt(((trial_gain(rows_g)-mean_gain_g)/std_gain_g).^2 + ((trial_r(rows_g)-mean_r_g)/std_r_g).^2);
    [~,rel] = min(z_dist);
    example_i(g) = rows_g(rel);
    fprintf('  representative trial: fly "%s", gain=%.3f, r=%.3f (z-dist to genotype mean=%.2f)\n', ...
        fly_short_label(fly_id{example_i(g)}),trial_gain(example_i(g)),trial_r(example_i(g)),z_dist(rel));

    [~,z_order] = sort(z_dist); % ascending: closest-to-mean first (all entries here are already non-NaN, so sort()'s NaN-ordering quirk doesn't apply)
    gallery{g}   = rows_g(z_order);
    gallery_z{g} = z_dist(z_order);
end
assert(all(isfinite(example_i)), 'could not find a valid closed-loop example trial in one or both genotypes')

%% ================= row 5 data: per-fly gain/r for the two DARK groups too =================
% same fly-level pooling as the CL groups above (fly_gain_corr_v2, one
% point per fly), just without any example-trial selection -- dark trials
% aren't used for the example columns, only for this group comparison.
geno_dark_rows = {empty_dark_rows, lpsp_dark_rows};
for g = 1:2
    rows_g = geno_dark_rows{g};
    flies_g = unique(fly_id(rows_g));
    fly_gain_g = nan(numel(flies_g),1);
    fly_r_g    = nan(numel(flies_g),1);
    for f = 1:numel(flies_g)
        trial_list = rows_g(strcmp(fly_id(rows_g),flies_g{f}));
        [fly_gain_g(f),fly_r_g(f)] = fly_gain_corr_v2(all_data,trial_list, ...
            mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt,vel_thresh,bump_thresh,rho_thresh,vel_max);
    end
    group_fly_gain{2+g} = fly_gain_g;
    group_fly_r{2+g}    = fly_r_g;
    fprintf('%s (dark): fly-level mean gain=%.3f, mean r=%.3f, n=%d flies\n', ...
        example_geno_label{g},mean(fly_gain_g,'omitnan'),mean(fly_r_g,'omitnan'),sum(~isnan(fly_gain_g)));
end

% manual override: force a specific gallery rank per genotype (see
% Figure_4_claude_scratch_gallery.png for the ranked candidates, rank 1 =
% the automatic closest-to-mean pick above) instead of that automatic
% pick (NaN = automatic). Per explicit request: rank 13 for empty>kir,
% rank 4 for lpsp>kir.
forced_rank = [13, 4];
for g = 1:2
    if ~isnan(forced_rank(g))
        example_i(g) = gallery{g}(forced_rank(g));
    end
end

%% ================= per-example: snippet + display traces =================
% within each representative trial, the displayed snippet_s-second window
% is the one whose discernible-bump samples span the WIDEST range of
% heading positions (pick_wide_heading_snippet, called inside
% build_example_display) -- see this file's header comment for why.
snippet_s      = 60; % s
snippet_step_s = 10; % s
min_window_n   = 50; % minimum valid samples within a candidate window

clear ex % ex(c) = E below needs to create ex fresh from E's own fields -- struct() would pre-allocate a 0-field struct that whole-struct assignment then rejects as "dissimilar structures"
for c = 1:2
    i = example_i(c);
    E = build_example_display(all_data(i),mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max,snippet_s,snippet_step_s,min_window_n);
    E.i = i;
    E.fly_short = fly_short_label(fly_id{i});
    E.trial_short = trial_short_label(all_data(i).meta);
    if isnan(is_dark(i))
        E.light_str = 'light condition unknown';
    elseif is_dark(i)
        E.light_str = 'dark';
    else
        E.light_str = 'closed loop';
    end
    ex(c) = E; %#ok<SAGROW>

    fprintf('\nexample %d (%s): fly "%s", trial "%s" (%s) -- r=%.3f, gain=%.3f, n_valid=%d, snippet [%.1f, %.1f] s\n', ...
        c,example_geno_label{c},E.fly_short,E.trial_short,E.light_str,trial_r(i),trial_gain(i),trial_n_valid(i),E.snippet_t0,E.snippet_t1);
end

%% ================= plot: 5 rows x 2 columns (rows 1-4: one column per example fly; row 5: genotype x light-condition group comparison) =================
bump_color    = [0 0.7 0.7]; % EPG blue (matches this codebase's own EPG-calcium/bump color, e.g. Figure_3_claude.m's gcamp_color)
heading_color = [0 0 0];     % black -- fly heading (position and velocity, cue-derived)
xl = [-vel_max,vel_max]; yl = [-bump_thresh,bump_thresh]; % shared scatter axis limits, both columns

fig = figure('color','w','Position',[50 50 1500 1850]); clf
t = tiledlayout(fig,5,2,'TileSpacing','loose','Padding','compact');

for c = 1:2
    E = ex(c);
    i = E.i;
    t_rel_xb = E.xb(E.xb_idx) - E.snippet_t0; % time relative to snippet start, so both columns share one x-axis
    t_rel_xf = E.xf(E.xf_idx) - E.snippet_t0;

    %% row 1: EB bump activity (im.z) heatmap for the snippet -- heatmap
    % ONLY, no position overlay (that's its own row below) -- ONE bump copy
    % across the EB's wedges (im.alpha), unlike a PB heatmap's two mirrored
    % hemispheres.
    ax1 = nexttile(t,c); hold(ax1,'on')
    z_snip = E.ex_im.z(:,E.xb_idx);
    z_clim = [prctile(min(z_snip,[],1),5), prctile(max(z_snip,[],1),95)];
    imagesc(ax1,t_rel_xb,unwrap(E.ex_im.alpha),z_snip,z_clim)
    colormap(ax1,white_to_color(bump_color))
    cb1 = colorbar(ax1); cb1.Label.String = 'z-score';
    xlim(ax1,[0,snippet_s]); ylim(ax1,[min(unwrap(E.ex_im.alpha)),max(unwrap(E.ex_im.alpha))])
    xlabel(ax1,'time (s)'); ylabel(ax1,'EB wedge angle (rad)')
    title(ax1,sprintf('%s example: fly %s, trial "%s" (%s, %.0fs snippet)\nEB bump activity (im.z)', ...
        example_geno_label{c},E.fly_short,E.trial_short,E.light_str,snippet_s),'Interpreter','none')

    %% row 2: bump position (mu, smoothed, EPG blue) and heading (-cue,
    % smoothed, black) traces, in their own row directly below the heatmap.
    mu_plot = wrap_to_pi(E.mu_smoothed_disp(E.xb_idx));
    mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan; % break the line at circular wraps, not connect across them
    heading_plot = wrap_to_pi(E.heading_smoothed_disp(E.xf_idx));
    heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;

    ax2 = nexttile(t,2+c); hold(ax2,'on')
    plot(ax2,t_rel_xb,mu_plot,'Color',bump_color,'LineWidth',1.2)
    plot(ax2,t_rel_xf,heading_plot,'Color',heading_color,'LineWidth',1.2)
    xlim(ax2,[0,snippet_s]); ylim(ax2,[-pi,pi])
    xlabel(ax2,'time (s)'); ylabel(ax2,'angle (rad)')
    legend(ax2,{'bump position (mu, smoothed)','fly heading (-cue, smoothed)'},'Location','best')
    title(ax2,'bump position vs. fly heading')

    %% row 3: heading velocity (black, cue-derived, "speed relative to the
    % VR environment") and bump velocity (EPG blue) traces over the same
    % snippet -- SMOOTHED but NOT lag-shifted for display (the lag found by
    % sweep 3/3 is applied only to the fit in row 4, so these traces are
    % shown on their own shared, real fictrac time axis rather than one
    % shifted relative to the other).
    ax3 = nexttile(t,4+c); hold(ax3,'on')
    plot(ax3,t_rel_xf,E.heading_vel_disp(E.xf_idx),'Color',heading_color,'LineWidth',1.2)
    plot(ax3,t_rel_xf,E.bump_vel_disp(E.xf_idx),'Color',bump_color,'LineWidth',1.2)
    yline(ax3,0,':','Color',[.6 .6 .6])
    xlim(ax3,[0,snippet_s])
    xlabel(ax3,'time (s)'); ylabel(ax3,'angular velocity (rad/s)')
    legend(ax3,{'fly heading velocity (cue-derived)','bump velocity'},'Location','best')
    title(ax3,sprintf('bump velocity lags heading velocity by %d frames (%.2fs) -- not shifted here for display', ...
        lag_frames_opt,lag_frames_opt*median(diff(E.xf))))

    %% row 4: scatter of fly heading velocity (cue-derived) vs. bump
    % velocity for the WHOLE EXAMPLE TRIAL (not just the snippet), using the
    % lag-shifted, validity-gated pair the winning pipeline actually fits --
    % i.e. exactly what ranked this trial among empty>kir closed-loop trials.
    ax4 = nexttile(t,6+c); hold(ax4,'on')
    scatter(ax4,E.hv_full(E.vd_full),E.bv_full(E.vd_full),10,heading_color,'filled','MarkerFaceAlpha',.2)
    plot(ax4,xl,[0,0],':','Color',[.6,.6,.6]); plot(ax4,[0,0],yl,':','Color',[.6,.6,.6])
    plot(ax4,xl,xl*trial_gain(i),'r','LineWidth',1.5)
    xlim(ax4,xl); ylim(ax4,yl); axis(ax4,'square')
    xlabel(ax4,'fly heading velocity, cue-derived (rad/s)'); ylabel(ax4,'bump velocity (rad/s)')
    title(ax4,sprintf('whole trial: gain=%.2f, r=%.2f, n=%d valid samples (of %d total)', ...
        trial_gain(i),trial_r(i),trial_n_valid(i),numel(E.hv_full)))
    text(ax4,xl(2),yl(1),sprintf('gain=%.2f\nr=%.2f',trial_gain(i),trial_r(i)), ...
        'HorizontalAlignment','right','VerticalAlignment','bottom')
end

%% row 5: genotype x light-condition group comparison, one point per fly --
% left panel: gain, right panel: correlation coefficient (r). Per explicit
% request, ordered condition-major (empty CL, lpsp CL, empty dark, lpsp
% dark) -- NOT lpsp_kir_claude.m's own genotype-major order (empty CL,
% empty dark, lpsp CL, lpsp dark) -- so the two closed-loop groups sit
% together, then the two dark groups. Genotype keeps a consistent color
% across both light conditions (gray=empty, magenta=lpsp -- lpsp's own
% color matches Figure_1_claude.m's "LPsP > syt7f" column convention);
% light condition is distinguished by grouped x-position, not color.
% group_fly_gain/group_fly_r (built above) are already in exactly this
% 1=emptyCL,2=lpspCL,3=emptyDark,4=lpspDark order, so no re-indexing is
% needed here.
group_cat_labels = {'empty (CL)','lpsp (CL)','empty (dark)','lpsp (dark)'};
empty_color = [.5 .5 .5]; lpsp_color = [.7 0 .7];
group_cat_colors = [empty_color; lpsp_color; empty_color; lpsp_color];

cat_x_fly    = [];
all_fly_gain = [];
all_fly_r    = [];
for g = 1:4
    n_g = numel(group_fly_gain{g});
    cat_x_fly    = [cat_x_fly; g*ones(n_g,1)]; %#ok<AGROW>
    all_fly_gain = [all_fly_gain; group_fly_gain{g}]; %#ok<AGROW>
    all_fly_r    = [all_fly_r; group_fly_r{g}]; %#ok<AGROW>
end

ax5a = nexttile(t,9); %#ok<NASGU>
groupplot(cat_x_fly,all_fly_gain,group_cat_labels,group_cat_colors)
ylabel('gain (bump vel. ~ 1 + heading vel.)')
title('gain, by genotype and light condition (one point per fly)')

ax5b = nexttile(t,10); %#ok<NASGU>
groupplot(cat_x_fly,all_fly_r,group_cat_labels,group_cat_colors)
ylabel('correlation coefficient (r)')
title('heading/bump velocity correlation, by genotype and light condition (one point per fly)')

sgtitle(t,sprintf('Figure 4: bump tracking examples, empty>kir vs. lpsp>kir (bump=%.2fs, heading=%.2fs, lag=%d frames smoothing pipeline)', ...
    mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt))

%% export
drawnow % force a full render pass before export -- headless (-nodisplay) batch runs have occasionally clipped a tick label's leading characters otherwise
if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig, fullfile(export_dir,'Figure_4_claude.png'), 'Resolution', 300)
if ~isfolder(all_figs_dir); mkdir(all_figs_dir); end
exportgraphics(fig, fullfile(all_figs_dir,'Fig4_V1.pdf'), 'ContentType', 'auto')

%% ===================== SCRATCH: gallery of ALL candidate trials per genotype =====================
% per explicit request, same spirit as Figure_3_claude.m's own "gallery of
% candidate example stims" scratch section: every closed-loop candidate
% trial for BOTH genotypes, laid out as one ranked column each (rank 1 =
% closest to that genotype's own fly-level mean gain/r, i.e. the automatic
% pick used in Figure 4 above; further down = progressively less
% representative), so a specific rank/trial can be picked by eye instead
% of always taking the automatic closest-to-mean pick. Each panel: that
% trial's own best 60s snippet (same selection as the main figure), heatmap
% with smoothed bump position (white) and heading (orange, chosen for
% contrast against both the teal colormap and the white mu line -- unlike
% the main figure's black, which would be hard to see here) overlaid
% directly on it, so tracking quality can be judged at a glance without
% needing 4 separate rows per candidate.
n_gallery_rows = max(cellfun(@numel,gallery));
fig_gallery = figure('color','w','Position',[50 50 1500 170*n_gallery_rows]); clf
tg = tiledlayout(fig_gallery,n_gallery_rows,2,'TileSpacing','compact','Padding','compact');

gallery_heading_color = [1 0.6 0]; % orange

for col = 1:2
    rows_ranked = gallery{col};
    z_ranked    = gallery_z{col};
    for r = 1:n_gallery_rows
        ax = nexttile(tg,(r-1)*2+col); hold(ax,'on')
        if r > numel(rows_ranked)
            axis(ax,'off'); continue
        end
        i = rows_ranked(r);
        G = build_example_display(all_data(i),mu_smooth_opt_s,heading_smooth_opt_s,lag_frames_opt, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max,snippet_s,snippet_step_s,min_window_n);

        z_snip = G.ex_im.z(:,G.xb_idx);
        z_clim = [prctile(min(z_snip,[],1),5), prctile(max(z_snip,[],1),95)];
        t_rel_xb = G.xb(G.xb_idx) - G.snippet_t0;
        t_rel_xf = G.xf(G.xf_idx) - G.snippet_t0;

        imagesc(ax,t_rel_xb,unwrap(G.ex_im.alpha),z_snip,z_clim)
        colormap(ax,white_to_color(bump_color))

        mu_plot = wrap_to_pi(G.mu_smoothed_disp(G.xb_idx));
        mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan;
        heading_plot = wrap_to_pi(G.heading_smoothed_disp(G.xf_idx));
        heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;
        plot(ax,t_rel_xb,mu_plot,'w','LineWidth',1)
        plot(ax,t_rel_xf,heading_plot,'Color',gallery_heading_color,'LineWidth',1)

        xlim(ax,[0,snippet_s]); ylim(ax,[min(unwrap(G.ex_im.alpha)),max(unwrap(G.ex_im.alpha))])
        pick_str = ''; if r==1; pick_str = ' <-- CURRENT PICK'; end
        title(ax,sprintf('rank %d: fly %s, trial "%s"\ngain=%.2f, r=%.2f, z-dist=%.2f%s', ...
            r,fly_short_label(fly_id{i}),trial_short_label(all_data(i).meta),trial_gain(i),trial_r(i),z_ranked(r),pick_str), ...
            'FontSize',7,'Interpreter','none')
        if col==1; ylabel(ax,'EB angle (rad)'); end
        if r==n_gallery_rows; xlabel(ax,'time (s)'); end
    end
end
title(tg,sprintf('candidate gallery for Figure 4 -- %s (left, n=%d) vs. %s (right, n=%d), ranked by closeness to that genotype''s own fly-level mean gain/r', ...
    example_geno_label{1},numel(gallery{1}),example_geno_label{2},numel(gallery{2})),'Interpreter','none')

exportgraphics(fig_gallery, fullfile(export_dir,'Figure_4_claude_scratch_gallery.png'), 'Resolution', 150)

%% ===================== functions =====================
% trial_folder_name/trial_fly_id/trial_short_label/trial_gain_vectors_v2/
% fly_gain_v2/groupplot are ported VERBATIM from data scripts/
% lpsp_kir_claude.m (see that script for the reasoning behind each design
% choice, summarized in this file's header above); get_xb/white_to_color/
% wrap_to_pi/circular_coverage are ported verbatim from Figure_1_claude.m /
% Figure_2_claude.m. build_example_display, fly_gain_corr_v2,
% lag_trim_time, and pick_wide_heading_snippet are new, written for this figure.

function n = trial_number(meta_path)
    % the integer right after the date-dash in the trial folder name, e.g.
    % "20231109-3_EPG_7f_empty_kir" -> 3 -- ported from
    % lpsp_kir_script_minimal.m's own light-condition fallback (regex out
    % every digit run in the trial folder name and take the second one:
    % {"20231109","3","7"} for that example). NaN if the folder name
    % doesn't match the expected "<8-digit date>-<trial num>..." pattern.
    tname = trial_folder_name(meta_path);
    m = regexp(tname,'^\d{8}-(\d+)','tokens','once');
    if isempty(m)
        n = nan;
    else
        n = str2double(m{1});
    end
end

function name = trial_folder_name(meta_path)
    % ported from lpsp_kir_claude.m, with ONE deliberate change: that
    % script splits meta_path on filesep, which is a safe assumption on
    % its own Windows-only workflow (meta stores literal "Z:\pablo\..."
    % paths) but silently fails on a Mac, where filesep is '/' and never
    % appears in these paths at all -- confirmed directly here (every
    % trial's fly grouping collapsed to one trial per "fly" the first time
    % this ran, because the un-split path never matched the '^fly\s*\d+$'
    % pattern below). split_path (same fix Figure_1_claude.m /
    % Figure_2_claude.m already use for this exact reason) splits on
    % either separator explicitly instead.
    parts = split_path(meta_path);
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        name = parts{end-1};
    else
        name = parts{end};
    end
end

function fid = trial_fly_id(meta_path)
    % ported from lpsp_kir_claude.m -- see trial_folder_name above for why
    % split_path (not strsplit(...,filesep)) is used here on a Mac.
    parts = split_path(meta_path);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
    if isempty(fly_part)
        fid = meta_path;
    else
        fid = strjoin(parts(1:fly_part(1)),filesep);
    end
end

function parts = split_path(p)
    % ported from Figure_1_claude.m / Figure_2_claude.m: meta paths are
    % literal Windows paths regardless of the host OS this script runs on
    % -- split on either separator explicitly.
    parts = strsplit(p,{'\\','/'});
end

function lbl = fly_short_label(fly_path)
    % ported from lpsp_kir_claude.m: "<date> fly N" from a fly_id path like
    % ...\<date folder>\fly N (split_path, not strsplit(...,filesep) --
    % same Mac-compatibility fix as trial_folder_name/trial_fly_id above).
    parts = split_path(fly_path);
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function lbl = trial_short_label(meta_path)
    % ported from lpsp_kir_claude.m
    tname = trial_folder_name(meta_path);
    m = regexp(tname,'^(\d{8}-\d+)','match','once');
    if isempty(m)
        lbl = tname;
    else
        lbl = m;
    end
end

function groupplot(cat_x, values, cat_labels, colors)
    % jittered per-point scatter + mean +/- SEM errorbar per category,
    % ported verbatim from lpsp_kir_claude.m (originally copied there from
    % lpsp_compartments_claude_script.m) -- plots into the CURRENT axes.
    hold on
    n_cat = numel(cat_labels);
    for c = 1:n_cat
        y = values(cat_x==c);
        y = y(~isnan(y));
        if isempty(y)
            continue
        end
        jitter = (rand(size(y))-.5)*.3;
        scatter(c+jitter,y,20,colors(c,:),'filled','MarkerFaceAlpha',.3)
        errorbar(c,mean(y),std(y)/sqrt(numel(y)),'o','Color',colors(c,:)*.6, ...
            'MarkerFaceColor',colors(c,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    xticks(1:n_cat); xticklabels(cat_labels)
    xlim([0.5,n_cat+0.5])
    y_lims = ylim;
    for c = 1:n_cat
        n = sum(cat_x==c & ~isnan(values));
        text(c,y_lims(1),sprintf('n=%d',n),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    end
    plot(xlim,[0,0],':k')
end

function E = build_example_display(trial,mu_smooth_s,heading_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max,snippet_s,snippet_step_s,min_window_n)
    % everything needed to display ONE trial as one example/gallery panel:
    % its own best snippet_s-second window (the one whose discernible-bump
    % samples span the WIDEST range of heading positions -- see
    % pick_wide_heading_snippet -- so the displayed bump is shown tracking
    % through many different headings, not just wherever the trial happens
    % to start) plus the smoothed display traces over the WHOLE trial
    % (cropped to the snippet by the caller, as needed for each panel).
    % Shared by both the two main example columns and the candidate
    % gallery below, so both use the exact same per-trial computation.
    ex_ft = trial.ft;
    ex_im = trial.im;
    xf = ex_ft.xf(:);
    dt_ft = median(diff(xf));

    xb = get_xb(ex_ft,size(ex_im.z,2));
    dt_im = median(diff(xb));
    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed_disp = smoothdata(unwrap(ex_im.mu(:)),'gaussian',win_im);
    win_hd = max(1,round(heading_smooth_s/dt_ft));
    heading_smoothed_disp = smoothdata(unwrap(-ex_ft.cue(:)),'gaussian',win_hd);
    heading_vel_disp = gradient(heading_smoothed_disp)/dt_ft;
    bump_vel_disp    = gradient(interp1(xb,mu_smoothed_disp,xf,'linear','extrap'))/dt_ft;

    [hv_full,bv_full,vd_full] = trial_gain_vectors_v2(trial,ex_im.mu,ex_im.rho, ...
        mu_smooth_s,heading_smooth_s,lag_frames,vel_thresh,bump_thresh,rho_thresh,vel_max);
    t_aligned = lag_trim_time(xf,lag_frames); % same length/pairing as hv_full/bv_full/vd_full
    heading_pos_aligned = lag_trim_time(wrap_to_pi(heading_smoothed_disp),lag_frames); % same trim convention as heading_vel, so it lines up index-for-index with vd_full

    [snippet_t0,snippet_t1] = pick_wide_heading_snippet(t_aligned,heading_pos_aligned,vd_full,snippet_s,snippet_step_s,min_window_n);

    E.ex_im = ex_im;
    E.xb = xb; E.xf = xf;
    E.xb_idx = xb>=snippet_t0 & xb<=snippet_t1;
    E.xf_idx = xf>=snippet_t0 & xf<=snippet_t1;
    E.snippet_t0 = snippet_t0; E.snippet_t1 = snippet_t1;
    E.mu_smoothed_disp = mu_smoothed_disp;
    E.heading_smoothed_disp = heading_smoothed_disp;
    E.heading_vel_disp = heading_vel_disp;
    E.bump_vel_disp = bump_vel_disp;
    E.hv_full = hv_full; E.bv_full = bv_full; E.vd_full = vd_full;
end

function [heading_vel,bump_vel,valid] = trial_gain_vectors_v2(trial,mu_raw,rho_raw,mu_smooth_s,heading_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % ported verbatim from lpsp_kir_claude.m -- see this file's header
    % comment (steps 1-5) for the full explanation of every stage here.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu_raw);
    xb = linspace(xf(1),xf(end),n_im)';
    dt_im = (xf(end)-xf(1)) / (n_im-1);

    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    bump_vel_full = gradient(interp1(xb,mu_smoothed,xf,'linear','extrap'))/dt;
    rho_full = interp1(xb,rho_raw(:),xf,'linear','extrap');

    win_hd = max(1,round(heading_smooth_s/dt));
    heading_smoothed = smoothdata(unwrap(-trial.ft.cue(:)),'gaussian',win_hd);
    heading_vel_full = gradient(heading_smoothed)/dt;

    if lag_frames == 0
        heading_vel = heading_vel_full; bump_vel = bump_vel_full; rho_i = rho_full;
    elseif lag_frames > 0
        heading_vel = heading_vel_full(1:end-lag_frames);
        bump_vel    = bump_vel_full(lag_frames+1:end);
        rho_i       = rho_full(lag_frames+1:end);
    else
        heading_vel = heading_vel_full(-lag_frames+1:end);
        bump_vel    = bump_vel_full(1:end+lag_frames);
        rho_i       = rho_full(1:end+lag_frames);
    end

    valid = abs(heading_vel) > vel_thresh & abs(heading_vel) < vel_max & abs(bump_vel) < bump_thresh & rho_i > rho_thresh;
end

function gain = fly_gain_v2(all_data,trial_list,mu_smooth_s,heading_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % ported verbatim from lpsp_kir_claude.m
    heading_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        trial = all_data(trial_list(k));
        [hv,bv,vd] = trial_gain_vectors_v2(trial,trial.im.mu,trial.im.rho,mu_smooth_s,heading_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        heading_vel = [heading_vel; hv]; %#ok<AGROW>
        bump_vel    = [bump_vel; bv]; %#ok<AGROW>
        valid       = [valid; vd]; %#ok<AGROW>
    end

    min_frac_moving = 0.01;
    n_valid = sum(valid);
    if n_valid < 50 || mean(abs(heading_vel) > vel_thresh) < min_frac_moving
        gain = nan;
        return
    end

    b = [ones(n_valid,1),heading_vel(valid)] \ bump_vel(valid);
    gain = b(2);
end

function [gain,r] = fly_gain_corr_v2(all_data,trial_list,mu_smooth_s,heading_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % same pooling/gate as fly_gain_v2 above, extended to also return r --
    % used for the per-fly genotype-mean reference (not the fly_gain_v2
    % smoothing/lag sweep itself, which only ever needed gain).
    heading_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        trial = all_data(trial_list(k));
        [hv,bv,vd] = trial_gain_vectors_v2(trial,trial.im.mu,trial.im.rho,mu_smooth_s,heading_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        heading_vel = [heading_vel; hv]; %#ok<AGROW>
        bump_vel    = [bump_vel; bv]; %#ok<AGROW>
        valid       = [valid; vd]; %#ok<AGROW>
    end

    min_frac_moving = 0.01;
    n_valid = sum(valid);
    if n_valid < 50 || mean(abs(heading_vel) > vel_thresh) < min_frac_moving
        gain = nan; r = nan;
        return
    end

    b = [ones(n_valid,1),heading_vel(valid)] \ bump_vel(valid);
    gain = b(2);
    r = corr(heading_vel(valid),bump_vel(valid));
end

function xb = get_xb(ft, n_im)
    % ported from Figure_1_claude.m / Figure_2_claude.m
    if isfield(ft,'xb') && numel(ft.xb) == n_im
        xb = ft.xb(:);
    else
        xb = linspace(ft.xf(1),ft.xf(end),n_im)';
    end
end

function cmap = white_to_color(max_color)
    % ported verbatim from Figure_1_claude.m / Figure_2_claude.m
    n = 256;
    cmap = [linspace(1,max_color(1),n)', linspace(1,max_color(2),n)', linspace(1,max_color(3),n)'];
end

function y = wrap_to_pi(x)
    % ported from Figure_1_claude.m (originally lpsp_compartments_claude_script.m)
    y = mod(x+pi,2*pi) - pi;
end

function t_trim = lag_trim_time(xf, lag_frames)
    % the real-time axis (a subset of xf) that trial_gain_vectors_v2's
    % shift-and-trim convention leaves heading_vel/bump_vel/valid aligned
    % to -- heading_vel(k) is "heading at this trimmed time", paired with
    % "bump velocity lag_frames later", so there is no single real instant
    % shared by both; xf(1:end-lag) (the HEADING side's own real time) is
    % used as the common index for plotting/windowing purposes.
    if lag_frames == 0
        t_trim = xf;
    elseif lag_frames > 0
        t_trim = xf(1:end-lag_frames);
    else
        t_trim = xf(-lag_frames+1:end);
    end
end

function [t0,t1] = pick_wide_heading_snippet(t_aligned, heading_pos, valid, snippet_s, snippet_step_s, min_window_n)
    % among snippet_s-second windows of t_aligned, drop any with fewer than
    % min_window_n valid (discernible-bump) samples, then take the window
    % whose valid samples span the WIDEST range of heading positions --
    % same "circular coverage" objective as Figure_1_claude.m's own
    % pick_snippet (2*pi minus the largest gap between sorted heading
    % samples around the circle), so the displayed snippet shows the bump
    % tracking through many different headings rather than oscillating
    % narrowly around one, or just happening to be wherever the trial
    % starts. Per explicit request, this REPLACES an earlier version of
    % this function that instead picked the highest-local-correlation
    % window -- that one tended to pick early, narrow-heading-range windows
    % whenever the trial's tracking happened to be locally cleanest there.
    if t_aligned(end) - t_aligned(1) <= snippet_s
        t0 = t_aligned(1); t1 = t_aligned(end);
        return
    end

    starts = t_aligned(1):snippet_step_s:(t_aligned(end)-snippet_s);
    cand_start = []; cand_cov = [];
    for w = starts
        idx = t_aligned>=w & t_aligned<w+snippet_s & valid;
        if sum(idx) < min_window_n
            continue
        end
        cand_start(end+1) = w; %#ok<AGROW>
        cand_cov(end+1)   = circular_coverage(heading_pos(idx)); %#ok<AGROW>
    end

    if isempty(cand_start)
        t0 = starts(1); t1 = t0 + snippet_s;
        return
    end

    [~,rel] = max(cand_cov);
    t0 = cand_start(rel);
    t1 = t0 + snippet_s;
end

function cov = circular_coverage(theta)
    % ported from Figure_1_claude.m: how much of the full circle (radians)
    % is spanned by these samples, measured as 2*pi minus the single
    % largest gap between consecutive angles once sorted around the circle
    % -- close to 2*pi means the samples are spread all the way around;
    % close to 0 means they're all clustered in one small arc.
    if isempty(theta)
        cov = 0;
        return
    end
    th = sort(mod(theta(:),2*pi));
    gaps = diff([th; th(1)+2*pi]);
    cov = 2*pi - max(gaps);
end
