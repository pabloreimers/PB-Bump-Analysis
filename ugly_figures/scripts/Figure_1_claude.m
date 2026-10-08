%% Figure_1_claude
% 3-column summary figure, one column per indicator/genotype:
%   col 1: EPG > GCaMP        (data/hackathon_20250729.mat, via hackathon_claude.m's recipe)
%   col 2: LPsP > syt7f       (data/lpsp_cl_data_20240206.mat, via lpsp_compartments_claude_script.m's recipe)
%   col 3: EPG > GRAB(DA2m)   (data/lpsp_cl_data_20240206.mat, via lpsp_compartments_claude_script.m's recipe)
%
% Rows:
%   1) example fly, example closed-loop trial: im.z heatmap with bump
%      position (im.mu, white) and fly heading (-ft.cue, colored) overlaid.
%   2) same example trial: im.mu and -ft.cue traces overlaid, line-only.
%   3) per-fly correlation summary, closed loop only, 2 x-tick categories:
%      circular position correlation (mu vs. -cue) and Pearson velocity
%      correlation (bump vel vs. fly rotation). One dot per fly. (Dark
%      trials are still processed below but no longer plotted -- same for
%      the row-5 co-imaging panels.)
%
% Loads the SAME three source .mat files combined by
% lpsp_compartments_claude_script.m (lpsp_cl, lpsp_cl_redo, epg_dlight),
% but only lpsp_cl actually contains syt7f/GRAB(DA2m) trials (lpsp_cl_redo
% is syt8m, epg_dlight is dLight) -- those other two are loaded for parity
% with that script's dataset set and then simply unused here, since neither
% of this figure's two lpsp columns wants them.
%
% ==THIS SCRIPT'S KEY SIMPLIFICATION vs. lpsp_compartments_claude_script.m==
% That script cross-checks closed-loop/dark against trialSettings.csv on
% the Z:\pablo raw-data drive (folder names disagree with the logged
% pattern for a real subset of lpsp_cl trials). That drive isn't mounted
% here, so this script falls back to the folder-name convention alone
% (meta contains '_dark' -> dark) for every lpsp_cl/lpsp_cl_redo trial --
% the same fallback that script itself uses when no CSV is found. This
% will mislabel whatever small subset of trials that script's CSV
% cross-check would have caught; revisit with drive access if that matters.
%
% ==ROW 3 METRICS AND LAGS==
% Ported from lpsp_compartments_claude_script.m section 7's own two per-fly
% metrics (figures 40 and 43 there): circular position correlation (mean
% of many 30s-window circ_corrcc(mu,-cue) values) and Pearson velocity
% correlation (corr(bump_vel, fly r_speed) over discernable-bump samples),
% one dot per fly, closed loop only. Fluorescence-vs-behavior lag is
% chosen SEPARATELY per indicator and per metric (6 lags total), as the
% lag in [-1,+1] s that maximizes the population-mean of that metric
% across flies -- see the "per-indicator, per-metric lag" section.
%
% ==HACKATHON (EPG>GCaMP) ADAPTATION==
% The hackathon dataset's berg4 trials change gain/luminance mid-trial (a
% staircase of ~7.5s constant-gain holds -- see hackathon_claude.m), so
% "the trial" isn't a well-defined single condition there. This script
% reuses hackathon_claude.m's own block-segmentation logic (find_runs on
% (gain, cue_brightness), with a dark block inheriting the most recently
% active closed-loop gain for that fly) to chunk every trial into
% constant-gain-or-dark blocks, then keeps only gain=0.8 blocks (both
% closed-loop and dark, per the user's request). Each kept block is one
% "unit" put through the same unit_metrics() as an lpsp trial: 30s sliding
% windows for the position correlation where the block is long enough,
% otherwise the whole block as a single window (most gain=0.8 berg4
% blocks are only ~7.5s long), plus one pool of (bump_vel, r_speed)
% samples for that fly's velocity correlation.
% Per hackathon_claude.m's own finding, the fly's rotation direction
% relative to bump motion is sign-FLIPPED on the berg4 rig relative to the
% older G4-pattern rig -- r_speed is sign-corrected for that (fly_sign),
% while -cue's convention is already consistent across both rigs (no
% correction needed there), exactly as documented in that script's
% trial_velocity_gain_signals().

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
addpath(fullfile(repo_root,'circ_stats'))
data_dir = fullfile(repo_root,'.data'); % processed datasets live in the gitignored .data/ (not data/, which only holds the small LED-stim summaries)
export_dir = fullfile(repo_root,'ugly_figures','exports');

%% shared analysis constants (ported from lpsp_compartments_claude_script.m / hackathon_claude.m)
flash_mad_thresh  = 8;    % detect_flash_frames: brightness-outlier MAD multiple
f0_pct            = 7;    % dF/F baseline percentile (hackathon bump recompute)
smooth_frames_f   = 5;    % hackathon: moving-average window on im.f before dF/F

rot_thresh_track  = 0.5;  % rad/s, "fly is rotating" (discernable-bump criterion)
rho_thresh_track  = 0.2;  % im.rho, "bump is discernable"
bump_vel_thresh   = 10;   % rad/s, exclude gradient/unwrap artifacts at wrap points
min_window_n      = 50;   % minimum discernable-bump samples to trust one circ_corrcc chunk
min_vel_n         = 100;  % minimum pooled samples to trust one fly's velocity Pearson r

window_s          = 30;   % s, lpsp sliding-window length for position correlation
window_step_s     = 10;   % s, lpsp sliding-window step

lag_grid_s        = -1:0.05:1; % s, candidate fluorescence-vs-behavior lags swept per indicator x metric (see the lag section below)

%% load the three lpsp_compartments_claude_script.m source datasets and combine
dataset_names = {'lpsp_cl','lpsp_cl_redo','epg_dlight'};
source_files  = {'lpsp_cl_data_20240206.mat','lpsp_cl_redo_data_20240306.mat','epg_dlight_20260415.mat'};

lpsp_sets = cell(1,numel(dataset_names));
for k = 1:numel(dataset_names)
    tmp = load(fullfile(data_dir,source_files{k}),'all_data');
    d = tmp.all_data(:);
    for i = 1:numel(d)
        d(i).dataset = dataset_names{k};
    end
    lpsp_sets{k} = d;
end
lpsp_data = combine_datasets(lpsp_sets{:});
fprintf('loaded %d combined lpsp trials (lpsp_cl + lpsp_cl_redo + epg_dlight)\n', numel(lpsp_data));

%% ===================== column 1 source (EPG calcium indicator) =====================
% col1_source selects which EPG calcium dataset fills the first column; the
% rest of the script only sees the C1 struct (trials, flash frames, fly ids,
% units) so the two sources are interchangeable downstream:
%   'hackathon' - .data/hackathon_20250729.mat: berg4 rig, EPG>syt8m,
%                 gain-staircase trials chunked into constant-gain blocks,
%                 gain=0.8 blocks kept (see prep_hackathon_column)
%   'epg7f'     - .data/epg_7f_20261007.mat: gain_change rig, 16 EPG>GCaMP7f
%                 flies, whole 300 s trials at closed-loop gain 0.7 (bar or
%                 dark), 100x256 px imaging (see prep_epg7f_column)
col1_source = 'epg7f';
switch col1_source
    case 'hackathon'
        C1 = prep_hackathon_column(data_dir, flash_mad_thresh, f0_pct, smooth_frames_f);
    case 'epg7f'
        C1 = prep_epg7f_column(data_dir, flash_mad_thresh);
    otherwise
        error('unknown col1_source "%s"', col1_source)
end
% downstream code keeps the historical hack_* names for column 1
hack_data        = C1.data;
hack_flash       = C1.flash;
hack_fly         = C1.fly;
hack_unit_trial  = C1.unit_trial;
hack_unit_s      = C1.unit_s;
hack_unit_e      = C1.unit_e;
hack_unit_fly    = C1.unit_fly;
hack_unit_isdark = C1.unit_isdark;

%% ===================== process lpsp (LPsP>syt7f, EPG>GRAB(DA2m)) =====================

n_lpsp = numel(lpsp_data);
lpsp_is_dark = false(n_lpsp,1);
lpsp_indicator = cell(n_lpsp,1);
lpsp_fly = cell(n_lpsp,1);
lpsp_flash = cell(n_lpsp,1);

for i = 1:n_lpsp
    meta = lpsp_data(i).meta;
    if isempty(meta) || ~ischar(meta)
        continue % leftover empty placeholder trial
    end

    if isfield(lpsp_data(i).ft,'pattern') && ~isempty(lpsp_data(i).ft.pattern)
        lpsp_is_dark(i) = contains(char(lpsp_data(i).ft.pattern),'background','IgnoreCase',true);
    else
        % no network-drive trialSettings.csv access here -- fall back to the
        % folder-name convention (see header comment for the caveat)
        lpsp_is_dark(i) = contains(meta,'_dark','IgnoreCase',true);
    end

    lpsp_indicator{i} = trial_indicator(meta);
    lpsp_fly{i}       = trial_fly_id(lpsp_data(i).dataset,meta);

    lpsp_flash{i} = detect_flash_frames(lpsp_data(i).im.f, flash_mad_thresh);
end

%% figure columns: only syt7f and GRAB(DA2m) are used from this combined dataset
lpsp_groups = {'syt7f','GRAB(DA2m)'};

% one "unit" per lpsp trial (the whole trial); which column it feeds
lpsp_unit_fly      = cell(n_lpsp,1);
lpsp_unit_isdark   = false(n_lpsp,1);
lpsp_unit_group    = zeros(n_lpsp,1); % 1=syt7f, 2=GRAB(DA2m), 0=neither (unused)
for i = 1:n_lpsp
    g = find(strcmp(lpsp_groups,lpsp_indicator{i}),1);
    if isempty(g)
        continue
    end
    lpsp_unit_group(i)  = g;
    lpsp_unit_fly{i}    = lpsp_fly{i};
    lpsp_unit_isdark(i) = lpsp_is_dark(i);
end
for g = 1:numel(lpsp_groups)
    fprintf('%s: n=%d trials\n', lpsp_groups{g}, sum(lpsp_unit_group==g));
end

%% ===================== assemble the 3 figure columns =====================
% per column: its units (trial + [start,end] index on ft.xf), fly ids,
% light condition, which flash-frame list / rotation-sign rule applies.
% The per-unit tracking metrics (wincorr, bumpvel, rspeed, offset,
% nbumpok) are filled in by the lag section right after this.
group_defs = struct( ...
    'title',    {C1.title,'LPsP > syt7f','EPG > GRAB(DA2m)'}, ...
    'unit_fly',    {hack_unit_fly,    lpsp_unit_fly(lpsp_unit_group==1),    lpsp_unit_fly(lpsp_unit_group==2)}, ...
    'unit_isdark', {hack_unit_isdark, lpsp_unit_isdark(lpsp_unit_group==1), lpsp_unit_isdark(lpsp_unit_group==2)}, ...
    'unit_trial',  {hack_unit_trial,  find(lpsp_unit_group==1),             find(lpsp_unit_group==2)}, ...
    'unit_s',      {hack_unit_s,      ones(sum(lpsp_unit_group==1),1),      ones(sum(lpsp_unit_group==2),1)}, ...
    'unit_e',      {hack_unit_e,      [],                                   []}, ...
    'all_data',    {hack_data,        lpsp_data,                            lpsp_data}, ...
    'flash',       {hack_flash,       lpsp_flash,                           lpsp_flash}, ...
    'use_berg4_sign', {C1.use_berg4_sign, false,                           false} ... % hackathon berg4 trials need r_speed sign-flipped (see header)
);
% lpsp units' block end = the whole trial (numel of that trial's own ft.xf)
for gi = 2:3
    rows = group_defs(gi).unit_trial;
    e_vals = zeros(numel(rows),1);
    for r = 1:numel(rows)
        e_vals(r) = numel(group_defs(gi).all_data(rows(r)).ft.xf);
    end
    group_defs(gi).unit_e = e_vals;
end

%% per-indicator, per-metric lag: maximize the POPULATION mean correlation
% Two lags per column, each picked by sweeping lag_grid_s (-1..+1 s) and
% taking the lag at which the mean across closed-loop flies of that metric
% is largest:
%   lag_pos_s(gi) - for the windowed circular position correlation
%   lag_vel_s(gi) - for the bump-velocity vs. rotation Pearson r
% Positive lag = fluorescence shifted EARLIER, i.e. fluorescence lags
% behavior (see lagged_track_signals). This replaces the earlier single
% per-indicator lag fit on |r_speed| vs. mean fluorescence -- a different
% signal from either metric, suspected of under-aligning the velocity
% metric for EPG>GCaMP. unit_metrics() computes both metrics identically
% for all three columns (sliding window_s windows for the position
% correlation; a unit shorter than one window, as most hackathon gain=0.8
% blocks are, is a single window).
P = struct('rot_thresh_track',rot_thresh_track,'rho_thresh_track',rho_thresh_track,'bump_vel_thresh',bump_vel_thresh, ...
           'min_window_n',min_window_n,'window_s',window_s,'window_step_s',window_step_s);
lag_pos_s = nan(1,3); lag_vel_s = nan(1,3);
lag_sweep_pos = nan(3,numel(lag_grid_s)); % population-mean position correlation at each candidate lag, per column
lag_sweep_vel = nan(3,numel(lag_grid_s)); % population-mean velocity correlation at each candidate lag, per column
for gi = 1:3
    G = group_defs(gi);
    for L = 1:numel(lag_grid_s)
        U = unit_metrics(G.all_data,G.flash,G.unit_trial,G.unit_s,G.unit_e,G.unit_isdark,lag_grid_s(L),G.use_berg4_sign,P);
        [pv,vv] = per_fly_cl(U,G.unit_fly,G.unit_isdark,min_vel_n);
        lag_sweep_pos(gi,L) = mean(pv,'omitnan');
        lag_sweep_vel(gi,L) = mean(vv,'omitnan');
    end
    [~,bp] = max(lag_sweep_pos(gi,:));
    [~,bv] = max(lag_sweep_vel(gi,:));
    lag_pos_s(gi) = lag_grid_s(bp);
    lag_vel_s(gi) = lag_grid_s(bv);
    fprintf('%s: position lag %+.2f s (pop. mean r=%.2f), velocity lag %+.2f s (pop. mean r=%.2f)\n', ...
        G.title, lag_pos_s(gi), lag_sweep_pos(gi,bp), lag_vel_s(gi), lag_sweep_vel(gi,bv));

    % final per-unit metrics, each at its own metric's lag
    U_pos = unit_metrics(G.all_data,G.flash,G.unit_trial,G.unit_s,G.unit_e,G.unit_isdark,lag_pos_s(gi),G.use_berg4_sign,P);
    U_vel = unit_metrics(G.all_data,G.flash,G.unit_trial,G.unit_s,G.unit_e,G.unit_isdark,lag_vel_s(gi),G.use_berg4_sign,P);
    group_defs(gi).unit_wincorr = U_pos.wincorr;
    group_defs(gi).unit_offset  = U_pos.offset;
    group_defs(gi).unit_nbumpok = U_pos.nbumpok;
    group_defs(gi).unit_bumpvel = U_vel.bumpvel;
    group_defs(gi).unit_rspeed  = U_vel.rspeed;
end

%% pick, per column, the example fly + a 120s snippet showing good bump tracking
% "good bump tracking" = the same criterion used throughout: high rotation
% (lots of discernable-bump samples) and low circular variance of the
% mu-(-cue) offset. That picks the FLY (pooling all its own closed-loop
% units) and, within that fly's own best closed-loop unit, WHICH 120s of it
% to show (the snippet_s-second sub-window with the lowest offset
% variance among discernable-bump samples) -- one consistent metric for
% both choices, not two different ones.
snippet_s      = 60;  % s, requested example-trace window length, all 3 columns
snippet_step_s = 10;  % s, sliding step when searching for the best snippet_s-second sub-window

% manual override: force a specific fly as the example for a column instead
% of the automatic pick_best_fly() choice below (leave '' for automatic).
% GRAB(DA2m)'s auto-picked fly (20221205_1) didn't look like the best
% example on inspection -- 20221122_1 (its 20221122-1_EPG_GRABDA(2m)_1
% trial, "Trial 1") was flagged as a more promising one instead. Likewise
% syt7f's auto-picked fly (20221114_2) was swapped for 20221123_2 (its
% 20221123-7_LPsP_syt7f_cl_2 trial, "Trial 7"). EPG>GCaMP's auto-picked
% fly (8_np_6) was on a menotaxis bout in its auto-picked snippet, and the
% first replacement (5_np_4 / trial 27 / t0 = 33.8 s, chosen purely on
% walking statistics) turned out to have the dimmest raw PB in the dataset
% (50x128 px, 10-frame SNR ~1.5). A second scan ranking trials by raw PB
% image quality (mean-image contrast and 10-frame SNR) as well as walking
% statistics picked 1_np_1 / hackathon trial 7 / t0 = 540 s: one of only
% two flies recorded at 100x256 px, with the sharpest PB in the dataset
% (mean-image contrast 8.4, 10-frame SNR 4.3) and walking statistics close
% to the other two columns' (full-circle coverage, circ. std 1.14, mean
% |r| 0.44 rad/s), though with some heading holds at the rear of the arena
% (see forced_trial / forced_snippet_t0).
forced_fly = {'1_np_1','20221123_2','20221122_1'};

% manual override: force a specific trial (index into that column's
% all_data) as the example, instead of the "fly's best closed-loop unit"
% choice below (NaN = automatic). Only meaningful with forced_fly set.
forced_trial = [7, nan, nan];

% manual override: force a specific snippet start time (s, on that trial's
% own ft.xf clock) for a column instead of pick_snippet()'s automatic
% choice below (NaN = automatic).
forced_snippet_t0 = [540.0, nan, nan];
if ~strcmp(col1_source,'hackathon')
    % the column-1 overrides above refer to hackathon flies/trials; for any
    % other source fall back to the automatic picks
    forced_fly{1} = ''; forced_trial(1) = nan; forced_snippet_t0(1) = nan;
end

for gi = 1:3
    G = group_defs(gi);
    rows_cl = find(~G.unit_isdark);
    fly_list_cl = unique(G.unit_fly(rows_cl));

    if ~isempty(forced_fly{gi})
        chosen_fly = forced_fly{gi};
        assert(any(strcmp(fly_list_cl,chosen_fly)), ...
            'forced_fly "%s" has no closed-loop units in column %d (%s)', chosen_fly, gi, G.title)
    else
        rot_amount = nan(numel(fly_list_cl),1);
        off_var    = nan(numel(fly_list_cl),1);
        for f = 1:numel(fly_list_cl)
            sel = rows_cl(strcmp(G.unit_fly(rows_cl),fly_list_cl{f}));
            rot_amount(f) = sum(G.unit_nbumpok(sel));
            off_pool = cat(1,G.unit_offset{sel});
            if ~isempty(off_pool)
                off_var(f) = circ_var(off_pool);
            end
        end
        chosen_fly = pick_best_fly(fly_list_cl,rot_amount,off_var);
    end

    % among that fly's own closed-loop units, prefer ones long enough to
    % hold a full snippet_s-second window (a hackathon gain=0.8 berg4 hold
    % can be as short as ~7.5s -- see header comment); fall back to the
    % single longest available unit if none reach snippet_s.
    sel = rows_cl(strcmp(G.unit_fly(rows_cl),chosen_fly));
    unit_dur_s = nan(numel(sel),1);
    for u = 1:numel(sel)
        trial_u = G.unit_trial(sel(u));
        unit_dur_s(u) = G.all_data(trial_u).ft.xf(G.unit_e(sel(u))) - G.all_data(trial_u).ft.xf(G.unit_s(sel(u)));
    end
    long_enough = unit_dur_s >= snippet_s;
    if any(long_enough)
        candidates = sel(long_enough);
    else
        [~,longest] = max(unit_dur_s);
        candidates = sel(longest);
    end
    if ~isnan(forced_trial(gi))
        % restrict to this fly's closed-loop units in the forced trial; if
        % a snippet start is also forced, take the unit that contains it
        candidates = sel(G.unit_trial(sel)==forced_trial(gi));
        assert(~isempty(candidates), 'forced_trial %d has no closed-loop units for fly "%s" in column %d (%s)', forced_trial(gi), chosen_fly, gi, G.title)
        if ~isnan(forced_snippet_t0(gi))
            xf_ft = G.all_data(forced_trial(gi)).ft.xf;
            contains_t0 = xf_ft(G.unit_s(candidates)) <= forced_snippet_t0(gi) & xf_ft(G.unit_e(candidates)) >= forced_snippet_t0(gi)+snippet_s;
            if any(contains_t0)
                candidates = candidates(contains_t0);
            end
        end
    end
    [~,best_rel] = max(G.unit_nbumpok(candidates));
    example_unit = candidates(best_rel);

    i  = G.unit_trial(example_unit);
    us = G.unit_s(example_unit);
    ue = G.unit_e(example_unit);
    ft = G.all_data(i).ft;

    % recompute the full per-sample lagged signals for this one example
    % trial (mu/cue/bump_ok at every timepoint, not just the pooled
    % bump_ok-only samples kept for the correlation stats above), so the
    % snippet search can slide a window and check each candidate's own
    % offset variance.
    xb = get_xb(ft,size(G.all_data(i).im.d,2));
    keep = ~G.flash{i}(:);
    if G.use_berg4_sign
        fly_sign = 1 - 2*strcmpi(ft.pattern,'berg4');
    else
        fly_sign = 1;
    end
    fs_i = 1/median(diff(ft.xf));
    [mu_l,rho_l,cue_l,r_speed_l,bump_vel_l,xf_l,orig_idx] = lagged_track_signals( ... % position lag: this is a position-tracking display
        ft.xf(:),xb,G.all_data(i).im.mu,G.all_data(i).im.rho,keep,ft.r_speed,ft.cue,round(lag_pos_s(gi)*fs_i),fly_sign);
    bump_ok = abs(r_speed_l)>rot_thresh_track & rho_l>rho_thresh_track & abs(bump_vel_l)<bump_vel_thresh & ...
              ~isnan(mu_l) & ~isnan(cue_l) & ~isnan(rho_l);
    in_unit = orig_idx>=us & orig_idx<=ue;

    [t0,t1] = pick_snippet(xf_l(in_unit),cue_l(in_unit),r_speed_l(in_unit),bump_ok(in_unit), ...
        ft.xf(us),ft.xf(ue),snippet_s,snippet_step_s,min_window_n,xb,keep);

    if ~isnan(forced_snippet_t0(gi))
        t0 = forced_snippet_t0(gi);
        t1 = t0 + snippet_s;
    end

    group_defs(gi).example_fly = chosen_fly;
    group_defs(gi).example_i   = i;
    group_defs(gi).example_t0  = t0;
    group_defs(gi).example_t1  = t1;
end

%% two strong-bump frames ~180 deg apart in each example snippet + their raw summed-z images
% For each column's example snippet: among non-flash frames whose bump
% strength (im.rho) is in the top (100-frame_pair_rho_pct)% of the snippet,
% take the pair of frames whose bump positions are at least
% frame_pair_min_sep_deg apart and, among those, has the largest
% min(rho_A, rho_B) -- i.e. two confidently-localised bumps on opposite
% sides of the PB. (Falls back to the most-separated pair if none reach
% the separation bar.) Then load those two frames' raw registered images
% (summed over z) straight from the trial's own imaging file on Z: so the
% heatmap row can be compared against actual PB images at the same two
% instants. Frames are sorted in time: A is earlier, B later.
frame_pair_min_sep_deg = 100;
frame_pair_rho_pct     = 80;
frame_avg_n            = 10;  % total frames averaged for each raw image, window [f-floor((n-1)/2), f+ceil((n-1)/2)] around A / B (a single z-summed frame is mostly shot noise, esp. the 50x128 px hackathon movies)
frame_img_smooth_sigma = 1;   % px, Gaussian spatial smoothing of the displayed raw image (0 = none)
frame_img_clim_pct     = [5 98]; % percentiles of all pixels (both images) for the raw-image color scale
% manual override of the selected frame TIMES (s, on that trial's own ft.xf
% clock) per column, [tA tB]; NaN = keep the automatic pick for that one.
% LPsP > syt7f's automatic B (t = 363.5 s) fell on a much faster turn than
% its A -- moved ~10 s earlier to a point with a similar rotational speed.
forced_frame_t = {[nan nan], [nan 353.5], [nan nan]};
if ~strcmp(col1_source,'hackathon'), forced_frame_t{1} = [nan nan]; end
for gi = 1:3
    G  = group_defs(gi);
    i  = G.example_i;
    im = G.all_data(i).im;
    ft = G.all_data(i).ft;
    xb = get_xb(ft,size(im.z,2));

    cand  = find(xb>=G.example_t0 & xb<=G.example_t1 & ~G.flash{i}(:));
    rho_c = im.rho(cand); mu_c = im.mu(cand);
    strong = rho_c >= prctile(rho_c,frame_pair_rho_pct);
    cand = cand(strong); rho_c = rho_c(strong); mu_c = mu_c(strong);

    [A,B]  = meshgrid(1:numel(cand));
    sep    = abs(circ_dist(mu_c(A),mu_c(B)));
    score  = min(rho_c(A),rho_c(B));
    upper  = A < B; % each unordered pair once
    ok     = upper & sep >= deg2rad(frame_pair_min_sep_deg);
    if any(ok(:))
        score(~ok) = -inf;
        [~,k] = max(score(:));
    else
        sep(~upper) = -inf;
        [~,k] = max(sep(:));
    end
    fr = sort([cand(A(k)), cand(B(k))]);
    for k = 1:2 % manual time override -> nearest imaging frame
        if ~isnan(forced_frame_t{gi}(k))
            [~,fr(k)] = min(abs(xb - forced_frame_t{gi}(k)));
        end
    end
    group_defs(gi).example_frames = fr;
    rs_at = interp1(ft.xf(:),abs(ft.r_speed(:)),xb(fr)); % fly rotational speed at each selected frame, for the log
    fprintf('%s: frame A %d (t=%.1f s, mu=%+.2f, rho=%.2f, |r|=%.2f rad/s), frame B %d (t=%.1f s, mu=%+.2f, rho=%.2f, |r|=%.2f rad/s), %.0f deg apart\n', ...
        G.title, fr(1), xb(fr(1)), im.mu(fr(1)), im.rho(fr(1)), rs_at(1), fr(2), xb(fr(2)), im.mu(fr(2)), im.rho(fr(2)), rs_at(2), ...
        rad2deg(abs(circ_dist(im.mu(fr(1)),im.mu(fr(2))))));

    im_mask = []; if isfield(im,'mask'), im_mask = im.mask; end % epg7f carries its PB mask in the dataset
    [imgs, mask] = load_summed_frames(G.all_data(i).meta, fr, size(im.z,2), frame_avg_n, im_mask);
    group_defs(gi).example_frame_imgs = imgs;
    group_defs(gi).example_mask       = mask;
end

%% per-fly correlation summary (row 3): closed-loop only, 2 categories:
% position, velocity -- each at its own metric's population-optimal lag
% (unit_wincorr was computed at lag_pos_s, unit_bumpvel/rspeed at
% lag_vel_s). Dark units are still computed above (and kept in group_defs)
% but deliberately not shown in this figure.
cat_labels = {'position','velocity'};
cat_colors = [0.20 0.45 0.85; 0.20 0.45 0.85];

for gi = 1:3
    G = group_defs(gi);
    U = struct('wincorr',{G.unit_wincorr},'bumpvel',{G.unit_bumpvel},'rspeed',{G.unit_rspeed});
    [pv,vv] = per_fly_cl(U,G.unit_fly,G.unit_isdark,min_vel_n);
    pv = pv(~isnan(pv)); vv = vv(~isnan(vv));
    group_defs(gi).cat_x = [ones(1,numel(pv)), 2*ones(1,numel(vv))];
    group_defs(gi).val   = [pv(:)', vv(:)'];
end

% shared y-limits across all 3 row-3 panels: 1 at the top, the minimum
% value actually present anywhere in the data at the bottom (rather than a
% fixed +/-1) so the panels use their full height for the real data range.
row3_ylim = [min([group_defs.val]), 1];

%% ===================== plot =====================
% each column gets its own white-to-color colormap (row 1 only) so its
% heatmap reads as "that indicator's own color" rather than all 3 sharing
% one generic colormap.
group_max_color = [0 .2 1;   % EPG > GCaMP    -> teal
                    .7 0 .7;  % LPsP > syt7f   -> magenta
                    0 .7 0];  % EPG > GRAB(DA2m) -> dark green

% manual override: fix a column's row-1 color scale to a literal [floor,
% ceiling] instead of the per-snippet percentile rule above (empty = use
% the percentile rule). Requested specifically for LPsP > syt7f's
% 20221123-7 trial.
forced_zclim = {[], [-2 8], []};

% row 2 bump-position trace style: 'line' (continuous mu trace, wrap
% breaks blanked) or 'scatter' (one dot per imaging frame, shown only
% where im.rho > row2_rho_min, so frames without a discernable bump drop
% out instead of being drawn). Heading (-cue) is always a line.
row2_mu_style = 'scatter';
row2_rho_min  = 0.1;
row2_dot_size = 6;

% figure width: every panel is made horizontally thinner by narrowing the
% whole page (1500 -> 1200 px, i.e. 0.8x), rather than popping each axes out
% of the tiledlayout individually -- same visual effect, keeps the layout
% intact. Height is unchanged.
fig_width_px = 1200;
figure('color','w','Position',[50 50 fig_width_px 2000]); clf
% 6 visual rows x 3 columns, but laid out on a finer grid so the rows can
% have different heights: row 1 (heatmaps) gets 3 grid units, every other
% row 2, so row 1 is 1.5x as tall as the rest. tile(r,c) gives the
% linear tile index of the top-left grid cell of visual row r, column c;
% each nexttile call spans [row_units(r) 1]. Visual rows:
%   1 heatmaps   2 raw summed-z frames A/B   3 mu vs. heading traces
%   4 per-fly correlations   5 co-imaging example   6 co-imaging summaries
% (code comments below still call the traces "row 2" and the correlation
% dots "row 3" from before the raw-frame row was inserted.)
row_units = [3 2 2 2 2 2];
row_start = cumsum([1, row_units(1:end-1)]); % first grid row of each visual row
tile = @(r,c) (row_start(r)-1)*3 + c;
t = tiledlayout(sum(row_units),3,'TileSpacing','loose','Padding','compact'); % 'loose' (not 'compact') leaves room for rows 1-2's outside-the-axes scale bars; rows 4-5 are the bonus co-imaging panels below
row2_axes = gobjects(1,3);
row3_axes = gobjects(1,3);

for gi = 1:3
    G  = group_defs(gi);
    i  = G.example_i;
    ft = G.all_data(i).ft;
    im = G.all_data(i).im;

    xb = get_xb(ft,size(im.z,2));
    t0 = G.example_t0; t1 = G.example_t1;
    xb_idx = find(xb>=t0 & xb<=t1);
    xf_idx = find(ft.xf>=t0 & ft.xf<=t1);

    mu_xb  = wrap_to_pi(im.mu(xb_idx));
    cue_xf = wrap_to_pi(-ft.cue(xf_idx));

    % color scale: floor and cap both set from THIS SNIPPET's own per-frame
    % extremes (min/max across wedges at each frame), not the literal
    % overall min/max -- ceiling = 95th percentile of the per-frame peak,
    % floor = 5th percentile of the per-frame trough (mirrors the ceiling's
    % rule at the other end; lowered from an earlier 50th-percentile floor,
    % which clipped too much of the low end to see the bump in the LPsP >
    % syt7f panel). Both are deliberately computed from the snippet itself
    % rather than the whole trial: these snippets were chosen FOR strong,
    % well-tracked bump activity, so their peaks/troughs are already biased
    % away from a typical trial-wide percentile -- using the whole trial's
    % own distribution clipped the very activity being highlighted. Flash
    % frames are already excluded by construction (pick_snippet only
    % accepts flash-free windows), so they can't distort either end.
    z_snip = im.z(:,xb_idx);
    z_clim = [prctile(min(z_snip,[],1),5), prctile(max(z_snip,[],1),95)];
    if ~isempty(forced_zclim{gi}) % manual override for this column's row-1 color scale
        z_clim = forced_zclim{gi};
    end

    % row 1: imagesc heatmap only (white-to-this-column's-color), no
    % bump/heading overlay -- see row 2 for those traces
    ax1 = nexttile(t,tile(1,gi),[row_units(1) 1]); hold(ax1,'on')
    imagesc(xb(xb_idx),unwrap(im.alpha),im.z(:,xb_idx),z_clim)
    colormap(ax1, white_to_color(group_max_color(gi,:)))
    cb = colorbar(ax1);
    cb.Ticks = z_clim;
    cb.TickLabels = compose('%.1f',z_clim);
    cb.Label.String = 'z-score';
    y_range1 = [min(unwrap(im.alpha)),max(unwrap(im.alpha))];
    xlim([t0,t1]); ylim(y_range1)
    xticks([])
    ylabel('PB angle (rad)')
    title(sprintf('%s\nexample fly %s (%.0fs closed-loop snippet)',G.title,G.example_fly,t1-t0),'Interpreter','none')
    add_time_scalebar(ax1,t0,t1,y_range1)
    % small black up-arrows just below the bottom edge at the two selected
    % frames (A earlier, B later), with the letter inside the axes right
    % above each arrow; their raw images are in the row directly below
    fr = G.example_frames;
    for k = 1:2
        plot(ax1,xb(fr(k)),y_range1(1)-0.035*diff(y_range1),'^k','MarkerFaceColor','k','MarkerSize',6,'Clipping','off')
        text(ax1,xb(fr(k)),y_range1(1)-0.07*diff(y_range1),char('A'+k-1),'HorizontalAlignment','center','VerticalAlignment','top', ...
            'FontWeight','bold','FontSize',9,'Clipping','off') % letter under the arrow, outside the heatmap
    end

    % new row 2: the raw registered PB image (summed over z, averaged over
    % frame_avg_n frames, lightly smoothed) at frames A and B,
    % side by side in one axes with a transparent gap between them. The
    % full field of view is shown (no mask clipping / outline); the shared
    % color scale is a percentile range of all pixels in both images. Same
    % white-to-color map as the heatmap above. A/B labels + times under
    % each image tie it back to its dashed line in the heatmap.
    ax1b = nexttile(t,tile(2,gi),[row_units(2) 1]); hold(ax1b,'on')
    imgs = G.example_frame_imgs;
    if frame_img_smooth_sigma > 0
        for k = 1:2
            imgs(:,:,k) = imgaussfilt(imgs(:,:,k),frame_img_smooth_sigma);
        end
    end
    % display rotated 180 deg (image and mask together, so the outline
    % stays registered)
    imgs = rot90(imgs,2);
    msk_disp = rot90(G.example_mask,2);
    
    % use only pixels inside the PB mask for the displayed color limits
    mask2 = rot90(G.example_mask, 2);

    [H,W,~] = size(imgs);
    gap  = max(4,round(W*0.04));
    comp = [imgs(:,:,1), nan(H,gap), imgs(:,:,2)];
    valid = ~isnan(comp) & [mask2,false(H,gap),mask2];
    pix = comp(valid);

    clim_f = prctile(pix, frame_img_clim_pct);
    h_img = imagesc(ax1b,comp,clim_f);
    h_img.AlphaData = ~isnan(comp);
    colormap(ax1b, white_to_color(group_max_color(gi,:)))
    axis(ax1b,'image'); axis(ax1b,'ij'); axis(ax1b,'off')
    % PB mask outline (dotted gray) on both images -- the image itself is
    % NOT clipped to the mask
    bnd = bwboundaries(msk_disp);
    for b = 1:numel(bnd)
        plot(ax1b,bnd{b}(:,2),bnd{b}(:,1),':','Color',[0.45 0.45 0.45],'LineWidth',0.8)
        plot(ax1b,bnd{b}(:,2)+W+gap,bnd{b}(:,1),':','Color',[0.45 0.45 0.45],'LineWidth',0.8)
    end
    text(ax1b,W/2,H+1,sprintf('A: t = %.1f s\nframe %d',xb(fr(1)),fr(1)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'FontWeight','bold')
    text(ax1b,1.5*W+gap,H+1,sprintf('B: t = %.1f s\nframe %d',xb(fr(2)),fr(2)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'FontWeight','bold')
    title(ax1b,sprintf('raw PB image at A / B (bumps %.0f%s apart)\nz-summed, mean of %d frames', ...
        rad2deg(abs(circ_dist(im.mu(fr(1)),im.mu(fr(2))))), char(176), frame_avg_n),'FontSize',9)

    % row 2: im.mu (this column's own row-1 color) and -ft.cue (black)
    % traces overlaid, line-only
    ax2 = nexttile(t,tile(3,gi),[row_units(3) 1]); hold(ax2,'on')
    if strcmp(row2_mu_style,'scatter')
        % one dot per imaging frame, only where the bump is discernable
        % (im.rho > row2_rho_min) -- frames with a weak/absent bump are
        % simply left blank rather than drawn as a noisy line
        rho_ok2 = im.rho(xb_idx) > row2_rho_min;
        hm2 = scatter(xb(xb_idx(rho_ok2)),mu_xb(rho_ok2),row2_dot_size,group_max_color(gi,:),'filled');
    else
        hm2 = plot(xb(xb_idx),mu_xb,'-','Color',group_max_color(gi,:),'LineWidth',1.2); hm2.YData(abs(diff(hm2.YData))>pi) = nan;
    end
    hc2 = plot(ft.xf(xf_idx),cue_xf,'-k','LineWidth',1.2); hc2.YData(abs(diff(hc2.YData))>pi) = nan;
    xlim([t0,t1]); ylim([-pi,pi])
    yticks([-pi,pi]); yticklabels({'-\pi','\pi'}) % angle axis in pi units (TeX labels)
    xticks([])
    ylabel('angle (rad)')
    add_time_scalebar(ax2,t0,t1,[-pi,pi])
    if gi == 1
        legend([hm2,hc2],{'bump position (mu)','fly heading (-cue)'},'Location','northeast') % 'northoutside' collided with row 1's new bottom-right-outside scale bar
    end
    row2_axes(gi) = ax2;

    % row 3: per-fly correlation summary, 4 categories
    ax3 = nexttile(t,tile(4,gi),[row_units(4) 1]);
    ylim(row3_ylim) % set before plotting so plot_fly_categories' "n=" labels land at the true bottom
    plot_fly_categories(G.cat_x,G.val,cat_labels,cat_colors,false) % no per-category n= labels: n is in the title
    ylabel('correlation')
    title(sprintf('n=%d flies (closed loop)\nlag: position %+.2f s, velocity %+.2f s',numel(unique(G.unit_fly(~G.unit_isdark))),lag_pos_s(gi),lag_vel_s(gi)))
    row3_axes(gi) = ax3;
end

%% bonus row 4: GRAB(DA2m)/jRGECo1a co-imaging example, 20260113-4
% a single trial where both channels were imaged simultaneously off the
% same PB (same frames, same mask -- see epg_coimaging_script.m /
% lpsp_compartments_claude_script.m section 8), so they're temporally
% aligned by construction (same xb). Plotted as two separate heatmaps
% (not the RGB-subtractive overlay used in that script's own figures)
% colored to match this figure's other columns: jRGECo1a (a calcium
% indicator, like GCaMP) in the EPG>GCaMP column's teal, GRAB(DA2m) in the
% EPG>GRAB(DA2m) column's green.
coimg_file = 'epg_coimage_20260120.mat';
tmp = load(fullfile(data_dir,coimg_file),'all_data');
coimg = tmp.all_data(:);

coimg_i = find(arrayfun(@(s) ischar(s.meta) && contains(s.meta,'20260113-4','IgnoreCase',true), coimg),1);
assert(~isempty(coimg_i), 'trial "20260113-4" not found in %s', coimg_file)

cft = coimg(coimg_i).ft;
cxb = get_xb(cft,size(coimg(coimg_i).grab.z,2));
cxf = cft.xf(:);

coimg_flash_grab = detect_flash_frames(coimg(coimg_i).grab.f, flash_mad_thresh);
coimg_flash_geco = detect_flash_frames(coimg(coimg_i).geco.f, flash_mad_thresh);
coimg_keep = ~(coimg_flash_grab(:) | coimg_flash_geco(:)); % flash on EITHER channel excludes the frame

% discernable-bump mask (both channels must show a bump), same
% rho/rotation gates used throughout this script -- no lag applied here
% (this panel is a visual time-aligned display, not a correlation
% analysis, so behavior/fluorescence alignment doesn't matter).
rho_g_t = interp1(cxb(~coimg_flash_grab(:)),coimg(coimg_i).grab.rho(~coimg_flash_grab(:)),cxf);
rho_r_t = interp1(cxb(~coimg_flash_geco(:)),coimg(coimg_i).geco.rho(~coimg_flash_geco(:)),cxf);
coimg_bump_ok = abs(cft.r_speed(:))>rot_thresh_track & rho_g_t>rho_thresh_track & rho_r_t>rho_thresh_track & ...
                ~isnan(rho_g_t) & ~isnan(rho_r_t);

[coimg_t0,coimg_t1] = pick_snippet(cxf,-cft.cue(:),cft.r_speed(:),coimg_bump_ok, ...
    cxf(1),cxf(end),snippet_s,snippet_step_s,min_window_n,cxb,coimg_keep);
coimg_xb_idx = find(cxb>=coimg_t0 & cxb<=coimg_t1);

coimg_channel = {coimg(coimg_i).geco, coimg(coimg_i).grab};
coimg_color   = [0 .7 .7; 0 .7 0]; % geco -> EPG>GCaMP teal, grab -> EPG>GRAB(DA2m) green
coimg_label   = {'jRGECo1a (co-imaged with GRAB(DA2m))','GRAB(DA2m) (co-imaged with jRGECo1a)'};

for k = 1:2
    ax = nexttile(t,tile(5,k),[row_units(5) 1]); hold(ax,'on')
    ch = coimg_channel{k};
    z_snip = ch.z(:,coimg_xb_idx);
    z_clim = [prctile(min(z_snip,[],1),5), prctile(max(z_snip,[],1),95)];
    imagesc(cxb(coimg_xb_idx),unwrap(ch.alpha),z_snip,z_clim)
    colormap(ax, white_to_color(coimg_color(k,:)))
    cb = colorbar(ax);
    cb.Ticks = z_clim;
    cb.TickLabels = compose('%.1f',z_clim);
    cb.Label.String = 'z-score';
    y_range = [min(unwrap(ch.alpha)),max(unwrap(ch.alpha))];
    xlim([coimg_t0,coimg_t1]); ylim(y_range)
    xticks([])
    ylabel('PB angle (rad)')
    title(sprintf('%s\n20260113-4, fly 2 (%.0fs snippet)',coimg_label{k},coimg_t1-coimg_t0),'Interpreter','none')
    add_time_scalebar(ax,coimg_t0,coimg_t1,y_range)
end

% third bonus panel: both channels' bump position (mu), overlaid, in their
% respective colors -- same time window/axis as the two heatmaps above
ax_mu = nexttile(t,tile(5,3),[row_units(5) 1]); hold(ax_mu,'on')
h_geco = plot(cxb(coimg_xb_idx),wrap_to_pi(coimg_channel{1}.mu(coimg_xb_idx)),'-','Color',coimg_color(1,:),'LineWidth',1.2);
h_geco.YData(abs(diff(h_geco.YData))>pi) = nan;
h_grab = plot(cxb(coimg_xb_idx),wrap_to_pi(coimg_channel{2}.mu(coimg_xb_idx)),'-','Color',coimg_color(2,:),'LineWidth',1.2);
h_grab.YData(abs(diff(h_grab.YData))>pi) = nan;
xlim([coimg_t0,coimg_t1]); ylim([-pi,pi])
xticks([])
ylabel('angle (rad)')
title(sprintf('bump position (mu), both channels\n20260113-4, fly 2 (%.0fs snippet)',coimg_t1-coimg_t0),'Interpreter','none')
add_time_scalebar(ax_mu,coimg_t0,coimg_t1,[-pi,pi])
lgd_mu = legend([h_geco,h_grab],{'jRGECo1a mu','GRAB(DA2m) mu'},'Location','southoutside','Orientation','horizontal'); % below the panel so it never covers the traces; created after the scale bar so its line isn't auto-added as a legend entry. Nudged further down at the very end (after the final layout pass) so it clears the scale bar's "10 s" label.

%% per-fly optimal lag: GRAB onto GECO, by light condition
% Step 1 of 2: find each fly's own best-aligning lag BEFORE computing any
% agreement statistic, so that statistic can be computed on lag-corrected
% (not raw, lag=0) traces below. For each trial, sweep candidate lags,
% circularly shifting GRAB (mu+rho together, so each shifted sample still
% carries its own valid discernible-bump gate) and recomputing mean
% |circ_dist| onto GECO at each one. SIGN CONVENTION (verified against a
% synthetic delayed copy): circshift(grab,+L) makes shifted_grab(t) =
% grab(t-L), so if GRAB truly LAGS GECO by d (grab(t) = geco(t-d)) the
% best alignment is at L = -d. I.e. a NEGATIVE circshift optimum means
% GRAB lags GECO; the plotted axis below is therefore -L ("GRAB delay"),
% so that right = GRAB lags GECO. A fly's optimal lag averages its own trials'
% metric-vs-lag curves (same per-trial-then-average convention as
% the row-3 lag sweep's population-mean-over-flies rule earlier in this script) and takes
% the argmin of that averaged curve. Pools ALL co-imaging trials, not just
% the 20260113-4 example above; discernable-bump restriction (both
% channels must show a bump, fly must be rotating) and flash-frame
% exclusion are the same criteria used throughout this script.
n_coimg = numel(coimg);
coimg_fly_id  = cell(n_coimg,1);
coimg_is_dark = false(n_coimg,1);
coimg_mu_g  = cell(n_coimg,1); % cached per-trial signals, reused below for the lag-corrected agreement/shuffle statistic
coimg_mu_r  = cell(n_coimg,1);
coimg_rho_g = cell(n_coimg,1);
coimg_rho_r = cell(n_coimg,1);
coimg_fs    = nan(n_coimg,1);

lag_grid_coimg_s = -5:0.1:5; % s, candidate lags for the grab-onto-geco lag sweep (widened from +/-2 s so the correlation curves' fall-off on both sides is visible)
coimg_lag_curve  = nan(n_coimg,numel(lag_grid_coimg_s)); % mean |circ_dist| at each candidate lag, per trial (used to PICK each fly's lag)
coimg_corr_curve = nan(n_coimg,numel(lag_grid_coimg_s)); % circ_corrcc(geco mu, shifted grab mu) at each candidate lag, per trial (what row 5 col 1 PLOTS)
% The correlation is computed in sliding window_s-second windows (same
% window_s/window_step_s/min_window_n as row 3's lpsp position correlation)
% and averaged per trial, NOT over the whole trial at once: circ_corrcc
% correlates deviations from each signal's own circular mean, and over a
% whole trial in which heading sweeps the full circle that mean is
% ill-defined and can differ arbitrarily between the two channels, which
% drags a near-perfect agreement (~0.9+ in every 60 s window of
% 20260113-4) down to ~0.55 whole-trial and even moves the peak lag. mean
% |circ_dist| (the lag-picking metric) doesn't have this problem.

for i = 1:n_coimg
    parts = split_path(coimg(i).meta);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
    coimg_fly_id{i}  = strjoin(parts(1:fly_part(1)),'/');
    coimg_is_dark(i) = contains(char(coimg(i).ft.pattern),'background','IgnoreCase',true);

    cflash_g = detect_flash_frames(coimg(i).grab.f, flash_mad_thresh);
    cflash_r = detect_flash_frames(coimg(i).geco.f, flash_mad_thresh);
    cxb_i = get_xb(coimg(i).ft,size(coimg(i).grab.z,2));
    cxf_i = coimg(i).ft.xf(:);
    keep_g = ~cflash_g(:); keep_r = ~cflash_r(:);

    mu_g_im = unwrap(coimg(i).grab.mu(:));
    mu_r_im = unwrap(coimg(i).geco.mu(:));
    mu_g_t  = interp1(cxb_i(keep_g),mu_g_im(keep_g),cxf_i);
    mu_r_t  = interp1(cxb_i(keep_r),mu_r_im(keep_r),cxf_i);
    rho_g_t = interp1(cxb_i(keep_g),coimg(i).grab.rho(keep_g),cxf_i);
    rho_r_t = interp1(cxb_i(keep_r),coimg(i).geco.rho(keep_r),cxf_i);
    coimg_mu_g{i} = mu_g_t; coimg_mu_r{i} = mu_r_t;
    coimg_rho_g{i} = rho_g_t; coimg_rho_r{i} = rho_r_t;
    coimg_fs(i) = 1/median(diff(cxf_i));

    lag_frames_grid_coimg = round(lag_grid_coimg_s*coimg_fs(i));
    coimg_win_starts = cxf_i(1):window_step_s:(cxf_i(end)-window_s);
    for Li = 1:numel(lag_frames_grid_coimg)
        lag = lag_frames_grid_coimg(Li);
        mu_g_shift  = circshift(mu_g_t,lag);
        rho_g_shift = circshift(rho_g_t,lag);
        ok_lag = abs(coimg(i).ft.r_speed(:))>rot_thresh_track & rho_g_shift>rho_thresh_track & rho_r_t>rho_thresh_track & ...
                 ~isnan(mu_g_shift) & ~isnan(mu_r_t);
        if sum(ok_lag) >= min_window_n
            coimg_lag_curve(i,Li) = mean(abs(circ_dist(mu_r_t(ok_lag),mu_g_shift(ok_lag))));
        end
        win_corrs = nan(numel(coimg_win_starts),1);
        for w = 1:numel(coimg_win_starts)
            idx = ok_lag & cxf_i>=coimg_win_starts(w) & cxf_i<coimg_win_starts(w)+window_s;
            if sum(idx) >= min_window_n
                win_corrs(w) = circ_corrcc(mu_r_t(idx),mu_g_shift(idx));
            end
        end
        coimg_corr_curve(i,Li) = mean(win_corrs,'omitnan');
    end
end

coimg_fly_lag = containers.Map('KeyType','char','ValueType','double'); % "flyid|0/1" -> that fly+condition's optimal lag (s)
coimg_fly_corr_curve = [];     % closed loop only: one row per fly, that fly's trial-averaged correlation-vs-lag curve
coimg_fly_corr_isexample = []; % ... and whether that row is the row-4 example fly (20260113 fly 2)
for c = 0 % 0 = closed loop only (1 = dark, excluded from the figure)
    rows = find(coimg_is_dark==logical(c));
    these_flies = unique(coimg_fly_id(rows));
    for f = 1:numel(these_flies)
        sel = rows(strcmp(coimg_fly_id(rows),these_flies{f}));
        mean_curve = mean(coimg_lag_curve(sel,:),1,'omitnan');
        if all(isnan(mean_curve))
            continue
        end
        [~,best_Li] = min(mean_curve);
        coimg_fly_lag(sprintf('%s|%d',these_flies{f},c)) = lag_grid_coimg_s(best_Li);
        coimg_fly_corr_curve(end+1,:)   = mean(coimg_corr_curve(sel,:),1,'omitnan'); %#ok<AGROW>
        coimg_fly_corr_isexample(end+1) = strcmp(these_flies{f},coimg_fly_id{coimg_i}); %#ok<AGROW>
    end
end

% row 5 col 1: per-fly correlation-vs-lag curves (faint gray; the row-4
% example fly in GRAB green), with the across-fly mean on top. Plotted
% against GRAB DELAY = -circshift lag (see sign-convention note above) so
% that right = GRAB lags GECO. Note the lag CORRECTION applied in the next
% section still uses each fly's argmin of mean |circ_dist| (coimg_fly_lag),
% not the argmax of these correlation curves -- the two are expected to
% agree closely but aren't forced to.
grab_delay_s = -lag_grid_coimg_s;
ax_coimg_lag = nexttile(t,tile(6,1),[row_units(6) 1]); hold(ax_coimg_lag,'on')
plot(grab_delay_s,coimg_fly_corr_curve(~coimg_fly_corr_isexample,:)','-','Color',[0 0 0 0.25],'LineWidth',0.8)
h_ex = plot(grab_delay_s,coimg_fly_corr_curve(logical(coimg_fly_corr_isexample),:)','-','Color',coimg_color(2,:),'LineWidth',0.8);
coimg_mean_corr_curve = mean(coimg_fly_corr_curve,1,'omitnan');
h_mn = plot(grab_delay_s,coimg_mean_corr_curve,'-','Color',[0 0 0],'LineWidth',2.5);
plot([0 0],[-1 1],':k')
[~,peak_Li] = max(coimg_mean_corr_curve);
coimg_pop_peak_delay_s = grab_delay_s(peak_Li);
plot(coimg_pop_peak_delay_s*[1 1],[-1 1],'--','Color',[0 0 0],'LineWidth',1.2) % population-mean peak (dashed, vs. the dotted line at 0)
text(coimg_pop_peak_delay_s,0.98,sprintf(' %+.1f s',coimg_pop_peak_delay_s),'Color',[0 0 0],'FontSize',8,'VerticalAlignment','top')
fprintf('GRAB-onto-GECO population correlation peak: GRAB delay = %+.1f s\n', coimg_pop_peak_delay_s);
xlim([min(grab_delay_s),max(grab_delay_s)]); ylim([-0.2 1])
xlabel(sprintf('GRAB delay relative to GECO (s)\n\\leftarrow GRAB leads GECO          GRAB lags GECO \\rightarrow'))
ylabel(sprintf('circular correlation\n(GRAB mu, GECO mu)'))
title(sprintf('GRAB(DA2m) vs. jRGECo1a:\nbump-position correlation vs. lag\n(one trace per fly, closed loop, n=%d)',size(coimg_fly_corr_curve,1)))
% legend lives in the otherwise-empty tile 15 (row 5, col 3) so it doesn't
% cover any of the per-fly traces
legend(ax_coimg_lag,[h_mn,h_ex(1)],{'mean across flies','example fly (row 4)'},'Location','southoutside','Orientation','horizontal') % below the panel (under the x-label)

%% per-fly summary: mean |circ_dist(grab mu, geco mu)|, by light condition
% Step 2 of 2: now apply EACH FLY'S OWN optimal lag (just found above) to
% GRAB before computing the agreement statistic -- so this reflects bump
% agreement once each fly's own sensor-kinetics-driven offset is corrected
% for, rather than the raw (lag=0) alignment. The shuffle-null control
% gets the SAME lag correction applied to GRAB first, then an ADDITIONAL
% large random circular shift (>= min_shuffle_lag_s seconds, with the
% usual near-0/near-trial-length exclusion band) applied to GECO -- so
% real and shuffled data are put through identical preprocessing and only
% differ in that extra decorrelating shift.
coimg_offset = cell(n_coimg,1);
min_shuffle_lag_s = 5;  % s, minimum |lag| for a circular-shuffle draw
n_shuffles        = 200; % shuffle draws per trial
coimg_shuffle_offset = cell(n_coimg,n_shuffles);
rng(1) % reproducible shuffles run-to-run

for i = 1:n_coimg
    mu_g_t = coimg_mu_g{i}; mu_r_t = coimg_mu_r{i};
    rho_g_t = coimg_rho_g{i}; rho_r_t = coimg_rho_r{i};
    fs_i = coimg_fs(i);
    r_speed_i = coimg(i).ft.r_speed(:);

    key = sprintf('%s|%d',coimg_fly_id{i},coimg_is_dark(i));
    if isKey(coimg_fly_lag,key)
        lag0 = round(coimg_fly_lag(key)*fs_i);
    else
        lag0 = 0; % this fly+condition had no usable lag curve -- fall back to uncorrected
    end
    mu_g_corr  = circshift(mu_g_t,lag0);
    rho_g_corr = circshift(rho_g_t,lag0);

    ok = abs(r_speed_i)>rot_thresh_track & rho_g_corr>rho_thresh_track & rho_r_t>rho_thresh_track & ...
         ~isnan(mu_g_corr) & ~isnan(mu_r_t);
    coimg_offset{i} = circ_dist(mu_g_corr(ok),mu_r_t(ok));

    n_i  = numel(mu_r_t);
    min_shift = round(min_shuffle_lag_s*fs_i);
    eligible_lags = (min_shift:(n_i-min_shift))';
    for sidx = 1:n_shuffles
        lag = eligible_lags(randi(numel(eligible_lags)));
        mu_r_shift  = circshift(mu_r_t,lag);
        rho_r_shift = circshift(rho_r_t,lag);
        ok_s = abs(r_speed_i)>rot_thresh_track & rho_g_corr>rho_thresh_track & rho_r_shift>rho_thresh_track & ...
               ~isnan(mu_g_corr) & ~isnan(mu_r_shift);
        coimg_shuffle_offset{i,sidx} = circ_dist(mu_g_corr(ok_s),mu_r_shift(ok_s));
    end
end

coimg_cat_labels = {'CL (real)','CL (shifted)'};
coimg_cat_colors = [0.20 0.45 0.85; 0.65 0.78 0.92];
coimg_cat_x = []; coimg_val = [];
for c = 0 % 0 = closed loop only (1 = dark, excluded from the figure)
    rows = find(coimg_is_dark==logical(c));
    these_flies = unique(coimg_fly_id(rows));
    for f = 1:numel(these_flies)
        sel = rows(strcmp(coimg_fly_id(rows),these_flies{f}));

        pooled = cat(1,coimg_offset{sel});
        if isempty(pooled)
            continue
        end
        coimg_cat_x(end+1) = 1 + 2*c; %#ok<AGROW>
        coimg_val(end+1)   = mean(abs(pooled)); %#ok<AGROW>

        % this fly's shuffle-null estimate: for each of the n_shuffles
        % independent draws, pool this fly's own trials' shuffled samples
        % and take the mean (mirroring the real-data statistic exactly),
        % then average that draw-level mean over all draws -- a
        % Monte-Carlo estimate of this fly's expected offset under the
        % null, directly comparable to its one real-data point.
        draw_means = nan(n_shuffles,1);
        for sidx = 1:n_shuffles
            pooled_shuf = cat(1,coimg_shuffle_offset{sel,sidx});
            if ~isempty(pooled_shuf)
                draw_means(sidx) = mean(abs(pooled_shuf));
            end
        end
        coimg_cat_x(end+1) = 2 + 2*c; %#ok<AGROW>
        coimg_val(end+1)   = mean(draw_means,'omitnan'); %#ok<AGROW>
    end
end

ax_coimg_summary = nexttile(t,tile(6,2),[row_units(6) 1]);
plot_fly_categories(coimg_cat_x,coimg_val,coimg_cat_labels,coimg_cat_colors)
ylim([0,pi])
yticks([0 pi/4 pi/2 3*pi/4 pi])
yticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
ylabel(sprintf('mean |circ\\_dist(GRAB mu, GECO mu)|\n(rad)'))
title(sprintf('GRAB(DA2m) vs. jRGECo1a bump agreement,\nper fly (lag-corrected) vs. %ds+\ncircularly-shifted null (n=%d draws/trial)',min_shuffle_lag_s,n_shuffles))

sgtitle('Figure 1: bump tracks fly heading across indicators/genotypes')

% hide rows 2-3's x-axis LINE only (keep row 3's category tick labels and
% its plotted y=0 dotted line) -- done as the very last step, right before
% export: XAxis.Axle is an internal ruler primitive that gets rebuilt on
% layout changes, so setting it mid-loop got silently undone by the
% subsequent nexttile/sgtitle calls re-laying out the tiledlayout.
drawnow
% scooch row 2 down toward row 3 to tighten the gap between them. The
% tiledlayout owns its children's Position, so pop each row-2 axes out of
% the layout first (reparenting to the figure keeps its on-screen
% Position), then shift it down. Row-2 axes carry no colorbar, and their
% scale bars are in-axes data-coordinate objects, so both move with them;
% the row-2 col-1 legend has a named Location and re-tracks its axes.
row2_shift = 0.025; % normalized figure units, downward
for gi = 1:3
    p = row2_axes(gi).Position;
    row2_axes(gi).Parent = gcf;
    row2_axes(gi).Position = [p(1), p(2)-row2_shift, p(3), p(4)];
end
drawnow
for gi = 1:3
    row2_axes(gi).XAxis.Axle.Visible = 'off';
    row3_axes(gi).XAxis.Axle.Visible = 'off';
end
% nudge the row-4 mu-overlay legend down so it clears the panel's "10 s"
% scale-bar label (which sits just below the axes box). Setting Position
% detaches the legend from the layout's auto-placement, which is fine
% here since this is the final layout.
lgd_mu.Position(2) = lgd_mu.Position(2) - 0.012; % normalized figure units
ax_mu.XAxis.Axle.Visible = 'off';
ax_coimg_summary.XAxis.Axle.Visible = 'off'; % (ax_coimg_lag keeps its x-axis: it's a real lag axis now, not categories)

%% export
if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(gcf, fullfile(export_dir,'Figure_1_claude.png'), 'Resolution', 300)
savefig(gcf, fullfile(export_dir,'Figure_1_claude.fig')) % so the PDF/PNG can be re-exported (openfig + exportgraphics) without the ~7 min recompute

all_figs_dir = fullfile(repo_root,'ugly_figures','all_figs');
if ~isfolder(all_figs_dir); mkdir(all_figs_dir); end
try
    exportgraphics(gcf, fullfile(all_figs_dir,'Fig1_V1.pdf'), 'ContentType', 'auto') % 'auto' rasterizes the dense heatmaps but keeps text/lines vector -- 'vector' would bloat the file turning each heatmap pixel into its own path
catch err
    % most likely the PDF is open in a viewer and locked -- don't fail the
    % whole (long) run over it; the PNG above has already been written
    warning('Figure1:pdfExport','could not write Fig1_V1.pdf (%s) -- is it open in a viewer?', err.message);
end

%% ===================== functions =====================

function combined = combine_datasets(varargin)
    % vertically concatenate struct arrays that don't all share the same
    % top-level fields, padding any fields missing from a given dataset
    % with [] so concatenation works. (ported from lpsp_compartments_claude_script.m)
    all_fields = {};
    for k = 1:numel(varargin)
        all_fields = union(all_fields,fieldnames(varargin{k}),'stable');
    end
    combined = [];
    for k = 1:numel(varargin)
        s = varargin{k}(:);
        missing = setdiff(all_fields,fieldnames(s));
        for m = 1:numel(missing)
            s(1).(missing{m}) = [];
        end
        s = orderfields(s,all_fields);
        if isempty(combined)
            combined = s;
        else
            combined = [combined; s]; %#ok<AGROW>
        end
    end
end

function ind = trial_indicator(meta_path)
    % ported from lpsp_compartments_claude_script.m
    if contains(meta_path,'syt7f','IgnoreCase',true)
        ind = 'syt7f';
    elseif contains(meta_path,'syt8m','IgnoreCase',true)
        ind = 'syt8m';
    elseif contains(meta_path,'dlight','IgnoreCase',true)
        ind = 'dLight';
    elseif contains(meta_path,'grab','IgnoreCase',true) || contains(meta_path,'DA2m','IgnoreCase',true)
        ind = 'GRAB(DA2m)';
    else
        ind = '';
    end
end

function fid = trial_fly_id(dataset_name, meta_path)
    % ported from lpsp_compartments_claude_script.m
    parts = split_path(meta_path);
    parts(cellfun(@isempty,parts)) = [];
    switch dataset_name
        case {'lpsp_cl_redo','epg_dlight'}
            fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
            if isempty(fly_part)
                fid = meta_path;
            else
                fid = strjoin(parts(1:fly_part(1)),filesep);
            end
        otherwise % lpsp_cl
            tname    = parts{end};
            date_str = regexp(tname,'^\d{8}','match','once');
            suffix   = regexp(tname,'_(\d+)$','tokens','once');
            if isempty(suffix)
                fid = date_str;
            else
                fid = [date_str,'_',suffix{1}];
            end
    end
end

function is_flash = detect_flash_frames(f, mad_k)
    % ported from lpsp_compartments_claude_script.m
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;
    med_m  = median(frame_mean);
    mad_m  = mad(frame_mean,1);
    med_cv = median(frame_cv);
    is_bright  = frame_mean > med_m + mad_k*mad_m;
    is_uniform = frame_cv <= med_cv;
    is_flash   = is_bright & is_uniform;
end

function xb = get_xb(ft, n_im)
    % real imaging timestamps when available and the right length, else a
    % linear placeholder spanning the trial (same fallback used throughout
    % lpsp_compartments_claude_script.m)
    if isfield(ft,'xb') && numel(ft.xb) == n_im
        xb = ft.xb(:);
    else
        xb = linspace(ft.xf(1),ft.xf(end),n_im)';
    end
end

function g = local_measured_gain(ft)
    % ported verbatim from hackathon_claude.m: estimate closed-loop VR gain
    % from behavior when no applied-gain field exists for this pattern.
    r  = ft.r_speed(1:end-2);
    c  = -gradient(ft.cue(3:end)) * 60;
    c(abs(c) > 50) = nan;
    xf = ft.xf(1:end-2);

    t = 1:round(max(ft.xf));
    g = nan(size(t));
    for j = 1:length(t)
        idx = xf > j-1 & xf < j+5 & abs(r) > .5 & abs(c) > .5;
        if sum(idx) > 0
            g(j) = r(idx) \ c(idx);
        end
    end
    g = g(abs(g) < 5);
end

function [run_starts,run_ends,run_vals] = find_runs(x)
    % ported verbatim from hackathon_claude.m
    x = x(:);
    n = numel(x);
    is_new = [true; x(2:end)~=x(1:end-1) | isnan(x(2:end)) | isnan(x(1:end-1))];
    run_starts = find(is_new);
    run_ends   = [run_starts(2:end)-1; n];
    run_vals   = x(run_starts);
end

function C = prep_hackathon_column(data_dir, flash_mad_thresh, f0_pct, smooth_frames_f)
    % column-1 source 'hackathon': load .data/hackathon_20250729.mat, recompute
    % the bump from smoothed im.f exactly as hackathon_claude.m, detect flash
    % frames, derive fly ids, and chunk every trial into constant-gain/dark
    % blocks keeping gain=0.8 (closed loop and dark). Returns the common
    % column struct C (data, flash, fly, unit_*, use_berg4_sign, title).
    tmp = load(fullfile(data_dir,'hackathon_20250729.mat'),'all_data');
    hack_data = tmp.all_data(:);
    fprintf('loaded %d hackathon trials\n', numel(hack_data));


    % -- recompute bump (mu,rho,z,d) from smoothed im.f, exactly as hackathon_claude.m
    for i = 1:numel(hack_data)
        f   = smoothdata(hack_data(i).im.f,2,'movmean',smooth_frames_f);
        f0  = prctile(f,f0_pct,2);
        dff = (f - f0) ./ f0;
        z   = zscore(dff,[],2);

        alpha    = hack_data(i).im.alpha;
        [x,y]    = pol2cart(alpha,z');
        [mu,rho] = cart2pol(mean(x,2),mean(y,2));
        mu = mod(mu,2*pi); mu(mu>pi) = mu(mu>pi)-2*pi;

        hack_data(i).im.mu  = mu;
        hack_data(i).im.rho = rho;
        hack_data(i).im.z   = z;
        hack_data(i).im.d   = dff;
    end

    % -- flash-frame detection + fly ID
    n_hack = numel(hack_data);
    hack_flash    = cell(n_hack,1);
    hack_fly      = cell(n_hack,1);

    for i = 1:n_hack
        hack_flash{i} = detect_flash_frames(hack_data(i).im.f, flash_mad_thresh);

        tok = regexp(hack_data(i).meta,'\d+_[a-z]{2}_\d+','match','once');
        if isempty(tok)
            parts = split_path(hack_data(i).meta);
            tok = parts{end-2};
        end
        hack_fly{i} = tok;
    end

    % -- segment every hackathon trial into constant-gain/dark blocks, gain=0.8 only
    % ported from hackathon_claude.m's compute_gain_blocks (block-splitting part
    % only -- this figure doesn't need its bump-mobility/gain_fly/gain_cue metrics)
    hack_flies = unique(hack_fly);

    hack_unit_trial  = [];
    hack_unit_s      = [];
    hack_unit_e      = [];
    hack_unit_fly    = {};
    hack_unit_isdark = [];

    for f = 1:numel(hack_flies)
        trial_idx = find(strcmp(hack_fly,hack_flies{f}));
        last_gain = nan;
        for i = trial_idx(:)'
            ft = hack_data(i).ft;
            xf = ft.xf(:);
            pattern = ft.pattern;

            if strcmpi(pattern,'berg4')
                cb    = ft.cue_brightness(:);
                state = round(ft.gain(:),2);
                state(cb==0) = -inf;
                [run_starts,run_ends,run_vals] = find_runs(state);
                run_isdark = isinf(run_vals);
            elseif contains(pattern,'background')
                run_starts = 1; run_ends = numel(xf); run_vals = nan; run_isdark = true;
            else
                run_starts = 1; run_ends = numel(xf); run_vals = median(local_measured_gain(ft)); run_isdark = false;
            end

            for k = 1:numel(run_starts)
                if run_isdark(k)
                    label_gain = last_gain;
                    is_dark_block = true;
                else
                    if ~isnan(run_vals(k))
                        last_gain = run_vals(k);
                    end
                    label_gain = run_vals(k);
                    is_dark_block = false;
                end

                if ~isnan(label_gain) && round(label_gain,1)==0.8
                    hack_unit_trial(end+1,1)  = i;   %#ok<AGROW>
                    hack_unit_s(end+1,1)      = run_starts(k); %#ok<AGROW>
                    hack_unit_e(end+1,1)      = run_ends(k);   %#ok<AGROW>
                    hack_unit_fly{end+1,1}    = hack_fly{i};   %#ok<AGROW>
                    hack_unit_isdark(end+1,1) = is_dark_block; %#ok<AGROW>
                end
            end
        end
    end
    n_hack_units = numel(hack_unit_trial);
    fprintf('%d gain=0.8 hackathon blocks (%d closed loop, %d dark)\n', n_hack_units, sum(~hack_unit_isdark), sum(hack_unit_isdark));

    C.data = hack_data; C.flash = hack_flash; C.fly = hack_fly;
    C.unit_trial = hack_unit_trial; C.unit_s = hack_unit_s; C.unit_e = hack_unit_e;
    C.unit_fly = hack_unit_fly; C.unit_isdark = logical(hack_unit_isdark);
    C.use_berg4_sign = true; % berg4 trials need r_speed sign-flipped (see header)
    C.title = 'EPG > GCaMP';
end

function C = prep_epg7f_column(data_dir, flash_mad_thresh)
    % column-1 source 'epg7f': load .data/epg_7f_20261007.mat (built by
    % data scripts/epg_7f_gain07_batch.m via claude/pb_gain_trial.m) and adapt
    % it to the lab all_data layout this script expects:
    %   - behavior at fictrac rate: that dataset stores cue/heading/r_speed
    %     RESAMPLED onto the volume clock in ft.cue etc. and keeps the 60 Hz
    %     originals as ft.*_raw -- the latter become ft.cue / ft.r_speed here
    %     (the script interpolates imaging onto ft.xf itself)
    %   - im.d (only used for its frame count) = im.f; im.mask kept for the
    %     raw-frame row; dropout volumes (im.badVol) are folded into the
    %     per-trial "flash" exclusion list
    %   - one unit per trial (whole trial), dark from ft.dark, fly id
    %     '<date>_flyN' (or '<date>' for dates without fly subfolders)
    % Sign conventions were checked on all closed-loop trials: bump position
    % correlates with -cue (as for lpsp) and bump velocity with +r_speed, so
    % no rotation sign flip is needed.
    tmp = load(fullfile(data_dir,'epg_7f_20261007.mat'),'all_data');
    src = tmp.all_data(:);
    n = numel(src);
    data  = struct('ft',cell(n,1),'im',cell(n,1),'meta',cell(n,1));
    flash = cell(n,1);
    fly   = cell(n,1);
    for i = 1:n
        s = src(i);
        ft = struct();
        ft.xf = s.ft.xf(:); ft.xb = s.ft.xb(:);
        ft.cue = s.ft.cue_raw(:); ft.r_speed = s.ft.r_speed_raw(:); ft.f_speed = s.ft.f_speed_raw(:);
        ft.pattern = char(s.ft.pattern); ft.dark = logical(s.ft.dark); ft.gain_empirical = s.ft.gain_empirical;
        im = struct('f',s.im.f,'z',s.im.z,'d',s.im.f,'mu',s.im.mu(:),'rho',s.im.rho(:),'alpha',s.im.alpha(:),'mask',logical(s.im.mask));
        data(i).ft = ft; data(i).im = im; data(i).meta = s.meta;
        bad = false(1,size(s.im.f,2));
        if isfield(s.im,'badVol') && ~isempty(s.im.badVol), bad = logical(s.im.badVol(:))'; end
        flash{i} = detect_flash_frames(s.im.f, flash_mad_thresh) | bad;
        parts = split_path(s.meta.flyDir); parts(cellfun(@isempty,parts)) = [];
        if startsWith(parts{end},'fly')
            fly{i} = sprintf('%s_%s', parts{end-1}, strrep(parts{end},' ',''));
        else
            fly{i} = parts{end};
        end
    end
    C.data = data; C.flash = flash; C.fly = fly;
    C.unit_trial  = (1:n)';
    C.unit_s      = ones(n,1);
    C.unit_e      = arrayfun(@(d) numel(d.ft.xf), data);
    C.unit_fly    = fly;
    C.unit_isdark = arrayfun(@(d) d.ft.dark, data);
    C.use_berg4_sign = false;
    C.title = 'EPG > GCaMP7f';
    fprintf('loaded %d epg_7f trials (%d flies; %d closed loop, %d dark)\n', n, numel(unique(fly)), sum(~C.unit_isdark), sum(C.unit_isdark));
end

function U = unit_metrics(all_data, flash, unit_trial, unit_s, unit_e, unit_isdark, lag_s, use_berg4_sign, P)
    % per-unit bump-tracking signals at ONE lag (lag_s, seconds; converted
    % to frames per trial from that trial's own ft.xf rate; positive =
    % fluorescence shifted earlier, see lagged_track_signals). A "unit" is a
    % (trial, start index, end index) chunk on the ft.xf timebase: one
    % gain=0.8 block for hackathon, one whole trial for lpsp. Returns:
    %   wincorr - circ_corrcc(mu,-cue) in sliding P.window_s windows (step
    %             P.window_step_s) over discernable-bump samples; a unit
    %             shorter than one window is treated as a single window
    %   bumpvel, rspeed - pooled discernable-bump samples (velocity Pearson r)
    %   offset  - circ_dist(mu,-cue) over discernable-bump samples
    %             (closed-loop units only)
    %   nbumpok - discernable-bump sample count
    % lagged_track_signals is run once per trial and cached across that
    % trial's units.
    n_units = numel(unit_trial);
    U.wincorr = cell(n_units,1);
    U.bumpvel = cell(n_units,1);
    U.rspeed  = cell(n_units,1);
    U.offset  = cell(n_units,1);
    U.nbumpok = zeros(n_units,1);
    cache = containers.Map('KeyType','double','ValueType','any');
    for u = 1:n_units
        i = unit_trial(u);
        if ~isKey(cache,i)
            ft = all_data(i).ft;
            xf = ft.xf(:);
            xb = get_xb(ft,size(all_data(i).im.d,2));
            keep = ~flash{i}(:);
            if use_berg4_sign
                fly_sign = 1 - 2*strcmpi(ft.pattern,'berg4'); % -1 for berg4, +1 otherwise
            else
                fly_sign = 1;
            end
            fs = 1/median(diff(xf));
            [mu_l,rho_l,cue_l,r_speed_l,bump_vel_l,xf_l,orig_idx] = lagged_track_signals( ...
                xf,xb,all_data(i).im.mu,all_data(i).im.rho,keep,ft.r_speed,ft.cue,round(lag_s*fs),fly_sign);
            bump_ok = abs(r_speed_l)>P.rot_thresh_track & rho_l>P.rho_thresh_track & abs(bump_vel_l)<P.bump_vel_thresh & ...
                      ~isnan(mu_l) & ~isnan(cue_l) & ~isnan(rho_l);
            cache(i) = struct('mu_l',mu_l,'cue_l',cue_l,'bump_vel_l',bump_vel_l,'r_speed_l',r_speed_l, ...
                              'xf_l',xf_l,'orig_idx',orig_idx,'bump_ok',bump_ok,'xf',xf);
        end
        c = cache(i);
        mask = c.bump_ok & c.orig_idx>=unit_s(u) & c.orig_idx<=unit_e(u);
        U.nbumpok(u) = sum(mask);
        U.bumpvel{u} = c.bump_vel_l(mask);
        U.rspeed{u}  = c.r_speed_l(mask);
        if ~unit_isdark(u)
            U.offset{u} = circ_dist(c.mu_l(mask),c.cue_l(mask));
        end

        t_lo = c.xf(unit_s(u)); t_hi = c.xf(unit_e(u));
        if t_hi - t_lo <= P.window_s
            win_starts = t_lo; win_len = (t_hi - t_lo) + 1e-9; % whole unit as one window (inclusive of t_hi)
        else
            win_starts = t_lo:P.window_step_s:(t_hi-P.window_s); win_len = P.window_s;
        end
        wc = [];
        for w = win_starts
            idx = mask & c.xf_l>=w & c.xf_l<w+win_len;
            if sum(idx) < P.min_window_n
                continue
            end
            cc = circ_corrcc(c.mu_l(idx),c.cue_l(idx));
            if ~isnan(cc)
                wc(end+1) = cc; %#ok<AGROW>
            end
        end
        U.wincorr{u} = wc;
    end
end

function [imgs, mask] = load_summed_frames(meta, frames, n_frames_expected, n_avg, im_mask)
    % raw PB images (summed over z, then averaged over n_avg frames
    % [f-floor((n_avg-1)/2), f+ceil((n_avg-1)/2)] around each requested
    % centre frame f -- for even n_avg the window is one frame longer after
    % f than before it) plus the
    % trial's PB mask, read from the trial's own registered imaging file on
    % the Z: drive. Two dataset layouts:
    %   hackathon_gain_change: <meta>\imagingData.mat, var imgData
    %       (y x frames; already z-summed), mask.mat one folder up
    %   lpsp_cl: <meta>\imagingData_reg_ch1_trial_001.mat (v7.3), var
    %       regProduct (y x z frames), partially loaded through matfile so
    %       the 1.4-1.7 GB file isn't read whole; mask.mat in the same folder
    % Frame k of these files is column k of im.z (process_im keeps frames
    % 1:1); n_frames_expected = size(im.z,2) is checked against the file,
    % and if they disagree the requested frames are rescaled proportionally
    % (with a warning) rather than erroring out of a long run.
    if isstruct(meta)
        % epg7f (gain_change rig): pb_gain_trial caches the registered,
        % z-summed movie as <trialDir>\claude_pipeline\imgData_sum_reg_trial_00k.mat
        % (var imgData_sum_reg, y x frames, single); the mask comes with the
        % dataset (im.mask), passed in as im_mask
        mf = matfile(fullfile(meta.trialDir,'claude_pipeline',sprintf('imgData_sum_reg_trial_%03d.mat',meta.tifIdx)));
        sz = size(mf,'imgData_sum_reg');
        n_file = sz(3);
        frames = remap_frames(frames,n_frames_expected,n_file,meta.trialDir);
        imgs = zeros(sz(1),sz(2),numel(frames));
        for k = 1:numel(frames)
            win = max(1,frames(k)-floor((n_avg-1)/2)):min(n_file,frames(k)+ceil((n_avg-1)/2));
            imgs(:,:,k) = mean(double(mf.imgData_sum_reg(:,:,win)),3);
        end
        M.mask = im_mask;
    elseif contains(meta,'hackathon_gain_change','IgnoreCase',true)
        S = load(fullfile(meta,'imagingData.mat'),'imgData');
        n_file = size(S.imgData,3);
        frames = remap_frames(frames,n_frames_expected,n_file,meta);
        imgs = zeros(size(S.imgData,1),size(S.imgData,2),numel(frames));
        for k = 1:numel(frames)
            win = max(1,frames(k)-floor((n_avg-1)/2)):min(n_file,frames(k)+ceil((n_avg-1)/2));
            imgs(:,:,k) = mean(double(S.imgData(:,:,win)),3);
        end
        M = load(fullfile(fileparts(meta),'mask.mat'),'mask');
    else
        mf = matfile(fullfile(meta,'imagingData_reg_ch1_trial_001.mat'));
        sz = size(mf,'regProduct');
        n_file = sz(4);
        frames = remap_frames(frames,n_frames_expected,n_file,meta);
        imgs = zeros(sz(1),sz(2),numel(frames));
        for k = 1:numel(frames)
            win = max(1,frames(k)-floor((n_avg-1)/2)):min(n_file,frames(k)+ceil((n_avg-1)/2));
            imgs(:,:,k) = mean(sum(double(mf.regProduct(:,:,:,win)),3),4);
        end
        M = load(fullfile(meta,'mask.mat'),'mask');
    end
    mask = logical(M.mask);
end

function frames = remap_frames(frames, n_expected, n_file, meta)
    if n_file ~= n_expected
        warning('Figure1:frameCountMismatch', ...
            '%s: imaging file has %d frames but im.z has %d columns -- rescaling requested frame indices proportionally', meta, n_file, n_expected);
        frames = max(1,min(n_file,round(frames * n_file / n_expected)));
    end
end

function [pos_vals, vel_vals] = per_fly_cl(U, unit_fly, unit_isdark, min_vel_n)
    % closed-loop per-fly summaries from unit_metrics output: mean windowed
    % position correlation (NaN if the fly has no usable windows) and
    % pooled velocity Pearson r (NaN if fewer than min_vel_n samples), one
    % entry per closed-loop fly, in unique(fly) order.
    rows = find(~unit_isdark);
    flies = unique(unit_fly(rows));
    pos_vals = nan(numel(flies),1);
    vel_vals = nan(numel(flies),1);
    for f = 1:numel(flies)
        sel = rows(strcmp(unit_fly(rows),flies{f}));
        wc = cat(2,U.wincorr{sel});
        if ~isempty(wc)
            pos_vals(f) = mean(wc,'omitnan');
        end
        bv = cat(1,U.bumpvel{sel}); rs = cat(1,U.rspeed{sel});
        if numel(bv) >= min_vel_n
            vel_vals(f) = corr(bv,rs);
        end
    end
end

function [mu_l,rho_l,cue_l,r_speed_l,bump_vel_l,xf_l,orig_idx] = lagged_track_signals( ...
        xf, xb, mu_im, rho_im, keep_frames, r_speed_raw, cue_raw, lag, r_sign)
    % shared per-trial signal prep for both the hackathon and lpsp columns:
    % interpolate bump position/strength onto the behavior (xf) timebase
    % (flash frames excluded first), apply this group's lag (fluorescence
    % lags behavior), and return everything needed for the discernable-bump
    % mask downstream. orig_idx maps each output sample back onto the
    % ORIGINAL xf frame index, so a caller can intersect with block
    % boundaries defined on that same original index (needed for
    % hackathon's gain=0.8 blocks; harmless/unused for lpsp's whole-trial
    % units). mu is unwrapped before interpolating/lagging (bump position is
    % cumulative); cue is used raw/wrapped, matching
    % lpsp_compartments_claude_script.m section 7's own convention -- both
    % circ_corrcc and circ_dist are wrap-safe, so cue never needs unwrapping.
    mu_im_u = unwrap(mu_im(:));
    mu_t  = interp1(xb(keep_frames),mu_im_u(keep_frames),xf);
    rho_t = interp1(xb(keep_frames),rho_im(keep_frames),xf);
    cue_t = -cue_raw(:);
    r_t   = r_sign * r_speed_raw(:);

    n = numel(xf);
    if lag == 0
        idx_lead = (1:n)'; idx_lag = (1:n)';
    elseif lag > 0
        idx_lead = (1:n-lag)'; idx_lag = (lag+1:n)';
    else
        idx_lead = (-lag+1:n)'; idx_lag = (1:n+lag)';
    end

    xf_l      = xf(idx_lead);
    cue_l     = cue_t(idx_lead);
    r_speed_l = r_t(idx_lead);
    mu_l      = mu_t(idx_lag);
    rho_l     = rho_t(idx_lag);
    orig_idx  = idx_lead;

    dt = mean(diff(xf));
    bump_vel_l = gradient(mu_l)/dt;
end

function [t0,t1] = pick_snippet(xf_l, cue_l, r_speed_l, bump_ok, t_lo, t_hi, snippet_s, step_s, min_n, xb, im_keep)
    % best snippet_s-second sub-window of [t_lo,t_hi]: among candidates
    % whose peak |r_speed| is at or below the MEDIAN peak speed across all
    % candidate windows in this unit (a data-driven "relatively low speed"
    % bar, rather than a fixed rad/s cutoff that wouldn't transfer across
    % this figure's very different rigs/indicators), pick the one with the
    % widest coverage of head directions (measured as 2*pi minus the
    % largest angular gap between its sampled headings -- close to 2*pi
    % means the fly faced nearly every direction at some point). This
    % favors a slow, wide-ranging meander over a fast, narrow-range spin:
    % a big turn's rapid mu/cue swings visually mask the bump (motivating
    % this speed cap), and "many different head directions" is what the
    % old plain-circ_var-of-heading objective was already reaching for,
    % just measured more directly here as actual angular coverage.
    %
    % peak |r_speed| is evaluated over EVERY sample in the window
    % (regardless of bump_ok) since a fast excursion looks messy in the
    % example plot whether or not it clears the discernable-bump gate;
    % coverage is evaluated only over bump_ok samples, consistent with
    % every other "good tracking" criterion used elsewhere in this script.
    %
    % min_n is a floor, not the real bar: a window where the fly barely
    % rotates has very few bump_ok samples, and any dispersion statistic
    % estimated from a handful of samples is unreliable (verified directly
    % with an earlier offset-variance objective -- a "perfect" window
    % turned out to have bump_ok samples over only ~14% of the snippet).
    % min_frac additionally requires bump_ok samples over a real fraction
    % of the window, so the winning snippet reflects sustained tracking.
    %
    % Candidate windows containing any flagged flash frame (xb/im_keep,
    % im_keep=true meaning "not a flash frame") are excluded outright.
    % Falls back to the whole range if it's already <= snippet_s (e.g. a
    % short hackathon gain=0.8 block) regardless of flash content -- no
    % shorter option to fall back to.
    min_frac = 0.2;
    if t_hi - t_lo <= snippet_s
        t0 = t_lo; t1 = t_hi;
        return
    end

    starts = t_lo:step_s:(t_hi-snippet_s);
    cand_start = []; cand_cov = []; cand_peak = [];
    for pass = 1:2 % pass 1: flash-free candidates only; pass 2 (only if pass 1 found none): flash allowed
        for w = starts
            if pass == 1 && any(~im_keep(xb>=w & xb<w+snippet_s))
                continue % a flash frame falls inside this candidate window
            end
            in_win = xf_l>=w & xf_l<w+snippet_s;
            idx = in_win & bump_ok;
            if sum(idx) < max(min_n, min_frac*sum(in_win))
                continue
            end
            cand_start(end+1) = w;                        %#ok<AGROW>
            cand_cov(end+1)   = circular_coverage(cue_l(idx)); %#ok<AGROW>
            cand_peak(end+1)  = max(abs(r_speed_l(in_win)));  %#ok<AGROW>
        end
        if ~isempty(cand_start)
            break
        end
    end

    if isempty(cand_start)
        % nothing met the bump_ok sample threshold at all -- just take the
        % first window rather than returning nothing
        t0 = starts(1); t1 = t0 + snippet_s;
        return
    end

    speed_cap = median(cand_peak);
    scores = cand_cov;
    scores(cand_peak > speed_cap) = -inf; % restrict the coverage argmax to the slower half of candidates
    [~,rel] = max(scores);
    t0 = cand_start(rel); t1 = t0 + snippet_s;
end

function cov = circular_coverage(theta)
    % how much of the full circle (radians) is spanned by these samples,
    % measured as 2*pi minus the single largest gap between consecutive
    % angles once sorted around the circle -- close to 2*pi means the
    % samples are spread all the way around; close to 0 means they're all
    % clustered in one small arc.
    if isempty(theta)
        cov = 0;
        return
    end
    th = sort(mod(theta(:),2*pi));
    gaps = diff([th; th(1)+2*pi]);
    cov = 2*pi - max(gaps);
end

function chosen_fly = pick_best_fly(fly_list, rot_amount, off_var)
    % "turned a lot" (rotation amount, i.e. discernable-bump sample count,
    % at or above the median across candidate flies) AND, among those, the
    % lowest circular variance of the mu-(-cue) offset.
    valid = ~isnan(rot_amount) & ~isnan(off_var);
    thresh = median(rot_amount(valid));
    candidates = valid & rot_amount >= thresh;
    if ~any(candidates)
        candidates = valid;
    end
    idx_list = find(candidates);
    [~,rel] = min(off_var(candidates));
    chosen_fly = fly_list{idx_list(rel)};
end

function cmap = white_to_color(max_color)
    % linear colormap from white (low) to max_color (high), for a single
    % axes -- imagesc's caxis/clim already maps the data's own [floor,
    % ceiling] onto this full range, so "bottoms out at white" just means
    % white sits at the first row.
    n = 256;
    cmap = [linspace(1,max_color(1),n)', linspace(1,max_color(2),n)', linspace(1,max_color(3),n)'];
end

function add_time_scalebar(ax, t0, t1, y_range)
    % 10 s horizontal scale bar below the bottom-right corner, OUTSIDE the
    % axes' own plotted data area (used instead of a numeric time axis on
    % rows 1-2 -- xticks/xlabel removed there) -- Clipping is turned off on
    % the bar/label themselves so they render in that outside margin
    % rather than being cut off at the axes box.
    bar_len = 10; % s
    x1 = t1 - 0.02*(t1-t0);
    x0 = x1 - bar_len;
    y0 = y_range(1) - 0.07*diff(y_range);
    line(ax,[x0,x1],[y0,y0],'Color','k','LineWidth',2.5,'Clipping','off')
    text(ax,x0+bar_len/2,y0,'10 s','Color','k','VerticalAlignment','top', ...
        'HorizontalAlignment','center','FontSize',9,'Clipping','off')
end

function parts = split_path(p)
    % meta paths are stored as literal Windows paths (Z:\pablo\...)
    % regardless of the host OS this script runs on, so splitting on
    % filesep (which is '/' on a Mac) would never actually split them --
    % split on either separator explicitly instead.
    parts = strsplit(p,{'\\','/'});
end

function y = wrap_to_pi(x)
    % ported from lpsp_compartments_claude_script.m
    y = mod(x+pi,2*pi) - pi;
end

function plot_fly_categories(cat_x, values, cat_labels, colors, show_n)
    % adapted from lpsp_compartments_claude_script.m's groupplot(): one
    % jittered dot per fly per category + mean+/-sem, with an n= label
    % (show_n=false suppresses the n= labels, e.g. when n is in the title).
    if nargin < 5
        show_n = true;
    end
    hold on
    n_cat = numel(cat_labels);
    for c = 1:n_cat
        y = values(cat_x==c);
        y = y(~isnan(y));
        if isempty(y)
            continue
        end
        jitter = (rand(size(y))-.5)*.3;
        scatter(c+jitter,y,20,colors(c,:),'filled','MarkerFaceAlpha',.4)
        errorbar(c,mean(y),std(y)/sqrt(numel(y)),'o','Color',colors(c,:)*.6, ...
            'MarkerFaceColor',colors(c,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    xticks(1:n_cat); xticklabels(cat_labels); xtickangle(20)
    xlim([0.5,n_cat+0.5])
    if show_n
        y_lims = ylim;
        for c = 1:n_cat
            n = sum(cat_x==c & ~isnan(values));
            text(c,y_lims(1),sprintf('n=%d',n),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
        end
    end
    plot(xlim,[0,0],':k')
end
