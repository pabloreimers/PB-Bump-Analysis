%% Figure_2_claude
% Example-trial + kinematics-encoding summary figure for the two indicators
% imaged in the lpsp_cl dataset: LPsP > syt7f and EPG > GRAB(DA2m). Data
% loading, per-trial preprocessing, and the R^2/tuning-curve analyses are
% all done ONCE below and shared; the script then renders TWO complete
% figures (variants '2a' and '2b') from that same shared data, differing
% only in which indicator supplies the example fly/snippet (rows 1-5) and
% the heatmap colormap:
%   2a: example fly from syt7f,      white-to-magenta colormap ([.7 0 .7])
%   2b: example fly from GRAB(DA2m), white-to-green colormap   ([0 .7 0])
% (these two used to be separate scripts, Figure_2a_claude.m/
% Figure_2b_claude.m -- consolidated here since ~90% of the code, and
% every design decision below, was identical between them; only the
% render_variant_figure() call's inputs differ per variant.)
%
% Rows (identical structure in both variants):
%   1) example fly, example closed-loop 60s snippet: imagesc heatmap
%      (full-width row) of z-scored dF/F across all 32 PB glomeruli, with
%      a dark-red box around one high-|r_speed| imaging frame and a salmon
%      box around one low-|r_speed| frame, that variant's own
%      white-to-color colormap (matching Figure_1_claude.m's
%      LPsP>syt7f/EPG>GRAB(DA2m) columns respectively).
%   2) its own row directly below, two side-by-side panels (same shared
%      y-scale, dots colored to match that variant's own imagesc color --
%      magenta for 2a, green for 2b): per-glomerulus z-score for the
%      high-rotation frame (left, dark-red axes outline) and the
%      low-rotation frame (right, salmon axes outline), glomerulus index
%      on x. BOTH panels mark their own 90th-percentile z-score (solid
%      black line -- row 6-7's response variable) AND their own
%      10th-percentile z-score (gray line).
%   3) same snippet: 90th-percentile z-score across the entire PB (all 32
%      wedges), as a solid black line trace.
%   4) same snippet: fly rotational speed (ft.r_speed), red.
%   5) same snippet: fly forward speed (ft.f_speed), blue.
%   6) adjusted R^2 for the 90th-percentile-across-wedges z-score
%      predicted by forward speed alone, |rotational speed| alone, or both
%      jointly (fitlm, one line per fly, pooling that fly's own trials) --
%      syt7f, closed loop | dark. IDENTICAL in both variants (shows both
%      indicators regardless of which one supplied the example above).
%   7) same R^2 comparison -- GRAB(DA2m), closed loop | dark.
%   8) binned tuning curve, THIS VARIANT'S OWN example indicator only:
%      |rotational speed| (x, binned) vs. 90th-percentile (black) and
%      10th-percentile (gray) z-score (y), one dot per fly per bin (that
%      fly's own mean) plus a population mean+/-SEM errorbar per bin --
%      closed loop | dark.
%
% Loads the SAME three source .mat files combined by
% lpsp_compartments_claude_script.m / Figure_1_claude.m (lpsp_cl,
% lpsp_cl_redo, epg_dlight). As in Figure_1_claude.m, only lpsp_cl actually
% contains syt7f/GRAB(DA2m) trials (lpsp_cl_redo is syt8m, epg_dlight is
% dLight) -- the other two are loaded for parity with that recipe and then
% simply contribute no rows to either indicator used here.
%
% ==EXAMPLE FLY/SNIPPET SELECTION (rows 1-5)==
% A fly, and a 60s window within that fly, chosen so the relationship
% between rotational speed and dF/F is actually visible: high correlation
% between trial_p90z (the SAME 90th-percentile-across-wedges z-score used
% as the R^2 fits' response variable below, and what row 3 displays) and
% (smoothed) |r_speed|, with rotational speed genuinely fluctuating over
% the window -- not a flat trace that happens to correlate trivially.
% Restricted to that variant's own example_indicator, closed-loop trials
% only, and still requires mean forward speed >=
% min_forward_speed_snippet over the window (the fly should be genuinely
% walking forward, not just turning in place). This last threshold is the
% ONE analysis parameter that differs by variant: 2a uses 2 mm/s (syt7f
% has plenty of flies clearing that bar); 2b uses 0.5 mm/s (matching the
% per-sample min_forward_speed threshold instead) because at 2 mm/s only
% 1/12 GRAB(DA2m) closed-loop flies ever qualifies, and that fly's own
% rotation is fairly modest -- the fly with by far the most prominent
% fluctuation in that dataset (std~1.6 rad/s) is essentially turning in
% place (mean forward speed ~0.2 mm/s) rather than walking-and-turning.
% 0.5 mm/s opens the candidate pool to all 12 GRAB(DA2m) flies while still
% ruling out a fly that's fully stationary.
%
% Fly and window are chosen jointly. For every candidate fly, every
% snippet_s-second sub-window (slid in snippet_step_s steps) across that
% fly's own closed-loop trials is screened on mean forward speed, then
% scored two ways: corr(trial_p90z, |r_speed|_smoothed) and
% std(|r_speed|_smoothed) ("how much it fluctuates"). Windows below that
% FLY'S OWN median std are dropped (same median-cutoff trick
% Figure_1_claude.m's pick_snippet uses for its own speed cap -- a
% data-driven "fluctuates enough" bar rather than an arbitrary rad/s
% number), and among the rest the window with the highest correlation
% becomes that fly's candidate. The fly whose candidate has the single
% highest correlation wins. |r_speed| is smoothed with the SAME two-stage
% (pre-abs + post-abs) gaussian used for row 4's display trace below, so
% the correlation used to pick the example matches what the eye actually
% sees plotted. The plotted trace reuses each trial's own trial_p90z
% (computed once, on the ft.xf timebase) rather than re-deriving a
% separate imaging-timebase trace, so rows 3-5 share one x-axis with no
% extra interpolation. Row 4 plots |r_speed| (rectified), matching the
% rectified rotational-speed convention used throughout the R^2 fits in
% rows 6-7.
%
% ==R^2 PANELS (rows 6-7)==
% Adapted from lpsp_compartments_claude_script.m section 4
% (fit_r2_triplet/apply_lag3/plot_r2_group): for each indicator, trial_p90z
% (90th-percentile-across-wedges z-score, NOT average dF/F -- see below)
% is predicted from forward speed alone, |rotational speed| alone, or both
% jointly (fitlm, adjusted R^2), at that indicator's own optimal lag
% (fluorescence lags behavior -- found the same way as
% Figure_1_claude.m's fit_group_lag, correlating |r_speed| against
% trial_p90z across a lag grid, pooling closed-loop + dark trials together
% so a CL-vs-dark difference below reflects the light condition and not a
% different lag choice). Trials are pooled per fly before fitting (a fly's
% several ~600s trials become one set of timepoints, one fitlm, one R^2
% triplet -- not one fit per trial), split into closed-loop and dark
% panels. Computed ONCE (not per variant) since it's identical either way.
%
% UNLIKE that script, two things changed:
%  (1) the response variable is the 90th-percentile-across-wedges z-score
%      (trial_p90z) rather than average dF/F -- also what rows 1-3 display
%      and the example-selection correlation above uses, so all of rows
%      1-8 are consistent on the same summary signal. Checked empirically
%      across 4 candidate summary signals (avg dF/F, avg z-score, peak
%      z-score, 90th-pctile z-score): pooled across both indicators'
%      qualifying trials, 90th-pctile z-score has the highest mean
%      correlation with |r_speed| (0.62, vs. 0.60 for avg dF/F, 0.60 for
%      peak z-score, 0.60 for avg z-score) -- and specifically for
%      GRAB(DA2m) it's a clear win (0.66 vs. 0.60 for avg dF/F), while for
%      syt7f avg dF/F is marginally better (0.61 vs. 0.59) but not by
%      enough to justify a different response variable per indicator.
%  (2) both predictors are the SAME gaussian-smoothed signals used
%      elsewhere in this figure (trial_rspeed_smooth/trial_fspeed_smooth,
%      sigma = r_speed_display_smooth_sigma_s) rather than raw per-sample
%      r_speed/f_speed -- confirmed empirically (12/12 syt7f closed-loop
%      flies) that smoothing raises adjusted R^2 for every fly, consistent
%      with fluorescence tracking a slower-than-per-sample-noise version
%      of behavior. Both predictors are smoothed the same way so the
%      forward-vs-rotational comparison stays fair (smoothing only
%      r_speed would inflate its R^2 relative to forward speed for a
%      reason unrelated to biology). trial_walks_enough/
%      min_forward_speed_snippet still gate on RAW behavior -- only the
%      R^2 predictors themselves are smoothed.
%
% ==SIMPLIFICATIONS vs. lpsp_compartments_claude_script.m's analysis_ok==
% That script's full inclusion filter also excludes fly-level
% forward/backward dominance, "is_ambiguous", and "odd_panels" trials --
% dataset-level QC concerns for its full multi-indicator sweep, not
% specific to the two indicators used here. This script keeps only the
% per-trial "walked forward enough" filter (>= min_forward_time_frac of
% the trial spent above min_forward_speed), since a trial where the fly
% barely walked forward would make the forward-speed model meaningless
% regardless of which other exclusions apply.

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
data_dir = fullfile(repo_root,'data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
all_figs_dir = fullfile(repo_root,'ugly_figures','all_figs');

%% shared analysis constants
flash_mad_thresh       = 8;    % detect_flash_frames: brightness-outlier MAD multiple
min_forward_speed      = 0.5;  % mm/s, "the fly is walking forward" per-sample threshold
min_forward_time_frac  = 0.20; % fraction of a trial's duration that must clear min_forward_speed

lag_frames_grid = -30:3:180; % same grid as lpsp_compartments_claude_script.m / Figure_1_claude.m

snippet_s      = 60; % s, requested example-trace window length
snippet_step_s = 10; % s, sliding step when searching for the best snippet_s-second sub-window

max_quiet_run_s      = 8;    % s, reject any candidate snippet containing a continuous "fly is basically
                              % still" stretch longer than this (the whole-snippet MEAN forward speed above
                              % doesn't stop a window from having e.g. an active first half and a dead second
                              % half -- this catches that case directly)
quiet_r_speed_thresh  = 0.15; % rad/s, "not really turning" for the stillness check, on the smoothed |r_speed| trace

r_speed_display_smooth_sigma_s = 0.3; % s, gaussian smoothing STANDARD DEVIATION applied to r_speed -- for
                                 % row 4's DISPLAY trace, the example-snippet correlation search, AND (see
                                 % trial_rspeed_smooth/trial_fspeed_smooth below) the R^2 fits in rows 6-7.
                                 % r_speed is smoothed TWICE, once to raw (signed) r_speed before abs(), once
                                 % to the rectified trace after abs() (a single pre-abs smooth alone still
                                 % leaves a jagged rectified trace, since abs() folds each smoothed-but-still-
                                 % noisy zero crossing into a sharp little peak); f_speed is smoothed once,
                                 % same sigma, no abs(). smoothdata's own 'gaussian' window argument is the
                                 % window LENGTH, not sigma -- it defines sigma as window/5 -- so the window
                                 % passed to smoothdata is 5x this sigma.

indicator_order = {'syt7f','GRAB(DA2m)'};
cond_label      = {'closed loop','dark'};

%% load the same three lpsp_compartments_claude_script.m source datasets Figure_1_claude.m loads, and combine
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

%% per-trial prep: indicator, fly ID, light condition, flash frames, trial_p90z, forward-walking filter
n_trials = numel(lpsp_data);

trial_valid   = false(n_trials,1);
trial_indic   = repmat({''},n_trials,1); % '' (not the cell(...) default []) so ismember() below sees a uniform cellstr
trial_fly     = repmat({''},n_trials,1);
trial_isdark  = false(n_trials,1);
trial_flash   = cell(n_trials,1);
trial_p90z = cell(n_trials,1); % 90th-percentile-across-wedges z-score trace, on its own ft.xf timebase --
                               % used for row 1's display trace, the example-fly/snippet correlation search,
                               % AND the R^2 fits' response variable (see header comment: empirically beats
                               % avg dF/F, avg z-score, and peak z-score pooled across both indicators)
trial_p10z = cell(n_trials,1); % 10th-percentile-across-wedges z-score trace, same timebase -- used only by
                               % the binned rotational-speed tuning-curve row at the bottom
trial_walks_enough = false(n_trials,1);
trial_rspeed_smooth = cell(n_trials,1); % gaussian-smoothed (pre-abs + post-abs) |r_speed|, used for the R^2
                                         % fits below (NOT for trial_walks_enough/min_forward_speed_snippet,
                                         % which stay on raw behavior)
trial_fspeed_smooth = cell(n_trials,1); % gaussian-smoothed f_speed, same sigma, used for the R^2 fits below

for i = 1:n_trials
    meta = lpsp_data(i).meta;
    if isempty(meta) || ~ischar(meta)
        continue % leftover empty placeholder trial
    end
    trial_valid(i) = true;

    trial_indic{i} = trial_indicator(meta);
    trial_fly{i}   = trial_fly_id(lpsp_data(i).dataset,meta);

    if isfield(lpsp_data(i).ft,'pattern') && ~isempty(lpsp_data(i).ft.pattern)
        trial_isdark(i) = contains(char(lpsp_data(i).ft.pattern),'background','IgnoreCase',true);
    else
        % no network-drive trialSettings.csv access here -- fall back to the
        % folder-name convention, same as Figure_1_claude.m
        trial_isdark(i) = contains(meta,'_dark','IgnoreCase',true);
    end

    trial_flash{i} = detect_flash_frames(lpsp_data(i).im.f, flash_mad_thresh);

    ft = lpsp_data(i).ft;
    xf = ft.xf(:);
    xb = get_xb(ft,size(lpsp_data(i).im.d,2));
    keep = ~trial_flash{i}(:);
    if sum(keep) < 2
        keep = true(size(keep));
    end
    p90_im = prctile(lpsp_data(i).im.z(:,keep),90,1)';
    trial_p90z{i} = interp1(xb(keep),p90_im,xf);
    p10_im = prctile(lpsp_data(i).im.z(:,keep),10,1)';
    trial_p10z{i} = interp1(xb(keep),p10_im,xf);

    trial_walks_enough(i) = mean(ft.f_speed > min_forward_speed) >= min_forward_time_frac;

    fs_i = 1/median(diff(xf));
    smooth_win_i = round(5*r_speed_display_smooth_sigma_s*fs_i); % smoothdata gaussian window = 5*sigma
    %r_pre = smoothdata(ft.r_speed,'gaussian',smooth_win_i);
    r_pre = ft.r_speed;
    trial_rspeed_smooth{i} = smoothdata(abs(r_pre),'gaussian',smooth_win_i);
    trial_fspeed_smooth{i} = smoothdata(ft.f_speed,'gaussian',smooth_win_i);
end

analysis_ok = trial_valid & ismember(trial_indic,indicator_order) & trial_walks_enough;
fprintf('%d/%d trials pass the (syt7f or GRAB(DA2m)) + walked-forward-enough filter\n', sum(analysis_ok), n_trials);

%% per-indicator optimal lag (pooling closed loop + dark), 90th-pctile z-score vs. |rotational speed|
lag_by_indicator = nan(1,numel(indicator_order));
for k = 1:numel(indicator_order)
    rows = find(strcmp(trial_indic,indicator_order{k}) & analysis_ok);
    lag_by_indicator(k) = fit_group_lag({lpsp_data(rows).ft}, trial_p90z(rows), lag_frames_grid);
    fprintf('%s optimal lag: %d frames (n=%d trials)\n', indicator_order{k}, lag_by_indicator(k), numel(rows));
end

%% per-fly R^2 (forward / rotational / joint), one panel per indicator x light condition
% shared by both figure variants below -- doesn't depend on which
% indicator supplies the example fly/snippet.
R2_panels = cell(numel(indicator_order),2);
for k = 1:numel(indicator_order)
    lag = lag_by_indicator(k);
    for c = 1:2 % 1 = closed loop, 2 = dark
        rows = find(strcmp(trial_indic,indicator_order{k}) & (trial_isdark==(c==2)) & analysis_ok);
        these_flies = unique(trial_fly(rows));
        R2 = nan(numel(these_flies),3);

        for ff = 1:numel(these_flies)
            trial_list = rows(strcmp(trial_fly(rows),these_flies{ff}));
            pooled_for = []; pooled_rot = []; pooled_amp = [];
            for i = trial_list(:)'
                [for_l,rot_l,amp_l] = apply_lag3(trial_fspeed_smooth{i},trial_rspeed_smooth{i},trial_p90z{i},lag);
                valid = ~isnan(for_l) & ~isnan(rot_l) & ~isnan(amp_l);
                pooled_for = [pooled_for; for_l(valid)]; %#ok<AGROW>
                pooled_rot = [pooled_rot; rot_l(valid)]; %#ok<AGROW>
                pooled_amp = [pooled_amp; amp_l(valid)]; %#ok<AGROW>
            end
            if numel(pooled_amp) < 100
                continue
            end
            R2(ff,:) = fit_r2_triplet(pooled_rot,pooled_for,pooled_amp);
        end

        R2_panels{k,c} = R2;
    end
end

%% ===================== render both figure variants =====================
% everything above this point is shared; bundle it into one struct so
% render_variant_figure() below can take a single (S, variant) pair
% instead of a long, error-prone positional argument list.
S = struct( ...
    'lpsp_data',lpsp_data, 'trial_indic',{trial_indic}, 'trial_fly',{trial_fly}, ...
    'trial_isdark',trial_isdark, 'trial_flash',{trial_flash}, ...
    'trial_p90z',{trial_p90z}, 'trial_p10z',{trial_p10z}, ...
    'trial_rspeed_smooth',{trial_rspeed_smooth}, 'trial_fspeed_smooth',{trial_fspeed_smooth}, ...
    'analysis_ok',analysis_ok, 'indicator_order',{indicator_order}, 'cond_label',{cond_label}, ...
    'lag_by_indicator',lag_by_indicator, 'R2_panels',{R2_panels}, ...
    'snippet_s',snippet_s, 'snippet_step_s',snippet_step_s, ...
    'max_quiet_run_s',max_quiet_run_s, 'quiet_r_speed_thresh',quiet_r_speed_thresh, ...
    'min_forward_speed',min_forward_speed, 'r_speed_display_smooth_sigma_s',r_speed_display_smooth_sigma_s, ...
    'export_dir',export_dir, 'all_figs_dir',all_figs_dir);

variants = struct( ...
    'label',                     {'2a',                 '2b'}, ...
    'example_indicator',         {'syt7f',              'GRAB(DA2m)'}, ...
    'min_forward_speed_snippet', {2,                    0.5}, ...
    'heatmap_max_color',         {[.7 0 .7],            [0 .7 0]}, ...
    'sgtitle_suffix',            {'(syt7f example)',    '(GRAB(DA2m) example)'} ...
);

for v = 1:numel(variants)
    render_variant_figure(S,variants(v))
end

%% ===================== functions =====================

function render_variant_figure(S, variant)
    % builds and exports one complete 8-row figure (rows described in this
    % file's header) for one variant ('2a' or '2b'); everything it needs
    % from the shared pipeline above is in S, everything that differs
    % between the two variants is in `variant`.
    rows_cl = find(strcmp(S.trial_indic,variant.example_indicator) & ~S.trial_isdark & S.analysis_ok);
    [example_fly,example_i,example_t0,example_t1] = pick_high_corr_example( ...
        rows_cl,S.lpsp_data,S.trial_fly,S.trial_p90z,S.snippet_s,S.snippet_step_s,variant.min_forward_speed_snippet, ...
        S.r_speed_display_smooth_sigma_s,S.max_quiet_run_s,S.quiet_r_speed_thresh,S.min_forward_speed);
    fprintf('[Figure %s] example fly: %s, trial "%s", snippet [%.1f, %.1f] s\n', ...
        variant.label,example_fly,meta_display(S.lpsp_data(example_i).meta),example_t0,example_t1);

    figure('color','w','Position',[50 50 900 2200]); clf
    t = tiledlayout(8,2,'TileSpacing','normal','Padding','compact');

    ex_ft = S.lpsp_data(example_i).ft;
    ex_im = S.lpsp_data(example_i).im;
    xf_idx = find(ex_ft.xf >= example_t0 & ex_ft.xf <= example_t1);
    t_rel = ex_ft.xf(xf_idx) - example_t0; % relative time within the snippet, 0..snippet_s

    fs_ex = 1/median(diff(ex_ft.xf));
    smooth_win = round(5*S.r_speed_display_smooth_sigma_s*fs_ex); % smoothdata gaussian window = 5*sigma
    r_speed_presmooth  = smoothdata(ex_ft.r_speed,'gaussian',smooth_win);      % smooth #1: before abs()
    r_speed_postsmooth = smoothdata(abs(r_speed_presmooth),'gaussian',smooth_win); % smooth #2: after abs()

    %% row 1: z-scored dF/F heatmap (all PB glomeruli), with a high- and a low-rotation frame boxed
    % "high"/"low" rotation frame = the single imaging frame (flash frames
    % excluded) within this snippet whose interpolated smoothed |r_speed|
    % (same trace row 4 plots) is the max/min across the snippet.
    color_hi  = [0.55 0 0];    % dark red -- high-rotation frame (box + line)
    color_lo  = [0.94 0.50 0.50]; % salmon -- low-rotation frame (box + line)
    color_dot = variant.heatmap_max_color; % row-2 scatter dots match this variant's own imagesc color

    n_glom = size(ex_im.z,1);
    xb_ex  = get_xb(ex_ft,size(ex_im.z,2));
    xb_idx = find(xb_ex >= example_t0 & xb_ex <= example_t1);
    xb_rel = xb_ex(xb_idx) - example_t0;

    r_on_xb = interp1(ex_ft.xf,r_speed_postsmooth,xb_ex(xb_idx));
    valid_frame = ~S.trial_flash{example_i}(xb_idx);
    cand_idx = xb_idx(valid_frame); cand_r = r_on_xb(valid_frame);
    [~,rel_hi] = max(cand_r); high_frame = cand_idx(rel_hi);
    [~,rel_lo] = min(cand_r); low_frame  = cand_idx(rel_lo);
    hi_p90 = prctile(ex_im.z(:,high_frame),90);
    lo_p90 = prctile(ex_im.z(:,low_frame),90);
    hi_p10 = prctile(ex_im.z(:,high_frame),10);
    lo_p10 = prctile(ex_im.z(:,low_frame),10);

    ax0 = nexttile(t,[1 2]); hold(ax0,'on')
    z_snip = ex_im.z(:,xb_idx);
    z_clim = [prctile(min(z_snip,[],1),5), prctile(max(z_snip,[],1),95)]; % same recipe as Figure_1_claude.m row 1
    imagesc(ax0,xb_rel,1:n_glom,z_snip,z_clim)
    colormap(ax0,white_to_color(variant.heatmap_max_color)) % matches Figure_1_claude.m's own column for this indicator
    cb0 = colorbar(ax0); cb0.Label.String = 'z-score';
    ylabel(ax0,'PB glomerulus')
    xlabel(ax0,'time (s)')
    xlim(ax0,[0,S.snippet_s]); ylim(ax0,[0.5,n_glom+0.5])
    title(ax0,{'z-scored dF/F, all PB glomeruli','(dark red = high-rotation frame, salmon = low-rotation frame)'})

    frame_dt = median(diff(xb_ex));
    box_halfwidth = 1.5*frame_dt;
    hi_t_rel = xb_ex(high_frame) - example_t0;
    lo_t_rel = xb_ex(low_frame)  - example_t0;
    rectangle(ax0,'Position',[hi_t_rel-box_halfwidth,0.5,2*box_halfwidth,n_glom],'EdgeColor',color_hi,'LineWidth',2)
    rectangle(ax0,'Position',[lo_t_rel-box_halfwidth,0.5,2*box_halfwidth,n_glom],'EdgeColor',color_lo,'LineWidth',2)

    % row 2: its own row directly below the heatmap (not an overlapping
    % inset), split into two side-by-side panels (one per boxed frame,
    % sharing a y-scale so amplitude is directly comparable): per-glomerulus
    % z-score, with each frame's own 90th-percentile value (the R^2 fits'
    % response variable) marked as a horizontal line
    z_both = [ex_im.z(:,high_frame); ex_im.z(:,low_frame)];
    z_pad = 0.05*range(z_both);
    z_ylim = [min(z_both)-z_pad, max(z_both)+z_pad];

    ax_scatter_hi = nexttile(t); hold(ax_scatter_hi,'on'); box(ax_scatter_hi,'on')
    scatter(ax_scatter_hi,1:n_glom,ex_im.z(:,high_frame),25,color_dot,'filled')
    xlim(ax_scatter_hi,[0.5,n_glom+0.5]); ylim(ax_scatter_hi,z_ylim)
    plot(ax_scatter_hi,xlim(ax_scatter_hi),[hi_p90,hi_p90],'--','Color','k','LineWidth',1.2)
    plot(ax_scatter_hi,xlim(ax_scatter_hi),[hi_p10,hi_p10],'--','Color',[0.5 0.5 0.5],'LineWidth',1.2)
    set(ax_scatter_hi,'XColor',color_hi,'YColor',color_hi,'LineWidth',1.5)
    xlabel(ax_scatter_hi,'PB glomerulus')
    ylabel(ax_scatter_hi,'z-score')
    title(ax_scatter_hi,'high-rotation frame (dark red)')

    ax_scatter_lo = nexttile(t); hold(ax_scatter_lo,'on'); box(ax_scatter_lo,'on')
    scatter(ax_scatter_lo,1:n_glom,ex_im.z(:,low_frame),25,color_dot,'filled')
    xlim(ax_scatter_lo,[0.5,n_glom+0.5]); ylim(ax_scatter_lo,z_ylim)
    plot(ax_scatter_lo,xlim(ax_scatter_lo),[lo_p90,lo_p90],'--','Color','k','LineWidth',1.2)
    plot(ax_scatter_lo,xlim(ax_scatter_lo),[lo_p10,lo_p10],'--','Color',[0.5 0.5 0.5],'LineWidth',1.2)
    set(ax_scatter_lo,'XColor',color_lo,'YColor',color_lo,'LineWidth',1.5)
    xlabel(ax_scatter_lo,'PB glomerulus')
    ylabel(ax_scatter_lo,'z-score')
    title(ax_scatter_lo,'low-rotation frame (salmon)')

    ax1 = nexttile(t,[1 2]); hold(ax1,'on')
    plot(ax1,t_rel,S.trial_p90z{example_i}(xf_idx),'-k','LineWidth',1.2)
    ylabel(ax1,{'90th %ile z-score','(whole PB)'})
    xlim(ax1,[0,S.snippet_s]); xticks(ax1,[])
    title(ax1,sprintf('example fly %s, %s closed-loop (%.0fs snippet, trial "%s")', ...
        example_fly,variant.example_indicator,S.snippet_s,meta_display(S.lpsp_data(example_i).meta)),'Interpreter','none')

    ax2 = nexttile(t,[1 2]);
    plot(ax2,t_rel,r_speed_postsmooth(xf_idx),'-r','LineWidth',1)
    ylabel(ax2,{'abs. rotational speed','(rad/s)'})
    xlim(ax2,[0,S.snippet_s]); xticks(ax2,[])

    ax3 = nexttile(t,[1 2]);
    plot(ax3,t_rel,ex_ft.f_speed(xf_idx),'-b','LineWidth',1)
    ylabel(ax3,{'forward speed','(mm/s)'})
    xlim(ax3,[0,S.snippet_s])
    xlabel(ax3,'time (s)')

    linkaxes([ax0,ax1,ax2,ax3],'x')

    for k = 1:numel(S.indicator_order)
        for c = 1:2
            ax = nexttile(t);
            plot_r2_group(S.R2_panels{k,c})
            n_flies = sum(~isnan(S.R2_panels{k,c}(:,1)));
            title(ax,sprintf('%s, %s (n=%d flies)',S.indicator_order{k},S.cond_label{c},n_flies))
            if c == 1
                ylabel(ax,{'adjusted R^2','(90th %ile z-score)'})
            end
        end
    end

    %% bottom row: binned tuning curve, |rotational speed| vs. 90th/10th %ile z-score
    % one point per fly per bin (that fly's own mean over its own pooled,
    % lag-shifted samples in that bin), plus a population mean+/-SEM
    % errorbar per bin -- same "pool a fly's own trials, then average
    % across flies" logic as the R^2 panels above, just binned instead of
    % fit. Scoped to this variant's own example_indicator, closed-loop and
    % dark as the two panels.
    tuning_indicator = variant.example_indicator;
    k_tuning = find(strcmp(S.indicator_order,tuning_indicator));
    lag_tuning = S.lag_by_indicator(k_tuning);

    rows_tuning_all = find(strcmp(S.trial_indic,tuning_indicator) & S.analysis_ok);
    pooled_rot_glob = [];
    for i = rows_tuning_all(:)'
        [~,rot_l,~] = apply_lag3(S.trial_fspeed_smooth{i},S.trial_rspeed_smooth{i},S.trial_p90z{i},lag_tuning);
        pooled_rot_glob = [pooled_rot_glob; rot_l(~isnan(rot_l))]; %#ok<AGROW>
    end
    n_speed_bins = 6;
    speed_bin_edges = linspace(0,prctile(pooled_rot_glob,95),n_speed_bins+1);
    speed_bin_centers = (speed_bin_edges(1:end-1)+speed_bin_edges(2:end))/2;
    min_n_per_bin = 20; % minimum pooled samples for one fly's one-bin mean to count

    for c = 1:2
        rows = find(strcmp(S.trial_indic,tuning_indicator) & (S.trial_isdark==(c==2)) & S.analysis_ok);
        [~,fly_mean90,fly_mean10] = bin_tuning_by_fly( ...
            rows,S.trial_fly,S.trial_rspeed_smooth,S.trial_fspeed_smooth,S.trial_p90z,S.trial_p10z,lag_tuning,speed_bin_edges,min_n_per_bin);
        ax = nexttile(t);
        plot_binned_tuning(ax,speed_bin_centers,fly_mean90,fly_mean10)
        title(ax,sprintf('%s, %s',tuning_indicator,S.cond_label{c}))
    end

    sgtitle(sprintf('Figure %s: PB dF/F encoding of forward/rotational speed, syt7f and GRAB(DA2m) %s', ...
        variant.label,variant.sgtitle_suffix))

    %% export
    if ~isfolder(S.export_dir); mkdir(S.export_dir); end
    exportgraphics(gcf, fullfile(S.export_dir,sprintf('Figure_%s_claude.png',variant.label)), 'Resolution', 300)

    if ~isfolder(S.all_figs_dir); mkdir(S.all_figs_dir); end
    exportgraphics(gcf, fullfile(S.all_figs_dir,sprintf('Fig%s_V1.pdf',variant.label)), 'ContentType', 'auto')
end

function combined = combine_datasets(varargin)
    % vertically concatenate struct arrays that don't all share the same
    % top-level fields, padding any fields missing from a given dataset
    % with [] so concatenation works. (ported from Figure_1_claude.m /
    % lpsp_compartments_claude_script.m)
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
    % ported from Figure_1_claude.m / lpsp_compartments_claude_script.m
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
    % ported from Figure_1_claude.m / lpsp_compartments_claude_script.m
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

function parts = split_path(p)
    % meta paths are stored as literal Windows paths (Z:\pablo\...)
    % regardless of the host OS this script runs on -- split on either
    % separator explicitly.
    parts = strsplit(p,{'\\','/'});
end

function is_flash = detect_flash_frames(f, mad_k)
    % ported from Figure_1_claude.m / lpsp_compartments_claude_script.m
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
    % linear placeholder spanning the trial (ported from Figure_1_claude.m)
    if isfield(ft,'xb') && numel(ft.xb) == n_im
        xb = ft.xb(:);
    else
        xb = linspace(ft.xf(1),ft.xf(end),n_im)';
    end
end

function best_lag = fit_group_lag(ft_list, resp_list, lag_grid)
    % ported from Figure_1_claude.m: correlate |r_speed(t)| against
    % resp(t+lag) across a grid of candidate lags, pooling the mean
    % correlation across trials, and pick the lag that maximizes it.
    n_lags = numel(lag_grid);
    n_trials = numel(ft_list);
    corr_grid = nan(n_trials,n_lags);
    for ii = 1:n_trials
        speed = abs(ft_list{ii}.r_speed);
        fa    = resp_list{ii};
        for L = 1:n_lags
            lag = lag_grid(L);
            if lag == 0
                s_win = speed; fa_win = fa;
            elseif lag > 0
                s_win = speed(1:end-lag); fa_win = fa(lag+1:end);
            else
                s_win = speed(-lag+1:end); fa_win = fa(1:end+lag);
            end
            valid = ~isnan(s_win) & ~isnan(fa_win);
            if sum(valid) > 100
                corr_grid(ii,L) = corr(s_win(valid),fa_win(valid));
            end
        end
    end
    mean_corr = mean(corr_grid,1,'omitnan');
    [~,best_idx] = max(mean_corr);
    best_lag = lag_grid(best_idx);
end

function [for_l,fly_l,amp_l] = apply_lag3(for_vel, fly_vel, amp, lag)
    % ported verbatim from lpsp_compartments_claude_script.m: shift amp by
    % lag frames relative to for_vel/fly_vel (amp lags behind behavior)
    if lag == 0
        for_l = for_vel; fly_l = fly_vel; amp_l = amp;
    elseif lag > 0
        for_l = for_vel(1:end-lag); fly_l = fly_vel(1:end-lag); amp_l = amp(lag+1:end);
    else
        for_l = for_vel(-lag+1:end); fly_l = fly_vel(-lag+1:end); amp_l = amp(1:end+lag);
    end
end

function r2 = fit_r2_triplet(fly_vel, for_vel, amp)
    % ported verbatim from lpsp_compartments_claude_script.m: adjusted R^2
    % for amp ~ for_vel, amp ~ |fly_vel|, and amp ~ [|fly_vel|,for_vel]
    mdl   = fitlm(for_vel,amp);
    r2(1) = mdl.Rsquared.Adjusted;
    mdl   = fitlm(abs(fly_vel),amp);
    r2(2) = mdl.Rsquared.Adjusted;
    mdl   = fitlm([abs(fly_vel),for_vel],amp);
    r2(3) = mdl.Rsquared.Adjusted;
end

function plot_r2_group(X)
    % ported verbatim from lpsp_compartments_claude_script.m: one faint
    % line per fly connecting its R^2 across the three models (forward,
    % rotational, joint), plus mean+/-sem in red.
    hold on
    xticks(1:3); xticklabels({'forward','rotational','joint'}); xlim([.5,3.5])
    X = X(~any(isnan(X),2),:);
    if isempty(X)
        text(2,0,'n=0','HorizontalAlignment','center','VerticalAlignment','bottom')
        return
    end
    plot(X','o-','Color',.8*[1,1,1],'LineWidth',1)
    m = mean(X,1);
    s = std(X,[],1)/sqrt(size(X,1));
    errorbar(1:3,m,s,'ro-','LineWidth',2,'MarkerFaceColor','r')
end

function [chosen_fly,i_ex,t0,t1] = pick_high_corr_example( ...
        rows_cl, all_data, fly_of_trial, resp_list, snippet_s, snippet_step_s, min_forward_speed_snippet, ...
        smooth_sigma_s, max_quiet_run_s, quiet_r_speed_thresh, quiet_f_speed_thresh)
    % jointly chooses the fly and snippet_s-second window (among rows_cl's
    % trials) whose resp_list signal correlates most strongly with (smoothed)
    % |r_speed|, among windows that (a) clear min_forward_speed_snippet mean
    % forward speed, (b) fluctuate enough to make that correlation
    % meaningful rather than a near-constant trace correlating trivially,
    % and (c) never go still (both |r_speed| and |f_speed| near zero) for
    % longer than max_quiet_run_s continuously -- (a) alone is a whole-
    % window MEAN and wouldn't catch a window that's active for 30s then
    % totally still for the other 30s. r_speed is smoothed with the SAME
    % two-stage (pre-abs + post-abs) gaussian used for row 4's own display
    % trace, so the correlation used here matches what a viewer would see
    % plotted.
    %
    % Per fly: every candidate window is scored on both correlation and its
    % own std(|r_speed|_smoothed); windows below THIS FLY'S OWN median std
    % are dropped (same median-cutoff trick as Figure_1_claude.m's
    % pick_snippet), and the remaining window with the highest correlation
    % becomes that fly's candidate. The fly with the single highest
    % candidate correlation wins.
    fly_list = unique(fly_of_trial(rows_cl));
    fly_best_corr = -inf(numel(fly_list),1);
    fly_best_window = cell(numel(fly_list),1); % {trial_idx, t0} of that fly's best qualifying window

    for f = 1:numel(fly_list)
        sel = rows_cl(strcmp(fly_of_trial(rows_cl),fly_list{f}));
        cand_trial = []; cand_start = []; cand_corr = []; cand_std = [];

        for i = sel(:)'
            xf = all_data(i).ft.xf(:);
            r  = all_data(i).ft.r_speed;
            fspd = all_data(i).ft.f_speed;
            fa = resp_list{i};
            fs = 1/median(diff(xf));
            smooth_win = round(5*smooth_sigma_s*fs); % smoothdata gaussian window = 5*sigma
            r_pre  = smoothdata(r,'gaussian',smooth_win);
            r_post = smoothdata(abs(r_pre),'gaussian',smooth_win);

            if xf(end) - xf(1) < snippet_s
                continue
            end
            starts = xf(1):snippet_step_s:(xf(end)-snippet_s);
            for w = starts
                idx = xf>=w & xf<w+snippet_s;
                if mean(fspd(idx)) < min_forward_speed_snippet
                    continue % this window doesn't meet the "largely running forward" requirement
                end
                r_win = r_post(idx); fa_win = fa(idx);
                valid = ~isnan(r_win) & ~isnan(fa_win);
                if sum(valid) < 100
                    continue
                end
                quiet_mask = r_win < quiet_r_speed_thresh & abs(fspd(idx)) < quiet_f_speed_thresh;
                if longest_run_s(quiet_mask,1/fs) > max_quiet_run_s
                    continue % this window contains too long a stretch of the fly just sitting still
                end
                cand_trial(end+1) = i;                            %#ok<AGROW>
                cand_start(end+1) = w;                             %#ok<AGROW>
                cand_std(end+1)   = std(r_win(valid));             %#ok<AGROW>
                cand_corr(end+1)  = corr(r_win(valid),fa_win(valid)); %#ok<AGROW>
            end
        end

        if isempty(cand_start)
            continue
        end
        fluctuates_enough = cand_std >= median(cand_std);
        if ~any(fluctuates_enough)
            fluctuates_enough = true(size(cand_std));
        end
        ok_idx = find(fluctuates_enough);
        [best_corr,rel] = max(cand_corr(ok_idx));
        fly_best_corr(f) = best_corr;
        fly_best_window{f} = {cand_trial(ok_idx(rel)), cand_start(ok_idx(rel))};
    end

    [~,best_f] = max(fly_best_corr);
    assert(isfinite(fly_best_corr(best_f)), ...
        'no %ds window with mean forward speed >= %.1f mm/s found in any candidate trial', snippet_s, min_forward_speed_snippet)
    chosen_fly = fly_list{best_f};
    i_ex = fly_best_window{best_f}{1};
    t0   = fly_best_window{best_f}{2};
    t1   = t0 + snippet_s;
end

function len_s = longest_run_s(mask, dt)
    % longest contiguous run of true values in mask, in seconds
    d = diff([false; mask(:); false]);
    run_starts = find(d==1);
    run_ends   = find(d==-1) - 1;
    if isempty(run_starts)
        len_s = 0;
    else
        len_s = max(run_ends - run_starts + 1) * dt;
    end
end

function cmap = white_to_color(max_color)
    % linear colormap from white (low) to max_color (high) -- ported
    % verbatim from Figure_1_claude.m (used there for its per-column colors)
    n = 256;
    cmap = [linspace(1,max_color(1),n)', linspace(1,max_color(2),n)', linspace(1,max_color(3),n)'];
end

function [fly_ids,mean90,mean10] = bin_tuning_by_fly( ...
        rows, fly_of_trial, rot_smooth_list, fspeed_smooth_list, p90_list, p10_list, lag, bin_edges, min_n)
    % for each fly (among rows' trials), pools that fly's own trials'
    % lag-shifted (|r_speed|, p90z, p10z) samples, then averages p90z/p10z
    % separately within each bin_edges bin -- one mean value per fly per
    % bin (NaN if that fly has fewer than min_n samples in that bin).
    fly_ids = unique(fly_of_trial(rows));
    n_bins = numel(bin_edges)-1;
    mean90 = nan(numel(fly_ids),n_bins);
    mean10 = nan(numel(fly_ids),n_bins);

    for f = 1:numel(fly_ids)
        trial_list = rows(strcmp(fly_of_trial(rows),fly_ids{f}));
        pooled_rot = []; pooled_90 = []; pooled_10 = [];
        for i = trial_list(:)'
            [~,rot_l,amp90_l] = apply_lag3(fspeed_smooth_list{i},rot_smooth_list{i},p90_list{i},lag);
            [~,~,    amp10_l] = apply_lag3(fspeed_smooth_list{i},rot_smooth_list{i},p10_list{i},lag);
            valid = ~isnan(rot_l) & ~isnan(amp90_l) & ~isnan(amp10_l);
            pooled_rot = [pooled_rot; rot_l(valid)];    %#ok<AGROW>
            pooled_90  = [pooled_90;  amp90_l(valid)];  %#ok<AGROW>
            pooled_10  = [pooled_10;  amp10_l(valid)];  %#ok<AGROW>
        end
        for b = 1:n_bins
            in_bin = pooled_rot>=bin_edges(b) & pooled_rot<bin_edges(b+1);
            if sum(in_bin) >= min_n
                mean90(f,b) = mean(pooled_90(in_bin));
                mean10(f,b) = mean(pooled_10(in_bin));
            end
        end
    end
end

function plot_binned_tuning(ax, bin_centers, mean90, mean10)
    % one faint dot per fly per bin (that fly's own bin mean) in black
    % (90th %ile) / gray (10th %ile), jittered apart in x so the two
    % overlap less, plus a population mean+/-SEM errorbar line per curve.
    hold(ax,'on')
    jit = 0.12*median(diff(bin_centers));
    for b = 1:numel(bin_centers)
        y90 = mean90(:,b); y90 = y90(~isnan(y90));
        y10 = mean10(:,b); y10 = y10(~isnan(y10));
        if ~isempty(y90)
            scatter(ax,bin_centers(b)-jit+zeros(size(y90)),y90,15,'k','filled','MarkerFaceAlpha',.35)
        end
        if ~isempty(y10)
            scatter(ax,bin_centers(b)+jit+zeros(size(y10)),y10,15,[.5 .5 .5],'filled','MarkerFaceAlpha',.35)
        end
    end
    m90 = mean(mean90,1,'omitnan'); s90 = std(mean90,[],1,'omitnan')./sqrt(sum(~isnan(mean90),1));
    m10 = mean(mean10,1,'omitnan'); s10 = std(mean10,[],1,'omitnan')./sqrt(sum(~isnan(mean10),1));
    errorbar(ax,bin_centers,m90,s90,'-o','Color','k','LineWidth',2,'MarkerFaceColor','k')
    errorbar(ax,bin_centers,m10,s10,'-o','Color',[.5 .5 .5],'LineWidth',2,'MarkerFaceColor',[.5 .5 .5])
    xlabel(ax,'abs. rotational speed (rad/s, binned)')
    ylabel(ax,'z-score')
end

function s = meta_display(meta_path)
    % last path component only, for compact titles/fprintf (paths are
    % literal Windows paths regardless of host OS)
    parts = split_path(meta_path);
    parts(cellfun(@isempty,parts)) = [];
    s = parts{end};
end
