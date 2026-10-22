%% Figure_3_claude
% Combines analyses from "data scripts/lpsp_p2x2_claude.m" (ATP
% iontophoresis onto the PB in LPsP>P2X2 vs empty>P2X2 flies -- per-stim
% ejection detection, bump-position response aligned to stim onset) and
% "data scripts/lpsp_p2x2_walking_script.m" (bump-mobility / gain-like
% path-length ratio before vs after the perturbation trial, split by
% genotype and light condition), following the Figure_1_claude.m /
% Figure_2_claude.m conventions in this folder (cell-mode sections, ported
% per-trial helper functions with header comments, one shared data-loading
% pass, exportgraphics at the end).
%
% SOURCE DATASET: lpsp_p2x2_reredo_20260518.mat only, per explicit request
% -- this is the dataset with a clean, complete ft.stims (ATP-ejection TTL)
% record (101/101 trials, 4 stims each) and an unambiguous "fly N _
% <empty|lpsp>" folder-based genotype tag (confirmed directly against all
% 101 trials when investigating this dataset for the prior comparison
% report). lpsp_p2x2_revisit_20250829.mat was set aside for now: it has no
% ft.stims at all, so its ejection timing would need a different,
% unvalidated peak-detector before it could support the same kind of panel.

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
addpath(fullfile(repo_root,'circ_stats'))
data_dir = fullfile(repo_root,'.data'); % Windows machine keeps datasets in .data/, the Mac in data/
if ~isfolder(data_dir); data_dir = fullfile(repo_root,'data'); end
export_dir = fullfile(repo_root,'ugly_figures','exports');
all_figs_dir = fullfile(repo_root,'ugly_figures','all_figs');

%% load lpsp_p2x2_reredo data
source_file = 'lpsp_p2x2_reredo_20260518.mat';
tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype + fly ID per trial, from the "fly N _ <genotype>" folder in all_data.meta
% ported from data scripts/lpsp_p2x2_claude.m -- see that script's header
% for why the tag sits one level above the trial folder here (unlike
% kir/tnt, where it's in the trial folder name itself).
genotype = cell(n_trials,1);
fly_id   = cell(n_trials,1);
for i = 1:n_trials
    fname = trial_fly_folder_name(all_data(i).meta);
    if contains(fname,'lpsp','IgnoreCase',true)
        genotype{i} = 'lpsp>p2x2';
    elseif contains(fname,'empty','IgnoreCase',true)
        genotype{i} = 'empty>p2x2';
    else
        genotype{i} = '';
    end
    fly_id{i} = trial_fly_id(all_data(i).meta);
end

%% precompute confirmed-ejection stim events for every trial -- shared by
% the example-trial selection, row 3, and the scratch diagnostic below.
% (This consolidates what used to be three separate copies of the same
% per-stim ejection-detection loop into one pass, which is also what made
% it easy to add the per-stim side classification below without having to
% keep three copies of that logic in sync.) Each confirmed ejection (same
% per-stim z-score detector as lpsp_p2x2_claude.m, peak_factor x baseline
% SD above baseline just after onset) gets: its onset time, its own
% peak-window (0-3s post-onset) per-wedge atp.d spatial profile (both the
% raw 32-element profile and the resulting right_frac scalar), and its
% bump-position (mu) / self-motion heading / im.rho traces, all
% interpolated onto one shared time grid and zeroed at that stim's own
% onset.
peri_win    = 5;
base_win    = [-peri_win,0];
peak_win    = [0,3];
peak_factor = 3;
t_grid = linspace(-5,30,351); % covers both the -5..15s window the example-trial selection scores below and the -2..30s window row 3/the scratch figure display

trial_events = cell(n_trials,1);
for i = 1:n_trials
    trial_events{i} = get_trial_stim_events(all_data(i), t_grid, peri_win, base_win, peak_win, peak_factor);
end
n_ejections_confirmed_all = cellfun(@numel,trial_events);
fprintf('\n%d/%d trials have >=1 confirmed ejection (%d confirmed ejections total)\n', ...
    sum(n_ejections_confirmed_all>0), n_trials, sum(n_ejections_confirmed_all));

%% which PB hemisphere each trial's ATP pipette targeted -- called from the
% MOST CONFIDENT individual ejection, then applied to the whole trial. Per
% explicit request: ejections within one trial are always the same
% physical side, so if e.g. the 3rd stim gives an unambiguous spatial
% profile, that call can be trusted for the whole trial even if another
% stim in the same trial looks ambiguous on its own.
%
% This replaces an earlier version that averaged atp.f (raw fluorescence)
% over the WHOLE TRIAL: stims are ~60s apart and only eject for a few
% seconds, so that whole-trial average pools in ~57s of inter-stim
% baseline for every ~3s of actual signal, diluting the side call for any
% trial where even one of the 4 stims was weak. Per stim here, right_frac
% instead comes from atp.d (dF/F, already baseline-subtracted per wedge)
% during just that ONE stim's own 0-3s peak window (computed once, in
% get_trial_stim_events); per trial, the confirmed-ejection stim with the
% largest |right_frac-0.5| (the clearest individual call) decides the
% whole trial's side. See the scratch figure exported below for a
% step-by-step look at this, including how often it changes the call
% relative to the old whole-trial method.
trial_is_right   = nan(n_trials,1); % NaN = no confirmed ejection this trial, no call possible
trial_right_frac = nan(n_trials,1); % the deciding stim's own right_frac
for i = 1:n_trials
    ev = trial_events{i};
    if isempty(ev)
        continue
    end
    rfs = [ev.right_frac];
    [~,most_confident] = max(abs(rfs-0.5));
    trial_right_frac(i) = rfs(most_confident);
    trial_is_right(i)   = trial_right_frac(i) > 0.5;
end

% comparison against the old whole-trial atp.f method, kept only to report
% how many trials the new per-stim method reclassifies (see scratch figure)
whole_trial_right_frac = nan(n_trials,1);
for i = 1:n_trials
    f_sum = sum(all_data(i).atp.f,2);
    whole_trial_right_frac(i) = mean(f_sum(1:16)) / (mean(f_sum(1:16)) + mean(f_sum(17:32)));
end
old_is_right  = whole_trial_right_frac > 0.5;
old_ambiguous = whole_trial_right_frac>0.4 & whole_trial_right_frac<0.6;
have_call = ~isnan(trial_is_right);
n_reclassified = sum(have_call & (trial_is_right ~= old_is_right));
fprintf('\nstim-side call: %d/%d trials got a call (the rest had 0 confirmed ejections)\n', sum(have_call), n_trials);
fprintf('  old whole-trial-atp.f method: %d/%d trials fell in the ambiguous 0.4-0.6 band\n', sum(old_ambiguous), n_trials);
fprintf('  new most-confident-stim method disagrees with the old method on %d/%d called trials\n', n_reclassified, sum(have_call));

%% pick two example LPsP (not empty) STIMS -- one per side -- that each
% show the bump MOVING after ejection AND have a clean, sizeable atp
% signal on that SAME stim (not just somewhere else in the same trial).
%
% Selection works per INDIVIDUAL STIM, not per trial: an earlier version
% scored/filtered at the trial level (a trial's BEST stim had to clear the
% amplitude bar), then always displayed that trial's stim #2 regardless of
% which stim had actually driven the trial's score -- for the picked
% opposite-side trial, stim #1 had a clean amp=2.64 but the displayed
% stim #2 was actually its weakest (amp=0.78), so the panel shown didn't
% match the trial's own qualification. Scoring and thresholding the exact
% stim that gets displayed avoids that mismatch.
%
% Bump movement score: for each confirmed-ejection stim, both mu and a
% self-motion-integrated heading (cumsum(r_speed) -- same "heading"
% construction lpsp_p2x2_revisit.m/lpsp_p2x2_walking_script.m use, robust
% in both closed-loop and dark trials, unlike -ft.cue which only tracks
% the fly in closed loop) are zeroed at that stim's own onset; under normal
% (non-perturbed) EPG tracking mu and heading move together 1:1, so a stim
% that visibly KICKS the bump off of where the fly's own rotation would
% predict shows up as a large |mu-heading divergence|, scored over the
% SETTLED 8-15s post-onset window (gated on mean im.rho>rho_gate so a
% noisy/lost bump can't win by accident, and using the later window rather
% than the first second or two, which is dominated by a broad,
% unstructured pan-PB calcium brightening that temporarily destabilizes
% the bump estimate).
% Amplitude filter: that SAME stim's own atp_peak_amp (whole-PB atp.d dF/F
% peak 0-5s post-onset minus its own -2-0s baseline, from get_trial_stim_events)
% must clear atp_peak_amp_thresh, so the displayed panel is guaranteed
% visually clean, not just movement-scored. Trials are also required to
% have >=2 confirmed ejections overall (not just this one stim), so the
% example isn't a one-off within an otherwise-failed trial.
settle_win = [8,15];
rho_gate   = 0.25;

lpsp_rows = find(strcmp(genotype,'lpsp>p2x2'));

% atp_peak_amp is now dF/F, not z-score (per explicit request -- see
% get_trial_stim_events), so the old z-score-calibrated threshold (1.2) no
% longer applies. Re-checked directly on this dataset: dF/F atp_peak_amp
% across all lpsp confirmed-ejection stims (n=200) has min=0.09, p25=1.81,
% median=2.64, p75=3.83, max=10.94 -- 2.0 sits just above the 25th
% percentile, comfortably separating a lower tail of weak/partial
% ejections from the bulk of clearly-clean ones.
atp_peak_amp_thresh = 2.0;
cand_trial = []; cand_stim = []; cand_score = []; cand_side = [];
for k = 1:numel(lpsp_rows)
    i = lpsp_rows(k);
    if n_ejections_confirmed_all(i) < 2; continue; end
    ev = trial_events{i};
    for s = 1:numel(ev)
        if ev(s).atp_peak_amp < atp_peak_amp_thresh; continue; end
        divergence = ev(s).mu_rel - ev(s).head_rel;
        settle_idx = t_grid>settle_win(1) & t_grid<=settle_win(2);
        if mean(ev(s).rho_grid(settle_idx),'omitnan') <= rho_gate; continue; end
        cand_trial(end+1) = i; %#ok<AGROW>
        cand_stim(end+1)  = s; %#ok<AGROW>
        cand_score(end+1) = abs(median(divergence(settle_idx),'omitnan')); %#ok<AGROW>
        cand_side(end+1)  = trial_is_right(i); %#ok<AGROW>
    end
end
fprintf('\n%d individual confirmed-ejection stims clear atp_peak_amp>=%.1f and the rho gate (candidate pool for both examples)\n', numel(cand_trial), atp_peak_amp_thresh);

[~,best] = max(cand_score);

% the reference side for row 3's flip convention below: whatever side THIS
% (top-scoring, before any manual rank override) stim actually was, per
% explicit request -- not a fixed "always call it right" choice.
reference_is_right = logical(cand_side(best));
if reference_is_right; ref_side_str = 'RIGHT'; else; ref_side_str = 'LEFT'; end

% rank every candidate within its own side by the same post-ejection
% bump-movement score used above, so a specific rank can be hand-picked
% below (from Figure_3_claude_scratch_example_gallery.png) instead of
% always taking the single top-scoring stim per side.
ref_pool = find(cand_side==reference_is_right);
[~,order_ref] = sort(cand_score(ref_pool),'descend');
ref_pool = ref_pool(order_ref);

% example pick: rank 3 on the reference/LEFT side, per explicit request --
% the top-ranked stim (rank 1) didn't look as clean in the gallery.
example_rank_ref = 3;

best = ref_pool(example_rank_ref);
example_i = cand_trial(best);
example_s = cand_stim(best);
best_score = cand_score(best);

fprintf('\nexample LPsP>P2X2 trial: fly %s, trial "%s", stim %d (rank %d of %d on the %s side)\n', fly_id{example_i}, meta_display(all_data(example_i).meta), example_s, example_rank_ref, numel(ref_pool), ref_side_str);
fprintf('  post-ejection bump displacement = %.2f rad, atp_peak_amp = %.2f\n', best_score, trial_events{example_i}(example_s).atp_peak_amp);
fprintf('  example stim side: %s -- rows 4-5 below are flipped to match this side\n', ref_side_str);

%% pick the SECOND example stim: an EMPTY>P2X2 control, SAME stim side as
% the LPsP example above (per explicit request -- was previously a second
% LPsP stim from the opposite hemisphere). Same candidate filters as the
% LPsP pool (>=2 confirmed ejections in the trial, atp_peak_amp >=
% atp_peak_amp_thresh, settled-window rho > rho_gate), but ranked by how
% closely its atp_peak_amp matches the LPsP example's own, NOT by the
% bump-movement score: ranking a control by "most bump movement" would
% pick the least representative empty stim, and ranking by "least
% movement" would cherry-pick the other way. Amplitude matching makes the
% two columns a like-for-like ejection, so any difference in the bump
% response is down to genotype.
% Per explicit request the pick should also have a very SMOOTH whole-PB
% atp trace (the amplitude-matched rank 1 was visibly wiggly), so: keep
% only stims whose atp_peak_amp is within +-emp_amp_tol of the LPsP
% example's, then rank those by trace roughness (atp_trace_roughness: SD
% of the -2..20s whole-PB atp dF/F trace around a 1s-gaussian-smoothed
% copy, as a fraction of its peak) -- smoothest first. Stims outside the
% amplitude band follow, still by amplitude match. example_rank_empty
% picks a rank by eye from column 2 of the gallery below: rank 3 (fly 4
% empty, 20260429, stim 4) has the closest peak (~9 dF/F) and slow decay
% to the LPsP example's own trace.
emp_amp_tol = 0.25;
example_rank_empty = 3;

ref_amp = trial_events{example_i}(example_s).atp_peak_amp;
empty_rows = find(strcmp(genotype,'empty>p2x2'));
emp_trial = []; emp_stim = []; emp_score = []; emp_amp_err = []; emp_rough = [];
for k = 1:numel(empty_rows)
    i = empty_rows(k);
    if n_ejections_confirmed_all(i) < 2 || trial_is_right(i)~=reference_is_right; continue; end
    ev = trial_events{i};
    for s = 1:numel(ev)
        if ev(s).atp_peak_amp < atp_peak_amp_thresh; continue; end
        settle_idx = t_grid>settle_win(1) & t_grid<=settle_win(2);
        if mean(ev(s).rho_grid(settle_idx),'omitnan') <= rho_gate; continue; end
        emp_trial(end+1)   = i; %#ok<AGROW>
        emp_stim(end+1)    = s; %#ok<AGROW>
        emp_score(end+1)   = abs(median(ev(s).mu_rel(settle_idx) - ev(s).head_rel(settle_idx),'omitnan')); %#ok<AGROW>
        emp_amp_err(end+1) = abs(ev(s).atp_peak_amp - ref_amp); %#ok<AGROW>
        emp_rough(end+1)   = atp_trace_roughness(all_data(i), ev(s).onset_t); %#ok<AGROW>
    end
end
emp_in_band = find(emp_amp_err <= emp_amp_tol*ref_amp);
emp_out_band = find(emp_amp_err > emp_amp_tol*ref_amp);
[~,o_in]  = sort(emp_rough(emp_in_band),'ascend');
[~,o_out] = sort(emp_amp_err(emp_out_band),'ascend');
emp_pool = [emp_in_band(o_in), emp_out_band(o_out)];
fprintf('\n%d/%d empty candidate stims within +-%.0f%% of the LPsP example''s atp amplitude (ranked by trace smoothness)\n', numel(emp_in_band), numel(emp_amp_err), 100*emp_amp_tol);

best2 = emp_pool(example_rank_empty);
example_i2 = emp_trial(best2);
example_s2 = emp_stim(best2);
best_score2 = emp_score(best2);

fprintf('\nsecond example EMPTY>P2X2 trial (same side, %s): fly %s, trial "%s", stim %d (rank %d of %d: amplitude-matched, then by atp trace smoothness)\n', ref_side_str, fly_id{example_i2}, meta_display(all_data(example_i2).meta), example_s2, example_rank_empty, numel(emp_pool));
fprintf('  post-ejection bump displacement = %.2f rad, atp_peak_amp = %.2f (LPsP example: %.2f), atp trace roughness = %.4f\n', best_score2, trial_events{example_i2}(example_s2).atp_peak_amp, ref_amp, emp_rough(best2));

%% ===================== SCRATCH: gallery of candidate example stims (row 1) =====================
% per explicit request ("show me some different options, I think we can
% find cleaner still") -- same candidate pool as the two picks above
% (cand_trial/cand_stim/cand_score/cand_side, already filtered on
% atp_peak_amp>=atp_peak_amp_thresh and the rho gate), just laid out as a
% ranked grid of EPG-calcium heatmaps (same zoom/clim convention as row 1
% of the main figure below) instead of always taking the single top-score
% stim. Column 1: LPsP candidates on the reference side, ranked by the
% post-ejection bump-movement score. Column 2: empty>P2X2 candidates on
% the SAME side, amplitude-matched to the LPsP pick, then ranked by atp trace smoothness. The
% current picks (example_rank_ref/example_rank_empty above) are marked.
% Reuses ref_pool/emp_pool computed above (not recomputed here).
gallery_n = 8;
gallery_zoom_pre_s  = 2;  % s before onset -- same convention as the main figure's zoom_pre_s below (defined later in the script, so not reused directly here)
gallery_zoom_post_s = 20; % s after onset -- same as zoom_post_s below
gallery_gcamp_color = [0 .2 1]; % same EPG blue as gcamp_color below (defined later in the script, so not reused directly here)

gallery_ref_pool = ref_pool(1:min(gallery_n,numel(ref_pool)));
gallery_emp_pool = emp_pool(1:min(gallery_n,numel(emp_pool)));

fig_gallery = figure('color','w','Position',[50 50 900 170*gallery_n]); clf
tg = tiledlayout(fig_gallery,gallery_n,2,'TileSpacing','compact','Padding','compact');

gallery_cols   = {gallery_ref_pool, gallery_emp_pool};
gallery_labels = {[ref_side_str ', LPsP'], [ref_side_str ', empty']};
gallery_pick   = [best, best2];
gallery_trial  = {cand_trial, emp_trial};
gallery_stim   = {cand_stim, emp_stim};
gallery_score  = {cand_score, emp_score};
for r = 1:gallery_n
    for col = 1:2
        ax = nexttile(tg); hold(ax,'on')
        pool = gallery_cols{col};
        if r > numel(pool)
            axis(ax,'off'); continue
        end
        k = pool(r);
        i = gallery_trial{col}(k); s = gallery_stim{col}(k);
        ex = all_data(i);
        xb = ex.ft.xb;
        stim_t = trial_events{i}(s).onset_t;
        zoom_t0 = stim_t - gallery_zoom_pre_s; zoom_t1 = stim_t + gallery_zoom_post_s;
        zoom_idx = xb>=zoom_t0 & xb<=zoom_t1;
        im_alpha = unwrap(ex.im.alpha);
        im_clim = [prctile(min(ex.im.z(:,zoom_idx),[],1),5), prctile(max(ex.im.z(:,zoom_idx),[],1),95)];

        imagesc(ax,xb,im_alpha,ex.im.z,im_clim)
        colormap(ax,white_to_color(gallery_gcamp_color))
        plot(ax,[stim_t,stim_t],[min(im_alpha),max(im_alpha)],':k','LineWidth',1)
        xlim(ax,[zoom_t0,zoom_t1]); ylim(ax,[min(im_alpha),max(im_alpha)])

        pick_str = ''; if k==gallery_pick(col); pick_str = ' <-- CURRENT PICK'; end
        title(ax,sprintf('rank %d (%s): fly %s, stim %d\nscore=%.2f rad, atp=%.2f%s', ...
            r, gallery_labels{col}, fly_id{i}, s, gallery_score{col}(k), trial_events{i}(s).atp_peak_amp, pick_str), ...
            'FontSize',7,'Interpreter','none')
        if col==1; ylabel(ax,'PB angle (rad)'); end
        if r==gallery_n; xlabel(ax,'time (s)'); end
    end
end
title(tg,sprintf('candidate gallery for Figure 3 rows 1-2 (top %d; LPsP ranked by post-ejection bump displacement, empty by atp amplitude match then trace smoothness)',gallery_n))

exportgraphics(fig_gallery, fullfile(export_dir,'Figure_3_claude_scratch_example_gallery.png'), 'Resolution', 200)

%% plot: EPG calcium (top), atp channel (middle), and this trial's own
% whole-PB average atp trace (bottom-of-the-example-rows), zoomed to the
% 2nd stim (-2s to +30s around its onset) -- one column per example trial,
% so the two examples (opposite stim sides) sit side by side. No bump/cue
% overlay -- just the raw heatmaps/trace. EPG calcium stays z-scored (im.z,
% per-trial-normalized, matches Figure_1_claude.m/Figure_2_claude.m's own
% convention for that channel); the atp channel uses dF/F (atp.d) instead
% of z-score, per explicit request -- z-scoring re-normalizes away the
% actual ejection amplitude, which is the thing worth comparing directly
% here. Color scale for each heatmap uses the same percentile rule as
% Figure_1_claude.m/Figure_2_claude.m (5th/95th percentile of the per-frame
% trough/peak), computed from THIS ZOOMED WINDOW's own data (not the whole
% trial) so contrast isn't diluted by the other 3 stims' worth of data
% outside the crop -- and computed SEPARATELY per example column, since the
% two trials' own signal ranges differ.
zoom_pre_s  = 2;  % s before the zoomed stim's onset
zoom_post_s = 20; % s after the zoomed stim's onset -- 20s per explicit request (was 30s)

gcamp_color = [0 .2 1];   % dark blue, matches Figure_1_claude.m's EPG calcium color (group_max_color(1,:)), per explicit request (was teal [0 .7 .7])
atp_color   = [.8 0 0];   % dark red, visually distinct from the gcamp row

% z-planes summed for the row-1 mean images, per explicit request: planes
% 1-2 (the most superficial, 4 um steps going deeper -- scan_params.json
% z_positions) hold bright cell bodies above the PB that distract from it
% in the full z-sum; the PB itself runs through planes 3-10 in both
% example flies. Display only -- im/atp signals still use all planes.
mean_img_planes = 3:10;

example_idx      = [example_i, example_i2];
example_stim_idx = [example_s, example_s2]; % which of that trial's confirmed-ejection stims was actually selected -- zoom into THIS one, not a fixed stim number
example_col_label = {ref_side_str, ref_side_str};
example_geno_label = {'LPsP>P2X2', 'empty>P2X2'};
example_disp = struct('ex',{},'xb',{},'stim_t',{},'zoom_t0',{},'zoom_t1',{},'im_alpha',{},'atp_alpha',{},'im_clim',{},'atp_clim',{},'atp_avg_full',{},'mean_rgb',{});
for c = 1:2
    i = example_idx(c);
    ex_c = all_data(i);
    xb_c = ex_c.ft.xb;

    stim_t_c = trial_events{i}(example_stim_idx(c)).onset_t; % the exact stim this trial/column was selected for
    zoom_t0_c = stim_t_c - zoom_pre_s;
    zoom_t1_c = stim_t_c + zoom_post_s;
    zoom_idx_c = xb_c >= zoom_t0_c & xb_c <= zoom_t1_c;

    example_disp(c).ex           = ex_c;
    example_disp(c).xb           = xb_c;
    example_disp(c).stim_t       = stim_t_c;
    example_disp(c).zoom_t0      = zoom_t0_c;
    example_disp(c).zoom_t1      = zoom_t1_c;
    example_disp(c).im_alpha     = unwrap(ex_c.im.alpha);
    example_disp(c).atp_alpha    = unwrap(ex_c.atp.alpha);
    example_disp(c).im_clim      = [prctile(min(ex_c.im.z(:,zoom_idx_c),[],1),5),  prctile(max(ex_c.im.z(:,zoom_idx_c),[],1),95)];
    example_disp(c).atp_clim     = [prctile(min(ex_c.atp.d(:,zoom_idx_c),[],1),5), prctile(max(ex_c.atp.d(:,zoom_idx_c),[],1),95)];
    example_disp(c).atp_avg_full = mean(ex_c.atp.d,1); % row 3: whole-PB average atp dF/F over time
    example_disp(c).mean_rgb     = trial_mean_rgb(ex_c.meta, gcamp_color, atp_color, mean_img_planes); % row 1 left: whole-trial mean image, gcamp blue + atp red, rotated so the PB runs vertically
end

%% row 4: average bump-position trace across ALL confirmed ejections, all
% flies, one line +/- SEM per genotype, -2s to +30s around stim onset.
% Adapted from lpsp_p2x2_claude.m's own "bump movement collapsed across
% stim side" figure (its figure 7): for every confirmed-ejection stim
% (already found in trial_events above), the bump position (mu_rel, zeroed
% at that stim's own onset) is NEGATED if the trial's side call
% (trial_is_right) doesn't match reference_is_right -- so every trace,
% regardless of which physical hemisphere actually got the pipette that
% day, is flipped onto the SAME convention as the example trial shown in
% rows 1-2 above (this trial happens to be a RIGHT-side stim, so rows 1-2
% and row 3 both read "positive = the bump moved the way a right-side stim
% would push it" -- had the example been a left-side stim, row 3 would
% instead flip everything onto a left-side convention, so its sign always
% visually matches what rows 1-2 show). Without any flip, right- and
% left-stim trials would tend to cancel each other out when pooled, since
% the PB's two hemispheres are mirror images.
% Pooling: each stim contributes one flipped trace; a fly's own qualifying
% stims (possibly from several trials) are averaged together first, THEN
% flies are averaged within a genotype for the plotted mean/SEM -- same
% one-point-per-fly convention used throughout this codebase (e.g.
% lpsp_kir_claude.m, Figure_1_claude.m's row 3) so a fly with more
% confirmed ejections doesn't get more weight than one with fewer. Both
% light conditions (closed loop + dark) are pooled together here since the
% request didn't ask for that split.
trace_pre_s  = -2; trace_post_s = 20; % 20s post-stim per explicit request (was 30s) -- rows 4-5
trace_idx = t_grid>=trace_pre_s & t_grid<=trace_post_s;
t_trace = t_grid(trace_idx);

geno_order = {'empty>p2x2','lpsp>p2x2'};
geno_color = [.5 .5 .5; gcamp_color]; % empty = gray (kept throughout rows 4-5), lpsp = the same EPG-calcium blue used in rows 1-3

% row 4 col 2 / row 5 col 2 (EPG panels) also use lpsp = the EPG-calcium
% blue; row 5 col 1 (atp panel) uses lpsp = the atp red (matching row 2/3)
% -- empty stays plain gray in every one of these panels (not a tinted
% shade of the channel color).
epg_geno_color = [.5 .5 .5; gcamp_color]; % empty = gray, lpsp = EPG-calcium blue
atp_geno_color = [.5 .5 .5; atp_color];   % empty = gray, lpsp = atp red

fly_trace_by_geno = cell(1,numel(geno_order)); % {genotype}: n_flies_g x numel(t_trace), one row per fly

for gi = 1:numel(geno_order)
    rows = find(strcmp(genotype,geno_order{gi}));
    these_flies = unique(fly_id(rows));
    fly_traces = nan(numel(these_flies),numel(t_trace));

    for f = 1:numel(these_flies)
        trial_list = rows(strcmp(fly_id(rows),these_flies{f}));
        stim_traces = []; % one row per confirmed-ejection stim, pooled across this fly's own trials

        for i = trial_list(:)'
            ev = trial_events{i};
            for s = 1:numel(ev)
                mu_rel = ev(s).mu_rel(trace_idx);
                if trial_is_right(i) ~= reference_is_right
                    mu_rel = -mu_rel; % flip this trial's side onto the example trial's own side
                end
                stim_traces(end+1,:) = mu_rel; %#ok<AGROW>
            end
        end

        if ~isempty(stim_traces)
            fly_traces(f,:) = mean(stim_traces,1,'omitnan');
        end
    end
    fly_trace_by_geno{gi} = fly_traces;
end

geno_mean = nan(numel(geno_order),numel(t_trace));
geno_sem  = nan(numel(geno_order),numel(t_trace));
geno_n    = zeros(1,numel(geno_order));
for gi = 1:numel(geno_order)
    M = fly_trace_by_geno{gi};
    have = ~all(isnan(M),2);
    M = M(have,:);
    geno_n(gi) = size(M,1);
    geno_mean(gi,:) = mean(M,1,'omitnan');
    geno_sem(gi,:)  = std(M,[],1,'omitnan') ./ sqrt(sum(~isnan(M),1));
end
fprintf('\nrow 4 (bump position vs. stim, flipped to match the example trial''s %s side): n=%d empty flies, n=%d lpsp flies\n', ref_side_str, geno_n(1), geno_n(2));

%% row 3 col 2: LPsP bump-position response split by stim number WITHIN a
% trial (1-4), one line per stim number -- per explicit request, replacing
% the whole-PB atp population panel. Ported from Figure_3_supplements.m
% (Figure 3 S1, top-right traces panel); see that script for the full
% split (both genotypes, session order, per-fly metrics) and the matching
% atp / EPG-calcium splits that rule out a bigger ejection or a bigger EPG
% response at later stims. Every trial has exactly 4 TTL stims and a trial
% with any confirmed ejection has all 4 confirmed, so ev(s) IS stim s; the
% match against the TTL onsets below just guards that assumption.
% Unlike the col 1 panel (zeroed at the onset sample), each trace here is
% baseline-subtracted with its own -2..0 s mean, as in the supplement.
% Same side flip and per-fly-then-across-flies pooling as col 1.
n_stim_per_trial = 4;
stimnum_base_ix = t_grid>=-2 & t_grid<0;
lpsp_trial_rows = find(strcmp(genotype,'lpsp>p2x2'));
lpsp_fly_list = unique(fly_id(lpsp_trial_rows));
stimnum_fly_traces = nan(numel(lpsp_fly_list),numel(t_trace),n_stim_per_trial);
for f = 1:numel(lpsp_fly_list)
    trial_list = lpsp_trial_rows(strcmp(fly_id(lpsp_trial_rows),lpsp_fly_list{f}));
    by_stim = cell(1,n_stim_per_trial);
    for i = trial_list(:)'
        st = logical(all_data(i).ft.stims(:));
        ttl_t = all_data(i).ft.xf(diff([false;st])==1);
        ev = trial_events{i};
        for e = 1:numel(ev)
            [dmin,sn] = min(abs(ttl_t - ev(e).onset_t));
            if dmin > 1e-6 || sn > n_stim_per_trial; continue; end
            mu_rel = ev(e).mu_rel - mean(ev(e).mu_rel(stimnum_base_ix),'omitnan');
            if trial_is_right(i) ~= reference_is_right; mu_rel = -mu_rel; end
            by_stim{sn}(end+1,:) = mu_rel(trace_idx);
        end
    end
    for sn = 1:n_stim_per_trial
        if ~isempty(by_stim{sn}); stimnum_fly_traces(f,:,sn) = mean(by_stim{sn},1,'omitnan'); end
    end
end
stimnum_mean = nan(n_stim_per_trial,numel(t_trace));
stimnum_sem  = nan(n_stim_per_trial,numel(t_trace));
stimnum_n    = zeros(1,n_stim_per_trial);
for sn = 1:n_stim_per_trial
    M = stimnum_fly_traces(:,:,sn); M = M(~all(isnan(M),2),:);
    stimnum_n(sn) = size(M,1);
    stimnum_mean(sn,:) = mean(M,1,'omitnan');
    stimnum_sem(sn,:)  = std(M,[],1,'omitnan') ./ sqrt(sum(~isnan(M),1));
end
fprintf('row 3 col 2 (lpsp bump position by stim number): n=%s flies per stim\n', mat2str(stimnum_n));

%% row 6: bump mobility (mov_ratio = bump path length / heading path
% length during walking bouts) before vs. after the ATP-ejection
% perturbation trial, DARK trials only -- ported from data scripts/
% lpsp_p2x2_walking_script.m (sections 2-4), which frames this exact
% pre/post comparison as its own candidate "Fig 3 new row" (see that
% script's section 4b header). Uses a SEPARATE dataset
% (lpsp_p2x2_walking_20260728, not lpsp_p2x2_reredo above) -- same two
% genotypes, but these are longer straight-walking sessions built for
% bout-level bump-mobility regression, not the ATP-response imaging
% trials rows 1-5 use. Per explicit request: plain gray dots/connecting
% lines (one line per fly, pre->post), no per-group color coding, split
% into two panels by genotype (empty left, lpsp right).
%
% Pipeline (ported as-is, mac path-split fixed the same way as
% trial_fly_folder_name/trial_fly_id above):
%  1) tag genotype/lighting per trial (contains '_empty_'/'_lpsp_' in
%     meta, contains 'background' in ft.pattern for dark)
%  2) detect which trial IS the perturbation trial (has_pulse): the
%     across-stim peri-stim average atp trace must clear peak_factor
%     baseline-SDs above its own pre-stim baseline (same detector shape as
%     get_trial_stim_events above, just applied to the trial's own
%     stim-averaged trace instead of per stim)
%  3) per trial, detect walking bouts (smoothed |r_speed| > turn_thresh,
%     small gaps closed / short bouts dropped) and each bout's bump path
%     length (sum|diff(smoothed mu)|) vs. heading path length (integral of
%     |r_speed|)
%  4) per fly, DARK trials only: pool every walking bout from all trials
%     strictly BEFORE vs. strictly AFTER that fly's own first detected
%     perturbation trial (excluding any other detected perturbation
%     trials), then fit ONE weighted linear regression (bump path length
%     ~ heading path length, weighted by sqrt(bout duration)) per pool --
%     this is mov_ratio, "how far the bump moves per unit the fly actually
%     turned." A pool needs >=walk_min_bouts pooled bouts to get a value
%     (else NaN, dropped from the plot).
% !!! IMPORTANT CAVEAT -- RHO AND MOBILITY ARE COUPLED (found 2026-10-14) !!!
% In lpsp_p2x2_walking_20260728_post0708, bump vector strength (im.rho) and
% the bump-mobility slope move together: across all dark trials with >=10
% bouts, Spearman r(trial-mean rho, mobility) = -0.40 (n=73, p<0.001), and
% within just the lpsp post-ATP dark trials r = -0.59 (n=11, p=0.06). The
% lpsp group also shows a post-ATP DROP in dark-trial rho (fly-level mean
% 0.437 -> 0.382, -12%, paired signrank p=0.09, n=14 flies) that the empty
% group does not (0.454 -> 0.478, p=0.77, n=11). So the post-ATP mobility
% increase comes with a weaker/broader bump. This is NOT necessarily an
% artifact -- a destabilized EPG representation would look exactly like
% this -- and the lpsp post-dark trials sit above pre-ATP mobility even at
% normal rho (median 2.5 for rho>=0.35 trials vs 1.2 for all pre-dark).
% But the two measures are not independent: path length is a total-
% variation metric that inflates mechanically when mu is noisy, and mu
% (one 60-sample Gaussian in walk_trial_bouts) is smoothed less than
% r_speed (30 in process_ft + 60 here). REPORT the rho change alongside
% any mobility result, and consider a rho-gated or rho-matched version as
% a control. Quantified in the Claude session that produced
% ugly_figures/scripts/fig_example_flies_pb_heatmaps.m.
walk_source_file = 'lpsp_p2x2_walking_20260728_post0708.mat'; % already-cached subset (on/after 7/8), per that script's own section 1
tmp = load(fullfile(data_dir,walk_source_file),'all_data');
walk_data = tmp.all_data(:);
walk_n_trials = numel(walk_data);
fprintf('\n\n===== Figure_3 row 6 (bump mobility, lpsp_p2x2_walking) =====\nloaded %d trials from %s\n', walk_n_trials, walk_source_file);

walk_is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), walk_data);
walk_is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), walk_data);
walk_group   = repmat({'empty>p2x2'},walk_n_trials,1);
walk_group(walk_is_lpsp) = {'lpsp>p2x2'};

walk_fly_id    = cell(walk_n_trials,1);
walk_trial_num = nan(walk_n_trials,1);
for i = 1:walk_n_trials
    [walk_fly_id{i}, walk_trial_num(i)] = walk_meta_fly_trial(walk_data(i).meta);
end

walk_peri_win = 5; walk_base_win = [-walk_peri_win,0]; walk_peak_win = [0,3]; walk_peak_factor = 3;
walk_has_pulse = false(walk_n_trials,1);
for i = 1:walk_n_trials
    walk_has_pulse(i) = walk_trial_has_pulse(walk_data(i), walk_peri_win, walk_base_win, walk_peak_win, walk_peak_factor);
end
fprintf('%d/%d trials flagged as a detected perturbation trial\n', sum(walk_has_pulse), walk_n_trials);

walk_smooth_window = 60; walk_turn_thresh = .25; walk_max_gap_frames = 30; walk_min_walking_frames = 30; % same values as lpsp_p2x2_walking_script.m section 2 (.5*60 samples at xf's ~60Hz)
walk_bout_mu  = cell(walk_n_trials,1);
walk_bout_cue = cell(walk_n_trials,1);
walk_bout_dur = cell(walk_n_trials,1);
walk_bout_rho = cell(walk_n_trials,1);
for i = 1:walk_n_trials
    [walk_bout_mu{i}, walk_bout_cue{i}, walk_bout_dur{i}, walk_bout_rho{i}] = walk_trial_bouts(walk_data(i), walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames);
end

[walk_fly_list,~,walk_fly_ix] = unique(walk_fly_id);
walk_n_flies = numel(walk_fly_list);

walk_min_bouts  = 15; % minimum pooled bout count per fly x pre/post cell to trust its mov_ratio estimate, same as lpsp_p2x2_walking_script.m section 4
walk_fly_group  = nan(walk_n_flies,1); % 1 = empty, 2 = lpsp
walk_fly_pre    = nan(walk_n_flies,1); % dark only
walk_fly_post   = nan(walk_n_flies,1);
% mean bump vector strength (im.rho) across the same pooled walking bouts,
% duration-weighted (= mean over every bout sample, pooled) -- the rho
% control the RHO AND MOBILITY ARE COUPLED caveat above asks to report
% alongside mobility. Same walk_min_bouts gate, so the same flies.
walk_fly_rho_pre  = nan(walk_n_flies,1);
walk_fly_rho_post = nan(walk_n_flies,1);

for f = 1:walk_n_flies
    trial_idx = find(walk_fly_ix==f);
    [~,order] = sort(walk_trial_num(trial_idx));
    trial_idx = trial_idx(order); % chronological order of this fly's trials

    walk_fly_group(f) = 1 + strcmp(walk_group{trial_idx(1)},'lpsp>p2x2');

    pulse_pos = find(walk_has_pulse(trial_idx));
    if isempty(pulse_pos); continue; end % no detected perturbation for this fly

    pre_idx  = trial_idx(1:pulse_pos(1)-1);
    post_idx = trial_idx(pulse_pos(1)+1:end);
    post_idx = post_idx(~walk_has_pulse(post_idx)); % exclude any other detected perturbation trials

    pre_idx  = pre_idx(walk_is_dark(pre_idx));
    post_idx = post_idx(walk_is_dark(post_idx));
    if isempty(pre_idx) || isempty(post_idx); continue; end

    pre_cue = vertcat(walk_bout_cue{pre_idx}); pre_mu = vertcat(walk_bout_mu{pre_idx}); pre_dur = vertcat(walk_bout_dur{pre_idx});
    if numel(pre_cue) >= walk_min_bouts
        w = sqrt(pre_dur);
        walk_fly_pre(f) = (w.*pre_cue) \ (w.*pre_mu);
        pre_rho = vertcat(walk_bout_rho{pre_idx});
        walk_fly_rho_pre(f) = sum(pre_rho.*pre_dur,'omitnan') / sum(pre_dur(~isnan(pre_rho)));
    end

    post_cue = vertcat(walk_bout_cue{post_idx}); post_mu = vertcat(walk_bout_mu{post_idx}); post_dur = vertcat(walk_bout_dur{post_idx});
    if numel(post_cue) >= walk_min_bouts
        w = sqrt(post_dur);
        walk_fly_post(f) = (w.*post_cue) \ (w.*post_mu);
        post_rho = vertcat(walk_bout_rho{post_idx});
        walk_fly_rho_post(f) = sum(post_rho.*post_dur,'omitnan') / sum(post_dur(~isnan(post_rho)));
    end
end

walk_geno_order = {'empty>p2x2','lpsp>p2x2'};
for g = 1:2
    walk_valid_g = walk_fly_group==g & ~isnan(walk_fly_pre) & ~isnan(walk_fly_post);
    fprintf('row 6 (%s, dark, pre/post bump mobility): n=%d flies\n', walk_geno_order{g}, sum(walk_valid_g));
end

%% mobility example row: one lpsp fly, dark trials pre- and post-ATP
% Picked from ugly_figures/scripts/fig_mobility_example_candidates.m (per
% explicit request, its top row): lpsp 20260708 fly 1 -- dark mobility
% 0.88 -> 2.20 (Figure 3's own pooled slopes, nearer the group than the
% +2.77 outlier 20260727 fly 4), bout rho stable (0.45 -> 0.42) so the
% change isn't the rho-coupling artifact. (The closed-loop panel was
% dropped per explicit request -- dark trials only.)
% Windows (win_s each), as scored there:
%   dark pre  : trial 4, best bump-vs-heading velocity correlation among windows where the fly turns >= 2pi (r=.53, 9 bouts)
%   dark post : trial 8, window path ratio nearest the trial slope among windows with >= 2pi of turning (2.02 vs 2.20, 5 bouts)
% Each trace panel sits next to that WHOLE TRIAL's bout scatter (bump vs
% heading path length per walking bout, walk_bout_mu/cue from row 6) with
% its sqrt(dur)-weighted through-origin slope -- the per-trial version of
% the pooled pre/post slopes in the panels below.
% CAVEAT (see that script): in the dark the bump does NOT follow the fly's
% large turns in any candidate pre-ATP -- pre-ATP dark mobility ~1 comes
% from bump jitter accumulating path length, so the pair reads as "small
% jitter -> large bump excursions while the fly barely turns", not "tracks
% heading -> stops tracking". Heading in the dark is cumsum(r_speed)*dt
% (sign/scale checked against -cue in closed loop: -cue = 0.73 x it, the
% closed-loop gain); a slow r_speed bias would show up as a heading ramp.
mob_ex_tag    = '20260708/fly 1';
mob_ex_trials = [4 8];     % dark pre, dark post (trial numbers in the folder name)
mob_ex_t0     = [116 179]; % window start (s, ft.xf clock)
mob_ex_win_s  = 40;
mob_ex_label  = {'dark, pre-ATP','dark, post-ATP'};
mob_ex_f = find(contains(walk_fly_list, mob_ex_tag)); assert(isscalar(mob_ex_f), 'fly %s not found', mob_ex_tag);
mob_ex = struct('i',{},'xb',{},'z',{},'alpha',{},'t',{},'head',{},'mu',{},'bout_t',{},'b_head',{},'b_bump',{},'b_dur',{},'slope',{});
for c = 1:2
    i = find(walk_fly_ix==mob_ex_f & walk_trial_num==mob_ex_trials(c)); assert(isscalar(i));
    assert(walk_is_dark(i), 'trial %d is not a dark trial', mob_ex_trials(c));
    tr = walk_data(i); xf = tr.ft.xf(:); dt = median(diff(xf)); xb = tr.ft.xb(:);
    t0 = mob_ex_t0(c); t1 = t0 + mob_ex_win_s;
    head = cumsum(tr.ft.r_speed(:))*dt;
    mu = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    kf = xf>=t0 & xf<=t1; z0 = xf>=t0 & xf<=t0+1; % re-zero both at the window's first second
    iw = walk_bout_mask(tr, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames);
    d = diff([false; iw(kf); false]); tk = xf(kf)-t0;
    kb = xb>=t0 & xb<=t1;
    mob_ex(c).i = i;
    mob_ex(c).xb = xb(kb)-t0; mob_ex(c).z = tr.im.z(:,kb); mob_ex(c).alpha = unwrap(tr.im.alpha);
    mob_ex(c).t = tk; mob_ex(c).head = head(kf) - mean(head(z0)); mob_ex(c).mu = mu(kf) - mean(mu(z0));
    mob_ex(c).bout_t = [tk(d(1:end-1)==1), tk(find(d==-1)-1)]; % [start end] per walking bout
    mob_ex(c).b_head = walk_bout_cue{i}; mob_ex(c).b_bump = walk_bout_mu{i}; mob_ex(c).b_dur = walk_bout_dur{i};
    w = sqrt(mob_ex(c).b_dur);
    mob_ex(c).slope = (w.*mob_ex(c).b_head) \ (w.*mob_ex(c).b_bump);
    fprintf('mobility example %s trial %d: %d bouts, trial slope %.2f\n', mob_ex_tag, mob_ex_trials(c), numel(w), mob_ex(c).slope);
end

%% plot: 2 example columns (opposite stim sides) for rows 1-3, then rows
% 4-5 (population, genotype-level) also split into 2 columns each, then
% row 6 (bump mobility, genotype-level, dark trials only)
% built invisible and only shown after export: a visible figure gets
% clamped to the screen height (e.g. 833 px on a laptop display), which
% squashes the export; an invisible one keeps the full requested size.
fig3 = figure('color','w','Visible','off','Position',[50 50 1700 1700]); clf
% underlying grid is 4 rows x 24 columns so that the bottom row can hold
% 4 equal-width panels (colspan 6 each), row 3 (the mobility example)
% trace/scatter pairs (colspan 7 + 5 -- the square scatters are width-
% limited, so they get the extra columns), and rows 1-2 two (colspan 12
% each), all inside one tiledlayout. Row 1 is a nested 5x1 layout per example column: a short
% whole-PB atp dF/F trace (scale bar instead of an axis, per explicit
% request -- replaces the old full-size row 3) sitting directly on top of
% the EPG calcium heatmap (4/5 of the height).
t = tiledlayout(4,24,'TileSpacing','loose','Padding','compact');

row1_axes = gobjects(1,2);
trace_axes = gobjects(1,2);
for c = 1:2
    d = example_disp(c);
    i = example_idx(c);
    zoom_idx = d.xb>=d.zoom_t0 & d.xb<=d.zoom_t1;

    % outer 1x6 layout per example column: the trial's mean image (col 1)
    % to the left of the trace+heatmap stack (cols 2-6). The image sits in
    % its own 5x1 layout, in the same bottom 4/5 as the heatmap, so the
    % rotated PB lines up vertically with the heatmap rows.
    t1 = tiledlayout(t,1,6,'TileSpacing','compact','Padding','none');
    t1.Layout.Tile = (c-1)*12+1; t1.Layout.TileSpan = [1 12];

    tm = tiledlayout(t1,5,1,'TileSpacing','none','Padding','none');
    tm.Layout.Tile = 1;
    axm = nexttile(tm,2,[4 1]);
    image(axm,d.mean_rgb); axis(axm,'image','off')
    text(axm,1,size(d.mean_rgb,1),' EPG','Color',gcamp_color,'FontWeight','bold','FontSize',8,'VerticalAlignment','bottom')
    text(axm,size(d.mean_rgb,2),size(d.mean_rgb,1),'ATP ','Color',atp_color,'FontWeight','bold','FontSize',8,'VerticalAlignment','bottom','HorizontalAlignment','right')

    t1b = tiledlayout(t1,5,1,'TileSpacing','none','Padding','none');
    t1b.Layout.Tile = 2; t1b.Layout.TileSpan = [1 5];

    % small atp trace on top: no axis, a vertical dF/F scale bar on the right
    axt = nexttile(t1b,1); hold(axt,'on')
    plot(axt,d.xb,d.atp_avg_full,'Color',atp_color,'LineWidth',1.2)
    trace_lim = [min(d.atp_avg_full(zoom_idx)), max(d.atp_avg_full(zoom_idx))];
    trace_lim = trace_lim + [-.05,.05]*diff(trace_lim);
    plot(axt,[d.stim_t,d.stim_t],trace_lim,':k')
    xlim(axt,[d.zoom_t0,d.zoom_t1]); ylim(axt,trace_lim)
    bar_len = nice_scale_len(diff(trace_lim)/2);
    bar_x = d.zoom_t1 + .01*(d.zoom_t1-d.zoom_t0);
    plot(axt,[bar_x,bar_x],trace_lim(1)+[0,bar_len],'k','LineWidth',2,'Clipping','off')
    text(axt,bar_x,trace_lim(1)+bar_len/2,sprintf('  %g dF/F',bar_len),'VerticalAlignment','middle','FontSize',9,'Clipping','off')
    axt.Visible = 'off';
    cb_pad = colorbar(axt); cb_pad.Visible = 'off'; % reserves the same right-hand gutter as the heatmap's colorbar so the two plot boxes line up
    trace_axes(c) = axt;

    axg = nexttile(t1b,2,[4 1]); hold(axg,'on')
    imagesc(axg,d.xb,d.im_alpha,d.ex.im.z,d.im_clim)
    colormap(axg,white_to_color(gcamp_color))
    cbg = colorbar(axg); cbg.Label.String = 'z-score';
    plot(axg,[d.stim_t,d.stim_t],[min(d.im_alpha),max(d.im_alpha)],':k','LineWidth',1)
    xlim(axg,[d.zoom_t0,d.zoom_t1]); ylim(axg,[min(d.im_alpha),max(d.im_alpha)])
    if c==1; ylabel(axg,'PB angle (rad)'); end
    title(t1,sprintf('%s-stim example (%s): whole-PB atp dF/F (top), EPG calcium im.z (bottom)\nfly %s, trial "%s", stim %d', example_col_label{c}, example_geno_label{c}, fly_id{i}, meta_display(d.ex.meta), example_stim_idx(c)),'Interpreter','none','FontSize',10,'FontWeight','bold')
    row1_axes(c) = axg;
end

% (former row 2, the atp-channel heatmaps under each example, moved to
% Figure 3 supplement S5 per explicit request -- see Figure_3_supplements.m)

% row 2 panels sit in the same 1x6 split as row 1 (col 1 left empty under
% the mean image, cols 2-6 for the plot) plus an invisible colorbar, per
% explicit request, so their plot boxes line up with the heatmaps above
row2_axes = gobjects(1,2);
for c = 1:2
    t2 = tiledlayout(t,1,6,'TileSpacing','compact','Padding','none');
    t2.Layout.Tile = 24+(c-1)*12+1; t2.Layout.TileSpan = [1 12];
    row2_axes(c) = nexttile(t2,2,[1 5]);
    cb_pad = colorbar(row2_axes(c)); cb_pad.Visible = 'off';
end
ax4a = row2_axes(1); hold(ax4a,'on')
h_geno = gobjects(1,numel(geno_order));
for gi = 1:numel(geno_order)
    patch(ax4a,[t_trace,fliplr(t_trace)],[geno_mean(gi,:)+geno_sem(gi,:),fliplr(geno_mean(gi,:)-geno_sem(gi,:))], ...
        geno_color(gi,:),'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off')
    h_geno(gi) = plot(ax4a,t_trace,geno_mean(gi,:),'Color',geno_color(gi,:),'LineWidth',2);
end
xline(ax4a,0,':k','HandleVisibility','off') % xline, not plot(ylim): spans whatever y-limits linkaxes settles on without driving them
plot(ax4a,[trace_pre_s,trace_post_s],[0,0],':k','HandleVisibility','off')
xlim(ax4a,[trace_pre_s,trace_post_s])
xlabel(ax4a,'time from stim onset (s)')
ylabel(ax4a,sprintf('bump pos. (rad, %s-flip)',ref_side_str))
legend(ax4a,h_geno,{sprintf('empty>P2X2 (n=%d flies)',geno_n(1)),sprintf('lpsp>P2X2 (n=%d flies)',geno_n(2))},'Location','best')
title(ax4a,sprintf('average bump-position response (flipped to match the %s-stim example''s side)',ref_side_str))

ax3b = row2_axes(2); hold(ax3b,'on') % row 3 col 2: lpsp bump position by stim number (replaces the whole-PB atp population panel)
stimnum_cmap = (1-linspace(.75,0,n_stim_per_trial)').*gcamp_color + linspace(.75,0,n_stim_per_trial)'.*[1 1 1]; % light -> dark with stim number, same ramp as Figure_3_supplements.m
h_stim = gobjects(1,n_stim_per_trial);
for sn = 1:n_stim_per_trial
    patch(ax3b,[t_trace,fliplr(t_trace)],[stimnum_mean(sn,:)+stimnum_sem(sn,:),fliplr(stimnum_mean(sn,:)-stimnum_sem(sn,:))], ...
        stimnum_cmap(sn,:),'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off')
    h_stim(sn) = plot(ax3b,t_trace,stimnum_mean(sn,:),'Color',stimnum_cmap(sn,:),'LineWidth',2);
end
xline(ax3b,0,':k','HandleVisibility','off')
plot(ax3b,[trace_pre_s,trace_post_s],[0,0],':k','HandleVisibility','off')
xlim(ax3b,[trace_pre_s,trace_post_s])
xlabel(ax3b,'time from stim onset (s)')
ylabel(ax3b,sprintf('bump pos. (rad, %s-flip, baseline-sub.)',ref_side_str))
legend(ax3b,h_stim,arrayfun(@(sn) sprintf('stim %d in trial (n=%d flies)',sn,stimnum_n(sn)),1:n_stim_per_trial,'UniformOutput',false),'Location','northwest')
title(ax3b,'lpsp>P2X2 bump-position response by stim number within trial')

% row 3: mobility example (one lpsp fly, dark pre | dark post), each a
% trace window (span 7: nested EPG heatmap on top, 2/5, traces below, 3/5)
% + that trial's bout scatter (span 5)
mob_trace_axes = gobjects(1,2); mob_scatter_axes = gobjects(1,2); mob_heat_axes = gobjects(1,2);
bout_shade = [.8 .8 .8];
for c = 1:2
    e = mob_ex(c);
    t3 = tiledlayout(t,5,1,'TileSpacing','none','Padding','none');
    t3.Layout.Tile = 48+(c-1)*12+1; t3.Layout.TileSpan = [1 7];
    axh = nexttile(t3,1,[2 1]); hold(axh,'on')
    imagesc(axh, e.xb, e.alpha, e.z, [prctile(min(e.z,[],1),20), prctile(max(e.z,[],1),95)]) % floor at the 20th (not 5th) percentile of per-frame troughs, per explicit request -- whitens the background so the bump stands out
    colormap(axh, white_to_color(gcamp_color))
    xlim(axh,[0 mob_ex_win_s]); ylim(axh,[min(e.alpha) max(e.alpha)]); xticks(axh,[])
    if c==1; ylabel(axh,'PB angle (rad)'); end
    mob_heat_axes(c) = axh;
    axt = nexttile(t3,3,[3 1]); hold(axt,'on')
    hb = gobjects(0);
    for b = 1:size(e.bout_t,1)
        hb(end+1) = patch(axt, e.bout_t(b,[1 2 2 1]), [-100 -100 100 100], bout_shade, 'EdgeColor','none'); %#ok<AGROW>
    end
    h1 = plot(axt, e.t, e.head, 'k-', 'LineWidth', 1.2);
    h2 = plot(axt, e.t, e.mu, '-', 'Color', gcamp_color, 'LineWidth', 1.2);
    xlim(axt,[0 mob_ex_win_s])
    if c==1; ylabel(axt,'angle (rad, re-zeroed)'); end
    lg = legend(axt,[h1 h2 hb(1)],{'heading (\int yaw velocity)','bump (im.mu)','walking bout'},'Location','northwest','FontSize',8);
    lg.Color = 'w'; lg.EdgeColor = 'w'; % opaque: over a grey bout the 'walking bout' swatch would vanish
    if c==1; fly_str = sprintf('lpsp>P2X2 example fly %s: ', mob_ex_tag); else; fly_str = ''; end
    title(t3, sprintf('%s%s (trial %d)', fly_str, mob_ex_label{c}, mob_ex_trials(c)), 'FontSize', 10, 'FontWeight', 'bold')
    mob_trace_axes(c) = axt;

    axs = nexttile(t,48+(c-1)*12+8,[1 5]); hold(axs,'on')
    scatter(axs, e.b_head, e.b_bump, 12+4*sqrt(e.b_dur), [.4 .4 .4], 'filled', 'MarkerFaceAlpha', .5)
    mob_scatter_axes(c) = axs;
end
% shared scales: both dark trace panels on one angle scale (explicit -- the
% +/-100 bout patches must not drive autoscaling), both scatters on one
% path-length scale so the two slopes compare directly
mob_ylim = max(cellfun(@(v) max(abs(v)), [{mob_ex.head},{mob_ex.mu}]))*1.1;
for c = 1:2; ylim(mob_trace_axes(c),[-mob_ylim mob_ylim]); end
linkaxes(mob_trace_axes,'xy'); linkaxes([mob_heat_axes, mob_trace_axes],'x'); set_pi_yticks(mob_trace_axes(1)); set_pi_yticks(mob_trace_axes(2))
sc_max = max([vertcat(mob_ex.b_head); vertcat(mob_ex.b_bump)])*1.05;
for c = 1:2
    axs = mob_scatter_axes(c); e = mob_ex(c);
    plot(axs, [0 sc_max], [0 sc_max], 'k:')
    plot(axs, [0 sc_max], e.slope*[0 sc_max], '-', 'Color', gcamp_color, 'LineWidth', 1.5)
    xlim(axs,[0 sc_max]); ylim(axs,[0 sc_max]); axis(axs,'square')
    xlabel(axs,'heading path / bout (rad)'); ylabel(axs,'bump path / bout (rad)')
    text(axs, .95*sc_max, .05*sc_max, sprintf('mobility = %.2f\n%d bouts', e.slope, numel(e.b_dur)), 'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',9,'Color',gcamp_color)
    title(axs, sprintf('trial %d, one dot per bout', mob_ex_trials(c)), 'FontSize', 10)
end

ax6a = nexttile(t,73,[1 6]); hold(ax6a,'on')
walk_valid = walk_fly_group==1 & ~isnan(walk_fly_pre) & ~isnan(walk_fly_post);
walk_pre_vals  = walk_fly_pre(walk_valid);
walk_post_vals = walk_fly_post(walk_valid);
plot(ax6a,[ones(sum(walk_valid),1),2*ones(sum(walk_valid),1)]',[walk_pre_vals,walk_post_vals]','Color',[.6,.6,.6,.5])
scatter(ax6a,ones(sum(walk_valid),1), walk_pre_vals, 'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
scatter(ax6a,2*ones(sum(walk_valid),1),walk_post_vals,'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
errorbar(ax6a,[1,2],[mean(walk_pre_vals,'omitnan'),mean(walk_post_vals,'omitnan')], ...
    [std(walk_pre_vals,'omitnan'),std(walk_post_vals,'omitnan')]/sqrt(sum(walk_valid)),'-ok','LineWidth',2)
xlim(ax6a,[.5,2.5]); xticks(ax6a,[1,2]); xticklabels(ax6a,{'pre','post'})
ylabel(ax6a,'bump path length / heading path length')
title(ax6a,sprintf('empty>P2X2, dark trials (n=%d flies)',sum(walk_valid)))

ax6b = nexttile(t,79,[1 6]); hold(ax6b,'on')
walk_valid = walk_fly_group==2 & ~isnan(walk_fly_pre) & ~isnan(walk_fly_post);
walk_pre_vals  = walk_fly_pre(walk_valid);
walk_post_vals = walk_fly_post(walk_valid);
plot(ax6b,[ones(sum(walk_valid),1),2*ones(sum(walk_valid),1)]',[walk_pre_vals,walk_post_vals]','Color',[.6,.6,.6,.5])
scatter(ax6b,ones(sum(walk_valid),1), walk_pre_vals, 'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
scatter(ax6b,2*ones(sum(walk_valid),1),walk_post_vals,'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
errorbar(ax6b,[1,2],[mean(walk_pre_vals,'omitnan'),mean(walk_post_vals,'omitnan')], ...
    [std(walk_pre_vals,'omitnan'),std(walk_post_vals,'omitnan')]/sqrt(sum(walk_valid)),'-ok','LineWidth',2)
xlim(ax6b,[.5,2.5]); xticks(ax6b,[1,2]); xticklabels(ax6b,{'pre','post'})
ylabel(ax6b,'bump path length / heading path length')
title(ax6b,sprintf('lpsp>P2X2, dark trials (n=%d flies)',sum(walk_valid)))

ax6c = nexttile(t,85,[1 6]); hold(ax6c,'on')
% delta in bump mobility (post - pre), dark trials, by genotype -- ported
% from lpsp_p2x2_walking_script.m section 6 (the dark-only column of its
% own "change in bump mobility after perturbation" figure), reusing
% walk_fly_pre/walk_fly_post/walk_fly_group computed above.
walk_delta = walk_fly_post - walk_fly_pre;
walk_delta_valid = ~isnan(walk_delta);
scatter(ax6c,walk_fly_group(walk_delta_valid),walk_delta(walk_delta_valid),'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
for g = 1:2
    gv = walk_delta_valid & walk_fly_group==g;
    errorbar(ax6c,g+.15,mean(walk_delta(gv),'omitnan'),std(walk_delta(gv),'omitnan')/sqrt(sum(gv)),'ok','LineWidth',2,'MarkerFaceColor','k')
end
plot(ax6c,[.5,2.5],[0,0],':k','HandleVisibility','off')
xlim(ax6c,[.5,2.5]); xticks(ax6c,[1,2]); xticklabels(ax6c,{'empty>P2X2','lpsp>P2X2'})
ylabel(ax6c,'\Delta bump mobility (post - pre)')
title(ax6c,sprintf('change in bump mobility after perturbation\ndark trials (n=%d empty, n=%d lpsp)', ...
    sum(walk_delta_valid&walk_fly_group==1), sum(walk_delta_valid&walk_fly_group==2)))
add_delta_pvals(ax6c, walk_delta(walk_delta_valid), walk_fly_group(walk_delta_valid))

ax6d = nexttile(t,91,[1 6]); hold(ax6d,'on')
% delta in mean bump vector strength (im.rho) across walking bouts (post -
% pre), dark trials, by genotype -- same flies, bouts and pre/post pools
% as ax6c (walk_fly_rho_pre/post above), so the two panels read as a pair:
% the rho control the RHO AND MOBILITY ARE COUPLED caveat asks for.
walk_rho_delta = walk_fly_rho_post - walk_fly_rho_pre;
walk_rho_delta_valid = walk_delta_valid & ~isnan(walk_rho_delta);
scatter(ax6d,walk_fly_group(walk_rho_delta_valid),walk_rho_delta(walk_rho_delta_valid),'filled','MarkerFaceColor',[.6,.6,.6],'MarkerFaceAlpha',.6)
for g = 1:2
    gv = walk_rho_delta_valid & walk_fly_group==g;
    errorbar(ax6d,g+.15,mean(walk_rho_delta(gv),'omitnan'),std(walk_rho_delta(gv),'omitnan')/sqrt(sum(gv)),'ok','LineWidth',2,'MarkerFaceColor','k')
    fprintf('row 6 rho delta (%s, dark): mean %+.3f, n=%d flies, signrank p=%.2f\n', walk_geno_order{g}, mean(walk_rho_delta(gv)), sum(gv), signrank(walk_rho_delta(gv)));
end
plot(ax6d,[.5,2.5],[0,0],':k','HandleVisibility','off')
xlim(ax6d,[.5,2.5]); xticks(ax6d,[1,2]); xticklabels(ax6d,{'empty>P2X2','lpsp>P2X2'})
ylabel(ax6d,'\Delta mean bout \rho (post - pre)')
title(ax6d,sprintf('change in mean bump \\rho across walking bouts\ndark trials (n=%d empty, n=%d lpsp)', ...
    sum(walk_rho_delta_valid&walk_fly_group==1), sum(walk_rho_delta_valid&walk_fly_group==2)))
add_delta_pvals(ax6d, walk_rho_delta(walk_rho_delta_valid), walk_fly_group(walk_rho_delta_valid))

linkaxes([trace_axes(1),row1_axes(1)],'x')
linkaxes([trace_axes(2),row1_axes(2)],'x')
for c = 1:2; add_time_scalebar(row1_axes(c), 10); add_time_scalebar(mob_trace_axes(c), 10); end % per explicit request: 10 s scale bar instead of a labeled time axis (rows 1 and 3)
linkaxes([ax4a,ax3b],'xy') % row 3: both bump position (rad) -- same scale so the by-stim lines read directly against the genotype average
set_pi_yticks(ax4a); set_pi_yticks(ax3b) % after linkaxes, so ticks match the shared limits
linkaxes([ax6a,ax6b],'y') % separate from rows 2-3 -- totally different y-scale (mobility ratio, not bump position/dF/F); ax6c (a delta, can be negative) is its own scale too
% bigger axis text, per explicit request: tick labels 12 pt, axis labels
% 14 pt, colorbars 11/13 pt; titles kept at their own sizes (setting an
% axes' FontSize would otherwise rescale them too). 16 pt labels were
% tried first -- the long row-4 y-labels overflowed their axes.
for ax_k = findall(fig3,'Type','axes')'
    if strcmp(ax_k.Visible,'off'); continue; end % the scale-bar-only trace axes
    title_fs = ax_k.Title.FontSize;
    ax_k.FontSize = 12;
    ax_k.Title.FontSize = title_fs;
    ax_k.XLabel.FontSize = 14;
    ax_k.YLabel.FontSize = 14;
end
for cb_k = findall(fig3,'Type','colorbar')'
    cb_k.FontSize = 11;
    cb_k.Label.FontSize = 13;
end

sgtitle('Figure 3 (draft). Supplements S1-S5: Figure_3_supplements.m','FontSize',10,'Interpreter','none')

%% export
if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig3, fullfile(export_dir,'Figure_3_claude.png'), 'Resolution', 300)
if ~isfolder(all_figs_dir); mkdir(all_figs_dir); end
exportgraphics(fig3, fullfile(all_figs_dir,'Fig3_V1.pdf'), 'ContentType', 'vector')
fig3.Visible = 'on';

%% ===================== SCRATCH: verbose walkthrough of the stim-side flip =====================
% exported separately from the main figure, purely to show step by step
% how the left/right call is made (now: most-confident-stim-decides-the-
% whole-trial, per explicit request) and what flipping does to real
% traces. Not part of the polished Figure 3 panels above.
fig_scratch = figure('color','w','Position',[50 50 1500 1000]); clf
ts = tiledlayout(fig_scratch,2,2,'TileSpacing','compact','Padding','compact');
side_color = [0.00,0.45,0.74; 0.85,0,0]; % RIGHT/LEFT diagnostic colors -- independent of geno_color (genotype), which these panels don't encode

% panel A: EVERY confirmed-ejection stim's own right_frac (not one value
% per trial), grouped by trial and sorted by that trial's deciding value --
% shows directly that stims within the same trial mostly agree with each
% other (supporting the assumption a single side applies to the whole
% trial), and marks with a black outline which stim was actually used to
% decide that trial's call.
axA = nexttile(ts); hold(axA,'on')
called_trials = find(have_call);
[~,ord] = sort(trial_right_frac(called_trials));
called_trials = called_trials(ord);
patch(axA,[0.4,0.6,0.6,0.4],[0,0,numel(called_trials)+1,numel(called_trials)+1],[.85,.85,.85],'FaceAlpha',.6,'EdgeColor','none')
for r = 1:numel(called_trials)
    i = called_trials(r);
    ev = trial_events{i};
    rfs = [ev.right_frac];
    if trial_is_right(i); c = side_color(1,:); else; c = side_color(2,:); end
    plot(axA,rfs,r*ones(size(rfs)),'-','Color',[c,.35])
    scatter(axA,rfs,r*ones(size(rfs)),16,c,'filled','MarkerFaceAlpha',.6)
    [~,mc] = max(abs(rfs-0.5));
    scatter(axA,rfs(mc),r,32,c,'filled','MarkerEdgeColor','k','LineWidth',1)
end
xline(axA,0.5,':k','LineWidth',1.5)
xlabel(axA,'right\_frac per individual confirmed-ejection stim (peak-window atp.d spatial profile)')
ylabel(axA,'trial rank (sorted by the DECIDING stim''s right\_frac)')
title(axA,sprintf('step 1: within-trial agreement (n=%d trials); outlined dot = deciding stim',numel(called_trials)))

% panel B: the raw peak-window spatial profile behind two example calls --
% the trial whose deciding stim was most confidently RIGHT, and the one
% most confidently LEFT.
[~,ex_right_i] = max(trial_right_frac);
[~,ex_left_i]  = min(trial_right_frac);
ev_r = trial_events{ex_right_i}; rfs_r = [ev_r.right_frac]; [~,mc_r] = max(abs(rfs_r-0.5));
ev_l = trial_events{ex_left_i};  rfs_l = [ev_l.right_frac]; [~,mc_l] = max(abs(rfs_l-0.5));
prof_right = ev_r(mc_r).wedge_peak; prof_right = prof_right/max(prof_right);
prof_left  = ev_l(mc_l).wedge_peak; prof_left  = prof_left/max(prof_left);

axB = nexttile(ts); hold(axB,'on')
plot(axB,1:32,prof_right,'-o','Color',side_color(1,:),'MarkerFaceColor',side_color(1,:))
plot(axB,1:32,prof_left, '-o','Color',side_color(2,:),'MarkerFaceColor',side_color(2,:))
xline(axB,16.5,':k','LineWidth',1.5)
legend(axB,{sprintf('example RIGHT trial''s deciding stim (right\\_frac=%.2f)',rfs_r(mc_r)), ...
            sprintf('example LEFT trial''s deciding stim (right\\_frac=%.2f)',rfs_l(mc_l))},'Location','south')
xlabel(axB,'atp wedge index (1-16 = first half, 17-32 = second half)')
ylabel(axB,'that stim''s own 0-3s peak-window atp.d, normalized to its own peak')
title(axB,'step 2: the (per-stim, peak-window) spatial profile behind the deciding call')

% panels C/D: a handful of individual confirmed-ejection stim traces
% (bump position, zeroed at onset) from BOTH sides, BEFORE (panel C) and
% AFTER (panel D) flipping onto the EXAMPLE TRIAL'S OWN side
% (reference_is_right) -- same logic as the row-3 pooling above, just
% shown per-stim instead of averaged away, so the effect of flipping is
% visible directly.
n_example_each_side = 4;
example_traces_raw = {}; example_is_right = [];
for want_right = [true,false]
    found = 0;
    for i = 1:n_trials
        if isnan(trial_is_right(i)) || trial_is_right(i) ~= want_right; continue; end
        ev = trial_events{i};
        for s = 1:numel(ev)
            if mean(ev(s).rho_grid(trace_idx),'omitnan') < rho_gate; continue; end % keep this diagnostic to reasonably-tracked stims
            example_traces_raw{end+1} = ev(s).mu_rel(trace_idx); %#ok<AGROW>
            example_is_right(end+1)   = want_right; %#ok<AGROW>
            found = found + 1;
            if found >= n_example_each_side; break; end
        end
        if found >= n_example_each_side; break; end
    end
end

axC = nexttile(ts); hold(axC,'on')
axD = nexttile(ts); hold(axD,'on')
for e = 1:numel(example_traces_raw)
    if example_is_right(e); c = side_color(1,:); else; c = side_color(2,:); end
    plot(axC,t_trace,example_traces_raw{e},'Color',c,'LineWidth',1)
    flipped = example_traces_raw{e};
    if example_is_right(e) ~= reference_is_right
        flipped = -flipped;
    end
    plot(axD,t_trace,flipped,'Color',c,'LineWidth',1)
end
for ax = [axC,axD]
    plot(ax,[0,0],ylim(ax),':k')
    plot(ax,[trace_pre_s,trace_post_s],[0,0],':k')
    xlim(ax,[trace_pre_s,trace_post_s])
    xlabel(ax,'time from stim onset (s)')
    ylabel(ax,'bump position rel. to onset (rad)')
end
h_r = plot(axC,nan,nan,'Color',side_color(1,:),'LineWidth',1.5);
h_l = plot(axC,nan,nan,'Color',side_color(2,:),'LineWidth',1.5);
legend(axC,[h_r,h_l],{'RIGHT-stim trial (unflipped)','LEFT-stim trial (unflipped)'},'Location','best')
title(axC,sprintf('step 3: %d example confirmed-ejection stims, BEFORE flip',numel(example_traces_raw)))
if reference_is_right
    lbl_r = 'RIGHT-stim trial (unchanged, matches example)';
    lbl_l = 'LEFT-stim trial (negated)';
else
    lbl_r = 'RIGHT-stim trial (negated)';
    lbl_l = 'LEFT-stim trial (unchanged, matches example)';
end
legend(axD,[h_r,h_l],{lbl_r,lbl_l},'Location','best')
title(axD,sprintf('step 4: same stims, AFTER flipping onto the example trial''s own side (%s)',ref_side_str))

sgtitle(ts,sprintf('scratch: stim-side flip walkthrough (side decided per-trial by its most confident stim, then flipped to the example trial''s %s side; %d/%d disagree with the old whole-trial method)', ...
    ref_side_str, n_reclassified, sum(have_call)))

if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig_scratch, fullfile(export_dir,'Figure_3_claude_scratch_stim_side_flip.png'), 'Resolution', 300)

%% ===================== Figure_3_dopamine: same pipeline, dopamine_ionto_redo dataset =====================
% Adapts everything above (ejection detection, two-opposite-hemisphere
% example trials, population bump-position/atp/EPG traces) to
% dopamine_ionto_redo_20260518.mat (empty>P2X2 driving ATP-triggered
% depolarization of dopaminergic input near the PB, EPG-syt8s imaging) --
% integrated here as its own figure per explicit request, after review as
% a standalone scratch script. See data scripts/dopamine_ionto_redo_claude.m's
% own header for the full investigation this is ported from.
%
% Three real differences from the p2x2 pipeline above, each isolated with
% a dop_ prefix / its own function rather than touched in the p2x2 code:
% 1) NO GENOTYPE SPLIT: 20/20 trials carry "empty_p2x2_da" in their trial
%    folder name, 0 carry "lpsp" -- confirmed directly, single group. Every
%    "one line per genotype" panel above is instead ONE line (all flies).
% 2) NO LIGHT-CONDITION SPLIT: all 20 trials are ft.pattern ==
%    "0003_4px_brightbar.mat" (closed loop); no dark trials.
% 3) HEMISPHERE CALL NEEDS ITS OWN FUNCTION (get_trial_stim_events_dopamine):
%    the per-stim peak-window atp.d profile needs clipping to positive
%    values only before splitting per hemisphere half here -- confirmed
%    necessary on this dataset (unlike p2x2) since a raw dip elsewhere is
%    not evidence of "ATP didn't land here" and was otherwise free to swing
%    the balance. Also fly folders here are bare "fly N" (no genotype
%    suffix), but trial_fly_id's regex on "fly\s*\d+" matches that fine.
%
% peri_win/base_win/peak_win/peak_factor/t_grid/settle_win/rho_gate are
% REUSED UNCHANGED from the p2x2 section above: dopamine_ionto_redo_claude.m
% confirmed the same 4-stims-per-trial, ~60s-apart timing on this dataset,
% so there's no reason to re-tune them. dop_atp_peak_amp_thresh is its own
% value (0.5, not 2.0): this dataset's raw dF/F signal is much smaller
% (checked directly: min=0.02, p25=0.11, median=0.29, p75=0.60, max=1.12
% vs. the p2x2 dataset's min=0.09/median=2.64/max=10.94).

%% load dopamine_ionto_redo data
dop_source_file = 'dopamine_ionto_redo_20260518.mat';
tmp = load(fullfile(data_dir,dop_source_file),'all_data');
dop_all_data = tmp.all_data(:);
dop_n_trials = numel(dop_all_data);
fprintf('\n\n===== Figure_3_dopamine =====\nloaded %d trials from %s\n', dop_n_trials, dop_source_file);

dop_fly_id = cell(dop_n_trials,1);
for i = 1:dop_n_trials
    dop_fly_id{i} = trial_fly_id(dop_all_data(i).meta);
end
fprintf('%d trials -> %d flies\n', dop_n_trials, numel(unique(dop_fly_id)));

%% precompute confirmed-ejection stim events (positive-clipped hemisphere metric)
dop_trial_events = cell(dop_n_trials,1);
for i = 1:dop_n_trials
    dop_trial_events{i} = get_trial_stim_events_dopamine(dop_all_data(i), t_grid, peri_win, base_win, peak_win, peak_factor);
end
dop_n_ejections_confirmed_all = cellfun(@numel,dop_trial_events);
fprintf('%d/%d trials have >=1 confirmed ejection (%d confirmed ejections total)\n', ...
    sum(dop_n_ejections_confirmed_all>0), dop_n_trials, sum(dop_n_ejections_confirmed_all));

%% which PB hemisphere each trial's ATP pipette targeted
dop_trial_is_right   = nan(dop_n_trials,1);
dop_trial_right_frac = nan(dop_n_trials,1);
for i = 1:dop_n_trials
    ev = dop_trial_events{i};
    if isempty(ev); continue; end
    rfs = [ev.right_frac];
    [~,most_confident] = max(abs(rfs-0.5));
    dop_trial_right_frac(i) = rfs(most_confident);
    dop_trial_is_right(i)   = dop_trial_right_frac(i) > 0.5;
end
dop_have_call = ~isnan(dop_trial_is_right);
fprintf('stim-side call: %d/%d trials got a call; side split: %d RIGHT, %d LEFT\n', ...
    sum(dop_have_call), dop_n_trials, sum(dop_trial_is_right(dop_have_call)==1), sum(dop_trial_is_right(dop_have_call)==0));

%% pick two example STIMS -- one per side -- with good post-ejection bump
% movement AND a clean atp_peak_amp on that same stim (identical selection
% logic to the p2x2 section above, minus the genotype restriction).
dop_atp_peak_amp_thresh = 0.5;

dop_cand_trial = []; dop_cand_stim = []; dop_cand_score = []; dop_cand_side = [];
for i = 1:dop_n_trials
    if dop_n_ejections_confirmed_all(i) < 2; continue; end
    ev = dop_trial_events{i};
    for s = 1:numel(ev)
        if ev(s).atp_peak_amp < dop_atp_peak_amp_thresh; continue; end
        divergence = ev(s).mu_rel - ev(s).head_rel;
        settle_idx = t_grid>settle_win(1) & t_grid<=settle_win(2);
        if mean(ev(s).rho_grid(settle_idx),'omitnan') <= rho_gate; continue; end
        dop_cand_trial(end+1) = i; %#ok<AGROW>
        dop_cand_stim(end+1)  = s; %#ok<AGROW>
        dop_cand_score(end+1) = abs(median(divergence(settle_idx),'omitnan')); %#ok<AGROW>
        dop_cand_side(end+1)  = dop_trial_is_right(i); %#ok<AGROW>
    end
end
fprintf('%d individual confirmed-ejection stims clear atp_peak_amp>=%.1f and the rho gate (%d RIGHT, %d LEFT)\n', ...
    numel(dop_cand_trial), dop_atp_peak_amp_thresh, sum(dop_cand_side==1), sum(dop_cand_side==0));

[dop_best_score,dop_best] = max(dop_cand_score);
dop_example_i = dop_cand_trial(dop_best);
dop_example_s = dop_cand_stim(dop_best);
dop_reference_is_right = logical(dop_cand_side(dop_best));
if dop_reference_is_right; dop_ref_side_str = 'RIGHT'; else; dop_ref_side_str = 'LEFT'; end

fprintf('example trial: fly %s, trial "%s", stim %d\n', dop_fly_id{dop_example_i}, meta_display(dop_all_data(dop_example_i).meta), dop_example_s);
fprintf('  post-ejection bump displacement = %.2f rad, atp_peak_amp = %.2f, side = %s\n', dop_best_score, dop_trial_events{dop_example_i}(dop_example_s).atp_peak_amp, dop_ref_side_str);

dop_opp_idx = find(dop_cand_side == ~dop_reference_is_right);
[dop_best_score2,dop_best2rel] = max(dop_cand_score(dop_opp_idx));
dop_best2 = dop_opp_idx(dop_best2rel);
dop_example_i2 = dop_cand_trial(dop_best2);
dop_example_s2 = dop_cand_stim(dop_best2);
if dop_reference_is_right; dop_opp_side_str = 'LEFT'; else; dop_opp_side_str = 'RIGHT'; end
fprintf('second example trial (opposite side, %s): fly %s, trial "%s", stim %d\n', dop_opp_side_str, dop_fly_id{dop_example_i2}, meta_display(dop_all_data(dop_example_i2).meta), dop_example_s2);
fprintf('  post-ejection bump displacement = %.2f rad, atp_peak_amp = %.2f\n', dop_best_score2, dop_trial_events{dop_example_i2}(dop_example_s2).atp_peak_amp);

%% example display prep (rows 1-3)
dop_example_idx       = [dop_example_i, dop_example_i2];
dop_example_stim_idx  = [dop_example_s, dop_example_s2];
dop_example_col_label = {dop_ref_side_str, dop_opp_side_str};
dop_example_disp = struct('ex',{},'xb',{},'stim_t',{},'zoom_t0',{},'zoom_t1',{},'im_alpha',{},'atp_alpha',{},'im_clim',{},'atp_clim',{},'atp_avg_full',{});
for c = 1:2
    i = dop_example_idx(c);
    ex_c = dop_all_data(i);
    xb_c = ex_c.ft.xb;

    stim_t_c = dop_trial_events{i}(dop_example_stim_idx(c)).onset_t;
    zoom_t0_c = stim_t_c - zoom_pre_s;
    zoom_t1_c = stim_t_c + zoom_post_s;
    zoom_idx_c = xb_c >= zoom_t0_c & xb_c <= zoom_t1_c;

    dop_example_disp(c).ex           = ex_c;
    dop_example_disp(c).xb           = xb_c;
    dop_example_disp(c).stim_t       = stim_t_c;
    dop_example_disp(c).zoom_t0      = zoom_t0_c;
    dop_example_disp(c).zoom_t1      = zoom_t1_c;
    dop_example_disp(c).im_alpha     = unwrap(ex_c.im.alpha);
    dop_example_disp(c).atp_alpha    = unwrap(ex_c.atp.alpha);
    dop_example_disp(c).im_clim      = [prctile(min(ex_c.im.z(:,zoom_idx_c),[],1),5),  prctile(max(ex_c.im.z(:,zoom_idx_c),[],1),95)];
    dop_example_disp(c).atp_clim     = [prctile(min(ex_c.atp.d(:,zoom_idx_c),[],1),5), prctile(max(ex_c.atp.d(:,zoom_idx_c),[],1),95)];
    dop_example_disp(c).atp_avg_full = mean(ex_c.atp.d,1);
end

%% population traces (rows 4-5): ONE line per panel (all flies pooled, no genotype dimension)
dop_bump_color = [.2 .2 .2]; % row 4 col 1: no channel identity (bump position, not a fluorescence signal) -- neutral dark gray

dop_fly_list = unique(dop_fly_id);
dop_n_flies  = numel(dop_fly_list);

[dop_bump_mean,dop_bump_sem,dop_bump_n] = pool_all_flies(@(i,s) sign_flip(dop_trial_events{i}(s).mu_rel(trace_idx), dop_trial_is_right(i), dop_reference_is_right), dop_fly_id, dop_fly_list, dop_trial_events, dop_n_flies, numel(t_trace));
[dop_atp_mean,dop_atp_sem,dop_atp_n]    = pool_all_flies(@(i,s) dop_trial_events{i}(s).atp_dff_rel(trace_idx), dop_fly_id, dop_fly_list, dop_trial_events, dop_n_flies, numel(t_trace));
[dop_stim_epg_mean,dop_stim_epg_sem,dop_stim_epg_n] = pool_all_flies(@(i,s) side_pick(dop_trial_events{i}(s).im_z_half1(trace_idx),dop_trial_events{i}(s).im_z_half2(trace_idx),dop_trial_is_right(i),true),  dop_fly_id, dop_fly_list, dop_trial_events, dop_n_flies, numel(t_trace));
[dop_nonstim_epg_mean,dop_nonstim_epg_sem,dop_nonstim_epg_n] = pool_all_flies(@(i,s) side_pick(dop_trial_events{i}(s).im_z_half1(trace_idx),dop_trial_events{i}(s).im_z_half2(trace_idx),dop_trial_is_right(i),false), dop_fly_id, dop_fly_list, dop_trial_events, dop_n_flies, numel(t_trace));

fprintf('row 4 col 1 (bump position, %s-flip): n=%d flies\n', dop_ref_side_str, dop_bump_n);
fprintf('row 5 col 1 (whole-PB atp dF/F): n=%d flies\n', dop_atp_n);
fprintf('row 4 col 2 (EPG, stim side): n=%d flies\n', dop_stim_epg_n);
fprintf('row 5 col 2 (EPG, non-stim side): n=%d flies\n', dop_nonstim_epg_n);

%% plot: 5 rows x 2 cols, same layout as Figure_3_claude -- rows 1-3 two
% example columns (opposite sides), rows 4-5 population (single-line, no genotype split)
figure('color','w','Position',[50 50 1500 1800]); clf
t = tiledlayout(5,2,'TileSpacing','loose','Padding','compact');

dop_row1_axes = gobjects(1,2);
for c = 1:2
    d = dop_example_disp(c);
    i = dop_example_idx(c);
    axg = nexttile(t); hold(axg,'on')
    imagesc(axg,d.xb,d.im_alpha,d.ex.im.z,d.im_clim)
    colormap(axg,white_to_color(gcamp_color))
    cbg = colorbar(axg); cbg.Label.String = 'z-score';
    plot(axg,[d.stim_t,d.stim_t],[min(d.im_alpha),max(d.im_alpha)],':k','LineWidth',1)
    xlim(axg,[d.zoom_t0,d.zoom_t1]); ylim(axg,[min(d.im_alpha),max(d.im_alpha)])
    if c==1; ylabel(axg,'PB angle (rad)'); end
    xticks(axg,[])
    title(axg,sprintf('%s-stim example: EPG calcium (im.z)\nfly %s, trial "%s", stim %d', dop_example_col_label{c}, dop_fly_id{i}, meta_display(d.ex.meta), dop_example_stim_idx(c)),'Interpreter','none')
    dop_row1_axes(c) = axg;
end

dop_row2_axes = gobjects(1,2);
for c = 1:2
    d = dop_example_disp(c);
    axa = nexttile(t); hold(axa,'on')
    imagesc(axa,d.xb,d.atp_alpha,d.ex.atp.d,d.atp_clim)
    colormap(axa,white_to_color(atp_color))
    cba = colorbar(axa); cba.Label.String = 'dF/F';
    plot(axa,[d.stim_t,d.stim_t],[min(d.atp_alpha),max(d.atp_alpha)],':k','LineWidth',1)
    xlim(axa,[d.zoom_t0,d.zoom_t1]); ylim(axa,[min(d.atp_alpha),max(d.atp_alpha)])
    if c==1; ylabel(axa,'PB angle (rad)'); end
    xticks(axa,[])
    title(axa,'atp channel (atp.d, dF/F) -- dotted line marks the stim onset')
    dop_row2_axes(c) = axa;
end

dop_row3_axes = gobjects(1,2);
for c = 1:2
    d = dop_example_disp(c);
    axt = nexttile(t); hold(axt,'on')
    plot(axt,d.xb,d.atp_avg_full,'Color',atp_color,'LineWidth',1.2)
    plot(axt,[d.stim_t,d.stim_t],ylim(axt),':k')
    xlim(axt,[d.zoom_t0,d.zoom_t1])
    if c==1; ylabel(axt,'mean atp dF/F'); end
    xlabel(axt,'time (s)')
    title(axt,'this trial''s own whole-PB average atp signal')
    dop_row3_axes(c) = axt;
end

dop_ax4a = nexttile(t); hold(dop_ax4a,'on')
patch(dop_ax4a,[t_trace,fliplr(t_trace)],[dop_bump_mean+dop_bump_sem,fliplr(dop_bump_mean-dop_bump_sem)],dop_bump_color,'FaceAlpha',.15,'EdgeColor','none')
plot(dop_ax4a,t_trace,dop_bump_mean,'Color',dop_bump_color,'LineWidth',2)
plot(dop_ax4a,[0,0],ylim(dop_ax4a),':k')
plot(dop_ax4a,[trace_pre_s,trace_post_s],[0,0],':k')
xlim(dop_ax4a,[trace_pre_s,trace_post_s])
xticks(dop_ax4a,[])
ylabel(dop_ax4a,sprintf('bump pos. (rad, %s-flip)',dop_ref_side_str))
title(dop_ax4a,sprintf('average bump-position response (flipped to match the %s-stim example''s side, n=%d flies)',dop_ref_side_str,dop_bump_n))

dop_ax4b = nexttile(t); hold(dop_ax4b,'on')
patch(dop_ax4b,[t_trace,fliplr(t_trace)],[dop_stim_epg_mean+dop_stim_epg_sem,fliplr(dop_stim_epg_mean-dop_stim_epg_sem)],gcamp_color,'FaceAlpha',.15,'EdgeColor','none')
plot(dop_ax4b,t_trace,dop_stim_epg_mean,'Color',gcamp_color,'LineWidth',2)
plot(dop_ax4b,[0,0],ylim(dop_ax4b),':k')
xlim(dop_ax4b,[trace_pre_s,trace_post_s])
xticks(dop_ax4b,[])
ylabel(dop_ax4b,'mean EPG z (stim side)')
title(dop_ax4b,sprintf('average EPG calcium (im.z), STIM-SIDE wedges only (n=%d flies)',dop_stim_epg_n))

dop_ax5a = nexttile(t); hold(dop_ax5a,'on')
patch(dop_ax5a,[t_trace,fliplr(t_trace)],[dop_atp_mean+dop_atp_sem,fliplr(dop_atp_mean-dop_atp_sem)],atp_color,'FaceAlpha',.15,'EdgeColor','none')
plot(dop_ax5a,t_trace,dop_atp_mean,'Color',atp_color,'LineWidth',2)
plot(dop_ax5a,[0,0],ylim(dop_ax5a),':k')
xlim(dop_ax5a,[trace_pre_s,trace_post_s])
xlabel(dop_ax5a,'time from stim onset (s)')
ylabel(dop_ax5a,'mean atp dF/F')
title(dop_ax5a,sprintf('average whole-PB atp response to ATP ejection (n=%d flies)',dop_atp_n))

dop_ax5b = nexttile(t); hold(dop_ax5b,'on')
patch(dop_ax5b,[t_trace,fliplr(t_trace)],[dop_nonstim_epg_mean+dop_nonstim_epg_sem,fliplr(dop_nonstim_epg_mean-dop_nonstim_epg_sem)],gcamp_color,'FaceAlpha',.15,'EdgeColor','none')
plot(dop_ax5b,t_trace,dop_nonstim_epg_mean,'Color',gcamp_color,'LineWidth',2)
plot(dop_ax5b,[0,0],ylim(dop_ax5b),':k')
xlim(dop_ax5b,[trace_pre_s,trace_post_s])
xlabel(dop_ax5b,'time from stim onset (s)')
ylabel(dop_ax5b,'mean EPG z (non-stim side)')
title(dop_ax5b,sprintf('average EPG calcium (im.z), NON-STIM-SIDE wedges only (n=%d flies)',dop_nonstim_epg_n))

linkaxes([dop_row1_axes(1),dop_row2_axes(1),dop_row3_axes(1)],'x')
linkaxes([dop_row1_axes(2),dop_row2_axes(2),dop_row3_axes(2)],'x')
linkaxes([dop_ax4a,dop_ax4b,dop_ax5a,dop_ax5b],'x')
sgtitle('Figure 3\_dopamine (draft): PB heatmap around a confirmed ejection, two example trials (opposite stim sides) + population bump-position, atp, and EPG (by side) response -- single group (dopamine\_ionto\_redo)')

exportgraphics(gcf, fullfile(export_dir,'Figure_3_dopamine.png'), 'Resolution', 300)
if ~isfolder(all_figs_dir); mkdir(all_figs_dir); end
exportgraphics(gcf, fullfile(all_figs_dir,'Fig3_dopamine_V1.pdf'), 'ContentType', 'vector')

%% ===================== functions =====================

function ev = get_trial_stim_events(trial, t_grid, peri_win, base_win, peak_win, peak_factor)
    % returns a struct array, one entry per CONFIRMED-ejection stim in this
    % trial (same per-stim z-score detector as lpsp_p2x2_claude.m: this
    % stim's own atp.d peak must clear peak_factor baseline-SDs above its
    % own pre-stim baseline), with:
    %   onset_t    - this stim's onset time (ft.xf clock)
    %   right_frac - this stim's OWN 0-3s peak-window, per-wedge atp.d
    %                spatial profile, first-16-wedges share of the total --
    %                computed from just the few seconds ATP is actually
    %                being delivered, NOT diluted by the ~57s of inter-stim
    %                baseline a whole-trial average would pool in
    %   wedge_peak - that same peak-window profile, all 32 wedges (kept for
    %                display -- right_frac is just its 1:16-vs-17:32 summary)
    %   mu_rel     - bump position (unwrapped im.mu), interpolated onto
    %                t_grid and zeroed at this stim's own onset
    %   head_rel   - self-motion-integrated heading (cumsum(r_speed)), same
    %                treatment -- used by the example-trial bump-movement
    %                score (mu_rel vs. head_rel divergence)
    %   rho_grid   - im.rho interpolated onto t_grid, for confidence gating
    %   atp_dff_rel - mean atp.d ACROSS ALL 32 WEDGES (whole-PB average
    %                dF/F, not baseline-zeroed) -- dF/F rather than z-score
    %                per explicit request, so this reflects the actual
    %                ejection amplitude rather than a per-trial-
    %                renormalized shape, interpolated onto t_grid
    %   atp_peak_amp - this stim's own whole-PB atp.d (dF/F) peak amplitude:
    %                  max atp_dff_rel over 0-5s post-onset minus its mean
    %                  over -2-0s pre-onset. A stim can clear the
    %                  peak_factor ejection-confirmation test above (which
    %                  only checks the peak is a few baseline-SDs above a
    %                  possibly tiny/quiet baseline) while still barely
    %                  registering in absolute dF/F terms -- this is a
    %                  size-of-effect measure for picking a visually clean
    %                  example, not part of the ejection call itself.
    %   im_z_half1 - mean EPG calcium (im.z, z-scored -- unlike the atp
    %                fields above, EPG calcium stays z-scored per explicit
    %                request) across wedges 1:16 only, interpolated onto
    %                t_grid (not baseline-zeroed)
    %   im_z_half2 - same, wedges 17:32. Which physical half is "stim side"
    %                vs "non-stim side" depends on this TRIAL's own
    %                trial_is_right call (not knowable inside this
    %                trial-only function) -- the caller picks half1 vs
    %                half2 accordingly.
    xb = trial.ft.xb; xf = trial.ft.xf;
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    atp_sig = sum(trial.atp.d,1);
    t_onset = xf(onsets);
    mu_unwrapped = unwrap(trial.im.mu);
    dt = median(diff(xf));
    heading = cumsum(trial.ft.r_speed)*dt; % ported from lpsp_p2x2_revisit.m/lpsp_p2x2_walking_script.m
    atp_dff_avg = mean(trial.atp.d,1); % whole-PB average atp dF/F, on the xb clock
    im_z_half1_avg = mean(trial.im.z(1:16,:),1);  % EPG calcium, wedges 1:16, on the xb clock
    im_z_half2_avg = mean(trial.im.z(17:32,:),1); % EPG calcium, wedges 17:32

    t_common_local = linspace(-peri_win,peri_win,101);
    ev = struct('onset_t',{},'right_frac',{},'wedge_peak',{},'mu_rel',{},'head_rel',{},'rho_grid',{},'atp_dff_rel',{},'atp_peak_amp',{},'im_z_half1',{},'im_z_half2',{});
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s);
        xf_rel = xf - t_onset(s);
        peri = interp1(xb_rel,atp_sig,t_common_local,'linear',nan);
        base_idx = t_common_local >= base_win(1) & t_common_local < base_win(2);
        peak_idx = t_common_local > peak_win(1) & t_common_local <= peak_win(2);
        base_mean = median(peri(base_idx),'omitnan');
        base_std  = std(peri(base_idx),'omitnan');
        peak_amp  = median(peri(peak_idx),'omitnan');
        if (peak_amp-base_mean) <= peak_factor*base_std
            continue % not a confirmed ejection -- skip this stim
        end

        peak_frames = xb_rel>peak_win(1) & xb_rel<=peak_win(2);
        wedge_peak = mean(trial.atp.d(:,peak_frames),2); % 32x1, this stim's own peak-window spatial profile
        rf = mean(wedge_peak(1:16)) / (mean(wedge_peak(1:16)) + mean(wedge_peak(17:32)));

        mu_peri   = interp1(xb_rel,mu_unwrapped,t_grid,'linear',nan);
        head_peri = interp1(xf_rel,heading,t_grid,'linear',nan);
        rho_peri  = interp1(xb_rel,trial.im.rho,t_grid,'linear',nan);

        k = numel(ev)+1;
        ev(k).onset_t    = t_onset(s);
        ev(k).right_frac = rf;
        ev(k).wedge_peak = wedge_peak;
        ev(k).mu_rel     = mu_peri   - interp1(t_grid,mu_peri,0);
        ev(k).head_rel   = head_peri - interp1(t_grid,head_peri,0);
        ev(k).rho_grid    = rho_peri;
        ev(k).atp_dff_rel = interp1(xb_rel,atp_dff_avg,t_grid,'linear',nan);

        pre_idx  = t_grid>=-2 & t_grid<0;
        post_idx = t_grid>=0  & t_grid<=5;
        ev(k).atp_peak_amp = max(ev(k).atp_dff_rel(post_idx)) - mean(ev(k).atp_dff_rel(pre_idx),'omitnan');

        ev(k).im_z_half1 = interp1(xb_rel,im_z_half1_avg,t_grid,'linear',nan);
        ev(k).im_z_half2 = interp1(xb_rel,im_z_half2_avg,t_grid,'linear',nan);
    end
end

function [fid, tnum] = walk_meta_fly_trial(meta_path)
    % fly ID + trial number from an lpsp_p2x2_walking .meta path (e.g.
    % "Z:\pablo\lpsp_p2x2_walking\20260708\fly 1\20260708-1_epg_8m_..."
    % ) -- ported from data scripts/lpsp_p2x2_walking_script.m section 3,
    % split on {'\','/'} rather than '\' alone so it also works on mac
    % (same fix as trial_fly_folder_name/trial_fly_id above). The [_-] in
    % the trial-number regex (not just '_') is needed because one trial's
    % folder uses a hyphen there instead of the usual underscore (see the
    % walking script's own comment on this).
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    % true if this trial's stim-triggered average atp trace shows a
    % defined peak just after stim onset -- ported from lpsp_p2x2_walking_
    % script.m section 3 (trial-level pulse detection: one call per trial,
    % from that trial's OWN stim-averaged trace, unlike get_trial_stim_events
    % above which confirms each stim individually)
    if ~(isfield(trial,'atp') && isfield(trial.atp,'d') && ~isempty(trial.atp.d))
        tf = false; return
    end
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    if isempty(onsets); tf = false; return; end

    xb = trial.ft.xb; xf = trial.ft.xf;
    t_onset = xf(onsets);
    atp_sig = sum(trial.atp.d,1);
    t_common = linspace(-peri_win,peri_win,101);

    atp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s);
        atp_aligned(s,:) = interp1(xb_rel,atp_sig,t_common,'linear',nan);
    end
    atp_peri = mean(atp_aligned,1,'omitnan');

    base_idx = t_common>=base_win(1) & t_common<base_win(2);
    peak_idx = t_common>peak_win(1) & t_common<=peak_win(2);
    base_mean = median(atp_peri(base_idx),'omitnan');
    base_std  = std(atp_peri(base_idx),'omitnan');
    peak_amp  = median(atp_peri(peak_idx),'omitnan');
    tf = (peak_amp-base_mean) > peak_factor*base_std;
end

function [mov_mu, mov_cue, dur, rho_mean] = walk_trial_bouts(trial, smooth_window, turn_thresh, max_gap_frames, min_walking_frames)
    % per-walking-bout bump path length (mov_mu) vs. heading path length
    % (mov_cue, despite the name -- built from ft.r_speed, not cue
    % position, so an experimenter-driven cue change while the fly isn't
    % walking doesn't show up as spurious "heading" movement) and bout
    % duration (dur) -- ported verbatim from lpsp_p2x2_walking_script.m
    % section 2. rho_mean (added for the row 6 rho-delta panel) is the
    % bout's mean UNSMOOTHED im.rho on the behavior clock -- same
    % construction as Figure_3_scratch.m section 4's ungated bout rho.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    mu = interp1(trial.ft.xb,unwrap(trial.im.mu),xf,'linear','extrap');
    rho = interp1(trial.ft.xb,trial.im.rho,xf,'linear','extrap');

    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    mov_mu  = nan(numel(bout_starts),1);
    mov_cue = nan(numel(bout_starts),1);
    dur     = nan(numel(bout_starts),1);
    rho_mean = nan(numel(bout_starts),1);
    for b = 1:numel(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(r_speed_smooth(bout_starts(b):bout_ends(b))),'omitnan')*dt;
        dur(b)     = (bout_ends(b)-bout_starts(b)+1)*dt;
        rho_mean(b) = mean(rho(bout_starts(b):bout_ends(b)),'omitnan');
    end
end

function name = trial_fly_folder_name(meta_path)
    % ported from data scripts/lpsp_p2x2_claude.m: the fly folder ("fly N
    % _ <genotype>") one level above the trial folder. Split on both '\'
    % and '/' (not just filesep) -- these .meta paths were saved from a
    % Windows machine (e.g. "Z:\pablo\...\fly 1 _ empty\..."), so on mac
    % filesep alone ('/') never splits them and every trial's genotype
    % silently comes back empty. Matches meta_display's own split below.
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        name = '';
    else
        name = parts{fly_part(1)};
    end
end

function fid = trial_fly_id(meta_path)
    % ported from data scripts/lpsp_p2x2_claude.m -- see trial_fly_folder_name
    % above for why this splits on {'\','/'} rather than filesep
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        fid = meta_path;
    else
        fid = strjoin(parts(1:fly_part(1)),'/');
    end
end

function s = meta_display(meta_path)
    % last path component, UNLESS it's the generic "registration" folder
    % that EVERY trial in both this figure's datasets ends in -- in that
    % case use the trial folder one level up instead, which is the only
    % place a real trial identifier (date-number, or in dopamine_ionto_redo's
    % case the only place a trial NUMBER at all) lives. For p2x2 this mostly
    % just swaps "registration" for that trial's own numbered folder name;
    % for dopamine_ionto_redo it's the difference between every title saying
    % the same uninformative "registration" and actually showing which of a
    % fly's several trials is pictured (confirmed directly: two different
    % example trials from the same fly both showed "registration" until
    % this fix was added).
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if strcmpi(parts{end},'registration') && numel(parts)>1
        s = parts{end-1};
    else
        s = parts{end};
    end
end

function iw = walk_bout_mask(trial, smooth_window, turn_thresh, max_gap_frames, min_walking_frames)
    % walking-bout mask on the ft.xf clock -- same definition and parameters
    % as walk_trial_bouts (Figure 3 row 6), for shading example traces
    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    iw = abs(r_speed_smooth) > turn_thresh;
    iw = imclose(iw,ones(max_gap_frames+1,1));
    iw = bwareaopen(iw,min_walking_frames);
end

function add_time_scalebar(ax, len_s)
    % replace a time axis with a len_s scale bar under the axes' bottom-right
    % corner. Call after the axes' limits are final (freezes them).
    xl = xlim(ax); yl = ylim(ax); xlim(ax,xl); ylim(ax,yl);
    ax.XColor = 'none';
    if strcmp(ax.YDir,'reverse'); y = yl(2) + .06*diff(yl); dy = .04*diff(yl); else; y = yl(1) - .06*diff(yl); dy = -.04*diff(yl); end
    x = xl(2) - [len_s 0];
    plot(ax, x, [y y], 'k-', 'LineWidth', 2, 'Clipping', 'off', 'HandleVisibility', 'off')
    text(ax, mean(x), y + dy, sprintf('%g s', len_s), 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'Clipping', 'off')
end

function set_pi_yticks(ax)
    % y ticks at multiples of pi/2 across the axis's current limits, labeled
    % in pi (TeX) -- same label style as Figure_1_claude.m's angle axes.
    % Freezes the y-limits first so the ticks can't drift if they rescale.
    % Step widens (pi/2 -> pi -> 2pi) so a large range keeps <= ~8 ticks.
    yl = ylim(ax); ylim(ax,yl);
    q = 1; while diff(yl)/(q*pi/2) > 8; q = q*2; end % step = q*pi/2
    k = q*(ceil(yl(1)/(q*pi/2)):floor(yl(2)/(q*pi/2)));
    yticks(ax,k*pi/2);
    yticklabels(ax,arrayfun(@pi_half_label,k,'UniformOutput',false));
end

function s = pi_half_label(k)
    % TeX label for k*pi/2: 0, \pi/2, \pi, 3\pi/2, 2\pi, -\pi/2, ...
    if k == 0; s = '0'; return; end
    if mod(k,2) == 0
        n = k/2; d = '';
    else
        n = k; d = '/2';
    end
    if n == 1; s = ['\pi' d]; elseif n == -1; s = ['-\pi' d]; else; s = sprintf('%d\\pi%s',n,d); end
end

function add_delta_pvals(ax, vals, grp)
    % p-values on a 2-group delta (post - pre) panel: two-sided signrank of
    % each group's per-fly deltas vs. 0 above that group, and a two-sided
    % ranksum between the groups on a bracket above both -- same tests as
    % Figure_3_scratch.m's tables.
    vals = vals(:); grp = grp(:);
    top = max(vals); span = max(vals) - min([vals; 0]);
    for g = 1:2
        p = signrank(vals(grp==g));
        text(ax, g, max(vals(grp==g)) + .08*span, sprintf('p=%.2g',p), 'HorizontalAlignment','center','FontSize',9)
    end
    p = ranksum(vals(grp==1), vals(grp==2));
    yb = top + .22*span;
    plot(ax, [1 1 2 2], yb - [.03 0 0 .03]*span, 'k', 'LineWidth', 1, 'HandleVisibility','off')
    text(ax, 1.5, yb + .03*span, sprintf('p=%.2g',p), 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9)
    yl = ylim(ax); ylim(ax, [yl(1), max(yl(2), yb + .15*span)])
end

function len = nice_scale_len(x)
    % largest 1/2/5 x 10^k value <= x -- a round number for a scale bar
    % about as long as x
    p = 10^floor(log10(x));
    m = [1 2 5 10];
    len = p*m(find(m*p <= x,1,'last'));
end

function rgb = trial_mean_rgb(meta, gcamp_color, atp_color, planes)
    % whole-trial mean image of one trial, as a two-channel RGB composite
    % on a white background, same white->color convention as the heatmaps
    % (white_to_color): each channel subtracts its complement from white,
    % so gcamp alone goes white->gcamp_color, atp alone white->atp_color,
    % and overlap darkens toward black.
    % Both channels come from the SAME file
    % (registration\imagingData_trial001.mat, img{1} = gcamp, img{2} = atp,
    % per-plane 100x256x10xframes int16, ~3 GB, ~20 s to load) so they're
    % guaranteed to share one registration -- its z-sum is identical to
    % registration\imagingData.mat, which only exists pre-summed, so the
    % per-plane file is needed to sum just the requested z-planes. Cropped
    % to the trial's mask.mat bounding box
    % (+pad) and rotated 90 deg counter-clockwise: process_im orders the
    % glomeruli from the mask skeleton's first (= leftmost, column-major
    % find) endpoint, so image-left = glomerulus 1 = bottom of the
    % heatmap (YDir normal), which is exactly where rot90 puts it. Each
    % channel is clipped to its own 1st-99.5th (gcamp: 95th) percentile over the PB mask
    % pixels only (per explicit request -- the atp pipette tip is ~4x
    % brighter than the ejection cloud, and stray bright structure outside
    % the PB would otherwise set the gcamp range).
    S = load(fullfile(meta,'imagingData_trial001.mat'),'img');
    g = sum(mean(S.img{1}(:,:,planes,:),4),3);
    a = sum(mean(S.img{2}(:,:,planes,:),4),3);
    clear S
    mk = load(fullfile(fileparts(meta),'mask.mat'),'mask');
    [ym,xm] = find(mk.mask);
    pad = 5;
    rows = max(1,min(ym)-pad):min(size(g,1),max(ym)+pad);
    cols = max(1,min(xm)-pad):min(size(g,2),max(xm)+pad);
    g = fliplr(rot90(g(rows,cols))); a = fliplr(rot90(a(rows,cols))); % fliplr: mirrored about the vertical axis, per explicit request (leaves the glomerulus order bottom-to-top unchanged)
    m = fliplr(rot90(mk.mask(rows,cols)));
    norm01 = @(x,ref,hi) min(max((x-prctile(ref,1))/(prctile(ref,hi)-prctile(ref,1)),0),1);
    % gcamp ceiling at the 95th percentile (not 99.5th), per explicit
    % request -- brings up the dimmer glomeruli
    rgb = max(1 - norm01(g,g(m),95).*reshape(1-gcamp_color,1,1,3) - norm01(a,a(m),99.5).*reshape(1-atp_color,1,1,3), 0);
end

function rough = atp_trace_roughness(ex, onset_t)
    % how wiggly one stim's whole-PB atp dF/F trace is, -2..20s around
    % onset_t (the row-1 display window): SD of the trace around a
    % 1s-gaussian-smoothed copy of itself, divided by the smoothed trace's
    % peak above its own -2..0s baseline -- so a big ejection isn't
    % penalized just for having bigger absolute frame-to-frame noise.
    xb = ex.ft.xb;
    idx = xb>=onset_t-2 & xb<=onset_t+20;
    tr = mean(ex.atp.d(:,idx),1);
    dt = median(diff(xb));
    sm = smoothdata(tr,'gaussian',round(1/dt));
    rough = std(tr-sm) / (max(sm) - mean(sm(1:round(2/dt))));
end

function cmap = white_to_color(max_color)
    % linear colormap from white (low) to max_color (high) -- ported
    % verbatim from Figure_1_claude.m/Figure_2_claude.m
    n = 256;
    cmap = [linspace(1,max_color(1),n)', linspace(1,max_color(2),n)', linspace(1,max_color(3),n)'];
end

function x = wrap_to_pi(x)
    % rewrap to (-pi,pi], same convention used throughout this codebase
    x = mod(x,2*pi);
    x(x>pi) = x(x>pi) - 2*pi;
end

function ev = get_trial_stim_events_dopamine(trial, t_grid, peri_win, base_win, peak_win, peak_factor)
    % same as get_trial_stim_events above, EXCEPT right_frac/wedge_peak
    % clip the peak-window atp.d response to positive values only before
    % splitting per hemisphere half -- confirmed necessary specifically on
    % dopamine_ionto_redo_20260518.mat (see the Figure_3_dopamine section's
    % header, point 3): a raw dip elsewhere in the PB is not evidence ATP
    % didn't land there, and left unclipped it was free to swing the
    % hemisphere balance on this dataset's smaller, noisier signal.
    xb = trial.ft.xb; xf = trial.ft.xf;
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    atp_sig = sum(trial.atp.d,1);
    t_onset = xf(onsets);
    mu_unwrapped = unwrap(trial.im.mu);
    dt = median(diff(xf));
    heading = cumsum(trial.ft.r_speed)*dt;
    atp_dff_avg = mean(trial.atp.d,1);
    im_z_half1_avg = mean(trial.im.z(1:16,:),1);
    im_z_half2_avg = mean(trial.im.z(17:32,:),1);

    t_common_local = linspace(-peri_win,peri_win,101);
    ev = struct('onset_t',{},'right_frac',{},'wedge_peak',{},'mu_rel',{},'head_rel',{},'rho_grid',{},'atp_dff_rel',{},'atp_peak_amp',{},'im_z_half1',{},'im_z_half2',{});
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s);
        xf_rel = xf - t_onset(s);
        peri = interp1(xb_rel,atp_sig,t_common_local,'linear',nan);
        base_idx = t_common_local >= base_win(1) & t_common_local < base_win(2);
        peak_idx = t_common_local > peak_win(1) & t_common_local <= peak_win(2);
        base_mean = median(peri(base_idx),'omitnan');
        base_std  = std(peri(base_idx),'omitnan');
        peak_amp  = median(peri(peak_idx),'omitnan');
        if (peak_amp-base_mean) <= peak_factor*base_std
            continue
        end

        peak_frames = xb_rel>peak_win(1) & xb_rel<=peak_win(2);
        wedge_peak = mean(trial.atp.d(:,peak_frames),2);
        wedge_peak_pos = wedge_peak; wedge_peak_pos(wedge_peak_pos<0) = 0; % positive-clipped, this dataset's fix
        rf = sum(wedge_peak_pos(1:16)) / sum(wedge_peak_pos);

        mu_peri   = interp1(xb_rel,mu_unwrapped,t_grid,'linear',nan);
        head_peri = interp1(xf_rel,heading,t_grid,'linear',nan);
        rho_peri  = interp1(xb_rel,trial.im.rho,t_grid,'linear',nan);

        k = numel(ev)+1;
        ev(k).onset_t     = t_onset(s);
        ev(k).right_frac  = rf;
        ev(k).wedge_peak  = wedge_peak;
        ev(k).mu_rel      = mu_peri   - interp1(t_grid,mu_peri,0);
        ev(k).head_rel    = head_peri - interp1(t_grid,head_peri,0);
        ev(k).rho_grid    = rho_peri;
        ev(k).atp_dff_rel = interp1(xb_rel,atp_dff_avg,t_grid,'linear',nan);

        pre_idx  = t_grid>=-2 & t_grid<0;
        post_idx = t_grid>=0  & t_grid<=5;
        ev(k).atp_peak_amp = max(ev(k).atp_dff_rel(post_idx)) - mean(ev(k).atp_dff_rel(pre_idx),'omitnan');

        ev(k).im_z_half1 = interp1(xb_rel,im_z_half1_avg,t_grid,'linear',nan);
        ev(k).im_z_half2 = interp1(xb_rel,im_z_half2_avg,t_grid,'linear',nan);
    end
end

function v = sign_flip(v, this_is_right, reference_is_right)
    % negate v unless this trial's own side call already matches the
    % reference side (used for Figure_3_dopamine's bump-position pooling,
    % same convention as the p2x2 section's inline trial_is_right~=reference_is_right flip)
    if this_is_right ~= reference_is_right
        v = -v;
    end
end

function v = side_pick(half1, half2, is_right, want_stim_side)
    % want_stim_side=true -> return the STIM-SIDE half for this trial;
    % false -> return the NON-STIM-SIDE half (used for Figure_3_dopamine's
    % EPG-by-side pooling)
    if is_right == want_stim_side
        v = half1;
    else
        v = half2;
    end
end

function [g_mean,g_sem,g_n] = pool_all_flies(get_trace_fn, fly_id, fly_list, trial_events, n_flies, t_len)
    % one value per fly first (that fly's own confirmed-ejection stims,
    % across all its trials, averaged together via get_trace_fn), then
    % averaged across flies -- same one-point-per-fly convention as the
    % p2x2 section's per-genotype pooling above, just with no genotype
    % dimension to loop over (used by Figure_3_dopamine, which has none).
    fly_traces = nan(n_flies,t_len);
    for f = 1:n_flies
        trial_list = find(strcmp(fly_id,fly_list{f}));
        stim_traces = [];
        for i = trial_list(:)'
            ev = trial_events{i};
            for s = 1:numel(ev)
                stim_traces(end+1,:) = get_trace_fn(i,s); %#ok<AGROW>
            end
        end
        if ~isempty(stim_traces)
            fly_traces(f,:) = mean(stim_traces,1,'omitnan');
        end
    end
    have = ~all(isnan(fly_traces),2);
    M = fly_traces(have,:);
    g_n = size(M,1);
    g_mean = mean(M,1,'omitnan');
    g_sem  = std(M,[],1,'omitnan') ./ sqrt(sum(~isnan(M),1));
end
