%% lpsp_p2x2_claude
% Loads the LPsP>P2X2 dataset (empty>P2X2 control vs lpsp>P2X2, ATP
% iontophoresis onto the PB to trigger P2X2-mediated depolarization) and
% reports how many flies of each genotype are in it, following the same
% genotype/fly-count pattern as lpsp_kir_claude.m and lpsp_tnt_claude.m.
%
% This dataset was investigated directly (not assumed to match kir/tnt)
% before writing this script. It differs from both in two ways:
%
% 1) GENOTYPE TAG LOCATION: unlike kir ("..._LPsP_kir"/"..._empty_kir")
%    and tnt ("..._lpsp_tnt"/"..._empty_tnt"), where the tag sits in the
%    TRIAL folder name itself, here the trial folder name carries no
%    genotype info at all (e.g. "20260428-1_epg_syt8s_blind_p2x2" --
%    confirmed directly, same suffix on every trial regardless of
%    genotype). The tag instead sits one level up, in the FLY folder name
%    ("fly 1 _ empty", "fly 2 _ lpsp", etc.), confirmed directly against
%    all 101 trials in lpsp_p2x2_reredo_20260518.mat with zero ambiguity
%    and zero within-fly genotype disagreement. Spacing around the tag is
%    inconsistent in the raw data ("fly 1 _ empty" vs "fly 7 _empty"), so
%    matching below is a simple case-insensitive contains() on the fly
%    folder string, not a fixed-format regex.
%
% 2) LIGHT CONDITION IS ALREADY STORED: unlike lpsp_kir_redo (0/74 trials
%    had ft.pattern), every trial here already has ft.pattern populated
%    (101/101, confirmed directly) as either "0001_background.mat" (dark)
%    or "0003_4px_brightbar.mat" (closed loop) -- no need for the
%    trialSettings.csv re-matching lpsp_kir_claude.m/lpsp_tnt_claude.m had
%    to do. Both light conditions are present for both genotypes in
%    roughly balanced numbers (confirmed directly: empty 26 cl/24 dark,
%    lpsp 26 cl/25 dark).
%
% Also note this dataset's raw tree sits at a "\blind\" level not present
% in lpsp_p2x2_reredo's other data pulls (lpsp_p2x2_reredo_20260331.mat,
% _20260501.mat) -- all_data.meta here points into
% Z:\pablo\lpsp_p2x2_reredo\blind\<date>\fly N _ <genotype>\<trial>\registration,
% i.e. one level ABOVE the trial folder is the genotype-tagged fly folder,
% and meta's own last path component is always the bare "registration"
% folder (not "registration_NNN" like some other datasets in this
% codebase), confirmed directly across all 101 trials.

%% load data
data_dir    = fullfile('.data');
source_file = 'lpsp_p2x2_reredo_20260518.mat';
base_dir    = 'Z:\pablo\lpsp_p2x2_reredo\blind\';
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype per trial, from the "fly N _ <genotype>" folder name in all_data.meta
% (NOT the trial folder name itself -- see header. Matched against the fly
% folder only, so a base path segment like "lpsp_p2x2_reredo" can never
% falsely match, same false-positive risk already handled for
% lpsp_kir_claude.m/lpsp_tnt_claude.m.)
genotype = cell(1,n_trials);
for i = 1:n_trials
    fname = trial_fly_folder_name(all_data(i).meta);
    if contains(fname,'lpsp','IgnoreCase',true)
        genotype{i} = 'lpsp>p2x2';
    elseif contains(fname,'empty','IgnoreCase',true)
        genotype{i} = 'empty>p2x2';
    else
        genotype{i} = '';
    end
end
n_unlabeled = sum(cellfun(@isempty,genotype));
if n_unlabeled > 0
    warning('%d/%d trials did not match either genotype pattern in their fly folder name', n_unlabeled, n_trials)
end

%% light condition per trial, read directly from the already-stored ft.pattern
has_pattern = arrayfun(@(s) isfield(s.ft,'pattern') && ~isempty(s.ft.pattern), all_data);
fprintf('\ntrials with ft.pattern stored: %d/%d\n', sum(has_pattern), n_trials);

is_dark = nan(1,n_trials);
for i = 1:n_trials
    if has_pattern(i)
        is_dark(i) = contains(all_data(i).ft.pattern,'background','IgnoreCase',true);
    else
        warning('trial %d (%s) has no ft.pattern', i, all_data(i).meta);
    end
end

fprintf('\n=== light condition summary ===\n');
fprintf('closed loop: %d\n', sum(is_dark==0));
fprintf('dark:        %d\n', sum(is_dark==1));
fprintf('unresolved:  %d\n', sum(isnan(is_dark)));

%% fly ID per trial (date folder + "fly N _ genotype"), and a genotype-consistency check
fly_id = cell(1,n_trials);
for i = 1:n_trials
    fly_id{i} = trial_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
fly_num = fly_num(:)'; % keep row-vector, matching is_dark/genotype
n_flies = numel(fly_list);

fly_genotype = cell(n_flies,1);
for f = 1:n_flies
    g = unique(genotype(fly_num==f));
    g(cellfun(@isempty,g)) = [];
    if numel(g) > 1
        warning('fly %s has inconsistent genotype labels across its trials: %s', fly_list{f}, strjoin(g,', '));
    end
    if ~isempty(g)
        fly_genotype{f} = g{1};
    else
        fly_genotype{f} = '';
    end
end
fprintf('\n%d trials -> %d flies\n', n_trials, n_flies);

%% flies per genotype, and trials per genotype x light condition
geno_order = {'empty>p2x2','lpsp>p2x2'};
fly_counts = zeros(1,numel(geno_order));
for k = 1:numel(geno_order)
    fly_counts(k) = sum(strcmp(fly_genotype,geno_order{k}));
end

fprintf('\n=== flies per genotype ===\n');
for k = 1:numel(geno_order)
    fprintf('  %-10s n=%d flies\n', geno_order{k}, fly_counts(k));
end
n_unlabeled_flies = sum(cellfun(@isempty,fly_genotype));
if n_unlabeled_flies > 0
    fprintf('  %-10s n=%d flies (no genotype match)\n', 'unlabeled', n_unlabeled_flies);
end

fprintf('\n=== trials per genotype x light condition ===\n');
cond_label = {'closed loop','dark'};
for k = 1:numel(geno_order)
    for c = 1:2
        n = sum(strcmp(genotype,geno_order{k}) & is_dark==(c-1));
        fprintf('  %-10s %-11s n=%d trials\n', geno_order{k}, cond_label{c}, n);
    end
end

%% plot: number of flies per genotype
figure(1); clf
set(gcf,'Name','flies per genotype')
bar(fly_counts)
set(gca,'XTickLabel',geno_order)
ylabel('number of flies')
title(sprintf('lpsp\\_p2x2\\_reredo (blind): flies per genotype (n=%d flies, %d trials)', n_flies, n_trials))
for k = 1:numel(geno_order)
    text(k, fly_counts(k)+0.1, num2str(fly_counts(k)), 'HorizontalAlignment','center')
end

%% detect ATP ejections per individual stim: ground truth (ft.stims rising edge) vs. observed atp response
% lpsp_p2x2_walking_script.m's own perturbation detector (its section "%%
% 3) label each trial") calls a whole TRIAL a "perturbation trial" if the
% stim-triggered AVERAGE atp trace (across however many stims that trial
% has) clears a baseline-SD threshold just after the stim onset. That
% collapses every stim in a trial into one decision. Here every
% INDIVIDUAL stim is scored on its own: for each rising edge of
% ft.stims (the ground-truth record that a stim command was sent,
% regardless of whether ATP was actually ejected), the atp channel
% (sum(atp.d,1), same channel-summing convention as elsewhere in this
% codebase) is pulled into a peri-stim window on its own relative-time
% axis, and that ONE stim's own pre-stim baseline (median/SD, not
% mean/max -- a single noisy sample shouldn't drive the call, same
% reasoning as the walking script) is compared against its own post-stim
% peak.
%
% Confirmed directly on this dataset before picking parameters/threshold:
% 101 trials all have ft.stims and a real (non-empty) atp field; every
% trial has exactly 4 stims (404 total), spaced ~60s apart, well clear of
% each trial's own start/end (first onset >=41s in, last onset leaves
% >=3s of trial after it) -- so no window here is truncated by a trial
% boundary. With peak_factor=3 (the walking script's own threshold,
% unchanged), 396/404 stims (98%) are called ejections with z-scores in
% the hundreds; the 8 that aren't all belong to exactly 2 trials where
% EVERY one of that trial's 4 stims failed (a real all-or-nothing
% pipette/ejection failure, not scattered borderline calls) -- see the
% per-trial failure list printed below.
peri_win    = 5;      % seconds before/after each stim's rising edge to examine
t_common    = linspace(-peri_win,peri_win,101); % common relative-time axis every stim gets interpolated onto
base_win    = [-peri_win,0]; % seconds before onset used as this stim's own baseline
peak_win    = [0,3];         % seconds after onset checked for a response peak
peak_factor = 3;             % post-stim peak must clear this many baseline-SDs above baseline to count as "an ejection"

stim_trial       = [];      % index into all_data this stim belongs to
stim_num         = [];      % 1-based order of this stim within its own trial
stim_onset_t     = [];      % onset time, ft.xf clock
stim_base_mean   = [];
stim_base_std    = [];
stim_peak_amp    = [];
stim_z           = [];       % (peak - baseline) / baseline_sd
stim_is_ejection = false(1,0);
stim_peri        = {};       % this stim's own peri-stim atp trace, on t_common (kept for plotting)
stim_im_peri     = {};       % this stim's own peri-stim EPG calcium (im.d) trace, on t_common, same clock as atp.d (confirmed: both 32 x n_frames on xb for every trial)
stim_im_base     = [];       % this stim's own pre-stim baseline of the calcium trace (its own median, not the atp baseline)

for i = 1:n_trials
    xb = all_data(i).ft.xb;
    xf = all_data(i).ft.xf;
    stims  = logical(all_data(i).ft.stims(:));
    onsets = find(diff([false;stims])==1); % rising edge of each stim command
    if isempty(onsets); continue; end

    atp_sig = sum(all_data(i).atp.d,1);
    im_sig  = sum(all_data(i).im.d,1);
    t_onset = xf(onsets);

    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s); % imaging clock, re-centered on this one stim's own onset
        peri    = interp1(xb_rel,atp_sig,t_common,'linear',nan);
        im_peri = interp1(xb_rel,im_sig, t_common,'linear',nan);

        base_idx = t_common >= base_win(1) & t_common <  base_win(2);
        peak_idx = t_common >  peak_win(1) & t_common <= peak_win(2);

        base_mean = median(peri(base_idx),'omitnan');
        base_std  = std(peri(base_idx),'omitnan');
        peak_amp  = median(peri(peak_idx),'omitnan');

        stim_trial(end+1)       = i; %#ok<AGROW>
        stim_num(end+1)         = s; %#ok<AGROW>
        stim_onset_t(end+1)     = t_onset(s); %#ok<AGROW>
        stim_base_mean(end+1)   = base_mean; %#ok<AGROW>
        stim_base_std(end+1)    = base_std; %#ok<AGROW>
        stim_peak_amp(end+1)    = peak_amp; %#ok<AGROW>
        stim_z(end+1)           = (peak_amp-base_mean)/base_std; %#ok<AGROW>
        stim_is_ejection(end+1) = (peak_amp-base_mean) > peak_factor*base_std; %#ok<AGROW>
        stim_peri{end+1}        = peri; %#ok<AGROW>
        stim_im_peri{end+1}     = im_peri; %#ok<AGROW>
        stim_im_base(end+1)     = median(im_peri(base_idx),'omitnan'); %#ok<AGROW>
    end
end

n_stims = numel(stim_trial);
trials_with_stims = unique(stim_trial);

fprintf('\n=== ATP ejection detection (per individual stim) ===\n');
fprintf('%d stim(s) found across %d/%d trials\n', n_stims, numel(trials_with_stims), n_trials);
fprintf('ejection detected: %d/%d stims (%.0f%%)\n', sum(stim_is_ejection), n_stims, 100*mean(stim_is_ejection));

stim_genotype = genotype(stim_trial);
stim_dark     = is_dark(stim_trial);
fprintf('\n=== stims per genotype x light condition ===\n');
for k = 1:numel(geno_order)
    for c = 1:2
        idx = strcmp(stim_genotype,geno_order{k}) & stim_dark==(c-1);
        if any(idx)
            fprintf('  %-10s %-11s %d/%d stims labeled ejection\n', geno_order{k}, cond_label{c}, sum(stim_is_ejection(idx)), sum(idx));
        end
    end
end

%% trials where EVERY stim failed to produce a detected ejection
full_fail_trials = trials_with_stims(arrayfun(@(t) all(~stim_is_ejection(stim_trial==t)), trials_with_stims));
fprintf('\n%d/%d trials had every one of their stims fail to produce a detected ejection:\n', numel(full_fail_trials), numel(trials_with_stims));
for t = full_fail_trials
    fprintf('  trial %d: %s\n', t, all_data(t).meta);
end

%% figure: overview of every individual stim's response, colored by ejection decision
[z_sorted,z_order] = sort(stim_z);
ej_sorted = stim_is_ejection(z_order);
ej_color  = [0.85,0.33,0.10]; % ejection
no_color  = [0.00,0.45,0.74]; % no ejection

figure(2); clf
set(gcf,'Name','atp ejection detection: overview','Position',[100,100,900,450])

subplot(1,2,1); hold on
scatter(find(~ej_sorted),z_sorted(~ej_sorted),25,no_color,'filled')
scatter(find(ej_sorted), z_sorted(ej_sorted), 25,ej_color, 'filled')
yline(peak_factor,':k')
xlabel('stim rank (sorted by z)')
ylabel('z = (peak - baseline) / baseline SD')
legend({'no ejection','ejection','threshold'},'Location','northwest')
title(sprintf('all %d stims (full range)',n_stims))

subplot(1,2,2); hold on
scatter(find(~ej_sorted),z_sorted(~ej_sorted),35,no_color,'filled')
scatter(find(ej_sorted), z_sorted(ej_sorted), 35,ej_color, 'filled')
yline(peak_factor,':k')
ylim([-20,20])
xlabel('stim rank (sorted by z)')
ylabel('z = (peak - baseline) / baseline SD')
title('zoomed to the decision boundary (|z|<20)')

sgtitle(sprintf('ejection detection margin: %d/%d stims labeled ejection (peak\\_factor=%g)',sum(stim_is_ejection),n_stims,peak_factor))

%% figure: one panel per fly, all of that fly's stims overlaid on the same axes
% condensed from an earlier version of this section (one subplot per
% STIM, up to 27 figures total) per explicit request -- this instead
% gives each FLY one panel, tiled together into a single figure, with all
% of that fly's peri-stim atp traces overlaid directly on top of each
% other, colored by ejection decision.
%
% each trace is shown as a deviation from ITS OWN baseline
% (peri - stim_base_mean), not raw atp units: raw baseline level varies
% trial-to-trial (e.g. trial 92's ~1.1-1.4 vs trial 1's ~0.4-1.0,
% confirmed directly) and would otherwise stack traces at vertical offsets
% that have nothing to do with the stim response itself -- subtracting
% each stim's own baseline lines every trace up at 0 pre-stim, so the
% post-stim deviation (the thing the ejection call is actually based on)
% is what's visually comparable across stims and across flies.
fly_order   = [find(strcmp(fly_genotype,'empty>p2x2'));find(strcmp(fly_genotype,'lpsp>p2x2'))]; % genotype-grouped tile order
n_fly_tiles = numel(fly_order);
cols_fly    = ceil(sqrt(n_fly_tiles));
rows_fly    = ceil(n_fly_tiles/cols_fly);
overlay_ylim = [-5,20]; % fixed across every tile so panels are directly comparable -- covers the ejected-response range with margin (checked directly against stim_peak_amp-stim_base_mean above)

figure(3); clf
set(gcf,'Name','atp ejection QC: per-fly overlay','Position',[50,50,220*cols_fly,170*rows_fly])
tl = tiledlayout(rows_fly,cols_fly,'TileSpacing','compact','Padding','compact');

for tile_ix = 1:n_fly_tiles
    f = fly_order(tile_ix);
    stim_idx_f = find(ismember(stim_trial,find(fly_num==f)));
    ax = nexttile(tl); hold(ax,'on')
    if isempty(stim_idx_f)
        axis(ax,'off')
        continue
    end

    patch(ax,peak_win([1,2,2,1]),overlay_ylim([1,1,2,2]),'g','FaceAlpha',.08,'EdgeColor','none')
    plot(ax,[t_common(1),t_common(end)],[0,0],':k')
    plot(ax,[0,0],overlay_ylim,':k')

    for j = stim_idx_f
        trace = stim_peri{j} - stim_base_mean(j); % this stim's own deviation from its own baseline
        if stim_is_ejection(j)
            c = [0.85,0.33,0.10,.7];
        else
            c = [0.00,0.45,0.74,.7];
        end
        plot(ax,t_common,trace,'Color',c,'LineWidth',1.1)
    end

    xlim(ax,[t_common(1),t_common(end)]); ylim(ax,overlay_ylim)
    n_ej_f = sum(stim_is_ejection(stim_idx_f));
    title(ax,sprintf('%s\n%d/%d ejected',fly_short_label(fly_list{f}),n_ej_f,numel(stim_idx_f)),'FontSize',7,'Interpreter','none')
end

h_ej = plot(nan,nan,'Color',[0.85,0.33,0.10],'LineWidth',1.5);
h_no = plot(nan,nan,'Color',[0.00,0.45,0.74],'LineWidth',1.5);
legend([h_ej,h_no],{'ejection','no ejection'},'Location','southoutside')
title(tl,'atp response per stim, overlaid per fly (deviation from that stim''s own baseline; green = peak window)','Interpreter','none')

%% figure: EPG calcium (im.d) response around confirmed ATP ejections, overlaid per fly
% same per-fly overlay style as figure 3, but on the imaging channel
% (im.d, summed across glomeruli, same convention as atp_sig above) and
% restricted to only the 396 stims already confirmed as real ejections
% (stim_is_ejection) -- the question here is whether EPG calcium changed
% given that ATP was actually delivered, not a repeat of the ejection call
% itself. im.d's own peri-stim trace (stim_im_peri) was pulled on the same
% xb-clock window as the atp trace in the same loop above (confirmed
% im.d/atp.d share the same 32 x n_frames shape/clock every trial), and is
% shown here as a deviation from ITS OWN pre-stim baseline (stim_im_base
% -- the calcium channel's own baseline, not atp's), for the same reason
% as figure 3: raw baseline fluorescence varies trial-to-trial and would
% otherwise stack traces at offsets unrelated to the stim response.
%
% checked directly before fixing the y-axis: across all 396 confirmed
% ejections, the calcium deviation ranges from about -4 to +19 (1st-99th
% pctile), so the same [-5,20] window used for the atp overlay also covers
% this comfortably and keeps the two figures on a directly comparable scale.
ej_idx = find(stim_is_ejection);

figure(4); clf
set(gcf,'Name','epg calcium QC: per-fly overlay (confirmed ejections only)','Position',[50,50,220*cols_fly,170*rows_fly])
tl2 = tiledlayout(rows_fly,cols_fly,'TileSpacing','compact','Padding','compact');

for tile_ix = 1:n_fly_tiles
    f = fly_order(tile_ix);
    stim_idx_f = ej_idx(ismember(stim_trial(ej_idx),find(fly_num==f))); % this fly's confirmed-ejection stims
    ax = nexttile(tl2); hold(ax,'on')
    if isempty(stim_idx_f)
        axis(ax,'off')
        title(ax,sprintf('%s\n(no confirmed ejections)',fly_short_label(fly_list{f})),'FontSize',7,'Interpreter','none')
        continue
    end

    patch(ax,peak_win([1,2,2,1]),overlay_ylim([1,1,2,2]),'g','FaceAlpha',.08,'EdgeColor','none')
    plot(ax,[t_common(1),t_common(end)],[0,0],':k')
    plot(ax,[0,0],overlay_ylim,':k')

    for j = stim_idx_f
        trace = stim_im_peri{j} - stim_im_base(j); % this stim's own calcium deviation from its own baseline
        plot(ax,t_common,trace,'Color',[0.13,0.55,0.13,.6],'LineWidth',1.1)
    end

    xlim(ax,[t_common(1),t_common(end)]); ylim(ax,overlay_ylim)
    title(ax,sprintf('%s\nn=%d ejections',fly_short_label(fly_list{f}),numel(stim_idx_f)),'FontSize',7,'Interpreter','none')
end

title(tl2,'EPG calcium (im.d) response, overlaid per fly, confirmed ATP ejections only (deviation from that stim''s own baseline; green = peak window)','Interpreter','none')

%% which PB hemisphere each trial's ATP pipette targeted, from the atp channel's own raw fluorescence spatial profile
% same first-half/second-half convention lpsp_p2x2_reredo_script.m's and
% lpsp_p2x2_script_minimal.m's own right_idx used (e.g. script_minimal.m:
% "sum(x(1:size(x,1)/2,:),'all') > sum(x(size(x,1)/2:end,:),'all')"), on
% atp.f (raw fluorescence, summed across the WHOLE trial) rather than a
% single stim's window -- confirmed directly that this dataset's 32-row
% atp channel is the same two-hemisphere-repeated-16-wedge layout those
% scripts assumed (im.alpha(1:16) == im.alpha(17:32) exactly), and that
% "right" (first half) vs "left" (second half) is a property of the
% TRIAL, not the fly: 24/25 flies got stims on BOTH sides across their own
% trials (confirmed directly), so side is computed and grouped per trial
% below, not assumed fixed per fly.
%
% this split isn't perfectly clean -- checked directly: the first-half
% fraction of total atp fluorescence ranges from 0.17 to 0.79 across
% trials, with 12/101 trials landing in an ambiguous 0.4-0.6 band. Flagged
% below, not excluded, since the source scripts this convention is copied
% from don't filter on it either.
% all_data is a COLUMN struct array, so arrayfun/cellfun over it return
% column vectors by default -- reshaped to rows below (.') to match every
% other per-trial array in this script (fly_num, is_dark, genotype), since
% otherwise indexing them with the same row-shaped index (e.g.
% stim_trial(ej_idx)) returns mismatched orientations: MATLAB preserves
% the INDEXED array's own orientation, not the index's, so a column
% trial_is_right(...) combined via & with a row fly_num(...) silently
% broadcasts into an NxN matrix instead of comparing elementwise -- caught
% directly (this is what errored at the line combining them below).
sum_atp_f  = arrayfun(@(x) sum(x.atp.f,2), all_data, 'UniformOutput',false);
right_frac = cellfun(@(v) mean(v(1:numel(v)/2)) / (mean(v(1:numel(v)/2))+mean(v(numel(v)/2+1:end))), sum_atp_f).';
trial_is_right = right_frac > 0.5;
n_ambiguous = sum(right_frac > 0.4 & right_frac < 0.6);
fprintf('\n%d/%d trials have an ambiguous (0.4-0.6) right/left atp fluorescence split\n', n_ambiguous, n_trials);

%% average bump position (im.mu), zeroed at the stim rising edge -- confirmed ejections only, split by genotype x stim side
% mu is unwrapped once per trial (im.mu is circular; unwrapping per-trial
% before windowing keeps continuity across the stim onset), then each
% confirmed-ejection stim's own peri-stim window is interpolated onto
% t_bump and shifted so its OWN value at t=0 is exactly 0 -- "bump
% position relative to where it was at the stim rising edge", per
% request, not an absolute heading.
%
% a wider window than the +/-5s atp/calcium windows above is used here
% (-5s to +15s): a bump-position effect of exciting/silencing LPsP is
% expected to play out more slowly than the fast atp/calcium kinetics,
% matching the wider pulse-triggered-average windows
% lpsp_p2x2_reredo_script.m/lpsp_p2x2_script_minimal.m used for their own
% bump-position analyses (win_start/win_end of -2/20 and -5/10
% respectively).
%
% traces are averaged ONE VALUE PER FLY first (a fly's own qualifying
% stims, within one genotype x side group, averaged together), THEN
% averaged across flies for the plotted mean +/- SEM -- same
% one-point-per-fly pooling convention used throughout this codebase
% (e.g. lpsp_kir_claude.m), so a fly with more confirmed ejections on one
% side doesn't get more weight than a fly with fewer.
bump_pre = -5; bump_post = 15;
t_bump   = linspace(bump_pre,bump_post,201);

mu_unwrapped = cell(1,n_trials);  % unwrap(im.mu), computed once per trial, lazily
cue_unwrapped = cell(1,n_trials); % unwrap(-ft.cue), same sign flip used throughout this codebase (e.g. lpsp_p2x2_walking_script.m, lpsp_kir_claude.m) so cue is directly comparable to mu
stim_mu_rel  = cell(size(stim_trial)); % filled only for confirmed-ejection stims (ej_idx)
stim_cue_rel = cell(size(stim_trial));

for jj = 1:numel(ej_idx)
    j = ej_idx(jj);
    i = stim_trial(j);
    if isempty(mu_unwrapped{i})
        mu_unwrapped{i} = unwrap(all_data(i).im.mu);
    end
    if isempty(cue_unwrapped{i})
        cue_unwrapped{i} = unwrap(-all_data(i).ft.cue);
    end

    xb_rel  = all_data(i).ft.xb - stim_onset_t(j);
    mu_peri = interp1(xb_rel,mu_unwrapped{i},t_bump,'linear',nan);
    stim_mu_rel{j} = mu_peri - interp1(t_bump,mu_peri,0); % zero this stim's own trace at its own onset

    xf_rel   = all_data(i).ft.xf - stim_onset_t(j); % cue lives on the fictrac clock (xf), not the imaging clock (xb)
    cue_peri = interp1(xf_rel,cue_unwrapped{i},t_bump,'linear',nan);
    stim_cue_rel{j} = cue_peri - interp1(t_bump,cue_peri,0);
end

side_order  = {'right','left'};
group_defs2 = struct('geno',{},'right',{},'label',{});
for gi = 1:numel(geno_order)
    for si = 1:2
        group_defs2(end+1) = struct('geno',geno_order{gi},'right',(si==1),'label',sprintf('%s, stim %s',geno_order{gi},side_order{si})); %#ok<SAGROW>
    end
end

[bump_mean,bump_sem,bump_n] = pool_group_traces(ej_idx,stim_mu_rel, group_defs2,fly_num,stim_trial,trial_is_right,fly_genotype,n_flies,numel(t_bump));
[cue_mean, cue_sem, ~     ] = pool_group_traces(ej_idx,stim_cue_rel,group_defs2,fly_num,stim_trial,trial_is_right,fly_genotype,n_flies,numel(t_bump)); % same stim pool as bump above, so cue_n == bump_n -- legend gating below uses bump_n only

fprintf('\n=== average bump position & cue, zeroed at stim onset (confirmed ejections only) ===\n');
for gd = 1:numel(group_defs2)
    fprintf('  %-22s n=%d flies\n', group_defs2(gd).label, bump_n(gd));
end

%% figure: the 4 bump traces above, overlaid on one set of axes, with the fly's own heading cue overlaid in gray
empty_color = [0.00,0.45,0.74];
lpsp_color  = [0.85,0.33,0.10];
side_style  = {'-','--'}; % right, left

figure(5); clf
set(gcf,'Name','average bump position aligned to stim onset','Position',[100,100,700,550])
plot_bump_cue_overlay(gca,t_bump,bump_mean,bump_sem,bump_n,cue_mean,cue_sem,group_defs2,empty_color,lpsp_color,side_style);
xlabel('time from stim onset (s)')
ylabel('position relative to stim onset (rad)')
title('average bump position (color) vs. heading cue (gray), zeroed at the atp stim rising edge (mean +/- SEM across flies, confirmed ejections only)')

%% figure: the same 4 bump+cue traces, split into two panels by light condition (closed loop vs. dark)
figure(6); clf
set(gcf,'Name','average bump position aligned to stim onset, by light condition','Position',[100,100,1300,550])
stim_is_dark = is_dark(stim_trial); % is_dark is per-trial; broadcast onto the per-stim arrays used to build ej_idx-restricted pools below

for c = 1:2 % 1 = closed loop, 2 = dark
    pool = ej_idx(stim_is_dark(ej_idx)==(c-1));
    [bm,bs,bn] = pool_group_traces(pool,stim_mu_rel, group_defs2,fly_num,stim_trial,trial_is_right,fly_genotype,n_flies,numel(t_bump));
    [cm,cs,~ ] = pool_group_traces(pool,stim_cue_rel,group_defs2,fly_num,stim_trial,trial_is_right,fly_genotype,n_flies,numel(t_bump));

    ax = subplot(1,2,c);
    plot_bump_cue_overlay(ax,t_bump,bm,bs,bn,cm,cs,group_defs2,empty_color,lpsp_color,side_style);
    xlabel(ax,'time from stim onset (s)')
    if c == 1
        ylabel(ax,'position relative to stim onset (rad)')
    end
    title(ax,cond_label{c})
end
sgtitle('average bump position (color) vs. heading cue (gray), zeroed at stim onset, split by light condition (confirmed ejections only)')

%% figure: bump movement collapsed across stim side, by negating mu (and cue) for left-side stims
% same light-condition-collapsed pool as figure 5, but now ALSO collapsing
% right-side and left-side stims onto one axis: since the two PB
% hemispheres are mirror images of each other, negating a left-side
% stim's own already-zeroed bump trace (stim_mu_rel) before pooling puts
% it on the same sign convention as a right-side stim's trace -- so
% "positive" now consistently means "bump moved in the direction expected
% from a right-side stim" regardless of which side actually got stimulated
% that trial. Heading cue is negated the same way, for the same reason.
% This doubles the stims contributing to each of the resulting 2 traces
% (empty, lpsp) relative to any one of figure 5's 4 side-specific traces.
stim_mu_rel_flip  = stim_mu_rel;
stim_cue_rel_flip = stim_cue_rel;
for jj = 1:numel(ej_idx)
    j = ej_idx(jj);
    if ~trial_is_right(stim_trial(j)) % left-side stim -- flip onto the right-side sign convention
        stim_mu_rel_flip{j}  = -stim_mu_rel{j};
        stim_cue_rel_flip{j} = -stim_cue_rel{j};
    end
end

flip_bump_mean = nan(2,numel(t_bump)); flip_bump_sem = nan(2,numel(t_bump)); flip_bump_n = zeros(1,2);
flip_cue_mean  = nan(2,numel(t_bump)); flip_cue_sem  = nan(2,numel(t_bump));

fprintf('\n=== average bump position, zeroed at stim onset, collapsed across stim side (confirmed ejections only) ===\n');
for gi = 1:numel(geno_order)
    fly_mu  = {};
    fly_cue = {};
    for f = 1:n_flies
        if ~strcmp(fly_genotype{f},geno_order{gi}); continue; end
        idx = ej_idx(fly_num(stim_trial(ej_idx))==f); % both sides now, no side filter -- that's the point of the flip above
        if isempty(idx); continue; end
        fly_mu{end+1}  = mean(cell2mat(stim_mu_rel_flip(idx)'), 1,'omitnan'); %#ok<SAGROW>
        fly_cue{end+1} = mean(cell2mat(stim_cue_rel_flip(idx)'),1,'omitnan'); %#ok<SAGROW>
    end
    M_mu  = cell2mat(fly_mu');
    M_cue = cell2mat(fly_cue');
    flip_bump_n(gi)      = size(M_mu,1);
    flip_bump_mean(gi,:) = mean(M_mu, 1,'omitnan');
    flip_bump_sem(gi,:)  = std(M_mu, 1,'omitnan')/sqrt(max(flip_bump_n(gi),1));
    flip_cue_mean(gi,:)  = mean(M_cue,1,'omitnan');
    flip_cue_sem(gi,:)   = std(M_cue,1,'omitnan')/sqrt(max(flip_bump_n(gi),1));
    fprintf('  %-22s n=%d flies\n', geno_order{gi}, flip_bump_n(gi));
end

% reuses plot_bump_cue_overlay: a 2-entry group_defs with 'right' fixed
% true for both just selects one consistent line style (there's no side
% left to distinguish once collapsed) without needing a separate plotting
% routine.
flip_group_defs = struct('geno',geno_order(:)','right',{true,true},'label',geno_order(:)');

figure(7); clf
set(gcf,'Name','bump movement collapsed across stim side','Position',[100,100,700,550])
plot_bump_cue_overlay(gca,t_bump,flip_bump_mean,flip_bump_sem,flip_bump_n,flip_cue_mean,flip_cue_sem,flip_group_defs,empty_color,lpsp_color,side_style);
xlabel('time from stim onset (s)')
ylabel('position relative to stim onset (rad), left-side stims sign-flipped')
title({'average bump position (color) vs. heading cue (gray), zeroed at stim onset,','collapsed across stim side by negating left-side stims (confirmed ejections only, both light conditions pooled)'})

%% figure: average glomerulus (im.d) heatmap around confirmed ATP ejections, by genotype x stim side
% mirrors lpsp_p2x2_walking_script.m's own section "21) same average
% glomerulus heatmap as figure 18, but keep left- and right-stimulated
% flies separate (no PB flip)" -- a 2x2 grid, rows = stim side, cols =
% genotype, each panel a 32-glomerulus x time heatmap of im.d (raw,
% unflipped -- physical hemisphere order preserved, that's the point of
% keeping side separate rather than mirrored), baseline-subtracted per
% row at its own t=0.
%
% two differences from that source section, both because this script
% already has the pieces built: aligned to the per-STIM confirmed
% ejections (ej_idx) rather than a per-TRIAL average (the walking script
% only detected ejections at the trial level), and side comes from
% trial_is_right (this trial's own atp fluorescence spatial profile)
% rather than a fixed per-fly flip sign, since most flies here get stims
% on both sides across their own trials (see the section that computed
% trial_is_right, above).
zero_idx = find(t_common==0,1);

heatmap_pool = cell(2,2); % {side, genotype}: side 1 = right-stim trials, 2 = left-stim trials; each accumulates a 32 x length(t_common) x n_stims stack
for jj = 1:numel(ej_idx)
    j = ej_idx(jj);
    i = stim_trial(j);

    xb_rel  = all_data(i).ft.xb - stim_onset_t(j);
    aligned = interp1(xb_rel,all_data(i).im.d',t_common,'linear',nan)'; % 32 x length(t_common)
    aligned = aligned - aligned(:,zero_idx); % baseline-subtract each glomerulus row to its own t=0 value

    side = 1 + ~trial_is_right(i);                    % 1 = right-stim trial, 2 = left-stim trial
    gi   = 1 + strcmp(genotype{i},'lpsp>p2x2');        % 1 = empty>p2x2, 2 = lpsp>p2x2 (matches geno_order)
    heatmap_pool{side,gi} = cat(3,heatmap_pool{side,gi},aligned);
end

mean_heatmap = cell(2,2);
n_heatmap    = zeros(2,2);
for side = 1:2
    for gi = 1:2
        n_heatmap(side,gi) = size(heatmap_pool{side,gi},3);
        if n_heatmap(side,gi) > 0
            mean_heatmap{side,gi} = mean(heatmap_pool{side,gi},3,'omitnan');
        end
    end
end
clim_max = max(cellfun(@(m) max(abs(m(:)),[],'omitnan'), mean_heatmap(~cellfun(@isempty,mean_heatmap))));

n_col = 256;
diverge_cmap = [linspace(0,1,n_col/2)',linspace(0,1,n_col/2)',ones(n_col/2,1); ... % blue -> white
                ones(n_col/2,1),linspace(1,0,n_col/2)',linspace(1,0,n_col/2)'];    % white -> red

side_labels = {'right-stim trials','left-stim trials'};
hemi_labels = {{'stimulated hemisphere','non-stim hemisphere'},{'non-stim hemisphere','stimulated hemisphere'}}; % {side}{wedge 1:16 label, 17:32 label}

figure(8); clf
set(gcf,'Name','glomerulus dF/F heatmap by genotype x stim side','Position',[100,100,900,700])
for side = 1:2
    for gi = 1:2
        subplot(2,2,(side-1)*2+gi); hold on
        if isempty(mean_heatmap{side,gi}); continue; end

        imagesc(t_common,1:32,mean_heatmap{side,gi})
        plot([t_common(1),t_common(end)],[16.5,16.5],'k:','LineWidth',1) % hemisphere boundary
        plot([0,0],[.5,32.5],'k:','LineWidth',1) % stim onset
        colormap(gca,diverge_cmap)
        caxis([-clim_max,clim_max])
        axis tight
        yticks([8.5,24.5]); yticklabels(hemi_labels{side})
        xlabel('time from stim onset (s)')
        c = colorbar; c.Label.String = 'dF/F (baseline-subtracted)';
        title(sprintf('%s, %s (n=%d stims)',side_labels{side},geno_order{gi},n_heatmap(side,gi)))
    end
end
sgtitle('average EPG calcium (im.d) heatmap around confirmed ATP ejections, by genotype x stim side (no PB flip applied)')

%% save all figures as PDF
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\p2x2_dangle';
if ~isfolder(fig_dir)
    mkdir(fig_dir)
end

fig_handles = findobj('Type','figure');
[~,order] = sort(arrayfun(@(f) f.Number, fig_handles));
fig_handles = fig_handles(order); % ascending by figure number, not creation/stacking order

save_failed = {};
for k = 1:numel(fig_handles)
    fig = fig_handles(k);
    fig_name = get(fig,'Name');
    if isempty(fig_name)
        fig_name = sprintf('figure_%d',fig.Number);
    end
    safe_name = regexprep(fig_name,'[^\w\-]+','_'); % filesystem-safe filename
    out_path  = fullfile(fig_dir,sprintf('fig%02d_%s.pdf',fig.Number,safe_name));
    try
        exportgraphics(fig,out_path)
        fprintf('saved %s\n', out_path);
    catch ME
        % e.g. the PDF is open in a viewer and locked for writing -- don't
        % let one such failure abort saving every other figure
        save_failed{end+1} = out_path; %#ok<AGROW>
        warning('could not save %s (%s) -- is it open in another program?', out_path, ME.message);
    end
end
if ~isempty(save_failed)
    fprintf('\n%d figure(s) failed to save -- close them in any viewer and re-run to update:\n', numel(save_failed));
    fprintf('  %s\n', save_failed{:});
end

%% functions
function [g_mean,g_sem,g_n] = pool_group_traces(stim_pool,trace_cell,group_defs,fly_num,stim_trial,trial_is_right,fly_genotype,n_flies,t_len)
    % one value per fly first (that fly's own stims within stim_pool that
    % belong to this genotype x side group, averaged together), then
    % averaged across flies for the returned mean/SEM -- see the comment
    % above the section that calls this, and lpsp_kir_claude.m for the
    % same pooling convention used elsewhere in this codebase.
    fly_trace = cell(n_flies,numel(group_defs));
    for gd = 1:numel(group_defs)
        for f = 1:n_flies
            if ~strcmp(fly_genotype{f},group_defs(gd).geno); continue; end
            idx = stim_pool(fly_num(stim_trial(stim_pool))==f & trial_is_right(stim_trial(stim_pool))==group_defs(gd).right);
            if isempty(idx); continue; end
            fly_trace{f,gd} = mean(cell2mat(trace_cell(idx)'),1,'omitnan');
        end
    end

    n_groups = numel(group_defs);
    g_mean = nan(n_groups,t_len);
    g_sem  = nan(n_groups,t_len);
    g_n    = zeros(1,n_groups);
    for gd = 1:n_groups
        have = ~cellfun(@isempty,fly_trace(:,gd));
        if ~any(have); continue; end % no fly had a qualifying stim in this pool for this group -- leave as NaN, handled by the plotting helper
        M = cell2mat(fly_trace(have,gd));
        g_n(gd) = size(M,1);
        g_mean(gd,:) = mean(M,1,'omitnan');
        g_sem(gd,:)  = std(M,1,'omitnan')/sqrt(max(g_n(gd),1));
    end
end

function plot_bump_cue_overlay(ax,t_axis,bump_mean,bump_sem,bump_n,cue_mean,cue_sem,group_defs,empty_color,lpsp_color,side_style)
    % draws every group's bump trace (colored by genotype, styled by
    % side) and its matching heading-cue trace (gray, same side style),
    % each as a mean line + shaded SEM band, into the given axes. groups
    % with no data (bump_n(gd)==0) are skipped rather than plotted as a
    % flat line of NaNs.
    axes(ax); hold(ax,'on')
    h = gobjects(1,numel(group_defs));
    h_cue = gobjects(0);
    for gd = 1:numel(group_defs)
        if bump_n(gd) == 0; continue; end
        if strcmp(group_defs(gd).geno,'empty>p2x2')
            c = empty_color;
        else
            c = lpsp_color;
        end
        ls = side_style{2-group_defs(gd).right};

        patch(ax,[t_axis,fliplr(t_axis)],[cue_mean(gd,:)+cue_sem(gd,:),fliplr(cue_mean(gd,:)-cue_sem(gd,:))], ...
            [.5,.5,.5],'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off')
        h_cue(end+1) = plot(ax,t_axis,cue_mean(gd,:),'Color',[.5,.5,.5],'LineStyle',ls,'LineWidth',1.5); %#ok<AGROW>

        patch(ax,[t_axis,fliplr(t_axis)],[bump_mean(gd,:)+bump_sem(gd,:),fliplr(bump_mean(gd,:)-bump_sem(gd,:))], ...
            c,'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off')
        h(gd) = plot(ax,t_axis,bump_mean(gd,:),'Color',c,'LineStyle',ls,'LineWidth',2);
    end
    plot(ax,[t_axis(1),t_axis(end)],[0,0],':k','HandleVisibility','off')
    plot(ax,[0,0],ylim(ax),':k','HandleVisibility','off')

    have = bump_n > 0;
    legend_h = h(have);
    legend_labels = arrayfun(@(gd) sprintf('%s (n=%d)',group_defs(gd).label,bump_n(gd)), find(have),'UniformOutput',false);
    if ~isempty(h_cue)
        legend_h = [legend_h,h_cue(1)];
        legend_labels = [legend_labels,{'heading cue (gray)'}];
    end
    legend(ax,legend_h,legend_labels,'Location','northwest','Interpreter','none','FontSize',7)
end

function lbl = fly_short_label(fly_path)
    % "<date> fly N _ genotype" from a fly_id path like ...\<date folder>\fly N _ genotype
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function name = trial_fly_folder_name(meta_path)
    % the fly folder ("fly N _ <genotype>") one level above the trial
    % folder, which is itself one level above the bare "registration"
    % folder that all_data.meta points to in this dataset (confirmed: no
    % "registration_NNN" variants here, unlike some other datasets in this
    % codebase -- see header).
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        name = ''; % shouldn't happen; caller treats this as unlabeled
    else
        name = parts{fly_part(1)};
    end
end

function fid = trial_fly_id(meta_path)
    % meta = ...\<date folder>\fly N _ <genotype>\<trial folder>\registration
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        fid = meta_path; % shouldn't happen; fall back to a unique id per trial
    else
        fid = strjoin(parts(1:fly_part(1)),filesep); % date folder + "fly N _ genotype"
    end
end
