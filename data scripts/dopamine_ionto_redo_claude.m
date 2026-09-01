%% dopamine_ionto_redo_claude
% Loads the dopamine_ionto_redo dataset (empty>P2X2, ATP iontophoresis
% used to trigger P2X2-mediated depolarization of dopaminergic input near
% the PB, EPG-syt8s imaging of the bump) and reproduces the same
% per-stim ATP-ejection detection / EPG-response analyses as
% lpsp_p2x2_claude.m -- minus the genotype comparison that script needed,
% since this dataset only has one experimental group.
%
% Investigated directly before writing this script. What differs from
% lpsp_p2x2_claude.m, and why some of that script's sections are dropped:
%
% 1) SINGLE GENOTYPE, NO SPLIT NEEDED: every one of the 20 trials' meta
%    path carries "empty_p2x2_da" in its trial folder name and none carry
%    "lpsp" (confirmed directly: 20/20 "empty", 0 "lpsp"). There's only
%    one experimental group here, so none of lpsp_p2x2_claude.m's
%    genotype bookkeeping (genotype{}, fly_genotype, geno_order,
%    group_defs2, ...) is needed -- every pooling step below is simply
%    per-fly-then-across-flies.
%
% 2) FLY FOLDER HAS NO GENOTYPE TAG: the fly folder here is the bare
%    "fly N" (no "_ <genotype>" suffix -- nothing to tag, confirmed
%    directly against all 20 trials, e.g.
%    "Z:\pablo\dopamine_ionto_redo\20260506\fly 1\20260506-1_epg_syt8s_empty_p2x2_da\registration").
%    trial_fly_id/fly_short_label below are unchanged from
%    lpsp_p2x2_claude.m (same regex on the "fly\s*\d+" folder), since a
%    bare "fly N" still matches that pattern fine. Fly numbering resets
%    per session date (e.g. "fly 1" exists on both 20260506 and 20260515
%    as two different physical flies), so fly identity is still (date
%    folder + fly folder), same convention as the reference script.
%
% 3) HEMISPHERE SPLIT NEEDED A DIFFERENT METRIC: the whole-trial right_frac
%    diagnostic lpsp_p2x2_claude.m used (atp.f, RAW fluorescence, summed
%    across the ENTIRE trial, first-half/second-half split) is far more
%    ambiguous here than in that dataset -- confirmed directly: 13/20
%    trials (65%) land in its 0.4-0.6 ambiguous band, vs. 12/101 (12%)
%    there. Investigated directly why: trial 1's own atp.f spatial
%    profile across all 32 rows is 407-492 (flat to within ~20%) -- raw
%    fluorescence here is dominated by a static/background dye level
%    present everywhere in the frame, not by where ATP actually got
%    ejected, so a whole-trial average mostly averages together noise.
%    Fixed below by using each trial's own STIM-EVOKED response instead
%    (peak-window, baseline-subtracted atp.d dF/F, positive part only,
%    summed per hemisphere half) -- confirmed directly this is
%    essentially unambiguous on this dataset (0/20 trials in the 0.4-0.6
%    band, 18/20 clear >0.94 or <0.05), so the heatmap and bump-position
%    analyses below ARE split by hemisphere, using that corrected metric
%    (trial_is_right), unlike an earlier version of this script.
%
% Also confirmed directly, informing what else is/isn't here:
%  - all 20 trials already have ft.pattern == "0003_4px_brightbar.mat"
%    (closed loop) -- no dark trials, so no light-condition split either.
%  - every trial has exactly 4 stims (80 total), at the same ~60s-apart
%    onsets (41/101/161/221s into a ~300s trial), well clear of trial
%    boundaries -- same stim timing as lpsp_p2x2_claude.m's dataset, so
%    that script's peri-stim window choices (peri_win/base_win/peak_win)
%    are reused unchanged.
%  - with the same peak_factor=3 ejection threshold (also unchanged),
%    65/80 stims (81%) are called ejections; the 15 that aren't fall
%    entirely within 2 trials (trial 2 and trial 10) where every one of
%    that trial's 4 stims failed -- same all-or-nothing pipette-failure
%    pattern as the reference dataset, not scattered borderline calls.

%% load data
data_dir    = fullfile('.data');
source_file = 'dopamine_ionto_redo_20260518.mat';
base_dir    = 'Z:\pablo\dopamine_ionto_redo\';
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% fly ID per trial (date folder + "fly N" folder), and fly count
fly_id = cell(1,n_trials);
for i = 1:n_trials
    fly_id{i} = trial_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
fly_num = fly_num(:)'; % row vector, matching every other per-trial array below (stim_trial, etc.)
n_flies = numel(fly_list);

trials_per_fly = arrayfun(@(f) sum(fly_num==f), 1:n_flies);
fprintf('\n%d trials -> %d flies\n', n_trials, n_flies);
fprintf('\n=== trials per fly ===\n');
for f = 1:n_flies
    fprintf('  %-30s n=%d trials\n', fly_short_label(fly_list{f}), trials_per_fly(f));
end

%% plot: number of trials per fly (answers "how many flies are in the dataset")
figure(1); clf
set(gcf,'Name','flies in dataset','Position',[100,100,700,400])
bar(trials_per_fly,'FaceColor',[0.00,0.45,0.74])
set(gca,'XTick',1:n_flies,'XTickLabel',arrayfun(@(f) fly_short_label(fly_list{f}),1:n_flies,'UniformOutput',false),'XTickLabelRotation',45)
ylabel('number of trials')
title(sprintf('dopamine\\_ionto\\_redo: n=%d flies, %d trials total',n_flies,n_trials))
for f = 1:n_flies
    text(f,trials_per_fly(f)+0.1,num2str(trials_per_fly(f)),'HorizontalAlignment','center')
end

%% detect ATP ejections per individual stim: ground truth (ft.stims rising edge) vs. observed atp response
% same per-stim detection as lpsp_p2x2_claude.m's own section of the same
% name -- see there for the full method rationale (each stim is scored on
% its own pre-stim baseline median/SD vs. its own post-stim peak, rather
% than collapsing a trial's stims into one average). Parameters reused
% unchanged (see header for why: same stim count/spacing/timing here).
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
stim_im_peri     = {};       % this stim's own peri-stim EPG calcium (im.d) trace, on t_common, same clock as atp.d
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

%% figure: per-trial timeline of the ft.stims command vs. the atp channel's own response
% directly answers "when do ejections happen": each trial's own ft.stims
% command trace (ground truth that a stim was SENT, on the fictrac clock
% xf) is drawn alongside the atp channel's own signal (whether ATP was
% actually EJECTED, on the imaging clock xb, normalized to its own trial
% max so every panel is on a comparable 0-1 scale), with a marker at each
% stim's own onset colored by its ejection call -- so a trial where the
% command fired but the channel never responded (trials 2 and 10, both
% full pipette failures per above) is visually obvious against a trial
% where the channel tracks the command.
n_col_t = ceil(sqrt(n_trials));
n_row_t = ceil(n_trials/n_col_t);

figure(3); clf
set(gcf,'Name','ft.stims command vs. atp channel response, per trial','Position',[50,50,220*n_col_t,150*n_row_t])
tl0 = tiledlayout(n_row_t,n_col_t,'TileSpacing','compact','Padding','compact');

for i = 1:n_trials
    ax = nexttile(tl0); hold(ax,'on')
    xf      = all_data(i).ft.xf;
    xb      = all_data(i).ft.xb;
    atp_sig = sum(all_data(i).atp.d,1);

    plot(ax,xf,all_data(i).ft.stims,'Color',[.7,.7,.7],'LineWidth',1)
    plot(ax,xb,atp_sig/max(atp_sig),'Color','k','LineWidth',1)

    idx     = stim_trial==i;
    onset_i = stim_onset_t(idx);
    ej_i    = stim_is_ejection(idx);
    for k = 1:numel(onset_i)
        if ej_i(k); c = ej_color; else; c = no_color; end
        plot(ax,onset_i(k),1.1,'v','Color',c,'MarkerFaceColor',c,'MarkerSize',5)
    end
    xlim(ax,[0,max(xf)]); ylim(ax,[-.1,1.25])
    title(ax,fly_short_label(all_data(i).meta),'FontSize',6,'Interpreter','none')
end
title(tl0,'ft.stims command (gray) vs. normalized atp channel response (black); triangle = per-stim ejection call (orange=ejection, blue=none)','Interpreter','none')

%% debug: per-trial peri-stim heatmaps (atp channel + EPG calcium), and a corrected hemisphere-targeting metric
% see header (point 3) for why the whole-trial atp.f right_frac metric
% failed here. This section builds the fix and the two diagnostic
% figures used to check it: for every trial, average that trial's own
% stims (ALL of them, including trials 2 and 10 where every stim failed
% the ejection call -- the question here is that trial's own spatial
% response, independent of whether it cleared the ejection threshold)
% into one 32-glomerulus x time heatmap, on both the atp.d channel (where
% did the ejection itself land) and the im.d channel (where did the EPG
% calcium respond), baseline-subtracted at each row's own t=0 -- the same
% per-stim computation the pooled heatmap section below uses, just kept
% per-trial instead of pooled across trials.
%
% the corrected hemisphere metric (evoked_right_frac) comes from the atp
% heatmap's own post-stim peak window (same peak_win as the ejection
% detector above): mean response per glomerulus row in that window,
% clipped to positive values only (an ejection makes dF/F go UP; only the
% increase should count toward "where did it respond"), summed per
% hemisphere half (rows 1-16 vs 17-32) and normalized to a fraction.
zero_idx      = find(t_common==0,1);
peak_idx_heat = t_common > peak_win(1) & t_common <= peak_win(2);

trial_atp_heat    = cell(1,n_trials); % this trial's own stims averaged into one 32 x length(t_common) heatmap, atp.d
trial_im_heat     = cell(1,n_trials); % same, im.d
evoked_right_frac = nan(1,n_trials);

for i = 1:n_trials
    idx = find(stim_trial==i);
    if isempty(idx); continue; end
    atp_stack = nan(32,numel(t_common),numel(idx));
    im_stack  = nan(32,numel(t_common),numel(idx));
    for k = 1:numel(idx)
        j = idx(k);
        xb_rel = all_data(i).ft.xb - stim_onset_t(j);
        a = interp1(xb_rel,all_data(i).atp.d',t_common,'linear',nan)';
        m = interp1(xb_rel,all_data(i).im.d', t_common,'linear',nan)';
        atp_stack(:,:,k) = a - a(:,zero_idx);
        im_stack(:,:,k)  = m - m(:,zero_idx);
    end
    trial_atp_heat{i} = mean(atp_stack,3,'omitnan');
    trial_im_heat{i}  = mean(im_stack, 3,'omitnan');

    resp = mean(trial_atp_heat{i}(:,peak_idx_heat),2); % 32 x 1, evoked atp response per glomerulus
    resp(resp<0) = 0;
    evoked_right_frac(i) = sum(resp(1:16)) / sum(resp);
end
trial_is_right = evoked_right_frac > 0.5;

sum_atp_f        = arrayfun(@(x) sum(x.atp.f,2), all_data, 'UniformOutput',false);
whole_right_frac = cellfun(@(v) mean(v(1:numel(v)/2)) / (mean(v(1:numel(v)/2))+mean(v(numel(v)/2+1:end))), sum_atp_f).';

fprintf('\n=== hemisphere targeting: whole-trial atp.f metric (old) vs. stim-evoked atp.d metric (fixed) ===\n');
fprintf('%-6s %-12s %-12s\n','trial','whole_frac','evoked_frac');
for i = 1:n_trials
    fprintf('%-6d %-12.3f %-12.3f\n', i, whole_right_frac(i), evoked_right_frac(i));
end
fprintf('whole-trial ambiguous (0.4-0.6): %d/%d\n', sum(whole_right_frac>0.4 & whole_right_frac<0.6), n_trials);
fprintf('stim-evoked ambiguous (0.4-0.6): %d/%d\n', sum(evoked_right_frac>0.4 & evoked_right_frac<0.6), n_trials);

n_col_cmap = 256;
diverge_cmap = [linspace(0,1,n_col_cmap/2)',linspace(0,1,n_col_cmap/2)',ones(n_col_cmap/2,1); ... % blue -> white
                ones(n_col_cmap/2,1),linspace(1,0,n_col_cmap/2)',linspace(1,0,n_col_cmap/2)'];    % white -> red

%% figure: per-trial ATP channel heatmap, aligned to stims, own color scale per panel
figure(4); clf
set(gcf,'Name','debug atp channel heatmap per trial','Position',[50,50,220*n_col_t,180*n_row_t])
tl4 = tiledlayout(n_row_t,n_col_t,'TileSpacing','compact','Padding','compact');
for i = 1:n_trials
    ax = nexttile(tl4); hold(ax,'on')
    h = trial_atp_heat{i};
    imagesc(ax,t_common,1:32,h)
    plot(ax,[t_common(1),t_common(end)],[16.5,16.5],'k:','LineWidth',1)
    plot(ax,[0,0],[.5,32.5],'k:','LineWidth',1)
    colormap(ax,diverge_cmap)
    clim_i = max(abs(h(:)),[],'omitnan'); if isempty(clim_i) || clim_i==0; clim_i = 1; end
    caxis(ax,[-clim_i,clim_i])
    axis(ax,'tight'); yticks(ax,[])
    if trial_is_right(i); side_lbl = 'R'; else; side_lbl = 'L'; end
    title(ax,sprintf('t%d %s\n(%s, %.2f)',i,fly_short_label(all_data(i).meta),side_lbl,evoked_right_frac(i)),'FontSize',6,'Interpreter','none')
end
title(tl4,'debug: per-trial ATP channel (atp.d) heatmap aligned to stims, own color scale per panel (top=rows17-32, bottom=rows1-16; label=evoked-response side call)','Interpreter','none')

%% figure: per-trial EPG calcium heatmap, aligned to stims, own color scale per panel
figure(5); clf
set(gcf,'Name','debug epg calcium heatmap per trial','Position',[50,50,220*n_col_t,180*n_row_t])
tl5 = tiledlayout(n_row_t,n_col_t,'TileSpacing','compact','Padding','compact');
for i = 1:n_trials
    ax = nexttile(tl5); hold(ax,'on')
    h = trial_im_heat{i};
    imagesc(ax,t_common,1:32,h)
    plot(ax,[t_common(1),t_common(end)],[16.5,16.5],'k:','LineWidth',1)
    plot(ax,[0,0],[.5,32.5],'k:','LineWidth',1)
    colormap(ax,diverge_cmap)
    clim_i = max(abs(h(:)),[],'omitnan'); if isempty(clim_i) || clim_i==0; clim_i = 1; end
    caxis(ax,[-clim_i,clim_i])
    axis(ax,'tight'); yticks(ax,[])
    if trial_is_right(i); side_lbl = 'R'; else; side_lbl = 'L'; end
    title(ax,sprintf('t%d %s (%s)',i,fly_short_label(all_data(i).meta),side_lbl),'FontSize',6,'Interpreter','none')
end
title(tl5,'debug: per-trial EPG calcium (im.d) heatmap aligned to stims, own color scale per panel (side label from the atp-channel evoked-response call)','Interpreter','none')

%% average EPG activity (im.d, summed across glomeruli) aligned to stim onset -- confirmed ejections only
% this dataset's atp/calcium peri-stim traces (stim_im_peri/stim_im_base)
% were already computed on the common t_common axis in the ejection-
% detection loop above. Pooled one value per fly first (that fly's own
% confirmed-ejection stims averaged together), then across flies for the
% plotted mean +/- SEM -- same one-point-per-fly convention used
% throughout this codebase (e.g. lpsp_kir_claude.m, lpsp_p2x2_claude.m),
% so a fly with more confirmed ejections doesn't get more weight than a
% fly with fewer.
ej_idx = find(stim_is_ejection);

fly_im_trace = cell(n_flies,1);
for f = 1:n_flies
    idx = ej_idx(fly_num(stim_trial(ej_idx))==f);
    if isempty(idx); continue; end
    fly_im_trace{f} = mean(cell2mat(stim_im_peri(idx)'),1,'omitnan');
end
have_im = ~cellfun(@isempty,fly_im_trace);
M_im    = cell2mat(fly_im_trace(have_im));
im_mean = mean(M_im,1,'omitnan');
im_sem  = std(M_im,1,'omitnan')/sqrt(size(M_im,1));
fprintf('\naverage EPG trace: pooled from %d/%d flies with >=1 confirmed ejection\n', size(M_im,1), n_flies);

figure(6); clf
set(gcf,'Name','average EPG activity aligned to stim onset','Position',[100,100,600,400])
hold on
patch([t_common,fliplr(t_common)],[im_mean+im_sem,fliplr(im_mean-im_sem)],[0.13,0.55,0.13],'FaceAlpha',.2,'EdgeColor','none')
plot(t_common,im_mean,'Color',[0.13,0.55,0.13],'LineWidth',2)
plot([0,0],ylim,':k')
xlabel('time from stim onset (s)')
ylabel('EPG activity (\Sigma dF/F across glomeruli)')
title(sprintf('average EPG (im.d) activity aligned to ATP stim onset (mean +/- SEM across %d flies, confirmed ejections only)',size(M_im,1)))

%% figure: average glomerulus (im.d) heatmap around confirmed ATP ejections, split by stim side
% same convention as lpsp_p2x2_claude.m's own heatmap section, minus the
% genotype split that script needed (single genotype here) -- but WITH
% the side split that script had, now that trial_is_right (the corrected
% stim-evoked hemisphere metric, computed in the debug section above) is
% unambiguous on this dataset. Each panel pools its own side's
% confirmed-ejection stims, baseline-subtracted per glomerulus row at its
% own t=0, no PB flip applied (raw, physical hemisphere order preserved)
% -- same reasoning as the reference script: the point of keeping side
% separate is seeing where each side's own stim actually landed.
heatmap_pool = cell(1,2); % {1}=right-stim trials, {2}=left-stim trials; each a 32 x length(t_common) x n_stims stack
for jj = 1:numel(ej_idx)
    j = ej_idx(jj);
    i = stim_trial(j);
    xb_rel  = all_data(i).ft.xb - stim_onset_t(j);
    aligned = interp1(xb_rel,all_data(i).im.d',t_common,'linear',nan)'; % 32 x length(t_common)
    aligned = aligned - aligned(:,zero_idx);
    side = 1 + ~trial_is_right(i); % 1 = right-stim trial, 2 = left-stim trial
    heatmap_pool{side} = cat(3,heatmap_pool{side},aligned);
end

mean_heatmap = cell(1,2);
n_heatmap    = zeros(1,2);
for side = 1:2
    n_heatmap(side) = size(heatmap_pool{side},3);
    if n_heatmap(side) > 0
        mean_heatmap{side} = mean(heatmap_pool{side},3,'omitnan');
    end
end
clim_max = max(cellfun(@(m) max(abs(m(:)),[],'omitnan'), mean_heatmap(~cellfun(@isempty,mean_heatmap))));

side_labels = {'right-stim trials','left-stim trials'};

figure(7); clf
set(gcf,'Name','glomerulus dF/F heatmap aligned to stims, by side','Position',[100,100,900,450])
for side = 1:2
    subplot(1,2,side); hold on
    if isempty(mean_heatmap{side}); continue; end
    imagesc(t_common,1:32,mean_heatmap{side})
    plot([t_common(1),t_common(end)],[16.5,16.5],'k:','LineWidth',1) % hemisphere boundary
    plot([0,0],[.5,32.5],'k:','LineWidth',1) % stim onset
    colormap(gca,diverge_cmap)
    caxis([-clim_max,clim_max])
    axis tight
    yticks([8.5,24.5]); yticklabels({'rows 1-16','rows 17-32'})
    xlabel('time from stim onset (s)')
    c = colorbar; c.Label.String = 'dF/F (baseline-subtracted)';
    title(sprintf('%s (n=%d stims)',side_labels{side},n_heatmap(side)))
end
sgtitle('average EPG calcium (im.d) heatmap around confirmed ATP ejections, by stim side (no PB flip applied)')

%% average bump position (im.mu), zeroed at the stim rising edge -- confirmed ejections only
% mu is unwrapped once per trial (im.mu is circular; unwrapping per-trial
% before windowing keeps continuity across the stim onset), then each
% confirmed-ejection stim's own peri-stim window is interpolated onto
% t_bump and shifted so its OWN value at t=0 is exactly 0 -- "bump
% position relative to where it was at the stim rising edge". A wider
% window than the +/-5s atp/calcium windows above is used here (-5s to
% +15s), same as lpsp_p2x2_claude.m, since a bump-position effect is
% expected to play out more slowly than the fast atp/calcium kinetics.
bump_pre = -5; bump_post = 15;
t_bump   = linspace(bump_pre,bump_post,201);

mu_unwrapped  = cell(1,n_trials);  % unwrap(im.mu), computed once per trial, lazily
cue_unwrapped = cell(1,n_trials);  % unwrap(-ft.cue), same sign flip used throughout this codebase so cue is directly comparable to mu
stim_mu_rel   = cell(size(stim_trial)); % filled only for confirmed-ejection stims (ej_idx)
stim_cue_rel  = cell(size(stim_trial));

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

%% average bump position by stim side -- confirmed ejections only
% pooled one value per fly first (that fly's own confirmed-ejection stims
% within one side, averaged together), then across flies for the plotted
% mean +/- SEM, same one-point-per-fly convention as everywhere else in
% this codebase. Two side-specific traces, using trial_is_right from the
% debug section above.
side_labels_bump = {'stim right','stim left'};
side_color       = [0.85,0.33,0.10]; % same for both sides; distinguished by line style below
side_style       = {'-','--'};

bump_mean = nan(2,numel(t_bump)); bump_sem = nan(2,numel(t_bump)); bump_n = zeros(1,2);
cue_mean  = nan(2,numel(t_bump)); cue_sem  = nan(2,numel(t_bump));
for side = 1:2
    [bump_mean(side,:),bump_sem(side,:),bump_n(side)] = pool_side_traces(ej_idx,stim_mu_rel, side==1,fly_num,stim_trial,trial_is_right,n_flies,numel(t_bump));
    [cue_mean(side,:), cue_sem(side,:), ~]             = pool_side_traces(ej_idx,stim_cue_rel,side==1,fly_num,stim_trial,trial_is_right,n_flies,numel(t_bump));
end
fprintf('\n=== average bump position by stim side, zeroed at stim onset (confirmed ejections only) ===\n');
for side = 1:2
    fprintf('  %-12s n=%d flies\n', side_labels_bump{side}, bump_n(side));
end

figure(8); clf
set(gcf,'Name','average bump position aligned to stim onset, by side','Position',[100,100,700,550])
hold on
h_bump = gobjects(1,2); h_cue = gobjects(1,2);
for side = 1:2
    if bump_n(side)==0; continue; end
    patch([t_bump,fliplr(t_bump)],[cue_mean(side,:)+cue_sem(side,:),fliplr(cue_mean(side,:)-cue_sem(side,:))], ...
        [.5,.5,.5],'FaceAlpha',.10,'EdgeColor','none','HandleVisibility','off')
    h_cue(side) = plot(t_bump,cue_mean(side,:),'Color',[.5,.5,.5],'LineStyle',side_style{side},'LineWidth',1.5);
    patch([t_bump,fliplr(t_bump)],[bump_mean(side,:)+bump_sem(side,:),fliplr(bump_mean(side,:)-bump_sem(side,:))], ...
        side_color,'FaceAlpha',.15,'EdgeColor','none','HandleVisibility','off')
    h_bump(side) = plot(t_bump,bump_mean(side,:),'Color',side_color,'LineStyle',side_style{side},'LineWidth',2);
end
plot([t_bump(1),t_bump(end)],[0,0],':k')
plot([0,0],ylim,':k')
xlabel('time from stim onset (s)')
ylabel('unwrapped position relative to stim onset (rad)')
have_side = bump_n > 0;
legend_h  = [h_bump(have_side),h_cue(find(have_side,1))];
legend_labels = [arrayfun(@(s) sprintf('%s (n=%d)',side_labels_bump{s},bump_n(s)),find(have_side),'UniformOutput',false),{'heading cue (gray)'}];
legend(legend_h,legend_labels,'Location','northwest')
title('average unwrapped bump position (color) vs. heading cue (gray), by stim side, zeroed at the ATP stim rising edge (mean +/- SEM across flies, confirmed ejections only)')

%% figure: bump movement collapsed across stim side, by negating mu (and cue) for left-side stims
% since the two PB hemispheres are mirror images of each other, negating
% a left-side stim's own already-zeroed bump trace (stim_mu_rel) before
% pooling puts it on the same sign convention as a right-side stim's
% trace -- so "positive" now consistently means "bump moved in the
% direction expected from a right-side stim" regardless of which side
% actually got stimulated that trial. Heading cue is negated the same
% way. This pools all confirmed-ejection stims from both sides into one
% higher-powered trace instead of figure 8's two side-specific ones.
stim_mu_rel_flip  = stim_mu_rel;
stim_cue_rel_flip = stim_cue_rel;
for jj = 1:numel(ej_idx)
    j = ej_idx(jj);
    if ~trial_is_right(stim_trial(j)) % left-side stim -- flip onto the right-side sign convention
        stim_mu_rel_flip{j}  = -stim_mu_rel{j};
        stim_cue_rel_flip{j} = -stim_cue_rel{j};
    end
end

fly_mu_flip  = cell(n_flies,1);
fly_cue_flip = cell(n_flies,1);
for f = 1:n_flies
    idx = ej_idx(fly_num(stim_trial(ej_idx))==f);
    if isempty(idx); continue; end
    fly_mu_flip{f}  = mean(cell2mat(stim_mu_rel_flip(idx)'), 1,'omitnan');
    fly_cue_flip{f} = mean(cell2mat(stim_cue_rel_flip(idx)'),1,'omitnan');
end
have_flip     = ~cellfun(@isempty,fly_mu_flip);
M_mu_flip     = cell2mat(fly_mu_flip(have_flip));
M_cue_flip    = cell2mat(fly_cue_flip(have_flip));
flip_mean     = mean(M_mu_flip,1,'omitnan');  flip_sem     = std(M_mu_flip,1,'omitnan')/sqrt(size(M_mu_flip,1));
flip_cue_mean = mean(M_cue_flip,1,'omitnan'); flip_cue_sem = std(M_cue_flip,1,'omitnan')/sqrt(size(M_cue_flip,1));
fprintf('\naverage bump position, collapsed across stim side: pooled from %d/%d flies with >=1 confirmed ejection\n', size(M_mu_flip,1), n_flies);

figure(9); clf
set(gcf,'Name','bump movement collapsed across stim side','Position',[100,100,600,450])
hold on
patch([t_bump,fliplr(t_bump)],[flip_cue_mean+flip_cue_sem,fliplr(flip_cue_mean-flip_cue_sem)],[.5,.5,.5],'FaceAlpha',.15,'EdgeColor','none')
h_cue2 = plot(t_bump,flip_cue_mean,'Color',[.5,.5,.5],'LineWidth',1.5);
patch([t_bump,fliplr(t_bump)],[flip_mean+flip_sem,fliplr(flip_mean-flip_sem)],[0.85,0.33,0.10],'FaceAlpha',.15,'EdgeColor','none')
h_bump2 = plot(t_bump,flip_mean,'Color',[0.85,0.33,0.10],'LineWidth',2);
plot([t_bump(1),t_bump(end)],[0,0],':k')
plot([0,0],ylim,':k')
xlabel('time from stim onset (s)')
ylabel('position relative to stim onset (rad), left-side stims sign-flipped')
legend([h_bump2,h_cue2],{sprintf('bump position (n=%d flies)',size(M_mu_flip,1)),'heading cue (gray)'},'Location','northwest')
title({'average unwrapped bump position (color) vs. heading cue (gray), zeroed at stim onset,','collapsed across stim side by negating left-side stims (confirmed ejections only)'})

%% save all figures as PDF
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\dopamine_ionto_redo';
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
function [g_mean,g_sem,g_n] = pool_side_traces(stim_pool,trace_cell,want_right,fly_num,stim_trial,trial_is_right,n_flies,t_len)
    % one value per fly first (that fly's own stims in stim_pool that
    % belong to the requested side, averaged together), then averaged
    % across flies -- see lpsp_p2x2_claude.m's pool_group_traces for the
    % same convention with a genotype dimension this dataset doesn't need.
    fly_trace = cell(n_flies,1);
    for f = 1:n_flies
        idx = stim_pool(fly_num(stim_trial(stim_pool))==f & trial_is_right(stim_trial(stim_pool))==want_right);
        if isempty(idx); continue; end
        fly_trace{f} = mean(cell2mat(trace_cell(idx)'),1,'omitnan');
    end
    have = ~cellfun(@isempty,fly_trace);
    g_mean = nan(1,t_len); g_sem = nan(1,t_len); g_n = sum(have);
    if g_n == 0; return; end
    M = cell2mat(fly_trace(have));
    g_mean = mean(M,1,'omitnan');
    g_sem  = std(M,1,'omitnan')/sqrt(g_n);
end

function lbl = fly_short_label(fly_path)
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N, or
    % from a full trial meta path -- strips down to the last 2 non-empty
    % path parts either way.
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        lbl = strjoin(parts(max(1,end-1):end),' ');
    else
        lbl = strjoin(parts(max(1,fly_part(1)-1):fly_part(1)),' ');
    end
end

function fid = trial_fly_id(meta_path)
    % meta = ...\<date folder>\fly N\<trial folder>\registration
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        fid = meta_path; % shouldn't happen; fall back to a unique id per trial
    else
        fid = strjoin(parts(1:fly_part(1)),filesep); % date folder + "fly N"
    end
end
