%% lpsp_rnai_turns_claude
% Turn-resolved look at the LPsP-RNAi dataset (lpsp_rnai_joint_nosmooth_reg_
% *.mat, same file and same genotype/fly/light-condition labeling as
% lpsp_rnai_claude_v2.m), for all five genotypes (lpsp>th, empty>th,
% lpsp>vglut, empty>vglut, lpsp>mcherry) and both light conditions (closed
% loop and dark, each in its own set of figures). Two questions, each with
% its own block of figures:
%
%   1) TURN-SIZE-RESOLVED BUMP MOTION. Every discrete turn the fly makes
%      is detected from its own rotational velocity (ft.r_speed), sized by
%      how far the fly actually rotated over the turn (ft.heading), and
%      binned by that size (15 deg, 30 deg, 45 deg, ... -- turn_bin_edges_deg
%      below). Within each bin, the fly's heading and the bump's position
%      are aligned to turn onset, sign-flipped so every turn is "positive",
%      and overlaid (mean +/- SEM), one row of subplots per fly. Group
%      overlays come in four variants: all turns vs. time; only turns with
%      a QUIET BASELINE (no counter-correction of a preceding opposite
%      turn); rotational VELOCITY vs. time; and displacement vs.
%      NORMALIZED time (0 = onset, 1 = offset). A summary figure then
%      plots how far the bump moved vs. how far the fly moved per
%      turn-size bin, per fly, by genotype.
%
%   2) GAIN vs. BUMP POSITION IN THE PB. Per fly, the velocity gain
%      (slope of bump velocity vs. the fly's own rotational velocity --
%      same regression, smoothing, lag and sample gating as
%      lpsp_rnai_claude_v2.m step 8 / gain_scratch_claude.m's settled
%      pipeline) is recomputed separately for samples in which the bump
%      sits in each PB glomerulus, and plotted against that position. A
%      turn-based version (per-turn bump displacement vs. fly displacement,
%      binned by where the bump was at turn onset) is shown alongside.
%
% "Bump position in the PB" here is the PVA angle im.mu mapped back onto
% the glomerulus it points at. im.alpha assigns the SAME n_centroid angles
% to the left and right halves of the PB (alpha = repmat(linspace(-pi,pi-
% 2pi/n,n),1,2), see lpsp_rnai.m's dataset-building code), so a PVA angle
% names one glomerulus in each hemisphere, not a unique one of the 2n
% clusters -- the position axis below is therefore "glomerulus within a
% hemisphere" (1..n_centroid), and a diagnostic figure (mean z-scored
% fluorescence across all 2n clusters as a function of that binned PVA
% angle) checks that the mapping really lands on a matching pair of
% clusters.
%
% Sign conventions (confirmed directly on the dataset before writing this):
% ft.heading integrates ft.r_speed (positive = one rotation direction),
% and the bump im.mu moves in the SAME sign as heading in closed loop (the
% closed-loop cue is -0.8 x heading, and mu tracks -cue), so bump and fly
% displacements are directly comparable with no sign flip; the expected
% closed-loop ratio is the rig's 0.8 cue gain, not 1.
%
% Run from the repo root, one %% section at a time (or top to bottom).
% Every figure is exported as a PDF to ugly_figures/rnai_turns/ (the
% per-fly overlay figures into its per_fly/ subfolder).
addpath(fullfile(pwd,'circ_stats'));
set(groot,'defaultAxesToolbarVisible','off'); % otherwise the axes toolbar can get baked into exported PNGs when run in batch mode

%% 0) parameters
data_dir    = fullfile('.data'); % this checkout keeps the .mat datasets in ".data" (see lpsp_rnai_claude_v2.m)
source_file = 'lpsp_rnai_joint_nosmooth_reg_20260831.mat';

geno_keep    = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'}; % which genotypes to analyze (all five in the dataset)
geno_colors  = [0.85,0.33,0.10;   % lpsp>th      orange
                0.00,0.45,0.74;   % empty>th     blue
                0.49,0.18,0.56;   % lpsp>vglut   purple
                0.47,0.67,0.19;   % empty>vglut  green
                0.93,0.69,0.13];  % lpsp>mcherry yellow
include_dark = true; % false = closed-loop trials only; true = also analyze dark trials as a second condition, in their own figures

% bump position estimate: settled velocity-gain pipeline from
% lpsp_rnai_claude_v2.m step 8 (im.f movmean 3 frames -> z-score -> PVA,
% then a 0.2 s Gaussian on the unwrapped mu), reused for BOTH questions
% so the per-turn traces and the per-position gains describe the same mu.
im_f_frames     = 3;    % frames, movmean on im.f before z-score + PVA
bump_smooth_s   = 0.20; % s, Gaussian on unwrapped mu (imaging timebase)
bump_lag_s      = 0.13; % s, calcium/bump signal lags behavior by ~8 fictrac frames (gain_scratch_claude.m); in seconds here since fictrac dt varies across trials

% turn detection (question 1): a turn is a contiguous stretch where the
% smoothed rotational velocity stays above turn_vel_thresh WITH ONE SIGN.
% The velocity smoothing here (0.5 s) is lighter than the 1.0 s the gain
% regression uses -- that window was tuned for a per-frame velocity
% regression, and would blur turn onsets/offsets by ~0.5 s each way,
% smearing the aligned traces below.
turn_vel_smooth_s = 0.50; % s, Gaussian on r_speed for turn detection
turn_vel_thresh   = 0.30; % rad/s, |smoothed r_speed| must exceed this throughout the turn
turn_max_gap_s    = 0.20; % s, sub-threshold dips shorter than this (same sign either side) don't split a turn
turn_min_dur_s    = 0.15; % s, shorter detections are discarded as noise
turn_rho_thresh   = 0.20; % mean bump vector strength over the turn window must exceed this (bump has to be there to be tracked)
turn_win_s        = [-1, 3]; % s relative to turn onset, for the aligned overlays
turn_grid_dt      = 0.02; % s, common time grid the overlays are interpolated onto (fictrac dt differs across trials)
% option 1 ("quiet baseline"): a turn additionally qualifies if |smoothed
% r_speed| stayed below turn_vel_thresh for turn_quiet_pre_s before onset,
% i.e. it was NOT a counter-correction of an immediately preceding turn in
% the other direction (those give the small-turn averages a descending
% pre-onset limb). Shown as its own group-overlay variant.
turn_quiet_pre_s  = 0.50; % s
% option 2 ("normalized time"): each turn's traces are also resampled on a
% grid of normalized time tau, tau=0 at onset and tau=1 at offset, so turns
% of different duration line up at both ends instead of smearing the mean.
turn_norm_tau     = (-0.5:0.02:1.5)';

% turn-size bins (deg of fly rotation over the turn): 15-deg-wide bins
% centered on 15,30,...,90, then coarser bins for the rarer big turns.
turn_bin_edges_deg = [7.5:15:97.5, 135, 180, 360];
min_turns_per_bin  = 3; % a fly x bin needs at least this many turns to contribute a mean trace / summary point

% velocity-gain regression (question 2): identical to lpsp_rnai_claude_v2.m
% step 8 (gain_fly, with intercept), except the fit is done per PB
% position bin.
gain_vel_smooth_s = 1.00; % s, Gaussian on r_speed (fictrac timebase)
gain_vel_thresh   = 0.2; gain_bump_thresh = 10; gain_rho_thresh = 0.2; gain_vel_max = 5; % rad/s gates, as in v2
min_samples_per_pos_bin = 200; % ~3 s of above-threshold turning at 60 Hz -- fewer and a per-bin slope is mostly noise
min_turns_per_pos_bin   = 5;   % turn-based version: per-bin through-origin slope needs at least this many turns
turn_min_size_for_pos_deg = 15; % turn-based version: ignore tiny turns (ratio of two small displacements is noise-dominated)

export_prefix    = 'turns';
export_dir       = fullfile('ugly_figures','rnai_turns'); % every group figure (PDF) goes here
per_fly_dir      = fullfile(export_dir,'per_fly');          % the ~350 per-fly overlay figures (fig2) go in their own subfolder
if ~exist(export_dir,'dir'), mkdir(export_dir); end
if ~exist(per_fly_dir,'dir'), mkdir(per_fly_dir); end
% optional PNG previews (quick to flip through) -- define preview_png_dir
% in the workspace BEFORE running the script to enable, e.g. when
% batch-testing; left empty here so the repo only gets the PDFs.
if ~exist('preview_png_dir','var'), preview_png_dir = ''; end
make_per_fly_figs = true; % one aligned-overlay figure per fly (question 1) -- can be many figures; set false to only make the group summaries

%% 1) load data + labels (same parsing as lpsp_rnai_claude_v2.m)
tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:); clear tmp
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

driver = cell(n_trials,1); target = cell(n_trials,1);
for i = 1:n_trials
    [fly_seg,trial_seg] = meta_fly_and_trial_seg(all_data(i).meta);
    combo = lower(regexprep([fly_seg,trial_seg],'[_\s]',''));
    if contains(combo,'lpsp'),      driver{i} = 'lpsp';
    elseif contains(combo,'empty'), driver{i} = 'empty';
    else, driver{i} = ''; warning('trial %d: could not determine driver from meta "%s"', i, all_data(i).meta);
    end
    if contains(combo,'thrnai'),      target{i} = 'th';
    elseif contains(combo,'vglut'),   target{i} = 'vglut';
    elseif contains(combo,'mcherry'), target{i} = 'mcherry';
    else,                             target{i} = 'th'; % bare "..._rnai" in the th project = implicit TH-RNAi
    end
end
genotype = strcat(driver,'>',target);

is_dark = false(n_trials,1);
for i = 1:n_trials
    is_dark(i) = contains(all_data(i).ft.pattern,'background');
end

fly_id = cell(n_trials,1);
for i = 1:n_trials
    fly_id{i} = meta_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
n_flies = numel(fly_list);

n_centroid = size(all_data(1).im.f,1)/2; % glomeruli per hemisphere (im.alpha repeats over the two halves)
alpha_hemi = all_data(1).im.alpha(1:n_centroid); alpha_hemi = alpha_hemi(:)';
fprintf('%d trials -> %d flies; %d glomeruli per hemisphere\n', n_trials, n_flies, n_centroid);

% timebase check: imaging (ft.xb) vs. fictrac (ft.xf) end times -- see bump_on_fictrac
xb_gap = nan(n_trials,1);
for i = 1:n_trials
    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == numel(all_data(i).im.mu)
        xb_gap(i) = all_data(i).ft.xf(end) - all_data(i).ft.xb(end);
    end
end
fprintf('imaging ends before fictrac by %.1f s (median; range %.1f-%.1f s) in %d/%d trials with a matching ft.xb\n', ...
    median(xb_gap,'omitnan'), min(xb_gap), max(xb_gap), sum(~isnan(xb_gap)), n_trials);

if include_dark, cond_list = [false,true]; else, cond_list = false; end
cond_name = {'closed loop','dark'};

%% 2) per-trial bump position (settled gain pipeline) -- computed once, reused by both questions
use_trials = find(ismember(genotype,geno_keep) & ismember(is_dark,cond_list));
pva_mu  = cell(n_trials,1); pva_rho = cell(n_trials,1); pva_fz = cell(n_trials,1);
for ti = use_trials'
    [pva_mu{ti},pva_rho{ti},pva_fz{ti}] = trial_pva_movmean(all_data(ti),im_f_frames);
end
fprintf('bump position estimated for %d trials (%s)\n', numel(use_trials), strjoin(geno_keep,', '));

%% 3) detect turns in every used trial, extract aligned traces + per-turn displacements
% one struct array entry per turn, with: trial index, fly, genotype, dark,
% size (deg), direction, aligned heading/bump traces on the common grid
% (deg, sign-flipped so the turn is positive), bump displacement over the
% turn (deg, lagged by bump_lag_s), and the bump's PB position at onset.
turn_t = (turn_win_s(1):turn_grid_dt:turn_win_s(2))';
turns = [];
for ti = use_trials'
    tr = trial_turns(all_data(ti),pva_mu{ti},pva_rho{ti},bump_smooth_s,bump_lag_s, ...
        turn_vel_smooth_s,turn_vel_thresh,turn_max_gap_s,turn_min_dur_s,turn_rho_thresh,turn_win_s,turn_t,turn_quiet_pre_s,turn_norm_tau);
    if isempty(tr), continue, end
    [tr.trial] = deal(ti); [tr.fly] = deal(fly_num(ti)); [tr.geno] = deal(genotype{ti}); [tr.dark] = deal(is_dark(ti));
    turns = [turns; tr(:)]; %#ok<AGROW>
end
if isempty(turns)
    error('no turns detected in any trial -- check turn_vel_thresh / turn_rho_thresh');
end
turn_size  = [turns.size_deg]';
turn_bin   = discretize(turn_size,turn_bin_edges_deg); % NaN = below smallest edge or above largest
n_bins     = numel(turn_bin_edges_deg)-1;
bin_labels = cell(1,n_bins);
for b = 1:n_bins
    if diff(turn_bin_edges_deg(b:b+1)) <= 15.01
        bin_labels{b} = sprintf('%g%s',mean(turn_bin_edges_deg(b:b+1)),char(176)); % 15-deg-wide bin: label by its center ("15deg", "30deg", ...)
    else
        bin_labels{b} = sprintf('%g-%g%s',turn_bin_edges_deg(b),turn_bin_edges_deg(b+1),char(176));
    end
end
fprintf('%d turns detected (%d within the size bins); per bin: %s\n', numel(turns), sum(~isnan(turn_bin)), ...
    strjoin(arrayfun(@(b) sprintf('%s=%d',bin_labels{b},sum(turn_bin==b)),1:n_bins,'UniformOutput',false),', '));
turn_quiet = [turns.quiet_pre]';
fprintf('%d of %d binned turns have a quiet baseline (|v| < %.2f rad/s for %.2fs before onset)\n', ...
    sum(turn_quiet & ~isnan(turn_bin)), sum(~isnan(turn_bin)), turn_vel_thresh, turn_quiet_pre_s);

%% figure: how turns are defined -- example trial with the most turns
turn_trials = [turns.trial]';
[~,ex_ti] = max(accumarray(turn_trials,1,[n_trials,1]));
trial = all_data(ex_ti); xf = trial.ft.xf(:); dt = median(diff(xf));
r_sm = smoothdata(trial.ft.r_speed(:),'gaussian',max(1,round(turn_vel_smooth_s/dt)));
heading_deg = rad2deg(trial_heading(trial));
[mu_f,~] = bump_on_fictrac(trial,pva_mu{ex_ti},pva_rho{ex_ti},bump_smooth_s);
these = turns(turn_trials==ex_ti);

figure(1); clf; set(gcf,'Name','turn definition, example trial','Position',[80,80,1300,600],'Color','w')
ax1 = subplot(2,1,1); hold on
plot(xf,r_sm,'k'); yline(turn_vel_thresh,':r'); yline(-turn_vel_thresh,':r'); yline(0,'-','Color',[.8,.8,.8])
for k = 1:numel(these)
    xs = these(k).onset_s + [0,these(k).dur_s];
    patch([xs,fliplr(xs)],[-5,-5,5,5],[1,.85,.6],'EdgeColor','none','FaceAlpha',.5)
end
ylim([-1,1]*min(5,max(abs(r_sm),[],'omitnan')*1.1)); ylabel(sprintf('r\\_speed, %.2fs Gaussian (rad/s)',turn_vel_smooth_s))
title(sprintf('trial %d (%s, %s): %d turns; shaded = detected turns (|v|>%.2f rad/s, same sign, >%.2fs)', ...
    ex_ti,genotype{ex_ti},cond_name{is_dark(ex_ti)+1},numel(these),turn_vel_thresh,turn_min_dur_s),'Interpreter','none')
ax2 = subplot(2,1,2); hold on
plot(xf,heading_deg-heading_deg(1),'k','DisplayName','fly heading')
mu0 = mu_f(find(~isnan(mu_f),1)); % first valid sample (fictrac can start before the first imaging frame)
plot(xf,rad2deg(mu_f-mu0),'Color',geno_colors(1,:),'DisplayName','bump (unwrapped PVA)')
xlabel('time (s)'); ylabel('cumulative rotation (deg)'); legend('Location','best')
linkaxes([ax1,ax2],'x'); xlim([xf(1),min(xf(end),xf(1)+120)])
export_fig(export_dir,preview_png_dir,sprintf('%s_fig1_turn_definition_example',export_prefix))

%% 4) QUESTION 1 -- per-fly aligned overlays: fly heading vs. bump, one subplot per turn-size bin
fly_geno = cell(n_flies,1);
for f = 1:n_flies
    g = unique(genotype(fly_num==f)); fly_geno{f} = g{1};
end
used_flies = unique([turns.fly]');

% per fly x condition x bin: mean aligned traces (fly, bump), n turns,
% and summary displacements -- stored for the group figures below
n_turns_fb      = zeros(n_flies,numel(cond_list),n_bins);
dfly_fb         = nan(n_flies,numel(cond_list),n_bins); % mean fly displacement over the turn (deg)
dbump_fb        = nan(n_flies,numel(cond_list),n_bins); % mean bump displacement over the turn (deg, lagged)
ratio_fb        = nan(n_flies,numel(cond_list),n_bins); % through-origin slope dbump ~ dfly across this fly's turns in the bin

for f = used_flies'
    gi = find(strcmp(geno_keep,fly_geno{f}));
    for ci = 1:numel(cond_list)
        sel_fc = [turns.fly]'==f & [turns.dark]'==cond_list(ci);
        if ~any(sel_fc), continue, end
        if make_per_fly_figs
            fig = figure(100+f); clf
            set(fig,'Name',sprintf('turn overlays, fly %d',f),'Position',[60,60,220*n_bins,330],'Color','w')
        end
        for b = 1:n_bins
            sel = find(sel_fc & turn_bin==b);
            n_turns_fb(f,ci,b) = numel(sel);
            if numel(sel) < min_turns_per_bin, continue, end
            FT = cat(2,turns(sel).fly_trace)';  % turns x time
            BT = cat(2,turns(sel).bump_trace)';
            dfly  = [turns(sel).dfly_deg]'; dbump = [turns(sel).dbump_deg]';
            dfly_fb(f,ci,b)  = mean(dfly,'omitnan');
            dbump_fb(f,ci,b) = mean(dbump,'omitnan');
            ratio_fb(f,ci,b) = dfly \ dbump;
            if make_per_fly_figs
                subplot(1,n_bins,b); hold on
                plot_sem(turn_t,FT,[0,0,0])
                plot_sem(turn_t,BT,geno_colors(gi,:))
                plot(turn_t,0.8*mean(FT,1,'omitnan'),'--','Color',[.5,.5,.5])
                xline(0,':k')
                title(sprintf('%s turns (n=%d)',bin_labels{b},numel(sel)))
                xlabel('time from turn onset (s)'); xlim(turn_win_s)
                if b == 1, ylabel('rotation (deg)'), end
            end
        end
        if make_per_fly_figs
            sgtitle(sprintf('%s, %s -- fly %d (%s): fly heading (black) vs. bump (color), mean +/- SEM; dashed = 0.8 x fly', ...
                fly_geno{f},cond_name{cond_list(ci)+1},f,fly_short_label(fly_list{f})),'Interpreter','none')
            export_fig(per_fly_dir,preview_png_dir,sprintf('%s_fig2_overlay_%s_%s_fly%02d',export_prefix,strrep(fly_geno{f},'>','-'),strrep(cond_name{cond_list(ci)+1},' ',''),f))
            close(fig)
        end
    end
end

%% figures: group-level aligned overlays (mean of per-fly means), one row per genotype -- four variants
% (a) all turns, displacement vs. time from onset (the original view)
% (b) OPTION 1: only turns with a quiet baseline (|v| below threshold for
%     turn_quiet_pre_s before onset) -- drops the counter-corrections that
%     give small-turn averages their descending pre-onset limb
% (c) OPTION 2a: all turns, rotational VELOCITY vs. time from onset (fly:
%     r_speed with the turn-detection smoothing; bump: d/dt of the unwrapped
%     PVA on the fictrac timebase, same smoothing, so the two are comparable)
% (d) OPTION 2b: all turns, displacement vs. NORMALIZED time (0 = onset,
%     1 = offset) so turns of different duration line up at both ends
turn_fly  = [turns.fly]';
turn_dark = [turns.dark]';
variants = struct( ...
    'label',  {'all turns, displacement', sprintf('quiet-baseline turns (%.2fs), displacement',turn_quiet_pre_s), 'all turns, rotational velocity', 'all turns, displacement vs. normalized time (0 = onset, 1 = offset)'}, ...
    'ffield', {'fly_trace','fly_trace','fly_vel_trace','fly_norm'}, ...
    'bfield', {'bump_trace','bump_trace','bump_vel_trace','bump_norm'}, ...
    't',      {turn_t, turn_t, turn_t, turn_norm_tau}, ...
    'quiet',  {false, true, false, false}, ...
    'xlab',   {'time from turn onset (s)','time from turn onset (s)','time from turn onset (s)','normalized time'}, ...
    'ylab',   {'rotation (deg)','rotation (deg)','rotational velocity (deg/s)','rotation (deg)'}, ...
    'suffix', {'fig3a_group_overlay','fig3b_group_overlay_quiet_baseline','fig3c_group_overlay_velocity','fig3d_group_overlay_normalized_time'});
for ci = 1:numel(cond_list)
    for v = 1:numel(variants)
        sel = turn_dark==cond_list(ci);
        if variants(v).quiet, sel = sel & turn_quiet; end
        FM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,variants(v).ffield,min_turns_per_bin);
        BM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,variants(v).bfield,min_turns_per_bin);
        figure(1000+10*ci+v); clf
        set(gcf,'Name',sprintf('group turn overlays (%s), %s',variants(v).label,cond_name{cond_list(ci)+1}),'Position',[60,60,220*n_bins,240*numel(geno_keep)],'Color','w')
        plot_group_overlay(FM,BM,variants(v).t,geno_keep,fly_geno,geno_colors,bin_labels,variants(v).xlab,variants(v).ylab,strcmp(variants(v).ffield,'fly_norm'))
        sgtitle(sprintf('%s, %s: fly heading (black) vs. bump (color), mean +/- SEM across flies (each fly = mean of its >= %d turns per bin); dashed = 0.8 x fly', ...
            cond_name{cond_list(ci)+1},variants(v).label,min_turns_per_bin),'Interpreter','none')
        export_fig(export_dir,preview_png_dir,sprintf('%s_%s_%s',export_prefix,variants(v).suffix,strrep(cond_name{cond_list(ci)+1},' ','')))
    end
end

%% figure: summary -- bump displacement vs. fly displacement per turn-size bin, per fly, by genotype
for ci = 1:numel(cond_list)
    figure(20+ci); clf
    set(gcf,'Name',sprintf('turn-size summary, %s',cond_name{cond_list(ci)+1}),'Position',[60,60,1500,500],'Color','w')
    subplot(1,3,1); hold on
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        X = squeeze(dfly_fb(flies_g,ci,:)); Y = squeeze(dbump_fb(flies_g,ci,:));
        if isvector(X), X = X(:)'; Y = Y(:)'; end
        plot(X',Y','-','Color',[geno_colors(gi,:),.25],'Marker','.','MarkerSize',8,'HandleVisibility','off')
        errorbar(mean(X,1,'omitnan'),mean(Y,1,'omitnan'),sem(Y,1),sem(Y,1),sem(X,1),sem(X,1),'o-','Color',geno_colors(gi,:), ...
            'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2,'DisplayName',sprintf('%s (n=%d flies)',geno_keep{gi},sum(any(~isnan(X),2))))
    end
    xl = [0,max([turn_bin_edges_deg(end-1); dfly_fb(:)],[],'omitnan')*1.05];
    plot(xl,xl,':k','DisplayName','1:1'); plot(xl,0.8*xl,'--','Color',[.5,.5,.5],'DisplayName','0.8:1 (cue gain)')
    xlabel('fly rotation over turn (deg)'); ylabel(sprintf('bump rotation over turn (deg, lag %.2fs)',bump_lag_s))
    legend('Location','southeast','Interpreter','none'); title('mean displacement per turn-size bin')
    axis square; xlim(xl); ylim([max(-50,min([-10; dbump_fb(:)],[],'omitnan')),xl(2)]) % floor clamped at -50: a lone fly x bin with a big negative mean shouldn't flatten the plot

    subplot(1,3,2); hold on
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        R = squeeze(ratio_fb(flies_g,ci,:)); if isvector(R), R = R(:)'; end
        xj = (1:n_bins) + (gi-(numel(geno_keep)+1)/2)*(0.8/numel(geno_keep));
        plot(repmat(xj,size(R,1),1)',R','.','Color',[geno_colors(gi,:),.35],'MarkerSize',10,'HandleVisibility','off')
        errorbar(xj,mean(R,1,'omitnan'),sem(R,1),'o-','Color',geno_colors(gi,:),'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2)
    end
    yline(0.8,'--','Color',[.5,.5,.5]); yline(1,':k'); yline(0,'-','Color',[.85,.85,.85])
    xticks(1:n_bins); xticklabels(bin_labels); xlim([0.5,n_bins+0.5])
    xlabel('turn size (fly rotation)'); ylabel('bump / fly displacement (per-fly through-origin slope)')
    title('gain per turn-size bin')

    subplot(1,3,3); hold on
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        N = squeeze(n_turns_fb(flies_g,ci,:)); if isvector(N), N = N(:)'; end
        xj = (1:n_bins) + (gi-(numel(geno_keep)+1)/2)*(0.8/numel(geno_keep));
        bar(xj,mean(N,1),0.7/numel(geno_keep),'FaceColor',geno_colors(gi,:),'EdgeColor','none','FaceAlpha',.7)
    end
    xticks(1:n_bins); xticklabels(bin_labels); xlim([0.5,n_bins+0.5])
    xlabel('turn size'); ylabel('mean turns per fly'); title('turn counts')
    sgtitle(sprintf('%s: turn-size-resolved bump motion, %s',cond_name{cond_list(ci)+1},strjoin(geno_keep,', ')),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig4_turn_size_summary_%s',export_prefix,strrep(cond_name{cond_list(ci)+1},' ','')))
end

%% 5) QUESTION 2 -- velocity gain as a function of the bump's PB position, per fly
% For each fly x condition, pool every trial's fictrac-timebase vectors
% (fly_vel, bump_vel, rho, wrapped mu -- all lagged/trimmed identically),
% keep the gain-gated samples, bin them by which glomerulus the bump's PVA
% angle points at, and fit bump_vel ~ 1 + fly_vel separately per bin.
pos_edges = [alpha_hemi - pi/n_centroid, alpha_hemi(end) + pi/n_centroid]; % bins centered on each glomerulus's PVA angle
gain_pos      = nan(n_flies,numel(cond_list),n_centroid); % per-bin gain
gain_all      = nan(n_flies,numel(cond_list));            % this fly's overall gain (all bins pooled), same gates
n_samp_pos    = zeros(n_flies,numel(cond_list),n_centroid);
occupancy_pos = nan(n_flies,numel(cond_list),n_centroid); % fraction of gated samples the bump spent in each bin
fz_by_pos     = nan(n_flies,numel(cond_list),2*n_centroid,n_centroid); % mean z-scored fluorescence of all 2n clusters, per PVA bin (mapping check)

for f = used_flies'
    for ci = 1:numel(cond_list)
        trial_list = find(fly_num==f & is_dark==cond_list(ci) & ismember(genotype,geno_keep));
        if isempty(trial_list), continue, end
        fly_vel = []; bump_vel = []; mu_w = []; valid = logical([]); fz_acc = zeros(2*n_centroid,n_centroid); fz_n = zeros(1,n_centroid);
        for ti = trial_list'
            [fv,bv,mw,vd,fz_bin] = trial_gain_by_position_vectors(all_data(ti),pva_mu{ti},pva_rho{ti},pva_fz{ti}, ...
                bump_smooth_s,gain_vel_smooth_s,bump_lag_s,gain_vel_thresh,gain_bump_thresh,gain_rho_thresh,gain_vel_max,pos_edges);
            fly_vel = [fly_vel; fv]; bump_vel = [bump_vel; bv]; mu_w = [mu_w; mw]; valid = [valid; vd]; %#ok<AGROW>
            fz_acc = fz_acc + fz_bin.sum; fz_n = fz_n + fz_bin.n;
        end
        fz_by_pos(f,ci,:,:) = fz_acc ./ fz_n;
        if sum(valid) < min_samples_per_pos_bin, continue, end
        b_all = [ones(sum(valid),1),fly_vel(valid)] \ bump_vel(valid);
        gain_all(f,ci) = b_all(2);
        pos_bin = discretize(wrap_to_range(mu_w,pos_edges(1)),pos_edges);
        for p = 1:n_centroid
            sel = valid & pos_bin==p;
            n_samp_pos(f,ci,p) = sum(sel);
            occupancy_pos(f,ci,p) = sum(sel)/sum(valid);
            if sum(sel) < min_samples_per_pos_bin || std(fly_vel(sel)) < 0.01, continue, end
            bb = [ones(sum(sel),1),fly_vel(sel)] \ bump_vel(sel);
            gain_pos(f,ci,p) = bb(2);
        end
    end
end

% turn-based version: per-turn displacement ratio, binned by bump position at turn onset
turn_gain_pos = nan(n_flies,numel(cond_list),n_centroid);
turn_mu_onset = [turns.mu_onset]';
turn_pos_bin  = discretize(wrap_to_range(turn_mu_onset,pos_edges(1)),pos_edges);
for f = used_flies'
    for ci = 1:numel(cond_list)
        for p = 1:n_centroid
            sel = [turns.fly]'==f & [turns.dark]'==cond_list(ci) & turn_pos_bin==p & turn_size >= turn_min_size_for_pos_deg;
            if sum(sel) < min_turns_per_pos_bin, continue, end
            turn_gain_pos(f,ci,p) = [turns(sel).dfly_deg]' \ [turns(sel).dbump_deg]';
        end
    end
end

%% figure: gain vs. PB position -- per fly (faint) + group mean (thick), by genotype
pos_x = 1:n_centroid;
for ci = 1:numel(cond_list)
    figure(30+ci); clf
    set(gcf,'Name',sprintf('gain vs. PB position, %s',cond_name{cond_list(ci)+1}),'Position',[60,60,1400,800],'Color','w')
    panels = {gain_pos,      'velocity gain (slope bump\_vel ~ 1 + fly\_vel)', 'velocity gain, per PB position (y clipped to [-1 3])';
              gain_pos ./ gain_all, 'gain / fly''s overall gain',          'velocity gain, normalized to each fly''s overall gain (y clipped to [-1 3])';
              turn_gain_pos, 'per-turn bump/fly displacement slope',           sprintf('turn-based gain (turns >= %g%s), by bump position at turn onset (y clipped to [-1 3])',turn_min_size_for_pos_deg,char(176));
              occupancy_pos, 'fraction of gated samples',                     'bump occupancy across PB positions (gated samples)'};
    for k = 1:size(panels,1)
        subplot(2,2,k); hold on
        V = panels{k,1};
        for gi = 1:numel(geno_keep)
            flies_g = find(strcmp(fly_geno,geno_keep{gi}));
            Y = squeeze(V(flies_g,ci,:)); if isvector(Y), Y = Y(:)'; end
            xj = pos_x + (gi-(numel(geno_keep)+1)/2)*(0.7/numel(geno_keep));
            plot(repmat(xj,size(Y,1),1)',Y','-','Color',[geno_colors(gi,:),.2],'Marker','.','MarkerSize',7,'HandleVisibility','off')
            errorbar(xj,mean(Y,1,'omitnan'),sem(Y,1),'o-','Color',geno_colors(gi,:),'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2, ...
                'DisplayName',sprintf('%s (n=%d flies)',geno_keep{gi},sum(any(~isnan(Y),2))))
        end
        if k <= 3
            yline(0.8,'--','Color',[.5,.5,.5],'HandleVisibility','off'); yline(1,':k','HandleVisibility','off'); yline(0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
            if k == 2, yline(1,'--','Color',[.5,.5,.5],'HandleVisibility','off'), end
            ylim([-1,3]) % a few flies have wild per-bin slopes (tiny fly_vel variance in that bin); clipped so the group trend stays readable
        else
            yline(1/n_centroid,':k')
        end
        xticks(pos_x); xticklabels(arrayfun(@(a) sprintf('%d (%.0f%s)',a,rad2deg(alpha_hemi(a)),char(176)),pos_x,'UniformOutput',false)); xtickangle(30)
        xlim([0.5,n_centroid+0.5]); xlabel('PB glomerulus within hemisphere (PVA angle)'); ylabel(panels{k,2})
        title(panels{k,3},'Interpreter','none')
        if k == 1, legend('Location','best','Interpreter','none'), end
    end
    sgtitle(sprintf('%s: gain (fly rotation -> bump rotation) vs. bump position in the PB; faint = flies, thick = mean +/- SEM', ...
        cond_name{cond_list(ci)+1}),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig5_gain_vs_pb_position_%s',export_prefix,strrep(cond_name{cond_list(ci)+1},' ','')))
end

%% figure: polar view of the same per-position gain (bump position is circular)
for ci = 1:numel(cond_list)
    figure(40+ci); clf
    set(gcf,'Name',sprintf('gain vs. PB position (polar), %s',cond_name{cond_list(ci)+1}),'Position',[60,60,380*numel(geno_keep),450],'Color','w')
    for gi = 1:numel(geno_keep)
        % polarplot needs a PolarAxes -- replace the cartesian subplot axes with one at the same position
        ax = subplot(1,numel(geno_keep),gi); pos = ax.Position; delete(ax);
        pax = polaraxes('Position',pos); hold(pax,'on')
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        Y = squeeze(gain_pos(flies_g,ci,:)); if isvector(Y), Y = Y(:)'; end
        th = [alpha_hemi, alpha_hemi(1)];
        for ff = 1:size(Y,1)
            polarplot(pax,th,max([Y(ff,:),Y(ff,1)],0),'-','Color',[geno_colors(gi,:),.2]);
        end
        m = mean(Y,1,'omitnan');
        polarplot(pax,th,max([m,m(1)],0),'-o','Color',geno_colors(gi,:),'LineWidth',2.5,'MarkerFaceColor',geno_colors(gi,:))
        polarplot(pax,linspace(-pi,pi,100),0.8*ones(1,100),'--','Color',[.5,.5,.5])
        rlim(pax,[0,max([1.5,min(3,max(m,[],'omitnan')*1.2)])])
        title(pax,sprintf('%s (n=%d flies)',geno_keep{gi},sum(any(~isnan(Y),2))),'Interpreter','none')
    end
    sgtitle(sprintf('%s: velocity gain by bump PVA angle (radius = gain, dashed = 0.8; negative gains clipped to 0)',cond_name{cond_list(ci)+1}),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig6_gain_vs_pb_position_polar_%s',export_prefix,strrep(cond_name{cond_list(ci)+1},' ','')))
end

%% figure: mapping check -- which of the 2n clusters light up for each PVA position bin
% rows = all 2*n_centroid clusters along the PB (left half then right
% half), columns = the PVA position bin. A correct mapping shows one
% bright cluster per hemisphere per column, on the diagonal in both halves.
for ci = 1:numel(cond_list)
    figure(50+ci); clf
    set(gcf,'Name',sprintf('PVA angle -> glomerulus mapping check, %s',cond_name{cond_list(ci)+1}),'Position',[60,60,380*numel(geno_keep),450],'Color','w')
    for gi = 1:numel(geno_keep)
        subplot(1,numel(geno_keep),gi)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        M = squeeze(mean(fz_by_pos(flies_g,ci,:,:),1,'omitnan'));
        imagesc(1:n_centroid,1:2*n_centroid,M); axis xy
        yline(n_centroid+0.5,'w','LineWidth',1.5)
        xlabel('PVA position bin'); ylabel(sprintf('PB cluster (1-%d left, %d-%d right)',n_centroid,n_centroid+1,2*n_centroid))
        xticks(2:2:n_centroid); yticks(4:4:2*n_centroid)
        cb = colorbar; ylabel(cb,'mean z-scored fluorescence')
        n_f = sum(any(~isnan(fz_by_pos(flies_g,ci,:,:)),[3,4]));
        title(sprintf('%s (n=%d flies)',geno_keep{gi},n_f),'Interpreter','none')
    end
    sgtitle(sprintf('%s -- mapping check: z-scored fluorescence of every PB cluster vs. the PVA position bin (glomerulus within hemisphere) it is assigned to',cond_name{cond_list(ci)+1}),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig7_pva_to_glomerulus_mapping_%s',export_prefix,strrep(cond_name{cond_list(ci)+1},' ','')))
end

%% report
for ci = 1:numel(cond_list)
    fprintf('\n=== per-fly summary (%s) ===\n',cond_name{cond_list(ci)+1});
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        fprintf('%s:\n',geno_keep{gi});
        for f = flies_g'
            if ~any(used_flies==f) || sum(n_turns_fb(f,ci,:)) == 0, continue, end
            fprintf('  fly %3d %-28s turns=%4d  overall gain=%5.2f  gain range across PB positions=[%5.2f %5.2f]\n', ...
                f,fly_short_label(fly_list{f}),sum(n_turns_fb(f,ci,:)),gain_all(f,ci),min(gain_pos(f,ci,:)),max(gain_pos(f,ci,:)));
        end
    end
end

%% functions
function [fly_seg,trial_seg] = meta_fly_and_trial_seg(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        trial_seg = parts{end-1}; fly_seg = parts{end-2};
    else
        trial_seg = parts{end};   fly_seg = parts{end-1};
    end
end

function fid = meta_fly_id(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        fly_part_end = numel(parts)-2;
    else
        fly_part_end = numel(parts)-1;
    end
    fid = strjoin(parts(1:fly_part_end),'\');
end

function lbl = fly_short_label(fly_path)
    parts = strsplit(fly_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function [mu_new,rho_new,f_z] = trial_pva_movmean(trial,smooth_frames)
    % same as lpsp_rnai_claude_v2.m: movmean over frames on raw im.f,
    % z-score per cluster, population vector against im.alpha.
    f_smooth = movmean(trial.im.f,smooth_frames,2);
    f_z = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);
    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp] = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function h = trial_heading(trial)
    % the fly's own cumulative rotation (rad, unwrapped). ft.heading is
    % essentially cumtrapz(r_speed) in this dataset (gain_scratch_claude.m);
    % fall back to integrating r_speed if a trial lacks it.
    if isfield(trial.ft,'heading') && numel(trial.ft.heading) == numel(trial.ft.xf)
        h = trial.ft.heading(:);
        h = fillmissing(h,'linear','EndValues','nearest');
    else
        r = trial.ft.r_speed(:); r(isnan(r)) = 0;
        h = cumtrapz(trial.ft.xf(:),r);
    end
end

function [mu_f,rho_f,xb] = bump_on_fictrac(trial,mu_raw,rho_raw,mu_smooth_s)
    % unwrapped, Gaussian-smoothed bump position + vector strength on the
    % fictrac timebase. Uses the dataset's own imaging timestamps (ft.xb)
    % rather than spreading the frames evenly over xf (as
    % lpsp_rnai_claude_v2.m does): in most trials the two end together,
    % but in a minority imaging stops up to ~15 s before fictrac (checked
    % directly on all 559 trials, see the timebase report in step 1), and
    % stretching would then drift the bump late relative to behavior by
    % that much at trial end. Fictrac samples after the last imaging frame
    % get NaN (no extrapolation) and are dropped downstream.
    xf = trial.ft.xf(:);
    n_im = numel(mu_raw);
    if isfield(trial.ft,'xb') && numel(trial.ft.xb) == n_im
        xb = trial.ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end
    dt_im = median(diff(xb));
    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_s = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    mu_f  = interp1(xb,mu_s,xf,'linear');
    rho_f = interp1(xb,rho_raw(:),xf,'linear');
end

function turns = trial_turns(trial,mu_raw,rho_raw,mu_smooth_s,lag_s,vel_smooth_s,vel_thresh,max_gap_s,min_dur_s,rho_thresh,win_s,t_grid,quiet_pre_s,tau_grid)
    % detect single-direction turns from smoothed r_speed and extract, per
    % turn: size (deg of fly rotation from onset to offset), aligned
    % heading/bump traces on t_grid (deg, relative to onset value,
    % sign-flipped so the turn is positive), the same on a normalized-time
    % grid tau_grid (0 = onset, 1 = offset), rotational-velocity traces on
    % t_grid (deg/s, both signals smoothed identically), whether the
    % quiet_pre_s before onset was sub-threshold ("quiet baseline"), and
    % the bump's displacement over the turn measured lag_s later than the
    % fly's (bump lags behavior). Traces are stored as single: with all
    % five genotypes and both light conditions there are ~10^5 turns x 6
    % traces, and double would be ~1 GB on top of the 3.7 GB dataset.
    xf = trial.ft.xf(:); dt = median(diff(xf));
    heading = trial_heading(trial);
    [mu_f,rho_f] = bump_on_fictrac(trial,mu_raw,rho_raw,mu_smooth_s);

    win = max(1,round(vel_smooth_s/dt));
    r_sm = smoothdata(trial.ft.r_speed(:),'gaussian',win);
    r_sm(isnan(r_sm)) = 0;
    mu_vel = smoothdata(gradient(mu_f)/dt,'gaussian',win); % bump velocity, same smoothing as the fly's -> comparable velocity overlays
    max_gap = max(1,round(max_gap_s/dt)); min_dur = max(1,round(min_dur_s/dt)); n_quiet = max(1,round(quiet_pre_s/dt));

    turns = struct('onset_s',{},'dur_s',{},'size_deg',{},'dir',{},'quiet_pre',{},'fly_trace',{},'bump_trace',{}, ...
                   'fly_vel_trace',{},'bump_vel_trace',{},'fly_norm',{},'bump_norm',{},'dfly_deg',{},'dbump_deg',{},'mu_onset',{},'rho_mean',{});
    for s = [1,-1]
        active = s*r_sm > vel_thresh;
        active = imclose(active,ones(max_gap+1,1)); % bridge brief dips, same sign only (opposite-sign samples are never "active" here)
        active = bwareaopen(active,min_dur);
        d = diff([false;active;false]);
        on = find(d==1); off = find(d==-1)-1;
        for k = 1:numel(on)
            t_on = xf(on(k)); t_off = xf(off(k));
            % need the whole overlay window and the lagged bump offset inside the trial
            if t_on + win_s(1) < xf(1) || max(t_on + win_s(2), t_off + lag_s) > xf(end), continue, end
            rng = on(k):off(k);
            if mean(rho_f(rng),'omitnan') < rho_thresh, continue, end
            dfly = heading(off(k)) - heading(on(k));
            if abs(dfly) < 1e-6, continue, end
            % onset value of the bump taken lag_s after fly onset too, so
            % displacement is measured over the same interval, shifted by the lag
            mu_on  = interp1(xf,mu_f,t_on + lag_s);
            mu_off = interp1(xf,mu_f,t_off + lag_s);
            tq = t_on + t_grid;
            bump_tr = interp1(xf,mu_f,tq);
            if isnan(mu_on) || isnan(mu_off) || any(isnan(bump_tr)), continue, end % turn runs past the last imaging frame
            tau_q = t_on + tau_grid*(t_off - t_on); % normalized-time sample points; NaN where they fall outside the trial/imaging
            mu_at_on = interp1(xf,mu_f,t_on);
            T.onset_s   = t_on;
            T.dur_s     = t_off - t_on;
            T.size_deg  = rad2deg(abs(dfly));
            T.dir       = s;
            T.quiet_pre = on(k) > n_quiet && all(abs(r_sm(on(k)-n_quiet:on(k)-1)) < vel_thresh);
            T.fly_trace  = single(s*rad2deg(interp1(xf,heading,tq) - heading(on(k))));
            T.bump_trace = single(s*rad2deg(bump_tr - mu_at_on));
            T.fly_vel_trace  = single(s*rad2deg(interp1(xf,r_sm,tq)));
            T.bump_vel_trace = single(s*rad2deg(interp1(xf,mu_vel,tq)));
            T.fly_norm  = single(s*rad2deg(interp1(xf,heading,tau_q) - heading(on(k))));
            T.bump_norm = single(s*rad2deg(interp1(xf,mu_f,tau_q) - mu_at_on));
            T.dfly_deg  = s*rad2deg(dfly);          % = +size_deg by construction
            T.dbump_deg = s*rad2deg(mu_off - mu_on);
            T.mu_onset  = wrap_pi(mu_at_on); % where the bump sat (PVA angle) when the turn began
            T.rho_mean  = mean(rho_f(rng),'omitnan');
            turns(end+1) = T; %#ok<AGROW>
        end
    end
    [~,order] = sort([turns.onset_s]); turns = turns(order);
end

function M = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,field,min_turns)
    % per fly x turn-size bin: mean of one aligned-trace field over that
    % fly's selected turns in the bin (NaN if fewer than min_turns);
    % n_flies x n_bins x n_timepoints
    n_t = numel(turns(1).(field));
    M = nan(n_flies,n_bins,n_t);
    for f = unique(turn_fly(sel))'
        for b = 1:n_bins
            idx = find(sel & turn_fly==f & turn_bin==b);
            if numel(idx) < min_turns, continue, end
            X = cat(2,turns(idx).(field))';
            M(f,b,:) = mean(X,1,'omitnan');
        end
    end
end

function plot_group_overlay(FM,BM,t,geno_keep,fly_geno,geno_colors,bin_labels,xlab,ylab,mark_offset)
    % one row per genotype, one column per turn-size bin: mean +/- SEM
    % across flies of the per-fly mean traces in FM (fly, black) and BM
    % (bump, genotype color), plus 0.8 x fly (dashed) and onset (dotted;
    % offset too when mark_offset, for the normalized-time variant)
    n_bins = size(FM,2);
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        for b = 1:n_bins
            subplot(numel(geno_keep),n_bins,(gi-1)*n_bins+b); hold on
            FT = reshape(FM(flies_g,b,:),numel(flies_g),[]);
            BT = reshape(BM(flies_g,b,:),numel(flies_g),[]);
            n_f = sum(~all(isnan(FT),2));
            if n_f >= 1
                plot_sem(t,FT,[0,0,0])
                plot_sem(t,BT,geno_colors(gi,:))
                plot(t,0.8*mean(FT,1,'omitnan'),'--','Color',[.5,.5,.5])
            end
            xline(0,':k'); if mark_offset, xline(1,':k'), end
            xlim([t(1),t(end)])
            title(sprintf('%s turns (n=%d flies)',bin_labels{b},n_f))
            if b == 1, ylabel(sprintf('%s\n%s',geno_keep{gi},ylab),'Interpreter','none'), end
            if gi == numel(geno_keep), xlabel(xlab), end
        end
    end
end

function [fly_vel,bump_vel,mu_w,valid,fz_bin] = trial_gain_by_position_vectors(trial,mu_raw,rho_raw,fz,mu_smooth_s,vel_smooth_s,lag_s,vel_thresh,bump_thresh,rho_thresh,vel_max,pos_edges)
    % fictrac-timebase fly_vel / bump_vel / rho / wrapped mu, with the
    % bump-side vectors shifted lag_s later than the fly-side ones (same
    % shift-and-trim as lpsp_rnai_claude_v2.m's trial_gain_vectors_v3,
    % with the lag in seconds -> this trial's frames). Also returns, for
    % the mapping check, the per-PVA-bin sum and count of z-scored
    % fluorescence across all 2n clusters (on the imaging timebase).
    xf = trial.ft.xf(:); dt = median(diff(xf));
    [mu_f,rho_f] = bump_on_fictrac(trial,mu_raw,rho_raw,mu_smooth_s);
    bump_vel_full = gradient(mu_f)/dt;
    fly_vel_full  = smoothdata(trial.ft.r_speed(:),'gaussian',max(1,round(vel_smooth_s/dt)));
    lag = round(lag_s/dt);
    if lag > 0
        fly_vel  = fly_vel_full(1:end-lag);
        bump_vel = bump_vel_full(lag+1:end);
        rho_i    = rho_f(lag+1:end);
        mu_w     = wrap_pi(mu_f(lag+1:end));
    else
        fly_vel = fly_vel_full; bump_vel = bump_vel_full; rho_i = rho_f; mu_w = wrap_pi(mu_f);
    end
    valid = abs(fly_vel) > vel_thresh & abs(fly_vel) < vel_max & abs(bump_vel) < bump_thresh & rho_i > rho_thresh & ~isnan(bump_vel);

    % mapping check: bin the (unsmoothed-in-time) PVA angle per imaging
    % frame and accumulate the z-scored fluorescence of every cluster
    n_pos = numel(pos_edges)-1;
    mu_im = wrap_to_range(mu_raw(:),pos_edges(1));
    rho_ok = rho_raw(:) > rho_thresh;
    pb = discretize(mu_im,pos_edges);
    fz_bin.sum = zeros(size(fz,1),n_pos); fz_bin.n = zeros(1,n_pos);
    for p = 1:n_pos
        sel = pb==p & rho_ok;
        fz_bin.sum(:,p) = sum(fz(:,sel),2,'omitnan');
        fz_bin.n(p) = sum(sel);
    end
end

function y = wrap_pi(x)
    y = mod(x+pi,2*pi)-pi;
end

function y = wrap_to_range(x,lo)
    % wrap angles into [lo, lo+2*pi)
    y = mod(x-lo,2*pi)+lo;
end

function s = sem(X,dim)
    s = std(X,0,dim,'omitnan') ./ sqrt(max(1,sum(~isnan(X),dim)));
end

function plot_sem(t,X,c)
    % mean +/- SEM shading across rows of X, line in color c
    t = t(:)';
    m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m);
    if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end

function export_fig(out_dir,png_dir,name)
    % PDF into out_dir (ugly_figures/rnai_turns), plus a PNG preview into
    % png_dir if one is given. Axes toolbars are hidden first: in batch
    % mode exportgraphics sometimes rasterizes the hover toolbar into the
    % image (it warns "Exported image displays axes toolbar"), so every
    % axes' toolbar is switched off and, if the warning still fires, the
    % export is repeated once as MATLAB's own message suggests. A locked
    % destination (file open in a viewer) only warns so the rest of the
    % script keeps running.
    for a = findall(gcf,'-property','Toolbar')'
        try, a.Toolbar.Visible = 'off'; catch, end
    end
    targets = {fullfile(out_dir,[name '.pdf'])};
    if ~isempty(png_dir)
        if ~exist(png_dir,'dir'), mkdir(png_dir); end
        targets{end+1} = fullfile(png_dir,[name '.png']);
    end
    for t = 1:numel(targets)
        for attempt = 1:2
            lastwarn('');
            try
                if endsWith(targets{t},'.png')
                    exportgraphics(gcf,targets{t},'Resolution',120)
                else
                    exportgraphics(gcf,targets{t})
                end
            catch ME
                warning('could not save %s (%s)', targets{t}, ME.message);
                break
            end
            if ~contains(lastwarn,'toolbar'), break, end
        end
    end
end
