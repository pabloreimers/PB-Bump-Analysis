%% lpsp_tnt_turns_claude
% Turn-resolved bump motion in the LPsP>TNT silencing dataset (empty>tnt
% control vs. lpsp>tnt, PB imaging), closed loop and dark in separate
% figures. A copy of lpsp_kir_turns_claude.m (EB imaging, LPsP>Kir) with
% only the dataset-specific lines changed, so the two silencing datasets
% are analyzed identically: same turn detection, thresholds, turn-size
% bins, lag, smoothing and figures. That script in turn follows question 1
% of lpsp_rnai_turns_claude.m.
%
% Data handling follows lpsp_tnt_claude.m (checked directly on
% lpsp_tnt_data_20240307.mat, 89 trials):
%   - genotype comes from the trial folder name ("..._lpsp_tnt" /
%     "..._empty_tnt", case varies in the raw data so matching is
%     case-insensitive), and light condition from each trial's own
%     trialSettings.csv under Z:\pablo\lpsp_tnt\ (no ft.pattern stored).
%   - PB vs. EB makes no difference here: im.z/im.alpha are hemisphere-
%     resolved (32 rows, alpha repeated twice), but this analysis only
%     uses the stored im.mu / im.rho, which are one value per frame.
%   - HEADING is -ft.cue (unwrapped), as in lpsp_tnt_claude.m. -ft.cue
%     moves at ~0.75x the integrated r_speed in most trials (median
%     velocity slope 0.75, 0.24-0.80, median corr 0.91; Kir: ~0.67), so
%     turns are sized in cue/arena degrees and the closed-loop reference is
%     bump = 1 x fly. The turn threshold is kept at Kir's 0.20 rad/s of cue
%     velocity (~0.27 rad/s of ball rotation here) so both silencing
%     datasets use the same detection rule. Set heading_source = 'r_speed'
%     to size turns by ball rotation instead (then expect ~0.75 in closed loop).
%   - the bump is the dataset's own im.mu / im.rho, smoothed with the same
%     0.2 s Gaussian on the unwrapped mu as the other turn scripts. im.mu
%     moves in the same sign as -ft.cue.
%   - there is no ft.xb, so imaging frames are spread evenly over ft.xf
%     (as lpsp_tnt_claude.m does).
%
% Run from the repo root, one %% section at a time (or top to bottom).
% Figures are exported as PDFs to ugly_figures/tnt_turns/ (per-fly
% overlays in its per_fly/ subfolder).
addpath(fullfile(pwd,'circ_stats'));
set(groot,'defaultAxesToolbarVisible','off'); % keeps the axes toolbar out of exported images in batch mode

%% 0) parameters
data_dir    = fullfile('.data');
source_file = 'lpsp_tnt_data_20240307.mat';
base_dir    = 'Z:\pablo\lpsp_tnt\'; % raw session tree holding each trial's trialSettings.csv (see lpsp_tnt_claude.m)

geno_keep   = {'empty>tnt','lpsp>tnt'};
geno_colors = [0.00,0.45,0.74;   % empty>tnt blue
               0.85,0.20,0.10];  % lpsp>tnt  red
include_dark = true;

heading_source = 'cue'; % 'cue' = unwrapped -ft.cue (lpsp_tnt_claude.m's heading); 'r_speed' = integrated ft.r_speed (ball rotation, ~1.33x larger)

% bump position: stored im.mu, Gaussian on the unwrapped trace
bump_smooth_s = 0.20; % s, same as lpsp_rnai_turns_claude.m
bump_lag_s    = 0.13; % s, bump lags behavior (same value as the RNAi turn analysis)

% turn detection: identical to lpsp_kir_turns_claude.m (velocity threshold
% in units of the heading signal -- 0.20 rad/s of cue velocity, ~0.27 rad/s
% of ball rotation in this dataset)
turn_vel_smooth_s = 0.50;
if strcmp(heading_source,'cue'), turn_vel_thresh = 0.20; else, turn_vel_thresh = 0.30; end
turn_max_gap_s    = 0.20;
turn_min_dur_s    = 0.15;
turn_rho_thresh   = 0.20;
turn_win_s        = [-1, 3];
turn_grid_dt      = 0.02;
turn_quiet_pre_s  = 0.50; % "quiet baseline" turns: |v| below threshold this long before onset
turn_norm_tau     = (-0.5:0.02:1.5)';

turn_bin_edges_deg = [7.5:15:97.5, 135, 180, 360];
min_turns_per_bin  = 3;
peak_settle_s      = 1.0; % peak window = 0..(median turn duration of the bin + this), as in lpsp_rnai_turns_claude.m section 6

if strcmp(heading_source,'cue'), cl_ref_gain = 1; else, cl_ref_gain = 0.75; end % expected closed-loop bump/fly ratio

export_prefix = 'tntturns';
export_dir    = fullfile('ugly_figures','tnt_turns');
per_fly_dir   = fullfile(export_dir,'per_fly');
if ~exist(export_dir,'dir'), mkdir(export_dir); end
if ~exist(per_fly_dir,'dir'), mkdir(per_fly_dir); end
if ~exist('preview_png_dir','var'), preview_png_dir = fullfile(export_dir,'pngs'); end % PNG copy of every PDF ('' = none)
per_fly_png_dir = ''; if ~isempty(preview_png_dir), per_fly_png_dir = fullfile(preview_png_dir,'per_fly'); end
if ~exist('make_per_fly_figs','var'), make_per_fly_figs = true; end

%% 1) load data + labels (same parsing as lpsp_tnt_claude.m)
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)
tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:); clear tmp
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

genotype = cell(n_trials,1);
for i = 1:n_trials
    tname = trial_folder_name(all_data(i).meta);
    if contains(tname,'_lpsp_tnt','IgnoreCase',true),      genotype{i} = 'lpsp>tnt';
    elseif contains(tname,'_empty_tnt','IgnoreCase',true), genotype{i} = 'empty>tnt';
    else,                                                  genotype{i} = '';
    end
end

d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
csv_by_name = containers.Map('KeyType','char','ValueType','char');
for k = 1:numel(d)
    [~,tname] = fileparts(fileparts(d(k).folder));
    csv_by_name(lower(tname)) = fullfile(d(k).folder,d(k).name);
end
is_dark = nan(n_trials,1);
for i = 1:n_trials
    tname = lower(trial_folder_name(all_data(i).meta));
    if ~isKey(csv_by_name,tname), warning('no trialSettings.csv for "%s"', tname); continue, end
    T = readtable(csv_by_name(tname),'TextType','string');
    is_dark(i) = contains(T.patternPath(1),'background','IgnoreCase',true);
end

fly_id = cell(n_trials,1);
for i = 1:n_trials, fly_id{i} = trial_fly_id(all_data(i).meta); end
[fly_list,~,fly_num] = unique(fly_id);
n_flies = numel(fly_list);
fly_geno = cell(n_flies,1);
for f = 1:n_flies
    g = unique(genotype(fly_num==f)); fly_geno{f} = g{1};
end

if include_dark, cond_list = [false,true]; else, cond_list = false; end
cond_name = {'closed loop','dark'};
for gi = 1:numel(geno_keep)
    fprintf('  %-10s %2d flies; trials: %d closed loop, %d dark\n', geno_keep{gi}, sum(strcmp(fly_geno,geno_keep{gi})), ...
        sum(strcmp(genotype,geno_keep{gi}) & is_dark==0), sum(strcmp(genotype,geno_keep{gi}) & is_dark==1));
end

%% 2) detect turns in every trial, extract aligned traces + per-turn displacements
use_trials = find(ismember(genotype,geno_keep) & ismember(is_dark,cond_list));
turn_t = (turn_win_s(1):turn_grid_dt:turn_win_s(2))';
turns = [];
for ti = use_trials'
    tr = trial_turns(all_data(ti),heading_source,bump_smooth_s,bump_lag_s,turn_vel_smooth_s,turn_vel_thresh, ...
        turn_max_gap_s,turn_min_dur_s,turn_rho_thresh,turn_win_s,turn_t,turn_quiet_pre_s,turn_norm_tau);
    if isempty(tr), continue, end
    [tr.trial] = deal(ti); [tr.fly] = deal(fly_num(ti)); [tr.geno] = deal(genotype{ti}); [tr.dark] = deal(logical(is_dark(ti)));
    turns = [turns; tr(:)]; %#ok<AGROW>
end
turn_size  = [turns.size_deg]';
turn_bin   = discretize(turn_size,turn_bin_edges_deg);
n_bins     = numel(turn_bin_edges_deg)-1;
bin_labels = cell(1,n_bins);
for b = 1:n_bins
    if diff(turn_bin_edges_deg(b:b+1)) <= 15.01
        bin_labels{b} = sprintf('%g%s',mean(turn_bin_edges_deg(b:b+1)),char(176));
    else
        bin_labels{b} = sprintf('%g-%g%s',turn_bin_edges_deg(b),turn_bin_edges_deg(b+1),char(176));
    end
end
turn_quiet = [turns.quiet_pre]';
turn_fly   = [turns.fly]';
turn_dark  = [turns.dark]';
turn_dur   = [turns.dur_s]';
turn_geno  = {turns.geno}';
fprintf('%d turns detected (%d within the size bins), heading = %s\n', numel(turns), sum(~isnan(turn_bin)), heading_source);
for ci = 1:numel(cond_list)
    for gi = 1:numel(geno_keep)
        sel = turn_dark==cond_list(ci) & strcmp(turn_geno,geno_keep{gi});
        fprintf('  %-11s %-10s %5d turns; per bin: %s\n', cond_name{cond_list(ci)+1}, geno_keep{gi}, sum(sel), ...
            strjoin(arrayfun(@(b) sprintf('%d',sum(sel & turn_bin==b)),1:n_bins,'UniformOutput',false),' / '));
    end
end
used_flies = unique(turn_fly);

%% figure: how turns are defined -- example trial with the most turns
turn_trials = [turns.trial]';
[~,ex_ti] = max(accumarray(turn_trials,1,[n_trials,1]));
trial = all_data(ex_ti); xf = trial.ft.xf(:); dt = median(diff(xf));
heading = trial_heading(trial,heading_source);
v_sm = gradient(smoothdata(heading,'gaussian',max(1,round(turn_vel_smooth_s/dt))))/dt;
mu_f = bump_on_fictrac(trial,bump_smooth_s);
these = turns(turn_trials==ex_ti);

figure(1); clf; set(gcf,'Name','turn definition, example trial','Position',[80,80,1300,600],'Color','w')
ax1 = subplot(2,1,1); hold on
plot(xf,v_sm,'k'); yline(turn_vel_thresh,':r'); yline(-turn_vel_thresh,':r'); yline(0,'-','Color',[.8,.8,.8])
for k = 1:numel(these)
    xs = these(k).onset_s + [0,these(k).dur_s];
    patch([xs,fliplr(xs)],[-5,-5,5,5],[1,.85,.6],'EdgeColor','none','FaceAlpha',.5)
end
ylim([-1,1]*min(5,max(abs(v_sm))*1.1)); ylabel(sprintf('heading velocity, %.2fs Gaussian (rad/s)',turn_vel_smooth_s))
title(sprintf('trial %d (%s, %s): %d turns; shaded = detected turns (|v|>%.2f rad/s, same sign, >%.2fs)', ...
    ex_ti,genotype{ex_ti},cond_name{is_dark(ex_ti)+1},numel(these),turn_vel_thresh,turn_min_dur_s),'Interpreter','none')
ax2 = subplot(2,1,2); hold on
plot(xf,rad2deg(heading-heading(1)),'k','DisplayName',sprintf('fly heading (%s)',heading_source))
plot(xf,rad2deg(mu_f-mu_f(1)),'Color',geno_colors(strcmp(geno_keep,genotype{ex_ti}),:),'DisplayName','bump (unwrapped im.mu)')
xlabel('time (s)'); ylabel('cumulative rotation (deg)'); legend('Location','best')
linkaxes([ax1,ax2],'x'); xlim([xf(1),min(xf(end),xf(1)+120)])
export_fig(export_dir,preview_png_dir,sprintf('%s_fig1_turn_definition_example',export_prefix))

%% 3) per fly x condition x turn-size bin: mean displacements (+ optional per-fly overlay figures)
n_turns_fb = zeros(n_flies,numel(cond_list),n_bins);
dfly_fb    = nan(n_flies,numel(cond_list),n_bins); % mean fly rotation over the turn (deg)
dbump_fb   = nan(n_flies,numel(cond_list),n_bins); % mean bump rotation over the turn (deg, lagged)
ratio_fb   = nan(n_flies,numel(cond_list),n_bins); % through-origin slope dbump ~ dfly across this fly's turns in the bin
for f = used_flies'
    gi = find(strcmp(geno_keep,fly_geno{f}));
    for ci = 1:numel(cond_list)
        sel_fc = turn_fly==f & turn_dark==cond_list(ci);
        if ~any(sel_fc), continue, end
        if make_per_fly_figs
            fig = figure(100+f); clf
            set(fig,'Name',sprintf('turn overlays, fly %d',f),'Position',[60,60,220*n_bins,330],'Color','w')
        end
        for b = 1:n_bins
            sel = find(sel_fc & turn_bin==b);
            n_turns_fb(f,ci,b) = numel(sel);
            if numel(sel) < min_turns_per_bin, continue, end
            dfly = [turns(sel).dfly_deg]'; dbump = [turns(sel).dbump_deg]';
            dfly_fb(f,ci,b)  = mean(dfly);
            dbump_fb(f,ci,b) = mean(dbump);
            ratio_fb(f,ci,b) = dfly \ dbump;
            if make_per_fly_figs
                subplot(1,n_bins,b); hold on
                plot_sem(turn_t,cat(2,turns(sel).fly_trace)',[0,0,0])
                plot_sem(turn_t,cat(2,turns(sel).bump_trace)',geno_colors(gi,:))
                xline(0,':k'); xlim(turn_win_s)
                title(sprintf('%s turns (n=%d)',bin_labels{b},numel(sel)))
                xlabel('time from turn onset (s)')
                if b == 1, ylabel('rotation (deg)'), end
            end
        end
        if make_per_fly_figs
            sgtitle(sprintf('%s, %s -- fly %d (%s): fly heading (black) vs. bump (color), mean +/- SEM', ...
                fly_geno{f},cond_name{cond_list(ci)+1},f,fly_short_label(fly_list{f})),'Interpreter','none')
            export_fig(per_fly_dir,per_fly_png_dir,sprintf('%s_fig2_overlay_%s_%s_fly%02d',export_prefix,strrep(fly_geno{f},'>','-'),strrep(cond_name{cond_list(ci)+1},' ',''),f))
            close(fig)
        end
    end
end

%% figures: group-level aligned overlays (mean of per-fly means), one row per genotype
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
        set(gcf,'Name',sprintf('group turn overlays (%s), %s',variants(v).label,cond_name{cond_list(ci)+1}),'Position',[60,60,220*n_bins,260*numel(geno_keep)],'Color','w')
        plot_group_overlay(FM,BM,variants(v).t,geno_keep,fly_geno,geno_colors,bin_labels,variants(v).xlab,variants(v).ylab,strcmp(variants(v).ffield,'fly_norm'))
        sgtitle(sprintf('%s, %s: fly heading (black) vs. bump (color), mean +/- SEM across flies (each fly = mean of its >= %d turns per bin)', ...
            cond_name{cond_list(ci)+1},variants(v).label,min_turns_per_bin),'Interpreter','none')
        export_fig(export_dir,preview_png_dir,sprintf('%s_%s_%s',export_prefix,variants(v).suffix,strrep(cond_name{cond_list(ci)+1},' ','')))
    end
end

%% figure: summary -- bump displacement vs. fly displacement per turn-size bin, per fly, by genotype
% per-bin genotype comparison: ranksum on the per-fly through-origin slopes
% (empty>tnt vs. lpsp>tnt), printed and marked above each bin
for ci = 1:numel(cond_list)
    cname = cond_name{cond_list(ci)+1};
    figure(20+ci); clf
    set(gcf,'Name',sprintf('turn-size summary, %s',cname),'Position',[60,60,1500,500],'Color','w')
    subplot(1,3,1); hold on
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        X = reshape(dfly_fb(flies_g,ci,:),numel(flies_g),[]); Y = reshape(dbump_fb(flies_g,ci,:),numel(flies_g),[]);
        plot(X',Y','-','Color',[geno_colors(gi,:),.25],'Marker','.','MarkerSize',8,'HandleVisibility','off')
        errorbar(mean(X,1,'omitnan'),mean(Y,1,'omitnan'),sem(Y,1),sem(Y,1),sem(X,1),sem(X,1),'o-','Color',geno_colors(gi,:), ...
            'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2,'DisplayName',sprintf('%s (n=%d flies)',geno_keep{gi},sum(any(~isnan(X),2))))
    end
    xl = [0,max([turn_bin_edges_deg(end-1); dfly_fb(:)],[],'omitnan')*1.05];
    plot(xl,xl,':k','DisplayName','1:1')
    if cl_ref_gain ~= 1, plot(xl,cl_ref_gain*xl,'--','Color',[.5,.5,.5],'DisplayName',sprintf('%.2f:1',cl_ref_gain)), end
    xlabel(sprintf('fly rotation over turn (deg, %s)',heading_source)); ylabel(sprintf('bump rotation over turn (deg, lag %.2fs)',bump_lag_s))
    legend('Location','northwest','Interpreter','none'); title('mean displacement per turn-size bin')
    axis square; xlim(xl); ylim([max(-50,min([-10; dbump_fb(:)],[],'omitnan')),xl(2)])

    subplot(1,3,2); hold on
    R = cell(1,numel(geno_keep));
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        R{gi} = reshape(ratio_fb(flies_g,ci,:),numel(flies_g),[]);
        xj = (1:n_bins) + (gi-(numel(geno_keep)+1)/2)*0.25;
        plot(repmat(xj,size(R{gi},1),1)',R{gi}','.','Color',[geno_colors(gi,:),.35],'MarkerSize',10,'HandleVisibility','off')
        errorbar(xj,mean(R{gi},1,'omitnan'),sem(R{gi},1),'o-','Color',geno_colors(gi,:),'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2)
    end
    yline(cl_ref_gain,'--','Color',[.5,.5,.5]); yline(0,'-','Color',[.85,.85,.85])
    ylim([-0.5,2]); yl = ylim;
    fprintf('\n=== %s: per-fly bump/fly slope per turn-size bin, %s vs. %s (ranksum) ===\n',cname,geno_keep{1},geno_keep{2});
    for b = 1:n_bins
        a = R{1}(:,b); a = a(~isnan(a)); c = R{2}(:,b); c = c(~isnan(c));
        if numel(a) < 3 || numel(c) < 3, continue, end
        p = ranksum(a,c);
        fprintf('  %-9s %s %.2f +/- %.2f (n=%d)   %s %.2f +/- %.2f (n=%d)   p=%.3g\n',bin_labels{b}, ...
            geno_keep{1},mean(a),std(a)/sqrt(numel(a)),numel(a),geno_keep{2},mean(c),std(c)/sqrt(numel(c)),numel(c),p);
        text(b,yl(2)-0.05*diff(yl),p_stars(p),'HorizontalAlignment','center','FontSize',9)
    end
    xticks(1:n_bins); xticklabels(bin_labels); xlim([0.5,n_bins+0.5])
    xlabel('turn size (fly rotation)'); ylabel('bump / fly displacement (per-fly through-origin slope)')
    title('gain per turn-size bin (ranksum between genotypes; y clipped)')

    subplot(1,3,3); hold on
    for gi = 1:numel(geno_keep)
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        N = reshape(n_turns_fb(flies_g,ci,:),numel(flies_g),[]);
        xj = (1:n_bins) + (gi-(numel(geno_keep)+1)/2)*0.35;
        bar(xj,mean(N,1),0.35,'FaceColor',geno_colors(gi,:),'EdgeColor','none','FaceAlpha',.7)
    end
    xticks(1:n_bins); xticklabels(bin_labels); xlim([0.5,n_bins+0.5])
    xlabel('turn size'); ylabel('mean turns per fly'); title('turn counts'); legend(geno_keep,'Interpreter','none')
    sgtitle(sprintf('%s: turn-size-resolved bump motion, %s',cname,strjoin(geno_keep,' vs. ')),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig4_turn_size_summary_%s',export_prefix,strrep(cname,' ','')))
end

%% figure: peak of each fly's mean trace vs. turn size, by genotype (all turns and quiet-baseline turns)
% same measure as lpsp_rnai_turns_claude.m section 6 (real-time version):
% per fly x bin, the max of its mean aligned trace over 0..(that bin's
% median turn duration + peak_settle_s), for the fly's heading and the bump.
% Unlike the displacement-at-offset slope above, this captures the bump's
% eventual excursion even though it lags the fly.
cat_x = nan(n_flies,1);
for gi = 1:numel(geno_keep), cat_x(strcmp(fly_geno,geno_keep{gi})) = gi; end
turn_sets = struct('tag',{'allturns','quiet'},'label',{'all turns','quiet-baseline turns'},'mask',{true(size(turn_quiet)),turn_quiet});
for ts = 1:numel(turn_sets)
for ci = 1:numel(cond_list)
    cname = cond_name{cond_list(ci)+1};
    sel = turn_dark==cond_list(ci) & turn_sets(ts).mask;
    FM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,'fly_trace',min_turns_per_bin);
    BM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,'bump_trace',min_turns_per_bin);
    fly_pk = nan(n_flies,n_bins); bump_pk = nan(n_flies,n_bins); win_end = nan(1,n_bins);
    for b = 1:n_bins
        win_end(b) = min(turn_win_s(2), median(turn_dur(sel & turn_bin==b)) + peak_settle_s);
        t_pk = turn_t >= 0 & turn_t <= win_end(b);
        fly_pk(:,b)  = max(reshape(FM(:,b,t_pk),n_flies,[]),[],2);
        bump_pk(:,b) = max(reshape(BM(:,b,t_pk),n_flies,[]),[],2);
    end
    ratio_pk = bump_pk ./ fly_pk;

    figure(40+10*ts+ci); clf
    set(gcf,'Name',sprintf('peak summary vs. turn size, %s, %s',turn_sets(ts).label,cname),'Position',[60,60,1500,480],'Color','w')
    panels = {fly_pk,'peak fly rotation (deg)','fly heading peak'; bump_pk,'peak bump rotation (deg)','bump peak'; ratio_pk,'bump peak / fly peak','bump/fly peak ratio'};
    fprintf('\n=== %s, %s: bump/fly PEAK ratio per turn-size bin (ranksum) ===\n',cname,turn_sets(ts).label);
    for k = 1:3
        subplot(1,3,k); hold on
        V = panels{k,1};
        for gi = 1:numel(geno_keep)
            Y = V(cat_x==gi,:);
            xj = (1:n_bins) + (gi-(numel(geno_keep)+1)/2)*0.25;
            plot(repmat(xj,size(Y,1),1)',Y','.','Color',[geno_colors(gi,:),.3],'MarkerSize',8,'HandleVisibility','off')
            errorbar(xj,mean(Y,1,'omitnan'),sem(Y,1),'o-','Color',geno_colors(gi,:),'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2, ...
                'DisplayName',sprintf('%s (n=%d flies)',geno_keep{gi},sum(any(~isnan(Y),2))))
        end
        if k == 3
            yline(cl_ref_gain,'--','Color',[.5,.5,.5],'HandleVisibility','off'); yline(0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
            ylim([-0.5,2]); yl = ylim;
            for b = 1:n_bins
                a = V(cat_x==1,b); a = a(~isnan(a)); c = V(cat_x==2,b); c = c(~isnan(c));
                if numel(a) < 3 || numel(c) < 3, continue, end
                p = ranksum(a,c);
                fprintf('  %-9s %s %.2f +/- %.2f (n=%d)   %s %.2f +/- %.2f (n=%d)   p=%.3g\n',bin_labels{b}, ...
                    geno_keep{1},mean(a),std(a)/sqrt(numel(a)),numel(a),geno_keep{2},mean(c),std(c)/sqrt(numel(c)),numel(c),p);
                text(b,yl(2)-0.05*diff(yl),p_stars(p),'HorizontalAlignment','center','FontSize',9)
            end
        end
        if k == 1, legend('Location','northwest','Interpreter','none'), end
        xticks(1:n_bins); xticklabels(bin_labels); xlim([0.5,n_bins+0.5])
        xlabel('turn size (fly rotation)'); ylabel(panels{k,2}); title(panels{k,3})
    end
    sgtitle(sprintf('%s, %s: peak of per-fly mean trace (max over 0..median turn duration + %.1f s) vs. turn size; faint = flies, thick = mean +/- SEM', ...
        cname,turn_sets(ts).label,peak_settle_s),'Interpreter','none')
    export_fig(export_dir,preview_png_dir,sprintf('%s_fig5_peak_%s_summary_%s',export_prefix,turn_sets(ts).tag,strrep(cname,' ','')))
end
end

%% functions
function name = trial_folder_name(meta_path)
    % trial folder = last path component of meta, or the one above a registration_NNN subfolder (as lpsp_tnt_claude.m)
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        name = parts{end-1};
    else
        name = parts{end};
    end
end

function fid = trial_fly_id(meta_path)
    % meta = ...\<date folder>\fly N\<trial folder>\registration_NNN -> "...\<date folder>\fly N"
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
    if isempty(fly_part)
        fid = meta_path;
    else
        fid = strjoin(parts(1:fly_part(1)),filesep);
    end
end

function lbl = fly_short_label(fly_path)
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function h = trial_heading(trial,source)
    % the fly's cumulative rotation (rad, unwrapped): -ft.cue (the arena
    % position, lpsp_tnt_claude.m's heading) or the integral of ft.r_speed
    if strcmp(source,'cue')
        h = unwrap(-trial.ft.cue(:));
    else
        r = trial.ft.r_speed(:); r(isnan(r)) = 0;
        h = cumtrapz(trial.ft.xf(:),r);
    end
end

function [mu_f,rho_f] = bump_on_fictrac(trial,mu_smooth_s)
    % unwrapped, Gaussian-smoothed im.mu and im.rho on the fictrac
    % timebase; no ft.xb in this dataset, so frames are spread evenly over
    % ft.xf (as lpsp_tnt_claude.m does)
    xf = trial.ft.xf(:);
    n_im = numel(trial.im.mu);
    xb = linspace(xf(1),xf(end),n_im)';
    win_im = max(1,round(mu_smooth_s/median(diff(xb))));
    mu_s = smoothdata(unwrap(trial.im.mu(:)),'gaussian',win_im);
    mu_f  = interp1(xb,mu_s,xf,'linear');
    rho_f = interp1(xb,trial.im.rho(:),xf,'linear');
end

function turns = trial_turns(trial,heading_source,mu_smooth_s,lag_s,vel_smooth_s,vel_thresh,max_gap_s,min_dur_s,rho_thresh,win_s,t_grid,quiet_pre_s,tau_grid)
    % same detection and per-turn extraction as lpsp_rnai_turns_claude.m's
    % trial_turns, except that the velocity used for detection is the
    % derivative of the Gaussian-smoothed heading (heading_source), so turn
    % size and turn detection come from the same signal
    xf = trial.ft.xf(:); dt = median(diff(xf));
    heading = trial_heading(trial,heading_source);
    [mu_f,rho_f] = bump_on_fictrac(trial,mu_smooth_s);

    win = max(1,round(vel_smooth_s/dt));
    v_sm = gradient(smoothdata(heading,'gaussian',win))/dt;
    mu_vel = smoothdata(gradient(mu_f)/dt,'gaussian',win);
    max_gap = max(1,round(max_gap_s/dt)); min_dur = max(1,round(min_dur_s/dt)); n_quiet = max(1,round(quiet_pre_s/dt));

    turns = struct('onset_s',{},'dur_s',{},'size_deg',{},'dir',{},'quiet_pre',{},'fly_trace',{},'bump_trace',{}, ...
                   'fly_vel_trace',{},'bump_vel_trace',{},'fly_norm',{},'bump_norm',{},'dfly_deg',{},'dbump_deg',{},'rho_mean',{});
    for s = [1,-1]
        active = s*v_sm > vel_thresh;
        active = imclose(active,ones(max_gap+1,1));
        active = bwareaopen(active,min_dur);
        d = diff([false;active;false]);
        on = find(d==1); off = find(d==-1)-1;
        for k = 1:numel(on)
            t_on = xf(on(k)); t_off = xf(off(k));
            if t_on + win_s(1) < xf(1) || max(t_on + win_s(2), t_off + lag_s) > xf(end), continue, end
            rng = on(k):off(k);
            if mean(rho_f(rng),'omitnan') < rho_thresh, continue, end
            dfly = heading(off(k)) - heading(on(k));
            if abs(dfly) < 1e-6, continue, end
            mu_on  = interp1(xf,mu_f,t_on + lag_s);
            mu_off = interp1(xf,mu_f,t_off + lag_s);
            tq = t_on + t_grid;
            bump_tr = interp1(xf,mu_f,tq);
            if isnan(mu_on) || isnan(mu_off) || any(isnan(bump_tr)), continue, end
            tau_q = t_on + tau_grid*(t_off - t_on);
            mu_at_on = interp1(xf,mu_f,t_on);
            T.onset_s   = t_on;
            T.dur_s     = t_off - t_on;
            T.size_deg  = rad2deg(abs(dfly));
            T.dir       = s;
            T.quiet_pre = on(k) > n_quiet && all(abs(v_sm(on(k)-n_quiet:on(k)-1)) < vel_thresh);
            T.fly_trace  = single(s*rad2deg(interp1(xf,heading,tq) - heading(on(k))));
            T.bump_trace = single(s*rad2deg(bump_tr - mu_at_on));
            T.fly_vel_trace  = single(s*rad2deg(interp1(xf,v_sm,tq)));
            T.bump_vel_trace = single(s*rad2deg(interp1(xf,mu_vel,tq)));
            T.fly_norm  = single(s*rad2deg(interp1(xf,heading,tau_q) - heading(on(k))));
            T.bump_norm = single(s*rad2deg(interp1(xf,mu_f,tau_q) - mu_at_on));
            T.dfly_deg  = s*rad2deg(dfly);
            T.dbump_deg = s*rad2deg(mu_off - mu_on);
            T.rho_mean  = mean(rho_f(rng),'omitnan');
            turns(end+1) = T; %#ok<AGROW>
        end
    end
    [~,order] = sort([turns.onset_s]); turns = turns(order);
end

function M = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,field,min_turns)
    % per fly x turn-size bin: mean of one aligned-trace field (NaN if < min_turns); n_flies x n_bins x n_t
    n_t = numel(turns(1).(field));
    M = nan(n_flies,n_bins,n_t);
    for f = unique(turn_fly(sel))'
        for b = 1:n_bins
            idx = find(sel & turn_fly==f & turn_bin==b);
            if numel(idx) < min_turns, continue, end
            M(f,b,:) = mean(cat(2,turns(idx).(field))',1,'omitnan');
        end
    end
end

function plot_group_overlay(FM,BM,t,geno_keep,fly_geno,geno_colors,bin_labels,xlab,ylab,mark_offset)
    % one row per genotype, one column per turn-size bin: mean +/- SEM
    % across flies of the per-fly mean traces (fly black, bump genotype color)
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
            end
            xline(0,':k'); if mark_offset, xline(1,':k'), end
            xlim([t(1),t(end)])
            title(sprintf('%s turns (n=%d flies)',bin_labels{b},n_f))
            if b == 1, ylabel(sprintf('%s\n%s',geno_keep{gi},ylab),'Interpreter','none'), end
            if gi == numel(geno_keep), xlabel(xlab), end
        end
    end
end

function s = p_stars(p)
    if p < 0.001, s = '***'; elseif p < 0.01, s = '**'; elseif p < 0.05, s = '*'; else, s = 'n.s.'; end
end

function s = sem(X,dim)
    s = std(X,0,dim,'omitnan') ./ sqrt(max(1,sum(~isnan(X),dim)));
end

function plot_sem(t,X,c)
    t = t(:)';
    m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m);
    if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end

function export_fig(out_dir,png_dir,name)
    % PDF into out_dir, plus a PNG copy into png_dir if given (same helper as lpsp_rnai_turns_claude.m)
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
                    exportgraphics(gcf,targets{t},'Resolution',200)
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
