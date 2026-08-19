%% lpsp_tnt_claude
% Loads the LPsP>TNT silencing dataset (empty>tnt control vs lpsp>tnt
% silencing) and repeats the same analysis pipeline built for the
% lpsp_kir_redo dataset in lpsp_kir_claude.m: genotype/light-condition
% labeling, PB activity + bump position visualization, velocity-gain
% scatter with a per-group optimal lag, sliding-window position-tracking
% accuracy, circular occupancy/"stickiness", and bump-mobility
% (path-length ratio) during walking bouts.
%
% This dataset was investigated directly (not assumed to match kir) before
% writing this script. Three real differences from lpsp_kir_claude.m:
%
% 1) GENOTYPE TAG: folder names end "..._empty_tnt" or "..._lpsp_tnt" (not
%    "_kir"), and case is inconsistent in the raw data itself -- both
%    "lpsp_tnt" and "LPsP_tnt" appear (e.g. trial 19-22 vs 27-28 in the
%    89-trial dataset) -- confirmed case-insensitive matching leaves 0/89
%    trials unmatched.
%
% 2) FILE LAYOUT IS SIMPLER: all_data.meta points directly into
%    Z:\pablo\lpsp_tnt\... and trialSettings.csv lives right under that
%    same tree (confirmed: 92 csv files on disk under this base_dir, all
%    89 all_data trials matched by name with zero ambiguity -- the 3 extra
%    csv names belong to raw trials that were never processed into
%    all_data). Unlike lpsp_kir_redo, there's no separate "\stacks\" tree
%    to worry about -- but the exact same name-indexed csv-matching code
%    from lpsp_kir_claude.m is reused as-is below, since it doesn't assume
%    a stacks split either way and this dataset's date folders aren't all
%    clean "YYYYMMDD" (one is "20231114_2"), which the fly-ID regex (keyed
%    off the "fly N" folder, not the date format) handles fine regardless.
%
% 3) im.alpha/im.z ARE HEMISPHERE-RESOLVED (32 rows), NOT PRE-COLLAPSED TO
%    16 LIKE lpsp_kir_redo: confirmed directly across multiple trials that
%    im.alpha repeats the exact same 16-value angular sequence twice
%    (wedges 1:16 and 17:32), the same structure
%    lpsp_compartments_claude_script.m documented for its own datasets.
%    Plotting all 32 rows directly against that repeating (non-monotonic)
%    angle vector would break imagesc's Y-axis (it assumes monotonic X/Y).
%    Per explicit user direction, the PB-activity heatmap figures below
%    show BOTH hemispheres stacked on a plain wedge-index y-axis (1-32,
%    with a divider line at the hemisphere boundary) rather than
%    collapsing to one hemisphere-averaged 16-row panel; bump position and
%    heading (both single circular values, unaffected by the 16-vs-32
%    wedge-count difference) are converted to a continuous wedge-index
%    coordinate and overlaid on BOTH hemisphere bands, since the same
%    physical bump position is visible in both. All other analyses here
%    (velocity gain, position tracking, occupancy, mobility) operate on
%    im.mu/im.rho, which are already single-valued per frame regardless of
%    how many raw wedges went into computing them -- so those sections are
%    unaffected by the hemisphere difference and are unchanged from
%    lpsp_kir_claude.m.
%
% Everything else -- thresholds (vel_thresh, rho_thresh, rot_thresh_track,
% etc.), the lag-sweep search range, and the walking-bout mobility
% parameters -- is reused unchanged from lpsp_kir_claude.m: fictrac rate
% (60Hz), trial duration (600s), and im.rho scale (median 0.37, vs. kir's
% ~0.4) are all close enough that re-tuning wasn't warranted.

%% load data
data_dir    = fullfile('.data');
source_file = 'lpsp_tnt_data_20240307.mat';
base_dir    = 'Z:\pablo\lpsp_tnt\';
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype per trial, from the raw folder name in all_data.meta
% matched against the trial folder name only (e.g. "20240109-1_epg_7f_lpsp_tnt"),
% not the full meta path -- the base data folder itself is named
% "lpsp_tnt", which would otherwise falsely match a plain
% contains(meta,'lpsp_tnt') check on every trial, including empty>tnt ones
% (same false-positive risk already found and fixed for lpsp_kir_claude.m).
genotype = cell(1,n_trials);
for i = 1:n_trials
    tname = trial_folder_name(all_data(i).meta);
    if contains(tname,'_lpsp_tnt','IgnoreCase',true)
        genotype{i} = 'lpsp>tnt';
    elseif contains(tname,'_empty_tnt','IgnoreCase',true)
        genotype{i} = 'empty>tnt';
    else
        genotype{i} = '';
    end
end
n_unlabeled = sum(cellfun(@isempty,genotype));
if n_unlabeled > 0
    warning('%d/%d trials did not match either genotype pattern in their folder name', n_unlabeled, n_trials)
end

%% light condition per trial: ft.pattern isn't stored here, read it from each trial's own trialSettings.csv
has_pattern = arrayfun(@(s) isfield(s.ft,'pattern') && ~isempty(s.ft.pattern), all_data);
fprintf('\ntrials with ft.pattern already stored: %d/%d\n', sum(has_pattern), n_trials);

is_dark     = nan(1,n_trials);
pattern_str = cell(1,n_trials);

if ~all(has_pattern)
    % one recursive dir() over the whole raw tree, keyed by trial folder
    % name, rather than one dir() call per trial
    d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
    csv_by_name = containers.Map('KeyType','char','ValueType','char');
    for k = 1:numel(d)
        [~,tname] = fileparts(fileparts(d(k).folder)); % strip trailing \csv
        if isKey(csv_by_name,lower(tname))
            warning('multiple trialSettings.csv share the trial folder name "%s"', tname);
        end
        csv_by_name(lower(tname)) = fullfile(d(k).folder,d(k).name);
    end
    fprintf('indexed %d trialSettings.csv under %s\n', numel(d), base_dir);
end

for i = 1:n_trials
    if has_pattern(i)
        pattern_str{i} = char(all_data(i).ft.pattern);
    else
        tname = lower(trial_folder_name(all_data(i).meta));
        if ~isKey(csv_by_name,tname)
            warning('no trialSettings.csv found for trial folder "%s"', tname);
            continue
        end
        T = readtable(csv_by_name(tname),'TextType','string');
        pattern_str{i} = char(T.patternPath(1));
    end
    is_dark(i) = contains(pattern_str{i},'background','IgnoreCase',true);
end

fprintf('\n=== light condition summary ===\n');
fprintf('closed loop: %d\n', sum(is_dark==0));
fprintf('dark:        %d\n', sum(is_dark==1));
fprintf('unresolved:  %d\n', sum(isnan(is_dark)));
for i = find(isnan(is_dark))
    fprintf('  %s\n', all_data(i).meta);
end

%% fly ID per trial (date folder + "fly N"), and a genotype-consistency check
fly_id = cell(1,n_trials);
for i = 1:n_trials
    fly_id{i} = trial_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
fly_num = fly_num(:)'; % keep row-vector, matching is_dark/genotype -- unique() returns fly_num as a column, and fly_num==f & is_dark==c-1 would otherwise implicitly broadcast to an n_trials x n_trials matrix
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

%% plot: number of flies per genotype
geno_order = {'empty>tnt','lpsp>tnt'};
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

figure(1); clf
set(gcf,'Name','flies per genotype')
bar(fly_counts)
set(gca,'XTickLabel',geno_order)
ylabel('number of flies')
title(sprintf('lpsp\\_tnt: flies per genotype (n=%d flies, %d trials)', n_flies, n_trials))
for k = 1:numel(geno_order)
    text(k, fly_counts(k)+0.1, num2str(fly_counts(k)), 'HorizontalAlignment','center')
end

%% figures: PB activity (im.z) with bump position (im.mu) and fly heading overlay, one row per fly, two columns (closed loop | dark)
% two figures, split by genotype (empty>tnt / lpsp>tnt) since a single
% ~40-row figure would be unreadable. within a fly's column, if more than
% one trial shares that light condition, those trials are concatenated
% frame-by-frame (not by real elapsed time -- trials aren't necessarily
% the same duration) with a dotted vertical line at each trial boundary,
% so every fly still gets exactly one panel per condition.
%
% UNLIKE lpsp_kir_claude.m: this dataset's im.z/im.alpha are
% hemisphere-resolved (32 rows; wedges 1:16 and 17:32 repeat the same
% angular sequence -- confirmed directly), not pre-collapsed to one
% 16-row hemisphere-equivalent. Plotting all 32 rows against that
% repeating (non-monotonic) angle vector would break imagesc's Y-axis. Per
% explicit user direction, both hemispheres are shown stacked on a plain
% wedge-index y-axis (1-32, divider line at the hemisphere boundary)
% instead of collapsing to one hemisphere-averaged panel. Bump position
% (mu) and heading are single circular values (unaffected by wedge count),
% so each is converted to a continuous wedge-index coordinate
% (theta_to_wedge_row, below) and drawn on BOTH hemisphere bands -- the
% same physical bump position is visible in both hemispheres at once.
%
% heading is -ft.cue (ft.cue is the closed-loop position signal driving the
% panels, which continues to track the fly's own rotation in dark trials
% too -- there's just no visible pattern -- so it's a valid heading proxy
% in both conditions), same convention as lpsp_kir_claude.m. heading is
% colored per that same convention: cyan for closed loop, magenta for dark.
cond_label = {'closed loop','dark'};
cond_color = {[0,1,1],[1,0,1]}; % cyan (closed loop), magenta (dark)
z_clim     = [-2,4]; % shared color scale across every panel (im.z 1st-99th pctile across this dataset is about [-3.4,3.2], similar order to kir's)
n_hemi     = numel(all_data(1).im.alpha)/2; % wedges per hemisphere (16) -- alpha repeats this many angular positions twice

for gi = 1:numel(geno_order)
    flies_g = find(strcmp(fly_genotype,geno_order{gi}));
    n_g = numel(flies_g);

    figure(1+gi); clf
    set(gcf,'Name',sprintf('%s: PB activity + bump position',geno_order{gi}),'Position',[50,50,900,max(170*n_g,300)])
    t = tiledlayout(n_g,2,'TileSpacing','tight','Padding','compact');

    for r = 1:n_g
        f = flies_g(r);
        for c = 1:2 % 1 = closed loop, 2 = dark
            trial_idx = find(fly_num==f & is_dark==(c-1))';
            ax = nexttile(t); hold(ax,'on')
            if isempty(trial_idx)
                axis(ax,'off')
                if r == 1
                    title(ax,cond_label{c})
                end
                continue
            end

            [x_cat,~,z_cat,mu_cat,heading_cat,bounds] = concat_trials_im(all_data,trial_idx);
            wedge_idx = (1:size(z_cat,1))'; % plain index, not real angle -- see header comment
            imagesc(ax,x_cat,wedge_idx,z_cat)
            set(ax,'YDir','normal','CLim',z_clim,'XLim',[x_cat(1),x_cat(end)],'YLim',[wedge_idx(1),wedge_idx(end)])
            yline(ax,n_hemi+0.5,':','Color',[.6,.6,.6]) % boundary between the two hemisphere bands

            mu_row = theta_to_wedge_row(mu_cat,n_hemi);
            plot(ax,x_cat,mu_row,'w','LineWidth',0.5)
            plot(ax,x_cat,mu_row+n_hemi,'w','LineWidth',0.5) % same bump position, shown in both hemisphere bands

            heading_row = theta_to_wedge_row(heading_cat,n_hemi);
            plot(ax,x_cat,heading_row,'Color',cond_color{c},'LineWidth',0.5)
            plot(ax,x_cat,heading_row+n_hemi,'Color',cond_color{c},'LineWidth',0.5)

            for b = bounds(1:end-1)
                xline(ax,b,':','Color',[.6,.6,.6]);
            end

            yticks(ax,[])
            if c == 1
                ylabel(ax,fly_short_label(fly_list{f}),'Rotation',0,'HorizontalAlignment','right', ...
                    'VerticalAlignment','middle','FontSize',7)
            end
            if r == 1
                title(ax,cond_label{c})
            end
            if r == n_g
                xlabel(ax,'frame')
            else
                xticks(ax,[])
            end
        end
    end

    cb = colorbar(ax);
    cb.Layout.Tile = 'east';
    ylabel(cb,'z-score')
    title(t,sprintf('%s -- PB activity (im.z, both hemispheres stacked), bump position (im.mu, white), heading (-ft.cue: cyan=CL, magenta=dark)',geno_order{gi}),'Interpreter','none')
end

%% figures: bump velocity vs. fly heading velocity (gain), one figure per genotype x light-condition group, one subplot per FLY
% same gain/regression analysis as lpsp_kir_claude.m: fly_vel = angular
% velocity of heading (gradient of -ft.cue), bump_vel = angular velocity
% of the fitted bump position (gradient of unwrapped im.mu, interpolated
% onto the fictrac timebase). bump_vel is offset by a lag (in frames)
% relative to fly_vel (the bump lags the fly's own turning), and a point
% is kept only where the fly is actually turning
% (vel_thresh<|fly_vel|<vel_max), the bump velocity estimate isn't
% degenerate (|bump_vel|<bump_thresh), and bump concentration is high
% enough to trust (rho>rho_thresh).
%
% a fly's own trials within one light condition are concatenated
% (fly_gain_vectors, below) into one combined set of points BEFORE
% fitting/correlating, so every fly gets exactly one panel and contributes
% exactly one point to every downstream comparison. A type-2 (x-on-y...
% here y-on-x) linear fit (bump_vel ~ 1 + fly_vel) is overlaid, with its
% slope ("gain") and the raw correlation printed per fly.
vel_thresh  = 0.2;
bump_thresh = 10;
rho_thresh  = 0.2;
vel_max     = 5;
scatter_xlim = [-10,10];
scatter_ylim = [-8,8];

group_defs = struct('geno',{},'dark',{},'label',{});
for gi = 1:numel(geno_order)
    for ci = 1:numel(cond_label)
        group_defs(end+1) = struct('geno',geno_order{gi},'dark',ci-1,'label',sprintf('%s (%s)',geno_order{gi},cond_label{ci})); %#ok<SAGROW>
    end
end

% which group (1..4) each TRIAL belongs to -- just a lookup used later to
% find a trial's own group's optimal lag (group_lag_frames(cat_x(i))) for
% the position-tracking section; not itself used for any group comparison,
% so it doesn't need fly-level de-duplication.
cat_x = nan(1,n_trials);
for gd = 1:numel(group_defs)
    cat_x(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark) = gd;
end

trial_dt = nan(1,n_trials);
for i = 1:n_trials
    trial_dt(i) = mean(diff(all_data(i).ft.xf));
end

% for each group, the list of its own flies and, per fly, that fly's own
% trial indices within this group -- built once here and reused by both
% the lag sweep and the per-fly gain panels below, so a fly's trials are
% always concatenated the same way.
group_fly_list   = cell(1,numel(group_defs));
group_fly_trials = cell(1,numel(group_defs));
for gd = 1:numel(group_defs)
    rows = find(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark);
    these_flies = unique(fly_num(rows));
    group_fly_list{gd}   = these_flies;
    group_fly_trials{gd} = arrayfun(@(f) rows(fly_num(rows)==f), these_flies, 'UniformOutput',false);
end

%% optimal lag, per genotype x light-condition group: sweep a lag grid and take whichever frame lag maximizes that group's own mean corr(fly heading vel, bump vel)
% same "sweep a lag grid, average correlation across a group's own flies,
% take the argmax" approach as lpsp_kir_claude.m, one value per fly (not
% per trial) at each candidate lag.
lag_frames_grid  = -10:1:40; % ~ -0.17s to +0.67s in ~1/60s steps at this rig's fictrac rate
n_lag_grid       = numel(lag_frames_grid);
lag_seconds_grid = lag_frames_grid * mean(trial_dt,'omitnan');

group_lag_frames = nan(1,numel(group_defs));
group_lag_corr   = cell(1,numel(group_defs)); % [n_flies_in_group x n_lag_grid], one row per fly
fprintf('\n=== optimal lag per genotype x light condition (maximizes mean corr(fly heading vel, bump vel), one value per fly) ===\n');
for gd = 1:numel(group_defs)
    these_flies = group_fly_list{gd};
    trial_lists = group_fly_trials{gd};
    n_f = numel(these_flies);

    corr_mat = nan(n_f,n_lag_grid);
    for ff = 1:n_f
        for L = 1:n_lag_grid
            [fv,bv,valid] = fly_gain_vectors(all_data,trial_lists{ff},lag_frames_grid(L),vel_thresh,bump_thresh,rho_thresh,vel_max);
            if sum(valid) > 100
                corr_mat(ff,L) = corr(fv(valid),bv(valid));
            end
        end
    end
    group_lag_corr{gd} = corr_mat;

    mean_corr_grp = mean(corr_mat,1,'omitnan');
    [best_corr,best_idx] = max(mean_corr_grp);
    group_lag_frames(gd) = lag_frames_grid(best_idx);
    fprintf('  %-24s %3d frames (%.3fs), mean r = %.3f (n=%d flies)\n', group_defs(gd).label, lag_frames_grid(best_idx), lag_seconds_grid(best_idx), best_corr, n_f);
end

figure(4); clf
set(gcf,'Name','optimal lag: fly heading velocity vs. bump velocity correlation','Position',[100,100,900,700])
axs_lag = gobjects(1,numel(group_defs));
for gd = 1:numel(group_defs)
    mean_corr_grp = mean(group_lag_corr{gd},1,'omitnan');
    [~,best_idx] = max(mean_corr_grp);

    axs_lag(gd) = subplot(2,2,gd); hold on
    plot(lag_seconds_grid,mean_corr_grp,'-k')
    plot(lag_seconds_grid(best_idx),mean_corr_grp(best_idx),'ro','MarkerFaceColor','r','MarkerSize',7)
    xline(0,':','Color',[.6,.6,.6])
    title(sprintf('%s (opt=%.2fs, %d frames, r=%.2f)',group_defs(gd).label,lag_seconds_grid(best_idx),lag_frames_grid(best_idx),mean_corr_grp(best_idx)))
    xlabel('lag: bump vel. relative to fly heading vel. (s)')
    ylabel('mean correlation (r), one value per fly')
end
linkaxes(axs_lag,'xy')

% per-fly gain/correlation and its group category, one point per fly
% (concatenating that fly's own trials in this group -- see fly_gain_vectors)
fly_gain     = [];
fly_corr_vel = [];
cat_x_fly    = [];

for gd = 1:numel(group_defs)
    these_flies = group_fly_list{gd};
    trial_lists = group_fly_trials{gd};
    n_f = numel(these_flies);
    lag = group_lag_frames(gd); % this group's own optimal lag, from the sweep above

    figure(4+gd); clf
    set(gcf,'Name',sprintf('gain: %s',group_defs(gd).label),'Position',[50,50,900,700])
    cols = max(ceil(sqrt(n_f)),1);
    rows_n = max(ceil(n_f/cols),1);
    tl = tiledlayout(rows_n,cols,'TileSpacing','compact','Padding','compact');

    for ff = 1:n_f
        nexttile(tl); hold on

        [fly_vel,bump_vel,valid] = fly_gain_vectors(all_data,trial_lists{ff},lag,vel_thresh,bump_thresh,rho_thresh,vel_max);

        scatter(fly_vel(valid),bump_vel(valid),8,[0,0.4470,0.7410],'filled','MarkerFaceAlpha',.15)
        plot(scatter_xlim,[0,0],':','Color',[.6,.6,.6])
        plot([0,0],scatter_ylim,':','Color',[.6,.6,.6])

        if sum(valid) > 1
            b = [ones(sum(valid),1),fly_vel(valid)] \ bump_vel(valid);
            r = corr(fly_vel(valid),bump_vel(valid));
            plot(scatter_xlim,scatter_xlim*b(2)+b(1),'r','LineWidth',1.5)
            text(scatter_xlim(2),scatter_ylim(1),sprintf('gain=%.2f\nr=%.2f',b(2),r), ...
                'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7)
            fly_gain(end+1)     = b(2); %#ok<AGROW>
            fly_corr_vel(end+1) = r; %#ok<AGROW>
            cat_x_fly(end+1)    = gd; %#ok<AGROW>
        end

        xlim(scatter_xlim); ylim(scatter_ylim)
        title(fly_short_label(fly_list{these_flies(ff)}),'FontSize',7,'Interpreter','none')
    end

    xlabel(tl,'fly heading velocity (rad/s)')
    ylabel(tl,'bump velocity (rad/s)')
    title(tl,sprintf('%s -- bump velocity vs. fly heading velocity (n=%d flies, lag=%d frames/%.2fs)',group_defs(gd).label,n_f,lag,lag*mean(trial_dt,'omitnan')),'Interpreter','none')
end

%% figure: gain and correlation coefficient by genotype x light condition, one point per fly
% same group-scatter style as lpsp_compartments_claude_script.m's own
% groupplot helper (jittered per-item points + mean +/- SEM errorbar).
% group_defs is ordered genotype-major (empty CL, empty dark, lpsp CL,
% lpsp dark), so each genotype's CL/dark pair shares one color.
cat_labels  = {group_defs.label};
base_colors = lines(numel(geno_order));
cat_colors  = zeros(numel(group_defs),3);
for k = 1:numel(geno_order)
    cat_colors(2*k-1,:) = base_colors(k,:);
    cat_colors(2*k,:)   = base_colors(k,:);
end

figure(9); clf
set(gcf,'Name','gain & correlation by genotype x light condition','Position',[100,100,900,700])
subplot(2,1,1)
groupplot(cat_x_fly,fly_gain,cat_labels,cat_colors)
ylabel('gain (bump vel / fly heading vel)')
title('bump velocity gain, by genotype and light condition (one point per fly)')

subplot(2,1,2)
groupplot(cat_x_fly,fly_corr_vel,cat_labels,cat_colors)
ylabel('correlation coefficient (r)')
title('bump vs. fly heading velocity correlation, by genotype and light condition (one point per fly)')

%% bump position tracking accuracy: sliding-window circular correlation and circ_var(circ_dist()) between bump position (im.mu) and fly heading (-ft.cue)
% same sliding-window analysis as lpsp_kir_claude.m / lpsp_compartments_claude_script.m's
% "bump position tracking fly heading" section (part 1: position, part 3:
% accuracy). within each window_s-second window (stepped every
% window_step_s seconds), restricted to "discernable bump" samples (fly
% rotating faster than rot_thresh_track, bump concentration im.rho above
% rho_thresh_track, and the bump's own angular velocity not a
% gradient/unwrap artifact), one circ_corrcc(mu,-cue) and one
% circ_var(circ_dist(mu,-cue)) is computed per window. bump signals (mu,
% rho) are shifted by that trial's own GROUP's optimal lag
% (group_lag_frames(cat_x(i)), from the sweep above) later than behavior
% (r_speed, cue).
rot_thresh_track = 0.5;  % rad/s -- "fly is rotating"
rho_thresh_track = 0.2;  % im.rho -- "bump is discernable"
bump_vel_thresh  = 10;   % rad/s -- exclude gradient/unwrap artifacts at wrap-around points
window_s         = 30;   % seconds -- sliding window length
window_step_s    = 10;   % seconds -- step between window starts
min_window_n     = 50;   % minimum discernable-bump samples within a window to trust it

chunk_window_corr       = cell(1,n_trials);
chunk_window_var        = cell(1,n_trials);
chunk_window_headingvar = cell(1,n_trials);

for i = find(~isnan(is_dark))
    [chunk_window_corr{i},chunk_window_var{i},chunk_window_headingvar{i}] = trial_window_tracking( ...
        all_data(i),group_lag_frames(cat_x(i)),rot_thresh_track,rho_thresh_track,bump_vel_thresh,window_s,window_step_s,min_window_n);
end

%% figure: bump position tracks fly heading -- circular correlation (window mean), one point per fly, by genotype x light condition
fly_pos_corr = [];
fly_cat_x_track = [];
for gd = 1:numel(group_defs)
    rows = find(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies
        wc_f = [chunk_window_corr{rows(fly_num(rows)==f)}];
        if isempty(wc_f)
            continue
        end
        fly_pos_corr(end+1)    = mean(wc_f); %#ok<AGROW>
        fly_cat_x_track(end+1) = gd; %#ok<AGROW>
    end
end

figure(10); clf
set(gcf,'Name','bump position tracks fly heading: circular correlation per fly','Position',[100,100,900,600])
groupplot(fly_cat_x_track,fly_pos_corr,cat_labels,cat_colors)
ylabel('circular correlation (mu vs. -cue, window mean)')
title(sprintf('bump position tracking accuracy, one point per fly (%ds windows, step %ds)',window_s,window_step_s))

%% figure: bump accuracy -- circular variance of (mu - (-cue)) offset, one point per fly, by genotype x light condition
% alongside each fly's offset variance (in that category's color), its own
% heading variance alone is plotted in gray as a reference.
fly_offset_var = [];
fly_heading_var = [];
fly_cat_x_var = [];
for gd = 1:numel(group_defs)
    rows = find(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies
        trial_list = rows(fly_num(rows)==f);
        wv_f = [chunk_window_var{trial_list}];
        wh_f = [chunk_window_headingvar{trial_list}];
        if isempty(wv_f)
            continue
        end
        fly_offset_var(end+1)  = mean(wv_f); %#ok<AGROW>
        fly_heading_var(end+1) = mean(wh_f); %#ok<AGROW>
        fly_cat_x_var(end+1)   = gd; %#ok<AGROW>
    end
end

figure(11); clf
set(gcf,'Name','bump accuracy: circular variance of offset, per fly','Position',[100,100,900,600])
hold on
gray = [.6,.6,.6];
for cIdx = 1:numel(cat_labels)
    y_off = fly_offset_var(fly_cat_x_var==cIdx);
    y_hdg = fly_heading_var(fly_cat_x_var==cIdx);
    if ~isempty(y_off)
        jit = (rand(size(y_off))-.5)*.25;
        scatter(cIdx-0.18+jit,y_off,20,cat_colors(cIdx,:),'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx-0.18,mean(y_off),std(y_off)/sqrt(numel(y_off)),'o','Color',cat_colors(cIdx,:)*.6, ...
            'MarkerFaceColor',cat_colors(cIdx,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    if ~isempty(y_hdg)
        jit = (rand(size(y_hdg))-.5)*.25;
        scatter(cIdx+0.18+jit,y_hdg,20,gray,'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx+0.18,mean(y_hdg),std(y_hdg)/sqrt(numel(y_hdg)),'o','Color',gray*.6, ...
            'MarkerFaceColor',gray*.6,'LineWidth',2,'MarkerSize',7)
    end
end
xticks(1:numel(cat_labels)); xticklabels(cat_labels)
xlim([0.5,numel(cat_labels)+0.5])
y_lims = ylim;
for cIdx = 1:numel(cat_labels)
    text(cIdx,y_lims(1),sprintf('n=%d',sum(fly_cat_x_var==cIdx)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
end
h_off = scatter(nan,nan,20,[0,0,0],'filled');
h_hdg = scatter(nan,nan,20,gray,'filled');
legend([h_off,h_hdg],{'circ.var(heading - bump)','circ.var(heading)'},'Location','eastoutside')
ylabel('circular variance (window mean)')
title('bump accuracy: offset variance (color) vs. heading variance alone (gray), one point per fly')

%% bump occupancy / "stickiness": does the bump cover as much of the circle as the fly's heading does, or does it cluster at a few preferred positions while heading roams more widely?
% two views of the same question, both restricted to each trial's own
% discernable-bump samples (same rot_thresh_track/rho_thresh_track/
% bump_vel_thresh criteria and per-group optimal lag as the
% position-tracking section above, but pooled raw -- not windowed):
%
% 1) a scalar stickiness index, one value per fly per condition:
%      stickiness = 1 - circ_var(mu) / circ_var(heading)
% 2) a 2D histogram of the mu-heading offset vs. fly heading, pooled
%    across every fly in a genotype x light-condition group.
chunk_mu_ok  = cell(1,n_trials);
chunk_cue_ok = cell(1,n_trials);
for i = find(~isnan(is_dark))
    [chunk_mu_ok{i},chunk_cue_ok{i}] = trial_bumpok_samples( ...
        all_data(i),group_lag_frames(cat_x(i)),rot_thresh_track,rho_thresh_track,bump_vel_thresh);
end

%% figure: bump occupancy stickiness index, one point per fly, by genotype x light condition
% each fly's own mu/cue samples are pooled here BEFORE computing any
% concentration/variance statistic, and every quantity below (stickiness
% ratio, concentration, histogram) is computed once per fly -- never by
% pooling raw angles across different flies (mu's zero-point is an
% arbitrary per-fly/per-session convention -- see lpsp_kir_claude.m for the
% two-flies-with-different-offsets argument for why pooling across flies
% here would be invalid).
occ_edges   = -pi:pi/16:pi; % 32 bins, matching the PB's own (per-hemisphere) wedge resolution
occ_centers = occ_edges(1:end-1) + diff(occ_edges)/2;

fly_stickiness   = [];
fly_mu_conc      = []; % per-fly circular concentration (1-circ_var) of mu alone
fly_heading_conc = []; % per-fly circular concentration of heading alone, for reference
fly_mu_hist      = []; % per-fly normalized histogram of wrapped mu (one row per fly), for the overlay figure below
fly_cat_x_stick  = [];
for gd = 1:numel(group_defs)
    rows = find(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies
        trial_list = rows(fly_num(rows)==f);
        mu_f  = cat(1,chunk_mu_ok{trial_list});
        cue_f = cat(1,chunk_cue_ok{trial_list});
        if numel(mu_f) < 100
            continue
        end
        mu_f_wrapped = mod(mu_f+pi,2*pi)-pi;
        fly_stickiness(end+1)   = 1 - circ_var(mu_f)/circ_var(cue_f); %#ok<AGROW>
        fly_mu_conc(end+1)      = 1 - circ_var(mu_f_wrapped); %#ok<AGROW>
        fly_heading_conc(end+1) = 1 - circ_var(cue_f); %#ok<AGROW>
        fly_mu_hist(end+1,:)    = histcounts(mu_f_wrapped,occ_edges,'Normalization','probability'); %#ok<AGROW>
        fly_cat_x_stick(end+1)  = gd; %#ok<AGROW>
    end
end

fprintf('\n=== bump occupancy / stickiness index (1 - circ.var(mu)/circ.var(heading)), per group ===\n');
for gd = 1:numel(group_defs)
    y = fly_stickiness(fly_cat_x_stick==gd);
    fprintf('  %-24s n=%2d flies   mean=%.3f\n', group_defs(gd).label, numel(y), mean(y,'omitnan'));
end

figure(12); clf
set(gcf,'Name','bump occupancy: stickiness index per fly','Position',[100,100,900,600])
groupplot(fly_cat_x_stick,fly_stickiness,cat_labels,cat_colors)
ylabel('stickiness = 1 - circ.var(mu) / circ.var(heading)')
title('bump occupancy: how much less of the circle the bump covers vs. heading, one point per fly')

%% figure: bump-heading offset (mu - heading) vs. fly heading, pooled 2D histogram per group
heading_edges = -pi:pi/16:pi;
offset_edges  = -pi:pi/16:pi;
figure(13); clf
set(gcf,'Name','bump-heading offset vs. fly heading, pooled 2D histogram','Position',[100,100,900,800])
axs_occ = gobjects(1,numel(group_defs));
for gd = 1:numel(group_defs)
    rows = strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark;
    mu_g     = cat(1,chunk_mu_ok{rows});
    cue_g    = cat(1,chunk_cue_ok{rows});
    offset_g = circ_dist(mu_g,cue_g);

    axs_occ(gd) = subplot(2,2,gd); hold on
    histogram2(cue_g,offset_g,heading_edges,offset_edges,'DisplayStyle','tile','Normalization','probability')
    plot([-pi,pi],[0,0],'w:','LineWidth',1) % reference: offset independent of heading
    axis square
    xlabel('fly heading (-cue, rad)')
    ylabel('mu - heading offset (rad)')
    title(sprintf('%s (n=%d samples)',group_defs(gd).label,numel(mu_g)))
end
linkaxes(axs_occ,'xy')

%% figure: circular occupancy -- how uniformly each FLY's own bump position (mu) and heading (-cue) cover the circle, by genotype x light condition
figure(14); clf
set(gcf,'Name','bump occupancy: circular concentration per fly','Position',[100,100,900,600])
hold on
gray = [.6,.6,.6];
for cIdx = 1:numel(cat_labels)
    y_mu  = fly_mu_conc(fly_cat_x_stick==cIdx);
    y_hdg = fly_heading_conc(fly_cat_x_stick==cIdx);
    if ~isempty(y_mu)
        jit = (rand(size(y_mu))-.5)*.25;
        scatter(cIdx-0.18+jit,y_mu,20,cat_colors(cIdx,:),'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx-0.18,mean(y_mu),std(y_mu)/sqrt(numel(y_mu)),'o','Color',cat_colors(cIdx,:)*.6, ...
            'MarkerFaceColor',cat_colors(cIdx,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    if ~isempty(y_hdg)
        jit = (rand(size(y_hdg))-.5)*.25;
        scatter(cIdx+0.18+jit,y_hdg,20,gray,'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx+0.18,mean(y_hdg),std(y_hdg)/sqrt(numel(y_hdg)),'o','Color',gray*.6, ...
            'MarkerFaceColor',gray*.6,'LineWidth',2,'MarkerSize',7)
    end
end
xticks(1:numel(cat_labels)); xticklabels(cat_labels)
xlim([0.5,numel(cat_labels)+0.5])
y_lims = ylim;
for cIdx = 1:numel(cat_labels)
    text(cIdx,y_lims(1),sprintf('n=%d',sum(fly_cat_x_stick==cIdx)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
end
h_mu  = scatter(nan,nan,20,[0,0,0],'filled');
h_hdg = scatter(nan,nan,20,gray,'filled');
legend([h_mu,h_hdg],{'bump concentration','heading concentration'},'Location','eastoutside')
ylabel('circular concentration, 1-circ.var (0=uniform, 1=one point)')
title('bump occupancy: circular concentration (color) vs. heading concentration alone (gray), one point per fly')

%% figure: per-fly bump-position histograms overlaid (not pooled), by genotype x light condition
figure(15); clf
set(gcf,'Name','bump position histograms, per fly (overlaid, not pooled)','Position',[100,100,900,800])
axs_ind = gobjects(1,numel(group_defs));
for gd = 1:numel(group_defs)
    rows_f = fly_cat_x_stick==gd;
    axs_ind(gd) = subplot(2,2,gd); hold on
    plot(occ_centers,fly_mu_hist(rows_f,:)','Color',[cat_colors(gd,:),.3],'LineWidth',1);
    plot(occ_centers,mean(fly_mu_hist(rows_f,:),1),'Color',cat_colors(gd,:)*.5,'LineWidth',2.5)
    xlim([-pi,pi])
    xlabel('bump position, mu (rad)')
    ylabel('fraction of samples (this fly)')
    title(sprintf('%s (n=%d flies)',group_defs(gd).label,sum(rows_f)))
end
linkaxes(axs_ind,'xy')

%% bump mobility: ratio of bump path length to fly heading path length, during walking bouts
% same "path-length gain" analysis as lpsp_p2x2_walking_script.m / lpsp_kir_claude.m:
% within each sustained-turning "walking bout" (fly_speed = |r_speed|,
% smoothed, above turn_thresh, with brief gaps bridged and short bouts
% dropped), the bump's path length (sum |diff(mu)|, smoothed) and the
% fly's own turning path length (integral of |r_speed| dt) are each summed
% per bout, then related by a weighted (by bout duration), through-origin
% least-squares fit: mov_ratio = bump path length per unit heading path
% length. r_speed (not -cue) is used for the fly's own turning speed, same
% reasoning as lpsp_kir_claude.m/lpsp_p2x2_walking_script.m.
mobility_smooth_window      = 60;   % samples (~1s at this rig's ~60Hz fictrac rate), gaussian smoothing window
mobility_turn_thresh        = 0.25; % rad/s, minimum heading speed to call a bout "walking"
mobility_max_gap_frames     = 30;   % frames of non-walking allowed within a bout before splitting it (0.5s)
mobility_min_walking_frames = 30;   % minimum bout length to keep (0.5s)

chunk_bout_mu    = cell(1,n_trials);
chunk_bout_speed = cell(1,n_trials);
chunk_bout_dur   = cell(1,n_trials);
for i = find(~isnan(is_dark))
    [chunk_bout_mu{i},chunk_bout_speed{i},chunk_bout_dur{i}] = trial_walking_bouts( ...
        all_data(i),mobility_smooth_window,mobility_turn_thresh,mobility_max_gap_frames,mobility_min_walking_frames);
end

n_bouts_per_fly = [];
for gd = 1:numel(group_defs)
    for ff = 1:numel(group_fly_list{gd})
        n_bouts_per_fly(end+1) = sum(cellfun(@numel,chunk_bout_mu(group_fly_trials{gd}{ff}))); %#ok<AGROW>
    end
end

figure(16); clf
set(gcf,'Name','bump mobility: bout-count diagnostic','Position',[100,100,600,450])
histogram(n_bouts_per_fly,'BinMethod','integers')
xlabel('number of walking bouts (this fly, this condition, all its trials pooled)')
ylabel('number of fly x condition entries')
title('how well-powered is each fly''s own mov\_ratio regression?')

%% figures: bump path length vs. heading path length per walking bout, one figure per genotype x light-condition group, one subplot per fly
fly_mob_ratio = [];
cat_x_mob     = [];

for gd = 1:numel(group_defs)
    these_flies = group_fly_list{gd};
    trial_lists = group_fly_trials{gd};
    n_f = numel(these_flies);

    figure(16+gd); clf
    set(gcf,'Name',sprintf('bump mobility: %s',group_defs(gd).label),'Position',[50,50,900,700])
    cols = max(ceil(sqrt(n_f)),1);
    rows_n = max(ceil(n_f/cols),1);
    tl = tiledlayout(rows_n,cols,'TileSpacing','compact','Padding','compact');

    for ff = 1:n_f
        trial_list = trial_lists{ff};
        mov_mu    = cat(1,chunk_bout_mu{trial_list});
        mov_speed = cat(1,chunk_bout_speed{trial_list});
        dur       = cat(1,chunk_bout_dur{trial_list});

        nexttile(tl); hold on
        if numel(mov_mu) >= 1
            w = sqrt(dur);
            ratio = (w.*mov_speed) \ (w.*mov_mu); % weighted, through-origin: bump path length per unit heading path length

            scatter(mov_speed,mov_mu,10+40*dur/max(dur),[0,0.4470,0.7410],'filled','MarkerFaceAlpha',.5)
            xl = [0,max([mov_speed;mov_mu],[],'omitnan')*1.05];
            if ~(diff(xl)>0)
                xl(2) = 1;
            end
            plot(xl,xl*ratio,'r','LineWidth',1.5)
            plot(xl,xl,':k') % reference: bump moves exactly as much as heading (ratio = 1)
            xlim(xl); ylim(xl); axis square
            text(xl(2),0,sprintf('ratio=%.2f',ratio),'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',7)

            fly_mob_ratio(end+1) = ratio; %#ok<AGROW>
            cat_x_mob(end+1)     = gd; %#ok<AGROW>
        else
            axis off
        end
        title(sprintf('%s (n=%d bouts)',fly_short_label(fly_list{these_flies(ff)}),numel(mov_mu)),'FontSize',7,'Interpreter','none')
    end

    xlabel(tl,'heading path length per bout (rad)')
    ylabel(tl,'bump path length per bout (rad)')
    title(tl,sprintf('%s -- bump vs. heading path length per walking bout (n=%d flies)',group_defs(gd).label,n_f),'Interpreter','none')
end

fprintf('\n=== bump mobility (bump path length / heading path length per walking bout), per group ===\n');
for gd = 1:numel(group_defs)
    y = fly_mob_ratio(cat_x_mob==gd);
    fprintf('  %-24s n=%2d flies   mean ratio=%.3f\n', group_defs(gd).label, numel(y), mean(y,'omitnan'));
end

%% figure: bump mobility (path-length ratio), one point per fly, by genotype x light condition
figure(21); clf
set(gcf,'Name','bump mobility: path-length ratio per fly','Position',[100,100,900,600])
groupplot(cat_x_mob,fly_mob_ratio,cat_labels,cat_colors)
hold on
plot(xlim,[1,1],':k') % reference: bump moves exactly as much as heading
ylabel('bump path length / heading path length (per fly, bout-duration-weighted)')
title('bump mobility relative to fly turning, by genotype and light condition')

%% save all figures as PDF
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\tnt';
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
function name = trial_folder_name(meta_path)
    % the trial folder is the last path component of meta, unless meta
    % points into a registration_NNN subfolder, in which case it's the
    % folder one level up.
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        name = parts{end-1};
    else
        name = parts{end};
    end
end

function fid = trial_fly_id(meta_path)
    % meta = ...\<date folder>\fly N\<trial folder>\registration_NNN
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
    if isempty(fly_part)
        fid = meta_path; % shouldn't happen; fall back to a unique id per trial
    else
        fid = strjoin(parts(1:fly_part(1)),filesep); % date folder + "fly N"
    end
end

function [x_cat,alpha,z_cat,mu_cat,heading_cat,bounds] = concat_trials_im(all_data,trial_idx)
    % concatenate the im.z/im.mu/heading traces of one or more trials (for
    % one fly, one light condition) end-to-end along the imaging-frame
    % axis, so a fly with more than one trial in a condition still gets a
    % single panel. concatenation is by frame count, not real elapsed time,
    % since trials aren't necessarily the same duration/frame rate -- this
    % is just a continuous strip to look at, not an aligned timebase.
    %
    % heading (-ft.cue, at fictrac's own sampling rate) is interpolated
    % onto this trial's imaging-frame timebase (xb, same construction used
    % throughout this codebase, e.g. lpsp_compartments_claude_script.m) so
    % it lines up frame-for-frame with im.z/im.mu on the shared x-axis.
    % "alpha" is returned as-is (informational only, e.g. for callers on
    % single-hemisphere datasets) -- callers on hemisphere-resolved
    % datasets (this one) build their own plain wedge-index y-axis instead.
    alpha = all_data(trial_idx(1)).im.alpha(:);
    z_cat = []; mu_cat = []; heading_cat = []; x_cat = []; bounds = [];
    x_offset = 0;
    for k = 1:numel(trial_idx)
        i    = trial_idx(k);
        z    = all_data(i).im.z;
        mu   = all_data(i).im.mu(:);
        n_im = size(z,2);
        x    = (0:n_im-1)' + x_offset;

        xf      = all_data(i).ft.xf;
        xb      = linspace(xf(1),xf(end),n_im)';
        heading = interp1(xf,unwrap(-all_data(i).ft.cue),xb);
        heading = mod(heading,2*pi);
        heading(heading > pi) = heading(heading > pi) - 2*pi;

        z_cat       = [z_cat, z]; %#ok<AGROW>
        mu_cat      = [mu_cat; mu]; %#ok<AGROW>
        heading_cat = [heading_cat; heading]; %#ok<AGROW>
        x_cat       = [x_cat; x]; %#ok<AGROW>
        x_offset = x(end) + 1;
        bounds(end+1) = x(end); %#ok<AGROW>
    end
end

function row = theta_to_wedge_row(theta,n_hemi)
    % maps a circular angle (radians, any wrapping) onto a continuous
    % [1,n_hemi+1) wedge-index coordinate matching how im.alpha's first
    % n_hemi entries are laid out (alpha(1)=-pi, evenly spaced up to just
    % under +pi) -- used to overlay a single circular value (bump position
    % or heading) on the plain wedge-index y-axis used for this dataset's
    % hemisphere-stacked PB heatmap panels (see header comment for why a
    % real-angle y-axis doesn't work here). the line is broken at circular
    % wraps (a jump of more than half the wedge range) so it doesn't draw
    % a spurious vertical connecting the top and bottom of a band.
    row = 1 + n_hemi*mod(theta+pi,2*pi)/(2*pi);
    row(find(abs(diff(row))>n_hemi/2)+1) = nan;
end

function lbl = fly_short_label(fly_path)
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function lbl = trial_short_label(meta_path)
    % "<date>-<trialnum>" recovered from the trial folder name, e.g.
    % "20240109-1_epg_7f_lpsp_tnt" -> "20240109-1"
    tname = trial_folder_name(meta_path);
    m = regexp(tname,'^(\d{8}-\d+)','match','once');
    if isempty(m)
        lbl = tname;
    else
        lbl = m;
    end
end

function [fly_vel,bump_vel,valid] = trial_gain_vectors(trial,lag,vel_thresh,bump_thresh,rho_thresh,vel_max)
    % fly_vel/bump_vel: angular velocity (rad/s) of the fly's heading
    % (-ft.cue) and of the fitted bump position (im.mu), both on the
    % fictrac timebase (ft.xf). bump_vel is shifted `lag` frames later than
    % fly_vel (positive lag: the bump lags behind the fly's own turning;
    % lag can also be zero or negative -- same shift-and-trim convention as
    % lpsp_compartments_claude_script.m's own lagged_corr, needed here so
    % the optimal-lag sweep can search both directions). valid flags
    % points where the fly is turning fast enough to be informative, the
    % bump velocity estimate isn't a degenerate outlier, and the bump is
    % concentrated enough (rho) to trust.
    xf   = trial.ft.xf;
    fr   = mean(diff(xf));
    n_im = numel(trial.im.mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    fly_vel_full  = gradient(-trial.ft.cue)/fr;
    bump_vel_full = gradient(interp1(xb,unwrap(trial.im.mu),xf))/fr;
    rho_full      = interp1(xb,trial.im.rho,xf);

    if lag == 0
        fly_vel  = fly_vel_full;
        bump_vel = bump_vel_full;
        rho      = rho_full;
    elseif lag > 0
        fly_vel  = fly_vel_full(1:end-lag);
        bump_vel = bump_vel_full(lag+1:end);
        rho      = rho_full(lag+1:end);
    else
        fly_vel  = fly_vel_full(-lag+1:end);
        bump_vel = bump_vel_full(1:end+lag);
        rho      = rho_full(1:end+lag);
    end

    valid = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh & abs(fly_vel) < vel_max;
end

function [fly_vel,bump_vel,valid] = fly_gain_vectors(all_data,trial_list,lag,vel_thresh,bump_thresh,rho_thresh,vel_max)
    % concatenates trial_gain_vectors across every trial in trial_list (one
    % fly's own trials within a single light condition), so a fly with
    % multiple trials in one condition contributes ONE combined set of
    % points instead of one set per trial -- which would otherwise let
    % that fly outweigh single-trial flies in any downstream fit,
    % correlation, or group comparison.
    fly_vel = []; bump_vel = []; valid = logical([]); % logical(.), not [] -- concatenating a logical with a plain double [] silently promotes the result to double, and indexing fv(valid) with the resulting 0/1 doubles (rather than a logical mask) errors on any 0
    for k = 1:numel(trial_list)
        [fv,bv,vd] = trial_gain_vectors(all_data(trial_list(k)),lag,vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end
end

function groupplot(cat_x, values, cat_labels, colors)
    % jittered per-point scatter + mean +/- SEM errorbar per category,
    % copied verbatim from lpsp_compartments_claude_script.m
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

function [win_corr,win_var,win_heading_var] = trial_window_tracking(trial,lag,rot_thresh,rho_thresh,bump_vel_thresh,window_s,window_step_s,min_window_n)
    % slides a window_s-second window (stepped every window_step_s
    % seconds) across one trial, and within each window restricted to
    % "discernable bump" samples returns one circ_corrcc(mu,-cue) and one
    % circ_var(circ_dist(mu,-cue)) -- see the comment above this function's
    % call site for what each of those two metrics captures. mu/rho are
    % shifted `lag` fictrac frames later than behavior (r_speed/cue), same
    % shift-and-trim convention (supporting zero/negative lag too) as
    % trial_gain_vectors above.
    xf   = trial.ft.xf;
    n_im = numel(trial.im.mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu_full  = interp1(xb,unwrap(trial.im.mu(:)),xf);
    rho_full = interp1(xb,trial.im.rho(:),xf);
    cue_full = -trial.ft.cue;

    if lag == 0
        r_speed_l = trial.ft.r_speed;
        cue_l     = cue_full;
        xf_l      = xf;
        mu_l      = mu_full;
        rho_l     = rho_full;
    elseif lag > 0
        r_speed_l = trial.ft.r_speed(1:end-lag);
        cue_l     = cue_full(1:end-lag);
        xf_l      = xf(1:end-lag);
        mu_l      = mu_full(lag+1:end);
        rho_l     = rho_full(lag+1:end);
    else
        r_speed_l = trial.ft.r_speed(-lag+1:end);
        cue_l     = cue_full(-lag+1:end);
        xf_l      = xf(-lag+1:end);
        mu_l      = mu_full(1:end+lag);
        rho_l     = rho_full(1:end+lag);
    end

    dt = mean(diff(xf_l));
    bump_vel_l = gradient(mu_l)/dt;

    bump_ok = abs(r_speed_l) > rot_thresh & rho_l > rho_thresh & abs(bump_vel_l) < bump_vel_thresh & ...
              ~isnan(mu_l) & ~isnan(cue_l) & ~isnan(rho_l);

    win_corr = []; win_var = []; win_heading_var = [];
    window_starts = xf_l(1):window_step_s:(xf_l(end)-window_s);
    for w = 1:numel(window_starts)
        in_win = xf_l >= window_starts(w) & xf_l < window_starts(w)+window_s;
        idx = in_win & bump_ok;
        if sum(idx) < min_window_n
            continue
        end
        % circ_corrcc's denominator is near-zero (giving NaN) when mu or
        % cue is ~constant relative to its own circular mean within this
        % window -- skip such degenerate windows rather than letting one
        % NaN poison the per-fly mean.
        c = circ_corrcc(mu_l(idx),cue_l(idx));
        if ~isnan(c)
            win_corr(end+1) = c; %#ok<AGROW>
        end
        win_var(end+1)         = circ_var(circ_dist(mu_l(idx),cue_l(idx))); %#ok<AGROW>
        win_heading_var(end+1) = circ_var(cue_l(idx)); %#ok<AGROW>
    end
end

function [mu_ok,cue_ok] = trial_bumpok_samples(trial,lag,rot_thresh,rho_thresh,bump_vel_thresh)
    % same discernable-bump restriction (fly rotating, bump concentrated,
    % bump velocity not a gradient/unwrap artifact) and lag shift-and-trim
    % as trial_window_tracking above, but returns the raw pooled samples
    % (mu, -cue) themselves rather than sliding-window summary statistics
    % -- used for the occupancy/stickiness analysis, where the full
    % distribution of positions visited (not just how well they correlate)
    % is what matters.
    xf   = trial.ft.xf;
    n_im = numel(trial.im.mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu_full  = interp1(xb,unwrap(trial.im.mu(:)),xf);
    rho_full = interp1(xb,trial.im.rho(:),xf);
    cue_full = -trial.ft.cue;

    if lag == 0
        r_speed_l = trial.ft.r_speed;
        cue_l     = cue_full;
        mu_l      = mu_full;
        rho_l     = rho_full;
    elseif lag > 0
        r_speed_l = trial.ft.r_speed(1:end-lag);
        cue_l     = cue_full(1:end-lag);
        mu_l      = mu_full(lag+1:end);
        rho_l     = rho_full(lag+1:end);
    else
        r_speed_l = trial.ft.r_speed(-lag+1:end);
        cue_l     = cue_full(-lag+1:end);
        mu_l      = mu_full(1:end+lag);
        rho_l     = rho_full(1:end+lag);
    end

    bump_vel_l = gradient(mu_l)/mean(diff(xf));

    bump_ok = abs(r_speed_l) > rot_thresh & rho_l > rho_thresh & abs(bump_vel_l) < bump_vel_thresh & ...
              ~isnan(mu_l) & ~isnan(cue_l) & ~isnan(rho_l);

    mu_ok  = mu_l(bump_ok);
    cue_ok = cue_l(bump_ok);
end

function [mov_mu,mov_speed,dur] = trial_walking_bouts(trial,smooth_window,turn_thresh,max_gap_frames,min_walking_frames)
    % detects "walking bouts" (sustained turning) in one trial and returns
    % this trial's own per-bout bump path length (mov_mu) and heading path
    % length (mov_speed, from r_speed), plus each bout's duration -- same
    % method as lpsp_p2x2_walking_script.m's own bump-mobility section.
    xf   = trial.ft.xf;
    dt   = median(diff(xf));
    n_im = numel(trial.im.mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu = interp1(xb,unwrap(trial.im.mu(:)),xf,'linear','extrap');

    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window); % force column, unlike cue, r_speed's native orientation isn't guaranteed
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    n_bouts   = numel(bout_starts);
    mov_mu    = nan(n_bouts,1);
    mov_speed = nan(n_bouts,1);
    dur       = nan(n_bouts,1);
    for b = 1:n_bouts
        rng          = bout_starts(b):bout_ends(b);
        mov_mu(b)    = sum(abs(diff(mu_smooth(rng))),'omitnan');
        mov_speed(b) = sum(abs(r_speed_smooth(rng)),'omitnan')*dt; % path length = integral of |speed| dt
        dur(b)       = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end
