%% lpsp_kir_claude
% Loads the LPsP>Kir silencing dataset (empty>kir control vs lpsp>kir
% silencing) and determines each trial's genotype and light condition
% (closed loop vs. dark), then reports how many flies of each genotype
% are in the dataset.
%
% Genotype is read directly from the raw folder name recorded in
% all_data.meta -- every trial folder is suffixed "..._LPsP_kir" or
% "..._empty_kir" (see e.g. 20231109-1_EPG_7f_LPsP_kir), the same string
% lpsp_kir_script_minimal.m checked for its own lpsp_idx (contains(meta,'_lpsp_')).
%
% Light condition is NOT stored in this dataset: unlike epg_dlight in
% lpsp_compartments_claude_script.m, no trial here has an all_data.ft.pattern
% field at all (checked directly against lpsp_kir_redo_data_20240206.mat --
% 0/74 trials have it). lpsp_kir_script_minimal.m worked around this with a
% trial-number-parity heuristic (even/odd trial number => dark/closed loop)
% -- fragile if a session ever skips or repeats a trial. Instead this reads
% each trial's own trialSettings.csv directly (same file/field
% lpsp_kir_script_minimal.m stored inline during processing, at
% all_data(i).ft.pattern = tmp2.patternPath{1}).
%
% Note: all_data.meta points into Z:\pablo\stacks\lpsp_kir_redo\... (the
% registered-stack tree), but trialSettings.csv lives under the raw session
% tree, Z:\pablo\lpsp_kir_redo\... -- confirmed these are separate trees,
% not just a path prefix difference (fileparts(meta) does not contain a csv
% folder). So trials are matched to their csv by trial FOLDER NAME against
% a recursive index of base_dir, not by editing meta's own path. This
% dataset has an exact 1:1 trial<->csv correspondence by name (74 trials,
% 74 trialSettings.csv, confirmed directly, no duplicates) -- no need for
% the fuzzy multi-candidate/ambiguity matching lpsp_compartments_claude_script.m
% uses for other datasets with renamed/duplicated raw folders.

%% load data
data_dir    = fullfile('.data');
source_file = 'lpsp_kir_redo_data_20240206.mat';
base_dir    = 'Z:\pablo\lpsp_kir_redo\';
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype per trial, from the raw folder name in all_data.meta
% matched against the trial folder name only (e.g. "20231109-1_EPG_7f_LPsP_kir"),
% not the full meta path -- the base data folder itself is named
% "lpsp_kir_redo", which would otherwise falsely match a plain
% contains(meta,'LPsP_kir') check on every trial, including empty>kir ones.
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
geno_order = {'empty>kir','lpsp>kir'};
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
title(sprintf('lpsp\\_kir\\_redo: flies per genotype (n=%d flies, %d trials)', n_flies, n_trials))
for k = 1:numel(geno_order)
    text(k, fly_counts(k)+0.1, num2str(fly_counts(k)), 'HorizontalAlignment','center')
end

%% figures: PB activity (im.z) with bump position (im.mu) and fly heading overlay, one row per fly, two columns (closed loop | dark)
% two figures, split by genotype (empty>kir / lpsp>kir) since a single
% 33-row figure would be unreadable. within a fly's column, if more than
% one trial shares that light condition (a few flies have 2 CL + 2 dark
% trials -- see the per-fly trial-count check above), those trials are
% concatenated frame-by-frame (not by real elapsed time -- trials aren't
% necessarily the same duration) with a dotted vertical line at each trial
% boundary, so every fly still gets exactly one panel per condition.
%
% heading is -ft.cue (ft.cue is the closed-loop position signal driving the
% panels, which continues to track the fly's own rotation in dark trials
% too -- there's just no visible pattern -- so it's a valid heading proxy
% in both conditions). the sign flip and use of cue as a heading stand-in
% match lpsp_kir_script_minimal.m's own diagnostic overlay
% (plot(ft.xf,-ft.cue,'m')); confirmed directly here too (corr(mu,-cue)=+0.64
% vs corr(mu,cue)=-0.64 on trial 1). heading is colored per the requested
% convention: cyan for closed loop, magenta for dark.
cond_label = {'closed loop','dark'};
cond_color = {[0,1,1],[1,0,1]}; % cyan (closed loop), magenta (dark)
z_clim     = [-2,4]; % shared color scale across every panel (im.z 1st-99th pctile across the dataset is about [-2.4,3.6])

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

            [x_cat,alpha,z_cat,mu_cat,heading_cat,bounds] = concat_trials_im(all_data,trial_idx);
            imagesc(ax,x_cat,alpha,z_cat)
            set(ax,'YDir','normal','CLim',z_clim,'XLim',[x_cat(1),x_cat(end)],'YLim',[alpha(1),alpha(end)])

            mu_plot = mu_cat;
            mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan; % break the line at circular wraps, not connect across them
            plot(ax,x_cat,mu_plot,'w','LineWidth',0.5)

            heading_plot = heading_cat;
            heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;
            plot(ax,x_cat,heading_plot,'Color',cond_color{c},'LineWidth',0.5)

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
    title(t,sprintf('%s -- PB activity (im.z), bump position (im.mu, white), heading (-ft.cue: cyan=CL, magenta=dark)',geno_order{gi}),'Interpreter','none')
end

%% figures: bump velocity vs. fly heading velocity (gain), one figure per genotype x light-condition group, one subplot per FLY
% same gain/regression analysis as lpsp_kir_script_minimal.m's own scatter
% section: fly_vel = angular velocity of heading (gradient of -ft.cue, the
% same heading proxy used for the overlay above), bump_vel = angular
% velocity of the fitted bump position (gradient of unwrapped im.mu,
% interpolated onto the fictrac timebase). bump_vel is offset by a lag
% (in frames) relative to fly_vel (the bump lags the fly's own turning),
% and a point is kept only where the fly is actually turning
% (vel_thresh<|fly_vel|<vel_max), the bump velocity estimate isn't
% degenerate (|bump_vel|<bump_thresh), and bump concentration is high
% enough to trust (rho>rho_thresh) -- identical thresholds to
% lpsp_kir_script_minimal.m.
%
% a handful of flies have 2 trials in one light condition rather than 1
% (see the per-fly trial-count check near the top of this script); those
% trials are concatenated (fly_gain_vectors, below) into one combined set
% of points per fly BEFORE fitting/correlating, so every fly gets exactly
% one panel and contributes exactly one point to every downstream
% comparison -- a fly with 2 trials was previously getting 2 panels and
% counted twice in the lag sweep and in figure 9's group comparison. A
% type-2 (x-on-y... here y-on-x) linear fit (bump_vel ~ 1 + fly_vel) is
% overlaid, with its slope ("gain") and the raw correlation printed per fly.
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
% same "sweep a lag grid, average correlation across a group's own trials,
% take the argmax" approach as lpsp_compartments_claude_script.m's own lag
% search (its section 2/2b: "find the lag that maximizes their
% correlation"), but applied to velocity correlation
% (corr(fly_vel,bump_vel) from fly_gain_vectors) rather than a
% fluorescence-vs-speed correlation, grouped by genotype x light condition
% rather than by indicator, and averaged ONE VALUE PER FLY (not per
% trial) at each candidate lag -- a fly with 2 trials in a condition
% contributes its own single (concatenated) correlation, matching every
% other fly's single contribution, instead of counting twice.
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
% groupplot helper (jittered per-item points + mean +/- SEM errorbar,
% reused verbatim below) -- applied here to the per-fly gain and
% correlation values computed in the section just above. group_defs is
% ordered genotype-major (empty CL, empty dark, lpsp CL, lpsp dark), so
% each genotype's CL/dark pair shares one color, same convention the
% compartments script used for its own indicator x condition groups.
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
% same sliding-window analysis as lpsp_compartments_claude_script.m's own
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
%
% part 1 (position): circ_corrcc measures how well bump position tracks
% heading moment-to-moment; many short windows (rather than one
% trial-long correlation) avoids being fooled by slow relative drift
% between two cumulative/integrated angular signals (mu and cue), even
% when short-timescale tracking is good.
% part 3 (accuracy): circ_var(circ_dist(mu,-cue)) is the circular variance
% of the offset between bump and heading within a window -- how tightly
% clustered that offset is, regardless of its mean value (which is
% arbitrary, since mu's own zero-point is arbitrary). circ_var(-cue) alone
% (the heading's own circular variance within the same window) is kept
% alongside it as a reference: a fly holding a near-constant heading makes
% a low offset variance unremarkable.
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
% heading variance alone is plotted in gray as a reference -- see comment above.
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
%    stickiness ~ 0 means the bump covers about as much of the circle as
%    heading does (normal tracking, circ_var scales together); stickiness
%    -> 1 means the bump occupies a much narrower range than heading did
%    (bump "stuck" near a limited set of positions while the fly's actual
%    heading kept changing -- an attractor-like state); stickiness < 0
%    would mean the bump moves around MORE than heading does (not expected
%    biologically, but not excluded by the formula, so worth noticing if it
%    ever shows up).
% 2) a 2D histogram of bump position (mu) vs. fly heading (-cue), pooled
%    across every fly in a genotype x light-condition group -- if the bump
%    tracks heading uniformly everywhere around the circle, density should
%    sit along the diagonal at every heading value; "stickiness" localized
%    to particular bump positions shows up as horizontal bands (mu
%    clustering near a fixed value across a wide range of heading) rather
%    than a clean diagonal, which the scalar index above can't distinguish
%    from uniformly-poor tracking.
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
% pooling raw angles across different flies. mu's zero-point is an
% arbitrary per-fly/per-session convention (there's no calibration tying
% it to a shared anatomical reference across flies), so a fly holding a
% rock-steady offset of pi and another holding a rock-steady offset of
% pi/2 are both maximally "stuck" individually, but pooling their raw mu
% values together would split that mass into two bins and make the
% pooled distribution look far LESS concentrated than either fly's own is
% -- exactly backwards. occ_edges/occ_centers (used for the per-fly
% histogram below) are defined here since this is the first place they're needed.
occ_edges   = -pi:pi/16:pi; % 32 bins, matching the PB's own wedge resolution
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

%% figure: bump-heading offset (mu - heading) vs. fly heading, pooled 2D histogram per group -- shows WHERE/how a lack of tracking shows up
% plotting mu directly against heading on fixed circular axes suffers a
% wraparound-seam artifact (values near +pi and -pi are the same physical
% angle but get split to opposite plot edges, producing spurious bright
% corners -- confirmed directly: an earlier version of this figure showed
% exactly that, not a clean diagonal). plotting the (mu-heading) OFFSET
% (circ_dist(mu,cue), so it's already correctly wrapped) against heading
% instead sidesteps that entirely.
%
% if the bump tracked heading with perfect gain everywhere, the offset
% would cluster near one constant value regardless of heading (a flat
% horizontal band). a negative-sloped band instead means the offset
% accumulates across the FULL heading excursion pooled into this plot --
% i.e. mu fails to keep pace with heading over the long run, even if
% (per the gain-scatter figures above) it tracks reasonably well
% instantaneously. that combination -- fine short-timescale velocity gain
% but poor long-timescale position-keeping -- is one signature of a bump
% that keeps drifting back toward a preferred set of positions rather
% than fully integrating every turn, so where along the heading axis a
% band sits/bends is worth inspecting even though this plot alone can't
% distinguish that from a uniformly-reduced gain that just happens to be
% smaller over long excursions than the short-window estimate above.
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
% per-fly circular concentration (1-circ_var, the mean resultant vector
% length: 0 = uniformly spread around the circle, 1 = all mass at one
% point), NOT a pooled-across-flies histogram concentration -- see the
% comment above the per-fly loop that computed fly_mu_conc/
% fly_heading_conc for why pooling raw angles across flies would be
% invalid here (different flies' mu zero-points/preferred offsets don't
% share a common reference). same paired-column style as the bump-accuracy
% figure above: bump concentration in that group's color, heading
% concentration alongside it in gray, one point per fly.
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
% each thin line is ONE fly's own normalized histogram of wrapped mu --
% overlaying rather than pooling means two flies stuck at different
% absolute positions each show up as their OWN peak, instead of being
% averaged into a flatter combined curve. a thick line shows the
% across-flies mean of these per-fly histograms (averaging already-peaky
% distributions, not raw angles, so it doesn't reintroduce the pooling
% problem -- though it can still look flatter than any one fly's own line
% if different flies peak at different positions, which is exactly the
% scenario this figure is designed to make visible rather than hide).
figure(15); clf
set(gcf,'Name','bump position histograms, per fly (overlaid, not pooled)','Position',[100,100,900,800])
axs_ind = gobjects(1,numel(group_defs));
for gd = 1:numel(group_defs)
    rows_f = fly_cat_x_stick==gd;
    axs_ind(gd) = subplot(2,2,gd); hold on
    h = plot(occ_centers,fly_mu_hist(rows_f,:)','Color',[cat_colors(gd,:),.3],'LineWidth',1);
    plot(occ_centers,mean(fly_mu_hist(rows_f,:),1),'Color',cat_colors(gd,:)*.5,'LineWidth',2.5)
    xlim([-pi,pi])
    xlabel('bump position, mu (rad)')
    ylabel('fraction of samples (this fly)')
    title(sprintf('%s (n=%d flies)',group_defs(gd).label,sum(rows_f)))
end
linkaxes(axs_ind,'xy')

%% bump mobility: ratio of bump path length to fly heading path length, during walking bouts
% same "path-length gain" analysis as lpsp_p2x2_walking_script.m's own
% section 2 ("how much does the bump (mu) move relative to how much the
% fly turns?"): within each sustained-turning "walking bout" (fly_speed =
% |r_speed|, smoothed, above turn_thresh, with brief gaps bridged and
% short bouts dropped), the bump's path length (sum |diff(mu)|, smoothed)
% and the fly's own turning path length (integral of |r_speed| dt) are
% each summed per bout, then related by a weighted (by bout duration),
% through-origin least-squares fit: mov_ratio = bump path length per unit
% heading path length. mov_ratio ~ 1 means the bump moves as much as the
% fly turns; < 1 means the bump moves less (a "sticky"/reluctant bump) --
% the same theme as the velocity-gain and occupancy sections above, but
% measured over sustained bouts of real turning rather than instantaneous
% velocity samples or raw circular spread.
%
% r_speed (not -cue) is used for the fly's own turning speed, matching
% lpsp_p2x2_walking_script.m's reasoning: that script's cue can change
% experimentally while the fly isn't walking (a stim-related confound
% specific to that dataset), which would otherwise look like spurious
% "heading" movement -- r_speed is the ball's actual measured rotation
% regardless of what the visual scene is doing, so it's used here for
% methodological consistency even though this dataset doesn't have that
% specific confound.
%
% same per-fly concatenation as the rest of this script: a fly's own
% bouts, pooled across however many trials it has in a light condition,
% feed ONE weighted regression -- a fly with 2 trials contributes exactly
% the same as a fly with 1, rather than 2 independent mov_ratio estimates.
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
% marker area scales with bout duration -- the weight each bout carries in
% that fly's own regression -- same convention as
% lpsp_p2x2_walking_script.m's own diagnostic scatter figure.
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
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\kir';
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

function lbl = fly_short_label(fly_path)
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function lbl = trial_short_label(meta_path)
    % "<date>-<trialnum>" recovered from the trial folder name, e.g.
    % "20231109-1_EPG_7f_LPsP_kir" -> "20231109-1"
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
    % fictrac timebase (ft.xf) -- same computation as
    % lpsp_kir_script_minimal.m's own gain-scatter section. bump_vel is
    % shifted `lag` frames later than fly_vel (positive lag: the bump lags
    % behind the fly's own turning; lag can also be zero or negative --
    % same shift-and-trim convention as lpsp_compartments_claude_script.m's
    % own lagged_corr, needed here so the optimal-lag sweep can search
    % both directions). valid flags points where the fly is turning fast
    % enough to be informative, the bump velocity estimate isn't a
    % degenerate outlier, and the bump is concentrated enough (rho) to trust.
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
