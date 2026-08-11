%% lpsp_compartments_claude_script
% Fresh analysis combining three previously-processed PB imaging datasets:
%   lpsp_cl      <- .data\lpsp_cl_data_20240206.mat       (Z:\pablo\lpsp_cl\)
%   lpsp_cl_redo <- .data\lpsp_cl_redo_data_20240306.mat  (Z:\pablo\lpsp_cl_redo\)
%   epg_dlight   <- .data\epg_dlight_20260417.mat         (Z:\pablo\epg_dlight\)
%
% Step 1: load all three into one combined all_data struct array, tagging
% every trial with which dataset/source file it came from.
%
% Step 2: check that every trial has the fields a standard PB-bump analysis
% needs, and determine closed-loop vs. dark for every trial. If
% all_data(i).ft.pattern is already stored (epg_dlight), that's used
% directly. Otherwise (lpsp_cl, lpsp_cl_redo) it's read from
% trialSettings.csv on the raw data drive, found via all_data(i).meta --
% not from folder names, which turn out to disagree with the logged
% pattern for a real subset of lpsp_cl trials (see report below).

%% load the three source datasets, tagging each trial with its origin
data_dir = fullfile('.data');

dataset_names = {'lpsp_cl','lpsp_cl_redo','epg_dlight'};
source_files  = {'lpsp_cl_data_20240206.mat','lpsp_cl_redo_data_20240306.mat','epg_dlight_20260417.mat'};
base_dirs     = {'Z:\pablo\lpsp_cl\','Z:\pablo\lpsp_cl_redo\','Z:\pablo\epg_dlight\'};

datasets = cell(1,numel(dataset_names));
for k = 1:numel(dataset_names)
    assert(isfolder(base_dirs{k}), 'cannot find %s -- check drive mapping', base_dirs{k})

    tmp = load(fullfile(data_dir,source_files{k}),'all_data');
    d = tmp.all_data(:);
    for i = 1:numel(d)
        d(i).dataset     = dataset_names{k};
        d(i).source_file = source_files{k};
    end
    datasets{k} = d;
    fprintf('loaded %-12s n=%3d trials from %s\n', dataset_names{k}, numel(d), source_files{k});
end

all_data = combine_datasets(datasets{:});
n_trials = numel(all_data);
fprintf('\ncombined all_data: %d trials total\n', n_trials);

%% index trialSettings.csv under each dataset's base dir, for trials that need it
% only lpsp_cl and lpsp_cl_redo trials ever fall back to this (epg_dlight
% trials already have ft.pattern stored -- see loop below), but the index
% is built for all three as a safety net in case a trial is missing it.
% one recursive dir() per dataset (not per trial) so this stays fast even
% though it walks the whole raw data tree on the network drive.
csv_indexes = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(dataset_names)
    fprintf('indexing trialSettings.csv under %s ...\n', base_dirs{k});
    csv_indexes(dataset_names{k}) = build_csv_index(base_dirs{k});
end

%% determine closed-loop vs. dark per trial, and check field completeness
ft_required = {'xf','f_speed','r_speed','cue'};
im_required = {'mu','rho','z','d','f','alpha'};

is_dark        = false(1,n_trials);
condition_src  = cell(1,n_trials);  % how is_dark was determined, for the report below
pattern_str    = cell(1,n_trials);
panels_mode    = cell(1,n_trials);
genotype       = cell(1,n_trials);
name_mismatch  = false(1,n_trials); % folder name and logged pattern disagree on dark/CL (lpsp_cl only)
odd_panels     = false(1,n_trials); % panelsMode isn't a clean "closed loop" (e.g. a pilot mixing modes)
is_ambiguous   = false(1,n_trials); % multiple trialSettings.csv candidates disagree -- needs manual review
is_empty_trial = false(1,n_trials); % trial entry has no meta/ft/im at all (leftover placeholder)
ambiguous_info = cell(1,n_trials);  % for ambiguous trials only: {path,pattern,is_dark} per candidate, for the report
missing_ft     = cell(1,n_trials);
missing_im     = cell(1,n_trials);

for i = 1:n_trials
    meta = all_data(i).meta;

    % a few trials in the source .mat files are leftover empty placeholders
    % (meta/ft/im all []) that were never cleaned out during processing --
    % catch those before doing any path parsing on meta.
    is_empty_trial(i) = isempty(meta) || ~(ischar(meta) || (isstring(meta) && isscalar(meta)));
    if is_empty_trial(i)
        condition_src{i} = 'empty trial (no meta/ft/im)';
        missing_ft{i}    = ft_required;
        missing_im{i}    = im_required;
        continue
    end

    folder_says_dark = contains(meta,'_dark','IgnoreCase',true);

    if isfield(all_data(i).ft,'pattern') && ~isempty(all_data(i).ft.pattern)
        % already known from processing -- no need to touch trialSettings.csv
        pattern_str{i}   = char(all_data(i).ft.pattern);
        is_dark(i)       = contains(pattern_str{i},'background','IgnoreCase',true);
        condition_src{i} = 'ft.pattern (already stored)';
        if strcmp(all_data(i).dataset,'lpsp_cl')
            name_mismatch(i) = is_dark(i) ~= folder_says_dark;
        end
    else
        tname = trial_folder_name(meta);
        res = resolve_csv_for_trial(meta, tname, csv_indexes(all_data(i).dataset));

        pattern_str{i} = res.info.pattern;
        panels_mode{i} = res.info.panels_mode;
        genotype{i}    = res.info.genotype;
        odd_panels(i)  = ~strcmp(res.status,'not_found') && ~isempty(res.info.panels_mode) && ...
                          ~contains(res.info.panels_mode,'closed loop','IgnoreCase',true);

        switch res.status
            case 'ambiguous'
                % candidates disagree on pattern/dark -- don't guess
                is_ambiguous(i)  = true;
                is_dark(i)       = folder_says_dark;
                condition_src{i} = 'AMBIGUOUS trialSettings.csv candidates -- used folder name, needs manual review';
                ambiguous_info{i} = struct('path',{res.candidates},'pattern',{cellfun(@(x)x.pattern,res.infos,'UniformOutput',false)});
            case 'not_found'
                is_dark(i)       = folder_says_dark;
                condition_src{i} = 'folder name fallback (no trialSettings.csv found)';
            otherwise % 'exact' or 'fallback' match, and candidates (if >1) agreed
                is_dark(i)       = res.info.is_dark;
                condition_src{i} = sprintf('trialSettings.csv (%s match)', res.status);
                if strcmp(all_data(i).dataset,'lpsp_cl')
                    name_mismatch(i) = is_dark(i) ~= folder_says_dark;
                end
        end
    end

    missing_ft{i} = find_missing_fields(all_data(i).ft, ft_required);
    missing_im{i} = find_missing_fields(all_data(i).im, im_required);
end

complete = cellfun(@isempty,missing_ft) & cellfun(@isempty,missing_im);
resolved = ~is_empty_trial & ~is_ambiguous & ~strcmp(condition_src,'folder name fallback (no trialSettings.csv found)');

%% summary report
fprintf('\n=== per-dataset summary ===\n');
for k = 1:numel(dataset_names)
    idx = strcmp({all_data.dataset},dataset_names{k});
    fprintf(['%-12s n=%3d   empty=%2d   closed loop=%3d   dark=%3d   ', ...
             'condition resolved=%3d/%-3d   ambiguous=%2d   incomplete fields=%3d   name/pattern mismatch=%3d\n'], ...
        dataset_names{k}, sum(idx), sum(idx & is_empty_trial), sum(idx & ~is_dark & ~is_empty_trial), sum(idx & is_dark & ~is_empty_trial), ...
        sum(idx & resolved), sum(idx & ~is_empty_trial), sum(idx & is_ambiguous), sum(idx & ~complete), sum(idx & name_mismatch));
end
fprintf(['%-12s n=%3d   empty=%2d   closed loop=%3d   dark=%3d   ', ...
         'condition resolved=%3d/%-3d   ambiguous=%2d   incomplete fields=%3d   name/pattern mismatch=%3d\n'], ...
    'TOTAL', n_trials, sum(is_empty_trial), sum(~is_dark & ~is_empty_trial), sum(is_dark & ~is_empty_trial), ...
    sum(resolved), sum(~is_empty_trial), sum(is_ambiguous), sum(~complete), sum(name_mismatch));

fprintf('\n=== genotype / sensor (from trialSettings.csv expName, when read) ===\n');
g = genotype(~cellfun(@isempty,genotype));
[ug,~,ug_idx] = unique(g);
counts = accumarray(ug_idx(:),1);
for k = 1:numel(ug)
    fprintf('  %-30s n=%d\n', ug{k}, counts(k));
end

fprintf('\n=== empty placeholder trials (no meta/ft/im -- likely leftover from processing) ===\n');
for i = find(is_empty_trial)
    fprintf('  [%-12s] trial index %d in %s\n', all_data(i).dataset, i, all_data(i).source_file);
end

fprintf('\n=== AMBIGUOUS trials: multiple trialSettings.csv candidates disagree -- needs manual review ===\n');
for i = find(is_ambiguous)
    fprintf('  [%-12s] %s\n', all_data(i).dataset, meta_display(all_data(i).meta));
    for c = 1:numel(ambiguous_info{i}.path)
        fprintf('      candidate: %s\n          pattern=%s\n', meta_display(ambiguous_info{i}.path{c}), ambiguous_info{i}.pattern{c});
    end
end
if ~any(is_ambiguous)
    fprintf('  (none)\n');
end

fprintf('\n=== trials with no trialSettings.csv found at all (fell back to folder-name for CL/dark) ===\n');
cl_dark_label = {'closed loop','dark'};
for i = find(~is_empty_trial & ~is_ambiguous & strcmp(condition_src,'folder name fallback (no trialSettings.csv found)'))
    fprintf('  [%-12s] %s  -> %s (folder name)\n', all_data(i).dataset, meta_display(all_data(i).meta), cl_dark_label{is_dark(i)+1});
end

fprintf('\n=== trials matched via ficTracData-filename fallback (trial folder name did not match a raw folder directly) ===\n');
for i = find(strcmp(condition_src,'trialSettings.csv (fallback match)'))
    fprintf('  [%-12s] %s\n', all_data(i).dataset, meta_display(all_data(i).meta));
end

fprintf('\n=== trials where folder name and logged pattern disagree on dark/closed-loop ===\n');
for i = find(name_mismatch)
    fprintf('  [%-12s] %s\n      pattern=%s -> using %s (folder name would have said %s)\n', ...
        all_data(i).dataset, meta_display(all_data(i).meta), pattern_str{i}, cl_dark_label{is_dark(i)+1}, cl_dark_label{~is_dark(i)+1});
end

fprintf('\n=== trials with an unexpected panelsMode (not a clean "closed loop") ===\n');
for i = find(odd_panels)
    fprintf('  [%-12s] %s  panelsMode="%s"\n', all_data(i).dataset, meta_display(all_data(i).meta), panels_mode{i});
end

fprintf('\n=== trials missing expected ft/im fields (excluding empty placeholder trials above) ===\n');
for i = find(~complete & ~is_empty_trial)
    fprintf('  [%-12s] %s\n      missing ft: {%s}   missing im: {%s}\n', ...
        all_data(i).dataset, meta_display(all_data(i).meta), strjoin(missing_ft{i},', '), strjoin(missing_im{i},', '));
end

%% indicator (sensor) label per trial, from the folder name
% every dataset spells out its sensor/indicator in the raw folder name, so
% this is derived independently of the trialSettings.csv/ft.pattern
% machinery above (and works uniformly for epg_dlight, which never touches
% trialSettings.csv at all now that ft.pattern is used directly).
indicator_order = {'syt7f','syt8m','GRAB(DA2m)','dLight'};
indicator = cell(1,n_trials);
for i = 1:n_trials
    if is_empty_trial(i)
        continue
    end
    indicator{i} = trial_indicator(all_data(i).meta);
end

fprintf('\n=== indicator counts (closed loop / dark) ===\n');
for k = 1:numel(indicator_order)
    idx = strcmp(indicator,indicator_order{k});
    fprintf('  %-12s n=%3d   closed loop=%3d   dark=%3d\n', indicator_order{k}, sum(idx), sum(idx & ~is_dark), sum(idx & is_dark));
end
n_unlabeled = sum(~is_empty_trial & cellfun(@isempty,indicator));
if n_unlabeled > 0
    fprintf('  %-12s n=%3d  (folder name did not match any known indicator)\n', 'unlabeled', n_unlabeled);
end

%% basic completeness mask, needed below for the fly-level walking check
% (the fuller inclusion masks that build on this are defined further down)
trace_ok = ~is_empty_trial & complete;

%% fly ID per trial, and forward- vs. backward-walking dominance per fly
% lpsp_cl_redo/epg_dlight meta paths have an explicit "fly N" folder, so
% that (plus its date folder, to disambiguate "fly 1" on different days)
% is used directly. lpsp_cl has no such folder; it falls back to the same
% date + trailing-suffix-digit heuristic already used for fly grouping in
% this codebase's lpsp_cl_script_minimal.m -- that heuristic occasionally
% mis-parses (e.g. a trial with no trailing "_N" at all falls back to
% date-only), so it's noted in the printed fly count below for a sanity check.
fly_id = cell(1,n_trials);
for i = find(~is_empty_trial)
    fly_id{i} = trial_fly_id(all_data(i).dataset,all_data(i).meta);
end

[fly_list,~,fly_num] = unique(fly_id(~is_empty_trial));
fly_num_full = nan(1,n_trials);
fly_num_full(~is_empty_trial) = fly_num;

fly_fwd_dominant = false(numel(fly_list),1);
fly_fwd_s = nan(numel(fly_list),1);
fly_bwd_s = nan(numel(fly_list),1);
fly_dataset = cell(numel(fly_list),1);
for f = 1:numel(fly_list)
    trials_f = find(fly_num_full==f & trace_ok);
    fwd_s = 0; bwd_s = 0;
    for i = trials_f
        dt = mean(diff(all_data(i).ft.xf));
        fwd_s = fwd_s + sum(all_data(i).ft.f_speed > 0)*dt;
        bwd_s = bwd_s + sum(all_data(i).ft.f_speed < 0)*dt;
    end
    fly_fwd_dominant(f) = fwd_s > bwd_s;
    fly_fwd_s(f) = fwd_s;
    fly_bwd_s(f) = bwd_s;
    if ~isempty(trials_f)
        fly_dataset{f} = all_data(trials_f(1)).dataset;
    end
end
fprintf('\nflies spending more time walking forward than backward: %d/%d (this is the fly-level inclusion criterion used by analysis_ok below)\n', sum(fly_fwd_dominant), numel(fly_list));
fprintf('excluded flies (more backward- than forward-walking time, summed across all of that fly''s own trials):\n');
for f = find(~fly_fwd_dominant)'
    fprintf('  [%-12s] %s   forward=%.1fs   backward=%.1fs\n', fly_dataset{f}, fly_list{f}, fly_fwd_s(f), fly_bwd_s(f));
end
if all(fly_fwd_dominant)
    fprintf('  (none)\n');
end

trial_fwd_dominant = false(1,n_trials);
trial_fwd_dominant(~is_empty_trial) = fly_fwd_dominant(fly_num);

%% per-trial walking-quality filter: exclude trials where the fly barely walked forward
% tuneable: a trial must spend at least min_forward_time_frac of its
% duration walking forward faster than min_forward_speed, or it's excluded
% (a trial-level filter, unlike the fly-level forward/backward dominance
% check above -- a fly can pass that check overall but still have one
% trial where it was mostly sitting still).
min_forward_speed     = 0.5;  % mm/s
min_forward_time_frac = 0.20; % fraction of the trial's duration

trial_forward_frac = nan(1,n_trials);
for i = find(trace_ok)
    trial_forward_frac(i) = mean(all_data(i).ft.f_speed > min_forward_speed);
end
trial_walks_enough = trial_forward_frac >= min_forward_time_frac;

fprintf('\ntrials spending at least %.0f%% of their duration walking forward faster than %.2g mm/s: %d/%d\n', ...
    min_forward_time_frac*100, min_forward_speed, sum(trial_walks_enough(trace_ok)), sum(trace_ok));

%% trial inclusion mask for the correlation/lag analyses below
% analysis_ok additionally excludes trials that shouldn't go into a
% closed-loop-vs-dark comparison (ambiguous condition, an odd panelsMode
% like the one pilot trial that mixed closed & open loop, a fly that spent
% more time walking backward than forward, or a trial where the fly barely
% walked forward at all)
analysis_ok = trace_ok & ~is_ambiguous & ~odd_panels & ~strcmp(indicator,'') & trial_fwd_dominant & trial_walks_enough;

%% 1) store average and peak PB fluorescence, interpolated onto the fictrac timebase
% average = mean dF/F across all 32 PB wedges (both hemispheres) at each
% imaging frame; peak = max dF/F across those same wedges. both are
% interpolated from the imaging timebase onto ft.xf so they line up
% sample-for-sample with r_speed/f_speed for the correlation analyses below.
for i = find(trace_ok)
    xf   = all_data(i).ft.xf;
    d_im = all_data(i).im.d;
    n_im = size(d_im,2);

    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
        xb = all_data(i).ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end

    avg_im  = mean(d_im,1)';
    peak_im = max(d_im,[],1)';

    all_data(i).ft.fluor_avg  = interp1(xb,avg_im,xf);
    all_data(i).ft.fluor_peak = interp1(xb,peak_im,xf);
end
fprintf('\nstored fluor_avg/fluor_peak traces for %d/%d trials\n', sum(trace_ok), n_trials);

%% 2) find the lag (fluorescence relative to rotational speed) that maximizes their correlation
% correlates |rotational speed| at time t against fluorescence at time
% t+lag, for a range of candidate lags -- matches the sign convention of
% the "abs(fly_vel) vs. amplitude" relationship used throughout this
% codebase's other analysis scripts (bump amplitude tracks turning speed
% regardless of turn direction).
lag_frames_grid = -30:3:180; % ~ -0.5s to +3s in ~50ms steps at this rig's ~60Hz fictrac rate
n_lags = numel(lag_frames_grid);

inc_idx = find(analysis_ok);
n_inc   = numel(inc_idx);
corr_avg  = nan(n_inc,n_lags);
corr_peak = nan(n_inc,n_lags);
trial_dt  = nan(n_inc,1);

for ii = 1:n_inc
    i = inc_idx(ii);
    speed = abs(all_data(i).ft.r_speed);
    fa    = all_data(i).ft.fluor_avg;
    fp    = all_data(i).ft.fluor_peak;
    trial_dt(ii) = mean(diff(all_data(i).ft.xf));

    for L = 1:n_lags
        [corr_avg(ii,L),corr_peak(ii,L)] = lagged_corr(speed,fa,fp,lag_frames_grid(L));
    end
end

dt_mean = mean(trial_dt);
lag_seconds_grid = lag_frames_grid * dt_mean;

mean_corr_avg  = mean(corr_avg,1,'omitnan');
mean_corr_peak = mean(corr_peak,1,'omitnan');
[~,best_idx]      = max(mean_corr_avg);
[~,best_idx_peak] = max(mean_corr_peak);
optimal_lag_frames = lag_frames_grid(best_idx);
optimal_lag_s      = lag_seconds_grid(best_idx);

fprintf('\noptimal lag (avg fluorescence, population mean correlation): %d frames (%.3f s)\n', optimal_lag_frames, optimal_lag_s);
fprintf('optimal lag (peak fluorescence, for comparison): %d frames (%.3f s)\n', lag_frames_grid(best_idx_peak), lag_seconds_grid(best_idx_peak));

figure(1); clf
set(gcf,'Name','optimal lag: fluorescence vs. rotational speed','Position',[100,100,900,700])

subplot(2,1,1); hold on
violin_lag_idx = 1:4:n_lags; % coarser subset of lags -- plotting all n_lags violins would be unreadable
cmap = parula(numel(violin_lag_idx));
for k = 1:numel(violin_lag_idx)
    violin_at(lag_seconds_grid(violin_lag_idx(k)), corr_avg(:,violin_lag_idx(k)), mean(diff(lag_seconds_grid))*1.5, cmap(k,:));
end
plot(xlim,[0,0],':k')
xline(optimal_lag_s,'r--','LineWidth',1.5);
xlabel('lag, fluorescence relative to |rotational speed| (s)')
ylabel({'correlation','(per trial)'})
title(sprintf('distribution across trials at each lag (n=%d trials)',n_inc))

subplot(2,1,2); hold on
stem(lag_seconds_grid,mean_corr_avg,'filled','Color','k')
stem(optimal_lag_s,mean_corr_avg(best_idx),'filled','Color','r','LineWidth',1.5)
plot(xlim,[0,0],':k')
xlabel('lag, fluorescence relative to |rotational speed| (s)')
ylabel({'mean correlation','across trials'})
title(sprintf('optimal lag = %.3f s (%d frames)',optimal_lag_s,optimal_lag_frames))

%% 2b) same lag analysis, broken up by indicator and by light condition
% one panel per indicator x condition group (reusing the same corr_avg
% matrix computed above, just averaged within each group's own trials
% instead of across all of them). violins aren't repeated here -- 8 panels
% of full violins would be unreadable -- so each panel shows the mean
% correlation-vs-lag curve (the bottom-panel style of the figure above)
% for that group only, with that group's own optimal lag marked.
group_indicator = cell(1,n_inc);
group_is_dark   = false(1,n_inc);
for ii = 1:n_inc
    group_indicator{ii} = indicator{inc_idx(ii)};
    group_is_dark(ii)   = is_dark(inc_idx(ii));
end

figure(2); clf
set(gcf,'Name','optimal lag, by indicator and light condition','Position',[100,100,1100,900])

cond_label = {'closed loop','dark'};
y_max = 0;
axs = gobjects(numel(indicator_order),2);
for r = 1:numel(indicator_order)
    for c = 1:2
        rows = strcmp(group_indicator,indicator_order{r}) & (group_is_dark==(c==2));
        axs(r,c) = subplot(numel(indicator_order),2,2*(r-1)+c); hold on
        n_grp = sum(rows);
        if n_grp == 0
            title(sprintf('%s, %s (n=0)',indicator_order{r},cond_label{c}))
            continue
        end
        mean_grp = mean(corr_avg(rows,:),1,'omitnan');
        [grp_peak,grp_best] = max(mean_grp);
        stem(lag_seconds_grid,mean_grp,'filled','Color','k')
        stem(lag_seconds_grid(grp_best),grp_peak,'filled','Color','r','LineWidth',1.5)
        plot(xlim,[0,0],':k')
        title(sprintf('%s, %s (n=%d, opt=%.2fs)',indicator_order{r},cond_label{c},n_grp,lag_seconds_grid(grp_best)))
        y_max = max(y_max,max(mean_grp));
        if r == numel(indicator_order)
            xlabel('lag (s)')
        end
        if c == 1
            ylabel('mean corr.')
        end
    end
end
linkaxes(axs(:),'x')
set(axs(:),'YLim',[-0.1,y_max*1.15])

%% indicator-level optimal lag (pooling closed loop + dark), for the group comparison below
% one lag per indicator, not per indicator x condition group, so that a
% CL-vs-dark difference in the plot below reflects the condition and not a
% different lag choice between its two bars.
indicator_best_idx = nan(1,numel(indicator_order));
for r = 1:numel(indicator_order)
    rows = strcmp(group_indicator,indicator_order{r});
    if ~any(rows)
        continue
    end
    m = mean(corr_avg(rows,:),1,'omitnan');
    [~,indicator_best_idx(r)] = max(m);
end

fprintf('\noptimal lag per indicator (avg fluorescence, pooling closed loop + dark):\n');
for r = 1:numel(indicator_order)
    if isnan(indicator_best_idx(r))
        continue
    end
    fprintf('  %-12s %.3f s (%d frames)\n', indicator_order{r}, lag_seconds_grid(indicator_best_idx(r)), lag_frames_grid(indicator_best_idx(r)));
end

%% 3) correlation between |rotational speed| and fluorescence at each indicator's own optimal lag,
%     by closed-loop vs. dark and by indicator
trial_corr_avg  = nan(1,n_trials);
trial_corr_peak = nan(1,n_trials);
for ii = 1:n_inc
    i = inc_idx(ii);
    r = find(strcmp(indicator_order,indicator{i}),1);
    if isempty(r) || isnan(indicator_best_idx(r))
        continue
    end
    trial_corr_avg(i)  = corr_avg(ii,indicator_best_idx(r));
    trial_corr_peak(i) = corr_peak(ii,indicator_best_idx(r));
end

cat_labels = {};
for k = 1:numel(indicator_order)
    cat_labels{end+1} = sprintf('%s (CL)',indicator_order{k});   %#ok<SAGROW>
    cat_labels{end+1} = sprintf('%s (dark)',indicator_order{k}); %#ok<SAGROW>
end
n_cat  = numel(cat_labels);
cat_x  = nan(1,n_trials);
for i = find(analysis_ok)
    ind_num = find(strcmp(indicator_order,indicator{i}));
    if isempty(ind_num); continue; end
    cat_x(i) = 2*(ind_num-1) + is_dark(i) + 1;
end
base_colors = lines(numel(indicator_order));
cat_colors  = zeros(n_cat,3);
for k = 1:numel(indicator_order) % each indicator's CL/dark pair shares one color
    cat_colors(2*k-1,:) = base_colors(k,:);
    cat_colors(2*k,:)   = base_colors(k,:);
end

lag_str_parts = {};
for r = 1:numel(indicator_order)
    if isnan(indicator_best_idx(r))
        continue
    end
    lag_str_parts{end+1} = sprintf('%s=%.2fs',indicator_order{r},lag_seconds_grid(indicator_best_idx(r))); %#ok<SAGROW>
end
lag_str = strjoin(lag_str_parts,', ');

figure(3); clf
set(gcf,'Name','fluorescence vs. rotational speed, by condition and indicator','Position',[100,100,900,700])

subplot(2,1,1)
groupplot(cat_x,trial_corr_avg,cat_labels,cat_colors);
ylabel({'correlation','(avg PB fluorescence vs. |speed|)'})
title({'average fluorescence, each indicator at its own optimal lag',lag_str})

subplot(2,1,2)
groupplot(cat_x,trial_corr_peak,cat_labels,cat_colors);
ylabel({'correlation','(peak PB fluorescence vs. |speed|)'})
title({'peak fluorescence, each indicator at its own optimal lag',lag_str})

%% 4) adjusted R^2 for fluorescence predicted by forward speed, rotational speed, or both -- one fly per line
% same three-model comparison as lpsp_cl_redo_script_minimal.m (forward
% speed alone, |rotational speed| alone, and both jointly, each fit
% against fluorescence with fitlm and compared by adjusted R^2), but
% pooling each fly's own trials together before fitting, so every fly
% contributes exactly one line, not one per trial.
%
% pooling raw timepoints (rather than averaging separately-fit per-trial
% R^2 values) is the right way to combine a fly's trials here: trial
% duration is ~600s for 199/201 analysis_ok trials (2 are 300s), so
% concatenating timepoints and fitting once automatically gives each
% trial equal weight in proportion to how much data it actually has --
% a 300s trial just contributes half as many points, exactly as it
% should, with no separate weighting step needed. averaging R^2 values
% directly wouldn't have a well-defined combined meaning and would need
% its own ad hoc weighting for those 2 shorter trials anyway.
r2_names     = {'average fluorescence','peak fluorescence'};
fluor_fields = {'fluor_avg','fluor_peak'};

for m = 1:2
    figure(3+m); clf
    set(gcf,'Name',sprintf('adjusted R^2: %s predicted by kinematics, one fly per line',r2_names{m}),'Position',[100,100,1100,900])

    axs = gobjects(numel(indicator_order),2);
    panel_R2 = cell(numel(indicator_order),2);

    for rI = 1:numel(indicator_order)
        if isnan(indicator_best_idx(rI))
            continue
        end
        lag = lag_frames_grid(indicator_best_idx(rI));

        for c = 1:2
            rows = find(strcmp(indicator,indicator_order{rI}) & (is_dark==(c==2)) & analysis_ok);
            these_flies = unique(fly_num_full(rows));
            R2 = nan(numel(these_flies),3);

            for ff = 1:numel(these_flies)
                trial_list = rows(fly_num_full(rows)==these_flies(ff));
                pooled_for = []; pooled_fly = []; pooled_amp = [];
                for i = trial_list(:)' % force a row so the loop steps through one trial at a time regardless of trial_list's orientation
                    [for_l,fly_l,amp_l] = apply_lag3(all_data(i).ft.f_speed,all_data(i).ft.r_speed,all_data(i).ft.(fluor_fields{m}),lag);
                    valid = ~isnan(for_l) & ~isnan(fly_l) & ~isnan(amp_l);
                    pooled_for = [pooled_for; for_l(valid)]; %#ok<AGROW>
                    pooled_fly = [pooled_fly; fly_l(valid)]; %#ok<AGROW>
                    pooled_amp = [pooled_amp; amp_l(valid)]; %#ok<AGROW>
                end
                if numel(pooled_amp) < 100
                    continue
                end
                R2(ff,:) = fit_r2_triplet(pooled_fly,pooled_for,pooled_amp);
            end

            panel_R2{rI,c} = R2;
            axs(rI,c) = subplot(numel(indicator_order),2,2*(rI-1)+c);
            plot_r2_group(R2);
            title(sprintf('%s, %s (n=%d flies)',indicator_order{rI},cond_label{c},sum(~isnan(R2(:,1)))))
            if c == 1
                ylabel('adjusted R^2')
            end
        end
    end

    % link/scale y per indicator (row), not globally across all indicators
    % -- see figure 3's comment above for why
    for rI = 1:numel(indicator_order)
        row_all = [panel_R2{rI,1}; panel_R2{rI,2}];
        if isempty(row_all)
            continue
        end
        row_max = max(row_all,[],'all');
        if isnan(row_max) || row_max <= 0
            continue
        end
        linkaxes(axs(rI,:),'y')
        set(axs(rI,:),'YLim',[min(0,min(row_all,[],'all')),row_max*1.1])
    end
end

%% Functions

function combined = combine_datasets(varargin)
    % vertically concatenate struct arrays that don't all share the same
    % top-level fields (e.g. only epg_dlight has a 'gain' field), padding
    % any fields missing from a given dataset with [] so concatenation works.
    all_fields = {};
    for k = 1:numel(varargin)
        all_fields = union(all_fields,fieldnames(varargin{k}),'stable');
    end

    combined = [];
    for k = 1:numel(varargin)
        s = varargin{k}(:);
        missing = setdiff(all_fields,fieldnames(s));
        for m = 1:numel(missing)
            s(1).(missing{m}) = []; % auto-adds this field (as []) to every element of s
        end
        s = orderfields(s,all_fields);
        if isempty(combined)
            combined = s;
        else
            combined = [combined; s]; %#ok<AGROW>
        end
    end
end

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

function idx = build_csv_index(base_dir)
    % index every trialSettings.csv under base_dir by its trial folder
    % name, in one recursive dir() call, so per-trial lookups are cheap.
    % this drive has plenty of byte-identical duplicate copies of the same
    % raw session (e.g. a "to do" staging copy alongside the final "2p
    % data" copy), so multiple hits per name are expected and are resolved
    % later by comparing their logged pattern/settings, not just counted.
    d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
    idx.names = cell(numel(d),1);
    idx.paths = cell(numel(d),1);
    for i = 1:numel(d)
        [~,trial_name] = fileparts(fileparts(d(i).folder)); % strip trailing \csv
        idx.names{i} = trial_name;
        idx.paths{i} = fullfile(d(i).folder,d(i).name);
    end

    % drop calibration/test recordings (e.g. "..._flash_test") -- these
    % share a date+trial-number with real experimental trials often enough
    % to falsely show up as a second, disagreeing candidate during fallback
    % matching, but they were never a real experimental condition
    is_test = cellfun(@(n) contains(n,'flash_test','IgnoreCase',true), idx.names);
    idx.names = idx.names(~is_test);
    idx.paths = idx.paths(~is_test);

    idx.by_name = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(idx.names)
        key = lower(idx.names{i});
        if isKey(idx.by_name,key)
            idx.by_name(key) = [idx.by_name(key), idx.paths(i)];
        else
            idx.by_name(key) = idx.paths(i);
        end
    end
    fprintf('  found %d trialSettings.csv files (%d distinct trial folder names)\n', numel(d), idx.by_name.Count);
end

function res = resolve_csv_for_trial(trial_dir, tname, idx)
    % match a trial to its trialSettings.csv (there can be more than one
    % identical copy on disk -- see build_csv_index). res.status is:
    %   'exact'      - trial folder name matched a raw folder directly
    %   'fallback'   - no direct name match; recovered via the trial's own
    %                  *ficTracData_DAQ* filename, which some "to analyze"
    %                  folders rename losing information the raw folder
    %                  still has (e.g. append an extra digit for a repeated
    %                  registration within one session, or drop a "_dark"
    %                  suffix)
    %   'ambiguous'  - multiple candidate trialSettings.csv were found and
    %                  they don't agree on pattern/dark -- needs a human
    %   'not_found'  - no candidate at all
    res.info = read_trial_condition('');
    res.status = 'not_found';
    res.candidates = {};
    res.infos = {};

    key = lower(tname);
    if isKey(idx.by_name,key)
        candidates = idx.by_name(key);
        used_fallback = false;
    else
        candidates = {};
        used_fallback = true; % only reached if the fallback below finds something
    end

    if isempty(candidates)
        base_id = fictrac_base_id(trial_dir);
        if ~isempty(base_id)
            hit = false(numel(idx.names),1);
            for i = 1:numel(idx.names)
                hit(i) = strcmp(idx.names{i},base_id) || startsWith(idx.names{i},[base_id,'_']);
            end
            candidates = idx.paths(hit);
        end
    end

    res.candidates = candidates;
    if isempty(candidates)
        res.status = 'not_found';
        return
    end

    infos = cell(numel(candidates),1);
    for i = 1:numel(candidates)
        infos{i} = read_trial_condition(candidates{i});
    end
    res.infos = infos;

    agree = true;
    for i = 2:numel(infos)
        if ~isequal(infos{i}.is_dark,infos{1}.is_dark) || ~strcmp(infos{i}.pattern,infos{1}.pattern)
            agree = false;
        end
    end

    res.info = infos{1};
    if ~agree
        res.status = 'ambiguous';
    elseif used_fallback
        res.status = 'fallback';
    else
        res.status = 'exact';
    end
end

function base_id = fictrac_base_id(trial_dir)
    % recover the true "<date>-<trialnum>" session id from the trial's own
    % *ficTracData_DAQ* filename -- this is what the original processing
    % pipeline itself used to find fictrac data for a trial, so it's a more
    % reliable anchor than the (sometimes renamed) "to analyze" folder name.
    base_id = '';
    d = dir(fullfile(trial_dir,'*ficTracData_DAQ*'));
    if isempty(d)
        return
    end
    m = regexp(d(1).name,'^(\d{8}-\d+)','match','once');
    if ~isempty(m)
        base_id = m;
    end
end

function info = read_trial_condition(csv_path)
    info.found       = false;
    info.pattern     = '';
    info.panels_mode = '';
    info.genotype    = '';
    info.is_dark     = nan;

    if isempty(csv_path) || ~isfile(csv_path)
        return
    end
    info.found = true;

    T = readtable(csv_path,'TextType','string');
    vars = T.Properties.VariableNames;

    if ismember('patternPath',vars) && height(T) > 0 && strlength(T.patternPath(1)) > 0
        info.pattern = char(T.patternPath(1));
        info.is_dark = contains(info.pattern,'background','IgnoreCase',true);
    end
    if ismember('panelsMode',vars) && height(T) > 0
        info.panels_mode = char(T.panelsMode(1));
    end
    if ismember('expName',vars) && height(T) > 0
        info.genotype = char(T.expName(1));
    end
end

function missing = find_missing_fields(s, required)
    missing = {};
    for k = 1:numel(required)
        f = required{k};
        if ~isfield(s,f) || isempty(s.(f)) || (isnumeric(s.(f)) && all(isnan(s.(f)(:))))
            missing{end+1} = f; %#ok<AGROW>
        end
    end
end

function str = meta_display(meta_path)
    % keep report lines short: show only the last 3 path components
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    str = strjoin(parts(max(1,end-2):end),filesep);
end

function ind = trial_indicator(meta_path)
    % every dataset spells out its sensor/indicator in the raw folder name
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
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];

    switch dataset_name
        case {'lpsp_cl_redo','epg_dlight'}
            % meta = ...\<date folder>\fly N\<trial folder>\registration_NNN
            fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
            if isempty(fly_part)
                fid = meta_path; % shouldn't happen; fall back to a unique id
            else
                fid = strjoin(parts(1:fly_part(1)),filesep); % date folder + "fly N"
            end
        otherwise % lpsp_cl: no "fly N" folder -- date + trailing suffix digit,
            % same heuristic as lpsp_cl_script_minimal.m's fly_id grouping
            tname    = trial_folder_name(meta_path);
            date_str = regexp(tname,'^\d{8}','match','once');
            suffix   = regexp(tname,'_(\d+)$','tokens','once');
            if isempty(suffix)
                fid = date_str;
            else
                fid = [date_str,'_',suffix{1}];
            end
    end
end

function [ca,cp] = lagged_corr(speed, fluor_avg, fluor_peak, lag)
    % correlate |speed(t)| against fluorescence at t+lag (frames), i.e.
    % fluorescence lagging speed by "lag" frames; same shift-and-trim
    % convention as this codebase's other lag/correlation code.
    if lag == 0
        s_win  = speed;
        fa_win = fluor_avg;
        fp_win = fluor_peak;
    elseif lag > 0
        s_win  = speed(1:end-lag);
        fa_win = fluor_avg(lag+1:end);
        fp_win = fluor_peak(lag+1:end);
    else
        s_win  = speed(-lag+1:end);
        fa_win = fluor_avg(1:end+lag);
        fp_win = fluor_peak(1:end+lag);
    end

    valid = ~isnan(s_win) & ~isnan(fa_win);
    ca = nan;
    if sum(valid) > 100
        ca = corr(s_win(valid),fa_win(valid));
    end

    valid = ~isnan(s_win) & ~isnan(fp_win);
    cp = nan;
    if sum(valid) > 100
        cp = corr(s_win(valid),fp_win(valid));
    end
end

function violin_at(x0, y, width, color)
    y = y(~isnan(y));
    if numel(y) < 5
        return
    end
    [f,yi] = ksdensity(y);
    f = f / max(f) * width;
    patch([x0-f,fliplr(x0+f)],[yi,fliplr(yi)],color,'FaceAlpha',.5,'EdgeColor','none');
    plot([x0-width,x0+width],[mean(y),mean(y)],'-','Color',color*.6,'LineWidth',2);
end

function groupplot(cat_x, values, cat_labels, colors)
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

function r2 = fit_r2_triplet(fly_vel, for_vel, amp)
    % adjusted R^2 for amp ~ for_vel, amp ~ |fly_vel|, and amp ~ [|fly_vel|,for_vel]
    mdl   = fitlm(for_vel,amp);
    r2(1) = mdl.Rsquared.Adjusted;
    mdl   = fitlm(abs(fly_vel),amp);
    r2(2) = mdl.Rsquared.Adjusted;
    mdl   = fitlm([abs(fly_vel),for_vel],amp);
    r2(3) = mdl.Rsquared.Adjusted;
end

function [for_l,fly_l,amp_l] = apply_lag3(for_vel, fly_vel, amp, lag)
    % shift amp by lag frames relative to for_vel/fly_vel (amp lags behind
    % behavior), same convention used throughout this script
    if lag == 0
        for_l = for_vel; fly_l = fly_vel; amp_l = amp;
    elseif lag > 0
        for_l = for_vel(1:end-lag); fly_l = fly_vel(1:end-lag); amp_l = amp(lag+1:end);
    else
        for_l = for_vel(-lag+1:end); fly_l = fly_vel(-lag+1:end); amp_l = amp(1:end+lag);
    end
end

function plot_r2_group(X)
    % one faint line per row (one trial, or one fly if the rows passed in
    % are already pooled across a fly's trials) connecting its R^2 across
    % the three models (forward, rotational, joint), plus mean+/-sem in
    % red -- matches the plotting style in lpsp_cl_redo_script_minimal.m
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
