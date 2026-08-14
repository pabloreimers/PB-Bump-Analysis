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

%% flash-frame detection: a bright light passing under the objective swamps the biological signal
% detected from raw fluorescence (im.f, before dF/F normalization), not
% dF/F itself, per the reasoning that a flash should show up most cleanly
% there. a frame is flagged as a flash if it is BOTH:
%   (a) a brightness outlier -- its across-wedge mean exceeds the trial's
%       own median wedge-mean by more than flash_mad_thresh times the
%       median absolute deviation (MAD) of that trial's wedge-means
%   (b) unusually spatially UNIFORM across the 32 wedges -- coefficient of
%       variation (SD/mean across wedges) at or below the trial's median.
%       a real PB bump is spatially localized (some wedges much brighter
%       than others: high CV); a flash adds roughly the same extra light
%       to every wedge, which raises the mean but *dilutes* the relative
%       spread, lowering CV. see the proof figures below for direct
%       evidence this is doing something sensible, not just guessing.
flash_mad_thresh = 8; % tuneable

is_flash       = cell(1,n_trials);
n_flash_frames = zeros(1,n_trials);
for i = find(trace_ok)
    is_flash{i}       = detect_flash_frames(all_data(i).im.f, flash_mad_thresh);
    n_flash_frames(i) = sum(is_flash{i});
end

n_im_total = sum(cellfun(@numel,is_flash(trace_ok)));
fprintf('\n=== flash-frame contamination (raw fluorescence outlier + spatially uniform) ===\n');
fprintf('%d/%d trials have at least one flagged flash frame\n', sum(n_flash_frames(trace_ok)>0), sum(trace_ok));
fprintf('total flagged frames: %d / %d imaging frames (%.2f%%)\n', sum(n_flash_frames(trace_ok)), n_im_total, 100*sum(n_flash_frames(trace_ok))/n_im_total);

fprintf('\nmost-contaminated trials (top 8 by flagged frame count):\n');
[~,order] = sort(n_flash_frames,'descend');
for k = 1:8
    i = order(k);
    if n_flash_frames(i) == 0
        break
    end
    fprintf('  [%-12s] %s   %d/%d frames flagged (%.1f%%)\n', ...
        all_data(i).dataset, meta_display(all_data(i).meta), n_flash_frames(i), numel(is_flash{i}), 100*n_flash_frames(i)/numel(is_flash{i}));
end

%% proof: show the detector is finding bright, spatially-uniform events, not just noise
[~,proof_trial] = max(n_flash_frames);
f_proof    = all_data(proof_trial).im.f;
frame_mean = mean(f_proof,1);
frame_cv   = std(f_proof,0,1) ./ frame_mean;
flagged    = is_flash{proof_trial};

figure(10); clf
set(gcf,'Name',sprintf('flash detection proof: %s',meta_display(all_data(proof_trial).meta)),'Position',[100,100,1100,850])

subplot(3,1,1); hold on
plot(frame_mean,'Color',[.6,.6,.6])
plot(find(flagged),frame_mean(flagged),'.r','MarkerSize',8)
ylabel('mean raw fluorescence (across 32 wedges)')
title(sprintf('%s -- %d/%d imaging frames flagged as flash',meta_display(all_data(proof_trial).meta),sum(flagged),numel(flagged)),'Interpreter','none')
legend({'all frames','flagged flash frames'},'Location','northwest')

subplot(3,1,2); hold on
plot(frame_cv,'Color',[.6,.6,.6])
plot(find(flagged),frame_cv(flagged),'.r','MarkerSize',8)
ylabel({'coefficient of variation','across wedges'})
xlabel('imaging frame')

subplot(3,1,3); hold on
scatter(frame_mean(~flagged),frame_cv(~flagged),8,[.6,.6,.6],'filled')
scatter(frame_mean(flagged),frame_cv(flagged),14,'r','filled')
xlabel('mean raw fluorescence'); ylabel('CV across wedges')
legend({'normal frames','flagged flash frames'})
title('flagged frames occupy the bright + relatively-uniform corner')

% per-wedge profile: one flagged flash frame vs. a normal frame nearby in time
flash_frames_list = find(flagged);
example_flash      = flash_frames_list(round(numel(flash_frames_list)/2));
window_start       = max(1,example_flash-50);
window_end         = min(numel(flagged),example_flash+50);
nearby_normal_rel  = find(~flagged(window_start:window_end),1);
nearby_normal      = nearby_normal_rel + window_start - 1;

figure(11); clf
set(gcf,'Name','flagged flash frame vs. a normal frame: per-wedge profile','Position',[100,100,700,420])
hold on
plot(f_proof(:,example_flash),'-or','LineWidth',1.5)
plot(f_proof(:,nearby_normal),'-ok','LineWidth',1.5)
xlabel('PB wedge (1-32)')
ylabel('raw fluorescence')
legend({sprintf('flagged flash frame %d',example_flash),sprintf('normal frame %d (nearby)',nearby_normal)})
title('a flash adds roughly uniform brightness across every wedge')

%% 1) store average and peak PB fluorescence, interpolated onto the fictrac timebase, excluding flash frames
% average = mean dF/F across all 32 PB wedges (both hemispheres) at each
% imaging frame; peak = max dF/F across those same wedges. flash frames
% are dropped before computing either, and before interpolating onto
% ft.xf, so the behavioral-timebase trace interpolates across a flash's
% duration using the surrounding good frames rather than being
% contaminated by it.
for i = find(trace_ok)
    xf   = all_data(i).ft.xf;
    d_im = all_data(i).im.d;
    n_im = size(d_im,2);

    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
        xb = all_data(i).ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end

    keep = ~is_flash{i};
    if sum(keep) < 2
        fprintf('  warning: trial %d has too few non-flash frames -- keeping all frames\n', i);
        keep = true(size(keep));
    end
    xb_keep = xb(keep);
    d_keep  = d_im(:,keep);

    avg_im  = mean(d_keep,1)';
    peak_im = max(d_keep,[],1)';

    all_data(i).ft.fluor_avg  = interp1(xb_keep,avg_im,xf);
    all_data(i).ft.fluor_peak = interp1(xb_keep,peak_im,xf);
end
fprintf('\nstored fluor_avg/fluor_peak traces for %d/%d trials (flash frames excluded)\n', sum(trace_ok), n_trials);

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

%% 5) diagnostic: find "flat" trials in a given indicator (fluorescence doesn't track rotational speed)
% tuneable: which indicator to inspect, and how low the ROTATIONAL-speed
% model's adjusted R^2 has to be (at that indicator's own optimal lag,
% per-trial -- not fly-pooled) to call a trial "flat". for each flagged
% trial, plots summed dFF, peak dFF, and |rotational speed| together over
% time so you can see what the fly/fluorescence were actually doing.
diagnostic_indicator = 'GRAB(DA2m)';
flat_r2_thresh        = 0.2;

r = find(strcmp(indicator_order,diagnostic_indicator),1);
assert(~isempty(r) && ~isnan(indicator_best_idx(r)), 'no optimal lag available for indicator "%s"', diagnostic_indicator)
lag = lag_frames_grid(indicator_best_idx(r));

diag_trials = find(strcmp(indicator,diagnostic_indicator) & analysis_ok);
diag_r2     = nan(numel(diag_trials),3);
for k = 1:numel(diag_trials)
    i = diag_trials(k);
    [for_l,fly_l,amp_l] = apply_lag3(all_data(i).ft.f_speed,all_data(i).ft.r_speed,all_data(i).ft.fluor_avg,lag);
    valid = ~isnan(for_l) & ~isnan(fly_l) & ~isnan(amp_l);
    if sum(valid) < 100
        continue
    end
    diag_r2(k,:) = fit_r2_triplet(fly_l(valid),for_l(valid),amp_l(valid));
end

flat_idx = diag_trials(diag_r2(:,2) < flat_r2_thresh);
fprintf('\n=== "flat" %s trials (rotational-speed R^2 < %.2f at lag=%.2fs): %d/%d ===\n', ...
    diagnostic_indicator, flat_r2_thresh, lag_seconds_grid(indicator_best_idx(r)), numel(flat_idx), numel(diag_trials));
for k = 1:numel(flat_idx)
    i = flat_idx(k);
    row = diag_trials==i;
    fprintf('  %s   [f=%.3f, r=%.3f, j=%.3f]\n', meta_display(all_data(i).meta), diag_r2(row,1), diag_r2(row,2), diag_r2(row,3));
end

for k = 1:numel(flat_idx)
    i = flat_idx(k);
    row = diag_trials==i;

    xf         = all_data(i).ft.xf;
    xb         = linspace(xf(1),xf(end),size(all_data(i).im.d,2));
    sum_dff_im = sum(all_data(i).im.d,1); % raw, at imaging-frame resolution -- for placing flash markers at their actual frame times
    sum_dff    = interp1(xb',sum_dff_im',xf);
    peak_dff   = all_data(i).ft.fluor_peak;
    r_speed    = abs(all_data(i).ft.r_speed);
    flash_mask = is_flash{i};

    figure(20+k); clf % offset well clear of figures 1-5 and the flash-detection proof figures (10-11) above
    set(gcf,'Name',sprintf('flat trial: %s',meta_display(all_data(i).meta)),'Position',[100,100,1000,400])
    subplot(2,1,1)
    yyaxis left
    plot(xf,sum_dff,'-','Color',[0,0.4470,0.7410]); hold on
    plot(xf,peak_dff,'-','Color',[0.4660,0.6740,0.1880])
    plot(xb(flash_mask),sum_dff_im(flash_mask),'.r','MarkerSize',10)
    ylabel('dFF')
    yyaxis right
    plot(xf,r_speed,'-','Color',[0.6,0.6,0.6])
    ylabel('|rotational speed|')
    xlabel('time (s)')
    legend({'summed dFF','peak dFF','flagged flash frames','|rotational speed|'},'Location','northoutside','Orientation','horizontal')
    title(sprintf('%s [%s]   R^2: forward=%.3f rotational=%.3f joint=%.3f', ...
        meta_display(all_data(i).meta), cond_label{is_dark(i)+1}, diag_r2(row,1), diag_r2(row,2), diag_r2(row,3)),'Interpreter','none')

    subplot(2,1,2)
    imagesc(xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)

    linkaxes(get(gcf,'Children'),'x')
    axis tight
end

%% 6) spatial profile of PB activity: hemisphere-averaged von Mises fit per frame
% averages the two PB hemispheres together -- wedges 1:16 and 17:32 share
% the exact same angular positions (confirmed directly: im.alpha repeats
% the same 16-value sequence twice) -- then fits each frame's resulting
% 16-point profile to a von Mises "bump":
%   y(theta) = baseline + amp*exp(kappa*(cos(theta-mu)-1))
% so y = baseline+amp at the peak (theta=mu) and y = baseline+amp*exp(-2*kappa)
% at the trough (theta=mu+pi).
%
% a per-frame nonlinear fit would be far too slow across ~1.6M frames, so
% mu and kappa are estimated analytically and in closed form, the same
% way this repo's own circ_stats toolbox estimates von Mises parameters:
% mu/kappa come from the resultant vector of the (rectified) activity
% treated as circular weights -- circ_r-style, including the standard
% bias correction for 16 evenly-spaced binned angles -- and kappa from R
% via the Fisher approximation used in circ_stats/circ_kappa.m. once mu
% and kappa are fixed, baseline and amp are an exact 2-parameter linear
% least-squares fit, fully vectorized across every frame of a trial at
% once (see fit_hemisphere_vm below). validated against real data before
% writing this in: fits track the actual profile shape well across a
% range of signal strengths (see the proof figure below).
%
% runs on im.z (z-scored per wedge across the trial), not im.d (dFF) --
% z-scoring puts every wedge on the same scale regardless of that wedge's
% own baseline brightness/expression level, which matters here since the
% fit and the rectification step both implicitly compare magnitudes
% across wedges. sections 1-5 (fluorescence traces, lag/correlation/R^2)
% are unaffected and still use dFF (im.d) -- that wasn't asked to change,
% and this way their already-validated numbers stay as they were.
weak_signal_frac  = 0.1; % tuneable: a frame is excluded if its rectified activity is below this fraction of the trial's own median (see fit_hemisphere_vm for why -- near-zero activity makes the resultant length numerically unstable)
min_active_wedges = 3;   % tuneable: a frame is also excluded unless at least this many of the 16 wedges have positive rectified activity -- guards against a single noisy wedge (the rest negative/zero) looking spuriously "concentrated"; verified directly against real data (a kappa=100 degenerate fit came from exactly 1 active wedge before this guard was added)

vm_mu             = cell(1,n_trials);
vm_kappa          = cell(1,n_trials);
vm_amp            = cell(1,n_trials);
vm_baseline       = cell(1,n_trials);
vm_valid          = cell(1,n_trials);
vm_actual_peak    = cell(1,n_trials);
vm_actual_trough  = cell(1,n_trials);

for i = find(trace_ok)
    [vm_mu{i},vm_kappa{i},vm_baseline{i},vm_amp{i},vm_valid{i},vm_actual_peak{i},vm_actual_trough{i}] = ...
        fit_hemisphere_vm(all_data(i).im.z, all_data(i).im.alpha, is_flash{i}, weak_signal_frac, min_active_wedges);
end

for i = find(trace_ok)
    xf   = all_data(i).ft.xf;
    n_im = numel(vm_valid{i});

    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
        xb = all_data(i).ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end

    keep = vm_valid{i};
    if sum(keep) < 2
        continue
    end
    all_data(i).ft.vm_baseline = interp1(xb(keep),vm_baseline{i}(keep)',xf);
    all_data(i).ft.vm_amp      = interp1(xb(keep),vm_amp{i}(keep)',xf);
    all_data(i).ft.vm_kappa    = interp1(xb(keep),vm_kappa{i}(keep)',xf);

    % actual (not fitted) peak/trough only need a flash to be excluded --
    % they're always well-defined from the data, unlike mu/kappa, which
    % need the extra weak-signal guard above
    keep_actual = ~is_flash{i};
    if sum(keep_actual) < 2
        continue
    end
    all_data(i).ft.vm_actual_peak   = interp1(xb(keep_actual),vm_actual_peak{i}(keep_actual)',xf);
    all_data(i).ft.vm_actual_trough = interp1(xb(keep_actual),vm_actual_trough{i}(keep_actual)',xf);
end
fprintf('\nfit hemisphere-averaged von Mises bump parameters for %d/%d trials (on z-scored activity, im.z)\n', sum(trace_ok), n_trials);

%% proof: the von Mises fit tracks the actual hemisphere-averaged profile shape
proof_i = find(trace_ok,1);
theta_proof  = all_data(proof_i).im.alpha(1:16); theta_proof = theta_proof(:);
z_hemi_proof = (all_data(proof_i).im.z(1:16,:) + all_data(proof_i).im.z(17:32,:))/2;

valid_idx = find(vm_valid{proof_i});
[~,order] = sort(vm_kappa{proof_i}(valid_idx));
example_frames = valid_idx(order(round([0.3,0.6,0.85,0.98]*numel(order))));

figure(12); clf
set(gcf,'Name','von Mises fit proof: hemisphere-averaged PB profile (z-scored)','Position',[100,100,1000,700])
for k = 1:4
    f = example_frames(k);
    subplot(2,2,k); hold on
    plot(theta_proof,z_hemi_proof(:,f),'ok','MarkerFaceColor','k')
    theta_fine = linspace(-pi,pi,200)';
    fit_fine = vm_baseline{proof_i}(f) + vm_amp{proof_i}(f)*exp(vm_kappa{proof_i}(f)*(cos(theta_fine-vm_mu{proof_i}(f))-1));
    plot(theta_fine,fit_fine,'-r','LineWidth',1.5)
    title(sprintf('frame %d: kappa=%.2f amp=%.2f baseline=%.2f',f,vm_kappa{proof_i}(f),vm_amp{proof_i}(f),vm_baseline{proof_i}(f)))
    xlabel('wedge angle (rad)'); ylabel('z-score (hemisphere-avg)')
    legend({'data','von Mises fit'},'Location','best')
end
sgtitle(sprintf('%s -- example frames spanning weak to strong signal',meta_display(all_data(proof_i).meta)),'Interpreter','none')

%% question 1: does activity scale by shifting the whole PB up, or by increasing bump amplitude?
speed_edges = 0:0.2:3;
speed_x     = speed_edges(1:end-1) + diff(speed_edges)/2;
min_bin_n   = 50; % minimum pooled samples required to trust a speed bin

has_vm = arrayfun(@(s) isfield(s.ft,'vm_baseline'), all_data);
inc_trials = find(analysis_ok(:) & has_vm(:));
chunk_speed = cell(numel(inc_trials),1); chunk_peak = cell(numel(inc_trials),1); chunk_trough = cell(numel(inc_trials),1);
chunk_baseline = cell(numel(inc_trials),1); chunk_amp = cell(numel(inc_trials),1); chunk_kappa = cell(numel(inc_trials),1);
for ii = 1:numel(inc_trials)
    i = inc_trials(ii);
    speed      = abs(all_data(i).ft.r_speed);
    baseline_t = all_data(i).ft.vm_baseline;
    amp_t      = all_data(i).ft.vm_amp;
    kappa_t    = all_data(i).ft.vm_kappa;
    peak_t     = baseline_t + amp_t;
    trough_t   = baseline_t + amp_t.*exp(-2*kappa_t);

    valid = ~isnan(peak_t);
    chunk_speed{ii}    = speed(valid);
    chunk_peak{ii}     = peak_t(valid);
    chunk_trough{ii}   = trough_t(valid);
    chunk_baseline{ii} = baseline_t(valid);
    chunk_amp{ii}      = amp_t(valid);
    chunk_kappa{ii}    = kappa_t(valid);
end
pooled_speed    = cat(1,chunk_speed{:});
pooled_peak     = cat(1,chunk_peak{:});
pooled_trough   = cat(1,chunk_trough{:});
pooled_baseline = cat(1,chunk_baseline{:});
pooled_amp      = cat(1,chunk_amp{:});
pooled_kappa    = cat(1,chunk_kappa{:});

binned_peak = nan(size(speed_x)); binned_trough = nan(size(speed_x));
binned_baseline = nan(size(speed_x)); binned_amp = nan(size(speed_x));
for j = 1:numel(speed_x)
    idx = pooled_speed>=speed_edges(j) & pooled_speed<speed_edges(j+1);
    if sum(idx) >= min_bin_n
        binned_peak(j)     = mean(pooled_peak(idx));
        binned_trough(j)   = mean(pooled_trough(idx));
        binned_baseline(j) = mean(pooled_baseline(idx));
        binned_amp(j)      = mean(pooled_amp(idx));
    end
end

% same pooling for the ACTUAL (not model-fitted) peak/trough -- the raw
% max/min across the 16 hemisphere-averaged wedges each frame. these only
% require a flash to be excluded (see fit_hemisphere_vm), so they can use
% a slightly larger set of frames than the fitted quantities above (which
% also drop weak-signal frames), and are pooled from their own field.
has_vm_actual = arrayfun(@(s) isfield(s.ft,'vm_actual_peak'), all_data);
inc_trials_actual = find(analysis_ok(:) & has_vm_actual(:));
chunk_speed_actual = cell(numel(inc_trials_actual),1);
chunk_peak_actual  = cell(numel(inc_trials_actual),1);
chunk_trough_actual = cell(numel(inc_trials_actual),1);
for ii = 1:numel(inc_trials_actual)
    i = inc_trials_actual(ii);
    speed    = abs(all_data(i).ft.r_speed);
    peak_a   = all_data(i).ft.vm_actual_peak;
    trough_a = all_data(i).ft.vm_actual_trough;
    valid_a  = ~isnan(peak_a);
    chunk_speed_actual{ii}  = speed(valid_a);
    chunk_peak_actual{ii}   = peak_a(valid_a);
    chunk_trough_actual{ii} = trough_a(valid_a);
end
pooled_speed_actual  = cat(1,chunk_speed_actual{:});
pooled_peak_actual   = cat(1,chunk_peak_actual{:});
pooled_trough_actual = cat(1,chunk_trough_actual{:});

binned_peak_actual = nan(size(speed_x)); binned_trough_actual = nan(size(speed_x));
for j = 1:numel(speed_x)
    idx = pooled_speed_actual>=speed_edges(j) & pooled_speed_actual<speed_edges(j+1);
    if sum(idx) >= min_bin_n
        binned_peak_actual(j)   = mean(pooled_peak_actual(idx));
        binned_trough_actual(j) = mean(pooled_trough_actual(idx));
    end
end

figure(13); clf
set(gcf,'Name','bump peak/trough vs. rotational speed','Position',[100,100,900,950])
subplot(3,1,1); hold on
plot(speed_x,binned_peak,'-o','Color',[0.8500,0.3250,0.0980],'LineWidth',2,'MarkerFaceColor',[0.8500,0.3250,0.0980])
plot(speed_x,binned_trough,'-o','Color',[0,0.4470,0.7410],'LineWidth',2,'MarkerFaceColor',[0,0.4470,0.7410])
legend({'peak (baseline+amp)','trough (baseline+amp*exp(-2kappa))'},'Location','best')
xlabel('|rotational speed| (rad/s)')
ylabel('z-score')
title({'FITTED peak/trough (von Mises model)','parallel lines = whole-PB mean shift.  only the peak line rising = amplitude increase.'})

subplot(3,1,2); hold on
plot(speed_x,binned_peak_actual,'-o','Color',[0.8500,0.3250,0.0980],'LineWidth',2,'MarkerFaceColor',[0.8500,0.3250,0.0980])
plot(speed_x,binned_trough_actual,'-o','Color',[0,0.4470,0.7410],'LineWidth',2,'MarkerFaceColor',[0,0.4470,0.7410])
legend({'actual peak (max of 16 wedges)','actual trough (min of 16 wedges)'},'Location','best')
xlabel('|rotational speed| (rad/s)')
ylabel('z-score')
title('ACTUAL (raw data) peak/trough -- max/min across the hemisphere-averaged wedges, no model fit')

subplot(3,1,3); hold on
plot(speed_x,binned_baseline,'-o','Color',[0.4940,0.1840,0.5560],'LineWidth',2,'MarkerFaceColor',[0.4940,0.1840,0.5560])
plot(speed_x,binned_amp,'-o','Color',[0.4660,0.6740,0.1880],'LineWidth',2,'MarkerFaceColor',[0.4660,0.6740,0.1880])
legend({'fitted baseline','fitted amplitude'},'Location','best')
xlabel('|rotational speed| (rad/s)')
ylabel('z-score')
title('fitted baseline vs. amplitude')

%% question 2: what fraction of the PB is within half-max of the peak?
% a large chunk of frames have kappa too low for the fit to ever drop to
% half-max anywhere around the circle -- half_max_width_frac reports
% these as 100% ("the whole PB counts"), which is mathematically correct
% but really means "no clearly localized bump that frame" (e.g. the fly
% is sitting still -- see figure 13: amplitude collapses at low speed).
% pooling those in with genuine sharp bumps inflates the width estimate,
% so this also reports a version restricted to |rotational speed|>1 rad/s,
% where figure 13 shows amplitude is clearly elevated above baseline.
half_max_frac = half_max_width_frac(pooled_kappa);
frac_no_peak  = mean(pooled_kappa <= log(2)/2);

fprintf('\n=== bump width (fraction of PB within half-max) ===\n');
fprintf('pooled across analysis_ok trials, n=%d frames: median=%.1f%%  mean=%.1f%%\n', ...
    sum(~isnan(half_max_frac)), 100*median(half_max_frac,'omitnan'), 100*mean(half_max_frac,'omitnan'));
fprintf('%.1f%% of those frames have no clearly localized bump at all (kappa too low for any half-max crossing)\n', 100*frac_no_peak);

high_speed_thresh = 1; % rad/s -- see figure 13, amplitude is clearly elevated above baseline past this point
is_high_speed = pooled_speed > high_speed_thresh;
half_max_frac_hs = half_max_width_frac(pooled_kappa(is_high_speed));
fprintf('restricted to |rotational speed| > %g rad/s (n=%d frames, a clear bump is present): median=%.1f%%  mean=%.1f%%\n', ...
    high_speed_thresh, sum(is_high_speed), 100*median(half_max_frac_hs,'omitnan'), 100*mean(half_max_frac_hs,'omitnan'));

trial_median_width = nan(1,n_trials);
for i = find(analysis_ok)
    if ~isfield(all_data(i).ft,'vm_kappa')
        continue
    end
    k = all_data(i).ft.vm_kappa;
    trial_median_width(i) = median(half_max_width_frac(k(~isnan(k))),'omitnan');
end
fprintf('\nbump width by indicator (median-of-trial-medians, %% of PB within half-max, all speeds):\n');
for r = 1:numel(indicator_order)
    idx = strcmp(indicator,indicator_order{r}) & ~isnan(trial_median_width);
    fprintf('  %-12s n=%3d trials   median=%.1f%%\n', indicator_order{r}, sum(idx), 100*median(trial_median_width(idx),'omitnan'));
end

% all-speeds width is dominated by low-speed/low-amplitude moments, where
% kappa is poorly defined for every indicator (amplitude collapses toward
% baseline at rest -- see figure 13), just more severely for some
% indicators than others depending on their typical resting concentration.
% restricting to timepoints where a bump is actually clearly present
% (|speed|>high_speed_thresh) is a fairer per-indicator comparison -- pooled
% across frames, not a median-of-per-trial-medians: |speed|>1 rad/s is a
% fairly rare state within any one 600s trial, so a per-trial median over
% that sparse a subset is itself too noisy to be a useful summary (this
% was tried first and gave unstable, saturated-looking results).
% reuses chunk_kappa/chunk_speed (indexed by inc_trials) from the pooling
% step above rather than re-reading every trial's fields again.
fprintf('\nbump width by indicator (pooled frames, restricted to |speed|>%g rad/s):\n', high_speed_thresh);
for r = 1:numel(indicator_order)
    ind_mask = strcmp(indicator(inc_trials),indicator_order{r});
    k_pool = cat(1,chunk_kappa{ind_mask});
    s_pool = cat(1,chunk_speed{ind_mask});
    hs = s_pool > high_speed_thresh;
    w_ind = half_max_width_frac(k_pool(hs));
    fprintf('  %-12s n=%d frames   median=%.1f%%\n', indicator_order{r}, sum(hs), 100*median(w_ind,'omitnan'));
end

figure(14); clf
set(gcf,'Name','bump width: fraction of PB within half-max','Position',[100,100,700,700])
subplot(2,1,1)
histogram(100*half_max_frac,0:5:100)
xlabel('% of PB within half-max')
ylabel('count (pooled frames, all speeds)')
title(sprintf('all frames: median = %.1f%% of the PB (%.0f%% have no clear bump)',100*median(half_max_frac,'omitnan'),100*frac_no_peak))

subplot(2,1,2)
histogram(100*half_max_frac_hs,0:5:100)
xlabel('% of PB within half-max')
ylabel(sprintf('count (|speed|>%g rad/s)',high_speed_thresh))
title(sprintf('restricted to |rotational speed|>%g rad/s: median = %.1f%% of the PB',high_speed_thresh,100*median(half_max_frac_hs,'omitnan')))

%% figure 14, broken up by indicator and light condition
% reuses chunk_kappa/chunk_speed (indexed by inc_trials) from the pooling
% step above, same 4x2 grid convention as figure 2's lag breakdown.
width_edges = 0:5:100;

figure(15); clf
set(gcf,'Name','bump width by indicator and light condition (all speeds)','Position',[100,100,1100,900])
axs15 = gobjects(numel(indicator_order),2);
for rI = 1:numel(indicator_order)
    for c = 1:2
        mask   = strcmp(indicator(inc_trials),indicator_order{rI}) & (is_dark(inc_trials)==(c==2));
        k_pool = cat(1,chunk_kappa{mask});
        w_pool = half_max_width_frac(k_pool);

        axs15(rI,c) = subplot(numel(indicator_order),2,2*(rI-1)+c);
        histogram(100*w_pool,width_edges,'Normalization','probability')
        title(sprintf('%s, %s (n=%d, median=%.0f%%)',indicator_order{rI},cond_label{c},numel(w_pool),100*median(w_pool,'omitnan')))
        if rI == numel(indicator_order)
            xlabel('% of PB within half-max')
        end
        if c == 1
            ylabel('%frames')
        end
    end
end
linkaxes(axs15(:),'x')

figure(16); clf
set(gcf,'Name',sprintf('bump width by indicator and light condition (|speed|>%g rad/s)',high_speed_thresh),'Position',[100,100,1100,900])
axs16 = gobjects(numel(indicator_order),2);
for rI = 1:numel(indicator_order)
    for c = 1:2
        mask   = strcmp(indicator(inc_trials),indicator_order{rI}) & (is_dark(inc_trials)==(c==2));
        k_pool = cat(1,chunk_kappa{mask});
        s_pool = cat(1,chunk_speed{mask});
        hs     = s_pool > high_speed_thresh;
        w_pool = half_max_width_frac(k_pool(hs));

        axs16(rI,c) = subplot(numel(indicator_order),2,2*(rI-1)+c);
        histogram(100*w_pool,width_edges,'Normalization','probability')
        title(sprintf('%s, %s (n=%d, median=%.0f%%)',indicator_order{rI},cond_label{c},numel(w_pool),100*median(w_pool,'omitnan')))
        if rI == numel(indicator_order)
            xlabel('% of PB within half-max')
        end
        if c == 1
            ylabel('% frames')
        end
    end
end
linkaxes(axs16(:),'x')

%% diagnostic: is the syt8m width result real, or is the width metric misbehaving?
% picks one syt8m trial (the one with the most |speed|>high_speed_thresh
% time, so there's plenty to look at) and plots its raw im.z heatmap next
% to the "% of PB below half-max" trace derived from the same per-frame
% kappa, plus rotational speed, all on a linked time axis -- so periods
% where the width metric says "no localized bump" can be checked directly
% against what the actual PB activity looked like right then.
diag_indicator_2 = 'syt8m';
ind_match2 = strcmp(indicator,diag_indicator_2);
has_kappa2 = arrayfun(@(s) isfield(s.ft,'vm_kappa'), all_data);
is_diag2 = analysis_ok(:) & ind_match2(:) & has_kappa2(:); % force column: mixing row/column here would silently broadcast into an NxN matrix instead of elementwise
n_hs = zeros(1,n_trials);
for i = find(is_diag2)' % transpose: find() on a column vector returns a column, and "for i = <column>" hands the whole column to i in one shot instead of looping element-by-element
    n_hs(i) = sum(abs(all_data(i).ft.r_speed) > high_speed_thresh);
end
[~,diag_trial_2] = max(n_hs);

xf2 = all_data(diag_trial_2).ft.xf;
if isfield(all_data(diag_trial_2).ft,'xb') && numel(all_data(diag_trial_2).ft.xb) == size(all_data(diag_trial_2).im.z,2)
    xb2 = all_data(diag_trial_2).ft.xb(:);
else
    xb2 = linspace(xf2(1),xf2(end),size(all_data(diag_trial_2).im.z,2))';
end

pct_below_hm = 100 - 100*half_max_width_frac(vm_kappa{diag_trial_2});
speed_im2    = interp1(xf2,abs(all_data(diag_trial_2).ft.r_speed),xb2);

figure(17); clf
set(gcf,'Name',sprintf('width diagnostic: %s',meta_display(all_data(diag_trial_2).meta)),'Position',[100,100,1100,800])

ax1 = subplot(3,1,1);
imagesc(xb2,unwrap(all_data(diag_trial_2).im.alpha),all_data(diag_trial_2).im.z)
ylabel('wedge angle (unwrapped, rad)')
title(sprintf('%s -- im.z heatmap',meta_display(all_data(diag_trial_2).meta)),'Interpreter','none')
pos=  get(gca,'Position');
colorbar
set(gca,'Position',pos)

ax2 = subplot(3,1,2);
plot(xb2,pct_below_hm,'-k')
ylabel({'% of PB','below half-max'})
ylim([-5,105])
title('derived from the same per-frame kappa used in figures 14-16')

ax3 = subplot(3,1,3);
plot(xb2,speed_im2,'-','Color',[0.6,0.6,0.6])
ylabel('|rotational speed| (rad/s)')
xlabel('time (s)')

linkaxes([ax1,ax2,ax3],'x')

%% 7) bump position tracking fly heading
% restricted throughout to "discernable bump" moments: the fly rotating
% faster than rot_thresh_track, the stored im.rho (pre-computed resultant
% vector length, from the original processing pipeline -- not the kappa
% from section 6) above rho_thresh_track, and the computed bump velocity
% not itself an unwrap/gradient artifact (|bump_vel|<bump_vel_thresh).
% same three-criteria style as the vel_thresh/rho_thresh/bump_thresh
% filter already used in lpsp_cl_script_minimal.m. flash frames are
% excluded from mu/rho before interpolating onto the fictrac timebase,
% same as sections 1 and 6.
%
% by convention, bump position (mu) is compared against -cue, not raw
% cue (this matches the sign flip already used when plotting cue against
% heading elsewhere in this codebase, e.g. "-unwrap(all_data(i).ft.cue)").
% "cue_t"/"cue_full"/"chunk_cue" below all hold the already-negated value.
%
% mu/rho are fluorescence-derived, so -- like the amplitude-vs-speed
% analyses in sections 2-4 -- they lag true behavior by a sensor-dependent
% kinetic delay. each trial is shifted by ITS OWN indicator's optimal lag
% (indicator_best_idx/lag_frames_grid, from section 2b) before anything
% else is computed: r_speed/cue_t/xf stay on the original clock, mu_t/rho_t
% are trimmed from the other end, same "behavior leads, fluorescence lags"
% convention as apply_lag3 elsewhere in this script.
rot_thresh_track = 0.5;  % rad/s -- "fly is rotating"
rho_thresh_track = 0.2;  % im.rho -- "bump is discernable"
bump_vel_thresh  = 10;   % rad/s -- exclude gradient/unwrap artifacts at wrap-around points

window_s      = 30; % seconds -- sliding window length for the position-correlation analysis (part 1)
window_step_s = 10; % seconds -- step between window starts (60 = non-overlapping)
min_window_n  = 50; % minimum discernable-bump (rotating) samples within a window to trust its correlation

has_cue = arrayfun(@(s) isfield(s.ft,'cue'), all_data);
inc_trials_track = find(analysis_ok(:) & has_cue(:));

chunk_mu         = cell(numel(inc_trials_track),1);
chunk_cue        = cell(numel(inc_trials_track),1);
chunk_bumpvel    = cell(numel(inc_trials_track),1);
chunk_rspeed     = cell(numel(inc_trials_track),1);
chunk_offset     = cell(numel(inc_trials_track),1);
chunk_window_corr = cell(numel(inc_trials_track),1); % one scalar circ_corrcc per valid 60s window in this trial
chunk_window_var  = cell(numel(inc_trials_track),1); % one scalar circ_var(circ_dist(mu,cue)) per valid window in this trial
chunk_window_headingvar = cell(numel(inc_trials_track),1); % one scalar circ_var(cue) per valid window -- reference for how much heading itself varied
n_bump_ok        = zeros(1,n_trials);

for ii = 1:numel(inc_trials_track)
    i = inc_trials_track(ii);
    xf   = all_data(i).ft.xf;
    n_im = size(all_data(i).im.d,2);
    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
        xb = all_data(i).ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end

    mu_im  = unwrap(all_data(i).im.mu(:));
    rho_im = all_data(i).im.rho(:);
    keep   = ~is_flash{i}(:);
    if sum(keep) < 2
        continue
    end

    mu_t  = interp1(xb(keep),mu_im(keep),xf);
    rho_t = interp1(xb(keep),rho_im(keep),xf);

    r_speed = all_data(i).ft.r_speed;
    cue_t   = -all_data(i).ft.cue; % by convention, bump position is compared to -cue, not raw cue (see comment above)

    r_ind = find(strcmp(indicator_order,indicator{i}),1);
    if isempty(r_ind) || isnan(indicator_best_idx(r_ind))
        continue
    end
    lag = lag_frames_grid(indicator_best_idx(r_ind));
    [behav,fluor] = apply_lag_multi({r_speed,cue_t,xf},{mu_t,rho_t},lag);
    r_speed_l = behav{1}; cue_l = behav{2}; xf_l = behav{3};
    mu_l      = fluor{1}; rho_l = fluor{2};

    dt = mean(diff(xf_l));
    bump_vel_l = gradient(mu_l)/dt;

    % ~isnan guards matter here: a NaN in im.mu/im.rho/ft.cue anywhere near
    % the "keep" window propagates through interp1 into a NaN stretch of
    % mu_t/rho_t (seen for real with epg_dlight trials), and neither
    % circ_corrcc nor circ_var skip NaNs on their own -- one NaN sample
    % left in the pool would otherwise silently turn the whole statistic
    % into NaN.
    bump_ok = abs(r_speed_l) > rot_thresh_track & rho_l > rho_thresh_track & abs(bump_vel_l) < bump_vel_thresh & ...
              ~isnan(mu_l) & ~isnan(cue_l) & ~isnan(rho_l);
    n_bump_ok(i) = sum(bump_ok);
    if ~any(bump_ok)
        continue
    end

    chunk_mu{ii}      = mu_l(bump_ok);
    chunk_cue{ii}     = cue_l(bump_ok);
    chunk_bumpvel{ii} = bump_vel_l(bump_ok);
    chunk_rspeed{ii}  = r_speed_l(bump_ok);
    chunk_offset{ii}  = circ_dist(mu_l(bump_ok),cue_l(bump_ok));

    % part 1 (position): slide a window_s-second window across the trial;
    % within each window, restrict to bump_ok (rotating + discernable-bump)
    % samples and compute one circ_corrcc if there are enough of them.
    window_starts = xf_l(1):window_step_s:(xf_l(end)-window_s);
    win_corrs = [];
    win_vars  = [];
    win_heading_vars = [];
    for w = 1:numel(window_starts)
        in_win = xf_l >= window_starts(w) & xf_l < window_starts(w)+window_s;
        idx    = in_win & bump_ok;
        if sum(idx) < min_window_n
            continue
        end
        % circ_corrcc's denominator is a product of two resultant-variance
        % sums; a window where mu or cue is (near-)constant relative to its
        % own circular mean makes that denominator ~0, giving NaN rather than
        % an error -- skip such degenerate windows instead of letting one NaN
        % poison the trial/indicator/pooled mean (same NaN-guard convention
        % used for individual samples elsewhere in this section).
        c = circ_corrcc(mu_l(idx),cue_l(idx));
        if ~isnan(c)
            win_corrs(end+1) = c; %#ok<AGROW>
        end
        % part 3 (accuracy), same windowing: circ_var of the mu-(-cue)
        % offset within this window (circ_var doesn't share circ_corrcc's
        % divide-by-near-zero failure mode, so no NaN guard needed here).
        win_vars(end+1) = circ_var(circ_dist(mu_l(idx),cue_l(idx))); %#ok<AGROW>
        % reference: circ_var of the heading itself within this window --
        % how much the fly's own heading varied, as a scale to compare the
        % offset variance above against (e.g. a fly holding a near-constant
        % heading makes a low offset variance unremarkable).
        win_heading_vars(end+1) = circ_var(cue_l(idx)); %#ok<AGROW>
    end
    chunk_window_corr{ii}       = win_corrs;
    chunk_window_var{ii}        = win_vars;
    chunk_window_headingvar{ii} = win_heading_vars;
end

pooled_mu      = cat(1,chunk_mu{:});
pooled_cue     = cat(1,chunk_cue{:});
pooled_bumpvel = cat(1,chunk_bumpvel{:});
pooled_rspeed  = cat(1,chunk_rspeed{:});
pooled_offset  = cat(1,chunk_offset{:});
all_window_corrs = cat(2,chunk_window_corr{:});

fprintf('\n=== bump position tracking fly heading (n=%d discernable-bump samples, pooled) ===\n', numel(pooled_mu));

% 1) bump position (mu) vs. fly heading (-cue): mean of many short
% (window_s-second) sliding-window circular correlations, rather than one
% circ_corrcc over an entire trial -- mu and cue are both cumulative/
% integrated angular signals, so comparing them over a long window is
% vulnerable to slow relative drift even when short-timescale tracking is
% good (checked directly: an example trial's own offset was tight and
% unimodal even though the whole-trial pooled offset variance was high).
fprintf('1) bump position (mu) vs. fly heading (-cue): mean %ds-window circular correlation = %.3f (n=%d windows)\n', ...
    window_s, mean(all_window_corrs), numel(all_window_corrs));

% 2) bump speed vs. fly rotational speed -- both already linear
% (rad/s, not wrapped), so ordinary Pearson correlation applies.
[vel_r,vel_p] = corr(pooled_bumpvel,pooled_rspeed);
fprintf('2) bump speed vs. fly rotational speed: r = %.3f (p=%.3g)\n', vel_r, vel_p);

% 3) bump accuracy: circular variance (and, for interpretability, circular
% SD in degrees) of the mu-(-cue) offset across discernable-bump moments.
offset_var = circ_var(pooled_offset);
offset_std = circ_std(pooled_offset);
fprintf('3) bump accuracy: circular variance of (mu-(-cue)) offset = %.3f (circular SD = %.1f deg)\n', offset_var, rad2deg(offset_std));

fprintf('\nby indicator:\n');
for r = 1:numel(indicator_order)
    mask   = strcmp(indicator(inc_trials_track),indicator_order{r});
    wc_r   = cat(2,chunk_window_corr{mask});
    bv_r   = cat(1,chunk_bumpvel{mask});
    rs_r   = cat(1,chunk_rspeed{mask});
    off_r  = cat(1,chunk_offset{mask});
    if numel(wc_r) < 5 || numel(bv_r) < 100
        fprintf('  %-12s (too few discernable-bump samples/windows)\n', indicator_order{r});
        continue
    end
    vel_r_r = corr(bv_r,rs_r);
    var_r   = circ_var(off_r);
    fprintf('  %-12s n_windows=%4d   position corr (window mean)=%.2f   speed corr=%.2f   offset circ.var=%.2f\n', ...
        indicator_order{r}, numel(wc_r), mean(wc_r), vel_r_r, var_r);
end

% cue is only a meaningful thing to track in closed loop -- in dark trials
% there's no visual heading reference, so pooling CL+dark together could
% dilute a real CL-only relationship. check that explicitly.
fprintf('\nby light condition:\n');
is_dark_track = is_dark(inc_trials_track);
for c = 1:2
    mask  = (is_dark_track==(c==2));
    wc_c  = cat(2,chunk_window_corr{mask});
    bv_c  = cat(1,chunk_bumpvel{mask});
    rs_c  = cat(1,chunk_rspeed{mask});
    off_c = cat(1,chunk_offset{mask});
    if numel(wc_c) < 5 || numel(bv_c) < 100
        fprintf('  %-12s (too few discernable-bump samples/windows)\n', cond_label{c});
        continue
    end
    vel_r_c = corr(bv_c,rs_c);
    var_c   = circ_var(off_c);
    fprintf('  %-12s n_windows=%4d   position corr (window mean)=%.2f   speed corr=%.2f   offset circ.var=%.2f\n', ...
        cond_label{c}, numel(wc_c), mean(wc_c), vel_r_c, var_c);
end

%% 4) summary figure, using the trial with the most discernable-bump samples as an example
[~,example_trial] = max(n_bump_ok);
example_ii = find(inc_trials_track==example_trial,1);

i    = example_trial;
xf   = all_data(i).ft.xf;
n_im = size(all_data(i).im.d,2);
if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
    xb = all_data(i).ft.xb(:);
else
    xb = linspace(xf(1),xf(end),n_im)';
end

mu_im_ex  = unwrap(all_data(i).im.mu(:));
rho_im_ex = all_data(i).im.rho(:);
keep_ex   = ~is_flash{i}(:);
mu_t_full  = interp1(xb(keep_ex),mu_im_ex(keep_ex),xf);
rho_t_full = interp1(xb(keep_ex),rho_im_ex(keep_ex),xf);
r_speed_full = all_data(i).ft.r_speed;
cue_full     = -all_data(i).ft.cue; % same -cue convention as the pooling loop above

r_ind_ex = find(strcmp(indicator_order,indicator{i}),1);
lag_ex   = lag_frames_grid(indicator_best_idx(r_ind_ex));
[behav_ex,fluor_ex] = apply_lag_multi({r_speed_full,cue_full,xf},{mu_t_full,rho_t_full},lag_ex);
r_speed_l_ex = behav_ex{1}; cue_l_ex = behav_ex{2}; xf_l_ex = behav_ex{3};
mu_l_ex      = fluor_ex{1}; rho_l_ex = fluor_ex{2};

dt_ex = mean(diff(xf_l_ex));
bump_vel_l_ex = gradient(mu_l_ex)/dt_ex;
bump_ok_full  = abs(r_speed_l_ex) > rot_thresh_track & rho_l_ex > rho_thresh_track & abs(bump_vel_l_ex) < bump_vel_thresh & ...
                ~isnan(mu_l_ex) & ~isnan(cue_l_ex) & ~isnan(rho_l_ex);

offset_full = circ_dist(mu_l_ex,cue_l_ex);
offset_full(~bump_ok_full) = nan;

% this trial's own offset/accuracy, NOT the pooled-across-all-trials
% version -- used below for the histogram in this same summary figure so
% it matches what's shown in the other panels
offset_example = offset_full(bump_ok_full);
offset_var_example = circ_var(offset_example);
offset_std_example = circ_std(offset_example);

% this trial's own mean window correlation for the position panel below,
% same metric now used everywhere else in this section (part 1)
pos_rho_example = mean(chunk_window_corr{example_ii});

theta_ex  = all_data(i).im.alpha(1:16); theta_ex = theta_ex(:);
z_hemi_ex = (all_data(i).im.z(1:16,:) + all_data(i).im.z(17:32,:))/2;

mu_plot  = nan_at_wrap(wrap_to_pi(mu_l_ex));
cue_plot = nan_at_wrap(wrap_to_pi(cue_l_ex));

figure(18); clf
set(gcf,'Name',sprintf('bump tracks fly heading: %s',meta_display(all_data(i).meta)),'Position',[100,100,1200,900])

subplot(3,2,[1,2]); hold on
imagesc(xb,theta_ex,z_hemi_ex)
set(gca,'YDir','normal')
h_mu  = plot(xf_l_ex,mu_plot,'-w','LineWidth',1.2);
h_cue = plot(xf_l_ex,cue_plot,'-r','LineWidth',1.2);
xlim([xf(1),xf(end)]); ylim([min(theta_ex),max(theta_ex)])
ylabel('angle (rad)'); xlabel('time (s)')
legend([h_mu,h_cue],{'bump position (mu)','fly heading (-cue)'},'TextColor','w','Location','eastoutside')
title(sprintf('%s (n=%d discernable-bump frames in this trial, lag=%.2fs)',meta_display(all_data(i).meta),n_bump_ok(i),lag_ex*dt_ex),'Interpreter','none')
pos = get(gca,'Position'); colorbar; set(gca,'Position',pos)

subplot(3,2,3)
scatter(wrap_to_pi(chunk_cue{example_ii}),wrap_to_pi(chunk_mu{example_ii}),8,'filled','MarkerFaceAlpha',.3)
xlabel('fly heading (-cue, rad)'); ylabel('bump position (mu, rad)')
axis equal; xlim([-pi,pi]); ylim([-pi,pi])
title(sprintf('position: mean %ds-window circ. corr = %.2f',window_s,pos_rho_example))

subplot(3,2,4)
scatter(chunk_rspeed{example_ii},chunk_bumpvel{example_ii},8,'filled','MarkerFaceAlpha',.3)
xlabel('fly rotational speed (rad/s)'); ylabel('bump speed (rad/s)')
title(sprintf('speed: r = %.2f',vel_r))

subplot(3,2,5)
plot(xf_l_ex,offset_full,'-k')
ylabel('mu - (-cue) offset (rad)'); xlabel('time (s)')
ylim([-pi,pi])
title('offset over time (this trial, discernable-bump periods only)')

subplot(3,2,6)
histogram(offset_example,-pi:pi/12:pi)
xlabel('mu - (-cue) offset (rad)'); ylabel('count (this trial only)')
title(sprintf('this trial: circ. var = %.2f (circ. SD = %.0f%s)',offset_var_example,rad2deg(offset_std_example),char(176)))

%% companion figures for steps 1-3: broken down by indicator and light condition
% same 4x2 grid convention as figures 2/15/16, reusing chunk_mu/chunk_cue/
% chunk_bumpvel/chunk_rspeed/chunk_offset (already restricted to
% discernable-bump samples) from the pooling loop above. scatters are
% subsampled with a fixed stride (not a random draw) purely so rendering
% stays fast for the largest groups (up to ~1M points) -- deterministic so
% the figure is identical every run.
max_scatter_pts = 4000;

%% 1) companion figure: circular correlation between bump position and fly heading, one point per FLY
% same grouped-scatter style as figure 3 (groupplot: one jittered point
% per item + mean+/-sem, across the 8 indicator x condition categories),
% but the "item" here is a FLY, not a trial -- a fly's own discernable-
% bump samples (within one light condition) are pooled together first,
% same per-fly pooling convention used for the R^2 analysis in figures 4-5,
% then ONE circ_corrcc is computed per fly per condition.
cat_labels_pos = {};
for k = 1:numel(indicator_order)
    cat_labels_pos{end+1} = sprintf('%s (CL)',indicator_order{k});   %#ok<SAGROW>
    cat_labels_pos{end+1} = sprintf('%s (dark)',indicator_order{k}); %#ok<SAGROW>
end
n_cat_pos = numel(cat_labels_pos);
base_colors_pos = lines(numel(indicator_order));
cat_colors_pos  = zeros(n_cat_pos,3);
for k = 1:numel(indicator_order) % each indicator's CL/dark pair shares one color
    cat_colors_pos(2*k-1,:) = base_colors_pos(k,:);
    cat_colors_pos(2*k,:)   = base_colors_pos(k,:);
end

fly_pos_corr = [];
fly_cat_x    = [];
for rI = 1:numel(indicator_order)
    for c = 1:2
        rows = find(strcmp(indicator,indicator_order{rI}) & (is_dark==(c==2)) & analysis_ok);
        rows = intersect(rows(:),inc_trials_track(:)); % only trials that made it into the tracking pool (had discernable-bump samples)
        these_flies = unique(fly_num_full(rows));
        for ff = 1:numel(these_flies)
            trial_list = rows(fly_num_full(rows)==these_flies(ff));
            iis = arrayfun(@(t) find(inc_trials_track==t,1), trial_list);
            wc_f = cat(2,chunk_window_corr{iis});
            if isempty(wc_f)
                continue
            end
            fly_pos_corr(end+1) = mean(wc_f); %#ok<AGROW>
            fly_cat_x(end+1)    = 2*(rI-1) + c; %#ok<AGROW>
        end
    end
end

figure(40); clf
set(gcf,'Name','1) bump position (mu vs. -cue): circular correlation per fly','Position',[100,100,900,600])
groupplot(fly_cat_x,fly_pos_corr,cat_labels_pos,cat_colors_pos)
ylabel('circular correlation (mu vs. -cue)')
title('one point per fly, pooled across that fly''s own discernable-bump samples')

%% 2) companion figure: correlation between bump speed and fly rotational speed, one point per FLY
% same per-fly pooling convention as the position companion figure above --
% a fly's own chunk_bumpvel/chunk_rspeed (already lag-shifted per that
% trial's indicator via apply_lag_multi in the pooling loop, and already
% restricted to discernable-bump/rotating/non-flash samples) are pooled
% across that fly's trials within one light condition, then ONE Pearson r
% is computed per fly per condition.
fly_speed_corr  = [];
fly_cat_x_speed = [];
for rI = 1:numel(indicator_order)
    for c = 1:2
        rows = find(strcmp(indicator,indicator_order{rI}) & (is_dark==(c==2)) & analysis_ok);
        rows = intersect(rows(:),inc_trials_track(:));
        these_flies = unique(fly_num_full(rows));
        for ff = 1:numel(these_flies)
            trial_list = rows(fly_num_full(rows)==these_flies(ff));
            iis = arrayfun(@(t) find(inc_trials_track==t,1), trial_list);
            bv_f = cat(1,chunk_bumpvel{iis});
            rs_f = cat(1,chunk_rspeed{iis});
            if numel(bv_f) < 100
                continue
            end
            fly_speed_corr(end+1)  = corr(bv_f,rs_f); %#ok<AGROW>
            fly_cat_x_speed(end+1) = 2*(rI-1) + c;     %#ok<AGROW>
        end
    end
end

figure(43); clf
set(gcf,'Name','2) bump speed vs. fly rotational speed: Pearson r per fly','Position',[100,100,900,600])
groupplot(fly_cat_x_speed,fly_speed_corr,cat_labels_pos,cat_colors_pos)
ylabel('correlation (bump speed vs. fly rotational speed)')
title('one point per fly, pooled across that fly''s own discernable-bump samples')

figure(41); clf
set(gcf,'Name','2) bump speed vs. fly rotational speed, by indicator and light condition','Position',[100,100,1100,1000])
axs41 = gobjects(numel(indicator_order),2);
for rI = 1:numel(indicator_order)
    for c = 1:2
        mask = strcmp(indicator(inc_trials_track),indicator_order{rI}) & (is_dark_track==(c==2));
        bv_g = cat(1,chunk_bumpvel{mask});
        rs_g = cat(1,chunk_rspeed{mask});
        axs41(rI,c) = subplot(numel(indicator_order),2,2*(rI-1)+c);
        if numel(bv_g) < 100
            title(sprintf('%s, %s (n=%d)',indicator_order{rI},cond_label{c},numel(bv_g)))
            continue
        end
        step = max(1,floor(numel(bv_g)/max_scatter_pts));
        show = 1:step:numel(bv_g);
        scatter(rs_g(show),bv_g(show),6,'filled','MarkerFaceAlpha',.3)
        r_g = corr(bv_g,rs_g);
        title(sprintf('%s, %s (n=%d, r=%.2f)',indicator_order{rI},cond_label{c},numel(bv_g),r_g))
        if rI == numel(indicator_order)
            xlabel('fly rotational speed (rad/s)')
        end
        if c == 1
            ylabel('bump speed (rad/s)')
        end
    end
end
linkaxes(axs41(:),'xy')

%% 3) companion figure: circular variance of (mu - (-cue)) offset, one point per FLY
% same windowed approach as the position companion figure (figure 40): each
% trial's chunk_window_var (per-window circ_var(circ_dist(mu,cue)), computed
% over the same sliding windows/discernable-bump restriction as the position
% correlation) is pooled across a fly's own trials within one light
% condition, then averaged into ONE value per fly per condition. alongside
% it, chunk_window_headingvar (per-window circ_var(cue), same windows) is
% averaged the same way as a gray reference -- how much the heading itself
% varied, to judge whether a low offset variance is actually meaningful.
fly_offset_var  = [];
fly_heading_var = [];
fly_cat_x_var   = [];
for rI = 1:numel(indicator_order)
    for c = 1:2
        rows = find(strcmp(indicator,indicator_order{rI}) & (is_dark==(c==2)) & analysis_ok);
        rows = intersect(rows(:),inc_trials_track(:));
        these_flies = unique(fly_num_full(rows));
        for ff = 1:numel(these_flies)
            trial_list = rows(fly_num_full(rows)==these_flies(ff));
            iis = arrayfun(@(t) find(inc_trials_track==t,1), trial_list);
            wv_f = cat(2,chunk_window_var{iis});
            wh_f = cat(2,chunk_window_headingvar{iis});
            if isempty(wv_f)
                continue
            end
            fly_offset_var(end+1)  = mean(wv_f); %#ok<AGROW>
            fly_heading_var(end+1) = mean(wh_f); %#ok<AGROW>
            fly_cat_x_var(end+1)   = 2*(rI-1) + c; %#ok<AGROW>
        end
    end
end

figure(42); clf
set(gcf,'Name','3) bump accuracy: mean windowed circ. var of (mu - (-cue)) offset, per fly','Position',[100,100,900,600])
hold on
n_cat_var = numel(cat_labels_pos);
gray = [.6,.6,.6];
for cIdx = 1:n_cat_var
    y_off = fly_offset_var(fly_cat_x_var==cIdx);
    y_hdg = fly_heading_var(fly_cat_x_var==cIdx);
    if ~isempty(y_off)
        jit = (rand(size(y_off))-.5)*.25;
        scatter(cIdx-0.18+jit,y_off,20,cat_colors_pos(cIdx,:),'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx-0.18,mean(y_off),std(y_off)/sqrt(numel(y_off)),'o','Color',cat_colors_pos(cIdx,:)*.6, ...
            'MarkerFaceColor',cat_colors_pos(cIdx,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    if ~isempty(y_hdg)
        jit = (rand(size(y_hdg))-.5)*.25;
        scatter(cIdx+0.18+jit,y_hdg,20,gray,'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx+0.18,mean(y_hdg),std(y_hdg)/sqrt(numel(y_hdg)),'o','Color',gray*.6, ...
            'MarkerFaceColor',gray*.6,'LineWidth',2,'MarkerSize',7)
    end
end
xticks(1:n_cat_var); xticklabels(cat_labels_pos)
xlim([0.5,n_cat_var+0.5])
y_lims = ylim;
for cIdx = 1:n_cat_var
    text(cIdx,y_lims(1),sprintf('n=%d',sum(fly_cat_x_var==cIdx)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
end
h_off = scatter(nan,nan,20,[0,0,0],'filled');
h_hdg = scatter(nan,nan,20,gray,'filled');
legend([h_off,h_hdg],{'circ.var(heading - bump)','circ.var(heading)'},'Location','eastoutside')
ylabel('circular variance (window mean)')
title('one point per fly: offset variance (color) vs. heading variance alone (gray)')

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

function is_flash = detect_flash_frames(f, mad_k)
    % f: [n_wedges x n_frames] raw fluorescence (before dF/F normalization).
    % flags frames that are both a brightness outlier relative to the rest
    % of the trial, and unusually uniform across wedges -- consistent with
    % a bright light passing under the objective, not genuine spatially-
    % localized PB bump activity (which would raise the mean but *increase*
    % the across-wedge spread, not shrink it).
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;

    med_m  = median(frame_mean);
    mad_m  = mad(frame_mean,1); % median absolute deviation
    med_cv = median(frame_cv);

    is_bright  = frame_mean > med_m + mad_k*mad_m;
    is_uniform = frame_cv <= med_cv;
    is_flash   = is_bright & is_uniform;
end

function [mu,kappa,baseline,amp,valid,actual_peak,actual_trough] = fit_hemisphere_vm(z_im, alpha, flash_mask, weak_signal_frac, min_active_wedges)
    % averages the two 16-wedge PB hemispheres (z_im rows 1:16 and 17:32
    % share identical angular positions) and fits each frame's combined
    % profile to a von Mises bump y = baseline + amp*exp(kappa*(cos(theta-mu)-1)).
    % z_im is expected to be z-scored activity (im.z), not dFF -- z-scoring
    % puts every wedge on the same scale regardless of its own baseline
    % brightness, which matters here since both the fit and the
    % rectification step compare magnitudes across wedges.
    %
    % mu/kappa come from the resultant vector of the RECTIFIED activity
    % (negative values clipped to 0) treated as circular weights, same
    % approach as circ_stats/circ_r.m + circ_kappa.m in this repo: R is
    % corrected for the 16 evenly-spaced angle bins (Zar's binned-data
    % correction), then kappa is Fisher's piecewise approximation from R
    % (the same formula circ_kappa.m uses). baseline/amp then follow from
    % an exact 2-parameter linear least-squares fit given that shape.
    %
    % a frame is marked invalid (NaN) for mu/kappa/baseline/amp if it's a
    % detected flash, if its total rectified activity is below
    % weak_signal_frac of the trial's own median (near-zero rectified
    % weight makes the resultant length numerically unstable -- it can
    % spuriously approach 1 from pure noise), OR if fewer than
    % min_active_wedges of the 16 wedges have any positive rectified
    % activity at all: a single noisy wedge with everything else
    % negative/zero can pass the sum_w check on total magnitude alone
    % while still being a degenerate, non-bump-shaped "spike" (verified
    % directly against real data: exactly this pattern produced a
    % kappa~100 degenerate fit before this guard was added). a genuine
    % bump should be spread across multiple adjacent wedges.
    % actual_peak/actual_trough (the data's own max/min across the 16
    % wedges, not model-derived) only need a flash excluded -- they're
    % always well-defined otherwise.
    n = 16;
    theta  = alpha(1:n); theta = theta(:);
    d_hemi = (z_im(1:n,:) + z_im(n+1:2*n,:)) / 2;

    actual_peak   = max(d_hemi,[],1);
    actual_trough = min(d_hemi,[],1);
    actual_peak(flash_mask)   = nan;
    actual_trough(flash_mask) = nan;

    w     = max(d_hemi,0);
    sum_w = sum(w,1);
    weak  = sum_w < weak_signal_frac*median(sum_w) | sum(w>0,1) < min_active_wedges;

    C  = sum(w.*cos(theta),1);
    S  = sum(w.*sin(theta),1);
    mu = atan2(S,C);

    d_bin    = 2*pi/n;
    bin_corr = d_bin/2/sin(d_bin/2); % Zar's correction for binned circular data
    R = min(sqrt(C.^2+S.^2)./sum_w * bin_corr, 0.995);

    kappa = nan(size(R));
    lo  = R < 0.53; mid = R>=0.53 & R<0.85; hi = R>=0.85;
    kappa(lo)  = 2*R(lo) + R(lo).^3 + 5*R(lo).^5/6;
    kappa(mid) = -.4 + 1.39*R(mid) + 0.43./(1-R(mid));
    kappa(hi)  = 1./(R(hi).^3 - 4*R(hi).^2 + 3*R(hi));

    basis = exp(kappa.*(cos(theta-mu)-1));
    Sx  = sum(basis,1);   Sxx = sum(basis.^2,1);
    Sy  = sum(d_hemi,1);  Sxy = sum(basis.*d_hemi,1);
    det_ = n*Sxx - Sx.^2;
    baseline = (Sxx.*Sy - Sx.*Sxy) ./ det_;
    amp      = (n*Sxy - Sx.*Sy) ./ det_;

    valid = ~weak & ~flash_mask;
    mu(~valid) = nan; kappa(~valid) = nan; baseline(~valid) = nan; amp(~valid) = nan;
end

function frac = half_max_width_frac(kappa)
    % fraction of the full PB (0-1) whose von Mises fit value is at or
    % above the half-max point (baseline + amp/2). below kappa=log(2)/2
    % the model never drops to half-max anywhere around the circle -- the
    % bump is broad enough that the whole PB counts as "within half-max".
    frac = ones(size(kappa));
    valid = kappa > log(2)/2;
    delta = acos(1 - log(2)./kappa(valid)); % half-width from the peak, radians
    frac(valid) = delta/pi; % (2*delta)/(2*pi)
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

function [leading_out, lagging_out] = apply_lag_multi(leading_in, lagging_in, lag)
    % cell-array generalization of apply_lag3: shifts every series in
    % lagging_in (fluorescence-derived) behind every series in leading_in
    % (behavior), same "fluorescence lags behavior by lag frames" convention.
    n = numel(leading_in{1});
    if lag == 0
        idx_lead = 1:n; idx_lag = 1:n;
    elseif lag > 0
        idx_lead = 1:n-lag; idx_lag = lag+1:n;
    else
        idx_lead = -lag+1:n; idx_lag = 1:n+lag;
    end
    leading_out = cellfun(@(x) x(idx_lead), leading_in, 'UniformOutput', false);
    lagging_out = cellfun(@(x) x(idx_lag), lagging_in, 'UniformOutput', false);
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

function y = wrap_to_pi(x)
    % wraps angles (radians) into (-pi,pi], without requiring the Mapping
    % Toolbox's wrapToPi
    y = mod(x+pi,2*pi) - pi;
end

function y = nan_at_wrap(x)
    % inserts NaN right after a >pi jump in a wrapped angle trace, so
    % plotting it doesn't draw a spurious line straight across the
    % wrap-around discontinuity
    y = x;
    jump = [false; abs(diff(x)) > pi];
    y(jump) = nan;
end
