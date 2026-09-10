%% grab_freeze_scratch
% Rotational speed vs. PB fluorescence (summed & peak), across sensors
% (GCaMP, GRAB(DA2m), dLight) and light/cue conditions (closed loop,
% frozen cue, dark).
%
% Step 1 (this installment): load the same combined dataset as
% lpsp_compartments_claude_script.m (lpsp_cl + lpsp_cl_redo + epg_dlight),
% using the same dataset-tagging/condition/indicator/fly-ID conventions,
% then detect freeze blocks in every closed-loop trial using the same
% detector as epg_dlight_claude_script.m ("frozen cue" = fly rotating but
% the visual cue isn't moving). Reports how many flies are in the combined
% dataset and how many freeze blocks were found per fly, as a sanity check
% before building the fluorescence-vs-speed comparison on top of it.
%
% syt7f and syt8m (synaptotagmin-tagged GCaMP calcium sensors) are kept as
% separate sensor groups here, alongside GRAB(DA2m) and dLight -- i.e. four
% sensor groups, not three, since the two GCaMP variants aren't pooled.

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
% trials already have ft.pattern stored)
csv_indexes = containers.Map('KeyType','char','ValueType','any');
for k = 1:numel(dataset_names)
    fprintf('indexing trialSettings.csv under %s ...\n', base_dirs{k});
    csv_indexes(dataset_names{k}) = build_csv_index(base_dirs{k});
end

%% determine closed-loop vs. dark per trial, indicator, and check field completeness
ft_required = {'xf','f_speed','r_speed','cue'};
im_required = {'mu','rho','z','d','f','alpha'};

is_dark        = false(1,n_trials);
condition_src  = cell(1,n_trials);
is_ambiguous   = false(1,n_trials);
odd_panels     = false(1,n_trials);
is_empty_trial = false(1,n_trials);
missing_ft     = cell(1,n_trials);
missing_im     = cell(1,n_trials);

for i = 1:n_trials
    meta = all_data(i).meta;

    is_empty_trial(i) = isempty(meta) || ~(ischar(meta) || (isstring(meta) && isscalar(meta)));
    if is_empty_trial(i)
        condition_src{i} = 'empty trial (no meta/ft/im)';
        missing_ft{i}    = ft_required;
        missing_im{i}    = im_required;
        continue
    end

    if isfield(all_data(i).ft,'pattern') && ~isempty(all_data(i).ft.pattern)
        pattern_str      = char(all_data(i).ft.pattern);
        is_dark(i)       = contains(pattern_str,'background','IgnoreCase',true);
        condition_src{i} = 'ft.pattern (already stored)';
    else
        tname = trial_folder_name(meta);
        res = resolve_csv_for_trial(meta, tname, csv_indexes(all_data(i).dataset));

        odd_panels(i) = ~strcmp(res.status,'not_found') && ~isempty(res.info.panels_mode) && ...
                        ~contains(res.info.panels_mode,'closed loop','IgnoreCase',true);

        switch res.status
            case 'ambiguous'
                is_ambiguous(i)  = true;
                is_dark(i)       = contains(meta,'_dark','IgnoreCase',true);
                condition_src{i} = 'AMBIGUOUS trialSettings.csv candidates -- used folder name, needs manual review';
            case 'not_found'
                is_dark(i)       = contains(meta,'_dark','IgnoreCase',true);
                condition_src{i} = 'folder name fallback (no trialSettings.csv found)';
            otherwise
                is_dark(i)       = res.info.is_dark;
                condition_src{i} = sprintf('trialSettings.csv (%s match)', res.status);
        end
    end

    missing_ft{i} = find_missing_fields(all_data(i).ft, ft_required);
    missing_im{i} = find_missing_fields(all_data(i).im, im_required);
end

complete  = cellfun(@isempty,missing_ft) & cellfun(@isempty,missing_im);
trace_ok  = ~is_empty_trial & complete;

%% indicator (sensor) label per trial, from the folder name
indicator_order = {'syt7f','syt8m','GRAB(DA2m)','dLight'};
indicator = cell(1,n_trials);
for i = 1:n_trials
    if is_empty_trial(i)
        continue
    end
    indicator{i} = trial_indicator(all_data(i).meta);
end

% sensor group labels used for the fluorescence-vs-speed comparison below --
% syt7f/syt8m kept separate (see header comment), GRAB(DA2m) relabeled for
% brevity, dLight unchanged
sensor_group = cell(1,n_trials);
sensor_group(strcmp(indicator,'syt7f'))      = {'syt7f'};
sensor_group(strcmp(indicator,'syt8m'))      = {'syt8m'};
sensor_group(strcmp(indicator,'GRAB(DA2m)')) = {'GRAB(DA)'};
sensor_group(strcmp(indicator,'dLight'))     = {'dLight'};

sensor_group_order = {'syt7f','syt8m','GRAB(DA)','dLight'};

fprintf('\n=== indicator counts (raw label / pooled group) ===\n');
for k = 1:numel(indicator_order)
    idx = strcmp(indicator,indicator_order{k});
    fprintf('  %-12s n=%3d   closed loop=%3d   dark=%3d\n', indicator_order{k}, sum(idx), sum(idx & ~is_dark), sum(idx & is_dark));
end
fprintf('  ---\n');
for k = 1:numel(sensor_group_order)
    idx = strcmp(sensor_group,sensor_group_order{k});
    fprintf('  %-12s n=%3d   closed loop=%3d   dark=%3d\n', sensor_group_order{k}, sum(idx), sum(idx & ~is_dark), sum(idx & is_dark));
end

%% fly ID per trial
fly_id = cell(1,n_trials);
for i = find(~is_empty_trial)
    fly_id{i} = trial_fly_id(all_data(i).dataset,all_data(i).meta);
end

[fly_list,~,fly_num] = unique(fly_id(~is_empty_trial));
fly_num_full = nan(1,n_trials);
fly_num_full(~is_empty_trial) = fly_num;
n_flies = numel(fly_list);

fprintf('\nn flies in combined dataset (lpsp_cl + lpsp_cl_redo + epg_dlight): %d\n', n_flies);
fprintf('n trials: %d (%d complete/usable, %d empty placeholder, %d ambiguous condition, %d odd panelsMode)\n', ...
    n_trials, sum(trace_ok), sum(is_empty_trial), sum(is_ambiguous), sum(odd_panels));

%% detect freeze blocks in each closed-loop trial
% a "freeze block" is a period where the fly is actively rotating (high
% |r_speed|) but the visual cue fails to move (low |d(cue)/dt|) -- same
% detector as epg_dlight_claude_script.m, applied here across all three
% combined datasets' closed-loop trials (dark trials have no cue to freeze).

r_thresh    = .1;      % rad/s, fly must be rotating at least this fast
cue_thresh  = 1e-2;    % rad/s, cue must be moving less than this to count as frozen
pre_smooth  = 90;      % smoothing window (frames) applied to the raw speeds before thresholding
post_smooth = 90;      % smoothing window (frames) used to accumulate recent speed into a slower "is moving" signal
min_block_s = .5;      % minimum duration (s) for a freeze period to count as a real block, not noise
max_gap_s   = 2;       % merge freeze periods separated by a gap shorter than this

freeze_idx_all = cell(1,n_trials); % per-sample: confidently-detectable frozen samples
block_mask_all = cell(1,n_trials); % per-sample, gap-merged version, for block boundaries
frac_freeze    = nan(1,n_trials);
n_blocks       = nan(1,n_trials);
block_dur      = cell(1,n_trials);
cr_all         = cell(1,n_trials); % smoothed |rotational speed|, kept for reuse in the fluorescence analysis
dt_all         = nan(1,n_trials);  % per-trial ft.xf sample interval, kept for reuse below

for i = find(trace_ok)
    dt = mean(diff(all_data(i).ft.xf));
    dt_all(i) = dt;

    r  = all_data(i).ft.r_speed;
    dr = smoothdata(r,'gaussian',pre_smooth);
    cr = smoothdata(abs(dr),'movmean',[post_smooth,0]);
    cr_all{i} = cr;

    if is_dark(i); continue; end

    min_block_frames = max(1,round(min_block_s/dt));
    max_gap_frames    = max(1,round(max_gap_s/dt));
    c  = all_data(i).ft.cue;

    dc = smoothdata([diff(unwrap(c))/dt;0],'gaussian',pre_smooth);
    cc = smoothdata(abs(dc),'movmean',[post_smooth,0]);

    freeze_idx = cr > r_thresh & bwareaopen(cc < cue_thresh, min_block_frames);
    freeze_idx_all{i} = freeze_idx;
    frac_freeze(i) = mean(freeze_idx);

    block_mask = ~bwareaopen(~freeze_idx, max_gap_frames);
    block_mask_all{i} = block_mask;

    labeled = bwlabel(block_mask);
    n_blocks(i) = max(labeled);
    dur = nan(n_blocks(i),1);
    for b = 1:n_blocks(i)
        dur(b) = sum(labeled==b) * dt;
    end
    block_dur{i} = dur;
end

has_freeze = trace_ok & ~is_dark & frac_freeze > .01 & n_blocks >= 1;
fprintf('\nclosed loop trials with detectable freeze blocks: %d / %d\n', sum(has_freeze), sum(trace_ok & ~is_dark));

%% per-fly freeze block summary
fly_n_blocks   = zeros(n_flies,1);
fly_cl_time_s  = zeros(n_flies,1); % total closed-loop time, this fly
fly_frz_time_s = zeros(n_flies,1); % total freeze time, this fly
fly_dataset    = cell(n_flies,1);
fly_sensors    = cell(n_flies,1);

for f = 1:n_flies
    trials_f = find(fly_num_full==f & trace_ok & ~is_dark);
    fly_n_blocks(f) = sum(n_blocks(trials_f),'omitnan');
    for i = trials_f
        dt = mean(diff(all_data(i).ft.xf));
        fly_cl_time_s(f)  = fly_cl_time_s(f)  + numel(all_data(i).ft.xf)*dt;
        fly_frz_time_s(f) = fly_frz_time_s(f) + sum(freeze_idx_all{i})*dt;
    end
    if ~isempty(trials_f)
        fly_dataset{f} = all_data(trials_f(1)).dataset;
    end
    fly_sensors{f} = strjoin(unique(sensor_group(fly_num_full==f & ~is_empty_trial)),'+');
end

fprintf('\n=== freeze blocks per fly (closed-loop trials only) ===\n');
for f = 1:n_flies
    fprintf('  [%-12s %-14s] %-45s   n_blocks=%3d   CL time=%6.0fs   freeze time=%6.0fs (%.1f%%)\n', ...
        fly_dataset{f}, fly_sensors{f}, fly_list{f}, fly_n_blocks(f), fly_cl_time_s(f), fly_frz_time_s(f), ...
        100*fly_frz_time_s(f)/max(fly_cl_time_s(f),eps));
end

fprintf('\ntotal freeze blocks across all flies: %d\n', sum(fly_n_blocks));
fprintf('flies with zero detected freeze blocks: %d/%d\n', sum(fly_n_blocks==0 & fly_cl_time_s>0), sum(fly_cl_time_s>0));

%% summary figure: flies and freeze blocks
figure(1); clf; set(gcf,'Name','grab_freeze_scratch: fly & freeze block counts','Position',[100,100,900,400])

subplot(1,2,1)
bar(categorical(sensor_group_order,sensor_group_order), ...
    arrayfun(@(k) numel(unique(fly_num_full(strcmp(sensor_group,sensor_group_order{k})))), 1:numel(sensor_group_order)))
ylabel('# flies (>=1 trial)')
title(sprintf('%d flies total',n_flies))

subplot(1,2,2)
histogram(fly_n_blocks,'BinMethod','integers')
xlabel('# freeze blocks per fly (closed-loop trials)')
ylabel('# flies')
title(sprintf('%d total blocks, %d flies with 0',sum(fly_n_blocks),sum(fly_n_blocks==0 & fly_cl_time_s>0)))

%% flash-frame detection, then summed & peak PB fluorescence per trial
% same detector as lpsp_compartments_claude_script.m: a frame is a "flash"
% (bright light passing under the objective, not a real bump) if it's both
% a brightness outlier and unusually spatially uniform across the 32 PB
% wedges. flagged frames are excluded before interpolating dF/F onto the
% behavior (ft.xf) timebase. runs on im.z (z-scored per wedge), matching
% the dff_avg/peak computation style in lpsp_compartments_claude_script.m.
flash_mad_thresh = 8;

is_flash   = cell(1,n_trials);
sum_cell   = cell(1,n_trials); % summed z-scored dFF, all 32 wedges
peak_cell  = cell(1,n_trials); % peak z-scored dFF, all 32 wedges (true max, no hemisphere averaging)
peak3_cell = cell(1,n_trials); % alternative peak: hemisphere-averaged (16 wedges), best 3 circularly-neighboring wedges averaged -- see section below

for i = find(trace_ok)
    is_flash{i} = detect_flash_frames(all_data(i).im.f, flash_mad_thresh);

    xf   = all_data(i).ft.xf;
    n_im = size(all_data(i).im.z,2);
    if isfield(all_data(i).ft,'xb') && numel(all_data(i).ft.xb) == n_im
        xb = all_data(i).ft.xb(:);
    else
        xb = linspace(xf(1),xf(end),n_im)';
    end

    keep = ~is_flash{i}(:);
    if sum(keep) < 2
        keep = true(size(keep)); % too few non-flash frames -- keep all rather than fail
    end

    dff = interp1(xb(keep), all_data(i).im.z(:,keep)', xf);
    sum_cell{i}  = sum(dff,2);
    peak_cell{i} = max(dff,[],2);

    peak3_im       = best3_neighbor_peak(all_data(i).im.z);
    peak3_cell{i}  = interp1(xb(keep), peak3_im(keep)', xf);
end
fprintf('\ncomputed summed/peak fluorescence for %d/%d trials (flash frames excluded)\n', sum(trace_ok), n_trials);

%% rotational speed vs. summed/peak PB fluorescence, per sensor x condition
% for each sensor group and each fly, pool that fly's own closed-loop
% samples (split into "closed loop" = cue moving vs. "frozen" = detected
% freeze block, via freeze_idx_all above) and dark samples, then bin by
% the smoothed |rotational speed| (cr_all). same per-fly-then-pooled
% approach and bin_by_speed helper as epg_dlight_claude_script.m.
speed_edges = 0:.25:3;
speed_x     = speed_edges(1:end-1) + diff(speed_edges)/2;
min_bin_s   = 5; % require this many seconds of pooled data per speed bin to trust it
nominal_dt  = mean(dt_all,'omitnan');

conditions = {'closed loop','frozen','dark'};
cond_rgb   = {[0,0,0],[1,0,0],[0,0,1]};

min_flies_per_bin = 3; % a speed bin's mean+/-sem line is only drawn where MORE than this many flies have data

binned_sum   = nan(numel(speed_x),3,n_flies,numel(sensor_group_order));
binned_peak  = nan(numel(speed_x),3,n_flies,numel(sensor_group_order));
binned_peak3 = nan(numel(speed_x),3,n_flies,numel(sensor_group_order));

for sg_idx = 1:numel(sensor_group_order)
    sg = sensor_group_order{sg_idx};
    for f = 1:n_flies
        cl_trials   = find(fly_num_full==f & trace_ok & ~is_dark & strcmp(sensor_group,sg));
        dark_trials = find(fly_num_full==f & trace_ok &  is_dark & strcmp(sensor_group,sg));

        cr_cl    = cat(1,cr_all{cl_trials});
        is_frz   = cat(1,freeze_idx_all{cl_trials});
        sum_cl   = cat(1,sum_cell{cl_trials});
        peak_cl  = cat(1,peak_cell{cl_trials});
        peak3_cl = cat(1,peak3_cell{cl_trials});

        cr_dk    = cat(1,cr_all{dark_trials});
        sum_dk   = cat(1,sum_cell{dark_trials});
        peak_dk  = cat(1,peak_cell{dark_trials});
        peak3_dk = cat(1,peak3_cell{dark_trials});

        binned_sum(:,1,f,sg_idx)   = bin_by_speed(cr_cl(~is_frz), sum_cl(~is_frz),   speed_edges, nominal_dt, min_bin_s);
        binned_sum(:,2,f,sg_idx)   = bin_by_speed(cr_cl(is_frz),  sum_cl(is_frz),    speed_edges, nominal_dt, min_bin_s);
        binned_sum(:,3,f,sg_idx)   = bin_by_speed(cr_dk,          sum_dk,            speed_edges, nominal_dt, min_bin_s);

        binned_peak(:,1,f,sg_idx)  = bin_by_speed(cr_cl(~is_frz), peak_cl(~is_frz),  speed_edges, nominal_dt, min_bin_s);
        binned_peak(:,2,f,sg_idx)  = bin_by_speed(cr_cl(is_frz),  peak_cl(is_frz),   speed_edges, nominal_dt, min_bin_s);
        binned_peak(:,3,f,sg_idx)  = bin_by_speed(cr_dk,          peak_dk,           speed_edges, nominal_dt, min_bin_s);

        binned_peak3(:,1,f,sg_idx) = bin_by_speed(cr_cl(~is_frz), peak3_cl(~is_frz), speed_edges, nominal_dt, min_bin_s);
        binned_peak3(:,2,f,sg_idx) = bin_by_speed(cr_cl(is_frz),  peak3_cl(is_frz),  speed_edges, nominal_dt, min_bin_s);
        binned_peak3(:,3,f,sg_idx) = bin_by_speed(cr_dk,          peak3_dk,          speed_edges, nominal_dt, min_bin_s);
    end
end

metric_names = {'summed dFF (z-scored, all 32 wedges)', ...
                 'peak dFF (max, z-scored, all 32 wedges)', ...
                 'peak dFF (best 3-neighbor avg, hemisphere-averaged 16 wedges)'};
metric_data  = {binned_sum, binned_peak, binned_peak3};

for m = 1:numel(metric_data)
    figure(2+m); clf; set(gcf,'Name',sprintf('%s vs. rotational speed, by sensor',metric_names{m}),'Position',[100,100,1000,800])
    axs = gobjects(1,numel(sensor_group_order));
    for sg_idx = 1:numel(sensor_group_order)
        axs(sg_idx) = subplot(2,2,sg_idx); hold on
        n_flies_sg = sum(any(any(~isnan(metric_data{m}(:,:,:,sg_idx)),1),2));

        for c = 1:3
            y = squeeze(metric_data{m}(:,c,:,sg_idx)); % speed bins x flies
            plot(speed_x,y,'Color',[cond_rgb{c},.15])
        end
        h = nan(3,1);
        for c = 1:3
            y = squeeze(metric_data{m}(:,c,:,sg_idx));
            n_ok   = sum(~isnan(y),2);
            mean_y = mean(y,2,'omitnan');
            sem_y  = std(y,0,2,'omitnan') ./ sqrt(n_ok);
            mean_y(n_ok <= min_flies_per_bin) = nan; % underpowered bins: don't draw a mean/sem line
            sem_y(n_ok <= min_flies_per_bin)  = nan;

            plot_mean_sem(speed_x,mean_y,sem_y,cond_rgb{c});
            h(c) = plot(speed_x,mean_y,'-o','Color',cond_rgb{c}, ...
                'MarkerFaceColor',cond_rgb{c},'MarkerSize',4,'linewidth',2);
        end
        legend(h,conditions,'Location','best')
        xlabel('rotational speed (rad/s)')
        ylabel(metric_names{m})
        title(sprintf('%s (n=%d flies)',sensor_group_order{sg_idx},n_flies_sg))
    end
    linkaxes(axs,'x')
end

%% save all figures as PDF
save_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'ugly_figures','grab_freeze_scratch');
if ~exist(save_dir,'dir'); mkdir(save_dir); end

for fig_num = 1:5
    fh = findobj('Type','figure','Number',fig_num);
    if isempty(fh); continue; end

    file_name = lower(get(fh,'Name'));
    file_name = regexprep(file_name,'[^a-z0-9]+','_');
    file_name = regexprep(file_name,'^_+|_+$','');

    try
        exportgraphics(fh, fullfile(save_dir,[file_name,'.pdf']))
        fprintf('saved figure %d -> %s.pdf\n', fig_num, file_name);
    catch err
        % most common cause: the PDF is still open (e.g. in a viewer),
        % which locks the file for writing -- don't let that kill the rest
        % of the run
        fprintf('WARNING: could not save figure %d (%s.pdf): %s\n', fig_num, file_name, err.message);
    end
end

%% Functions

function combined = combine_datasets(varargin)
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

function name = trial_folder_name(meta_path)
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        name = parts{end-1};
    else
        name = parts{end};
    end
end

function idx = build_csv_index(base_dir)
    d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
    idx.names = cell(numel(d),1);
    idx.paths = cell(numel(d),1);
    for i = 1:numel(d)
        [~,trial_name] = fileparts(fileparts(d(i).folder));
        idx.names{i} = trial_name;
        idx.paths{i} = fullfile(d(i).folder,d(i).name);
    end

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
        used_fallback = true;
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

function ind = trial_indicator(meta_path)
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
            fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
            if isempty(fly_part)
                fid = meta_path;
            else
                fid = strjoin(parts(1:fly_part(1)),filesep);
            end
        otherwise % lpsp_cl: no "fly N" folder -- date + trailing suffix digit
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

function is_flash = detect_flash_frames(f, mad_k)
    % f: [n_wedges x n_frames] raw fluorescence (before dF/F normalization).
    % flags frames that are both a brightness outlier relative to the rest
    % of the trial, and unusually uniform across wedges -- consistent with
    % a bright light passing under the objective, not genuine spatially-
    % localized PB bump activity.
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;

    med_m  = median(frame_mean);
    mad_m  = mad(frame_mean,1);
    med_cv = median(frame_cv);

    is_bright  = frame_mean > med_m + mad_k*mad_m;
    is_uniform = frame_cv <= med_cv;
    is_flash   = is_bright & is_uniform;
end

function m = bin_by_speed(cr, y, edges, dt, min_sec)
    m = nan(length(edges)-1,1);
    for j = 1:(length(edges)-1)
        idx = cr >= edges(j) & cr < edges(j+1);
        if sum(idx)*dt > min_sec
            m(j) = mean(y(idx),'omitnan');
        end
    end
end

function h = plot_mean_sem(t,m,s,c)
    % draws a mean+/-sem patch from precomputed mean (m) and sem (s) vectors
    % (as opposed to plotsem's old behavior of computing them internally from
    % raw data) -- lets the caller NaN out underpowered bins beforehand so
    % the patch and the mean line it accompanies always agree on which bins
    % are shown.
    t = reshape(t,1,[]); m = reshape(m,1,[]); s = reshape(s,1,[]);

    valid = ~isnan(m) & ~isnan(s);
    t = t(valid); m = m(valid); s = s(valid);

    h = patch([t,fliplr(t)],[m+s,fliplr(m-s)],c,'FaceAlpha',.2,'EdgeColor','none');
end

function peak3 = best3_neighbor_peak(z_im)
    % z_im: [n_wedges x n_frames] z-scored PB activity (im.z), full 32-wedge
    % ring. averages the two hemispheres wedge-for-wedge (wedge k with wedge
    % k+n/2, e.g. 1 with 17) into one 16-wedge ring, then for every frame
    % finds the 3 circularly-adjacent wedges (wrapping 16->1) whose average
    % is highest, and returns that average as the frame's peak estimate --
    % an alternative to the single-wedge max, less sensitive to one noisy
    % wedge spiking on its own.
    n_full = size(z_im,1);
    n_hemi = n_full/2;
    z_hemi = (z_im(1:n_hemi,:) + z_im(n_hemi+1:end,:)) / 2; % [16 x n_frames]

    avg3 = zeros(n_hemi,size(z_hemi,2));
    for w = 0:2
        avg3 = avg3 + circshift(z_hemi,-w,1);
    end
    avg3 = avg3/3; % avg3(k,:) = mean of wedges k,k+1,k+2 (circular)

    peak3 = max(avg3,[],1);
end
