%% load in data
base_dir = 'Z:\pablo\epg_dlight\';
data_file = fullfile(fileparts(fileparts(mfilename('fullpath'))),'.data','epg_dlight_20260415.mat');
load(data_file); % loads all_data

n_trials = length(all_data);

%% classify trials: dark vs closed loop, and assign each trial to a fly
is_dark  = false(n_trials,1);
fly_str  = cell(n_trials,1);

for i = 1:n_trials
    is_dark(i) = contains(all_data(i).ft.pattern,'background');

    tmp = strsplit(all_data(i).meta,'\');
    fly_folder = tmp(cellfun(@(x)(startsWith(x,'fly ')),tmp));
    if isempty(fly_folder)
        fly_str{i} = all_data(i).meta;
    else
        fly_str{i} = strjoin(tmp(1:find(strcmp(tmp,fly_folder{1}))),'\');
    end
end

[fly_list,~,fly_id] = unique(fly_str);
n_flies = length(fly_list);

fprintf('n flies: %d\n', n_flies);
fprintf('n trials: %d (dark: %d, closed loop: %d)\n', n_trials, sum(is_dark), sum(~is_dark));

%% detect freeze blocks in each closed-loop trial
% a "freeze block" is a period where the fly is actively rotating (high
% |r_speed|) but the visual cue fails to move (low |d(cue)/dt|), i.e. the
% closed-loop feedback has been decoupled from behavior. dark trials have
% no meaningful cue, so freeze blocks are only assessed for closed-loop
% (bar pattern) trials.

r_thresh    = .1;      % rad/s, fly must be rotating at least this fast
cue_thresh  = 1e-2;    % rad/s, cue must be moving less than this to count as frozen
pre_smooth  = 90;      % smoothing window (frames) applied to the raw speeds before thresholding
post_smooth = 90;      % smoothing window (frames) used to accumulate recent speed into a slower "is moving" signal
min_block_s = .5;      % minimum duration (s) for a freeze period to count as a real block, not noise

freeze_idx_all  = cell(n_trials,1);
frac_freeze     = nan(n_trials,1);
n_blocks        = nan(n_trials,1);
block_dur       = cell(n_trials,1);
cr_all          = cell(n_trials,1); % smoothed |rotational speed|, kept for reuse below

for i = 1:n_trials
    dt = mean(diff(all_data(i).ft.xf));

    r  = all_data(i).ft.r_speed;
    dr = smoothdata(r,'gaussian',pre_smooth);
    cr = smoothdata(abs(dr),'movmean',[post_smooth,0]);
    cr_all{i} = cr;

    if is_dark(i); continue; end

    min_block_frames = max(1,round(min_block_s/dt));
    c  = all_data(i).ft.cue;

    dc = smoothdata([diff(unwrap(c))/dt;0],'gaussian',pre_smooth);
    cc = smoothdata(abs(dc),'movmean',[post_smooth,0]);

    freeze_idx = cr > r_thresh & bwareaopen(cc < cue_thresh, min_block_frames);
    freeze_idx_all{i} = freeze_idx;
    frac_freeze(i) = mean(freeze_idx);

    labeled = bwlabel(freeze_idx);
    n_blocks(i) = max(labeled);
    dur = nan(n_blocks(i),1);
    for b = 1:n_blocks(i)
        dur(b) = sum(labeled==b) * dt;
    end
    block_dur{i} = dur;
end

has_freeze = ~is_dark & frac_freeze > .01 & n_blocks >= 1; % closed-loop trials with a non-trivial amount of freeze time

fprintf('closed loop trials with detectable freeze blocks: %d / %d\n', sum(has_freeze), sum(~is_dark));

%% summary figure: fly/trial counts
figure('Name','epg_dlight summary'); clf

subplot(2,2,1)
trial_type = repmat({'closed loop'},n_trials,1);
trial_type(is_dark)  = {'dark'};
trial_type(has_freeze) = {'closed loop + freeze'};
[types,~,type_id] = unique(trial_type);
counts = histcounts(type_id,.5:1:(length(types)+.5));
bar(counts)
set(gca,'XTickLabel',types)
ylabel('# trials')
title(sprintf('%d trials, %d flies',n_trials,n_flies))

subplot(2,2,2)
trials_per_fly = histcounts(fly_id,.5:1:(n_flies+.5));
histogram(trials_per_fly,'BinMethod','integers')
xlabel('trials per fly')
ylabel('# flies')

subplot(2,2,3)
histogram(n_blocks(~is_dark),'BinMethod','integers')
xlabel('# freeze blocks per closed-loop trial')
ylabel('# trials')

subplot(2,2,4)
all_dur = cat(1,block_dur{~is_dark});
histogram(all_dur,0:1:60)
xlabel('freeze block duration (s)')
ylabel('# blocks')

%% example trial: fly heading, cue, and detected freeze blocks
i = find(has_freeze,1);

figure('Name','example freeze detection'); clf
a1 = subplot(2,1,1); hold on
plot(all_data(i).ft.xf, all_data(i).ft.heading,'k')
plot(all_data(i).ft.xf, -unwrap(all_data(i).ft.cue),'c')
y = ylim;
patch([all_data(i).ft.xf;flipud(all_data(i).ft.xf)], ...
      [1e6*freeze_idx_all{i};-1e6*flipud(freeze_idx_all{i})], ...
      'r','FaceAlpha',.2,'EdgeColor','none')
ylim(y)
legend('fly heading','cue position','freeze block')
title(all_data(i).meta,'Interpreter','none')
ylabel('position (rad, unwrapped)')

a2 = subplot(2,1,2); hold on
plot(all_data(i).ft.xf, all_data(i).ft.r_speed,'k')
y = ylim;
patch([all_data(i).ft.xf;flipud(all_data(i).ft.xf)], ...
      [1e6*freeze_idx_all{i};-1e6*flipud(freeze_idx_all{i})], ...
      'r','FaceAlpha',.2,'EdgeColor','none')
ylim(y)
xlabel('time (s)')
ylabel('r speed (rad/s)')

linkaxes([a1,a2],'x')

%% fluorescence vs. rotational speed: closed loop vs. frozen vs. dark
% for each trial, interpolate the (z-scored) per-cluster fluorescence onto
% the behavior timescale, then collapse across clusters into three
% summary traces. each timepoint is labeled closed loop (moving cue),
% frozen (freeze block), or dark, and binned by the smoothed rotational
% speed computed above (cr_all).

speed_edges = 0:.25:3;                 % rad/s bins for cr
speed_x     = speed_edges(1:end-1) + diff(speed_edges)/2;
min_bin_s   = 5;                       % require this many seconds of data per bin to trust it
nominal_dt  = mean(diff(all_data(1).ft.xf));

sum_cell   = cell(n_trials,1);
peak_cell  = cell(n_trials,1);
range_cell = cell(n_trials,1);

for i = 1:n_trials
    dff = interp1(all_data(i).ft.xb, all_data(i).im.z', all_data(i).ft.xf);
    sum_cell{i}   = sum(dff,2);
    peak_cell{i}  = max(dff,[],2);
    range_cell{i} = max(dff,[],2) - min(dff,[],2);
end

conditions = {'closed loop','frozen','dark'};
cond_color = {'k','r','b'};

binned_sum   = nan(length(speed_x),3,n_flies);
binned_peak  = nan(length(speed_x),3,n_flies);
binned_range = nan(length(speed_x),3,n_flies);

for f = 1:n_flies
    cl_trials   = find(fly_id==f & ~is_dark);
    dark_trials = find(fly_id==f & is_dark);

    % closed loop (moving) and frozen timepoints, pooled across this fly's closed-loop trials
    cr_cl   = cat(1,cr_all{cl_trials});
    is_frz  = cat(1,freeze_idx_all{cl_trials});
    sum_cl  = cat(1,sum_cell{cl_trials});
    peak_cl = cat(1,peak_cell{cl_trials});
    rng_cl  = cat(1,range_cell{cl_trials});

    % dark timepoints, pooled across this fly's dark trials
    cr_dk   = cat(1,cr_all{dark_trials});
    sum_dk  = cat(1,sum_cell{dark_trials});
    peak_dk = cat(1,peak_cell{dark_trials});
    rng_dk  = cat(1,range_cell{dark_trials});

    binned_sum(:,1,f)   = bin_by_speed(cr_cl(~is_frz), sum_cl(~is_frz),  speed_edges, nominal_dt, min_bin_s);
    binned_sum(:,2,f)   = bin_by_speed(cr_cl(is_frz),  sum_cl(is_frz),   speed_edges, nominal_dt, min_bin_s);
    binned_sum(:,3,f)   = bin_by_speed(cr_dk,          sum_dk,           speed_edges, nominal_dt, min_bin_s);

    binned_peak(:,1,f)  = bin_by_speed(cr_cl(~is_frz), peak_cl(~is_frz), speed_edges, nominal_dt, min_bin_s);
    binned_peak(:,2,f)  = bin_by_speed(cr_cl(is_frz),  peak_cl(is_frz),  speed_edges, nominal_dt, min_bin_s);
    binned_peak(:,3,f)  = bin_by_speed(cr_dk,          peak_dk,          speed_edges, nominal_dt, min_bin_s);

    binned_range(:,1,f) = bin_by_speed(cr_cl(~is_frz), rng_cl(~is_frz),  speed_edges, nominal_dt, min_bin_s);
    binned_range(:,2,f) = bin_by_speed(cr_cl(is_frz),  rng_cl(is_frz),   speed_edges, nominal_dt, min_bin_s);
    binned_range(:,3,f) = bin_by_speed(cr_dk,          rng_dk,           speed_edges, nominal_dt, min_bin_s);
end

metric_names = {'summed dFF (z-scored, all clusters)','peak dFF (z-scored, all clusters)','peak amplitude, max-min dFF (z-scored, all clusters)'};
metric_data  = {binned_sum, binned_peak, binned_range};

for m = 1:3
    figure('Name',metric_names{m}); clf; hold on
    h = nan(3,1);
    for c = 1:3
        y = squeeze(metric_data{m}(:,c,:)); % speed bins x flies
        plotsem(speed_x,y',cond_color{c});
        h(c) = plot(speed_x,mean(y,2,'omitnan'),cond_color{c},'linewidth',2);
    end
    legend(h,conditions)
    xlabel('rotational speed (rad/s)')
    ylabel(metric_names{m})
end

%% Functions

function m = bin_by_speed(cr, y, edges, dt, min_sec)
    m = nan(length(edges)-1,1);
    for j = 1:(length(edges)-1)
        idx = cr >= edges(j) & cr < edges(j+1);
        if sum(idx)*dt > min_sec
            m(j) = mean(y(idx),'omitnan');
        end
    end
end

function h = plotsem(t,x,c)
    m = mean(x,1,'omitnan');
    s = std(x,[],1,'omitnan') ./ sqrt(sum(~isnan(x),1));
    t = reshape(t,1,[]);

    h = patch([t,fliplr(t)],[m+s,m-s],c,'FaceAlpha',.2,'EdgeColor','none');
end
