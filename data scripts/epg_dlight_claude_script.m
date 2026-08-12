%% load in data
base_dir = 'Z:\pablo\epg_dlight\';
data_file = fullfile(fileparts(fileparts(mfilename('fullpath'))),'data','epg_dlight_20260415.mat');
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
max_gap_s   = 2;       % merge freeze periods separated by a gap shorter than this: a brief dip in |r_speed| below r_thresh looks like a gap in freeze_idx even though the cue never actually resumed moving, and would otherwise split one physical freeze block into several spurious ones

freeze_idx_all  = cell(n_trials,1); % per-sample: confidently-detectable frozen samples (cue frozen AND fly rotating). used for all the per-timepoint binning above.
block_mask_all  = cell(n_trials,1); % per-sample, gap-merged version of freeze_idx: used only to find genuine block boundaries/onsets below, so a block's start isn't falsely reset by the fly briefly stopping mid-freeze.
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
    max_gap_frames    = max(1,round(max_gap_s/dt));
    c  = all_data(i).ft.cue;

    dc = smoothdata([diff(unwrap(c))/dt;0],'gaussian',pre_smooth);
    cc = smoothdata(abs(dc),'movmean',[post_smooth,0]);

    freeze_idx = cr > r_thresh & bwareaopen(cc < cue_thresh, min_block_frames);
    freeze_idx_all{i} = freeze_idx;
    frac_freeze(i) = mean(freeze_idx);

    block_mask = ~bwareaopen(~freeze_idx, max_gap_frames); % fill in gaps shorter than max_gap_frames
    block_mask_all{i} = block_mask;

    labeled = bwlabel(block_mask);
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
figure(1); clf; set(gcf,'Name','epg_dlight summary')

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

figure(2); clf; set(gcf,'Name','example freeze detection')
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

sum_cell        = cell(n_trials,1);
peak_cell       = cell(n_trials,1); % peak fluorescence: max of the hemisphere-averaged clusters (no fit) -- drives the speed/consistency plots below
peak_fit_cell   = cell(n_trials,1); % von Mises fit peak, hemisphere-averaged, kept for comparison below
peak_raw_cell   = cell(n_trials,1); % previous peak definition: max(dFF,[],2) over all raw clusters, kept for comparison below
peak_2hemi_cell = cell(n_trials,1); % alternative: no hemisphere averaging, two independent per-hemisphere fits, kept for comparison below
range_cell      = cell(n_trials,1);
dff_avg_cell    = cell(n_trials,1); % hemisphere-averaged (n_cluster/2) profile, kept below for the example fit figure
bump_fit        = struct('mu',[],'kappa',[],'amp',[],'offset',[]);
bump_fit(n_trials) = bump_fit(1);
bump_fit_2hemi  = struct('muL',[],'kappaL',[],'ampL',[],'offsetL',[], ...
                          'muR',[],'kappaR',[],'ampR',[],'offsetR',[]);
bump_fit_2hemi(n_trials) = bump_fit_2hemi(1);

for i = 1:n_trials
    dff = interp1(all_data(i).ft.xb, all_data(i).im.z', all_data(i).ft.xf);
    sum_cell{i}      = sum(dff,2);
    range_cell{i}    = max(dff,[],2) - min(dff,[],2);
    peak_raw_cell{i} = max(dff,[],2);

    % average the two PB hemispheres cluster-for-cluster (cluster c with
    % cluster c + n_cluster/2, e.g. cluster 1 with 17 of 32) before fitting
    % the bump, so the fit sees one clean profile instead of two noisy,
    % partially-redundant ones
    n_cluster  = size(dff,2);
    half       = n_cluster/2;
    dff_avg    = (dff(:,1:half) + dff(:,half+1:end)) / 2;
    alpha_half = linspace(-pi,pi,half);

    [mu,kappa,amp,offset,peak] = fit_bump_vonmises(dff_avg,alpha_half);

    dff_avg_cell{i}    = dff_avg;
    bump_fit(i).mu     = mu;
    bump_fit(i).kappa  = kappa;
    bump_fit(i).amp    = amp;
    bump_fit(i).offset = offset;

    peak_fit_cell{i} = peak; % von Mises fit peak dFF, kept for comparison below
    peak_cell{i}     = max(dff_avg,[],2); % peak fluorescence: max of the hemisphere-averaged clusters, no fit

    % alternative: skip the averaging and fit each hemisphere's 16 clusters
    % independently (same alpha_half domain for both, since it's the same
    % local position around the bridge on each side) -- this is two
    % separate, independently-located bumps instead of one, so unlike the
    % averaged fit above, the profile across all 32 raw clusters can show
    % two distinct peaks (one per hemisphere) rather than being forced into one
    [muL,kappaL,ampL,offsetL,peakL] = fit_bump_vonmises(dff(:,1:half),      alpha_half);
    [muR,kappaR,ampR,offsetR,peakR] = fit_bump_vonmises(dff(:,half+1:end),  alpha_half);

    bump_fit_2hemi(i).muL     = muL;
    bump_fit_2hemi(i).kappaL  = kappaL;
    bump_fit_2hemi(i).ampL    = ampL;
    bump_fit_2hemi(i).offsetL = offsetL;
    bump_fit_2hemi(i).muR     = muR;
    bump_fit_2hemi(i).kappaR  = kappaR;
    bump_fit_2hemi(i).ampR    = ampR;
    bump_fit_2hemi(i).offsetR = offsetR;

    peak_2hemi_cell{i} = (peakL + peakR) / 2; % combine the two independent hemisphere peaks the same way averaging would have
end

%% example timepoint: hemisphere-averaged bump and its von Mises fit
% pick the frame with the strongest population-vector bump, using im.rho
% (computed independently during preprocessing -- see e.g. cshl_poster_figures.m
% for the same pattern: rho = interp1(xb,all_data(i).im.rho,xf)) rather than
% our own fit's amplitude, so the example isn't chosen by the very method
% it's meant to illustrate

ex_rho = -inf; ex_i = 1; ex_t = 1;
for i = 1:n_trials
    rho_i = interp1(all_data(i).ft.xb, all_data(i).im.rho, all_data(i).ft.xf);
    [r,t] = max(rho_i);
    if r > ex_rho
        ex_rho = r; ex_i = i; ex_t = t;
    end
end

ex_alpha  = linspace(-pi,pi,size(dff_avg_cell{ex_i},2));
theta_grid = linspace(-pi,pi,200);
ex_fit = bump_fit(ex_i).offset(ex_t) + bump_fit(ex_i).amp(ex_t) * ...
    exp(bump_fit(ex_i).kappa(ex_t) * (cos(theta_grid - bump_fit(ex_i).mu(ex_t)) - 1));

ex_fit_peak = max(ex_fit); % the "max value of the fit" for this example timepoint
fprintf('example timepoint (trial %d, frame %d): fit peak dFF = %.3f\n', ex_i, ex_t, ex_fit_peak);

figure(3); clf; set(gcf,'Name','example bump fit'); hold on
plot(ex_alpha, dff_avg_cell{ex_i}(ex_t,:), 'ko', 'MarkerFaceColor','k', 'MarkerSize',6)
plot(theta_grid, ex_fit, 'r-', 'linewidth',1.5)
plot(bump_fit(ex_i).mu(ex_t), ex_fit_peak, 'rv', 'MarkerFaceColor','r')
xlabel('cluster angle around PB (rad)')
ylabel('hemisphere-averaged dFF (z-scored)')
legend('hemisphere-averaged clusters','von Mises fit','fit peak','Location','best')
title(sprintf('%s, frame %d',all_data(ex_i).meta,ex_t),'Interpreter','none')

%% same example timepoint, fit without hemisphere averaging (two independent bumps)
% same trial/frame as above, but now each hemisphere's 16 raw clusters are
% fit on their own (bump_fit_2hemi, computed in the loop above) instead of
% being averaged together first. plotted on the shared local-angle axis so
% you can see directly whether the two hemispheres' bumps agree or not.

half_ex = numel(ex_alpha);
dff_ex  = interp1(all_data(ex_i).ft.xb, all_data(ex_i).im.z', all_data(ex_i).ft.xf(ex_t)); % raw (unaveraged) clusters, this one frame

fitL = bump_fit_2hemi(ex_i).offsetL(ex_t) + bump_fit_2hemi(ex_i).ampL(ex_t) * ...
    exp(bump_fit_2hemi(ex_i).kappaL(ex_t) * (cos(theta_grid - bump_fit_2hemi(ex_i).muL(ex_t)) - 1));
fitR = bump_fit_2hemi(ex_i).offsetR(ex_t) + bump_fit_2hemi(ex_i).ampR(ex_t) * ...
    exp(bump_fit_2hemi(ex_i).kappaR(ex_t) * (cos(theta_grid - bump_fit_2hemi(ex_i).muR(ex_t)) - 1));
peakL_ex = max(fitL);
peakR_ex = max(fitR);

figure(4); clf; set(gcf,'Name','example bump fit, no hemisphere averaging'); hold on
plot(ex_alpha, dff_ex(1:half_ex),         'o', 'Color',[0,.45,.74],  'MarkerFaceColor',[0,.45,.74])
plot(ex_alpha, dff_ex(half_ex+1:end),     'o', 'Color',[.85,.33,.1], 'MarkerFaceColor',[.85,.33,.1])
plot(theta_grid, fitL, '-', 'Color',[0,.45,.74],  'linewidth',1.5)
plot(theta_grid, fitR, '-', 'Color',[.85,.33,.1], 'linewidth',1.5)
plot(bump_fit_2hemi(ex_i).muL(ex_t), peakL_ex, 'v', 'Color',[0,.45,.74],  'MarkerFaceColor',[0,.45,.74])
plot(bump_fit_2hemi(ex_i).muR(ex_t), peakR_ex, 'v', 'Color',[.85,.33,.1], 'MarkerFaceColor',[.85,.33,.1])
xlabel('cluster angle around PB (rad)')
ylabel('raw per-cluster dFF (z-scored)')
legend('hemisphere L clusters','hemisphere R clusters','hemisphere L fit','hemisphere R fit','Location','best')
title(sprintf('%s, frame %d -- no hemisphere averaging',all_data(ex_i).meta,ex_t),'Interpreter','none')

%% compare the von Mises peak to the old max-cluster peak, and to summed dFF
% "old" peak = max(dFF,[],2) over all raw clusters (the previous
% definition); "fit" peak = the von Mises fit peak, computed above from the
% hemisphere-averaged clusters. pooled across every behavioral sample from
% every trial; subsampled just for the scatter plots since that pool is far
% too dense to render usefully as individual points.

all_peak_fit = cat(1,peak_fit_cell{:});
all_peak_raw = cat(1,peak_raw_cell{:});
all_sum      = cat(1,sum_cell{:});

ok = ~isnan(all_peak_fit) & ~isnan(all_peak_raw) & ~isnan(all_sum);

r_fit_vs_raw = corr(all_peak_fit(ok),all_peak_raw(ok));
r_fit_vs_sum = corr(all_peak_fit(ok),all_sum(ok));

fprintf('von Mises peak vs. old max-cluster peak: r = %.3f (n = %d samples)\n', r_fit_vs_raw, sum(ok));
fprintf('von Mises peak vs. summed dFF:           r = %.3f (n = %d samples)\n', r_fit_vs_sum, sum(ok));

n_plot = min(20000,sum(ok));
plot_idx = randsample(find(ok),n_plot);

figure(5); clf; set(gcf,'Name','von Mises peak vs. old peak & summed dFF','Position',[100,100,900,420])

subplot(1,2,1); hold on
scatter(all_peak_raw(plot_idx),all_peak_fit(plot_idx),4,'filled','MarkerFaceAlpha',.15,'MarkerFaceColor','k')
axis equal
lims = [min([xlim,ylim]),max([xlim,ylim])];
plot(lims,lims,':k')
xlim(lims); ylim(lims)
xlabel('old peak: max(dFF), all raw clusters')
ylabel('von Mises fit peak')
title(sprintf('r = %.2f, n = %d',r_fit_vs_raw,sum(ok)))

subplot(1,2,2); hold on
scatter(all_sum(plot_idx),all_peak_fit(plot_idx),4,'filled','MarkerFaceAlpha',.15,'MarkerFaceColor','k')
xlabel('summed dFF, all raw clusters')
ylabel('von Mises fit peak')
title(sprintf('r = %.2f, n = %d',r_fit_vs_sum,sum(ok)))

%% compare the no-averaging (two-hemisphere) peak to the other three peak definitions
% peak_2hemi = mean of the two independent per-hemisphere fit peaks (no
% pre-averaging of the raw clusters). compared against: the hemisphere-
% averaged fit peak used above, the old max-cluster peak, and summed dFF.

all_peak_2hemi = cat(1,peak_2hemi_cell{:});
ok2 = ok & ~isnan(all_peak_2hemi);

r_2hemi_vs_fit = corr(all_peak_2hemi(ok2),all_peak_fit(ok2));
r_2hemi_vs_raw = corr(all_peak_2hemi(ok2),all_peak_raw(ok2));
r_2hemi_vs_sum = corr(all_peak_2hemi(ok2),all_sum(ok2));

fprintf('no-averaging peak vs. hemisphere-averaged fit peak: r = %.3f (n = %d samples)\n', r_2hemi_vs_fit, sum(ok2));
fprintf('no-averaging peak vs. old max-cluster peak:         r = %.3f (n = %d samples)\n', r_2hemi_vs_raw, sum(ok2));
fprintf('no-averaging peak vs. summed dFF:                   r = %.3f (n = %d samples)\n', r_2hemi_vs_sum, sum(ok2));

plot_idx2 = randsample(find(ok2),min(20000,sum(ok2)));

figure(6); clf; set(gcf,'Name','no-averaging peak vs. other peak definitions','Position',[100,100,1300,420])

subplot(1,3,1); hold on
scatter(all_peak_fit(plot_idx2),all_peak_2hemi(plot_idx2),4,'filled','MarkerFaceAlpha',.15,'MarkerFaceColor','k')
axis equal
lims2 = [min([xlim,ylim]),max([xlim,ylim])];
plot(lims2,lims2,':k')
xlim(lims2); ylim(lims2)
xlabel('fit peak, hemisphere-averaged')
ylabel('fit peak, no averaging (mean of 2 hemispheres)')
title(sprintf('r = %.2f, n = %d',r_2hemi_vs_fit,sum(ok2)))

subplot(1,3,2); hold on
scatter(all_peak_raw(plot_idx2),all_peak_2hemi(plot_idx2),4,'filled','MarkerFaceAlpha',.15,'MarkerFaceColor','k')
xlabel('old peak: max(dFF), all raw clusters')
ylabel('fit peak, no averaging')
title(sprintf('r = %.2f, n = %d',r_2hemi_vs_raw,sum(ok2)))

subplot(1,3,3); hold on
scatter(all_sum(plot_idx2),all_peak_2hemi(plot_idx2),4,'filled','MarkerFaceAlpha',.15,'MarkerFaceColor','k')
xlabel('summed dFF, all raw clusters')
ylabel('fit peak, no averaging')
title(sprintf('r = %.2f, n = %d',r_2hemi_vs_sum,sum(ok2)))

%% fluorescence vs. rotational speed: closed loop vs. frozen vs. dark (continued)

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

metric_names = {'summed dFF (z-scored, all clusters)','peak dFF (max, hemisphere-averaged)','peak amplitude, max-min dFF (z-scored, all clusters)'};
metric_data  = {binned_sum, binned_peak, binned_range};
cond_rgb     = {[0,0,0],[1,0,0],[0,0,1]};

for m = 1:3
    figure(6+m); clf; set(gcf,'Name',metric_names{m}); hold on

    % faint individual-fly traces, one per fly per condition, drawn first so the
    % mean +/- sem traces (with bin markers) sit on top and stay legible
    for c = 1:3
        y = squeeze(metric_data{m}(:,c,:)); % speed bins x flies
        plot(speed_x,y,'Color',[cond_rgb{c},.15])
    end

    h = nan(3,1);
    for c = 1:3
        y = squeeze(metric_data{m}(:,c,:));
        plotsem(speed_x,y',cond_color{c});
        h(c) = plot(speed_x,mean(y,2,'omitnan'),'-o','Color',cond_rgb{c}, ...
            'MarkerFaceColor',cond_rgb{c},'MarkerSize',4,'linewidth',2);
    end
    legend(h,conditions)
    xlabel('rotational speed (rad/s)')
    ylabel(metric_names{m})
    title(sprintf('faint lines = individual flies (n=%d)',n_flies))
end

%% per-fly consistency: paired differences in summed & peak dFF vs. rotational speed
% one figure per pairwise comparison, per metric. each fly contributes one
% faint gray difference trace (its own closed-loop/frozen/dark binned
% curves, subtracted bin-by-bin), overlaid with the mean +/- sem across flies.

diff_pairs = {[1,3],[1,2],[2,3]}; % {closed loop vs dark, closed loop vs frozen, frozen vs dark}
diff_names = {'closed loop - dark','closed loop - frozen','frozen - dark'};

consistency_metrics = {'summed dFF','peak dFF'};
consistency_data    = {binned_sum, binned_peak};

for q = 1:length(consistency_data)
    for p = 1:length(diff_pairs)
        a = diff_pairs{p}(1);
        b = diff_pairs{p}(2);
        d = squeeze(consistency_data{q}(:,a,:) - consistency_data{q}(:,b,:)); % speed bins x flies

        figure(9 + (q-1)*length(diff_pairs) + p); clf; set(gcf,'Name',sprintf('%s: %s',consistency_metrics{q},diff_names{p})); hold on
        plot(speed_x,d,'Color',[.5,.5,.5,.5])
        plotsem(speed_x,d','k');
        plot(speed_x,mean(d,2,'omitnan'),'-ok','linewidth',2,'MarkerFaceColor','k','MarkerSize',4)
        plot(xlim,[0,0],':k','linewidth',1)
        xlabel('rotational speed (rad/s)')
        ylabel(sprintf('%s, %s',consistency_metrics{q},diff_names{p}))
        title(sprintf('faint lines = individual flies (n=%d)',n_flies))
    end
end

%% fluorescence deficit over the course of a freeze block, aligned to onset
% for every detected freeze block, align a window of data to the block's
% starting edge (t=0). at each aligned sample we know the instantaneous
% rotational speed, so we look up "what summed dFF would be expected if
% the cue were still moving at that speed" and compare it to what was
% actually observed. one trace is plotted per rotational speed bin.
%
% the "expected if moving" value for a given block uses THAT FLY's own
% closed-loop curve (binned_sum(:,1,fly)), not the population-average
% curve: flies differ a lot in their overall closed-loop gain (see the
% spread of faint lines in the sum/peak/range figures above), and flies
% with more/longer freeze blocks are not a random sample of flies, so
% comparing against the population average leaves a persistent
% non-zero offset even before the freeze starts. matching each block to
% its own fly's reference removes that. blocks are also first averaged
% within each fly before averaging across flies (equal per-fly weight,
% consistent with how the reference curves themselves were built),
% rather than pooling every block equally.

rel_t       = -15:.25:25;  % time relative to freeze onset (s)
min_flies   = 5;           % require at least this many flies contributing at a given (time, speed bin) to trust it

block_list = zeros(0,2); % [trial_idx, block_label]
for i = 1:n_trials
    if is_dark(i) || isnan(n_blocks(i)) || n_blocks(i)==0; continue; end
    for b = 1:n_blocks(i)
        block_list = [block_list; i, b]; %#ok<AGROW>
    end
end

n_align      = size(block_list,1);
diff_mat     = nan(n_align,length(rel_t));
speed_mat    = nan(n_align,length(rel_t));
fly_of_block = fly_id(block_list(:,1));

for k = 1:n_align
    i = block_list(k,1);
    b = block_list(k,2);
    ref_fly = squeeze(binned_sum(:,1,fly_of_block(k))); % this fly's own closed-loop reference curve

    % onset is found from block_mask_all (gap-merged), not freeze_idx_all,
    % so a momentary dip in rotation speed mid-freeze doesn't get mistaken
    % for a return to closed loop followed by a fresh onset
    labeled    = bwlabel(block_mask_all{i});
    onset_idx  = find(labeled==b,1,'first');
    onset_time = all_data(i).ft.xf(onset_idx);
    xq = onset_time + rel_t;

    cr_q    = interp1(all_data(i).ft.xf, cr_all{i},          xq);
    amp_q   = interp1(all_data(i).ft.xf, sum_cell{i},        xq);
    label_q = interp1(all_data(i).ft.xf, double(labeled), xq, 'nearest');

    % keep only samples that are genuinely on the correct side of this
    % specific freeze block: before t=0, only keep samples where the cue
    % is not frozen at all (label==0); at/after t=0, only keep samples
    % still inside THIS block's own frozen span (label==b) -- e.g. if the
    % cue re-starts moving 2s into a block, later timepoints get dropped
    % rather than silently counted as "observed while frozen"
    keep = false(1,length(rel_t));
    keep(rel_t<0)  = label_q(rel_t<0)==0;
    keep(rel_t>=0) = label_q(rel_t>=0)==b;

    speed_mat(k,:) = cr_q;
    diff_mat(k,:)  = amp_q - interp1(speed_x, ref_fly, cr_q); % observed - expected
    diff_mat(k,~keep)  = nan;
    speed_mat(k,~keep) = nan;
end

fprintf('n freeze blocks aligned: %d\n', n_align);

figure(16); clf; set(gcf,'Name','freeze-onset aligned fluorescence deficit','Position',[100,100,900,550]); hold on
cmap = parula(length(speed_x));
h = []; leg_labels = {};
smooth_win = 5; % ~1.25s of smoothing (rel_t step is .25s), for display only

for j = 1:(length(speed_edges)-1)
    bin_idx = speed_mat >= speed_edges(j) & speed_mat < speed_edges(j+1);
    trace   = nan(1,length(rel_t));
    for t = 1:length(rel_t)
        col_idx = bin_idx(:,t);
        flies_here = unique(fly_of_block(col_idx));
        fly_vals = nan(length(flies_here),1);
        for ff = 1:length(flies_here)
            rows = col_idx & fly_of_block==flies_here(ff);
            fly_vals(ff) = mean(diff_mat(rows,t),'omitnan');
        end
        if sum(~isnan(fly_vals)) >= min_flies
            trace(t) = mean(fly_vals,'omitnan');
        end
    end
    if all(isnan(trace)); continue; end

    nan_mask = isnan(trace);
    trace = smoothdata(trace,'movmean',smooth_win,'omitnan');
    trace(nan_mask) = nan; % don't let smoothing bridge over bins that failed the min_flies cutoff

    h(end+1)          = plot(rel_t,trace,'Color',cmap(j,:),'linewidth',1.5); %#ok<SAGROW>
    leg_labels{end+1} = sprintf('%.2f-%.2f rad/s',speed_edges(j),speed_edges(j+1)); %#ok<SAGROW>
end

y = ylim;
plot([0,0],y,':k','linewidth',1)
ylim(y)
plot(xlim,[0,0],':k','linewidth',1)
legend(h,leg_labels,'Location','eastoutside')
xlabel('time from freeze onset (s)')
ylabel('summed dFF: observed - closed-loop expectation')
title('freeze onset alignment (negative = below closed-loop expectation)')

%% save all figures as PDF
% one PDF per figure, named after that figure's 'Name' (set above via
% set(gcf,'Name',...)), lowercased with non-alphanumeric runs collapsed to
% a single underscore -- e.g. 'peak dFF: closed loop - frozen' ->
% peak_dff_closed_loop_frozen.pdf

save_dir = fullfile(fileparts(fileparts(mfilename('fullpath'))),'ugly_figures','epg_dlight');
if ~exist(save_dir,'dir'); mkdir(save_dir); end

for fig_num = 1:16
    fh = findobj('Type','figure','Number',fig_num);
    if isempty(fh); continue; end

    file_name = lower(get(fh,'Name'));
    file_name = regexprep(file_name,'[^a-z0-9]+','_');
    file_name = regexprep(file_name,'^_+|_+$','');

    exportgraphics(fh, fullfile(save_dir,[file_name,'.pdf']))
    fprintf('saved figure %d -> %s.pdf\n', fig_num, file_name);
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

    valid = ~isnan(m) & ~isnan(s); % drop bins with no data so a trailing NaN doesn't blank the whole patch
    t = t(valid); m = m(valid); s = s(valid);

    h = patch([t,fliplr(t)],[m+s,fliplr(m-s)],c,'FaceAlpha',.2,'EdgeColor','none');
end

function [mu,kappa,amp,offset,peak] = fit_bump_vonmises(y,alpha)
% fit_bump_vonmises  per-row (per-timepoint) von Mises "bump" fit:
%   y(theta) = offset + amp * exp(kappa*(cos(theta-mu)-1))
% where kappa is the concentration parameter (kappa = 1/width: larger kappa
% -> narrower bump) and amp/offset set its height/baseline.
%
%   y      [n_time x n_cluster] signal at each of n_cluster fixed angles
%   alpha  [1 x n_cluster] angle (rad) of each column of y
%
% mu and kappa come from the population vector (weighted circular mean) of
% the non-negative part of y, converted to kappa with the same point
% estimate circ_kappa.m uses (Fisher, 1993, eq. p.88) -- the resultant-length
% -> kappa formula assumes non-negative "counts", so only the negative part
% of the (z-scored) signal is clipped for this step; amp/offset below are
% still fit to the real, signed data. Given mu/kappa, the model is linear in
% [amp,offset], so that fit is closed-form (exact 2-parameter least squares)
% rather than an iterative optimizer run at every timepoint of a whole
% experiment's worth of data.

n     = numel(alpha);
alpha = reshape(alpha,1,n);

w  = max(y,0);
sw = sum(w,2);
cx = w*cos(alpha)';
cy = w*sin(alpha)';

mu = atan2(cy,cx);
R  = zeros(size(sw));
R(sw>0) = hypot(cx(sw>0),cy(sw>0)) ./ sw(sw>0);

kappa = kappa_from_R(R,n);

basis = exp(kappa.*(cos(alpha-mu)-1)); % n_time x n_cluster

Sb  = sum(basis,2);
Sbb = sum(basis.^2,2);
Sy  = sum(y,2);
Sby = sum(basis.*y,2);

det  = Sbb.*n - Sb.^2;
flat = abs(det) < 1e-9; % basis is ~constant (kappa~0, no clear direction) -> amp/offset aren't separable

amp    = nan(size(det));
offset = nan(size(det));
amp(~flat)    = (Sby(~flat).*n - Sb(~flat).*Sy(~flat)) ./ det(~flat);
offset(~flat) = (Sbb(~flat).*Sy(~flat) - Sb(~flat).*Sby(~flat)) ./ det(~flat);

amp(flat)    = 0;
offset(flat) = mean(y(flat,:),2);

theta_grid = linspace(-pi,pi,64);
y_grid = offset + amp.*exp(kappa.*(cos(theta_grid-mu)-1));
peak   = max(y_grid,[],2);
end

function kappa = kappa_from_R(R,N)
% vectorized version of circ_kappa.m's point estimate (Fisher, 1993, eq.
% p.88) for a fixed sample size N shared across every row of R

R = min(R,1-1e-4); % avoid the R->1 singularity in the R>=.85 branch below

kappa = zeros(size(R));
lo  = R < .53;
mid = R >= .53 & R < .85;
hi  = R >= .85;

kappa(lo)  = 2*R(lo) + R(lo).^3 + 5*R(lo).^5/6;
kappa(mid) = -.4 + 1.39*R(mid) + .43./(1-R(mid));
kappa(hi)  = 1./(R(hi).^3 - 4*R(hi).^2 + 3*R(hi));

if N < 15
    small = kappa < 2;
    kappa(small)  = max(kappa(small) - 2./(N*kappa(small)), 0);
    kappa(~small) = (N-1)^3 .* kappa(~small) ./ (N^3+N);
end
end
