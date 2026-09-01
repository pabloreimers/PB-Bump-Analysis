%% lpsp_rnai_claude_v2
% Parameter-optimization pass on the LPsP-RNAi bump-mobility pipeline
% (lpsp_rnai_claude.m). Four steps, each with its own figure(s):
%   1) load the re-registered dataset (lpsp_rnai_joint_nosmooth_reg_*.mat)
%   2) re-estimate bump position from a 5-FRAME MOVING AVERAGE of raw
%      per-glomerulus fluorescence (im.f), z-scored, then PVA'd
%   3) sweep candidate heading-smoothing windows and pick whichever
%      MINIMIZES offset variability (circ_var(circ_dist(mu,heading))) in
%      the empty-driver controls (empty>th, empty>vglut), closed-loop
%      trials only
%   4) recompute walking-bout bump mobility with that pipeline, then
%      compare weighting the per-bout regression by fly path length
%      (mov_speed) instead of bout duration (dur), picking whichever
%      gives the better goodness-of-fit, and use the winner for a final
%      group plot
%
% addpath is explicit (not assumed already on the caller's MATLABPATH) so
% this script is self-contained when run non-interactively.
addpath(fullfile(pwd,'circ_stats'));

%% 1) load data
data_dir    = fullfile('data');
source_file = 'lpsp_rnai_joint_nosmooth_reg_20260831.mat';

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype (driver x target) per trial, from all_data.meta
% identical parsing to lpsp_rnai_claude.m -- see that script's header
% comment for the full rationale (blinded-trial handling, implicit
% TH-RNAi default, etc.); repeated verbatim here rather than shared, since
% the two scripts are meant to stand alone.
driver = cell(n_trials,1);
target = cell(n_trials,1);

for i = 1:n_trials
    [fly_seg,trial_seg] = meta_fly_and_trial_seg(all_data(i).meta);
    combo = lower(regexprep([fly_seg,trial_seg],'[_\s]',''));

    if contains(combo,'lpsp')
        driver{i} = 'lpsp';
    elseif contains(combo,'empty')
        driver{i} = 'empty';
    else
        driver{i} = '';
        warning('trial %d: could not determine driver from meta "%s"', i, all_data(i).meta);
    end

    if contains(combo,'thrnai')
        target{i} = 'th';
    elseif contains(combo,'vglut')
        target{i} = 'vglut';
    elseif contains(combo,'mcherry')
        target{i} = 'mcherry';
    else
        target{i} = 'th';
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
fprintf('%d trials -> %d flies\n', n_trials, n_flies);

geno_order = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};
cond_label = {'closed loop','dark'};

%% 2) re-estimate bump position: 5-frame moving average of im.f, then zscore, then PVA
% movmean(...,5) is a plain 5-SAMPLE centered moving average (2 frames
% either side) on the imaging-frame axis -- deliberately simpler than
% lpsp_rnai_claude.m's gaussian/time-converted smoothing, per this
% script's brief. Because it's specified as a frame count (not seconds),
% it is NOT converted per-trial via each trial's own imaging frame period
% -- unlike every time-based smoothing window elsewhere in this repo's
% claude scripts, 5 frames here means exactly 5 samples on every trial,
% regardless of that trial's own frame rate.
bump_smooth_frames = 5;

chunk_mu_new  = cell(n_trials,1);
chunk_rho_new = cell(n_trials,1);
chunk_fz      = cell(n_trials,1);
for i = 1:n_trials
    [chunk_mu_new{i},chunk_rho_new{i},chunk_fz{i}] = trial_pva_movmean(all_data(i),bump_smooth_frames);
end

%% figure: sanity check the new bump-position estimate against the upstream im.mu fit, one example trial
% picks the trial with the highest mean upstream rho (im.rho) as a
% "should look clean" example -- not cherry-picked for agreement, just
% for having an actual bump to look at.
mean_rho_upstream = arrayfun(@(t) mean(t.im.rho,'omitnan'), all_data);
[~,example_i] = max(mean_rho_upstream);

trial = all_data(example_i);
xf = trial.ft.xf;
n_im = size(trial.im.f,2);
alpha_disp = unwrap(trial.im.alpha(:));

mu_old = unwrap(trial.im.mu(:));

figure(1); clf
set(gcf,'Name','sanity check: new bump position vs upstream im.mu','Position',[100,100,1000,400])
hold on
imagesc(1:n_im,alpha_disp,chunk_fz{example_i})
set(gca,'YDir','normal','XLim',[1,n_im],'YLim',[alpha_disp(1),alpha_disp(end)])
mu_old_plot = wrap_break(mu_old);
plot(1:n_im,mu_old_plot,'c','LineWidth',1)
plot(1:n_im,mu_old_plot+2*pi,'c','LineWidth',1)
mu_new_plot = wrap_break(chunk_mu_new{example_i});
plot(1:n_im,mu_new_plot,'w','LineWidth',1)
plot(1:n_im,mu_new_plot+2*pi,'w','LineWidth',1)
cb = colorbar; ylabel(cb,'z-score (5-frame moving average)')
xlabel('imaging frame'); ylabel('PB angle (rad)')
title(sprintf('trial %d (highest upstream mean rho=%.2f): upstream im.mu (cyan) vs 5-frame-movmean-PVA mu (white)',example_i,mean_rho_upstream(example_i)))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig1_bump_sanity_check.png'),'Resolution',150)

%% 3) sweep heading-smoothing window, minimize offset variability in empty controls, closed loop only
% offset variability = circ_var(circ_dist(mu,heading)), after subtracting
% each fly's own mean offset (a fly can have a fixed, arbitrary rotational
% bias between mu's zero and the visual cue's zero -- that's not "noise",
% only deviation AROUND that fixed offset is) -- same construction as
% lpsp_rnai_joint.m's own "extract the offset variability" section
% (circ_dist(-cue,mu), demean via atan2(mean(sin),mean(cos)), restrict to
% abs(r_speed)>r_thresh and rho>rho_thresh), with the SAME r_thresh=.5,
% rho_thresh=.5 that script used for this specific metric.
%
% mu is fixed at this point (the 5-frame-movmean estimate from step 2) --
% only the heading (-ft.cue) smoothing window is swept. Heading is
% smoothed with a plain moving average too (movmean), for the same
% "simple, literal" reason as step 2, over a window given in SECONDS and
% converted to samples via each trial's own fictrac sample period
% (median(diff(xf))) -- fictrac sampling isn't at a fixed rate across
% every trial in this dataset (confirmed directly: 37800 vs 36006 samples
% for nominally similar trial durations), so a shared raw sample count
% would not mean the same real-time window on every trial.
r_thresh   = .5;
rho_thresh = .5;
heading_smooth_candidates_s = [0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5];

opt_geno   = {'empty>th','empty>vglut'};
opt_trials = find(ismember(genotype,opt_geno) & ~is_dark);
opt_flies  = unique(fly_num(opt_trials));
fprintf('\noptimizing heading smoothing on %d empty-control flies (%d closed-loop trials)\n', numel(opt_flies), numel(opt_trials));

n_cand = numel(heading_smooth_candidates_s);
fly_offset_var = nan(numel(opt_flies),n_cand);

for c = 1:n_cand
    win_s = heading_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        f = opt_flies(k);
        trial_list = opt_trials(fly_num(opt_trials)==f);
        pooled_offset = [];
        for ti = trial_list'
            trial = all_data(ti);
            xf = trial.ft.xf;
            dt = median(diff(xf));
            win_samples = max(1,round(win_s/dt));

            heading_smooth = movmean(unwrap(-trial.ft.cue(:)),win_samples);

            n_im = numel(chunk_mu_new{ti});
            xb = linspace(xf(1),xf(end),n_im)';
            mu_i  = interp1(xb,unwrap(chunk_mu_new{ti}(:)),xf,'linear','extrap');
            rho_i = interp1(xb,chunk_rho_new{ti}(:),xf,'linear','extrap');

            offset = circ_dist(heading_smooth,mu_i);
            valid  = abs(trial.ft.r_speed(:)) > r_thresh & rho_i > rho_thresh;
            pooled_offset = [pooled_offset; offset(valid)]; %#ok<AGROW>
        end
        if numel(pooled_offset) < 10
            continue % not enough valid samples for this fly to trust a circ_var estimate
        end
        offset_mean = atan2(mean(sin(pooled_offset)),mean(cos(pooled_offset)));
        offset_demeaned = circ_dist(pooled_offset,offset_mean);
        fly_offset_var(k,c) = circ_var(offset_demeaned);
    end
end

mean_offset_var = mean(fly_offset_var,1,'omitnan');
sem_offset_var  = std(fly_offset_var,0,1,'omitnan') ./ sqrt(sum(~isnan(fly_offset_var),1));
[~,opt_c] = min(mean_offset_var);
heading_smooth_opt_s = heading_smooth_candidates_s(opt_c);
fprintf('optimal heading-smoothing window: %.2fs (mean offset circ_var=%.4f)\n', heading_smooth_opt_s, mean_offset_var(opt_c));

figure(2); clf
set(gcf,'Name','heading-smoothing optimization','Position',[100,100,700,500])
hold on
errorbar(heading_smooth_candidates_s,mean_offset_var,sem_offset_var,'-ok','MarkerFaceColor','k')
plot(heading_smooth_opt_s,mean_offset_var(opt_c),'o','MarkerSize',12,'Color','r','LineWidth',2)
xlabel('heading moving-average window (s)')
ylabel('offset variability: circ\_var(circ\_dist(mu,heading))')
title(sprintf('heading-smoothing sweep, empty>th + empty>vglut, closed loop (n=%d flies) -- optimum = %.2fs',numel(opt_flies),heading_smooth_opt_s))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig2_heading_smoothing_sweep.png'),'Resolution',150)

% bout-detection constants used by both the mu-smoothing sweep below and
% the full-dataset bout recomputation in step 4 -- declared once, here,
% since both need them. heading is smoothed at heading_smooth_opt_s
% (step 3, just above), applied to r_speed (a directly-measured velocity)
% rather than to cue position -- same reasoning as lpsp_rnai_claude.m's
% own bout functions: r_speed doesn't need an unwrap (it's not circular),
% so it can be smoothed and used directly for both walking-bout detection
% and the path-length integral.
turn_thresh        = 0.25; % rad/s
max_gap_s          = 0.5;
min_walking_s      = 0.5;
bout_rho_thresh    = 0.4; % a bout's own mean bump vector strength must exceed this to be kept

%% 3b) sweep an ADDITIONAL mu-smoothing window, target mobility ratio == 1 in empty controls, closed loop only
% step 2 fixed the bump position at a 5-frame movmean of raw fluorescence
% (chunk_mu_new) -- this sweeps a SEPARATE, additional smoothing stage
% applied directly to that already-computed bump POSITION trace (not a
% re-sweep of bump_smooth_frames itself), on the theory that a healthy
% control fly (empty>th, empty>vglut) with no manipulation should have its
% bump track heading close to 1:1 in closed loop -- so whichever
% additional mu smoothing brings the CONTROL flies' own mobility ratio
% closest to 1 is treated as correcting residual bump-tracking noise, not
% real biology, and gets applied dataset-wide below.
%
% Uses the SAME opt_flies/opt_trials as step 3 (empty>th + empty>vglut,
% closed loop) and the SAME heading_smooth_opt_s already found. The
% per-fly regression here is deliberately UNWEIGHTED, independent of step
% 5's later (separate) comparison of weighting schemes -- evaluating this
% sweep would otherwise depend on which weighting scheme wins, and that
% comparison itself runs on bouts built from whichever mu-smoothing wins
% here, a circular dependency; unweighted also happens to win step 5's
% comparison on every run of this pipeline so far.
mu_smooth_candidates_s = [0, 0.05, 0.1, 0.2, 0.35, 0.5, 0.75, 1, 1.5, 2, 3, 5];
min_bouts_for_sweep = 5;

n_mu_cand = numel(mu_smooth_candidates_s);
fly_ratio_sweep = nan(numel(opt_flies),n_mu_cand);

for c = 1:n_mu_cand
    win_s = mu_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        f = opt_flies(k);
        trial_list = opt_trials(fly_num(opt_trials)==f);
        mov_mu_f = []; mov_speed_f = [];
        for ti = trial_list'
            trial = all_data(ti);
            xf = trial.ft.xf;
            n_im = numel(chunk_mu_new{ti});
            dt_im = (xf(end)-xf(1)) / (n_im-1);
            win_samples_im = max(1,round(win_s/dt_im));
            mu_smoothed = movmean(unwrap(chunk_mu_new{ti}(:)),win_samples_im);

            [mov_mu_ti,mov_speed_ti] = trial_walking_bouts_v2( ...
                trial,mu_smoothed,chunk_rho_new{ti},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);
            mov_mu_f    = [mov_mu_f; mov_mu_ti]; %#ok<AGROW>
            mov_speed_f = [mov_speed_f; mov_speed_ti]; %#ok<AGROW>
        end
        if numel(mov_mu_f) < min_bouts_for_sweep
            continue
        end
        fly_ratio_sweep(k,c) = mov_speed_f \ mov_mu_f; % unweighted, through-origin
    end
end

mean_ratio_sweep = mean(fly_ratio_sweep,1,'omitnan');
sem_ratio_sweep  = std(fly_ratio_sweep,0,1,'omitnan') ./ sqrt(sum(~isnan(fly_ratio_sweep),1));
[~,opt_mu_c] = min(abs(mean_ratio_sweep-1));
mu_smooth_opt_s = mu_smooth_candidates_s(opt_mu_c);
fprintf('optimal additional mu-smoothing window: %.2fs (mean control-fly ratio=%.3f)\n', mu_smooth_opt_s, mean_ratio_sweep(opt_mu_c));

figure(10); clf
set(gcf,'Name','mu-smoothing optimization','Position',[100,100,900,500])
hold on
errorbar(mu_smooth_candidates_s,mean_ratio_sweep,sem_ratio_sweep,'-ok','MarkerFaceColor','k')
plot(mu_smooth_opt_s,mean_ratio_sweep(opt_mu_c),'o','MarkerSize',12,'Color','r','LineWidth',2)
yline(1,':k')
xlabel('additional mu moving-average window (s)')
ylabel('mobility ratio (unweighted, mov\_mu ~ mov\_speed)')
title(sprintf('mu-smoothing sweep, empty>th + empty>vglut, closed loop (n=%d flies) -- optimum = %.2fs (ratio=%.3f)', ...
    numel(opt_flies),mu_smooth_opt_s,mean_ratio_sweep(opt_mu_c)))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig10_mu_smoothing_sweep.png'),'Resolution',150)

% apply the winning additional smoothing dataset-wide -- everything from
% here on (step 4's bout recomputation, and every downstream figure) uses
% chunk_mu_final, not the raw step-2 chunk_mu_new.
chunk_mu_final = cell(n_trials,1);
for i = 1:n_trials
    xf = all_data(i).ft.xf;
    n_im = numel(chunk_mu_new{i});
    dt_im = (xf(end)-xf(1)) / (n_im-1);
    win_samples_im = max(1,round(mu_smooth_opt_s/dt_im));
    chunk_mu_final{i} = movmean(unwrap(chunk_mu_new{i}(:)),win_samples_im);
end

%% 4) recompute walking-bout bump mobility with the optimized pipeline
chunk_bout_mu    = cell(n_trials,1);
chunk_bout_speed = cell(n_trials,1);
chunk_bout_dur   = cell(n_trials,1);
for i = 1:n_trials
    [chunk_bout_mu{i},chunk_bout_speed{i},chunk_bout_dur{i}] = trial_walking_bouts_v2( ...
        all_data(i),chunk_mu_final{i},chunk_rho_new{i},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);
end

%% 5) optimize the per-bout regression weighting: duration vs. fly path length
% pools EVERY bout from EVERY fly x light-condition instance with at
% least min_bouts_for_r2 bouts (no genotype/light restriction here --
% unlike step 3, this is about the general shape of the bump/heading
% path-length relationship, not genotype comparison) and, per such
% instance, fits mov_mu ~ mov_speed through the origin under each
% candidate weighting scheme, then reports each scheme's own weighted R^2
% (1 - weighted SS_resid / weighted SS_total around the weighted mean) --
% averaged across instances as the scheme's overall goodness-of-fit.
min_bouts_for_r2 = 10;

weight_schemes = struct('name',{'unweighted','duration','path length'}, ...
                         'fun',{@(dur,speed) ones(size(dur)), @(dur,speed) sqrt(dur), @(dur,speed) sqrt(speed)});

group_defs_all = struct('geno',{},'dark',{},'label',{});
for gi = 1:numel(geno_order)
    for ci = 1:numel(cond_label)
        group_defs_all(end+1) = struct('geno',geno_order{gi},'dark',ci-1,'label',sprintf('%s (%s)',geno_order{gi},cond_label{ci})); %#ok<SAGROW>
    end
end

r2_by_scheme = cell(1,numel(weight_schemes));
for s = 1:numel(weight_schemes)
    r2_by_scheme{s} = [];
end

for gd = 1:numel(group_defs_all)
    rows = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mov_mu    = cat(1,chunk_bout_mu{trial_list});
        mov_speed = cat(1,chunk_bout_speed{trial_list});
        dur       = cat(1,chunk_bout_dur{trial_list});
        if numel(mov_mu) < min_bouts_for_r2
            continue
        end
        for s = 1:numel(weight_schemes)
            w = weight_schemes(s).fun(dur,mov_speed);
            ratio = (w.*mov_speed) \ (w.*mov_mu);
            pred  = ratio*mov_speed;
            wmean = sum(w.*mov_mu)/sum(w);
            ss_res = sum(w.*(mov_mu-pred).^2);
            ss_tot = sum(w.*(mov_mu-wmean).^2);
            r2 = 1 - ss_res/ss_tot;
            r2_by_scheme{s}(end+1) = r2; %#ok<AGROW>
        end
    end
end

fprintf('\n=== per-bout regression weighting comparison (n=%d fly x light-condition instances, >=%d bouts each) ===\n', numel(r2_by_scheme{1}), min_bouts_for_r2);
mean_r2 = nan(1,numel(weight_schemes));
for s = 1:numel(weight_schemes)
    mean_r2(s) = mean(r2_by_scheme{s},'omitnan');
    fprintf('  %-15s mean R^2 = %.3f\n', weight_schemes(s).name, mean_r2(s));
end
[~,best_s] = max(mean_r2);
fprintf('winning weighting scheme: %s\n', weight_schemes(best_s).name);

figure(3); clf
set(gcf,'Name','regression weighting comparison','Position',[100,100,700,500])
hold on
for s = 1:numel(weight_schemes)
    y = r2_by_scheme{s};
    x = s*ones(size(y));
    jitter = (rand(size(y))-.5)*.3;
    scatter(x+jitter,y,15,'k','filled','MarkerFaceAlpha',.2)
    errorbar(s,mean(y,'omitnan'),std(y,'omitnan')/sqrt(numel(y)),'or','MarkerFaceColor','r','LineWidth',2)
end
xticks(1:numel(weight_schemes)); xticklabels({weight_schemes.name})
ylabel('per-fly-x-condition weighted R^2 (mov\_mu ~ mov\_speed, through origin)')
title(sprintf('regression weighting comparison -- winner: %s',weight_schemes(best_s).name))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig3_weighting_comparison.png'),'Resolution',150)

best_weight_fun = weight_schemes(best_s).fun;

%% figure: how walking bouts are defined, a couple of example trials, with each one's own regression fit
% left column: same trace-vs-threshold-vs-detected-bout view as before
% (re-running trial_walking_bouts_v2 with its r_speed_smooth/is_walking
% outputs kept, purely for display). Right column: THIS TRIAL's own bouts
% (chunk_bout_mu/chunk_bout_speed{ti}, computed once already, back in
% step 4) as a scatter, with the through-origin fit under the WINNING
% weighting scheme from step 5 (best_weight_fun) overlaid, plus its R^2
% and adjusted R^2 (standard 1-predictor correction: adj = 1-(1-R^2)*(n-1)/(n-p-1),
% p=1 -- the same DOF penalty regardless of how many bouts happen to be
% in this one trial, so a fit from only a handful of bouts is
% appropriately discounted relative to one from many). Example 1 is the
% same "highest upstream rho" trial as figure 1's bump sanity check;
% example 2 is whichever trial ended up with the most detected bouts, as
% a busier contrasting example.
[~,busiest_i] = max(cellfun(@numel,chunk_bout_mu));
example_trials = [example_i, busiest_i];

figure(5); clf
set(gcf,'Name','how walking bouts are defined','Position',[100,100,1400,600])
for e = 1:numel(example_trials)
    ti = example_trials(e);
    [~,~,~,r_speed_smooth_ex,is_walking_ex] = trial_walking_bouts_v2( ...
        all_data(ti),chunk_mu_final{ti},chunk_rho_new{ti},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);
    t = all_data(ti).ft.xf - all_data(ti).ft.xf(1);

    subplot(2,2,(e-1)*2+1); hold on
    bout_edges = find(diff([0;is_walking_ex;0])~=0);
    for b = 1:2:numel(bout_edges)-1
        patch(t([bout_edges(b),bout_edges(b+1)-1,bout_edges(b+1)-1,bout_edges(b)]), ...
              [-1,-1,1,1]*max(abs(all_data(ti).ft.r_speed)),[1,1,0],'FaceAlpha',.25,'EdgeColor','none')
    end
    plot(t,all_data(ti).ft.r_speed,'Color',[.7,.7,.7],'LineWidth',.5)
    plot(t,r_speed_smooth_ex,'k','LineWidth',1)
    yline(turn_thresh,'--r'); yline(-turn_thresh,'--r')
    xlim([t(1),t(end)])
    xlabel('time (s)'); ylabel('r\_speed (rad/s)')
    n_bouts_ex = numel(chunk_bout_mu{ti});
    title(sprintf('trial %d, %s -- %d bouts (raw=gray, smoothed=black, thresh=\\pm%.2f)', ...
        ti,fly_short_label(fly_list{fly_num(ti)}),n_bouts_ex,turn_thresh),'FontSize',9)

    mov_mu_ex    = chunk_bout_mu{ti};
    mov_speed_ex = chunk_bout_speed{ti};
    dur_ex       = chunk_bout_dur{ti};
    w_ex     = best_weight_fun(dur_ex,mov_speed_ex);
    ratio_ex = (w_ex.*mov_speed_ex) \ (w_ex.*mov_mu_ex);
    pred_ex  = ratio_ex*mov_speed_ex;
    wmean_ex = sum(w_ex.*mov_mu_ex)/sum(w_ex);
    ss_res   = sum(w_ex.*(mov_mu_ex-pred_ex).^2);
    ss_tot   = sum(w_ex.*(mov_mu_ex-wmean_ex).^2);
    r2_ex    = 1 - ss_res/ss_tot;
    n_ex     = numel(mov_mu_ex);
    r2_adj_ex = 1 - (1-r2_ex)*(n_ex-1)/(n_ex-1-1); % p=1 (one slope parameter)

    subplot(2,2,(e-1)*2+2); hold on
    scatter(mov_speed_ex,mov_mu_ex,25,dur_ex,'filled')
    xl = xlim; xl(1) = 0;
    plot(xl,ratio_ex*xl,'-k','LineWidth',1.5)
    plot(xl,xl,':','Color',[.6,.6,.6])
    xlim(xl)
    cb = colorbar; ylabel(cb,'bout duration (s)')
    xlabel('heading path length (rad)'); ylabel('bump path length (rad)')
    title(sprintf('%s fit: slope=%.2f, R^2=%.2f, adj R^2=%.2f',weight_schemes(best_s).name,ratio_ex,r2_ex,r2_adj_ex),'FontSize',9)
end
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig5_bout_definition.png'),'Resolution',150)

%% figure: imagesc + bump + heading overlay, for the high-slope ("busiest") example trial
% same imagesc/overlay convention as lpsp_rnai_claude.m's diagnostic
% figures: alpha unwrapped to the continuous -pi:~3pi range (im.alpha
% repeats the same -pi:pi sequence once per PB hemisphere), bump/heading
% each drawn twice (as-is, and shifted +2*pi) so the overlay tracks both
% hemisphere bands, not just one. Bump = chunk_mu_final{busiest_i} (step
% 2's PVA estimate plus the additional mu-smoothing from step 3b).
% Heading = -ft.cue,
% smoothed at heading_smooth_opt_s (step 3's optimum) on the fictrac
% timebase, THEN interpolated onto the image timebase for display --
% smoothing must happen before interpolation/wrapping, never after, or it
% would smear across the circular variable's false -pi/pi discontinuities.
% slope/R^2 recomputed HERE, fresh, from this trial's OWN current
% chunk_bout_mu/chunk_bout_speed/chunk_bout_dur (post rho-filter, under
% whichever weighting scheme is currently winning) -- NOT copy-pasted
% from whatever figure 5 happened to show on an earlier run. Those bout
% arrays are pipeline outputs that change whenever an upstream parameter
% (like bout_rho_thresh) changes, and busiest_i itself is picked by
% max bout count, so both the identity of "the busiest trial" and its own
% fit can silently shift between runs -- recomputing here instead of
% hardcoding a remembered number is what keeps this title honest.
mov_mu_busiest    = chunk_bout_mu{busiest_i};
mov_speed_busiest = chunk_bout_speed{busiest_i};
dur_busiest       = chunk_bout_dur{busiest_i};
w_busiest     = best_weight_fun(dur_busiest,mov_speed_busiest);
ratio_busiest = (w_busiest.*mov_speed_busiest) \ (w_busiest.*mov_mu_busiest);
pred_busiest  = ratio_busiest*mov_speed_busiest;
wmean_busiest = sum(w_busiest.*mov_mu_busiest)/sum(w_busiest);
r2_busiest = 1 - sum(w_busiest.*(mov_mu_busiest-pred_busiest).^2)/sum(w_busiest.*(mov_mu_busiest-wmean_busiest).^2);

plot_bump_heading_overlay(6,all_data(busiest_i),chunk_mu_final{busiest_i},chunk_fz{busiest_i},heading_smooth_opt_s, ...
    sprintf('trial %d, %s -- bump position (white) vs %.2fs-smoothed heading (-ft.cue, cyan) -- this is the busiest fly in figure 5 (slope=%.2f, R^2=%.2f)', ...
        busiest_i,fly_short_label(fly_list{fly_num(busiest_i)}),heading_smooth_opt_s,ratio_busiest,r2_busiest), ...
    fullfile('ugly_figures','exports','v2_fig6_busiest_trial_overlay.png'))

%% 6) final bump-mobility group plot, using the optimized pipeline + winning weighting scheme
min_bouts_per_fly = 20;
fly_mob_ratio = [];
fly_adj_r2    = []; % this fly's own goodness-of-fit (adjusted R^2, same 1-parameter correction as figure 5/6), not just its slope
cat_x_mob     = [];
fly_mob_fly   = []; % fly_num index for this entry -- lets later code trace an entry back to an actual fly/trial
for gd = 1:numel(group_defs_all)
    rows = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mov_mu    = cat(1,chunk_bout_mu{trial_list});
        mov_speed = cat(1,chunk_bout_speed{trial_list});
        dur       = cat(1,chunk_bout_dur{trial_list});
        if numel(mov_mu) >= min_bouts_per_fly
            w = best_weight_fun(dur,mov_speed);
            ratio = (w.*mov_speed) \ (w.*mov_mu);
            pred  = ratio*mov_speed;
            wmean = sum(w.*mov_mu)/sum(w);
            r2    = 1 - sum(w.*(mov_mu-pred).^2)/sum(w.*(mov_mu-wmean).^2);
            n     = numel(mov_mu);
            adj_r2 = 1 - (1-r2)*(n-1)/(n-1-1); % p=1 (one slope parameter), same correction as figures 5/6

            fly_mob_ratio(end+1) = ratio; %#ok<AGROW>
            fly_adj_r2(end+1)    = adj_r2; %#ok<AGROW>
            cat_x_mob(end+1)     = gd; %#ok<AGROW>
            fly_mob_fly(end+1)   = f; %#ok<AGROW>
        end
    end
end

keep_group = ismember(1:numel(group_defs_all),cat_x_mob)';
kept_idx = find(keep_group);
gd_to_plot_idx = zeros(numel(group_defs_all),1);
gd_to_plot_idx(kept_idx) = 1:numel(kept_idx);

base_colors = lines(numel(geno_order));
all_cat_colors = zeros(numel(group_defs_all),3);
for k = 1:numel(geno_order)
    all_cat_colors(2*k-1,:) = base_colors(k,:);
    all_cat_colors(2*k,:)   = base_colors(k,:);
end

cat_x_mob_plot = gd_to_plot_idx(cat_x_mob);
cat_labels     = {group_defs_all(kept_idx).label};
cat_colors     = all_cat_colors(kept_idx,:);

% figure 4 ALSO requires this fly's own adjusted R^2 to clear min_adj_r2
% -- a separate filter from min_bouts_per_fly (a fly can clear the bout
% count and still have a poorly-fit regression, per figure 8) -- computed
% as its own subset/remap so figure 8 (below) can keep showing the FULL,
% unfiltered distribution of adjusted R^2 across every bout-count-
% qualifying fly, including the ones this filter drops.
min_adj_r2 = 0.6;
r2_ok = fly_adj_r2 > min_adj_r2;

keep_group_r2 = ismember(1:numel(group_defs_all),cat_x_mob(r2_ok))';
kept_idx_r2 = find(keep_group_r2);
gd_to_plot_idx_r2 = zeros(numel(group_defs_all),1);
gd_to_plot_idx_r2(kept_idx_r2) = 1:numel(kept_idx_r2);

cat_x_mob_r2_plot   = gd_to_plot_idx_r2(cat_x_mob(r2_ok));
fly_mob_ratio_r2ok  = fly_mob_ratio(r2_ok);
cat_labels_r2       = {group_defs_all(kept_idx_r2).label};
cat_colors_r2       = all_cat_colors(kept_idx_r2,:);

fprintf('\n=== final pipeline: bump_smooth_frames=%d, heading_smooth=%.2fs, weighting=%s, min_bouts_per_fly=%d, min_adj_r2=%.2f ===\n', ...
    bump_smooth_frames,heading_smooth_opt_s,weight_schemes(best_s).name,min_bouts_per_fly,min_adj_r2);
for gd = 1:numel(group_defs_all)
    n_bouts_ok = sum(cat_x_mob==gd);
    n_r2_ok    = sum(cat_x_mob==gd & r2_ok);
    fprintf('  %-28s %d/%d flies also clear adj R^2 > %.2f\n', group_defs_all(gd).label, n_r2_ok, n_bouts_ok, min_adj_r2);
end

figure(4); clf
set(gcf,'Name','bump mobility (fully optimized pipeline)','Position',[100,100,max(700,120*numel(cat_labels_r2)),600])
groupplot(cat_x_mob_r2_plot,fly_mob_ratio_r2ok,cat_labels_r2,cat_colors_r2)
hold on
plot(xlim,[1,1],':k')
ylabel('bump path length / heading path length (per fly, optimized weighting)')
title(sprintf('bump mobility -- %d-frame bump smoothing, %.2fs heading smoothing, %s-weighted regression, adj R^2 > %.2f', ...
    bump_smooth_frames,heading_smooth_opt_s,weight_schemes(best_s).name,min_adj_r2))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig4_final_group_plot.png'),'Resolution',150)

%% figure: per-fly goodness of fit (adjusted R^2), by genotype and light condition
% same flies/categories as figure 4 (fly_adj_r2 and fly_mob_ratio come out
% of the SAME loop, same min_bouts_per_fly qualification), so this is
% directly comparable to figure 4 point-for-point: a genotype/condition
% whose mobility ratio (figure 4) looks unusual is worth checking here too
% -- a low adjusted R^2 means that ratio rests on a noisy/poorly-fit
% regression, not just an unusual slope.
figure(8); clf
set(gcf,'Name','bump mobility goodness of fit (adjusted R^2)','Position',[100,100,max(700,120*numel(cat_labels)),600])
groupplot(cat_x_mob_plot,fly_adj_r2,cat_labels,cat_colors)
hold on
plot(xlim,[1,1],':k')
ylabel('adjusted R^2 (mov\_mu ~ mov\_speed, per fly, optimized weighting)')
title(sprintf('bump mobility goodness of fit -- %d-frame bump smoothing, %.2fs heading smoothing, %s-weighted regression', ...
    bump_smooth_frames,heading_smooth_opt_s,weight_schemes(best_s).name))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig8_adj_r2_by_group.png'),'Resolution',150)

%% 7) figure: average bump vector strength (rho), by genotype and light condition
% independent QC check, not gated by min_bouts_per_fly or by walking at
% all -- every fly with at least one trial in a genotype x light-condition
% group gets one point here, using chunk_rho_new (the new PVA's own rho)
% averaged across every frame of every trial that fly has in that
% condition. Lets low imaging/bump quality in a genotype or condition be
% seen directly, rather than only inferred indirectly from a low or noisy
% mobility ratio.
fly_mean_rho = [];
cat_x_rho    = [];
for gd = 1:numel(group_defs_all)
    rows = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        fly_mean_rho(end+1) = mean(cat(1,chunk_rho_new{trial_list}),'omitnan'); %#ok<AGROW>
        cat_x_rho(end+1)    = gd; %#ok<AGROW>
    end
end

keep_group_rho = ismember(1:numel(group_defs_all),cat_x_rho)';
kept_idx_rho = find(keep_group_rho);
gd_to_plot_idx_rho = zeros(numel(group_defs_all),1);
gd_to_plot_idx_rho(kept_idx_rho) = 1:numel(kept_idx_rho);

cat_x_rho_plot = gd_to_plot_idx_rho(cat_x_rho);
cat_labels_rho = {group_defs_all(kept_idx_rho).label};
cat_colors_rho = all_cat_colors(kept_idx_rho,:);

figure(7); clf
set(gcf,'Name','average bump vector strength (rho)','Position',[100,100,max(700,120*numel(cat_labels_rho)),600])
groupplot(cat_x_rho_plot,fly_mean_rho,cat_labels_rho,cat_colors_rho)
hold on
yline(bout_rho_thresh,':k')
ylabel('mean rho (bump vector strength), per fly, all frames')
title(sprintf('average bump vector strength (rho), by genotype and light condition (dotted line = rho threshold %.2f)',bout_rho_thresh))
exportgraphics(gcf,fullfile('ugly_figures','exports','v2_fig7_mean_rho_by_group.png'),'Resolution',150)

%% figures: example empty>th closed-loop flies with mobility ratio > example_ratio_thresh
% picks flies from the SAME set figure 4 actually plots (bout-count AND
% adj-R^2 qualified, so these are "real" high-ratio flies, not ones that
% would already have been filtered out as noise) -- sorted by ratio
% descending and capped at example_max_n, so this reports the most
% extreme qualifying examples rather than an arbitrary subset. For each,
% picks that fly's own busiest trial in this condition (most bouts) as
% the one to display, same selection rule as figure 5/6's busiest_i.
example_geno        = 'empty>th';
example_dark        = 0; % closed loop
example_ratio_thresh = 1.5;
example_max_n        = 3;

example_pool = find(strcmp({group_defs_all(cat_x_mob).geno},example_geno) & ...
                     [group_defs_all(cat_x_mob).dark]==example_dark & ...
                     fly_mob_ratio > example_ratio_thresh & r2_ok);
[~,order] = sort(fly_mob_ratio(example_pool),'descend');
example_pool = example_pool(order(1:min(example_max_n,numel(example_pool))));

fprintf('\n=== %d %s (%s) example fly(s) with ratio > %.2f (bout-count and adj-R^2 qualified) ===\n', ...
    numel(example_pool),example_geno,cond_label{example_dark+1},example_ratio_thresh);

for k = 1:numel(example_pool)
    f = fly_mob_fly(example_pool(k));
    ratio_k  = fly_mob_ratio(example_pool(k));
    r2_k     = fly_adj_r2(example_pool(k));
    trial_list = find(fly_num==f & strcmp(genotype,example_geno) & is_dark==example_dark);
    [~,busiest_rel] = max(cellfun(@(c) numel(c),chunk_bout_mu(trial_list)));
    ti = trial_list(busiest_rel);

    fprintf('  %s: ratio=%.2f, adj R^2=%.2f, showing trial %d (%d bouts)\n', ...
        fly_short_label(fly_list{f}),ratio_k,r2_k,ti,numel(chunk_bout_mu{ti}));

    plot_bump_heading_overlay(50+k,all_data(ti),chunk_mu_final{ti},chunk_fz{ti},heading_smooth_opt_s, ...
        sprintf('%s (%s), fly %s, trial %d -- fly ratio=%.2f, adj R^2=%.2f -- bump (white) vs %.2fs-smoothed heading (cyan)', ...
            example_geno,cond_label{example_dark+1},fly_short_label(fly_list{f}),ti,ratio_k,r2_k,heading_smooth_opt_s), ...
        fullfile('ugly_figures','exports',sprintf('v2_fig9_example%d_%s_%s.png', ...
            k,strrep(example_geno,'>','-'),strrep(cond_label{example_dark+1},' ','_'))))
end

%% functions
function [fly_seg,trial_seg] = meta_fly_and_trial_seg(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        trial_seg = parts{end-1};
        fly_seg   = parts{end-2};
    else
        trial_seg = parts{end};
        fly_seg   = parts{end-1};
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
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N.
    parts = strsplit(fly_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function [mu_new,rho_new,f_z] = trial_pva_movmean(trial,smooth_frames)
    % bump position/strength via a plain SAMPLE-COUNT moving average
    % (movmean) of raw fluorescence (im.f), z-scored per glomerulus, then
    % population-vector-averaged against each glomerulus's own angular
    % position (im.alpha) -- same PVA formula as lpsp_rnai_claude.m's
    % trial_smoothed_bump_pva, but movmean(.,smooth_frames) in place of a
    % gaussian/time-converted smooth.
    f_smooth = movmean(trial.im.f,smooth_frames,2);
    f_z = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);

    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp] = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function y = wrap_break(x)
    % unwrap-safe display helper: wraps an unwrapped angle to -pi:pi and
    % breaks the line at the resulting circular jumps (NaN-inserted) so a
    % plot doesn't connect across wraps.
    y = mod(x,2*pi);
    y(y>pi) = y(y>pi) - 2*pi;
    y(find(abs(diff(y))>pi)+1) = nan;
end

function plot_bump_heading_overlay(fig_num,trial,mu,fz,heading_smooth_s,title_str,export_path)
    % imagesc + bump + heading overlay for one trial -- shared by figure 6
    % and any other single-trial diagnostic that wants the same view.
    % alpha unwrapped to the continuous -pi:~3pi range (im.alpha repeats
    % the same -pi:pi sequence once per PB hemisphere), bump/heading each
    % drawn twice (as-is, and shifted +2*pi) so the overlay tracks both
    % hemisphere bands, not just one. mu is used AS-IS (already smoothed
    % upstream, via trial_pva_movmean, on the image timebase). heading =
    % -ft.cue, smoothed at heading_smooth_s on the fictrac timebase, THEN
    % interpolated onto the image timebase -- smoothing must happen
    % before interpolation/wrapping, never after, or it would smear
    % across the circular variable's false -pi/pi discontinuities.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    heading_win = max(1,round(heading_smooth_s/dt));
    heading_smooth_full = movmean(unwrap(-trial.ft.cue(:)),heading_win);

    n_im = numel(mu);
    xb = linspace(xf(1),xf(end),n_im)';
    heading_on_xb = interp1(xf,heading_smooth_full,xb,'linear','extrap');

    alpha_disp = unwrap(trial.im.alpha(:));
    mu_plot      = wrap_break(mu);
    heading_plot = wrap_break(heading_on_xb);

    figure(fig_num); clf
    set(gcf,'Name','PB activity + bump + heading overlay','Position',[100,100,1100,400])
    hold on
    imagesc(1:n_im,alpha_disp,fz)
    set(gca,'YDir','normal','XLim',[1,n_im],'YLim',[alpha_disp(1),alpha_disp(end)])
    plot(1:n_im,mu_plot,'w','LineWidth',1)
    plot(1:n_im,mu_plot+2*pi,'w','LineWidth',1)
    plot(1:n_im,heading_plot,'Color',[0,1,1],'LineWidth',1)
    plot(1:n_im,heading_plot+2*pi,'Color',[0,1,1],'LineWidth',1)
    cb = colorbar; ylabel(cb,'z-score (5-frame moving average)')
    xlabel('imaging frame'); ylabel('PB angle (rad)')
    title(title_str,'Interpreter','none')
    exportgraphics(gcf,export_path,'Resolution',150)
end

function [mov_mu,mov_speed,dur,r_speed_smooth,is_walking] = trial_walking_bouts_v2(trial,mu,rho,heading_smooth_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh)
    % walking-bout detection and per-bout bump/heading path lengths, with
    % ALL smoothing windows given in seconds and converted to samples via
    % THIS trial's own fictrac sample period (dt = median(diff(xf))) --
    % fictrac sampling rate isn't identical across every trial in this
    % dataset, so a shared raw sample count would not represent the same
    % real-time window on every trial. mu (already smoothed upstream, via
    % trial_pva_movmean) is interpolated onto xf and used AS-IS -- no
    % additional bout-level smoothing is applied on top of it, unlike
    % lpsp_rnai_claude.m's pipeline, since this script's bump smoothing is
    % meant to be exactly the one fixed choice from step 2.
    %
    % A bout is only kept if its OWN mean rho (bump vector strength,
    % interpolated the same way as mu) exceeds bout_rho_thresh -- a bout
    % detected purely from the fly's own turning can still coincide with a
    % stretch where the bump signal itself is weak/unreliable (see figure
    % 6's "busiest" example trial, where the bump looks like noise rather
    % than a tracked bump), and such a bout's path length shouldn't be
    % trusted just because the fly was moving.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu);
    xb = linspace(xf(1),xf(end),n_im)';

    mu_f  = interp1(xb,unwrap(mu(:)),xf,'linear','extrap');
    rho_f = interp1(xb,rho(:),xf,'linear','extrap');

    heading_win = max(1,round(heading_smooth_s/dt));
    r_speed_smooth = movmean(trial.ft.r_speed(:),heading_win);
    fly_speed = abs(r_speed_smooth);

    max_gap_frames     = max(1,round(max_gap_s/dt));
    min_walking_frames = max(1,round(min_walking_s/dt));

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    n_bouts = numel(bout_starts);
    mov_mu    = nan(n_bouts,1);
    mov_speed = nan(n_bouts,1);
    dur       = nan(n_bouts,1);
    bout_rho  = nan(n_bouts,1);
    for b = 1:n_bouts
        rng = bout_starts(b):bout_ends(b);
        mov_mu(b)    = sum(abs(diff(mu_f(rng))),'omitnan');
        mov_speed(b) = sum(abs(r_speed_smooth(rng)),'omitnan')*dt;
        dur(b)       = (bout_ends(b)-bout_starts(b)+1)*dt;
        bout_rho(b)  = mean(rho_f(rng),'omitnan');
    end

    keep = bout_rho > bout_rho_thresh;
    mov_mu    = mov_mu(keep);
    mov_speed = mov_speed(keep);
    dur       = dur(keep);
end

function groupplot(cat_x, values, cat_labels, colors)
    % jittered per-point scatter + mean +/- SEM errorbar per category, n
    % embedded in the xtick label itself (see lpsp_rnai_claude.m's own
    % groupplot for why: a floating text() object for "n=" can visually
    % collide with data/neighboring categories once there are more than a
    % handful of categories, and xticklabels() splits any label string
    % containing a literal newline into separate entries rather than
    % rendering 2 lines -- so no embedded '\n' either).
    hold on
    n_cat = numel(cat_labels);
    labels_with_n = cell(1,n_cat);
    for c = 1:n_cat
        y = values(cat_x==c);
        y = y(~isnan(y));
        labels_with_n{c} = sprintf('%s (n=%d)',cat_labels{c},numel(y));
        if isempty(y)
            continue
        end
        jitter = (rand(size(y))-.5)*.3;
        scatter(c+jitter,y,20,colors(c,:),'filled','MarkerFaceAlpha',.3)
        errorbar(c,mean(y),std(y)/sqrt(numel(y)),'o','Color',colors(c,:)*.6, ...
            'MarkerFaceColor',colors(c,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    xticks(1:n_cat); xticklabels(labels_with_n); xtickangle(15)
    xlim([0.5,n_cat+0.5])
end
