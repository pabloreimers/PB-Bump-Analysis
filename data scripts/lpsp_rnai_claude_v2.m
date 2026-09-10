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
data_dir    = fullfile('.data'); % matches lpsp_kir_claude.m -- this repo checkout has ".data", not "data" (confirmed while testing the gain integration in step 8)
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
export_fig('v2_fig1_bump_sanity_check')

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
        % strip NaN BEFORE the mean()/circ_var() calls below, not after --
        % mean() and circ_var() don't ignore NaN by default, so a single
        % dropped-frame NaN anywhere in cue/mu for this fly would otherwise
        % silently poison its entire offset_mean/circ_var result to NaN,
        % indistinguishable from "not enough data" (confirmed directly:
        % this was happening for nearly every fly before this fix, despite
        % pooled_offset routinely holding 10,000+ samples).
        pooled_offset = pooled_offset(~isnan(pooled_offset));
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
export_fig('v2_fig2_heading_smoothing_sweep')

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
export_fig('v2_fig10_mu_smoothing_sweep')

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
export_fig('v2_fig3_weighting_comparison')

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
export_fig('v2_fig5_bout_definition')

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
    'v2_fig6_busiest_trial_overlay')

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
min_adj_r2 = 0.5;
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
export_fig('v2_fig4_final_group_plot')

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
export_fig('v2_fig8_adj_r2_by_group')

%% figure: offset variability (circ_var(circ_dist(mu,heading))), for the SAME flies/trials as figure 4
% restricted to exactly the fly x light-condition entries figure 4 plots
% (r2_ok -- bout-count AND adj-R^2 qualified), not recomputed on a
% different set.
%
% Smoothing AND sample restriction both match this pipeline's actual
% mobility computation, not step 3's separate (and, tried first here,
% much stricter) r_thresh/rho_thresh=.5 instantaneous masking -- that
% combination was tuned for step 3's own sweep on healthy empty-control
% closed-loop flies only, and confirmed directly to leave nearly every
% fly in most OTHER genotype x light-condition groups with fewer than 10
% valid samples (e.g. lpsp>th closed loop dropped from 35 r2_ok flies to
% n=1) -- far too strict once applied dataset-wide. Instead, this reuses
% trial_walking_bouts_v2's own is_walking mask (SAME turn_thresh,
% max_gap_s, min_walking_s, bout_rho_thresh as step 4), so the offset
% variability is computed over EXACTLY the walking-bout samples that
% actually feed the mobility ratio itself. mu = chunk_mu_final (step 2 +
% step 3b, used AS-IS). heading = movmean(unwrap(-ft.cue)) at
% heading_smooth_opt_s -- the SAME window trial_walking_bouts_v2 applies
% to r_speed, just applied here to the heading POSITION instead of its
% derivative (smoothing a signal by window W and differentiating is the
% same real-time-scale operation as smoothing the derivative by W
% directly, for a moving-average kernel). Demeaned via
% atan2(mean(sin),mean(cos)) before circ_var, same as step 3 -- a fixed
% rotational offset between mu's zero and the cue's zero isn't noise,
% only deviation AROUND it is.
fly_offset_var_final = nan(size(fly_mob_ratio));
included = find(r2_ok);
for m = 1:numel(included)
    idx = included(m);
    gd  = cat_x_mob(idx);
    f   = fly_mob_fly(idx);
    trial_list = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark & fly_num==f);

    pooled_offset = [];
    for ti = trial_list'
        trial = all_data(ti);
        xf = trial.ft.xf;
        dt = median(diff(xf));
        heading_win = max(1,round(heading_smooth_opt_s/dt));
        heading_smooth_full = movmean(unwrap(-trial.ft.cue(:)),heading_win);

        n_im = numel(chunk_mu_final{ti});
        xb = linspace(xf(1),xf(end),n_im)';
        mu_i = interp1(xb,unwrap(chunk_mu_final{ti}(:)),xf,'linear','extrap');

        [~,~,~,~,is_walking_ti] = trial_walking_bouts_v2( ...
            trial,chunk_mu_final{ti},chunk_rho_new{ti},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);

        offset = circ_dist(heading_smooth_full,mu_i);
        pooled_offset = [pooled_offset; offset(is_walking_ti)]; %#ok<AGROW>
    end
    % strip NaN BEFORE mean()/circ_var() -- see the identical fix (and the
    % comment explaining why) in step 3's own offset-variability loop above.
    pooled_offset = pooled_offset(~isnan(pooled_offset));
    if numel(pooled_offset) < 10
        continue % not enough valid samples for this fly to trust a circ_var estimate
    end
    offset_mean = atan2(mean(sin(pooled_offset)),mean(cos(pooled_offset)));
    offset_demeaned = circ_dist(pooled_offset,offset_mean);
    fly_offset_var_final(idx) = circ_var(offset_demeaned);
end

figure(13); clf
set(gcf,'Name','offset variability (mu vs heading), figure-4 flies','Position',[100,100,max(700,120*numel(cat_labels_r2)),600])
groupplot(cat_x_mob_r2_plot,fly_offset_var_final(r2_ok),cat_labels_r2,cat_colors_r2)
ylabel('offset variability: circ\_var(circ\_dist(mu,heading))')
title(sprintf('offset variability -- same flies as figure 4 -- mu: %d-frame + step-3b smoothing, heading: %.2fs movmean', ...
    bump_smooth_frames,heading_smooth_opt_s))
export_fig('v2_fig13_offset_variability_by_group')

%% figure: offset variability over a 30s SLIDING WINDOW, mean per fly, same flies/trials as figure 4
% figure 13 pools every walking-bout sample across a fly's whole trial(s)
% into ONE circ_var, demeaned by that fly's OWN OVERALL mean offset -- so
% slow drift in the offset over the course of a trial (a fixed bias
% between mu's zero and the cue's zero can itself wander slowly) gets
% counted as "variability" right alongside genuine moment-to-moment
% tracking noise. This instead slides a 30s window across each trial,
% computes circ_var within each window (demeaned by THAT WINDOW's own
% mean offset, not the whole trial's), and averages the per-window values
% for one fly -- so slow drift across windows no longer inflates the
% result, only within-window jitter does. Same is_walking mask, mu, and
% heading smoothing as figure 13 (and the same NaN-before-mean/circ_var
% fix); windows are skipped if they don't have at least min_window_n
% walking samples of their own to trust a circ_var estimate from.
offset_window_s      = 30;
offset_window_step_s = 5;
min_window_n         = 30;

fly_offset_var_window = nan(size(fly_mob_ratio));
for m = 1:numel(included)
    idx = included(m);
    gd  = cat_x_mob(idx);
    f   = fly_mob_fly(idx);
    trial_list = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark & fly_num==f);

    win_vars = [];
    for ti = trial_list'
        trial = all_data(ti);
        xf = trial.ft.xf;
        dt = median(diff(xf));
        heading_win = max(1,round(heading_smooth_opt_s/dt));
        heading_smooth_full = movmean(unwrap(-trial.ft.cue(:)),heading_win);

        n_im = numel(chunk_mu_final{ti});
        xb = linspace(xf(1),xf(end),n_im)';
        mu_i = interp1(xb,unwrap(chunk_mu_final{ti}(:)),xf,'linear','extrap');

        [~,~,~,~,is_walking_ti] = trial_walking_bouts_v2( ...
            trial,chunk_mu_final{ti},chunk_rho_new{ti},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);

        offset = circ_dist(heading_smooth_full,mu_i);

        window_starts = xf(1):offset_window_step_s:(xf(end)-offset_window_s);
        for w = 1:numel(window_starts)
            in_win = xf >= window_starts(w) & xf < window_starts(w)+offset_window_s;
            win_offset = offset(in_win & is_walking_ti);
            win_offset = win_offset(~isnan(win_offset));
            if numel(win_offset) < min_window_n
                continue
            end
            win_mean = atan2(mean(sin(win_offset)),mean(cos(win_offset)));
            win_vars(end+1) = circ_var(circ_dist(win_offset,win_mean)); %#ok<AGROW>
        end
    end
    if numel(win_vars) < 1
        continue % no window for this fly had enough walking samples to trust
    end
    fly_offset_var_window(idx) = mean(win_vars);
end

figure(14); clf
set(gcf,'Name','offset variability, 30s sliding window','Position',[100,100,max(700,120*numel(cat_labels_r2)),600])
groupplot(cat_x_mob_r2_plot,fly_offset_var_window(r2_ok),cat_labels_r2,cat_colors_r2)
ylabel('mean offset variability across 30s windows: mean(circ\_var(circ\_dist(mu,heading)))')
title(sprintf('offset variability, %ds sliding window (step %ds) -- same flies as figure 4 -- mu: %d-frame + step-3b smoothing, heading: %.2fs movmean', ...
    offset_window_s,offset_window_step_s,bump_smooth_frames,heading_smooth_opt_s))
export_fig('v2_fig14_offset_variability_sliding_window')

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
export_fig('v2_fig7_mean_rho_by_group')

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
        sprintf('v2_fig9_example%d_%s_%s', ...
            k,strrep(example_geno,'>','-'),strrep(cond_label{example_dark+1},' ','_')))
end

%% figures: every closed-loop trial, one figure per genotype (lpsp>th, empty>th), one subplot per trial
% every trial gets its own subplot (not one per fly -- a fly with 2
% closed-loop trials gets 2 subplots here), using the same
% plot_bump_heading_overlay_ax drawing as figures 6/9, just laid out in a
% grid instead of one trial per standalone figure. No per-subplot
% colorbar (one shared colorbar for the whole figure instead) and a
% terser per-subplot title (just fly + trial), since there isn't room for
% the full single-trial title at this scale.
grid_genos = {'lpsp>th','empty>th'};

for g = 1:numel(grid_genos)
    geno_g = grid_genos{g};
    trial_list = find(strcmp(genotype,geno_g) & ~is_dark);

    n_t = numel(trial_list);
    ncols = ceil(sqrt(n_t));
    nrows = ceil(n_t/ncols);

    figure(20+g); clf
    set(gcf,'Name',sprintf('%s: every closed-loop trial',geno_g),'Position',[50,50,min(1800,280*ncols),min(1000,200*nrows)])
    t = tiledlayout(nrows,ncols,'TileSpacing','compact','Padding','compact');
    for k = 1:n_t
        ti = trial_list(k);
        nexttile(t);
        plot_bump_heading_overlay_ax(all_data(ti),chunk_mu_final{ti},chunk_fz{ti},heading_smooth_opt_s, ...
            sprintf('%s, trial %d',fly_short_label(fly_list{fly_num(ti)}),ti))
        xticks([]); yticks([]); xlabel(''); ylabel('')
    end
    cb = colorbar; cb.Layout.Tile = 'east'; ylabel(cb,'z-score')
    title(t,sprintf('%s, closed loop -- every trial (n=%d) -- bump (white) vs %.2fs-smoothed heading (cyan)',geno_g,n_t,heading_smooth_opt_s))
    export_fig(sprintf('v2_fig12_%s_all_closed_loop_trials',strrep(geno_g,'>','-')))
end

%% 8) velocity gain (bump vs. visual cue, bump vs. fly's own rotation), every genotype x light-condition group
% A SEPARATE gain metric from this script's own bump-mobility pipeline
% above (steps 1-7): gain_cue/gain_fly are regression slopes between
% bump velocity and cue/fly-rotation velocity, answering "does the bump
% move at the right RATE relative to what it's tracking", rather than
% mobility's "how much total path length does the bump cover relative to
% heading, over a whole walking bout". Developed and parameter-swept
% separately in gain_scratch_claude.m (not reproduced here) -- ported in
% as the SETTLED result of that search, reusing this script's own
% trial_pva_movmean/export_fig/groupplot rather than duplicating them.
%
% Why the smoothing constants below differ from this script's own
% heading_smooth_opt_s/mu_smooth_opt_s (steps 2-3): those were tuned for
% the mobility-ratio metric (integrated path length over whole bouts,
% forgiving of frame-level velocity noise); gain_cue/gain_fly are
% instantaneous velocity regressions, far more sensitive to per-frame
% noise, and were sequentially/greedily swept on their OWN criterion
% (per-fly combined MSE from gain_cue's target of 1 and gain_fly's target
% of 0.8, pooled over empty-control flies, closed loop only) in
% gain_scratch_claude.m's own section 10. bump/cue/fly_vel smoothing use
% smoothdata(...,'gaussian',...), NOT this script's movmean, and cue/fly
% rotation get INDEPENDENT windows (no reason a visual-scene signal and a
% fictrac ball-rotation signal need the same amount of smoothing).
%
% gain_im_f_frames intentionally does NOT reuse bump_smooth_frames (this
% script's own step-2 constant, 5 frames) -- the gain sweep found its own
% winner (3 frames) independently and there's no reason the two should
% coincide.
gain_im_f_frames    = 3;    % frames, movmean on im.f before z-score+PVA (via trial_pva_movmean, same function as this script's own step 2)
gain_bump_smooth_s  = 0.20; % s, single Gaussian pass on unwrapped mu (image timebase)
gain_cue_smooth_s   = 0.75; % s, single Gaussian pass on unwrapped cue (fictrac timebase)
gain_vel_smooth_s   = 1.00; % s, single Gaussian pass on r_speed directly (fictrac timebase)
gain_lag_frames     = 8;    % fictrac frames (~0.13s): cue_vel/fly_vel led bump_vel by this much on the trial used to find it (gain_scratch_claude.m)
gain_vel_thresh     = 0.2; gain_bump_thresh = 10; gain_rho_thresh = 0.2; gain_vel_max = 5; % rad/s, inherited starting point from lpsp_kir_claude.m's own gain-scatter section

% light-condition-major group ordering (every closed-loop group, then
% every dark group) rather than this script's own genotype-major
% group_defs_all (step 5) -- easier to compare across genotypes within
% one light condition at a glance. Dark groups are included deliberately:
% gain_cue should drop toward 0 there (no visual pattern for the panel to
% show, so cue_vel carries little real signal), a useful built-in
% negative control rather than something to exclude.
gain_group_defs = struct('geno',{},'dark',{},'label',{});
for ci = 1:numel(cond_label)
    for gi = 1:numel(geno_order)
        gain_group_defs(end+1) = struct('geno',geno_order{gi},'dark',ci-1,'label',sprintf('%s (%s)',geno_order{gi},cond_label{ci})); %#ok<SAGROW>
    end
end

gain_fly_list   = cell(1,numel(gain_group_defs));
gain_fly_trials = cell(1,numel(gain_group_defs));
for gd = 1:numel(gain_group_defs)
    rows = find(strcmp(genotype,gain_group_defs(gd).geno) & is_dark==gain_group_defs(gd).dark);
    these_flies = unique(fly_num(rows));
    gain_fly_list{gd}   = these_flies;
    gain_fly_trials{gd} = arrayfun(@(f) rows(fly_num(rows)==f), these_flies, 'UniformOutput',false);
end

% PVA cache at gain_im_f_frames, computed once per trial that's actually
% used by any group above (every trial belonging to one of geno_order,
% both light conditions) -- trial_pva_movmean is this script's own step-2
% function, called here with a DIFFERENT frame count than chunk_mu_new
% (which used bump_smooth_frames=5) since the gain sweep found its own
% winner independently.
gain_pva_mu  = cell(n_trials,1);
gain_pva_rho = cell(n_trials,1);
gain_all_trials = find(ismember(genotype,geno_order));
for ti = gain_all_trials'
    [gain_pva_mu{ti},gain_pva_rho{ti}] = trial_pva_movmean(all_data(ti),gain_im_f_frames);
end

base_colors_gain = lines(numel(geno_order));
cat_colors_gain  = zeros(numel(gain_group_defs),3);
for gd = 1:numel(gain_group_defs)
    gi = find(strcmp(geno_order,gain_group_defs(gd).geno));
    cat_colors_gain(gd,:) = base_colors_gain(gi,:);
end
cat_labels_gain = {gain_group_defs.label};

fly_gain_cue_bygroup  = []; fly_gain_fly_bygroup  = []; cat_x_gain  = [];
fly_gain_cue_noint    = []; fly_gain_fly_noint    = []; cat_x_gain_noint = [];
fprintf('\n=== velocity gain (bump vs. cue, bump vs. fly rotation), every genotype x light-condition group ===\n');
for gd = 1:numel(gain_group_defs)
    these_flies = gain_fly_list{gd};
    trial_lists = gain_fly_trials{gd};
    gc  = nan(numel(these_flies),1); gf  = nan(numel(these_flies),1); % affine (with intercept)
    gc0 = nan(numel(these_flies),1); gf0 = nan(numel(these_flies),1); % through-origin
    for ff = 1:numel(these_flies)
        [gc(ff),gf(ff)]   = fly_gain_v3(all_data,trial_lists{ff},gain_pva_mu,gain_pva_rho, ...
            gain_bump_smooth_s,gain_cue_smooth_s,gain_vel_smooth_s,gain_lag_frames, ...
            gain_vel_thresh,gain_bump_thresh,gain_rho_thresh,gain_vel_max,true);
        [gc0(ff),gf0(ff)] = fly_gain_v3(all_data,trial_lists{ff},gain_pva_mu,gain_pva_rho, ...
            gain_bump_smooth_s,gain_cue_smooth_s,gain_vel_smooth_s,gain_lag_frames, ...
            gain_vel_thresh,gain_bump_thresh,gain_rho_thresh,gain_vel_max,false);
    end
    fprintf('  %-24s gain_cue = %.3f +/- %.3f, gain_fly = %.3f +/- %.3f (n=%d flies)\n', ...
        gain_group_defs(gd).label, mean(gc,'omitnan'), std(gc,'omitnan')/sqrt(sum(~isnan(gc))), ...
        mean(gf,'omitnan'), std(gf,'omitnan')/sqrt(sum(~isnan(gf))), sum(~isnan(gc)));

    fly_gain_cue_bygroup = [fly_gain_cue_bygroup; gc]; %#ok<AGROW>
    fly_gain_fly_bygroup = [fly_gain_fly_bygroup; gf]; %#ok<AGROW>
    cat_x_gain           = [cat_x_gain; gd*ones(numel(these_flies),1)]; %#ok<AGROW>
    fly_gain_cue_noint    = [fly_gain_cue_noint; gc0]; %#ok<AGROW>
    fly_gain_fly_noint    = [fly_gain_fly_noint; gf0]; %#ok<AGROW>
    cat_x_gain_noint      = [cat_x_gain_noint; gd*ones(numel(these_flies),1)]; %#ok<AGROW>
end

figure(30); clf
set(gcf,'Name','velocity gain, every genotype x light-condition group','Position',[100,100,1200,800])
subplot(2,1,1)
groupplot(cat_x_gain,fly_gain_cue_bygroup,cat_labels_gain,cat_colors_gain)
hold on; plot(xlim,[1,1],':k'); plot(xlim,[0,0],'-','Color',[.85,.85,.85])
ylabel('gain\_cue = slope(bump\_vel ~ 1 + cue\_vel)')
title('bump vs. visual cue (target=1 in closed loop, ~0 expected in dark)')
subplot(2,1,2)
groupplot(cat_x_gain,fly_gain_fly_bygroup,cat_labels_gain,cat_colors_gain)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ 1 + fly\_vel)')
title('bump vs. fly rotation (target=0.8 in closed loop)')
sgtitle(sprintf('velocity gain: im.f=%d frames, bump=%.2fs, cue=%.2fs, fly\\_vel=%.2fs, lag=%d frames', ...
    gain_im_f_frames,gain_bump_smooth_s,gain_cue_smooth_s,gain_vel_smooth_s,gain_lag_frames),'Interpreter','none')
export_fig('v2_fig30_velocity_gain_by_group')

figure(31); clf
set(gcf,'Name','velocity gain (through-origin fit), every genotype x light-condition group','Position',[100,100,1200,800])
subplot(2,1,1)
groupplot(cat_x_gain_noint,fly_gain_cue_noint,cat_labels_gain,cat_colors_gain)
hold on; plot(xlim,[1,1],':k'); plot(xlim,[0,0],'-','Color',[.85,.85,.85])
ylabel('gain\_cue = slope(bump\_vel ~ cue\_vel)')
title('bump vs. visual cue (target=1 in closed loop, ~0 expected in dark)')
subplot(2,1,2)
groupplot(cat_x_gain_noint,fly_gain_fly_noint,cat_labels_gain,cat_colors_gain)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ fly\_vel)')
title('bump vs. fly rotation (target=0.8 in closed loop)')
sgtitle(sprintf('THROUGH-ORIGIN velocity gain: im.f=%d frames, bump=%.2fs, cue=%.2fs, fly\\_vel=%.2fs, lag=%d frames', ...
    gain_im_f_frames,gain_bump_smooth_s,gain_cue_smooth_s,gain_vel_smooth_s,gain_lag_frames),'Interpreter','none')
export_fig('v2_fig31_velocity_gain_by_group_through_origin')

%% 9) bump amplitude vs. rotational speed: individual flies (faint) + group mean +/- SEM (thick), by genotype
% same "does bump amplitude scale with rotational speed" analysis as
% lpsp_kir_claude.m / lpsp_tnt_claude.m: amplitude is the "actual" peak
% z-score across glomeruli each frame (not a model fit), |rotational
% speed| is ft.r_speed, and amplitude is shifted `gain_lag_frames` (step
% 8's already-established "behavior leads bump" lag, ~0.13s) later than
% speed -- reused here rather than running a separate amplitude-specific
% lag sweep, since that wasn't asked for.
%
% amplitude is built from chunk_fz (step 2's 5-frame-movmean z-scored
% fluorescence, the same signal this script's own bump-position PVA is
% fit from), NOT a fresh smoothing choice. im.f/chunk_fz is
% hemisphere-resolved -- confirmed directly (im.alpha's two 16-element
% halves are bit-for-bit identical) -- so the two hemispheres are averaged
% together first in trial_speed_amp_v2, matching the "actual peak (max of
% 16 wedges)" convention lpsp_compartments_claude_script.m used for this
% same hemisphere-resolved shape, before taking the peak across wedges
% (same approach as lpsp_tnt_claude.m's own trial_speed_amp).
%
% TWO genotype groupings, per explicit request, since mcherry is the
% shared lpsp-driver-only control for both RNAi lines: {empty>th, lpsp>th,
% lpsp>mcherry} and {empty>vglut, lpsp>vglut, lpsp>mcherry}. Colors are
% consistent across both plots: empty>X (wild-type-like control) = black,
% lpsp>X (the RNAi being tested) = red, lpsp>mcherry (driver-only control,
% shared by both groupings) = blue. THREE figures per grouping (closed
% loop only, dark only, closed loop + dark pooled per fly), same
% three-figure layout as lpsp_kir_claude.m/lpsp_tnt_claude.m -- six
% figures total.
amp_speed_edges = 0:0.2:3; % rad/s, matching lpsp_compartments_claude_script.m
amp_speed_x     = amp_speed_edges(1:end-1) + diff(amp_speed_edges)/2;
min_bin_n_fly   = 10; % minimum per-fly pooled samples required to trust a bin

chunk_speed_amp = cell(n_trials,1);
chunk_peak_amp  = cell(n_trials,1);
for i = 1:n_trials
    [chunk_speed_amp{i},chunk_peak_amp{i}] = trial_speed_amp_v2(all_data(i),chunk_fz{i},gain_lag_frames);
end

amp_groupings = struct('name',{'th','vglut'}, ...
    'genos',{{'empty>th','lpsp>th','lpsp>mcherry'},{'empty>vglut','lpsp>vglut','lpsp>mcherry'}});
amp_colors = [0,0,0; 1,0,0; 0,0.4470,0.7410]; % empty>X = black, lpsp>X = red, lpsp>mcherry = blue

for gi = 1:numel(amp_groupings)
    genos_g = amp_groupings(gi).genos;
    n_g = numel(genos_g);

    binned_cl   = cell(1,n_g);
    binned_dark = cell(1,n_g);
    binned_all  = cell(1,n_g);
    for k = 1:n_g
        binned_cl{k}   = fly_binned_amp(genos_g{k},0,genotype,is_dark,fly_num,chunk_speed_amp,chunk_peak_amp,amp_speed_edges,amp_speed_x,min_bin_n_fly);
        binned_dark{k} = fly_binned_amp(genos_g{k},1,genotype,is_dark,fly_num,chunk_speed_amp,chunk_peak_amp,amp_speed_edges,amp_speed_x,min_bin_n_fly);
        binned_all{k}  = fly_binned_amp(genos_g{k},[],genotype,is_dark,fly_num,chunk_speed_amp,chunk_peak_amp,amp_speed_edges,amp_speed_x,min_bin_n_fly); % [] = both light conditions pooled per fly
    end

    fig_base = 40 + 3*(gi-1);

    figure(fig_base); clf
    set(gcf,'Name',sprintf('bump amplitude vs speed: %s group, closed loop only',amp_groupings(gi).name),'Position',[100,100,700,550])
    plot_amp_by_genotype(binned_cl,genos_g,amp_colors,amp_speed_x)
    title(sprintf('bump amplitude vs. rotational speed -- %s genotypes, closed loop only',amp_groupings(gi).name),'Interpreter','none')
    export_fig(sprintf('v2_fig%d_amp_speed_%s_closedloop',fig_base,amp_groupings(gi).name))

    figure(fig_base+1); clf
    set(gcf,'Name',sprintf('bump amplitude vs speed: %s group, dark only',amp_groupings(gi).name),'Position',[100,100,700,550])
    plot_amp_by_genotype(binned_dark,genos_g,amp_colors,amp_speed_x)
    title(sprintf('bump amplitude vs. rotational speed -- %s genotypes, dark only',amp_groupings(gi).name),'Interpreter','none')
    export_fig(sprintf('v2_fig%d_amp_speed_%s_dark',fig_base+1,amp_groupings(gi).name))

    figure(fig_base+2); clf
    set(gcf,'Name',sprintf('bump amplitude vs speed: %s group, CL+dark pooled',amp_groupings(gi).name),'Position',[100,100,700,550])
    plot_amp_by_genotype(binned_all,genos_g,amp_colors,amp_speed_x)
    title(sprintf('bump amplitude vs. rotational speed -- %s genotypes, closed loop + dark pooled within each fly',amp_groupings(gi).name),'Interpreter','none')
    export_fig(sprintf('v2_fig%d_amp_speed_%s_combined',fig_base+2,amp_groupings(gi).name))
end

%% 10) bump occupancy: does the bump cover as much of the circle as heading does, or does it cluster at a few preferred positions?
% same "bump occupancy / stickiness" analysis as lpsp_kir_claude.m /
% lpsp_tnt_claude.m, ported onto this script's own optimized pipeline
% rather than copied verbatim:
%   mu      = chunk_mu_final (step 2's PVA + step 3b's additional smoothing)
%   heading = movmean(unwrap(-ft.cue)) at heading_smooth_opt_s (step 3's optimum)
%   samples restricted to trial_walking_bouts_v2's own is_walking mask --
%     the SAME "discernable bump" convention already used by this script's
%     own offset-variability figures (13/14), NOT step 3's stricter
%     instantaneous r_thresh/rho_thresh=.5 masking, which that step's own
%     header comment already found "far too strict" once applied
%     dataset-wide (e.g. lpsp>th closed loop dropped to n=1 fly).
%   NO lag shift between mu and heading -- matching figures 13/14's own
%     convention for this exact comparison (unlike lpsp_kir_claude.m/
%     lpsp_tnt_claude.m, which shift bump signals by a velocity-gain-
%     derived lag): heading_smooth_opt_s's 1.5s window was itself chosen
%     in step 3 specifically to minimize this same offset's variability,
%     so adding a separate lag on top would be re-optimizing with an
%     extra free parameter this script never swept for.
%
% grouped by ALL 10 genotype x light-condition categories (group_defs_all,
% the same convention as this script's own figures 4/7/8/13/14), not the
% 2-genotype-plus-mcherry line-plot grouping used for the amplitude-vs-
% speed figures above -- a groupplot's categorical x-axis handles 10
% categories fine, unlike overlaying 10 sets of individual-fly lines on
% one plot.
%
% three comparisons, all per fly:
%   1) stickiness = 1 - circ_var(mu)/circ_var(heading)
%   2) circular concentration (1-circ_var), bump vs. heading side by side
%   3) circ_var(mu) - circ_var(heading) -- the same comparison as a plain
%      subtraction instead of a ratio (explicitly requested alongside kir/tnt)
chunk_mu_ok  = cell(n_trials,1);
chunk_cue_ok = cell(n_trials,1);
for i = 1:n_trials
    [chunk_mu_ok{i},chunk_cue_ok{i}] = trial_bumpok_samples_v2( ...
        all_data(i),chunk_mu_final{i},chunk_rho_new{i},heading_smooth_opt_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);
end

% entropy-based counterparts to the circ_var-based quantities below (see
% the header comment on the entropy figures, further down, for why
% circ_var alone can be misleading here): H = -sum(p.*log(p)) of each
% fly's own histogram over occ_edges, via the circ_entropy helper. Unlike
% circ_var, entropy only depends on how mass is spread across bins, not
% their angular position, so it isn't fooled by a symmetric multimodal
% (e.g. antipodal two-peaked) distribution the way a resultant-vector-
% based statistic can be. occ_edges wasn't previously needed in this
% script (unlike lpsp_kir_claude.m/lpsp_tnt_claude.m, this step never
% built a per-fly histogram-overlay figure), so it's defined here for the
% first time, matching those scripts' own 32-bin PB-wedge resolution.
occ_edges   = -pi:pi/16:pi; % 32 bins, matching the PB's own wedge resolution
n_occ_bins  = numel(occ_edges)-1;
max_entropy = log(n_occ_bins); % entropy of a perfectly uniform distribution over n_occ_bins bins -- normalizer for the "concentration" figures below

min_samples_stick = 100;
fly_stickiness   = [];
fly_mu_conc      = [];
fly_heading_conc = [];
fly_var_diff     = [];
fly_stickiness_ent   = []; % entropy-based stickiness: 1 - H(mu)/H(heading)
fly_mu_conc_ent      = []; % entropy-based concentration of mu: 1 - H(mu)/max_entropy
fly_heading_conc_ent = []; % entropy-based concentration of heading, for reference
fly_var_diff_ent     = []; % H(mu) - H(heading)
cat_x_stick      = [];
for gd = 1:numel(group_defs_all)
    rows = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark);
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mu_f  = cat(1,chunk_mu_ok{trial_list});
        cue_f = cat(1,chunk_cue_ok{trial_list});
        if numel(mu_f) < min_samples_stick
            continue
        end
        cue_var = circ_var(cue_f);
        if cue_var < 1e-3
            % this fly's heading barely varied at all within its own
            % walking-bout-restricted samples (near-delta-function
            % clustering) -- dividing by ~0 blows the ratio up to +/-Inf
            % (confirmed directly: this happened for a real empty>th dark
            % fly, poisoning that group's entire mean to -Inf), so the
            % ratio is skipped (left NaN, same as groupplot already
            % excludes) rather than trusted here. fly_mu_conc/
            % fly_heading_conc/fly_var_diff don't divide by cue_var, so
            % they're unaffected and still computed.
            fly_stickiness(end+1) = nan; %#ok<AGROW>
        else
            fly_stickiness(end+1) = 1 - circ_var(mu_f)/cue_var; %#ok<AGROW>
        end
        fly_mu_conc(end+1)      = 1 - circ_var(mu_f); %#ok<AGROW>
        fly_heading_conc(end+1) = 1 - cue_var; %#ok<AGROW>
        fly_var_diff(end+1)     = circ_var(mu_f) - cue_var; %#ok<AGROW>

        % unlike lpsp_kir_claude.m/lpsp_tnt_claude.m (whose cue_f is
        % already confined to [-pi,pi], coming straight from -ft.cue with
        % no unwrap), this script's own trial_bumpok_samples_v2 unwraps
        % BOTH mu and heading before smoothing (needed to smooth correctly
        % across wrap boundaries), so both mu_f and cue_f here are
        % cumulative/unbounded -- histcounts with fixed occ_edges doesn't
        % wrap modularly, so both need to be wrapped back into [-pi,pi]
        % before binning, unlike the circ_var-based quantities above
        % (which don't care about wrapping either way).
        mu_f_wrapped  = mod(mu_f+pi,2*pi)-pi;
        cue_f_wrapped = mod(cue_f+pi,2*pi)-pi;
        Hm = circ_entropy(mu_f_wrapped,occ_edges);
        Hh = circ_entropy(cue_f_wrapped,occ_edges);
        if Hh < 1e-3
            % same near-zero-denominator guard as the circ_var ratio above
            fly_stickiness_ent(end+1) = nan; %#ok<AGROW>
        else
            fly_stickiness_ent(end+1) = 1 - Hm/Hh; %#ok<AGROW>
        end
        fly_mu_conc_ent(end+1)      = 1 - Hm/max_entropy; %#ok<AGROW>
        fly_heading_conc_ent(end+1) = 1 - Hh/max_entropy; %#ok<AGROW>
        fly_var_diff_ent(end+1)     = Hm - Hh; %#ok<AGROW>

        cat_x_stick(end+1)      = gd; %#ok<AGROW>
    end
end

keep_group_stick = ismember(1:numel(group_defs_all),cat_x_stick)';
kept_idx_stick = find(keep_group_stick);
gd_to_plot_idx_stick = zeros(numel(group_defs_all),1);
gd_to_plot_idx_stick(kept_idx_stick) = 1:numel(kept_idx_stick);

cat_x_stick_plot = gd_to_plot_idx_stick(cat_x_stick);
cat_labels_stick = {group_defs_all(kept_idx_stick).label};
cat_colors_stick = all_cat_colors(kept_idx_stick,:);

fprintf('\n=== bump occupancy / stickiness index (1 - circ.var(mu)/circ.var(heading)), per group ===\n');
for gd = 1:numel(group_defs_all)
    y = fly_stickiness(cat_x_stick==gd);
    y = y(~isnan(y)); % same NaN exclusion as groupplot's own embedded n= label, so the count printed here matches what's actually plotted
    if isempty(y); continue; end
    fprintf('  %-28s n=%2d flies   mean=%.3f\n', group_defs_all(gd).label, numel(y), mean(y));
end

figure(46); clf
set(gcf,'Name','bump occupancy: stickiness index per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
groupplot(cat_x_stick_plot,fly_stickiness,cat_labels_stick,cat_colors_stick)
ylabel('stickiness = 1 - circ\_var(mu) / circ\_var(heading)')
title('bump occupancy: how much less of the circle the bump covers vs. heading, one point per fly')
export_fig('v2_fig46_occupancy_stickiness')

%% figure: circular occupancy -- bump concentration vs. heading concentration, per fly
figure(47); clf
set(gcf,'Name','bump occupancy: circular concentration per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
hold on
gray = [.6,.6,.6];
labels_conc_n = cell(1,numel(cat_labels_stick)); % n embedded in xtick label, same rationale as this script's own groupplot
for cIdx = 1:numel(cat_labels_stick)
    y_mu  = fly_mu_conc(cat_x_stick_plot==cIdx);
    y_hdg = fly_heading_conc(cat_x_stick_plot==cIdx);
    labels_conc_n{cIdx} = sprintf('%s (n=%d)',cat_labels_stick{cIdx},numel(y_mu));
    if ~isempty(y_mu)
        jit = (rand(size(y_mu))-.5)*.25;
        scatter(cIdx-0.18+jit,y_mu,20,cat_colors_stick(cIdx,:),'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx-0.18,mean(y_mu),std(y_mu)/sqrt(numel(y_mu)),'o','Color',cat_colors_stick(cIdx,:)*.6, ...
            'MarkerFaceColor',cat_colors_stick(cIdx,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    if ~isempty(y_hdg)
        jit = (rand(size(y_hdg))-.5)*.25;
        scatter(cIdx+0.18+jit,y_hdg,20,gray,'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx+0.18,mean(y_hdg),std(y_hdg)/sqrt(numel(y_hdg)),'o','Color',gray*.6, ...
            'MarkerFaceColor',gray*.6,'LineWidth',2,'MarkerSize',7)
    end
end
xticks(1:numel(cat_labels_stick)); xticklabels(labels_conc_n); xtickangle(15)
xlim([0.5,numel(cat_labels_stick)+0.5])
h_mu  = scatter(nan,nan,20,[0,0,0],'filled');
h_hdg = scatter(nan,nan,20,gray,'filled');
legend([h_mu,h_hdg],{'bump concentration','heading concentration'},'Location','eastoutside')
ylabel('circular concentration, 1-circ\_var (0=uniform, 1=one point)')
title('bump occupancy: circular concentration (color) vs. heading concentration alone (gray), one point per fly')
export_fig('v2_fig47_occupancy_concentration')

%% figure: circ_var(mu) - circ_var(heading), per fly
% same comparison as figure 46's ratio, but a plain subtraction: positive
% means the bump varies MORE (is LESS concentrated/stuck) than heading --
% the opposite direction from "stickiness", which is framed so higher = more stuck.
figure(48); clf
set(gcf,'Name','bump occupancy: circ_var(mu) - circ_var(heading) per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
groupplot(cat_x_stick_plot,fly_var_diff,cat_labels_stick,cat_colors_stick)
ylabel('circ\_var(mu) - circ\_var(heading)')
title('bump occupancy: circ\_var(mu) - circ\_var(heading), by genotype and light condition (one point per fly)')
export_fig('v2_fig48_occupancy_var_diff')

%% figure: bump occupancy, entropy-based stickiness (same comparison as figure 46, using Shannon entropy instead of circ_var)
figure(49); clf
set(gcf,'Name','bump occupancy: entropy stickiness index per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
groupplot(cat_x_stick_plot,fly_stickiness_ent,cat_labels_stick,cat_colors_stick)
ylabel('entropy stickiness = 1 - H(mu) / H(heading)')
title('bump occupancy (entropy-based): how much less of the circle the bump covers vs. heading, one point per fly')
export_fig('v2_fig49_occupancy_entropy_stickiness')

%% figure: circular occupancy -- entropy-based concentration, bump vs. heading, per fly
figure(50); clf
set(gcf,'Name','bump occupancy: entropy concentration per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
hold on
gray = [.6,.6,.6];
labels_conc_ent_n = cell(1,numel(cat_labels_stick));
for cIdx = 1:numel(cat_labels_stick)
    y_mu  = fly_mu_conc_ent(cat_x_stick_plot==cIdx);
    y_hdg = fly_heading_conc_ent(cat_x_stick_plot==cIdx);
    labels_conc_ent_n{cIdx} = sprintf('%s (n=%d)',cat_labels_stick{cIdx},numel(y_mu));
    if ~isempty(y_mu)
        jit = (rand(size(y_mu))-.5)*.25;
        scatter(cIdx-0.18+jit,y_mu,20,cat_colors_stick(cIdx,:),'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx-0.18,mean(y_mu),std(y_mu)/sqrt(numel(y_mu)),'o','Color',cat_colors_stick(cIdx,:)*.6, ...
            'MarkerFaceColor',cat_colors_stick(cIdx,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    if ~isempty(y_hdg)
        jit = (rand(size(y_hdg))-.5)*.25;
        scatter(cIdx+0.18+jit,y_hdg,20,gray,'filled','MarkerFaceAlpha',.3)
        errorbar(cIdx+0.18,mean(y_hdg),std(y_hdg)/sqrt(numel(y_hdg)),'o','Color',gray*.6, ...
            'MarkerFaceColor',gray*.6,'LineWidth',2,'MarkerSize',7)
    end
end
xticks(1:numel(cat_labels_stick)); xticklabels(labels_conc_ent_n); xtickangle(15)
xlim([0.5,numel(cat_labels_stick)+0.5])
h_mu  = scatter(nan,nan,20,[0,0,0],'filled');
h_hdg = scatter(nan,nan,20,gray,'filled');
legend([h_mu,h_hdg],{'bump concentration','heading concentration'},'Location','eastoutside')
ylabel('entropy concentration, 1 - H/max entropy (0=uniform, 1=one bin)')
title('bump occupancy (entropy-based): concentration (color) vs. heading concentration alone (gray), one point per fly')
export_fig('v2_fig50_occupancy_entropy_concentration')

%% figure: H(mu) - H(heading), per fly
% same comparison as figure 49's ratio, but a plain subtraction: positive
% means the bump is MORE spread out (LESS concentrated) than heading --
% the opposite direction from "stickiness", which is framed so higher = more stuck.
figure(51); clf
set(gcf,'Name','bump occupancy: H(mu) - H(heading) per fly','Position',[100,100,max(700,120*numel(cat_labels_stick)),600])
groupplot(cat_x_stick_plot,fly_var_diff_ent,cat_labels_stick,cat_colors_stick)
ylabel('H(mu) - H(heading) (nats)')
title('bump occupancy (entropy-based): H(mu) - H(heading), by genotype and light condition (one point per fly)')
export_fig('v2_fig51_occupancy_entropy_H_diff')

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

function x_filled = fill_nan_gaps_pi(x)
    % ft.cue's scattered NaN dropouts happen almost exclusively at the
    % wrapped +/-pi boundary (confirmed in gain_scratch_claude.m: 98.1% of
    % gaps have a bracketing value within 0.3 rad of +/-pi, across every
    % trial in this dataset) -- linear interpolation on the raw wrapped
    % signal can walk the WRONG way around the circle across such a gap
    % (a spurious jump of up to 2*pi, confirmed directly: 26.7% of all
    % gaps produced a >1 rad linear jump despite a true circular distance
    % <0.5 rad), so dropped samples are set to exactly pi -- landing them
    % at the same boundary their neighbors already sit at, so unwrap()'s
    % own +/-2*pi correction lines them up correctly -- rather than
    % interpolated.
    x_filled = x(:);
    x_filled(isnan(x_filled)) = pi;
end

function [fly_vel,cue_vel,bump_vel,valid] = trial_gain_vectors_v3(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % all three vectors on the fictrac timebase (ft.xf), GAUSSIAN
    % smoothing (not this script's own movmean) with cue_smooth_s and
    % vel_smooth_s as INDEPENDENT windows -- see step 8's header comment.
    % lag_frames shifts cue_vel/fly_vel earlier relative to bump_vel (the
    % bump/calcium signal lags behind behavior), same shift-and-trim
    % convention as lpsp_kir_claude.m's own trial_gain_vectors.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu_raw);
    xb = linspace(xf(1),xf(end),n_im)';
    dt_im = (xf(end)-xf(1)) / (n_im-1);

    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    bump_vel_full = gradient(interp1(xb,mu_smoothed,xf,'linear','extrap'))/dt;
    rho_full = interp1(xb,rho_raw(:),xf,'linear','extrap');

    win_vel = max(1,round(vel_smooth_s/dt));
    fly_vel_full = smoothdata(trial.ft.r_speed(:),'gaussian',win_vel);

    win_cue = max(1,round(cue_smooth_s/dt));
    cue_filled = fill_nan_gaps_pi(trial.ft.cue);
    cue_smoothed = smoothdata(unwrap(-cue_filled),'gaussian',win_cue);
    cue_vel_full = gradient(cue_smoothed)/dt;

    if lag_frames == 0
        fly_vel = fly_vel_full; cue_vel = cue_vel_full; bump_vel = bump_vel_full; rho_i = rho_full;
    elseif lag_frames > 0
        fly_vel  = fly_vel_full(1:end-lag_frames);
        cue_vel  = cue_vel_full(1:end-lag_frames);
        bump_vel = bump_vel_full(lag_frames+1:end);
        rho_i    = rho_full(lag_frames+1:end);
    else
        fly_vel  = fly_vel_full(-lag_frames+1:end);
        cue_vel  = cue_vel_full(-lag_frames+1:end);
        bump_vel = bump_vel_full(1:end+lag_frames);
        rho_i    = rho_full(1:end+lag_frames);
    end

    valid = abs(fly_vel) > vel_thresh & abs(fly_vel) < vel_max & abs(bump_vel) < bump_thresh & rho_i > rho_thresh;
end

function [gain_cue,gain_fly] = fly_gain_v3(all_data,trial_list,mu_cache,rho_cache,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max,use_intercept)
    % pools trial_gain_vectors_v3 across every trial in trial_list (one
    % fly's own trials within a single genotype x light-condition group),
    % then fits bump_vel ~ predictor_vel either with an intercept
    % (use_intercept=true, matching lpsp_kir_claude.m's own gain-scatter
    % convention -- an intercept absorbs any constant offset between the
    % two, e.g. from a residual lag/calibration mismatch, so the slope
    % alone isn't biased by it) or through the origin (use_intercept=false
    % -- a different question, "bump velocity per unit predictor velocity
    % with no allowance for a constant offset", more sensitive to such an
    % offset since nothing else can absorb it).
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        ti = trial_list(k);
        trial = all_data(ti);
        [fv,cv,bv,vd] = trial_gain_vectors_v3(trial,mu_cache{ti},rho_cache{ti},mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    % conservative activity floor: this fly's own pooled trials (in this
    % light condition) must show real turning (|fly_vel|>vel_thresh) for
    % at least min_frac_moving of ALL samples, not just the "valid"-gated
    % ones -- confirmed directly (gain_scratch_claude.m) on a fly moving
    % above threshold for only 0.25% of a 600s dark trial (91/36006
    % samples, confined to two brief blips): its "gain" was a regression
    % fit to a handful of noise-dominated samples (slope -3.245), not a
    % real tracking relationship. 1% is a deliberately loose floor -- this
    % is about excluding near-total quiescence, not borderline cases.
    min_frac_moving = 0.01;
    n_valid = sum(valid);
    if n_valid < 50 || mean(abs(fly_vel) > vel_thresh) < min_frac_moving
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    if use_intercept
        b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid);
        gain_cue = b_cue(2);
        b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid);
        gain_fly = b_fly(2);
    else
        % a through-origin slope is sum(x.*y)/sum(x.^2) -- with no
        % intercept to absorb it, a predictor whose variance collapses
        % toward zero (possible for cue_vel specifically, since "valid"
        % above constrains fly_vel/bump_vel/rho but not cue_vel itself)
        % drives sum(x.^2) toward zero and the slope toward +/-Inf.
        min_predictor_std = 0.01; % rad/s
        if std(cue_vel(valid)) < min_predictor_std || std(fly_vel(valid)) < min_predictor_std
            gain_cue = nan;
            gain_fly = nan;
            return
        end
        gain_cue = cue_vel(valid) \ bump_vel(valid);
        gain_fly = fly_vel(valid) \ bump_vel(valid);
    end
end

function export_fig(name)
    % saves the CURRENT figure both as a quick-preview PNG
    % (ugly_figures/exports/<name>.png) and as a PDF
    % (ugly_figures/rnai/<name>.pdf) -- one call site for every figure in
    % this script, so adding/changing an output format only needs editing
    % here once. Each write is wrapped so a locked destination (e.g. the
    % PDF open in a viewer) only warns and skips that one file instead of
    % throwing and aborting the rest of this multi-minute script -- a real
    % failure mode hit repeatedly in practice, since these exported files
    % get reopened often to check them.
    try
        exportgraphics(gcf,fullfile('ugly_figures','exports',[name '.png']),'Resolution',150)
    catch ME
        warning('could not save %s.png (%s) -- is it open in another program?', name, ME.message);
    end
    try
        exportgraphics(gcf,fullfile('ugly_figures','rnai',[name '.pdf']))
    catch ME
        warning('could not save %s.pdf (%s) -- is it open in another program?', name, ME.message);
    end
end

function y = wrap_break(x)
    % unwrap-safe display helper: wraps an unwrapped angle to -pi:pi and
    % breaks the line at the resulting circular jumps (NaN-inserted) so a
    % plot doesn't connect across wraps.
    y = mod(x,2*pi);
    y(y>pi) = y(y>pi) - 2*pi;
    y(find(abs(diff(y))>pi)+1) = nan;
end

function plot_bump_heading_overlay(fig_num,trial,mu,fz,heading_smooth_s,title_str,export_name)
    % imagesc + bump + heading overlay for one trial, as its own standalone
    % figure -- shared by figure 6 and any other single-trial diagnostic
    % that wants the same view. See plot_bump_heading_overlay_ax for the
    % actual drawing logic (used here, and reused as-is for multi-trial
    % subplot grids).
    figure(fig_num); clf
    set(gcf,'Name','PB activity + bump + heading overlay','Position',[100,100,1100,400])
    plot_bump_heading_overlay_ax(trial,mu,fz,heading_smooth_s,title_str)
    cb = colorbar; ylabel(cb,'z-score (5-frame moving average)')
    export_fig(export_name)
end

function plot_bump_heading_overlay_ax(trial,mu,fz,heading_smooth_s,title_str)
    % draws imagesc + bump + heading overlay into the CURRENT axes (no
    % figure()/clf(), no colorbar/export) -- call subplot(...) or
    % nexttile(...) first to pick where this lands. alpha unwrapped to the
    % continuous -pi:~3pi range (im.alpha repeats the same -pi:pi sequence
    % once per PB hemisphere), bump/heading each drawn twice (as-is, and
    % shifted +2*pi) so the overlay tracks both hemisphere bands, not just
    % one. mu is used AS-IS (already smoothed upstream, via
    % trial_pva_movmean, on the image timebase). heading = -ft.cue,
    % smoothed at heading_smooth_s on the fictrac timebase, THEN
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

    hold on
    imagesc(1:n_im,alpha_disp,fz)
    set(gca,'YDir','normal','XLim',[1,n_im],'YLim',[alpha_disp(1),alpha_disp(end)])
    plot(1:n_im,mu_plot,'w','LineWidth',1)
    plot(1:n_im,mu_plot+2*pi,'w','LineWidth',1)
    plot(1:n_im,heading_plot,'Color',[0,1,1],'LineWidth',1)
    plot(1:n_im,heading_plot+2*pi,'Color',[0,1,1],'LineWidth',1)
    xlabel('imaging frame'); ylabel('PB angle (rad)')
    title(title_str,'Interpreter','none')
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

function [speed_l,amp_l] = trial_speed_amp_v2(trial,fz,lag)
    % |rotational speed| (behavior) and peak bump amplitude (max z-score
    % across glomeruli each imaging frame, not a model fit), both on the
    % fictrac timebase, with amplitude shifted `lag` frames later than
    % speed -- same shift-and-trim convention as trial_gain_vectors_v3
    % above. fz is hemisphere-resolved (confirmed: im.alpha's two
    % 16-element halves are identical), so the two hemispheres are
    % averaged together first, matching lpsp_tnt_claude.m's own
    % trial_speed_amp for this same data shape.
    xf   = trial.ft.xf;
    n_im = size(fz,2);
    xb   = linspace(xf(1),xf(end),n_im)';

    n_hemi  = size(fz,1)/2;
    z_hemi  = (fz(1:n_hemi,:) + fz(n_hemi+1:end,:))/2;
    peak_im = max(z_hemi,[],1)';

    speed_full = abs(trial.ft.r_speed);
    amp_full   = interp1(xb,peak_im,xf,'linear','extrap');

    if lag == 0
        speed_l = speed_full;
        amp_l   = amp_full;
    elseif lag > 0
        speed_l = speed_full(1:end-lag);
        amp_l   = amp_full(lag+1:end);
    else
        speed_l = speed_full(-lag+1:end);
        amp_l   = amp_full(1:end+lag);
    end
end

function binned = fly_binned_amp(geno,dark_val,genotype,is_dark,fly_num,chunk_speed,chunk_amp,speed_edges,speed_x,min_bin_n_fly)
    % per-fly binned amplitude-vs-|speed| tuning curve for ONE genotype,
    % restricted to one light condition (dark_val = 0 or 1) or pooling
    % both if dark_val is []. one row per fly; a bin is left NaN for a fly
    % whose own pooled samples don't clear min_bin_n_fly there.
    if isempty(dark_val)
        rows = find(strcmp(genotype,geno));
    else
        rows = find(strcmp(genotype,geno) & is_dark==dark_val);
    end
    these_flies = unique(fly_num(rows));
    n_f = numel(these_flies);
    binned = nan(n_f,numel(speed_x));
    for ff = 1:n_f
        trial_list = rows(fly_num(rows)==these_flies(ff));
        speed_f = []; amp_f = [];
        for ti = trial_list(:)'
            speed_f = [speed_f; chunk_speed{ti}]; %#ok<AGROW>
            amp_f   = [amp_f; chunk_amp{ti}]; %#ok<AGROW>
        end
        for j = 1:numel(speed_x)
            idx = speed_f>=speed_edges(j) & speed_f<speed_edges(j+1);
            if sum(idx) >= min_bin_n_fly
                binned(ff,j) = mean(amp_f(idx),'omitnan');
            end
        end
    end
end

function plot_amp_by_genotype(binned_by_geno, geno_labels, geno_colors, speed_x)
    % binned_by_geno{k}: [n_flies x numel(speed_x)] binned amplitude
    % curves for genotype k. draws every fly's own curve faintly in that
    % genotype's color (not a distinct color per fly), plus a thick mean
    % +/- SEM line per genotype on top, all on the current axes -- one
    % shared plot per genotype comparison, not one subplot per genotype.
    % n is embedded in the legend label (same rationale as this script's
    % own groupplot: a floating text() n= label risks colliding with data).
    hold on
    h = gobjects(1,numel(geno_labels));
    legend_labels = cell(1,numel(geno_labels));
    for k = 1:numel(geno_labels)
        binned = binned_by_geno{k};
        for ff = 1:size(binned,1)
            plot(speed_x,binned(ff,:),'-','Color',[geno_colors(k,:),0.25],'LineWidth',0.75)
        end
        m = mean(binned,1,'omitnan');
        s = std(binned,0,1,'omitnan') ./ sqrt(sum(~isnan(binned),1));
        h(k) = errorbar(speed_x,m,s,'-o','Color',geno_colors(k,:),'LineWidth',2.5, ...
            'MarkerFaceColor',geno_colors(k,:),'MarkerSize',5);
        legend_labels{k} = sprintf('%s (n=%d)',geno_labels{k},size(binned,1));
    end
    legend(h,legend_labels,'Location','best','Interpreter','none')
    xlabel('|rotational speed| (rad/s)')
    ylabel('peak bump amplitude')
end

function [mu_ok,cue_ok] = trial_bumpok_samples_v2(trial,mu,rho,heading_smooth_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh)
    % pooled (mu, heading) samples for one trial, restricted to
    % trial_walking_bouts_v2's own is_walking mask -- the SAME
    % "discernable bump" restriction this script's own offset-variability
    % figures (13/14) already use, reused here rather than inventing a
    % separate instantaneous threshold (see this function's call site for
    % why). mu is interpolated onto the fictrac timebase and used AS-IS
    % (already smoothed upstream); heading is -ft.cue smoothed at
    % heading_smooth_s via movmean, same construction as everywhere else
    % in this script that needs a smoothed heading trace. No lag shift
    % between the two (see call site).
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu);
    xb = linspace(xf(1),xf(end),n_im)';

    mu_i = interp1(xb,unwrap(mu(:)),xf,'linear','extrap');

    heading_win = max(1,round(heading_smooth_s/dt));
    heading_smooth_full = movmean(unwrap(-trial.ft.cue(:)),heading_win);

    [~,~,~,~,is_walking] = trial_walking_bouts_v2(trial,mu,rho,heading_smooth_s,turn_thresh,max_gap_s,min_walking_s,bout_rho_thresh);

    valid = is_walking & ~isnan(mu_i) & ~isnan(heading_smooth_full);
    mu_ok  = mu_i(valid);
    cue_ok = heading_smooth_full(valid);
end

function H = circ_entropy(x, edges)
    % Shannon entropy (natural log, nats) of a wrapped circular variable's
    % histogram -- a modality-agnostic "how peaky is this distribution"
    % measure, unlike circ_var (1 - resultant vector length), which is
    % based on VECTOR AVERAGING and can be fooled by a symmetric
    % multimodal distribution: two sharp peaks exactly opposite each other
    % on the circle sum to a near-zero resultant vector (circ_var near 1,
    % "looks nearly uniform") even though the distribution is actually
    % tightly concentrated in two spots, not spread out at all. Entropy
    % only looks at how mass is distributed across bins, not their angular
    % position, so it doesn't have this blind spot.
    p = histcounts(x,edges,'Normalization','probability');
    p = p(p>0); % 0*log(0) is defined as 0 by convention -- just drop empty bins rather than computing 0*(-Inf) = NaN
    H = -sum(p.*log(p));
end
