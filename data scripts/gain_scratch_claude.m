%% gain_scratch_claude
% Scratch pipeline for a NEW gain metric: how tightly the bump (im.mu)
% tracks the visual cue / the fly's own rotation, in closed loop.
%
% This is a fresh gain definition, not a re-run of the velocity-gain
% analysis in lpsp_kir_claude.m (bump_vel ~ fly heading vel, via -ft.cue
% as a heading PROXY, because that dataset has no independent heading
% signal) -- kept here only as a reference for the regression convention
% (through-intercept fit, bump_vel(lag) ~ 1 + predictor_vel, same
% vel/bump/rho thresholds as a starting point) and the general "sweep a
% smoothing window, pick the value that optimizes some criterion on
% control flies" pattern used again below.
%
% This dataset (lpsp_rnai_joint_nosmooth_reg_20260831.mat, see
% lpsp_rnai_claude_v2.m for the full labeling convention) has BOTH signals
% separately, confirmed directly on trial 1:
%   ft.heading -- the fly's own rotation, essentially cumtrapz(r_speed)
%                 (slope of gradient(heading) ~ r_speed is 0.99, corr .99)
%   ft.cue     -- the closed-loop visual scene position. NOT well
%                 correlated with heading instantaneously in this dataset
%                 (corr(diff(heading),diff(unwrap(-cue))) ~ 0.02 on trial
%                 1) -- cue is the scene position under a fixed 0.8
%                 closed-loop gain relative to the fly's own rotation,
%                 not a clean heading proxy the way -ft.cue was used in
%                 lpsp_kir_claude.m's dataset (which had no ft.heading at
%                 all).
% So there are two candidate "gain" definitions here:
%   gain_cue = slope(bump_vel ~ 1 + cue_vel)  -- bump vs. the visual scene
%              it's actually tracking; expected ~1 for perfect tracking.
%   gain_fly = slope(bump_vel ~ 1 + fly_vel)  -- bump vs. the fly's own
%              rotation (r_speed); expected ~0.8, since the SCENE itself
%              only moves at 0.8x the fly's own rotation in closed loop,
%              so a bump that tracks the scene perfectly only tracks the
%              fly's own turning at 0.8x.
% gain_cue is used as the OPTIMIZATION criterion below (cue is the more
% direct ground truth for "is the bump tracking what it's supposed to
% track", and isn't itself a fictrac measurement subject to the same
% sensor noise as r_speed); gain_fly is reported alongside it at the end
% for comparison, at the SAME optimized parameters.
%
% ft.cue has a small (~0.1-0.4%, confirmed on the first 20 trials) rate of
% scattered NaN dropouts mid-trial (not an edge-truncation artifact) --
% linearly interpolated over (fill_nan_gaps, below) before any
% unwrap/smooth/diff, since unwrap+gradient would otherwise propagate a
% single dropped sample into a much wider corrupted stretch.
%
% Three smoothing stages are swept, SEQUENTIALLY/GREEDILY (fix stage
% N+1..3 at a neutral default, sweep stage N, fix it at its winner, sweep
% stage N+1, ...) rather than as a full 3-D grid, both for tractability
% and because this mirrors lpsp_rnai_claude_v2.m's own step 3 -> step 3b
% chaining (mu fixed while heading is swept, then heading's winner reused
% while an additional mu smoothing is swept):
%   1) im.f smoothing (frames, movmean along time) before z-score + PVA
%      -> re-derives mu/rho from scratch each candidate (same "plain
%      sample-count movmean" convention as lpsp_rnai_claude_v2.m's
%      trial_pva_movmean, deliberately simpler than a gaussian/time
%      converted smooth).
%   2) additional mu smoothing (seconds, movmean on the unwrapped bump
%      POSITION trace, on the imaging timebase, converted to samples via
%      that trial's own image dt) -- same construction as
%      lpsp_rnai_claude_v2.m step 3b.
%   3) heading smoothing (seconds, applied identically to r_speed
%      directly -- already a velocity -- and to the unwrapped cue
%      POSITION before differentiating -- smoothing a signal by window W
%      then differentiating is the same real-time operation as smoothing
%      its derivative by W directly, same reasoning
%      lpsp_rnai_claude_v2.m used for its own heading smoothing) --
%      converted to samples via that trial's own fictrac dt
%      (median(diff(xf)), which is NOT constant across every trial in
%      this dataset -- same reasoning as lpsp_rnai_claude_v2.m).
%
% At each candidate, the OPTIMIZATION CRITERION is the per-fly gain_cue's
% mean squared error from the target value of 1, pooled across BOTH empty
% controls (empty>th, empty>vglut), closed loop only:
%   crit = mean( (fly_gain_cue - 1).^2 )   -- one value per fly, pooled.
% MSE-from-target decomposes into bias^2 + variance, so minimizing it
% favors whichever smoothing choice makes the control flies' own gain
% both closer to 1 on average AND more tightly clustered there -- exactly
% "optimally tight around a value of 1", not just "unbiased" or just
% "low-variance" alone.
addpath(fullfile(pwd,'circ_stats'));

%% 1) load data
data_dir    = fullfile('.data'); % matches lpsp_kir_claude.m -- this repo checkout has ".data", not "data" (lpsp_rnai_claude_v2.m's "data" reference is stale)
source_file = 'lpsp_rnai_joint_nosmooth_reg_20260831.mat';

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% 2) genotype (driver x target) and light condition per trial -- identical parsing to lpsp_rnai_claude_v2.m
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
[fly_list,~,fly_num] = unique(fly_id); %#ok<ASGLU>
fprintf('%d trials -> %d flies\n', n_trials, numel(fly_list));

%% 3) restrict to empty controls (empty>th, empty>vglut), closed loop only
opt_geno   = {'empty>th','empty>vglut'};
opt_trials = find(ismember(genotype,opt_geno) & ~is_dark);
opt_flies  = unique(fly_num(opt_trials));
fprintf('\noptimizing on %d empty-control flies (%d closed-loop trials: %d empty>th, %d empty>vglut)\n', ...
    numel(opt_flies), numel(opt_trials), sum(strcmp(genotype(opt_trials),'empty>th')), sum(strcmp(genotype(opt_trials),'empty>vglut')));

opt_fly_trials = arrayfun(@(f) opt_trials(fly_num(opt_trials)==f), opt_flies, 'UniformOutput',false);

%% 4) fixed thresholds for the gain regression (inherited starting point from lpsp_kir_claude.m's own gain-scatter section, not re-tuned here)
vel_thresh  = 0.2; % rad/s -- exclude near-stationary samples where velocity noise dominates the ratio
bump_thresh = 10;  % rad/s -- exclude gradient/unwrap artifacts in bump velocity
rho_thresh  = 0.2; % im.rho (re-derived per candidate) -- exclude samples where the bump is too weak to trust
vel_max     = 5;   % rad/s -- exclude fly_vel outliers (tracking glitches)

%% 5) sweep 1/3: im.f smoothing (frames), holding mu/heading smoothing at neutral defaults
f_smooth_candidates   = [1,3,5,7,9,11,15,21,31];
default_mu_smooth_s      = 0;
default_heading_smooth_s = 0.1;

n_f_cand = numel(f_smooth_candidates);
fly_crit_f = nan(numel(opt_flies),n_f_cand);
fly_gain_cue_f = nan(numel(opt_flies),n_f_cand);

fprintf('\n=== sweep 1/3: im.f smoothing (frames) ===\n');
for c = 1:n_f_cand
    f_smooth = f_smooth_candidates(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gain_cue,gain_fly] = fly_gain(all_data,trial_list,f_smooth,default_mu_smooth_s,default_heading_smooth_s, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max); %#ok<ASGLU>
        fly_gain_cue_f(k,c) = gain_cue;
    end
    fly_crit_f(:,c) = (fly_gain_cue_f(:,c)-1).^2;
    fprintf('  f_smooth=%2d frames: mean gain_cue=%.3f, mean MSE-from-1=%.4f (n=%d flies)\n', ...
        f_smooth, mean(fly_gain_cue_f(:,c),'omitnan'), mean(fly_crit_f(:,c),'omitnan'), sum(~isnan(fly_gain_cue_f(:,c))));
end
mean_crit_f = mean(fly_crit_f,1,'omitnan');
[~,best_fc] = min(mean_crit_f);
f_smooth_opt = f_smooth_candidates(best_fc);
fprintf('optimal im.f smoothing: %d frames (mean MSE-from-1=%.4f)\n', f_smooth_opt, mean_crit_f(best_fc));

figure(1); clf
set(gcf,'Name','sweep 1/3: im.f smoothing','Position',[100,100,700,500])
subplot(2,1,1); hold on
plot(f_smooth_candidates,mean(fly_gain_cue_f,1,'omitnan'),'-ok','MarkerFaceColor','k')
yline(1,':k')
xline(f_smooth_opt,'--r')
xlabel('im.f smoothing window (frames)'); ylabel('mean gain\_cue (target=1)')
title('empty>th + empty>vglut, closed loop')
subplot(2,1,2); hold on
plot(f_smooth_candidates,mean_crit_f,'-ok','MarkerFaceColor','k')
plot(f_smooth_opt,mean_crit_f(best_fc),'o','MarkerSize',12,'Color','r','LineWidth',2)
xlabel('im.f smoothing window (frames)'); ylabel('mean squared error from gain\_cue=1')
title(sprintf('optimum = %d frames',f_smooth_opt))

%% 6) sweep 2/3: additional mu smoothing (seconds), holding im.f smoothing at its winner
mu_smooth_candidates_s = [0,0.05,0.1,0.2,0.35,0.5,0.75,1,1.5,2];

n_mu_cand = numel(mu_smooth_candidates_s);
fly_crit_mu = nan(numel(opt_flies),n_mu_cand);
fly_gain_cue_mu = nan(numel(opt_flies),n_mu_cand);

fprintf('\n=== sweep 2/3: additional mu smoothing (s), im.f smoothing fixed at %d frames ===\n',f_smooth_opt);
for c = 1:n_mu_cand
    mu_smooth_s = mu_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gain_cue,gain_fly] = fly_gain(all_data,trial_list,f_smooth_opt,mu_smooth_s,default_heading_smooth_s, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max); %#ok<ASGLU>
        fly_gain_cue_mu(k,c) = gain_cue;
    end
    fly_crit_mu(:,c) = (fly_gain_cue_mu(:,c)-1).^2;
    fprintf('  mu_smooth=%.2fs: mean gain_cue=%.3f, mean MSE-from-1=%.4f (n=%d flies)\n', ...
        mu_smooth_s, mean(fly_gain_cue_mu(:,c),'omitnan'), mean(fly_crit_mu(:,c),'omitnan'), sum(~isnan(fly_gain_cue_mu(:,c))));
end
mean_crit_mu = mean(fly_crit_mu,1,'omitnan');
[~,best_muc] = min(mean_crit_mu);
mu_smooth_opt_s = mu_smooth_candidates_s(best_muc);
fprintf('optimal additional mu smoothing: %.2fs (mean MSE-from-1=%.4f)\n', mu_smooth_opt_s, mean_crit_mu(best_muc));

figure(2); clf
set(gcf,'Name','sweep 2/3: mu smoothing','Position',[100,100,700,500])
subplot(2,1,1); hold on
plot(mu_smooth_candidates_s,mean(fly_gain_cue_mu,1,'omitnan'),'-ok','MarkerFaceColor','k')
yline(1,':k')
xline(mu_smooth_opt_s,'--r')
xlabel('additional mu smoothing window (s)'); ylabel('mean gain\_cue (target=1)')
title('empty>th + empty>vglut, closed loop')
subplot(2,1,2); hold on
plot(mu_smooth_candidates_s,mean_crit_mu,'-ok','MarkerFaceColor','k')
plot(mu_smooth_opt_s,mean_crit_mu(best_muc),'o','MarkerSize',12,'Color','r','LineWidth',2)
xlabel('additional mu smoothing window (s)'); ylabel('mean squared error from gain\_cue=1')
title(sprintf('optimum = %.2fs',mu_smooth_opt_s))

%% 7) sweep 3/3: heading smoothing (seconds), holding im.f + mu smoothing at their winners
heading_smooth_candidates_s = [0,0.05,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];

n_hd_cand = numel(heading_smooth_candidates_s);
fly_crit_hd = nan(numel(opt_flies),n_hd_cand);
fly_gain_cue_hd = nan(numel(opt_flies),n_hd_cand);
fly_gain_fly_hd = nan(numel(opt_flies),n_hd_cand);

fprintf('\n=== sweep 3/3: heading smoothing (s), im.f=%d frames, mu smoothing=%.2fs fixed ===\n',f_smooth_opt,mu_smooth_opt_s);
for c = 1:n_hd_cand
    heading_smooth_s = heading_smooth_candidates_s(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gain_cue,gain_fly] = fly_gain(all_data,trial_list,f_smooth_opt,mu_smooth_opt_s,heading_smooth_s, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_gain_cue_hd(k,c) = gain_cue;
        fly_gain_fly_hd(k,c) = gain_fly;
    end
    fly_crit_hd(:,c) = (fly_gain_cue_hd(:,c)-1).^2;
    fprintf('  heading_smooth=%.2fs: mean gain_cue=%.3f, mean gain_fly=%.3f, mean MSE-from-1=%.4f (n=%d flies)\n', ...
        heading_smooth_s, mean(fly_gain_cue_hd(:,c),'omitnan'), mean(fly_gain_fly_hd(:,c),'omitnan'), ...
        mean(fly_crit_hd(:,c),'omitnan'), sum(~isnan(fly_gain_cue_hd(:,c))));
end
mean_crit_hd = mean(fly_crit_hd,1,'omitnan');
[~,best_hdc] = min(mean_crit_hd);
heading_smooth_opt_s = heading_smooth_candidates_s(best_hdc);
fprintf('optimal heading smoothing: %.2fs (mean MSE-from-1=%.4f)\n', heading_smooth_opt_s, mean_crit_hd(best_hdc));

figure(3); clf
set(gcf,'Name','sweep 3/3: heading smoothing','Position',[100,100,700,500])
subplot(2,1,1); hold on
plot(heading_smooth_candidates_s,mean(fly_gain_cue_hd,1,'omitnan'),'-ok','MarkerFaceColor','k')
plot(heading_smooth_candidates_s,mean(fly_gain_fly_hd,1,'omitnan'),'-o','Color',[0,0.4470,0.7410],'MarkerFaceColor',[0,0.4470,0.7410])
yline(1,':k'); yline(0.8,':','Color',[0,0.4470,0.7410])
xline(heading_smooth_opt_s,'--r')
legend({'gain\_cue (target=1)','gain\_fly (target=0.8)'},'Location','best')
xlabel('heading smoothing window (s)'); ylabel('mean gain')
title('empty>th + empty>vglut, closed loop')
subplot(2,1,2); hold on
plot(heading_smooth_candidates_s,mean_crit_hd,'-ok','MarkerFaceColor','k')
plot(heading_smooth_opt_s,mean_crit_hd(best_hdc),'o','MarkerSize',12,'Color','r','LineWidth',2)
xlabel('heading smoothing window (s)'); ylabel('mean squared error from gain\_cue=1')
title(sprintf('optimum = %.2fs',heading_smooth_opt_s))

fprintf('\n=== optimized pipeline ===\n');
fprintf('  im.f smoothing:       %d frames\n', f_smooth_opt);
fprintf('  additional mu smooth: %.2fs\n', mu_smooth_opt_s);
fprintf('  heading smoothing:    %.2fs\n', heading_smooth_opt_s);

%% 8) SUPERSEDED by section 9 below: the movmean sweep above (sections 5-7) never actually converged --
% heading_smooth_opt_s kept climbing to whatever the edge of its own
% candidate grid happened to be (section 7's printed sweep), rather than
% settling on an interior optimum, and the resulting gain_cue never got
% close to 1. Diagnosing this on trial 535 (an empty>th trial hand-picked
% for high rotational speed AND the lowest bump-cue POSITION offset
% variability among fast-turning empty-control trials -- i.e. a trial
% where the bump visibly tracks the cue well) showed why: differentiating
% (gradient()) a movmean-smoothed position trace to get velocity
% amplifies frame-to-frame noise enough to wash out a real tracking
% relationship, and there was also a real (if small) LAG between cue and
% bump that this section never accounted for. A follow-up 3-parameter
% grid on trial 535 alone (GAUSSIAN bump-smoothing width x GAUSSIAN
% cue-smoothing width x lag, sweeping bump/cue width via smoothdata(...,
% 'gaussian',...) instead of movmean) found a real interior optimum
% (bump=cue~0.75-1.0s Gaussian width, lag~8 frames/0.13s, cue leading
% bump) instead of another edge effect. Section 9 applies that settled
% choice (GAUSSIAN, not movmean; PLUS a fixed lag) across every
% empty-control fly, not just trial 535.
%
% (Also confirmed directly, on this same cue trace: iterating a Gaussian
% smooth N times with window W is empirically -- not just in the
% idealized-infinite-kernel theory -- equivalent to ONE pass at window
% W*sqrt(N) here, corr>0.9999 in every case checked; edge/truncation
% effects are negligible at these window sizes relative to trial length.
% So "1s single pass" below is the settled width itself, not an arbitrary
% stand-in for some larger multi-pass search.)

%% 9) settled pipeline: 1s Gaussian smoothing (single pass) on both bump (mu) and cue/fly_vel, lag=8 frames -- applied across every empty-control fly
% im.f smoothing stays at section 5's own winner (f_smooth_opt) -- that
% stage wasn't revisited by the trial-535 diagnostic.
bump_smooth_gaussian_s = 1; % s, single Gaussian pass on unwrapped mu (image timebase)
cue_smooth_gaussian_s  = 1; % s, single Gaussian pass on unwrapped cue AND on r_speed directly (fictrac timebase) -- same window for both, matching section 7's "one heading window for both predictors" convention
lag_frames = 8;             % fictrac frames (~0.13s): cue_vel/fly_vel led bump_vel by this much on trial 535, stable across that trial's whole bump/cue-width grid

fly_gain_cue_final = nan(numel(opt_flies),1);
fly_gain_fly_final = nan(numel(opt_flies),1);
fly_n_valid_final  = nan(numel(opt_flies),1);
fly_geno_final     = cell(numel(opt_flies),1);
for k = 1:numel(opt_flies)
    trial_list = opt_fly_trials{k};
    [gain_cue,gain_fly,n_valid] = fly_gain_gauss(all_data,trial_list,f_smooth_opt,bump_smooth_gaussian_s,cue_smooth_gaussian_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max);
    fly_gain_cue_final(k) = gain_cue;
    fly_gain_fly_final(k) = gain_fly;
    fly_n_valid_final(k)  = n_valid;
    fly_geno_final{k}     = genotype{trial_list(1)};
end

cat_x = ones(numel(opt_flies),1) + strcmp(fly_geno_final,'empty>vglut'); % 1=empty>th, 2=empty>vglut
cat_labels = {'empty>th','empty>vglut'};
cat_colors = lines(2);

figure(4); clf
set(gcf,'Name','final gain, settled Gaussian+lag pipeline','Position',[100,100,900,450])
subplot(1,2,1)
groupplot(cat_x,fly_gain_cue_final,cat_labels,cat_colors)
hold on; plot(xlim,[1,1],':k')
ylabel('gain\_cue = slope(bump\_vel ~ 1 + cue\_vel)')
title(sprintf('bump vs. visual cue (target=1), n=%d flies',sum(~isnan(fly_gain_cue_final))))

subplot(1,2,2)
groupplot(cat_x,fly_gain_fly_final,cat_labels,cat_colors)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ 1 + fly\_vel)')
title(sprintf('bump vs. fly rotation (target=0.8), n=%d flies',sum(~isnan(fly_gain_fly_final))))

sgtitle(sprintf('settled pipeline: im.f=%d frames, bump/cue Gaussian=%.1fs, lag=%d frames', ...
    f_smooth_opt,bump_smooth_gaussian_s,lag_frames),'Interpreter','none')

fprintf('\n=== final gain, settled Gaussian+lag pipeline (mean +/- SEM across flies) ===\n');
for g = 1:2
    y_cue = fly_gain_cue_final(cat_x==g);
    y_fly = fly_gain_fly_final(cat_x==g);
    fprintf('  %-12s gain_cue = %.3f +/- %.3f, gain_fly = %.3f +/- %.3f (n=%d flies)\n', ...
        cat_labels{g}, mean(y_cue,'omitnan'), std(y_cue,'omitnan')/sqrt(sum(~isnan(y_cue))), ...
        mean(y_fly,'omitnan'), std(y_fly,'omitnan')/sqrt(sum(~isnan(y_fly))), sum(~isnan(y_cue)));
end

%% 10) comprehensive 4-parameter sweep, now that the NaN-fill bug (fixed above) no longer swamps the search
% Section 9 applied one FIXED choice (1s Gaussian on everything, lag=8)
% found from a single-trial diagnostic. This sweeps all four smoothing
% stages properly, SEQUENTIALLY/GREEDILY (same pattern as sections 5-7 /
% lpsp_rnai_claude_v2.m's own step chaining): fix the later stages at a
% neutral default, sweep the current stage, fix it at its winner, move
% on. cue and fly_vel smoothing -- tied together as one "heading_smooth_s"
% in every earlier sweep in this script -- are swept INDEPENDENTLY here,
% since there's no a priori reason a visual-scene position signal and a
% fictrac ball-rotation signal need the same amount of smoothing. lag
% stays fixed at lag_frames (section 9, 8 frames/0.13s) -- not re-swept.
%
% Criterion: combined per-fly MSE-from-target, summing gain_cue's MSE
% from its target of 1 and gain_fly's MSE from its target of 0.8 -- so a
% smoothing choice only wins if it moves BOTH gains toward their own
% targets (not one at the other's expense).
combined_crit = @(gc,gf) (gc-1).^2 + (gf-0.8).^2;
default_smooth_s = 1; % neutral hold for whichever of mu/cue/vel isn't being swept yet, matching section 9's own settled value

%% 10a) stage 1/4: im.f smoothing (frames) -- re-derives mu/rho from scratch each candidate
f_smooth_candidates2 = [1,3,5,7,9,11,15,21];
n_c1 = numel(f_smooth_candidates2);
fly_gc1 = nan(numel(opt_flies),n_c1); fly_gf1 = nan(numel(opt_flies),n_c1);

fprintf('\n=== 4-param sweep, stage 1/4: im.f smoothing (frames) ===\n');
for c = 1:n_c1
    fsm = f_smooth_candidates2(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gc,gf] = fly_gain_gauss2(all_data,trial_list,fsm,default_smooth_s,default_smooth_s,default_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_gc1(k,c) = gc; fly_gf1(k,c) = gf;
    end
    crit_c = mean(combined_crit(fly_gc1(:,c),fly_gf1(:,c)),'omitnan');
    fprintf('  f_smooth=%2d frames: mean gain_cue=%.3f, mean gain_fly=%.3f, combined MSE=%.4f\n', ...
        fsm, mean(fly_gc1(:,c),'omitnan'), mean(fly_gf1(:,c),'omitnan'), crit_c);
end
crit1 = mean(combined_crit(fly_gc1,fly_gf1),1,'omitnan');
[~,best_c1] = min(crit1);
f_smooth_opt2 = f_smooth_candidates2(best_c1);
fprintf('winner: im.f smoothing = %d frames (combined MSE=%.4f)\n', f_smooth_opt2, crit1(best_c1));

%% 10b) precompute the PVA (mu_raw, rho_raw) cache at f_smooth_opt2 for every opt_trial -- reused by stages 2-4 so they don't redo this every candidate
pva_mu_cache  = cell(n_trials,1);
pva_rho_cache = cell(n_trials,1);
for ti = opt_trials'
    [pva_mu_cache{ti},pva_rho_cache{ti}] = trial_pva(all_data(ti),f_smooth_opt2);
end

%% 10c) stage 2/4: bump (mu) Gaussian smoothing window (s), im.f fixed at its winner
mu_smooth_candidates2_s = [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];
n_c2 = numel(mu_smooth_candidates2_s);
fly_gc2 = nan(numel(opt_flies),n_c2); fly_gf2 = nan(numel(opt_flies),n_c2);

fprintf('\n=== 4-param sweep, stage 2/4: bump (mu) Gaussian smoothing (s), im.f=%d frames fixed ===\n',f_smooth_opt2);
for c = 1:n_c2
    msm = mu_smooth_candidates2_s(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gc,gf] = fly_gain_cached(all_data,trial_list,pva_mu_cache,pva_rho_cache,msm,default_smooth_s,default_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_gc2(k,c) = gc; fly_gf2(k,c) = gf;
    end
    crit_c = mean(combined_crit(fly_gc2(:,c),fly_gf2(:,c)),'omitnan');
    fprintf('  mu_smooth=%.2fs: mean gain_cue=%.3f, mean gain_fly=%.3f, combined MSE=%.4f\n', ...
        msm, mean(fly_gc2(:,c),'omitnan'), mean(fly_gf2(:,c),'omitnan'), crit_c);
end
crit2 = mean(combined_crit(fly_gc2,fly_gf2),1,'omitnan');
[~,best_c2] = min(crit2);
mu_smooth_opt2 = mu_smooth_candidates2_s(best_c2);
fprintf('winner: bump (mu) smoothing = %.2fs (combined MSE=%.4f)\n', mu_smooth_opt2, crit2(best_c2));

%% 10d) stage 3/4: cue Gaussian smoothing window (s), im.f + mu fixed at their winners
cue_smooth_candidates2_s = [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];
n_c3 = numel(cue_smooth_candidates2_s);
fly_gc3 = nan(numel(opt_flies),n_c3); fly_gf3 = nan(numel(opt_flies),n_c3);

fprintf('\n=== 4-param sweep, stage 3/4: cue Gaussian smoothing (s), im.f=%d frames, mu=%.2fs fixed ===\n',f_smooth_opt2,mu_smooth_opt2);
for c = 1:n_c3
    csm = cue_smooth_candidates2_s(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gc,gf] = fly_gain_cached(all_data,trial_list,pva_mu_cache,pva_rho_cache,mu_smooth_opt2,csm,default_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_gc3(k,c) = gc; fly_gf3(k,c) = gf;
    end
    crit_c = mean(combined_crit(fly_gc3(:,c),fly_gf3(:,c)),'omitnan');
    fprintf('  cue_smooth=%.2fs: mean gain_cue=%.3f, mean gain_fly=%.3f, combined MSE=%.4f\n', ...
        csm, mean(fly_gc3(:,c),'omitnan'), mean(fly_gf3(:,c),'omitnan'), crit_c);
end
crit3 = mean(combined_crit(fly_gc3,fly_gf3),1,'omitnan');
[~,best_c3] = min(crit3);
cue_smooth_opt2 = cue_smooth_candidates2_s(best_c3);
fprintf('winner: cue smoothing = %.2fs (combined MSE=%.4f)\n', cue_smooth_opt2, crit3(best_c3));

%% 10e) stage 4/4: fly_vel (r_speed) Gaussian smoothing window (s), im.f + mu + cue fixed at their winners
vel_smooth_candidates2_s = [0,0.1,0.2,0.35,0.5,0.75,1,1.5,2,3];
n_c4 = numel(vel_smooth_candidates2_s);
fly_gc4 = nan(numel(opt_flies),n_c4); fly_gf4 = nan(numel(opt_flies),n_c4);

fprintf('\n=== 4-param sweep, stage 4/4: fly_vel (r_speed) Gaussian smoothing (s), im.f=%d frames, mu=%.2fs, cue=%.2fs fixed ===\n', ...
    f_smooth_opt2,mu_smooth_opt2,cue_smooth_opt2);
for c = 1:n_c4
    vsm = vel_smooth_candidates2_s(c);
    for k = 1:numel(opt_flies)
        trial_list = opt_fly_trials{k};
        [gc,gf] = fly_gain_cached(all_data,trial_list,pva_mu_cache,pva_rho_cache,mu_smooth_opt2,cue_smooth_opt2,vsm,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_gc4(k,c) = gc; fly_gf4(k,c) = gf;
    end
    crit_c = mean(combined_crit(fly_gc4(:,c),fly_gf4(:,c)),'omitnan');
    fprintf('  vel_smooth=%.2fs: mean gain_cue=%.3f, mean gain_fly=%.3f, combined MSE=%.4f\n', ...
        vsm, mean(fly_gc4(:,c),'omitnan'), mean(fly_gf4(:,c),'omitnan'), crit_c);
end
crit4 = mean(combined_crit(fly_gc4,fly_gf4),1,'omitnan');
[~,best_c4] = min(crit4);
vel_smooth_opt2 = vel_smooth_candidates2_s(best_c4);
fprintf('winner: fly_vel smoothing = %.2fs (combined MSE=%.4f)\n', vel_smooth_opt2, crit4(best_c4));

fprintf('\n=== 4-param sweep: winning pipeline ===\n');
fprintf('  im.f smoothing:   %d frames\n', f_smooth_opt2);
fprintf('  bump (mu) smooth: %.2fs\n', mu_smooth_opt2);
fprintf('  cue smooth:       %.2fs\n', cue_smooth_opt2);
fprintf('  fly_vel smooth:   %.2fs\n', vel_smooth_opt2);
fprintf('  lag (fixed):      %d frames\n', lag_frames);

figure(5); clf
set(gcf,'Name','4-param sweep: stage-by-stage curves','Position',[100,100,1200,350])
tl5 = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
nexttile(tl5); hold on
plot(f_smooth_candidates2,crit1,'-ok','MarkerFaceColor','k')
plot(f_smooth_opt2,crit1(best_c1),'o','MarkerSize',10,'Color','r','LineWidth',2)
xlabel('im.f smoothing (frames)'); ylabel('combined MSE'); title('stage 1: im.f')
nexttile(tl5); hold on
plot(mu_smooth_candidates2_s,crit2,'-ok','MarkerFaceColor','k')
plot(mu_smooth_opt2,crit2(best_c2),'o','MarkerSize',10,'Color','r','LineWidth',2)
xlabel('bump smoothing (s)'); title('stage 2: bump (mu)')
nexttile(tl5); hold on
plot(cue_smooth_candidates2_s,crit3,'-ok','MarkerFaceColor','k')
plot(cue_smooth_opt2,crit3(best_c3),'o','MarkerSize',10,'Color','r','LineWidth',2)
xlabel('cue smoothing (s)'); title('stage 3: cue')
nexttile(tl5); hold on
plot(vel_smooth_candidates2_s,crit4,'-ok','MarkerFaceColor','k')
plot(vel_smooth_opt2,crit4(best_c4),'o','MarkerSize',10,'Color','r','LineWidth',2)
xlabel('fly\_vel smoothing (s)'); title('stage 4: fly\_vel')
title(tl5,'4-param sweep: combined MSE-from-target at each stage (holding later stages at a 1s default, earlier stages at their own winners)','Interpreter','none')

%% 10f) final per-fly gain at the winning 4-parameter combo, group plot
fly_gain_cue_final2 = nan(numel(opt_flies),1);
fly_gain_fly_final2 = nan(numel(opt_flies),1);
fly_geno_final2     = cell(numel(opt_flies),1);
for k = 1:numel(opt_flies)
    trial_list = opt_fly_trials{k};
    [gc,gf] = fly_gain_cached(all_data,trial_list,pva_mu_cache,pva_rho_cache,mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max);
    fly_gain_cue_final2(k) = gc;
    fly_gain_fly_final2(k) = gf;
    fly_geno_final2{k}     = genotype{trial_list(1)};
end

cat_x2 = ones(numel(opt_flies),1) + strcmp(fly_geno_final2,'empty>vglut');

figure(6); clf
set(gcf,'Name','final gain, 4-param sweep winner','Position',[100,100,900,450])
subplot(1,2,1)
groupplot(cat_x2,fly_gain_cue_final2,cat_labels,cat_colors)
hold on; plot(xlim,[1,1],':k')
ylabel('gain\_cue = slope(bump\_vel ~ 1 + cue\_vel)')
title(sprintf('bump vs. visual cue (target=1), n=%d flies',sum(~isnan(fly_gain_cue_final2))))

subplot(1,2,2)
groupplot(cat_x2,fly_gain_fly_final2,cat_labels,cat_colors)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ 1 + fly\_vel)')
title(sprintf('bump vs. fly rotation (target=0.8), n=%d flies',sum(~isnan(fly_gain_fly_final2))))

sgtitle(sprintf('4-param sweep winner: im.f=%d frames, bump=%.2fs, cue=%.2fs, fly\\_vel=%.2fs, lag=%d frames', ...
    f_smooth_opt2,mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames),'Interpreter','none')

fprintf('\n=== final gain, 4-param sweep winner (mean +/- SEM across flies) ===\n');
for g = 1:2
    y_cue = fly_gain_cue_final2(cat_x2==g);
    y_fly = fly_gain_fly_final2(cat_x2==g);
    fprintf('  %-12s gain_cue = %.3f +/- %.3f, gain_fly = %.3f +/- %.3f (n=%d flies)\n', ...
        cat_labels{g}, mean(y_cue,'omitnan'), std(y_cue,'omitnan')/sqrt(sum(~isnan(y_cue))), ...
        mean(y_fly,'omitnan'), std(y_fly,'omitnan')/sqrt(sum(~isnan(y_fly))), sum(~isnan(y_cue)));
end

%% 11) apply the winning 4-parameter pipeline to EVERY genotype x light-condition group in the dataset
% section 10 only optimized/reported on empty controls, closed loop
% (opt_geno/opt_trials). This applies that SAME winning pipeline
% (f_smooth_opt2, mu_smooth_opt2, cue_smooth_opt2, vel_smooth_opt2,
% lag_frames -- none re-tuned here) to every driver x target genotype
% (including lpsp>th/lpsp>vglut/lpsp>mcherry, not just the empty
% controls) in BOTH light conditions (closed loop AND dark) -- same
% group_defs/cat_colors convention as lpsp_kir_claude.m's own
% genotype-major grouping (each genotype's CL/dark pair shares one
% color). Dark trials are included deliberately: gain_cue should be
% close to 0 there (no visual cue driving the panel in the dark, so
% cue_vel carries no real signal for the bump to track), which is a
% useful built-in negative control on this pipeline rather than
% something to filter out.
geno_order_all = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};
cond_label_all = {'closed loop','dark'};

group_defs_all = struct('geno',{},'dark',{},'label',{});
% light-condition-major ordering (all closed-loop groups first, then all
% dark groups) so the two blocks can be compared side by side -- color
% (below) still ties each genotype's two entries together across the
% resulting gap, since it's assigned by genotype, not by position.
for ci = 1:numel(cond_label_all)
    for gi = 1:numel(geno_order_all)
        group_defs_all(end+1) = struct('geno',geno_order_all{gi},'dark',ci-1,'label',sprintf('%s (%s)',geno_order_all{gi},cond_label_all{ci})); %#ok<SAGROW>
    end
end

group_fly_list_all   = cell(1,numel(group_defs_all));
group_fly_trials_all = cell(1,numel(group_defs_all));
for gd = 1:numel(group_defs_all)
    rows = find(strcmp(genotype,group_defs_all(gd).geno) & is_dark==group_defs_all(gd).dark);
    these_flies = unique(fly_num(rows));
    group_fly_list_all{gd}   = these_flies;
    group_fly_trials_all{gd} = arrayfun(@(f) rows(fly_num(rows)==f), these_flies, 'UniformOutput',false);
end

% extend the section 10b PVA cache (built only for opt_trials) to cover
% every trial actually needed by any group here -- f_smooth_opt2 is
% fixed, so a trial already cached (an empty-control closed-loop trial)
% is skipped rather than recomputed.
all_group_trials = find(ismember(genotype,geno_order_all));
for ti = all_group_trials'
    if isempty(pva_mu_cache{ti})
        [pva_mu_cache{ti},pva_rho_cache{ti}] = trial_pva(all_data(ti),f_smooth_opt2);
    end
end

fly_gain_cue_allgroups = [];
fly_gain_fly_allgroups = [];
cat_x_allgroups        = [];
fprintf('\n=== gain, winning pipeline, every genotype x light-condition group ===\n');
for gd = 1:numel(group_defs_all)
    these_flies = group_fly_list_all{gd};
    trial_lists = group_fly_trials_all{gd};
    gc_g = nan(numel(these_flies),1);
    gf_g = nan(numel(these_flies),1);
    for ff = 1:numel(these_flies)
        [gc,gf] = fly_gain_cached(all_data,trial_lists{ff},pva_mu_cache,pva_rho_cache, ...
            mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        gc_g(ff) = gc;
        gf_g(ff) = gf;
    end
    fprintf('  %-24s gain_cue = %.3f +/- %.3f, gain_fly = %.3f +/- %.3f (n=%d flies)\n', ...
        group_defs_all(gd).label, mean(gc_g,'omitnan'), std(gc_g,'omitnan')/sqrt(sum(~isnan(gc_g))), ...
        mean(gf_g,'omitnan'), std(gf_g,'omitnan')/sqrt(sum(~isnan(gf_g))), sum(~isnan(gc_g)));

    fly_gain_cue_allgroups = [fly_gain_cue_allgroups; gc_g]; %#ok<AGROW>
    fly_gain_fly_allgroups = [fly_gain_fly_allgroups; gf_g]; %#ok<AGROW>
    cat_x_allgroups        = [cat_x_allgroups; gd*ones(numel(these_flies),1)]; %#ok<AGROW>
end

cat_labels_all = {group_defs_all.label};
base_colors_all = lines(numel(geno_order_all));
cat_colors_all = zeros(numel(group_defs_all),3);
for gd = 1:numel(group_defs_all)
    gi = find(strcmp(geno_order_all,group_defs_all(gd).geno));
    cat_colors_all(gd,:) = base_colors_all(gi,:);
end

figure(7); clf
set(gcf,'Name','gain, every genotype x light-condition group','Position',[100,100,1200,800])
subplot(2,1,1)
groupplot(cat_x_allgroups,fly_gain_cue_allgroups,cat_labels_all,cat_colors_all)
hold on; plot(xlim,[1,1],':k'); plot(xlim,[0,0],'-','Color',[.85,.85,.85])
ylabel('gain\_cue = slope(bump\_vel ~ 1 + cue\_vel)')
title('bump vs. visual cue (target=1 in closed loop, ~0 expected in dark)')

subplot(2,1,2)
groupplot(cat_x_allgroups,fly_gain_fly_allgroups,cat_labels_all,cat_colors_all)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ 1 + fly\_vel)')
title('bump vs. fly rotation (target=0.8 in closed loop)')

sgtitle(sprintf('winning pipeline applied to every group: im.f=%d frames, bump=%.2fs, cue=%.2fs, fly\\_vel=%.2fs, lag=%d frames', ...
    f_smooth_opt2,mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames),'Interpreter','none')

%% 12) same groups, same winning pipeline, but gain as a THROUGH-ORIGIN fit (no intercept) instead of section 11's affine fit
% section 11's gain_cue/gain_fly come from bump_vel ~ 1 + predictor_vel
% (an intercept absorbs any constant offset between the two, e.g. from a
% residual lag/calibration mismatch, so the SLOPE alone isn't biased by
% it). Forcing the fit through the origin (bump_vel ~ predictor_vel,
% dropping the "1") answers a different question -- "how much bump
% velocity per unit predictor velocity, with no allowance for a constant
% offset" -- and is more sensitive to any such offset, since it can no
% longer be absorbed separately. Reuses the SAME pva cache and group
% definitions as section 11 (fly_gain_cached_noint below is
% fly_gain_cached with the intercept column removed from both \ solves).
fly_gain_cue_allgroups_noint = [];
fly_gain_fly_allgroups_noint = [];
cat_x_allgroups_noint        = [];
fprintf('\n=== gain (through-origin fit, no intercept), winning pipeline, every group ===\n');
for gd = 1:numel(group_defs_all)
    these_flies = group_fly_list_all{gd};
    trial_lists = group_fly_trials_all{gd};
    gc_g = nan(numel(these_flies),1);
    gf_g = nan(numel(these_flies),1);
    for ff = 1:numel(these_flies)
        [gc,gf] = fly_gain_cached_noint(all_data,trial_lists{ff},pva_mu_cache,pva_rho_cache, ...
            mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        gc_g(ff) = gc;
        gf_g(ff) = gf;
    end
    fprintf('  %-24s gain_cue = %.3f +/- %.3f, gain_fly = %.3f +/- %.3f (n=%d flies)\n', ...
        group_defs_all(gd).label, mean(gc_g,'omitnan'), std(gc_g,'omitnan')/sqrt(sum(~isnan(gc_g))), ...
        mean(gf_g,'omitnan'), std(gf_g,'omitnan')/sqrt(sum(~isnan(gf_g))), sum(~isnan(gc_g)));

    fly_gain_cue_allgroups_noint = [fly_gain_cue_allgroups_noint; gc_g]; %#ok<AGROW>
    fly_gain_fly_allgroups_noint = [fly_gain_fly_allgroups_noint; gf_g]; %#ok<AGROW>
    cat_x_allgroups_noint        = [cat_x_allgroups_noint; gd*ones(numel(these_flies),1)]; %#ok<AGROW>
end

figure(8); clf
set(gcf,'Name','gain (through-origin fit), every genotype x light-condition group','Position',[100,100,1200,800])
subplot(2,1,1)
groupplot(cat_x_allgroups_noint,fly_gain_cue_allgroups_noint,cat_labels_all,cat_colors_all)
hold on; plot(xlim,[1,1],':k'); plot(xlim,[0,0],'-','Color',[.85,.85,.85])
ylabel('gain\_cue = slope(bump\_vel ~ cue\_vel)')
title('bump vs. visual cue (target=1 in closed loop, ~0 expected in dark)')

subplot(2,1,2)
groupplot(cat_x_allgroups_noint,fly_gain_fly_allgroups_noint,cat_labels_all,cat_colors_all)
hold on; plot(xlim,[0.8,0.8],':k')
ylabel('gain\_fly = slope(bump\_vel ~ fly\_vel)')
title('bump vs. fly rotation (target=0.8 in closed loop)')

sgtitle(sprintf('THROUGH-ORIGIN fit, winning pipeline: im.f=%d frames, bump=%.2fs, cue=%.2fs, fly\\_vel=%.2fs, lag=%d frames', ...
    f_smooth_opt2,mu_smooth_opt2,cue_smooth_opt2,vel_smooth_opt2,lag_frames),'Interpreter','none')

%% 13) save every figure as a PDF into ugly_figures/rnai -- same convention as lpsp_kir_claude.m's own "save all figures" section
fig_dir = fullfile('ugly_figures','rnai');
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
    out_path  = fullfile(fig_dir,sprintf('gain_scratch_fig%02d_%s.pdf',fig.Number,safe_name));
    try
        exportgraphics(fig,out_path)
        fprintf('saved %s\n', out_path);
    catch ME
        save_failed{end+1} = out_path; %#ok<AGROW>
        warning('could not save %s (%s) -- is it open in another program?', out_path, ME.message);
    end
end
if ~isempty(save_failed)
    fprintf('\n%d figure(s) failed to save -- close them in any viewer and re-run to update:\n', numel(save_failed));
    fprintf('  %s\n', save_failed{:});
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

function x_filled = fill_nan_gaps(x)
    % ft.cue's scattered NaN dropouts (~0.1-0.4% of samples, confirmed in
    % conversation) happen almost exclusively right at the wrapped +/-pi
    % boundary (98.1% of gaps have a bracketing value within 0.3 rad of
    % +/-pi, checked across every trial in the dataset) -- so plain LINEAR
    % interpolation on the raw wrapped signal (the previous version of
    % this function) ignores the wrap and can walk the wrong way around
    % the circle across a gap, producing a spurious jump of up to +/-2*pi
    % (confirmed: 26.7% of all gaps produced a >1 rad linear-interpolated
    % jump despite a true circular distance <0.5 rad; max jump seen was
    % exactly 2*pi) -- which then becomes a large, spurious velocity spike
    % once unwrapped, smoothed, and differentiated downstream (a real
    % contributor to the "ringing" tendrils seen in cue_vel scatter plots).
    % Since the bracketing values already sit at the wrap boundary,
    % setting dropped samples to exactly pi (rather than interpolating)
    % lands them at that same boundary, so unwrap()'s own +/-2*pi
    % correction lines them up with their neighbors correctly instead of
    % introducing a fake jump.
    x_filled = x(:);
    bad = isnan(x_filled);
    x_filled(bad) = pi;
end

function [mu_new,rho_new] = trial_pva(trial,f_smooth_frames)
    % bump position/strength via a plain sample-count moving average
    % (movmean) of raw fluorescence (im.f), z-scored per glomerulus, then
    % population-vector-averaged against each glomerulus's own angular
    % position (im.alpha) -- same construction as
    % lpsp_rnai_claude_v2.m's trial_pva_movmean.
    f_smooth = movmean(trial.im.f,f_smooth_frames,2);
    f_z = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);

    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp] = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function [fly_vel,cue_vel,bump_vel,valid] = trial_gain_vectors(trial,mu_raw,rho_raw,mu_smooth_s,heading_smooth_s, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % all three vectors on the fictrac timebase (ft.xf):
    %   fly_vel  = heading_smooth_s-smoothed r_speed (already a velocity,
    %              smoothed directly -- see header comment on why that's
    %              equivalent to smoothing-then-differentiating a position)
    %   cue_vel  = d/dt of heading_smooth_s-smoothed unwrap(-cue), cue's
    %              scattered NaN dropouts filled first (fill_nan_gaps)
    %   bump_vel = d/dt of mu_smooth_s-smoothed unwrap(mu_raw) (image
    %              timebase), interpolated onto the fictrac timebase
    % valid flags samples where the fly is actually turning fast enough to
    % be informative (vel_thresh<|fly_vel|<vel_max), the bump velocity
    % isn't a degenerate/unwrap-artifact outlier (|bump_vel|<bump_thresh),
    % and the bump is concentrated enough to trust (rho>rho_thresh) --
    % same thresholds/logic as lpsp_kir_claude.m's own trial_gain_vectors,
    % now gating both candidate predictors (cue_vel, fly_vel) at once
    % since they share the same "is the fly turning" mask.
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu_raw);
    xb = linspace(xf(1),xf(end),n_im)';

    dt_im = (xf(end)-xf(1)) / (n_im-1);
    mu_win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed = movmean(unwrap(mu_raw(:)),mu_win_im);
    bump_vel = gradient(interp1(xb,mu_smoothed,xf,'linear','extrap'))/dt;
    rho_i    = interp1(xb,rho_raw(:),xf,'linear','extrap');

    heading_win = max(1,round(heading_smooth_s/dt));
    fly_vel = movmean(trial.ft.r_speed(:),heading_win);

    cue_filled = fill_nan_gaps(trial.ft.cue);
    cue_smoothed = movmean(unwrap(-cue_filled),heading_win);
    cue_vel = gradient(cue_smoothed)/dt;

    valid = abs(fly_vel) > vel_thresh & abs(fly_vel) < vel_max & abs(bump_vel) < bump_thresh & rho_i > rho_thresh;
end

function [gain_cue,gain_fly,n_valid] = fly_gain(all_data,trial_list,f_smooth_frames,mu_smooth_s,heading_smooth_s, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % pools trial_gain_vectors across every trial in trial_list (one fly's
    % own closed-loop trials), then fits bump_vel ~ 1 + predictor_vel
    % (through-intercept, matching lpsp_kir_claude.m's own gain-scatter
    % convention) for both candidate predictors -- so a fly with multiple
    % qualifying trials contributes one combined pair of gain estimates,
    % not one pair per trial.
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        trial = all_data(trial_list(k));
        [mu_raw,rho_raw] = trial_pva(trial,f_smooth_frames);
        [fv,cv,bv,vd] = trial_gain_vectors(trial,mu_raw,rho_raw,mu_smooth_s,heading_smooth_s, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    n_valid = sum(valid);
    if n_valid < 50
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid);
    gain_cue = b_cue(2);
    b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid);
    gain_fly = b_fly(2);
end

function [fly_vel,cue_vel,bump_vel,valid] = trial_gain_vectors_gauss(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % settled pipeline (see section 9's header comment): GAUSSIAN, single-
    % pass smoothing (smoothdata(...,'gaussian',...), not movmean) on both
    % the bump position trace and cue/fly_vel, PLUS a fixed lag shifting
    % cue_vel/fly_vel earlier relative to bump_vel (the bump/calcium
    % signal lags behind behavior) -- otherwise the same construction as
    % trial_gain_vectors above (window sizes converted to samples via each
    % trial's own image/fictrac dt, fly_vel smoothed directly since it's
    % already a velocity, cue's NaN dropouts filled first).
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu_raw);
    xb = linspace(xf(1),xf(end),n_im)';
    dt_im = (xf(end)-xf(1)) / (n_im-1);

    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    bump_vel_full = gradient(interp1(xb,mu_smoothed,xf,'linear','extrap'))/dt;
    rho_full = interp1(xb,rho_raw(:),xf,'linear','extrap');

    win_fictrac = max(1,round(cue_smooth_s/dt));
    fly_vel_full = smoothdata(trial.ft.r_speed(:),'gaussian',win_fictrac);

    cue_filled = fill_nan_gaps(trial.ft.cue);
    cue_smoothed = smoothdata(unwrap(-cue_filled),'gaussian',win_fictrac);
    cue_vel_full = gradient(cue_smoothed)/dt;

    % shift-and-trim: positive lag_frames means cue_vel/fly_vel are EARLIER
    % than bump_vel by that many fictrac frames (the bump follows behavior
    % with a delay), same convention as lpsp_kir_claude.m's own
    % trial_gain_vectors (there applied to shift bump relative to fly_vel
    % directly; here applied to shift cue_vel/fly_vel relative to bump_vel,
    % per this section's own "lag frames to lag cue to match bump" framing).
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

function [gain_cue,gain_fly,n_valid] = fly_gain_gauss(all_data,trial_list,f_smooth_frames,mu_smooth_s,cue_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % same per-fly pooling as fly_gain above, calling
    % trial_gain_vectors_gauss instead of trial_gain_vectors.
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        trial = all_data(trial_list(k));
        [mu_raw,rho_raw] = trial_pva(trial,f_smooth_frames);
        [fv,cv,bv,vd] = trial_gain_vectors_gauss(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    n_valid = sum(valid);
    if n_valid < 50
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid);
    gain_cue = b_cue(2);
    b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid);
    gain_fly = b_fly(2);
end

function [fly_vel,cue_vel,bump_vel,valid] = trial_gain_vectors_gauss2(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % same construction as trial_gain_vectors_gauss, but cue_smooth_s and
    % vel_smooth_s are INDEPENDENT windows (section 9/trial_gain_vectors_gauss
    % used one shared window, cue_smooth_s, for both) -- see section 10's
    % header comment for why. lag_frames shifts cue_vel/fly_vel earlier
    % relative to bump_vel, same convention as trial_gain_vectors_gauss.
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
    cue_filled = fill_nan_gaps(trial.ft.cue);
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

function [gain_cue,gain_fly,n_valid] = fly_gain_gauss2(all_data,trial_list,f_smooth_frames,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % same per-fly pooling as fly_gain_gauss, but re-derives mu_raw/rho_raw
    % from scratch each call (f_smooth_frames is being swept) -- used only
    % by section 10's stage 1; stages 2-4 reuse a precomputed PVA cache via
    % fly_gain_cached instead, since f_smooth is fixed by then.
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        trial = all_data(trial_list(k));
        [mu_raw,rho_raw] = trial_pva(trial,f_smooth_frames);
        [fv,cv,bv,vd] = trial_gain_vectors_gauss2(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    n_valid = sum(valid);
    if n_valid < 50
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid);
    gain_cue = b_cue(2);
    b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid);
    gain_fly = b_fly(2);
end

function [gain_cue,gain_fly,n_valid] = fly_gain_cached(all_data,trial_list,mu_cache,rho_cache,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % same per-fly pooling as fly_gain_gauss2, but pulls mu_raw/rho_raw
    % from a precomputed cache (indexed by trial number, filled once in
    % section 10b) instead of recomputing the PVA every call -- section
    % 10's stages 2-4 sweep mu/cue/vel smoothing only, so the underlying
    % PVA (which only depends on im.f smoothing) never changes and doesn't
    % need to be redone for every candidate.
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        ti = trial_list(k);
        trial = all_data(ti);
        [fv,cv,bv,vd] = trial_gain_vectors_gauss2(trial,mu_cache{ti},rho_cache{ti},mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    % conservative activity floor: this fly's own pooled trials (in this
    % light condition) must show real turning (|fly_vel|>vel_thresh) for
    % at least min_frac_moving of ALL samples, not just the "valid"-gated
    % ones -- confirmed on trial 533 (empty>th, dark, gain_cue=-3.245):
    % that fly was moving above threshold for only 0.25% of a 600s trial
    % (91/36006 samples, confined to two brief blips), so its "gain" was
    % really a regression fit to a handful of noise-dominated samples,
    % not a real tracking relationship. 1% is a deliberately loose floor
    % -- this is about excluding near-total quiescence, not borderline cases.
    min_frac_moving = 0.01;
    n_valid = sum(valid);
    if n_valid < 50 || mean(abs(fly_vel) > vel_thresh) < min_frac_moving
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid);
    gain_cue = b_cue(2);
    b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid);
    gain_fly = b_fly(2);
end

function [gain_cue,gain_fly,n_valid] = fly_gain_cached_noint(all_data,trial_list,mu_cache,rho_cache,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    % identical pooling to fly_gain_cached, but fits bump_vel ~ predictor_vel
    % THROUGH THE ORIGIN (no intercept column) -- see section 12's header
    % comment for why that's a different question than the affine fit.
    % Unlike the affine fit, a through-origin slope is just
    % sum(x.*y)/sum(x.^2) -- with no intercept to absorb it, a predictor
    % whose variance collapses toward zero (confirmed to happen for
    % cue_vel on at least one dark-condition fly, where "valid" only
    % constrains fly_vel/bump_vel/rho, not cue_vel itself) drives sum(x.^2)
    % toward zero and the slope toward +/-Inf -- not a real gain, a
    % divide-by-near-zero artifact. min_predictor_std guards against that.
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        ti = trial_list(k);
        trial = all_data(ti);
        [fv,cv,bv,vd] = trial_gain_vectors_gauss2(trial,mu_cache{ti},rho_cache{ti},mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel  = [fly_vel; fv]; %#ok<AGROW>
        cue_vel  = [cue_vel; cv]; %#ok<AGROW>
        bump_vel = [bump_vel; bv]; %#ok<AGROW>
        valid    = [valid; vd]; %#ok<AGROW>
    end

    % same conservative activity floor as fly_gain_cached -- see that
    % function's comment for the trial-533 case this is guarding against.
    min_frac_moving = 0.01;
    n_valid = sum(valid);
    min_predictor_std = 0.01; % rad/s
    if n_valid < 50 || std(cue_vel(valid)) < min_predictor_std || std(fly_vel(valid)) < min_predictor_std ...
            || mean(abs(fly_vel) > vel_thresh) < min_frac_moving
        gain_cue = nan;
        gain_fly = nan;
        return
    end

    gain_cue = cue_vel(valid) \ bump_vel(valid);
    gain_fly = fly_vel(valid) \ bump_vel(valid);
end

function groupplot(cat_x, values, cat_labels, colors)
    % jittered per-point scatter + mean +/- SEM errorbar per category, n
    % embedded in the xtick label -- same convention as
    % lpsp_rnai_claude_v2.m's own groupplot.
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
