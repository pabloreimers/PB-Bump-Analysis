%% hackathon_claude
% Summarize how many flies and how many trials per fly are in
% hackathon_20250729.mat, for each trial whether the visual cue was in
% closed loop, dark, or a mix of both, and (where the trial has a
% closed-loop portion) what VR gain was applied to the cue.
%
% Fly identity isn't a field on its own -- it has to be parsed out of
% all_data(i).meta (the raw data folder path). Most sessions name a
% per-fly subfolder like "15_pr_5" (fly 15, experimenter initials, 5th fly
% that experimenter ran); a couple of early 20250727 trials predate that
% convention, so those fall back to their session folder name instead.
%
% Closed-loop/dark comes from all_data(i).ft.pattern:
%   - pattern 'berg4': the arena reports a literal cue brightness in
%     all_data(i).ft.cue_brightness (0 = dark, >0 = cue on), and a single
%     berg4 trial can switch between the two mid-trial, so it can come out
%     "mixed". These trials also carry the actual applied gain directly in
%     all_data(i).ft.gain, so we just read it off during the cue-on samples.
%   - any other pattern (a G4 pattern filename): the fly sees either a
%     bar/grating pattern for the whole trial (closed loop) or a uniform
%     'background' pattern (dark) -- there's no per-frame brightness trace
%     or applied-gain field for these, so closed loop/dark comes from the
%     pattern name, and the gain has to be estimated from behavior: regress
%     the cue's angular velocity against the fly's yaw velocity in a
%     sliding window, same approach as the all_data(i).gain.vr calculation
%     in hackathon_script.m.

%% load data
clear all
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(repo_root) % for natsort.m
load(fullfile(repo_root,'.data','hackathon_20250729.mat')) % loads variable "all_data"

%% recompute the bump (mu, rho) from a 5-frame moving-average-smoothed im.f
% im.f is the raw per-cluster fluorescence, one value per imaging frame --
% noisy enough at the single-frame level that it visibly roughens the bump
% trace. Smooth it before recomputing dF/F, z-score, and the population
% vector average, following the same recipe as process_im() in
% hackathon_script.m, just starting one step later (from f rather than
% from raw imgData, since that's all that's saved in this .mat file).
f0_pct = 7; % matches process_im's default in hackathon_script.m

for i = 1:numel(all_data)
    f   = smoothdata(all_data(i).im.f,2,'movmean',5); % smooth along time (dim 2): im.f is [clusters x frames]
    f0  = prctile(f,f0_pct,2);
    dff = (f - f0) ./ f0;
    z   = zscore(dff,[],2);

    alpha    = all_data(i).im.alpha;
    [x,y]    = pol2cart(alpha,z');
    [mu,rho] = cart2pol(mean(x,2),mean(y,2));
    mu = mod(mu,2*pi); mu(mu>pi) = mu(mu>pi)-2*pi;

    all_data(i).im.mu  = mu;
    all_data(i).im.rho = rho;
    all_data(i).im.z   = z;
    all_data(i).im.d   = dff;
end

%% assign each trial to a fly, a closed-loop/dark/mixed category, and a closed-loop gain label
n_trials = numel(all_data);
fly_id   = cell(n_trials,1);
category = nan(n_trials,1); % 0 = closed loop, 1 = dark, 2 = mixed
gain_txt = cell(n_trials,1); % text label for the closed-loop gain, '' if the trial has no closed-loop portion (or no usable estimate)

for i = 1:n_trials
    tok = regexp(all_data(i).meta,'\d+_[a-z]{2}_\d+','match','once');
    if isempty(tok)
        parts = strsplit(all_data(i).meta,filesep);
        tok = parts{end-2}; %fall back to the session folder for trials that predate the fly-naming convention
    end
    fly_id{i} = tok;

    ft = all_data(i).ft;
    pattern = ft.pattern;

    if strcmpi(pattern,'berg4')
        cb = ft.cue_brightness;
        if all(cb==0)
            category(i) = 1;
        elseif all(cb>0)
            category(i) = 0;
        else
            category(i) = 2;
        end
        g = ft.gain(cb>0); % applied gain is known directly -- just read it during the cue-on samples
    elseif contains(pattern,'background')
        category(i) = 1;
        g = [];
    else
        category(i) = 0;
        g = local_measured_gain(ft); % no applied-gain field -- estimate it from behavior
    end

    gain_txt{i} = format_gain(g);
end

flies   = natsort(unique(fly_id));
n_flies = numel(flies);

trials_per_fly = zeros(n_flies,1);
for f = 1:n_flies
    trials_per_fly(f) = sum(strcmp(fly_id,flies{f}));
end
n_cols = max(trials_per_fly);

fprintf('%d flies, %d trials total, %d trials for the fly with the most\n',n_flies,n_trials,n_cols)

%% build an RGB image: rows = flies, columns = trials, color = category
colors = [0.20 0.45 0.85;  % 0 closed loop
          0.10 0.10 0.10;  % 1 dark
          0.90 0.55 0.10]; % 2 mixed
text_colors = [1 1 1; 1 1 1; 0 0 0]; % contrasting text color per category

img    = ones(n_flies,n_cols,3); % white = no trial (fly had fewer trials than the max)
labels = cell(n_flies,n_cols);
label_colors = cell(n_flies,n_cols);
for f = 1:n_flies
    idx = find(strcmp(fly_id,flies{f}));
    for t = 1:numel(idx)
        img(f,t,:) = colors(category(idx(t))+1,:);
        labels{f,t} = gain_txt{idx(t)};
        label_colors{f,t} = text_colors(category(idx(t))+1,:);
    end
end

%% plot
figure('color','w','Position',[100 100 900 800]); clf
image(img)
set(gca,'YTick',1:n_flies,'YTickLabel',flies,'XTick',1:n_cols,'TickLength',[0 0],'TickLabelInterpreter','none')
xlabel('trial number')
ylabel('fly')
title(sprintf('%d flies, %d trials',n_flies,n_trials))
box on

hold on
for f = 1:n_flies
    for t = 1:n_cols
        if ~isempty(labels{f,t})
            text(t,f,labels{f,t},'HorizontalAlignment','center','FontSize',7,'Color',label_colors{f,t})
        end
    end
end

h = gobjects(3,1);
for c = 1:3
    h(c) = patch(nan,nan,colors(c,:));
end
legend(h,{'closed loop','dark','mixed'},'Location','eastoutside')
fontsize(gcf,14,'points')

%% export
exportgraphics(gcf, fullfile(repo_root,'ugly_figures','exports','hackathon_claude.png'), 'Resolution', 300)

%% pick the heading-side (r_speed) smoothing that makes closed-loop gain=0.8 mobility track a 1:1 slope
% Bump mobility (below) regresses bump path length (from im.mu) against
% heading path length (from r_speed) across walking bouts, using a gaussian
% smoothing window on each side. The bump side is already de-noised by the
% 5-frame f smoothing above, so it keeps a fixed, short window; here we
% sweep the heading (r_speed) side's window and pick whichever value pools
% every closed-loop, gain=0.8 walking bout (across all flies/trials/blocks)
% closest to the y=x line -- i.e. the fly's own turning should predict the
% bump 1:1 when gain is close to nominal, and this calibrates the smoothing
% mismatch between the imaging and behavior clocks that would otherwise bias
% that ratio away from 1.

smooth_window_mu   = 60;           %samples (~1s @ ~60Hz), gaussian smoothing window for the bump (im.mu) side -- unchanged
turn_thresh        = .25;          %rad/s, minimum heading speed to call a bout "walking"
max_gap_frames     = round(.5*60); %frames of non-walking allowed within a bout before splitting it
min_walking_frames = round(.5*60); %minimum bout length to keep
min_block_dur      = 1;            %seconds of usable (in-bout) data required to trust a block's mobility estimate

candidate_windows = [6 15 30 60 90 120 180 300 450 600]; % samples, ~0.1s to 10s @ 60Hz
sweep_slope = nan(size(candidate_windows));
sweep_r2    = nan(size(candidate_windows));

for wi = 1:numel(candidate_windows)
    [~,~,~,~,~,~,pooled] = compute_gain_blocks(all_data,fly_id,flies,n_flies,smooth_window_mu,candidate_windows(wi),turn_thresh,max_gap_frames,min_walking_frames,min_block_dur);

    keep = ~pooled.dark & round(pooled.gain,1)==0.8;
    x = pooled.cue(keep); y = pooled.mu(keep); w = sqrt(pooled.dur(keep));

    sweep_slope(wi) = (w.*x) \ (w.*y); % pooled best-fit slope, for reference
    ss_res = sum(w.*(y-x).^2);         % goodness of fit specifically to the y=x (slope=1) line, not to the best-fit slope
    ss_tot = sum(w.*(x-sum(w.*x)/sum(w)).^2);
    sweep_r2(wi) = 1 - ss_res/ss_tot;

    fprintf('heading smoothing = %3d samples (%.2fs): pooled slope=%.2f, R^2 to y=x=%.2f (n=%d bouts)\n', ...
        candidate_windows(wi),candidate_windows(wi)/60,sweep_slope(wi),sweep_r2(wi),sum(keep))
end

[~,best_i] = min(abs(sweep_slope-1));
smooth_window_r = candidate_windows(best_i);
fprintf('-> using %d samples (%.2fs) of heading smoothing (pooled slope closest to 1)\n',smooth_window_r,smooth_window_r/60)

figure('color','w','Position',[100 100 700 450]); clf
yyaxis left
plot(candidate_windows/60,sweep_slope,'-o'); hold on; plot(xlim,[1,1],':k')
ylabel('pooled slope (bump / heading path length)')
yyaxis right
plot(candidate_windows/60,sweep_r2,'-s')
ylabel('R^2 to the y=x line')
xlabel('r_{speed} smoothing window (s)')
title('heading-smoothing sweep, closed-loop gain=0.8 walking bouts')
xline(smooth_window_r/60,':k','chosen')
fontsize(gcf,14,'points')
exportgraphics(gcf, fullfile(repo_root,'ugly_figures','exports','hackathon_claude_smoothing_sweep.png'), 'Resolution', 300)

%% quantify bump mobility for closed-loop (constant-gain) blocks and dark blocks, using the chosen smoothing
% Bump mobility is the path-length gain metric from lpsp_p2x2_walking_script.m:
% detect walking bouts from r_speed, then run a duration-weighted regression of
% bump path length (from im.mu) against heading path length (from r_speed)
% across bouts -> mov_ratio. Here it's computed per BLOCK rather than per
% trial, where a block is a contiguous stretch of constant closed-loop gain,
% or a contiguous dark stretch. A berg4 "ramp" is really a staircase of many
% short constant-gain holds (~7.5s each), so it naturally becomes many
% separate gain conditions instead of being thrown out for "not being
% constant". Each closed-loop block is labeled with its own gain (applied,
% or for non-berg4 patterns, the behaviorally-estimated gain from
% local_measured_gain above); each dark block is labeled with the most
% recently active closed-loop gain for that fly, which can carry over from
% an earlier trial (e.g. a whole-trial "background" dark trial).

[block_dark,block_gain,block_mov,block_trial,block_start,block_end] = compute_gain_blocks( ...
    all_data,fly_id,flies,n_flies,smooth_window_mu,smooth_window_r,turn_thresh,max_gap_frames,min_walking_frames,min_block_dur);

fprintf('%d closed-loop blocks, %d dark blocks with a usable mobility estimate\n',sum(~block_dark),sum(block_dark))

%% group plot: bump mobility vs. closed-loop gain (own gain if closed loop, preceding gain if dark)
gain_group = round(block_gain,1); % bin to the experiment's ~0.1 gain-step resolution
groups     = unique(gain_group);
colors2    = [0.20 0.45 0.85; 0.10 0.10 0.10]; % closed loop, dark
jit        = 0.12;

figure('color','w','Position',[100 100 1000 500]); clf; hold on
for gi = 1:numel(groups)
    for is_dark = 0:1
        idx = gain_group==groups(gi) & block_dark==is_dark;
        if ~any(idx); continue; end
        x = gi + jit*(2*is_dark-1); % closed loop shifted left, dark shifted right
        swarmchart(x*ones(sum(idx),1),block_mov(idx),15,colors2(is_dark+1,:),'filled','MarkerFaceAlpha',.5,'XJitterWidth',jit)
        errorbar(x,mean(block_mov(idx)),std(block_mov(idx))/sqrt(sum(idx)),'o','Color',colors2(is_dark+1,:)*.5,'LineWidth',1.5,'MarkerFaceColor',colors2(is_dark+1,:)*.5)
    end
end
plot(xlim,[1,1],':k')
xticks(1:numel(groups)); xticklabels(compose('%.1f',groups))
xlabel('closed-loop gain (own gain if closed loop, preceding closed-loop gain if dark)')
ylabel('bump path length / heading path length')
title('bump mobility by closed-loop gain')

h = gobjects(2,1);
h(1) = scatter(nan,nan,30,colors2(1,:),'filled');
h(2) = scatter(nan,nan,30,colors2(2,:),'filled');
legend(h,{'closed loop','dark'},'Location','best')
fontsize(gcf,14,'points')

exportgraphics(gcf, fullfile(repo_root,'ugly_figures','exports','hackathon_claude_mobility.png'), 'Resolution', 300)

%% example plots: imagesc + bump + heading overlay, for gain 0.8 & 1.6, closed loop vs dark
% For each of the 4 (gain x lighting) combinations, pick the block (from the
% same block list behind the group plot above) with the most walking data and
% re-plot its raw imaging heatmap (im.z) with the bump (mu, white) and the
% fly's own heading (integrated from r_speed, zeroed at the start of the
% window) overlaid -- same style as the example figure in hackathon_script.m.
% Heading comes from r_speed rather than the visual cue so that it means the
% same thing in both closed loop and dark (the cue is gain-scaled and only
% actually seen by the fly when the trial is in closed loop).

example_gains   = [0.8 1.6];
lighting_labels = {'closed loop','dark'};
line_colors     = {'c','m'}; % closed loop = cyan, dark = magenta, matching hackathon_script.m's convention

figure('color','w','Position',[100 100 1400 900]); clf
for gi = 1:numel(example_gains)
    for is_dark = 0:1
        idx = find(gain_group==example_gains(gi) & block_dark==is_dark);
        if isempty(idx)
            continue
        end
        [~,best] = max(block_end(idx)-block_start(idx)); % longest example of this condition
        bi = idx(best);

        i  = block_trial(bi);
        s  = block_start(bi);
        e  = block_end(bi);
        ft = all_data(i).ft;
        im = all_data(i).im;

        t0 = ft.xf(s); t1 = ft.xf(e);
        xb_idx = find(ft.xb>=t0 & ft.xb<=t1);
        xf_idx = find(ft.xf>=t0 & ft.xf<=t1);

        heading = cumsum(ft.r_speed(xf_idx))/60; % physical heading, zeroed at the start of this window (not the gain-scaled visual cue)
        heading = mod(heading,2*pi); heading(heading>pi) = heading(heading>pi)-2*pi;

        subplot(2,2,(gi-1)*2+is_dark+1)
        imagesc(ft.xb(xb_idx),unwrap(im.alpha),im.z(:,xb_idx))
        hold on
        hh = plot(ft.xf(xf_idx),heading,line_colors{is_dark+1},'LineWidth',1);
        hh.YData(abs(diff(hh.YData))>pi) = nan;
        hm = plot(ft.xb(xb_idx),im.mu(xb_idx),'w','LineWidth',1);
        hm.YData(abs(diff(hm.YData))>pi) = nan;

        xlabel('time (s)'); ylabel('PB angle (rad)')
        title(sprintf('gain %.1f, %s\nfly %s, trial %i (%.0fs window)',example_gains(gi),lighting_labels{is_dark+1},fly_id{i},i,t1-t0),'Interpreter','none')
    end
end
sgtitle('example bump (white) vs heading (color) traces, by closed-loop gain condition')
fontsize(gcf,12,'points')

exportgraphics(gcf, fullfile(repo_root,'ugly_figures','exports','hackathon_claude_examples.png'), 'Resolution', 300)

%% functions
function g = local_measured_gain(ft)
    % Estimate the closed-loop VR gain from behavior: regress the cue's
    % angular velocity against the fly's yaw velocity in a sliding window.
    % There's no ground-truth gain field for these patterns, so this is
    % the best available substitute.
    r  = ft.r_speed(1:end-2);
    c  = -gradient(ft.cue(3:end)) * 60;
    c(abs(c) > 50) = nan;
    xf = ft.xf(1:end-2);

    t = 1:round(max(ft.xf));
    g = nan(size(t));
    for j = 1:length(t)
        idx = xf > j-1 & xf < j+5 & abs(r) > .5 & abs(c) > .5;
        if sum(idx) > 0
            g(j) = r(idx) \ c(idx);
        end
    end
    g = g(abs(g) < 5); % drop unstable few-sample estimates
end

function s = format_gain(g)
    g = g(~isnan(g));
    if isempty(g)
        s = '';
        return
    end
    g_lo = min(g);
    g_hi = max(g);
    if g_hi - g_lo < 0.05 % effectively constant across the closed-loop portion
        s = sprintf('%.1f',median(g));
    else
        s = sprintf('%.1f-%.1f',g_lo,g_hi);
    end
end

function [block_dark,block_gain,block_mov,block_trial,block_start,block_end,pooled] = compute_gain_blocks( ...
        all_data,fly_id,flies,n_flies,smooth_window_mu,smooth_window_r,turn_thresh,max_gap_frames,min_walking_frames,min_block_dur)
    % Segment every trial into closed-loop (constant-gain) / dark blocks and
    % compute the path-length bump-mobility ratio for each, clipping walking
    % bouts to block boundaries. Also returns "pooled", the raw per-bout
    % (mov_cue,mov_mu,dur) triples with their block's gain/dark label
    % attached, so a caller can pool bouts across blocks (e.g. to fit a
    % single slope for "all closed-loop gain=0.8 bouts") rather than being
    % limited to one regression per block.
    block_dark  = [];
    block_gain  = [];
    block_mov   = [];
    block_trial = [];
    block_start = [];
    block_end   = [];

    pooled_mu   = [];
    pooled_cue  = [];
    pooled_dur  = [];
    pooled_gain = [];
    pooled_dark = [];

    for f = 1:n_flies
        trial_idx = find(strcmp(fly_id,flies{f}));
        last_gain = nan; % most recently active closed-loop gain for this fly, carried across trials

        for i = trial_idx'
            ft = all_data(i).ft;
            xf = ft.xf(:);
            dt = median(diff(xf));

            % walking-bout detection over the whole trial, exactly as in lpsp_p2x2_walking_script.m
            mu             = interp1(ft.xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');
            r_speed_smooth = smoothdata(ft.r_speed(:),'gaussian',smooth_window_r);
            mu_smooth      = smoothdata(mu,'gaussian',smooth_window_mu);

            is_walking = abs(r_speed_smooth) > turn_thresh;
            is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
            is_walking = bwareaopen(is_walking,min_walking_frames);

            d = diff([false;is_walking;false]);
            bout_starts = find(d==1);
            bout_ends   = find(d==-1)-1;

            % split this trial into closed-loop (constant-gain) / dark blocks
            pattern = ft.pattern;
            if strcmpi(pattern,'berg4')
                cb    = ft.cue_brightness(:);
                state = round(ft.gain(:),2);
                state(cb==0) = -inf; %flag dark samples distinctly from any real gain value
                [run_starts,run_ends,run_vals] = find_runs(state);
                run_isdark = isinf(run_vals);
            elseif contains(pattern,'background')
                run_starts = 1; run_ends = numel(xf); run_vals = nan; run_isdark = true;
            else
                run_starts = 1; run_ends = numel(xf); run_vals = median(local_measured_gain(ft)); run_isdark = false;
            end

            for k = 1:numel(run_starts)
                s = run_starts(k); e = run_ends(k);

                % clip this trial's walking bouts to the block, and rerun the same weighted regression on the clipped bouts
                mov_mu = []; mov_cue = []; dur = [];
                for b = 1:numel(bout_starts)
                    cs = max(bout_starts(b),s);
                    ce = min(bout_ends(b),e);
                    if ce > cs
                        mov_mu(end+1,1)  = sum(abs(diff(mu_smooth(cs:ce))),'omitnan');
                        mov_cue(end+1,1) = sum(abs(r_speed_smooth(cs:ce)),'omitnan')*dt;
                        dur(end+1,1)     = (ce-cs+1)*dt;
                    end
                end

                if sum(dur) < min_block_dur
                    mov_ratio = nan;
                else
                    w = sqrt(dur);
                    mov_ratio = (w.*mov_cue) \ (w.*mov_mu);
                end

                if run_isdark(k)
                    label_gain    = last_gain;
                    is_dark_block = true;
                else
                    if ~isnan(run_vals(k))
                        last_gain = run_vals(k); %update regardless of whether this block itself has a usable mobility estimate
                    end
                    label_gain    = run_vals(k);
                    is_dark_block = false;
                end

                if ~isnan(label_gain) && ~isnan(mov_ratio)
                    block_dark(end+1,1)  = is_dark_block;
                    block_gain(end+1,1)  = label_gain;
                    block_mov(end+1,1)   = mov_ratio;
                    block_trial(end+1,1) = i;
                    block_start(end+1,1) = s;
                    block_end(end+1,1)   = e;

                    n_b = numel(mov_mu);
                    pooled_mu   = [pooled_mu;   mov_mu];   %#ok<AGROW>
                    pooled_cue  = [pooled_cue;  mov_cue];  %#ok<AGROW>
                    pooled_dur  = [pooled_dur;  dur];      %#ok<AGROW>
                    pooled_gain = [pooled_gain; repmat(label_gain,n_b,1)];    %#ok<AGROW>
                    pooled_dark = [pooled_dark; repmat(is_dark_block,n_b,1)]; %#ok<AGROW>
                end
            end
        end
    end

    pooled = struct('mu',pooled_mu,'cue',pooled_cue,'dur',pooled_dur,'gain',pooled_gain,'dark',pooled_dark);
end

function [run_starts,run_ends,run_vals] = find_runs(x)
    % Split x into maximal runs of consecutive equal values (NaN-safe: two
    % consecutive NaNs are NOT considered equal, so NaN samples each end up
    % as their own singleton run rather than being silently merged).
    x = x(:);
    n = numel(x);
    is_new = [true; x(2:end)~=x(1:end-1) | isnan(x(2:end)) | isnan(x(1:end-1))];
    run_starts = find(is_new);
    run_ends   = [run_starts(2:end)-1; n];
    run_vals   = x(run_starts);
end
