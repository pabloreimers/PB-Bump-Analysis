%% fig_example_flies_pb_heatmaps
% One example fly per genotype from the Figure 3 bump-mobility dataset
% (.data/lpsp_p2x2_walking_20260728_post0708.mat). Per fly, one row:
%   col 1: mean PB image from imgData_reg (frame-average of the registered
%          movie that im.mu / im.z were extracted from -- confirmed in
%          data scripts/lpsp_p2x2_script_minimal.m lines 71-82: it loads
%          imgData_reg.mat and calls process_im(imgData_reg,...); the
%          imgData_smooth_reg / imgData_sum variants are NOT what was used)
%          with the ROI mask outline and the 32 glomerulus centroids.
%   col 2-4: z-scored PB activity (im.z, 32 wedges x time) for one
%          closed-loop trial, one dark trial before the first ATP pulse
%          trial, one dark trial after it. Bump position (im.mu) overlaid
%          in white; cue (-ft.cue) in cyan for the closed-loop trial.
%   col 5: bump mobility on the pre vs post dark trial: per walking bout,
%          bump path length vs heading path length (same bout definition
%          and sqrt(duration)-weighted through-origin slope as Figure 3 /
%          lpsp_p2x2_walking_script.m sections 2+4).
%
% EXAMPLE SELECTION (section "score flies" prints the full table):
%   tracking quality = circular correlation between -ft.cue and im.mu on
%   each pre-ATP closed-loop trial, restricted to samples with rho > 0.2
%   AND cue angular speed > 0.2 rad/s (so frozen-cue stretches, where the
%   bump legitimately drifts, don't count against the fly), averaged over
%   the fly's pre-ATP closed-loop trials. Mobility = Figure 3's pooled
%   dark pre / post slope.
%   * lpsp>p2x2: 20260727 fly 4 -- tracking r = 0.88 on both of its first
%     two closed-loop trials (30, 31; the fly-level mean is only 0.29
%     because its later pre-ATP closed-loop trials 33 and 35 scored -0.66
%     and 0.06 -- the heatmap column uses trial 30), pre-dark mobility ~1
%     (trials 32, 34: 1.31, 0.88), post-dark 3.82 (trial 40). The biggest
%     post-ATP increase in the dataset among flies with good tracking on
%     at least some closed-loop trials. NOTE: in trial 40 the fly barely
%     turned (heading path per bout mostly < 3 rad), so the 3.82 slope
%     rests on small bouts -- see the col 5 scatter.
%   * empty>p2x2: 20260728 fly 12 -- best tracking in the dataset
%     (r = 0.80, 0.94), pre-dark 0.93 / 0.86 (trials 97, 99), post-dark
%     1.11 (trial 103), i.e. essentially unchanged.
%   Runners up: lpsp 20260708 fly 1 (r 0.18/0.73, 0.88 -> 2.20), empty
%   20260728 fly 17 (r 0.86/0.87, but post-dark rose to 1.61).

% !!! IMPORTANT CAVEAT -- RHO AND MOBILITY ARE COUPLED (found 2026-10-14) !!!
% In lpsp_p2x2_walking_20260728_post0708, bump vector strength (im.rho) and
% the bump-mobility slope move together: across all dark trials with >=10
% bouts, Spearman r(trial-mean rho, mobility) = -0.40 (n=73, p<0.001), and
% within just the lpsp post-ATP dark trials r = -0.59 (n=11, p=0.06). The
% lpsp group also shows a post-ATP DROP in dark-trial rho (fly-level mean
% 0.437 -> 0.382, -12%, paired signrank p=0.09, n=14 flies) that the empty
% group does not (0.454 -> 0.478, p=0.77, n=11). So the post-ATP mobility
% increase comes with a weaker/broader bump. This is NOT necessarily an
% artifact -- a destabilized EPG representation would look exactly like
% this -- and the lpsp post-dark trials sit above pre-ATP mobility even at
% normal rho (median 2.5 for rho>=0.35 trials vs 1.2 for all pre-dark).
% But the two measures are not independent: path length is a total-
% variation metric that inflates mechanically when mu is noisy, and mu
% (one 60-sample Gaussian in walk_trial_bouts) is smoothed less than
% r_speed (30 in process_ft + 60 here). REPORT the rho change alongside
% any mobility result, and consider a rho-gated or rho-matched version as
% a control. Quantified in the Claude session that produced
% ugly_figures/scripts/fig_example_flies_pb_heatmaps.m.

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root,'circ_stats'))
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
if ~isfolder(export_dir); mkdir(export_dir); end

%% load
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% label trials
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id);
n_flies = numel(fly_list);

%% bouts per trial (Figure 3 settings)
bout_mu = cell(n,1); bout_cue = bout_mu; bout_dur = bout_mu;
for i = 1:n
    [bout_mu{i}, bout_cue{i}, bout_dur{i}] = walk_trial_bouts(all_data(i), 60, .25, 30, 30);
end

%% score flies: tracking quality (pre-ATP closed loop) + dark mobility pre/post
rho_thr = .2; cue_speed_thr = .2; min_bouts = 15;
fly_track = nan(n_flies,1); fly_pre = fly_track; fly_post = fly_track; fly_lpsp = false(n_flies,1);
fly_trials = cell(n_flies,1); % chronological trial indices
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o); fly_trials{f} = idx;
    fly_lpsp(f) = any(is_lpsp(idx));
    p = find(has_pulse(idx),1); if isempty(p); continue; end
    pre = idx(1:p-1); post = idx(p+1:end); post = post(~has_pulse(post));
    r = nan(numel(pre),1);
    for k = 1:numel(pre)
        if is_dark(pre(k)); continue; end
        r(k) = trial_tracking_r(all_data(pre(k)), rho_thr, cue_speed_thr);
    end
    fly_track(f) = mean(r,'omitnan');
    fly_pre(f)  = pooled_slope(bout_cue(pre(is_dark(pre))),   bout_mu(pre(is_dark(pre))),   bout_dur(pre(is_dark(pre))),   min_bouts);
    fly_post(f) = pooled_slope(bout_cue(post(is_dark(post))), bout_mu(post(is_dark(post))), bout_dur(post(is_dark(post))), min_bouts);
end
[~,o] = sort(fly_track,'descend');
fprintf('%-32s %-5s %-7s %-6s %-6s %-6s\n','fly','geno','trackR','pre','post','delta');
for f = o'
    if isnan(fly_track(f)); continue; end
    g = 'empty'; if fly_lpsp(f); g = 'lpsp'; end
    fprintf('%-32s %-5s %-7.2f %-6.2f %-6.2f %-6.2f\n', strrep(fly_list{f},'Z:/pablo/lpsp_p2x2_walking/',''), g, fly_track(f), fly_pre(f), fly_post(f), fly_post(f)-fly_pre(f));
end

%% example picks (see header) -- fly folder tag + trial numbers
clear ex % a stale ex (e.g. with snip_t0 set by an earlier run) would otherwise survive the field-by-field assignment below
ex(1).tag = '20260727/fly 4';  ex(1).label = 'lpsp>p2x2';  ex(1).cl = 30; ex(1).dark_pre = 34; ex(1).dark_post = 40;
ex(2).tag = '20260728/fly 12'; ex(2).label = 'empty>p2x2'; ex(2).cl = 98; ex(2).dark_pre = 99; ex(2).dark_post = 103;
for e = 1:2
    f = find(contains(fly_list, ex(e).tag)); assert(isscalar(f), 'fly tag %s not unique/found', ex(e).tag)
    ex(e).f = f;
    ex(e).i_cl   = fly_trials{f}(trial_num(fly_trials{f})==ex(e).cl);
    ex(e).i_pre  = fly_trials{f}(trial_num(fly_trials{f})==ex(e).dark_pre);
    ex(e).i_post = fly_trials{f}(trial_num(fly_trials{f})==ex(e).dark_post);
    assert(~is_dark(ex(e).i_cl) && is_dark(ex(e).i_pre) && is_dark(ex(e).i_post), 'lighting mismatch for %s', ex(e).tag)
    % mean PB image from the registered movie the bump was extracted from
    reg_file = fullfile(all_data(ex(e).i_cl).meta, 'imgData_reg.mat');
    fprintf('loading %s ...\n', reg_file)
    tmp = load(reg_file, 'imgData_reg');
    ex(e).mean_img = mean(tmp.imgData_reg, 3); clear tmp
    fprintf('  %s: tracking r = %.2f, dark mobility pre %.2f -> post %.2f\n', ex(e).tag, fly_track(f), fly_pre(f), fly_post(f))
end

%% pick a representative 30 s snippet per trial
% "representative" = among the 30 s windows with the highest walking-bout
% coverage (fraction of samples inside a detected walking bout, same bout
% definition as the mobility metric; windows within 90% of the best
% coverage count), the one with the MEDIAN mean rho -- so the snippet shows
% the bump during active walking without cherry-picking its cleanest or
% worst stretch. Override per trial via ex(e).snip_t0 = [t_cl, t_pre, t_post].
snip_len = 30;
for e = 1:2
    trial_ii = [ex(e).i_cl, ex(e).i_pre, ex(e).i_post];
    if ~isfield(ex,'snip_t0') || isempty(ex(e).snip_t0)
        ex(e).snip_t0 = nan(1,3);
        for c = 1:3
            i = trial_ii(c); xf = all_data(i).ft.xf; xb = all_data(i).ft.xb; dt = median(diff(xf));
            w = round(snip_len/dt);
            is_walking = walk_trial_is_walking(all_data(i), 60, .25, 30, 30);
            cover = movmean(double(is_walking), w, 'Endpoints','discard'); % bout coverage of each w-sample window starting at k
            rho_f = interp1(xb, all_data(i).im.rho, xf, 'linear', 'extrap');
            rho_w = movmean(rho_f, w, 'Endpoints','discard');
            cand = find(cover >= .9*max(cover));
            [~,k] = min(abs(rho_w(cand) - median(rho_w(cand)))); k = cand(k);
            ex(e).snip_t0(c) = xf(k); % window = [xf(k), xf(k)+snip_len]
        end
    end
    fprintf('%s snippet starts (s): CL %.0f, dark pre %.0f, dark post %.0f\n', ex(e).tag, ex(e).snip_t0);
end

%% figure
% layout: 4 rows x 5 cols. Each fly takes two rows: heatmap snippets on
% the upper row, heading + bump traces (not overlaid on the heatmap) on
% the lower row. Mean image and mobility scatter span both rows.
fig = figure(32); clf
set(fig,'Color','w','Units','pixels','Position',[40 40 1900 900])
tl = tiledlayout(4,5,'TileSpacing','compact','Padding','compact');
heat_clim = [-1.5 3];
for e = 1:2
    r0 = (e-1)*2; % first row index (0-based) for this fly
    % --- col 1: mean PB image + mask outline + centroids (spans 2 rows)
    ax = nexttile(tl, r0*5+1, [2 1]);
    im = ex(e).mean_img; lo = prctile(im(:),1); hi = prctile(im(:),99.5);
    imagesc(ax, im, [lo hi]); colormap(ax, bone); axis(ax,'image'); hold(ax,'on')
    set(ax,'XTick',[],'YTick',[])
    msk = all_data(ex(e).i_cl).im.mask;
    B = bwboundaries(msk); for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), 'y-', 'LineWidth', .7); end
    c = all_data(ex(e).i_cl).im.centroids; plot(ax, c(:,2), c(:,1), 'r.', 'MarkerSize', 7)
    title(ax, sprintf('%s  (%s)\nmean imgData_reg, trial %d', ex(e).label, strrep(ex(e).tag,'/',' '), ex(e).cl), 'FontSize', 9, 'Interpreter','none')

    % --- cols 2-4: 30 s z-score heatmap snippet (upper row) + traces (lower row)
    trial_ii = [ex(e).i_cl, ex(e).i_pre, ex(e).i_post];
    ttl = {sprintf('closed loop, trial %d', ex(e).cl), sprintf('dark PRE atp, trial %d', ex(e).dark_pre), sprintf('dark POST atp, trial %d', ex(e).dark_post)};
    for c = 1:3
        i = trial_ii(c); t0 = ex(e).snip_t0(c); t1 = t0 + snip_len;
        xb = all_data(i).ft.xb; xf = all_data(i).ft.xf; dt = median(diff(xf));
        kb = xb >= t0 & xb <= t1; kf = xf >= t0 & xf <= t1;

        % heatmap
        ax = nexttile(tl, r0*5+1+c);
        imagesc(ax, xb(kb)-t0, 1:32, all_data(i).im.z(:,kb), heat_clim); colormap(ax, parula); hold(ax,'on')
        yline(ax, 16.5, 'w:');
        set(ax,'YDir','normal','YTick',[1 16 17 32],'XTickLabel',[],'FontSize',8)
        if c == 1, ylabel(ax,'PB glomerulus'); end
        nb = numel(bout_mu{i}); sl = nan; if nb >= 5, sl = pooled_slope(bout_cue(i), bout_mu(i), bout_dur(i), 5); end
        title(ax, sprintf('%s  (mobility %.2f, %d bouts)\n%d s snippet from t = %.0f s, mean rho %.2f', ttl{c}, sl, nb, snip_len, t0, mean(all_data(i).im.rho(kb))), 'FontSize', 8)
        xlim(ax,[0 snip_len])

        % traces: bump (mu, unwrapped) and heading, both re-zeroed at the snippet start
        ax = nexttile(tl, (r0+1)*5+1+c); hold(ax,'on')
        mu_u = unwrap(all_data(i).im.mu(kb)); mu_u = mu_u - mu_u(1);
        if is_dark(i)
            head = cumsum(all_data(i).ft.r_speed(kf))*dt; head_lab = 'heading (\int r\_speed)';
        else
            head = unwrap(-all_data(i).ft.cue(kf)); head = head - head(1); head_lab = 'heading (-cue)';
        end
        plot(ax, xf(kf)-t0, head, '-', 'Color', [0 .6 .8], 'LineWidth', 1.2)
        plot(ax, xb(kb)-t0, mu_u, 'k-', 'LineWidth', 1.2)
        xlim(ax,[0 snip_len]); xlabel(ax,'time in snippet (s)')
        yl = max([abs(head); abs(mu_u); pi]); ylim(ax, [-yl yl]*1.05)
        step = pi*max(1, ceil(yl/(2*pi))); kk = (-floor(yl/step):floor(yl/step)); % pi-multiple ticks, coarser when the range is large
        yticks(ax, kk*step); yticklabels(ax, arrayfun(@(k) pi_label(k*step/pi), kk, 'UniformOutput', false))
        if c == 1, ylabel(ax,'angle (rad, re-zeroed)'); end
        legend(ax, {head_lab, 'bump (im.mu)'}, 'Location','best','FontSize',7,'Box','off')
        set(ax,'FontSize',8)
    end

    % --- col 5: bump mobility, pre vs post dark trial (spans 2 rows)
    ax = nexttile(tl, r0*5+5, [2 1]); hold(ax,'on')
    cols = {[.3 .3 .3], [.85 .2 .2]}; lab = {'dark PRE','dark POST'};
    for c = 1:2
        i = trial_ii(c+1);
        scatter(ax, bout_cue{i}, bout_mu{i}, 20 + 4*sqrt(bout_dur{i}), cols{c}, 'filled', 'MarkerFaceAlpha', .6)
        sl = pooled_slope(bout_cue(i), bout_mu(i), bout_dur(i), 5);
        xl = [0 max([bout_cue{i}; 1])]; plot(ax, xl, sl*xl, '-', 'Color', cols{c}, 'LineWidth', 1.5)
        lab{c} = sprintf('%s (trial %d): slope %.2f', lab{c}, trial_num(i), sl);
    end
    plot(ax, [0 max(xlim(ax))], [0 max(xlim(ax))], 'k:')
    axis(ax,'square'); xlabel(ax,'heading path length per bout (rad)'); ylabel(ax,'bump path length per bout (rad)')
    legend(ax, lab, 'Location','northwest','FontSize',7,'Box','off')
    title(ax, sprintf('bump mobility, dark trials\n(fly pooled: pre %.2f, post %.2f)', fly_pre(ex(e).f), fly_post(ex(e).f)), 'FontSize', 9)
end
cb = colorbar(nexttile(tl,4)); cb.Label.String = 'z-scored dF/F'; cb.Layout.Tile = 'east';
title(tl, sprintf('example flies: PB anatomy, %d s z-scored activity snippets with heading + bump traces below, and dark-trial bump mobility before vs after ATP', snip_len), 'FontSize', 11)

%% export
exportgraphics(fig, fullfile(export_dir,'fig_example_flies_pb_heatmaps.png'), 'Resolution', 150)
fprintf('exported to %s\n', fullfile(export_dir,'fig_example_flies_pb_heatmaps.png'))

%% helpers
function lab = pi_label(m)
    % '-2\pi', '0', '\pi' etc. for an integer multiple m of pi
    if m == 0, lab = '0'; elseif m == 1, lab = '\pi'; elseif m == -1, lab = '-\pi'; else, lab = [num2str(m) '\pi']; end
end

function is_walking = walk_trial_is_walking(trial, smooth_window, turn_thresh, max_gap_frames, min_walking_frames)
    % the walking-bout mask from walk_trial_bouts, on the behavior clock
    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    is_walking = abs(r_speed_smooth) > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);
end

function r = trial_tracking_r(trial, rho_thr, cue_speed_thr)
    % circ corr between cue and mu, only on high-rho samples while the cue is moving
    xb = trial.ft.xb; xf = trial.ft.xf;
    cue_u = interp1(xf, unwrap(-trial.ft.cue), xb, 'linear', 'extrap');
    cue_sp = abs(gradient(smoothdata(cue_u,'movmean',5), xb));
    good = trial.im.rho > rho_thr & cue_sp > cue_speed_thr;
    r = nan;
    if sum(good) > 100
        r = circ_corrcc(mod(cue_u(good),2*pi), mod(trial.im.mu(good),2*pi));
    end
end

function s = pooled_slope(cue_c, mu_c, dur_c, min_bouts)
    % sqrt(duration)-weighted through-origin slope of bump vs heading path length over pooled bouts
    c = vertcat(cue_c{:}); m = vertcat(mu_c{:}); d = vertcat(dur_c{:});
    s = nan;
    if numel(c) >= min_bouts, w = sqrt(d); s = (w.*c) \ (w.*m); end
end

function [fid, tnum] = walk_meta_fly_trial(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    if ~(isfield(trial,'atp') && isfield(trial.atp,'d') && ~isempty(trial.atp.d)); tf = false; return; end
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    if isempty(onsets); tf = false; return; end
    xb = trial.ft.xb; xf = trial.ft.xf;
    t_onset = xf(onsets);
    atp_sig = sum(trial.atp.d,1);
    t_common = linspace(-peri_win,peri_win,101);
    atp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets)
        atp_aligned(s,:) = interp1(xb - t_onset(s),atp_sig,t_common,'linear',nan);
    end
    atp_peri = mean(atp_aligned,1,'omitnan');
    base_idx = t_common>=base_win(1) & t_common<base_win(2);
    peak_idx = t_common>peak_win(1) & t_common<=peak_win(2);
    tf = (median(atp_peri(peak_idx),'omitnan')-median(atp_peri(base_idx),'omitnan')) > peak_factor*std(atp_peri(base_idx),'omitnan');
end

function [mov_mu, mov_cue, dur] = walk_trial_bouts(trial, smooth_window, turn_thresh, max_gap_frames, min_walking_frames)
    % ported verbatim from Figure_3_claude.m / lpsp_p2x2_walking_script.m section 2
    xf = trial.ft.xf;
    dt = median(diff(xf));
    mu = interp1(trial.ft.xb,unwrap(trial.im.mu),xf,'linear','extrap');
    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);
    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);
    d = diff([false;is_walking;false]);
    bout_starts = find(d==1); bout_ends = find(d==-1)-1;
    mov_mu = nan(numel(bout_starts),1); mov_cue = mov_mu; dur = mov_mu;
    for b = 1:numel(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(r_speed_smooth(bout_starts(b):bout_ends(b))),'omitnan')*dt;
        dur(b)     = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end
