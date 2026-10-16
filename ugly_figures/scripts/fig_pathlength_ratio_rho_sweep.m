%% fig_pathlength_ratio_rho_sweep
% Whole-trial bump mobility as a PATH-LENGTH RATIO, with a parametric sweep
% of the rho threshold used to discount poorly-estimated bump samples.
% Built for one example fly first (the empty>p2x2 control fly picked in
% fig_example_flies_pb_heatmaps.m, 20260728 fly 12 -- best closed-loop
% tracking in the dataset, dark mobility 0.90 -> 1.11 by Figure 3's metric).
%
% WHAT IS COMPUTED, per trial and per rho threshold:
%   * Everything lives on the behavior clock (ft.xf, ~60 Hz). The bump
%     position im.mu is unwrapped and linearly interpolated from the imaging
%     clock (ft.xb, ~10.5 Hz) onto xf, and im.rho likewise.
%   * Smoothing: SAME as Figure 3's mobility metric (walk_trial_bouts) so the
%     numbers are comparable -- a 60-sample (~1 s) Gaussian on the unwrapped,
%     upsampled mu and on ft.r_speed. (Figure 3 additionally inherits the
%     dataset's own smoothing: 10-frame Gaussian on the movie, 30-sample
%     Gaussian on r_speed -- see the smoothing note in the walking script.)
%   * fly path length  = sum over kept samples of |r_speed_smooth| * dt
%     (rad of heading rotation, integrated yaw -- the same "heading" Figure 3
%     uses; cue is not used so dark trials are handled identically)
%   * bump path length = sum of |diff(mu_smooth)| over pairs of CONSECUTIVE
%     kept samples (a gap in the mask breaks the path: we do not count the
%     jump across a discarded stretch as bump movement)
%   * ratio = bump path / fly path on the SAME kept samples. Masking applies
%     to both signals: if a sample is dropped for low rho its fly rotation is
%     dropped too, otherwise the ratio would fall just because the mask
%     removes bump movement but keeps fly movement.
%   * kept samples = rho_f > threshold, for threshold = 0, 0.05, ..., 0.70.
%     threshold 0 keeps everything (the raw whole-trial ratio).
%   * Unlike Figure 3's metric there is NO walking-bout segmentation: the
%     whole trial counts, standing included. So bump drift while the fly
%     stands still adds to the numerator with nothing in the denominator.
%     For context the same sweep is also run on walking samples only (same
%     bout definition as Figure 3, |r_speed_smooth| > 0.25 rad/s with gap
%     closing / minimum duration), shown as a second set of lines.
%
% Trials are labeled pre / PULSE / post relative to the fly's first
% detected ATP-pulse trial (same detector as Figure 3), and closed loop vs
% dark from ft.pattern.

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root,'circ_stats'))
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% parameters
example_tag   = '20260728/fly 12';   % empty>p2x2 control example (fig_example_flies_pb_heatmaps.m)
rho_thr_list  = 0:0.05:0.70;
smooth_window = 60;                   % samples on xf, as in Figure 3's walk_trial_bouts
turn_thresh   = .25; max_gap_frames = 30; min_walking_frames = 30; % walking mask, as in Figure 3

%% label trials of this fly
fly_id = cell(n,1); trial_num = nan(n,1);
for i = 1:n, [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta); end
idx = find(contains(fly_id, example_tag)); [~,o] = sort(trial_num(idx)); idx = idx(o);
n_tr = numel(idx);
is_dark   = arrayfun(@(i) contains(all_data(i).ft.pattern,'background'), idx);
has_pulse = arrayfun(@(i) walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3), idx);
p1 = find(has_pulse,1);
phase = repmat({'pre'},n_tr,1); phase(has_pulse) = {'PULSE'}; phase(p1+1:end) = {'post'}; phase(has_pulse) = {'PULSE'};
geno = 'empty>p2x2'; if contains(all_data(idx(1)).meta,'_lpsp_'), geno = 'lpsp>p2x2'; end
fprintf('%s (%s): %d trials\n', example_tag, geno, n_tr);

%% sweep
n_thr = numel(rho_thr_list);
S = struct('fly_path',nan(n_tr,n_thr),'bump_path',nan(n_tr,n_thr),'ratio',nan(n_tr,n_thr),'frac_kept',nan(n_tr,n_thr), ...
           'fly_path_w',nan(n_tr,n_thr),'bump_path_w',nan(n_tr,n_thr),'ratio_w',nan(n_tr,n_thr),'frac_kept_w',nan(n_tr,n_thr));
dur_s = nan(n_tr,1); frac_walking = nan(n_tr,1); rho_all = cell(n_tr,1); mean_rho = nan(n_tr,1);
for k = 1:n_tr
    tr = all_data(idx(k));
    xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf));
    mu_f  = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    rho_f = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    mu_s  = smoothdata(mu_f, 'gaussian', smooth_window);
    r_s   = smoothdata(tr.ft.r_speed(:), 'gaussian', smooth_window);
    is_walking = abs(r_s) > turn_thresh;
    is_walking = imclose(is_walking, ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking, min_walking_frames);
    dur_s(k) = xf(end)-xf(1); frac_walking(k) = mean(is_walking); rho_all{k} = rho_f; mean_rho(k) = mean(rho_f);
    for t = 1:n_thr
        keep = rho_f > rho_thr_list(t);
        [S.fly_path(k,t), S.bump_path(k,t), S.frac_kept(k,t)] = masked_paths(mu_s, r_s, dt, keep);
        [S.fly_path_w(k,t), S.bump_path_w(k,t), S.frac_kept_w(k,t)] = masked_paths(mu_s, r_s, dt, keep & is_walking);
    end
end
S.ratio   = S.bump_path   ./ S.fly_path;
S.ratio_w = S.bump_path_w ./ S.fly_path_w;

%% print per-trial table (whole trial)
fprintf('\n=== %s: whole-trial path-length ratio (bump / fly) vs rho threshold ===\n', example_tag);
fprintf('%-6s %-5s %-6s %-6s %-6s %-6s |', 'trial','light','phase','dur_s','walk%','rho');
fprintf(' %5.2f', rho_thr_list); fprintf('\n');
for k = 1:n_tr
    lt = 'CL'; if is_dark(k), lt = 'dark'; end
    fprintf('%-6d %-5s %-6s %-6.0f %-6.0f %-6.2f |', trial_num(idx(k)), lt, phase{k}, dur_s(k), 100*frac_walking(k), mean_rho(k));
    fprintf(' %5.2f', S.ratio(k,:)); fprintf('\n');
end
fprintf('\n--- fraction of samples kept at each threshold ---\n');
for k = 1:n_tr
    fprintf('%-6d %-5s %-6s %20s|', trial_num(idx(k)), '', phase{k}, '');
    fprintf(' %5.2f', S.frac_kept(k,:)); fprintf('\n');
end
fprintf('\n--- fly path (rad) / bump path (rad) at threshold 0 and 0.30 ---\n');
t0 = 1; t3 = find(abs(rho_thr_list-0.30)<1e-9);
for k = 1:n_tr
    fprintf('%-6d %-6s  thr0: fly %6.1f bump %6.1f ratio %.2f  |  thr0.30: fly %6.1f bump %6.1f ratio %.2f (kept %.0f%%)\n', ...
        trial_num(idx(k)), phase{k}, S.fly_path(k,t0), S.bump_path(k,t0), S.ratio(k,t0), S.fly_path(k,t3), S.bump_path(k,t3), S.ratio(k,t3), 100*S.frac_kept(k,t3));
end
fprintf('\n--- walking-only ratio at threshold 0 (for comparison with Figure 3''s bout metric) ---\n');
for k = 1:n_tr, fprintf('%-6d %-6s %.2f\n', trial_num(idx(k)), phase{k}, S.ratio_w(k,t0)); end

%% figure
fig = figure(34); clf
set(fig,'Color','w','Units','pixels','Position',[40 40 1700 950])
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
% colors: phase; line style: light
ph_col = containers.Map({'pre','PULSE','post'}, {[.3 .3 .3],[.85 .5 .1],[.85 .2 .2]});
ls = {'-','--'}; % CL, dark
lab = cell(n_tr,1);
for k = 1:n_tr, lt = 'CL'; if is_dark(k), lt = 'dark'; end; lab{k} = sprintf('tr %d %s %s', trial_num(idx(k)), lt, phase{k}); end

% (a) ratio vs threshold, whole trial
ax = nexttile(tl,1); hold(ax,'on')
for k = 1:n_tr, plot(ax, rho_thr_list, S.ratio(k,:), ls{is_dark(k)+1}, 'Color', ph_col(phase{k}), 'LineWidth', 1.5, 'Marker','.', 'DisplayName', lab{k}); end
yline(ax,1,':k','HandleVisibility','off'); xlabel(ax,'rho threshold (samples with rho <= thr discarded)'); ylabel(ax,'bump path / fly path (same kept samples)')
title(ax,'whole trial: path-length ratio vs rho threshold'); legend(ax,'Location','northwest','FontSize',7,'Box','off','NumColumns',2); ylim(ax,[0 3]); grid(ax,'on')

% (b) ratio vs threshold, walking only
ax = nexttile(tl,2); hold(ax,'on')
for k = 1:n_tr, plot(ax, rho_thr_list, S.ratio_w(k,:), ls{is_dark(k)+1}, 'Color', ph_col(phase{k}), 'LineWidth', 1.5, 'Marker','.'); end
yline(ax,1,':k'); xlabel(ax,'rho threshold'); ylabel(ax,'bump path / fly path, walking samples only')
title(ax,'walking samples only (Figure 3 bout definition): ratio vs rho threshold'); ylim(ax,[0 3]); grid(ax,'on')

% (c) fraction of samples kept
ax = nexttile(tl,3); hold(ax,'on')
for k = 1:n_tr, plot(ax, rho_thr_list, S.frac_kept(k,:), ls{is_dark(k)+1}, 'Color', ph_col(phase{k}), 'LineWidth', 1.5, 'Marker','.'); end
xlabel(ax,'rho threshold'); ylabel(ax,'fraction of trial samples kept'); ylim(ax,[0 1]); grid(ax,'on')
title(ax,'how much of each trial survives the threshold')

% (d) fly path and bump path vs threshold (whole trial), normalized to their threshold-0 value
ax = nexttile(tl,4); hold(ax,'on')
for k = 1:n_tr
    plot(ax, rho_thr_list, S.fly_path(k,:)/S.fly_path(k,1), ls{is_dark(k)+1}, 'Color', [0 .55 .8], 'LineWidth', 1.2)
    plot(ax, rho_thr_list, S.bump_path(k,:)/S.bump_path(k,1), ls{is_dark(k)+1}, 'Color', ph_col(phase{k}), 'LineWidth', 1.2)
end
xlabel(ax,'rho threshold'); ylabel(ax,'path length retained (fraction of threshold-0 value)'); ylim(ax,[0 1]); grid(ax,'on')
title(ax,'what the mask removes: fly path (blue) vs bump path (phase color), per trial')

% (e) rho distributions per trial
ax = nexttile(tl,5); hold(ax,'on')
for k = 1:n_tr
    [f,x] = ecdf(rho_all{k}); plot(ax, x, f, ls{is_dark(k)+1}, 'Color', ph_col(phase{k}), 'LineWidth', 1.5, 'DisplayName', lab{k});
end
xlabel(ax,'rho (interpolated onto behavior clock)'); ylabel(ax,'cumulative fraction of samples'); grid(ax,'on')
title(ax,'rho distribution per trial (ecdf)'); legend(ax,'Location','southeast','FontSize',7,'Box','off','NumColumns',2)

% (f) summary: ratio at a few thresholds, trial order on x
ax = nexttile(tl,6); hold(ax,'on')
show_thr = [0 .2 .3 .4 .5]; cm = copper(numel(show_thr)+1);
for j = 1:numel(show_thr)
    t = find(abs(rho_thr_list-show_thr(j))<1e-9);
    plot(ax, 1:n_tr, S.ratio(:,t), 'o-', 'Color', cm(j,:), 'MarkerFaceColor', cm(j,:), 'LineWidth', 1.3, 'DisplayName', sprintf('thr %.2f', show_thr(j)));
end
yline(ax,1,':k','HandleVisibility','off')
if ~isempty(p1), xline(ax, p1, '-', 'Color', ph_col('PULSE'), 'LineWidth', 2, 'HandleVisibility','off'); end
xticks(ax,1:n_tr); xticklabels(ax, lab); xtickangle(ax,40); ylabel(ax,'bump path / fly path (whole trial)'); ylim(ax,[0 3]); grid(ax,'on')
title(ax,'ratio per trial in chronological order (orange line = first pulse trial)'); legend(ax,'Location','northwest','FontSize',7,'Box','off')

title(tl, sprintf('%s (%s): whole-trial path-length mobility with rho-threshold sweep -- solid = closed loop, dashed = dark; gray = pre, orange = pulse trial, red = post', example_tag, geno), 'Interpreter','none')
exportgraphics(fig, fullfile(export_dir,'fig_pathlength_ratio_rho_sweep_fly12.png'), 'Resolution',130)
fprintf('\nexported to %s\n', fullfile(export_dir,'fig_pathlength_ratio_rho_sweep_fly12.png'))

%% helpers
function [fly_path, bump_path, frac_kept] = masked_paths(mu_s, r_s, dt, keep)
    % path lengths over kept samples only. bump path counts |diff| only
    % between consecutive kept samples (both ends kept), so the jump across
    % a discarded stretch is not counted as bump movement.
    frac_kept = mean(keep);
    fly_path  = sum(abs(r_s(keep))) * dt;
    both = keep(1:end-1) & keep(2:end);
    d = abs(diff(mu_s));
    bump_path = sum(d(both));
end

function [fid, tnum] = walk_meta_fly_trial(meta_path)
    parts = strsplit(meta_path,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    if ~(isfield(trial,'atp') && isfield(trial.atp,'d') && ~isempty(trial.atp.d)); tf = false; return; end
    stims = logical(trial.ft.stims(:)); onsets = find(diff([false;stims])==1);
    if isempty(onsets); tf = false; return; end
    xb = trial.ft.xb; xf = trial.ft.xf; t_onset = xf(onsets); atp_sig = sum(trial.atp.d,1);
    t_common = linspace(-peri_win,peri_win,101);
    atp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets), atp_aligned(s,:) = interp1(xb - t_onset(s),atp_sig,t_common,'linear',nan); end
    atp_peri = mean(atp_aligned,1,'omitnan');
    base_idx = t_common>=base_win(1) & t_common<base_win(2); peak_idx = t_common>peak_win(1) & t_common<=peak_win(2);
    tf = (median(atp_peri(peak_idx),'omitnan')-median(atp_peri(base_idx),'omitnan')) > peak_factor*std(atp_peri(base_idx),'omitnan');
end
