%% fig_pathlength_ratio_smooth_sweep_allflies
% Whole-trial bump/fly PATH-LENGTH RATIO per fly, pooled by category, as in
% fig_pathlength_ratio_rho_sweep_allflies.m, but now sweeping the Gaussian
% smoothing window (columns) as well as the rho threshold (rows), with two
% changes to how smoothing is done:
%
%   1) EQUIVALENT SMOOTHING OF BOTH TRACES. The bump trace (unwrapped im.mu,
%      upsampled onto the behavior clock) and the fly's heading (cumtrapz of
%      ft.r_speed, i.e. the same "heading" Figure 3 uses) are both smoothed
%      with the same Gaussian window, and both path lengths are then
%      sum|diff| of the smoothed trace over consecutive kept samples. The
%      earlier scripts smoothed r_speed and summed |r_speed|*dt, which is
%      the same operation on the heading, so this is a cleanup rather than
%      a change -- but it makes the symmetry explicit. NOTE the dataset
%      itself is not symmetric: im.mu inherits a 10-frame (~0.95 s)
%      Gaussian on the movie, r_speed a 30-sample (~0.5 s) Gaussian in
%      process_ft (lpsp_p2x2_script_minimal.m). That cannot be undone here.
%   2) MASK BEFORE SMOOTHING. Samples with rho <= threshold are set to NaN
%      in BOTH traces before smoothing (smoothdata omits NaNs, so each kept
%      segment is smoothed using only its own samples -- a discarded
%      stretch never leaks into its neighbours), then re-masked after
%      smoothing so the path only counts pairs of consecutive kept samples.
%      The earlier scripts smoothed first and masked afterwards.
%
% Window = 1 sample means no additional smoothing beyond what the dataset
% carries. Everything else (categories, pooling by summing paths, the
% minimum total fly rotation per cell, median/IQR on a log axis) is as in
% fig_pathlength_ratio_rho_sweep_allflies.m.
%
% !!! RHO / MOBILITY CAVEAT (Figure_3_claude.m) still applies: the mask
% removes bump path faster than fly path at every threshold.

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
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
rho_thr_list = [0.2 0.3 0.4];          % rows
win_list     = [1 30 60 120 300];       % columns, samples on the ~60 Hz behavior clock
win_s_label  = @(w,dt) sprintf('%.2g s', w*dt);
min_fly_path = 10;                      % rad total masked fly rotation per fly x category cell
cat_names = {'CL pre','dark pre','CL pulse 1','CL post','CL pulse 2','dark post'}; n_cat = numel(cat_names);
geno_keep   = {'empty>p2x2','lpsp>p2x2'};
geno_colors = [0.00,0.45,0.74; 0.85,0.33,0.10];

%% labels (as in the rho-sweep script)
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);
fly_geno = cell(n_flies,1); trial_cat = nan(n,1);
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = geno_keep{1 + is_lpsp(idx(1))};
    pp = find(has_pulse(idx)); if isempty(pp), continue, end
    pre = idx(1:pp(1)-1); post = idx(pp(1)+1:end); post = post(~has_pulse(post));
    trial_cat(pre(~is_dark(pre))) = 1; trial_cat(pre(is_dark(pre))) = 2; trial_cat(idx(pp(1))) = 3;
    trial_cat(post(~is_dark(post))) = 4; trial_cat(post(is_dark(post))) = 6;
    if numel(pp) >= 2, trial_cat(idx(pp(2))) = 5; end
end

%% per-trial paths for every (threshold, window)
n_thr = numel(rho_thr_list); n_win = numel(win_list);
fly_path = nan(n,n_thr,n_win); bump_path = nan(n,n_thr,n_win); dt_all = nan(n,1);
for i = find(~isnan(trial_cat))'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf)); dt_all(i) = dt;
    mu_f    = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    rho_f   = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    r = tr.ft.r_speed(:); r(isnan(r)) = 0; heading = cumtrapz(xf, r);
    for t = 1:n_thr
        keep = rho_f > rho_thr_list(t);
        m = mu_f; m(~keep) = nan; h = heading; h(~keep) = nan;
        for w = 1:n_win
            if win_list(w) > 1
                ms = smoothdata(m, 'gaussian', win_list(w)); hs = smoothdata(h, 'gaussian', win_list(w)); % NaNs omitted -> each kept segment smoothed on its own
            else
                ms = m; hs = h;
            end
            ms(~keep) = nan; hs(~keep) = nan;
            both = keep(1:end-1) & keep(2:end);
            dm = abs(diff(ms)); dh = abs(diff(hs));
            bump_path(i,t,w) = sum(dm(both)); fly_path(i,t,w) = sum(dh(both));
        end
    end
end
dt_med = median(dt_all,'omitnan');

%% pool per fly x category
R = nan(n_flies,n_cat,n_thr,n_win); FP = R; NT = zeros(n_flies,n_cat);
for f = 1:n_flies
    for c = 1:n_cat
        ii = find(fly_ix==f & trial_cat==c); NT(f,c) = numel(ii);
        if isempty(ii), continue, end
        fp = squeeze(sum(fly_path(ii,:,:),1)); bp = squeeze(sum(bump_path(ii,:,:),1));
        r = bp ./ fp; r(fp < min_fly_path) = nan;
        R(f,c,:,:) = r; FP(f,c,:,:) = fp;
    end
end
used = ~all(isnan(R(:,:,1,1)),2);
fprintf('%d flies (%d empty, %d lpsp); windows = %s samples = %s\n', sum(used), sum(used & strcmp(fly_geno,geno_keep{1})), sum(used & strcmp(fly_geno,geno_keep{2})), ...
    mat2str(win_list), strjoin(arrayfun(@(w) win_s_label(w,dt_med), win_list, 'UniformOutput', false), ', '));

%% report: medians per category + paired dark post - pre, per (thr, win)
for t = 1:n_thr
    for w = 1:n_win
        fprintf('\n== rho thr %.2f, window %d samples (%s): median [IQR] (n) ==\n', rho_thr_list(t), win_list(w), win_s_label(win_list(w),dt_med));
        fprintf('%-12s', 'geno'); fprintf(' %-20s', cat_names{:}); fprintf('\n');
        for gi = 1:2
            fg = used & strcmp(fly_geno,geno_keep{gi}); fprintf('%-12s', geno_keep{gi});
            for c = 1:n_cat
                v = R(fg,c,t,w); v = v(~isnan(v));
                fprintf(' %5.2f [%4.2f-%4.2f] (%2d)', median(v), prctile(v,25), prctile(v,75), numel(v));
            end
            fprintf('\n');
        end
        for gi = 1:2
            fg = used & strcmp(fly_geno,geno_keep{gi});
            d = R(fg,6,t,w) - R(fg,2,t,w); d = d(~isnan(d)); p = nan; if numel(d) >= 5, p = signrank(d); end
            fprintf('   %-10s dark post - dark pre: median %+.2f, n=%d, signrank p=%.3f\n', geno_keep{gi}, median(d), numel(d), p);
        end
    end
end

%% figure: rows = rho threshold, cols = smoothing window
fig = figure(36); clf
set(fig,'Color','w','Units','pixels','Position',[20 20 1900 1000])
tl = tiledlayout(n_thr, n_win, 'TileSpacing','compact','Padding','compact');
xoff = [-0.17, 0.17];
for t = 1:n_thr
    for w = 1:n_win
        ax = nexttile(tl, (t-1)*n_win + w); hold(ax,'on')
        for gi = 1:2
            fg = find(used & strcmp(fly_geno,geno_keep{gi}));
            X = (1:n_cat) + xoff(gi); Y = reshape(R(fg,:,t,w), numel(fg), n_cat);
            for k = 1:numel(fg), ok = ~isnan(Y(k,:)); plot(ax, X(ok), Y(k,ok), '-', 'Color', [geno_colors(gi,:) .18], 'LineWidth', .7, 'HandleVisibility','off'); end
            jit = (rand(numel(fg),1)-.5)*0.10;
            scatter(ax, repmat(X,numel(fg),1)+jit, Y, 12, geno_colors(gi,:), 'filled', 'MarkerFaceAlpha', .5, 'HandleVisibility','off')
            m = median(Y,1,'omitnan'); lo = m - prctile(Y,25,1); hi = prctile(Y,75,1) - m;
            errorbar(ax, X, m, lo, hi, 'o', 'Color', geno_colors(gi,:), 'MarkerFaceColor','w', 'MarkerSize', 6, 'LineWidth', 1.6, 'CapSize', 5, ...
                'DisplayName', sprintf('%s (n=%d)', geno_keep{gi}, numel(fg)))
        end
        yline(ax, 1, ':k', 'HandleVisibility','off')
        set(ax, 'YScale','log'); ylim(ax, [0.3 10]); yticks(ax, [0.3 0.5 1 2 5 10]); grid(ax,'on')
        xticks(ax, 1:n_cat); xlim(ax, [.4 n_cat+.6])
        if t == n_thr, xticklabels(ax, cat_names); xtickangle(ax, 35); else, xticklabels(ax, []); end
        if w == 1, ylabel(ax, sprintf('rho thr %.1f\nbump / fly path', rho_thr_list(t))); end
        if t == 1, title(ax, sprintf('Gaussian window %d samples (%s)', win_list(w), win_s_label(win_list(w),dt_med)), 'FontSize', 9); end
        if t == 1 && w == 1, legend(ax, 'Location','northwest', 'FontSize', 7, 'Box','off', 'Interpreter','none'); end
        set(ax, 'FontSize', 8)
    end
end
title(tl, {'whole-trial bump/fly path-length ratio per fly, pooled by category -- rows: rho threshold, columns: Gaussian smoothing window (same window on bump and heading, applied within kept segments only)', ...
    sprintf('dot = fly, open circle = genotype median with IQR (log y); mask applied before smoothing and to both traces; cells with < %g rad total fly rotation dropped; window 1 = dataset smoothing only', min_fly_path)}, 'FontSize', 10)
exportgraphics(fig, fullfile(export_dir,'fig_pathlength_ratio_smooth_sweep_allflies.png'), 'Resolution', 130)
fprintf('\nexported to %s\n', fullfile(export_dir,'fig_pathlength_ratio_smooth_sweep_allflies.png'))

%% helpers
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
