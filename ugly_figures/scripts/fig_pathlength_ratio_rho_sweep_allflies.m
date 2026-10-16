%% fig_pathlength_ratio_rho_sweep_allflies
% Population version of fig_pathlength_ratio_rho_sweep.m (whole-trial
% bump/fly PATH-LENGTH RATIO with a rho-threshold mask applied to both
% signals -- see that script's header for the exact definition). Here:
%
%   1) per fly, trials are pooled into six categories, in chronological
%      order of the experiment:  CL pre | dark pre | CL pulse 1 | CL post |
%      CL pulse 2 | dark post.  "pre" = before the first detected ATP-pulse
%      trial, "post" = after it, excluding pulse trials (so "CL post" also
%      includes the closed-loop trial(s) between the two pulse trials, as
%      Figure 3 does). Pooling = summing bump path and fly path over the
%      category's trials before taking the ratio, which is identical to
%      concatenating the trials (the bump path never counts a jump across
%      a trial boundary because the mask-aware path breaks there anyway).
%   2) each fly is one dot per category; the two genotypes are offset side
%      by side; thick markers = genotype mean +/- SEM over flies.
%   3) one subpanel per rho threshold.
%
% Smoothing / clocks / mask conventions are identical to the single-fly
% script (1 s Gaussian on upsampled mu and on r_speed, mask = rho_f > thr
% applied to both numerator and denominator, no bout segmentation).
%
% !!! RHO / MOBILITY CAVEAT (Figure_3_claude.m): rho and bump movement are
% coupled, and the single-fly sweep showed the ratio falls monotonically
% with threshold in EVERY trial because low-rho samples carry
% disproportionate bump path. There is no threshold-independent answer;
% what this figure tests is whether the pre -> post ORDERING is stable
% across thresholds.

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
rho_thr_list  = [0 0.2 0.3 0.4 0.5 0.6]; % one subpanel each
smooth_window = 60;
cat_names = {'CL pre','dark pre','CL pulse 1','CL post','CL pulse 2','dark post'};
min_fly_path = 10; % rad of total (masked) fly rotation a fly x category cell needs before its ratio is trusted (~1.6 turns); below this the ratio blows up
n_cat = numel(cat_names);
geno_keep   = {'empty>p2x2','lpsp>p2x2'};
geno_colors = [0.00,0.45,0.74; 0.85,0.33,0.10];

%% labels
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);
fly_geno = cell(n_flies,1);

% category per trial (NaN = fly without pulse trial, or 3rd+ pulse trial)
trial_cat = nan(n,1);
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = geno_keep{1 + is_lpsp(idx(1))};
    pp = find(has_pulse(idx)); if isempty(pp), continue, end
    pre = idx(1:pp(1)-1); post = idx(pp(1)+1:end); post = post(~has_pulse(post));
    trial_cat(pre(~is_dark(pre)))   = 1; trial_cat(pre(is_dark(pre)))   = 2;
    trial_cat(idx(pp(1)))           = 3;
    trial_cat(post(~is_dark(post))) = 4; trial_cat(post(is_dark(post))) = 6;
    if numel(pp) >= 2, trial_cat(idx(pp(2))) = 5; end
    if numel(pp) >= 3, fprintf('note: %s has %d pulse trials; only the first two are shown as pulse categories\n', fly_list{f}, numel(pp)); end
end
dark_pulse = find(has_pulse & is_dark(:) & ~isnan(trial_cat));
if ~isempty(dark_pulse), fprintf('WARNING: %d pulse trial(s) are dark trials but are placed in the "CL pulse" categories\n', numel(dark_pulse)); end

%% per-trial sweep: fly path and bump path for each threshold
n_thr = numel(rho_thr_list);
fly_path = nan(n,n_thr); bump_path = nan(n,n_thr); n_kept = nan(n,n_thr); n_samp = nan(n,1);
for i = find(~isnan(trial_cat))'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf));
    mu_s  = smoothdata(interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap'), 'gaussian', smooth_window);
    rho_f = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    r_s   = smoothdata(tr.ft.r_speed(:), 'gaussian', smooth_window);
    n_samp(i) = numel(xf);
    for t = 1:n_thr
        keep = rho_f > rho_thr_list(t);
        [fly_path(i,t), bump_path(i,t), fk] = masked_paths(mu_s, r_s, dt, keep);
        n_kept(i,t) = fk*numel(xf);
    end
end

%% pool per fly x category (sum paths, then ratio)
R = nan(n_flies,n_cat,n_thr); FK = nan(n_flies,n_cat,n_thr); FP = nan(n_flies,n_cat,n_thr); NT = zeros(n_flies,n_cat);
for f = 1:n_flies
    for c = 1:n_cat
        ii = find(fly_ix==f & trial_cat==c); NT(f,c) = numel(ii);
        if isempty(ii), continue, end
        fp = sum(fly_path(ii,:),1); bp = sum(bump_path(ii,:),1);
        r = bp ./ fp; r(fp < min_fly_path) = nan; % too little fly rotation -> ratio unstable
        R(f,c,:)  = r; FP(f,c,:) = fp;
        FK(f,c,:) = sum(n_kept(ii,:),1) / sum(n_samp(ii));
    end
end
used = ~all(isnan(R(:,:,1)),2);
fprintf('%d flies with a pulse trial (%d empty, %d lpsp)\n', sum(used), sum(used & strcmp(fly_geno,geno_keep{1})), sum(used & strcmp(fly_geno,geno_keep{2})));
fprintf('trials per category (mean over flies): %s\n', strjoin(arrayfun(@(c) sprintf('%s=%.1f', cat_names{c}, mean(NT(used,c))), 1:n_cat, 'UniformOutput', false), ', '));
% cells dropped for too little fly rotation, and large-ratio cells, at the lowest threshold
fprintf('\nfly x category cells with total fly rotation < %g rad (dropped) or ratio > 5 at rho thr %.2f:\n', min_fly_path, rho_thr_list(1));
for f = find(used)'
    for c = 1:n_cat
        if NT(f,c) == 0, continue, end
        if FP(f,c,1) < min_fly_path || R(f,c,1) > 5
            fprintf('  %-22s %-11s %-10s fly path %6.1f rad, bump path %7.1f rad, ratio %6.2f, %d trial(s)\n', strrep(fly_list{f},'Z:/pablo/lpsp_p2x2_walking/',''), fly_geno{f}, cat_names{c}, FP(f,c,1), R(f,c,1)*FP(f,c,1), R(f,c,1), NT(f,c));
        end
    end
end
n_drop = squeeze(sum(sum(~isnan(FP) & FP < min_fly_path,1),2));
fprintf('cells dropped per threshold (%s): %s\n', mat2str(rho_thr_list), mat2str(n_drop(:)'));

%% print table: genotype mean (SEM) per category per threshold
for t = 1:n_thr
    fprintf('\n== rho thr %.2f: bump/fly path ratio, median [IQR] (n flies) ==\n', rho_thr_list(t));
    fprintf('%-12s', 'geno'); fprintf(' %-18s', cat_names{:}); fprintf('\n');
    for gi = 1:2
        fg = used & strcmp(fly_geno,geno_keep{gi});
        fprintf('%-12s', geno_keep{gi});
        for c = 1:n_cat
            v = R(fg,c,t); v = v(~isnan(v));
            fprintf(' %5.2f [%4.2f-%4.2f] (%2d)', median(v), prctile(v,25), prctile(v,75), numel(v));
        end
        fprintf('\n');
    end
    % paired dark post - dark pre per genotype
    for gi = 1:2
        fg = used & strcmp(fly_geno,geno_keep{gi});
        d = R(fg,6,t) - R(fg,2,t); d = d(~isnan(d));
        p = nan; if numel(d) >= 5, p = signrank(d); end
        fprintf('   %-10s dark post - dark pre: median %+.2f (mean %+.2f), n=%d paired, signrank p=%.3f | median kept frac pre %.2f, post %.2f\n', ...
            geno_keep{gi}, median(d), mean(d), numel(d), p, median(FK(fg,2,t),'omitnan'), median(FK(fg,6,t),'omitnan'));
    end
end

%% figure
fig = figure(35); clf
set(fig,'Color','w','Units','pixels','Position',[40 40 1750 900])
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
xoff = [-0.17, 0.17];
for t = 1:n_thr
    ax = nexttile(tl,t); hold(ax,'on')
    for gi = 1:2
        fg = find(used & strcmp(fly_geno,geno_keep{gi}));
        X = (1:n_cat) + xoff(gi);
        Y = reshape(R(fg,:,t), numel(fg), n_cat);
        % faint per-fly lines across the categories each fly has
        for k = 1:numel(fg)
            ok = ~isnan(Y(k,:));
            plot(ax, X(ok), Y(k,ok), '-', 'Color', [geno_colors(gi,:) .18], 'LineWidth', .8, 'HandleVisibility', 'off')
        end
        % per-fly dots with a little jitter
        jit = (rand(numel(fg),1)-.5)*0.10;
        scatter(ax, repmat(X,numel(fg),1) + jit, Y, 18, geno_colors(gi,:), 'filled', 'MarkerFaceAlpha', .55, 'HandleVisibility', 'off')
        % genotype median with IQR (means are dominated by a few near-stationary cells; see the printed outlier list)
        m = median(Y,1,'omitnan'); lo = m - prctile(Y,25,1); hi = prctile(Y,75,1) - m;
        errorbar(ax, X, m, lo, hi, 'o', 'Color', geno_colors(gi,:), 'MarkerFaceColor', 'w', 'MarkerSize', 7, 'LineWidth', 1.8, 'CapSize', 6, ...
            'DisplayName', sprintf('%s (n=%d flies)', geno_keep{gi}, numel(fg)))
    end
    yline(ax, 1, ':k', 'HandleVisibility', 'off')
    xticks(ax, 1:n_cat); xticklabels(ax, cat_names); xtickangle(ax, 30); xlim(ax, [.4 n_cat+.6])
    set(ax, 'YScale', 'log'); ylim(ax, [0.3 30]); yticks(ax, [0.3 0.5 1 2 5 10 30]); grid(ax, 'on')
    ylabel(ax, 'bump path / fly path (pooled trials)')
    kept_str = sprintf('median kept: %.0f%%', 100*median(FK(used,:,t),'all','omitnan'));
    title(ax, sprintf('rho threshold %.2f  (%s)', rho_thr_list(t), kept_str))
    if t == 1, legend(ax, 'Location', 'northwest', 'FontSize', 8, 'Box', 'off', 'Interpreter', 'none'), end
end
title(tl, {'whole-trial bump/fly path-length ratio per fly, trials pooled by category (chronological order left to right)', ...
    sprintf('dot = fly, faint line = same fly, open circle = genotype median with IQR (log y); mask (rho > thr) applied to bump AND fly path; no bout segmentation; cells with < %g rad total fly rotation dropped', min_fly_path)}, 'FontSize', 11)
exportgraphics(fig, fullfile(export_dir,'fig_pathlength_ratio_rho_sweep_allflies.png'), 'Resolution', 130)
fprintf('\nexported to %s\n', fullfile(export_dir,'fig_pathlength_ratio_rho_sweep_allflies.png'))

%% helpers
function [fly_path, bump_path, frac_kept] = masked_paths(mu_s, r_s, dt, keep)
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
