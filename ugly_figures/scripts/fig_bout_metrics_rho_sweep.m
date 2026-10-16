%% fig_bout_metrics_rho_sweep
% Bout-based bump mobility (Figure 3's walking-bout framework) under three
% definitions, with bump samples below a rho threshold NaN'd out, swept
% over the threshold. PRE-stim trials only, both genotypes.
%
% Per trial:
%   * bump (unwrap im.mu) and rho interpolated onto the 60 Hz behavior
%     clock; heading = cumtrapz of r_speed.
%   * MASK FIRST: bump samples with rho <= thr -> NaN. Then the 1 s
%     Gaussian (60 samples, Figure 3's window) on bump and on r_speed;
%     smoothdata omits NaNs so each kept stretch is smoothed on its own;
%     masked samples are re-NaN'd after smoothing. Heading is integrated
%     from the smoothed r_speed (unmasked -- the mask is applied when the
%     per-bout numbers are taken, below).
%   * walking bouts from |smoothed r_speed| > 0.25 rad/s, gaps <= 0.5 s
%     closed, bouts < 0.5 s dropped (Figure 3).
%   * per bout, on the KEPT samples only:
%       path:   bump = sum|diff(bump)| over consecutive kept pairs;
%               heading = sum|r_speed| dt over kept samples
%       |net|:  |bump(last kept) - bump(first kept)| vs same for heading
%       signed: the same differences with sign
%     A bout is dropped if < 50% of its samples are kept or < 2 are kept.
%   * per fly x category: sqrt(duration)-weighted through-origin slope
%     over pooled bouts (>= 15 bouts), duration = kept samples * dt.
% thr = 0 keeps everything and reproduces the earlier ungated numbers.

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
rho_thr_list = [0 0.2 0.3 0.4 0.5];
smooth_window = 60; turn_thresh = .25; max_gap = 30; min_walk = 30; min_bouts = 15; min_kept_frac = .5;
cat_names = {'CL pre','dark pre'}; n_cat = 2;
geno_keep = {'empty>p2x2','lpsp>p2x2'}; gc = [0 .45 .74; .85 .33 .1];
def_names = {'path-length slope','|net displacement| slope','signed net displacement slope'};

%% labels
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);
trial_cat = nan(n,1); fly_geno = cell(n_flies,1);
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = geno_keep{1+is_lpsp(idx(1))};
    p = find(has_pulse(idx),1); if isempty(p), continue, end
    pre = idx(1:p-1); trial_cat(pre(~is_dark(pre))) = 1; trial_cat(pre(is_dark(pre))) = 2;
end

%% per trial x threshold: per-bout quantities
n_thr = numel(rho_thr_list);
B = cell(n, n_thr); % each: struct with fields path_b, path_h, net_b, net_h, dur
for i = find(~isnan(trial_cat))'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf));
    mu_f  = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    rho_f = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    r_s   = smoothdata(tr.ft.r_speed(:), 'gaussian', smooth_window); r_s(isnan(r_s)) = 0;
    head  = cumtrapz(xf, r_s);
    w = abs(r_s) > turn_thresh; w = imclose(w, ones(max_gap+1,1)); w = bwareaopen(w, min_walk);
    d = diff([false; w; false]); b0 = find(d==1); b1 = find(d==-1)-1; nb = numel(b0);
    for t = 1:n_thr
        keep = rho_f > rho_thr_list(t);
        m = mu_f; m(~keep) = nan;
        if smooth_window > 1, ms = smoothdata(m, 'gaussian', smooth_window); else, ms = m; end
        ms(~keep) = nan;
        S = struct('path_b', nan(nb,1), 'path_h', nan(nb,1), 'net_b', nan(nb,1), 'net_h', nan(nb,1), 'dur', nan(nb,1));
        for b = 1:nb
            seg = (b0(b):b1(b))'; k = keep(seg);
            if sum(k) < 2 || mean(k) < min_kept_frac, continue, end
            segk = seg(k);
            both = k(1:end-1) & k(2:end); dm = abs(diff(ms(seg)));
            S.path_b(b) = sum(dm(both)); S.path_h(b) = sum(abs(r_s(segk)))*dt;
            S.net_b(b) = ms(segk(end)) - ms(segk(1)); S.net_h(b) = head(segk(end)) - head(segk(1));
            S.dur(b) = numel(segk)*dt;
        end
        ok = ~isnan(S.dur); S = structfun(@(v) v(ok), S, 'UniformOutput', false);
        B{i,t} = S;
    end
end

%% per fly x category x threshold x definition
R = nan(n_flies, n_cat, n_thr, 3); NB = zeros(n_flies, n_cat, n_thr);
sl = @(x,y,d) (sqrt(d).*x) \ (sqrt(d).*y);
for f = 1:n_flies
    for c = 1:n_cat
        ii = find(fly_ix==f & trial_cat==c); if isempty(ii), continue, end
        for t = 1:n_thr
            S = [B{ii,t}]; if isempty(S), continue, end
            pb = vertcat(S.path_b); ph = vertcat(S.path_h); nb_ = vertcat(S.net_b); nh = vertcat(S.net_h); du = vertcat(S.dur);
            NB(f,c,t) = numel(du); if numel(du) < min_bouts, continue, end
            R(f,c,t,1) = sl(ph, pb, du); R(f,c,t,2) = sl(abs(nh), abs(nb_), du); R(f,c,t,3) = sl(nh, nb_, du);
        end
    end
end
used = any(~isnan(R(:,:,1,1)),2);

%% report
for t = 1:n_thr
    fprintf('\n== rho thr %.2f ==  (median [IQR] (n flies); bouts kept: median per fly x cat)\n%-10s %-10s', rho_thr_list(t), 'category','geno');
    fprintf(' %-26s', def_names{:}); fprintf(' %s\n', 'nbouts');
    for c = 1:n_cat, for gi = 1:2
        sel = used & strcmp(fly_geno, geno_keep{gi});
        fprintf('%-10s %-10s', cat_names{c}, geno_keep{gi});
        for k = 1:3, v = R(sel,c,t,k); v = v(~isnan(v)); fprintf(' %5.2f [%5.2f-%5.2f] (%2d)   ', median(v), prctile(v,25), prctile(v,75), numel(v)); end
        fprintf(' %d\n', round(median(NB(sel,c,t))));
    end, end
end

%% figure: rows = threshold, cols = definition
fig = figure(83); clf; set(fig,'Color','w','Units','pixels','Position',[30 30 1300 240*n_thr])
tl = tiledlayout(n_thr, 3, 'TileSpacing','compact','Padding','compact');
xo = [-.12 .12];
for t = 1:n_thr
    for k = 1:3
        ax = nexttile(tl, (t-1)*3+k); hold(ax,'on')
        for gi = 1:2
            sel = find(used & strcmp(fly_geno, geno_keep{gi}));
            Y = reshape(R(sel,:,t,k), numel(sel), n_cat); X = (1:n_cat) + xo(gi);
            for f = 1:size(Y,1), plot(ax, X, Y(f,:), '-o', 'Color', [gc(gi,:) .35], 'MarkerFaceColor', gc(gi,:), 'MarkerSize', 3.5, 'HandleVisibility','off'); end
            m = median(Y,1,'omitnan'); lo = m - prctile(Y,25,1); hi = prctile(Y,75,1) - m;
            errorbar(ax, X, m, lo, hi, 'o', 'Color', gc(gi,:), 'MarkerFaceColor','w', 'LineWidth', 1.6, 'CapSize', 7, 'DisplayName', sprintf('%s (n=%d)', geno_keep{gi}, sum(any(~isnan(Y),2))))
        end
        yline(ax, 1, ':k', 'HandleVisibility','off'); yline(ax, .73, '--', 'Color', [.5 .5 .5], 'HandleVisibility','off'); if k == 3, yline(ax, 0, '-k', 'HandleVisibility','off'); end
        xticks(ax, 1:n_cat); xlim(ax, [.5 n_cat+.5]); ylim(ax, [-.5 3.2]); grid(ax,'on'); set(ax,'FontSize',8)
        if t == n_thr, xticklabels(ax, cat_names); else, xticklabels(ax, []); end
        if k == 1, ylabel(ax, sprintf('rho thr %.1f\nbump / heading per bout', rho_thr_list(t))); end
        if t == 1, title(ax, def_names{k}); end
        if t == 1 && k == 1, legend(ax, 'Location','northwest','Box','off','Interpreter','none','FontSize',7); end
    end
end
title(tl, {'PRE-stim bout-based mobility, bump samples with rho <= thr NaN''d before smoothing; per-bout numbers on kept samples only; bouts with < 50% kept dropped', ...
    'dotted = 1, dashed = 0.73 cue gain; dot = fly, circle = median with IQR'}, 'FontSize', 10)
exportgraphics(fig, fullfile(export_dir,'fig_bout_metrics_rho_sweep.png'), 'Resolution', 130)
fprintf('\nexported to %s\n', fullfile(export_dir,'fig_bout_metrics_rho_sweep.png'))

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
