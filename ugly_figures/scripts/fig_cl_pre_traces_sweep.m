%% fig_cl_pre_traces_sweep
% Sanity check of the path-length mobility metric on the cleanest case:
% empty>p2x2 control flies, closed-loop trials BEFORE any ATP. For each fly
% a 3 x 5 grid (rows = rho threshold 0.2/0.3/0.4, cols = Gaussian window
% 0/0.5/1/2/5 s) shows a 60 s snippet of the fly's heading (black), the cue
% (cyan) and the bump (orange), all unwrapped, after the mask-then-smooth
% treatment of fig_pathlength_ratio_smooth_sweep_allflies.m. Raw unsmoothed
% bump is drawn faintly underneath; discarded (rho <= thr) stretches are
% shaded. Each panel's title gives that fly's pooled CL-pre ratios:
%   bump/heading -- what the mobility metric reports. Expected ~cue gain
%                   (0.73 on this rig, measured in lpsp_p2x2_walking_turns_
%                   claude.m) for a bump that tracks the cue perfectly.
%   bump/cue     -- expected ~1 for perfect tracking. This is the direct
%                   "is the metric near 1 for controls in closed loop" test.
% Cue path uses unwrap(-ft.cue) with the same mask and smoothing as the
% other two traces. Heading = cumtrapz(r_speed).
%
% A summary figure then shows, per (thr, window), the per-fly CL-pre
% bump/cue and bump/heading ratios across all empty flies.

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports','cl_pre_traces');
if ~isfolder(export_dir), mkdir(export_dir), end
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% parameters
rho_thr_list = [0.2 0.3 0.4];
win_list     = [1 30 60 120 300];
snip_len     = 60; % s
geno_show    = 'empty>p2x2';

%% labels
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);
cl_pre = cell(n_flies,1); fly_geno = cell(n_flies,1);
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = 'empty>p2x2'; if is_lpsp(idx(1)), fly_geno{f} = 'lpsp>p2x2'; end
    p = find(has_pulse(idx),1); if isempty(p), continue, end
    pre = idx(1:p-1); cl_pre{f} = pre(~is_dark(pre));
end
flies = find(strcmp(fly_geno,geno_show) & ~cellfun(@isempty,cl_pre));
fprintf('%d %s flies with closed-loop pre-ATP trials\n', numel(flies), geno_show);

%% per trial: traces on behavior clock
n_thr = numel(rho_thr_list); n_win = numel(win_list);
T = struct(); % per trial cache
for i = unique(vertcat(cl_pre{flies}))'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:);
    T(i).xf = xf; T(i).dt = median(diff(xf));
    T(i).mu   = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    T(i).rho  = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    r = tr.ft.r_speed(:); r(isnan(r)) = 0; T(i).head = cumtrapz(xf, r);
    c = -tr.ft.cue(:); c(isnan(c)) = 0; T(i).cue = unwrap(c);
    r_s = smoothdata(tr.ft.r_speed(:),'gaussian',60); w = abs(r_s) > .25; w = imclose(w,ones(31,1)); w = bwareaopen(w,30); T(i).walk = w;
end

%% pooled CL-pre ratios per fly x (thr, win), plus per-fly figures
ratio_h = nan(numel(flies),n_thr,n_win); ratio_c = ratio_h; kept = ratio_h;
for k = 1:numel(flies)
    f = flies(k); ii = cl_pre{f};
    % snippet: from the fly's first CL-pre trial, the 60 s window with most walking coverage
    i0 = ii(1); dt = T(i0).dt; wlen = round(snip_len/dt);
    cover = movmean(double(T(i0).walk), wlen, 'Endpoints','discard'); [~,k0] = max(cover);
    seg = k0:min(k0+wlen-1, numel(T(i0).xf)); t_seg = T(i0).xf(seg) - T(i0).xf(seg(1));

    fig = figure(40); clf; set(fig,'Color','w','Units','pixels','Position',[20 20 1900 950])
    tl = tiledlayout(n_thr, n_win, 'TileSpacing','compact','Padding','compact');
    for t = 1:n_thr
        for w = 1:n_win
            % pooled ratios over all CL-pre trials
            bp = 0; hp = 0; cp = 0; nk = 0; ns = 0;
            for i = ii'
                keep = T(i).rho > rho_thr_list(t);
                [ms, hs, cs] = mask_smooth(T(i).mu, T(i).head, T(i).cue, keep, win_list(w));
                both = keep(1:end-1) & keep(2:end);
                dm = abs(diff(ms)); dh = abs(diff(hs)); dc = abs(diff(cs));
                bp = bp + sum(dm(both)); hp = hp + sum(dh(both)); cp = cp + sum(dc(both));
                nk = nk + sum(keep); ns = ns + numel(keep);
            end
            ratio_h(k,t,w) = bp/hp; ratio_c(k,t,w) = bp/cp; kept(k,t,w) = nk/ns;

            % snippet traces from the first trial
            keep = T(i0).rho > rho_thr_list(t);
            [ms, hs, cs] = mask_smooth(T(i0).mu, T(i0).head, T(i0).cue, keep, win_list(w));
            ax = nexttile(tl,(t-1)*n_win+w); hold(ax,'on')
            % re-zero every trace at the first kept sample in the snippet
            k1 = seg(find(keep(seg),1)); if isempty(k1), k1 = seg(1); end
            raw = T(i0).mu(seg) - T(i0).mu(k1);
            h0 = hs(k1); c0 = cs(k1); m0 = ms(k1);
            if isnan(h0), h0 = T(i0).head(k1); end; if isnan(c0), c0 = T(i0).cue(k1); end; if isnan(m0), m0 = T(i0).mu(k1); end
            % shade discarded stretches
            yl = [-1 1]*max([abs(raw); abs(hs(seg)-h0); abs(cs(seg)-c0); 2*pi],[],'omitnan')*1.1;
            d = diff([false; ~keep(seg); false]); on = find(d==1); off = find(d==-1)-1;
            for g = 1:numel(on), patch(ax, t_seg([on(g) off(g) off(g) on(g)]), yl([1 1 2 2]), [.92 .92 .92], 'EdgeColor','none'); end
            plot(ax, t_seg, raw, '-', 'Color', [.75 .75 .75], 'LineWidth', .6)
            plot(ax, t_seg, cs(seg)-c0, '-', 'Color', [0 .7 .85], 'LineWidth', 1.1)
            plot(ax, t_seg, hs(seg)-h0, 'k-', 'LineWidth', 1.1)
            plot(ax, t_seg, ms(seg)-m0, '-', 'Color', [.85 .33 .1], 'LineWidth', 1.1)
            ylim(ax, yl); xlim(ax, [0 snip_len]); set(ax,'FontSize',7)
            title(ax, sprintf('thr %.1f, win %s | kept %.0f%% | bump/head %.2f, bump/cue %.2f', rho_thr_list(t), win_label(win_list(w),dt), 100*kept(k,t,w), ratio_h(k,t,w), ratio_c(k,t,w)), 'FontSize', 7.5)
            if w == 1, ylabel(ax, sprintf('rho thr %.1f\nunwrapped angle (rad)', rho_thr_list(t))); end
            if t == n_thr, xlabel(ax, 'time in snippet (s)'); end
        end
    end
    short = strrep(fly_list{f},'Z:/pablo/lpsp_p2x2_walking/','');
    title(tl, {sprintf('%s (%s), closed-loop pre-ATP: %d trial(s) pooled for ratios; snippet = trial %d from t = %.0f s', short, geno_show, numel(ii), trial_num(i0), T(i0).xf(seg(1))), ...
        'black = heading, cyan = cue (-ft.cue), orange = bump (masked, then smoothed within kept segments), gray = raw bump; shaded = rho <= thr (discarded)'}, 'FontSize', 10)
    exportgraphics(fig, fullfile(export_dir, sprintf('cl_pre_traces_%s.png', regexprep(short,'[/ ]','_'))), 'Resolution', 110)
    fprintf('fly %-22s  bump/cue at thr .3 by window: %s   bump/head: %s\n', short, mat2str(round(squeeze(ratio_c(k,2,:))',2)), mat2str(round(squeeze(ratio_h(k,2,:))',2)));
end

%% summary: per-fly CL-pre ratios across the sweep
fig = figure(41); clf; set(fig,'Color','w','Units','pixels','Position',[40 40 1300 520])
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
cols = lines(n_thr);
for which = 1:2
    ax = nexttile(tl,which); hold(ax,'on')
    if which == 1, Rr = ratio_c; lab = 'bump path / cue path'; ref = 1; else, Rr = ratio_h; lab = 'bump path / heading path'; ref = 0.73; end
    for t = 1:n_thr
        xj = (1:n_win) + (t-2)*0.2;
        Y = squeeze(Rr(:,t,:)); % flies x win
        plot(ax, repmat(xj,size(Y,1),1)', Y', '.', 'Color', [cols(t,:) .4], 'MarkerSize', 9, 'HandleVisibility','off')
        m = median(Y,1,'omitnan'); lo = m - prctile(Y,25,1); hi = prctile(Y,75,1) - m;
        errorbar(ax, xj, m, lo, hi, 'o-', 'Color', cols(t,:), 'MarkerFaceColor','w', 'LineWidth', 1.6, 'DisplayName', sprintf('rho thr %.1f', rho_thr_list(t)))
    end
    yline(ax, ref, '--k', 'HandleVisibility','off')
    xticks(ax, 1:n_win); xticklabels(ax, arrayfun(@(w) win_label(w,1/60), win_list, 'UniformOutput', false)); xlabel(ax, 'Gaussian window')
    ylabel(ax, lab); ylim(ax, [0 2.5]); grid(ax,'on'); legend(ax,'Location','northwest','Box','off')
    title(ax, sprintf('%s, empty>p2x2 closed-loop pre-ATP (n=%d flies); dashed = %.2f', lab, numel(flies), ref))
end
title(tl, 'does the metric read ~1 for controls in closed loop? dot = fly (pooled CL-pre trials), circle = median with IQR')
exportgraphics(fig, fullfile(export_dir,'cl_pre_ratio_summary.png'), 'Resolution', 130)

fprintf('\n== empty>p2x2 CL pre, median over flies [IQR] ==\n%-8s', 'thr\win'); fprintf(' %-22s', win_label(win_list(1),1/60), win_label(win_list(2),1/60), win_label(win_list(3),1/60), win_label(win_list(4),1/60), win_label(win_list(5),1/60)); fprintf('\n');
for t = 1:n_thr
    fprintf('%-8.1f', rho_thr_list(t));
    for w = 1:n_win, v = ratio_c(:,t,w); fprintf(' cue %.2f [%.2f-%.2f]  ', median(v), prctile(v,25), prctile(v,75)); end; fprintf('\n%-8s', '');
    for w = 1:n_win, v = ratio_h(:,t,w); fprintf(' head %.2f [%.2f-%.2f] ', median(v), prctile(v,25), prctile(v,75)); end; fprintf('\n');
end
fprintf('exported per-fly figures + summary to %s\n', export_dir);

%% helpers
function [ms, hs, cs] = mask_smooth(mu, head, cue, keep, win)
    m = mu; m(~keep) = nan; h = head; h(~keep) = nan; c = cue; c(~keep) = nan;
    if win > 1
        ms = smoothdata(m,'gaussian',win); hs = smoothdata(h,'gaussian',win); cs = smoothdata(c,'gaussian',win);
    else
        ms = m; hs = h; cs = c;
    end
    ms(~keep) = nan; hs(~keep) = nan; cs(~keep) = nan;
end

function s = win_label(w, dt)
    if w <= 1, s = 'none'; else, s = sprintf('%.2g s', w*dt); end
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
