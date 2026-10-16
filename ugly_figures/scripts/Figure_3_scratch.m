%% Figure_3_scratch
% Consolidated pickup point for the Figure 3 bump-mobility re-analysis
% (sessions of 2026-10-14 to 2026-10-16). Runs top to bottom on
% .data/lpsp_p2x2_walking_20260728_post0708.mat and reproduces the final
% figures of that work. Read the WHERE THINGS STAND block before changing
% anything in Figure_3_claude.m's mobility row.
%
% ===================== WHERE THINGS STAND =====================
% 1) BUMP EXTRACTION. im.mu / im.rho / im.z were extracted from
%    imgData_reg (NOT imgData_sum_reg): lpsp_p2x2_script_minimal.m:71-82
%    loads imgData_reg.mat and calls process_im on it. process_im applies a
%    10-frame Gaussian to the movie before the wedge PVA, and essentially no
%    smoothing to mu itself (im_win = {10,1}). ft.cue / ft.r_speed carry a
%    30-sample Gaussian from process_ft. So the dataset's baked-in
%    smoothing is UNEQUAL: ~0.95 s on the bump, ~0.5 s on behavior.
% 2) RHO AND MOBILITY ARE COUPLED (flag also in Figure_3_claude.m and
%    lpsp_p2x2_walking_script.m). Across dark trials, trial-mean rho vs
%    Figure 3's path-length mobility: Spearman -0.40 (n=73, p<.001). lpsp
%    whole-trial dark rho drops post-ATP (0.44 -> 0.38, p=.09), empty does
%    not. Path length is a total-variation measure and inflates when mu is
%    noisy. Report any mobility change alongside rho.
% 3) THE UNGATED FIGURE 3 RESULT DEPENDS ON LOW-RHO SAMPLES. With bump
%    samples below a rho threshold NaN'd before smoothing (sections 4-5
%    below), the lpsp dark path-length delta is +0.44 (p=.01 within, .06
%    between) at thr 0 and 0.1, then +0.14 at 0.2 and +0.20 at 0.3 (p>.1,
%    indistinguishable from empty's +0.07..+0.11). The samples that carry
%    it sit in the 0.1-0.2 rho band. Bout-averaged rho itself does NOT
%    change post-ATP in either genotype (|delta| < 0.02, p>.5), so the gate
%    is not removing a general weakening -- it removes a jittery minority.
% 4) NET DISPLACEMENT IS THE BETTER-BEHAVED METRIC. Per bout, end-minus-
%    start of bump vs heading (sqrt-duration-weighted through-origin
%    slope, same bouts as Figure 3): reads the known closed-loop cue gain
%    (0.73, measured in lpsp_p2x2_walking_turns_claude.m) for controls
%    with tight scatter (0.68 [0.62-0.74]) and is almost threshold-
%    independent. |net| shows an lpsp dark increase only WITH the gate:
%    +0.30 at thr .3 (p=.08 within, .10 between, 8/10 flies up), +0.64 at
%    thr .4 (p=.04/.01, but only 8 lpsp flies survive). SIGNED net never
%    changes -> post-ATP bump is less anchored, not higher-gain. Path
%    length and |net| are sensitive to opposite ends of data quality.
% 5) CONTROLS IN THE DARK ARE NOT AT 1. Pre-ATP empty dark path-length
%    slope median 1.07 [0.88-1.46]; fly 18 (20260728) is at 2.3 with a
%    clean bump (fig_bout_diagnostic.m): the bump genuinely sweeps ~2x the
%    fly in darkness, not an artifact. Dark baseline is stable across the
%    two pre-ATP dark trials (first vs later: p=.47 empty, .58 lpsp).
% 6) FLIES WORTH KNOWING ABOUT (controls, from fig_cl_pre_traces_sweep.m):
%    20260727 fly 3 -- diffuse low-rho bump, poor tracking even at high
%    rho; candidate for exclusion. 20260721 fly 2 -- bump mirrors the cue
%    for ~500 s (sign flip; hemisphere-ambiguous PVA or real inversion);
%    passes every rho gate. 20260728 fly 18 -- clean bump that drifts
%    while standing. Several flies flip the sign of cue-bump correlation
%    between closed-loop trials (e.g. 20260721 fly 1: .95, -.80, .64).
% 7) TWO-SIDED tests throughout (signrank within genotype, ranksum
%    between). lpsp>p2x2 is bimodal pre-stim on the signed-net metric
%    (half track tightly in the dark, half ~0), so medians hide structure.
%    The four high-baseline lpsp flies (path slope > 1.8) are all from
%    20260728; a day effect has not been excluded.
% 8) EXAMPLE FLIES for fig_example_flies_pb_heatmaps.m: current picks are
%    empty 20260728 fly 12 and lpsp 20260727 fly 4; the gallery
%    (fig_example_fly_gallery.m) was built so a different pick can be made.
%    Figure 3's lateralized turn analysis lives in data scripts/
%    lpsp_p2x2_walking_turns_claude.m (flips right-stim flies).
%
% Related scripts written in these sessions (all in ugly_figures/scripts
% unless noted): fig_bump_vs_cue_pre_atp, fig_example_flies_pb_heatmaps,
% fig_example_fly_gallery, fig_pathlength_ratio_rho_sweep(_allflies),
% fig_pathlength_ratio_smooth_sweep_allflies, fig_cl_pre_traces_sweep,
% fig_flagged_flies_heatmaps, fig_bout_diagnostic, fig_bout_metrics_rho_sweep,
% data scripts/lpsp_p2x2_walking_turns_claude.
% ===============================================================
%
% What this script does, in order:
%   1 load + labels (fly, trial order, pulse detection, 4 categories)
%   2 per-trial bout cache (Figure 3 bout definition)
%   3 gated bout metrics for a rho-threshold sweep, 3 definitions
%   4 bout rho (ungated, and over kept samples)
%   5 figure A: pre-stim, 3 definitions x thresholds, both genotypes
%   6 figure B: per-fly deltas, 3 definitions + bout rho, by threshold
%   7 figure C: path-length delta + bout rho delta, thresholds 0..0.3
%   8 printed tables for A-C
%   9 figure D: candidate example-fly gallery (one per genotype)
%  10 figure E: the example-fly figure (mean PB image, heatmap snippets + traces, bout mobility)

%% 1) load + labels
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root,'circ_stats'))
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports','figure3_scratch');
if ~isfolder(export_dir), mkdir(export_dir), end
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);
geno_keep = {'empty>p2x2','lpsp>p2x2'}; gc = [0 .45 .74; .85 .33 .1];
cat_names = {'CL pre','dark pre','CL post','dark post'};
fly_geno = cell(n_flies,1); tc = nan(n,1); % trial category 1..4; NaN = pulse trial or fly without pulse
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = geno_keep{1+is_lpsp(idx(1))};
    p = find(has_pulse(idx),1); if isempty(p), continue, end
    pre = idx(1:p-1); post = idx(p+1:end); post = post(~has_pulse(post));
    tc(pre(~is_dark(pre))) = 1; tc(pre(is_dark(pre))) = 2; tc(post(~is_dark(post))) = 3; tc(post(is_dark(post))) = 4;
end
fprintf('%d trials, %d flies (%d empty, %d lpsp); %d usable trials; %d flies without a pulse trial dropped\n', n, n_flies, ...
    sum(strcmp(fly_geno,geno_keep{1})), sum(strcmp(fly_geno,geno_keep{2})), sum(~isnan(tc)), sum(arrayfun(@(f) ~any(has_pulse(fly_ix==f)), 1:n_flies)));

%% 2) per-trial bout cache (Figure 3 bout definition)
sm = 60; turn_thresh = .25; max_gap = 30; min_walk = 30;
C = cell(n,1);
for i = find(~isnan(tc))'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:);
    q = struct();
    q.dt = median(diff(xf));
    q.mu_f = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear','extrap');   % bump, unwrapped, on behavior clock, UNsmoothed
    q.rho_f = interp1(xb, tr.im.rho(:), xf, 'linear','extrap');
    r_s = smoothdata(tr.ft.r_speed(:),'gaussian',sm); r_s(isnan(r_s)) = 0;
    q.r_s = r_s; q.head = cumtrapz(xf, r_s);                             % heading from smoothed yaw (Figure 3)
    w = abs(r_s) > turn_thresh; w = imclose(w, ones(max_gap+1,1)); w = bwareaopen(w, min_walk);
    d = diff([false;w;false]); q.b0 = find(d==1); q.b1 = find(d==-1)-1;
    C{i} = q;
end

%% 3) gated bout metrics: rho threshold sweep x 3 definitions
% mask-then-smooth: bump samples with rho <= thr -> NaN, 1 s Gaussian
% (smoothdata omits NaN), re-NaN; per bout on kept samples only; drop a
% bout if < 50% kept or < 2 kept; per fly x category sqrt(dur)-weighted
% through-origin slope over pooled bouts (>= 15 bouts).
thr_list = [0 0.1 0.2 0.3 0.4 0.5]; n_thr = numel(thr_list);
def_names = {'path-length slope','|net displacement| slope','signed net displacement slope'};
min_bouts = 15; min_kept_frac = .5;
sl = @(x,y,d) (sqrt(d).*x) \ (sqrt(d).*y);
R = nan(n_flies,4,3,n_thr); NB = zeros(n_flies,4,n_thr);
for t = 1:n_thr
    thr = thr_list(t); Bt = cell(n,1);
    for i = find(~isnan(tc))'
        q = C{i}; nb = numel(q.b0); keep = q.rho_f > thr;
        m = q.mu_f; m(~keep) = nan; ms = smoothdata(m,'gaussian',sm); ms(~keep) = nan;
        S = struct('path_b',nan(nb,1),'path_h',nan(nb,1),'net_b',nan(nb,1),'net_h',nan(nb,1),'dur',nan(nb,1));
        for b = 1:nb
            seg = (q.b0(b):q.b1(b))'; k = keep(seg);
            if sum(k) < 2 || mean(k) < min_kept_frac, continue, end
            segk = seg(k); both = k(1:end-1) & k(2:end); dm = abs(diff(ms(seg)));
            S.path_b(b) = sum(dm(both)); S.path_h(b) = sum(abs(q.r_s(segk)))*q.dt;
            S.net_b(b) = ms(segk(end)) - ms(segk(1)); S.net_h(b) = q.head(segk(end)) - q.head(segk(1));
            S.dur(b) = numel(segk)*q.dt;
        end
        ok = ~isnan(S.dur); Bt{i} = structfun(@(v) v(ok), S, 'UniformOutput', false);
    end
    for f = 1:n_flies
        for cc = 1:4
            ii = find(fly_ix==f & tc==cc); if isempty(ii), continue, end
            S = [Bt{ii}]; pb = vertcat(S.path_b); ph = vertcat(S.path_h); nbv = vertcat(S.net_b); nh = vertcat(S.net_h); du = vertcat(S.dur);
            NB(f,cc,t) = numel(du); if numel(du) < min_bouts, continue, end
            R(f,cc,1,t) = sl(ph,pb,du); R(f,cc,2,t) = sl(abs(nh),abs(nbv),du); R(f,cc,3,t) = sl(nh,nbv,du);
        end
    end
end
D = cat(2, R(:,3,:,:)-R(:,1,:,:), R(:,4,:,:)-R(:,2,:,:)); % flies x [CL dark] x def x thr : post - pre
dn = {'CL (post - pre)','dark (post - pre)'};

%% 4) bout rho: ungated mean over bout samples, and mean over kept samples per threshold
rho_u = nan(n_flies,4); rho_g = nan(n_flies,4,n_thr);
for f = 1:n_flies
    for cc = 1:4
        ii = find(fly_ix==f & tc==cc); if isempty(ii), continue, end
        acc = [];
        for i = ii', q = C{i}; for b = 1:numel(q.b0), acc = [acc; q.rho_f(q.b0(b):q.b1(b))]; end, end %#ok<AGROW>
        if isempty(acc), continue, end
        rho_u(f,cc) = mean(acc);
        for t = 1:n_thr, rho_g(f,cc,t) = mean(acc(acc > thr_list(t))); end
    end
end
Drho_u = [rho_u(:,3)-rho_u(:,1), rho_u(:,4)-rho_u(:,2)];
Drho_g = cat(2, rho_g(:,3,:)-rho_g(:,1,:), rho_g(:,4,:)-rho_g(:,2,:));

%% 5) figure A: pre-stim, both genotypes, 3 definitions, rows = threshold
used = any(~isnan(R(:,:,1,1)),2); xo = [-.12 .12];
fig = figure(90); clf; set(fig,'Color','w','Units','pixels','Position',[30 30 1300 210*n_thr])
tl = tiledlayout(n_thr,3,'TileSpacing','compact','Padding','compact');
for t = 1:n_thr
    for k = 1:3
        ax = nexttile(tl,(t-1)*3+k); hold(ax,'on')
        for gi = 1:2
            sel = find(used & strcmp(fly_geno,geno_keep{gi})); Y = reshape(R(sel,1:2,k,t),numel(sel),2); X = (1:2)+xo(gi);
            for f = 1:size(Y,1), plot(ax, X, Y(f,:), '-o', 'Color',[gc(gi,:) .35], 'MarkerFaceColor',gc(gi,:), 'MarkerSize',3.5, 'HandleVisibility','off'); end
            m = median(Y,1,'omitnan'); lo = m-prctile(Y,25,1); hi = prctile(Y,75,1)-m;
            errorbar(ax, X, m, lo, hi, 'o', 'Color',gc(gi,:), 'MarkerFaceColor','w', 'LineWidth',1.6, 'CapSize',7, 'DisplayName',sprintf('%s (n=%d)',geno_keep{gi},sum(any(~isnan(Y),2))))
        end
        yline(ax,1,':k','HandleVisibility','off'); yline(ax,.73,'--','Color',[.5 .5 .5],'HandleVisibility','off'); if k==3, yline(ax,0,'-k','HandleVisibility','off'); end
        xticks(ax,1:2); xlim(ax,[.5 2.5]); ylim(ax,[-.5 3.2]); grid(ax,'on'); set(ax,'FontSize',8)
        if t==n_thr, xticklabels(ax,cat_names(1:2)); else, xticklabels(ax,[]); end
        if k==1, ylabel(ax,sprintf('rho thr %.1f\nbump / heading per bout',thr_list(t))); end
        if t==1, title(ax,def_names{k}); end
        if t==1 && k==1, legend(ax,'Location','northwest','Box','off','Interpreter','none','FontSize',7); end
    end
end
title(tl,{'A. PRE-stim bout mobility, three definitions, rows = rho threshold (bump samples <= thr NaN''d before smoothing)','dotted = 1, dashed = 0.73 cue gain; dot = fly, circle = median with IQR'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'A_prestim_definitions_rho_sweep.png'), 'Resolution',130)

%% 6) figure B: per-fly deltas, 3 definitions + bout rho (ungated), rows = threshold
fig = figure(91); clf; set(fig,'Color','w','Units','pixels','Position',[30 30 1650 215*n_thr])
tl = tiledlayout(n_thr,4,'TileSpacing','compact','Padding','compact');
for t = 1:n_thr
    for k = 1:4
        ax = nexttile(tl,(t-1)*4+k); hold(ax,'on')
        for cc = 1:2
            ab = cell(1,2);
            for gi = 1:2
                if k <= 3, v = D(strcmp(fly_geno,geno_keep{gi}),cc,k,t); else, v = Drho_u(strcmp(fly_geno,geno_keep{gi}),cc); end
                v = v(~isnan(v)); ab{gi} = v; x = cc + xo(gi);
                scatter(ax, x + (rand(numel(v),1)-.5)*.08, v, 14, gc(gi,:), 'filled', 'MarkerFaceAlpha',.55, 'HandleVisibility','off')
                errorbar(ax, x, median(v), median(v)-prctile(v,25), prctile(v,75)-median(v), 'o', 'Color',gc(gi,:), 'MarkerFaceColor','w', 'LineWidth',1.5, 'CapSize',5, 'DisplayName',geno_keep{gi})
                ytxt = 1.75; if k == 4, ytxt = .22; end
                text(ax, x, ytxt, sprintf('p=%.2f',signrank(v)), 'Color',gc(gi,:), 'FontSize',6.5, 'HorizontalAlignment','center')
            end
            ytxt = 2.0; if k == 4, ytxt = .26; end
            text(ax, cc, ytxt, sprintf('between p=%.2f',ranksum(ab{1},ab{2})), 'FontSize',6.5, 'HorizontalAlignment','center')
        end
        yline(ax,0,'-k','HandleVisibility','off'); xticks(ax,1:2); xlim(ax,[.5 2.5]); grid(ax,'on'); set(ax,'FontSize',8)
        if k <= 3, ylim(ax,[-1.3 2.2]); else, ylim(ax,[-.3 .3]); end
        if t==n_thr, xticklabels(ax,dn); else, xticklabels(ax,[]); end
        if k==1, ylabel(ax,sprintf('rho thr %.1f\npost - pre (per fly)',thr_list(t))); end
        if t==1, if k<=3, title(ax,def_names{k}); else, title(ax,'mean rho over walking-bout samples (ungated)'); end, end
        if t==1 && k==1, legend(ax,'Location','southwest','Box','off','Interpreter','none','FontSize',7); end
    end
end
title(tl,{'B. per-fly delta (post - pre ATP), rows = rho threshold, cols = mobility definition + bout rho','dot = fly, circle = median with IQR; two-sided signrank vs 0 per genotype, two-sided ranksum between. Rho column ungated (identical across rows).'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'B_delta_definitions_rho_sweep.png'), 'Resolution',130)

%% 7) figure C: path-length delta + bout rho over kept samples, thresholds 0..0.3
thr_sel = [0 0.1 0.2 0.3]; ti = arrayfun(@(x) find(abs(thr_list-x)<1e-9), thr_sel);
fig = figure(92); clf; set(fig,'Color','w','Units','pixels','Position',[40 40 900 215*numel(thr_sel)])
tl = tiledlayout(numel(thr_sel),2,'TileSpacing','compact','Padding','compact');
for t = 1:numel(thr_sel)
    for k = 1:2
        ax = nexttile(tl,(t-1)*2+k); hold(ax,'on')
        for cc = 1:2
            ab = cell(1,2);
            for gi = 1:2
                if k == 1, v = D(strcmp(fly_geno,geno_keep{gi}),cc,1,ti(t)); else, v = Drho_g(strcmp(fly_geno,geno_keep{gi}),cc,ti(t)); end
                v = v(~isnan(v)); ab{gi} = v; x = cc + xo(gi);
                scatter(ax, x + (rand(numel(v),1)-.5)*.08, v, 18, gc(gi,:), 'filled', 'MarkerFaceAlpha',.6, 'HandleVisibility','off')
                errorbar(ax, x, median(v), median(v)-prctile(v,25), prctile(v,75)-median(v), 'o', 'Color',gc(gi,:), 'MarkerFaceColor','w', 'LineWidth',1.6, 'CapSize',6, 'DisplayName',sprintf('%s (n=%d)',geno_keep{gi},numel(v)))
                ytxt = 1.55; if k == 2, ytxt = .22; end
                text(ax, x, ytxt, sprintf('p=%.2f',signrank(v)), 'Color',gc(gi,:), 'FontSize',7, 'HorizontalAlignment','center')
            end
            ytxt = 1.8; if k == 2, ytxt = .26; end
            text(ax, cc, ytxt, sprintf('between p=%.2f',ranksum(ab{1},ab{2})), 'FontSize',7, 'HorizontalAlignment','center')
        end
        yline(ax,0,'-k','HandleVisibility','off'); xticks(ax,1:2); xlim(ax,[.5 2.5]); grid(ax,'on'); set(ax,'FontSize',8)
        if k == 1, ylim(ax,[-1 2]); ylabel(ax,sprintf('rho thr %.1f\npost - pre (per fly)',thr_sel(t))); else, ylim(ax,[-.3 .3]); end
        if t == numel(thr_sel), xticklabels(ax,dn); else, xticklabels(ax,[]); end
        if t == 1, if k == 1, title(ax,'path-length slope (Figure 3 metric), bump samples <= thr NaN''d'); else, title(ax,'mean rho over walking-bout samples with rho > thr'); end, end
        if t == 1 && k == 1, legend(ax,'Location','southwest','Box','off','Interpreter','none','FontSize',7); end
    end
end
title(tl,{'C. per-fly delta (post - pre ATP): path-length mobility and bout rho, rows = rho threshold','dot = fly, circle = median with IQR; two-sided signrank vs 0 per genotype, two-sided ranksum between'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'C_pathlength_and_rho_delta.png'), 'Resolution',140)

%% 8) tables
fprintf('\n== per-fly delta (post - pre), median (n) | signrank p | ranksum p between ==\n');
for k = 1:3
    fprintf('\n%s\n%-5s %-45s %-45s\n', def_names{k}, 'thr', 'CL: empty | lpsp | between', 'dark: empty | lpsp | between');
    for t = 1:n_thr
        fprintf('%-5.1f', thr_list(t));
        for cc = 1:2
            a = D(strcmp(fly_geno,geno_keep{1}),cc,k,t); a = a(~isnan(a)); b = D(strcmp(fly_geno,geno_keep{2}),cc,k,t); b = b(~isnan(b));
            fprintf(' %+.2f(%2d) p%.2f | %+.2f(%2d) p%.2f | p%.2f      ', median(a), numel(a), signrank(a), median(b), numel(b), signrank(b), ranksum(a,b));
        end
        fprintf('\n');
    end
end
fprintf('\nbout rho delta (ungated): ');
for cc = 1:2, for gi = 1:2, v = Drho_u(strcmp(fly_geno,geno_keep{gi}),cc); v = v(~isnan(v)); fprintf('%s %s %+.3f (p%.2f)  ', dn{cc}(1:4), geno_keep{gi}, median(v), signrank(v)); end, end; fprintf('\n');
fprintf('\nexported A-C to %s\n', export_dir);

%% 9) example-fly gallery: score every fly, one panel per fly per genotype
% Scores: r_CL = circ corr of -cue vs mu on pre-ATP closed-loop trials
% (rho > .2 and cue moving > .2 rad/s), mean over trials; mob = ungated
% path-length slope dark pre -> post (R(:,:,1,1)); rho = mean rho over
% dark trials pre -> post. Panels sorted by r_CL. Pick from these.
track_r = nan(n,1); trial_rho = nan(n,1);
for i = 1:n
    trial_rho(i) = mean(all_data(i).im.rho);
    if ~is_dark(i), track_r(i) = trial_tracking_r(all_data(i), .2, .2); end
end
G = struct('f',{},'short',{},'geno',{},'r',{},'r_list',{},'pre_idx',{},'post_idx',{},'cl_pre_idx',{},'mob_pre',{},'mob_post',{},'rho_pre',{},'rho_post',{});
for f = 1:n_flies
    if ~any(tc(fly_ix==f) >= 1), continue, end
    g = struct(); g.f = f; g.short = strrep(fly_list{f},'Z:/pablo/lpsp_p2x2_walking/',''); g.geno = fly_geno{f};
    g.cl_pre_idx = find(fly_ix==f & tc==1); g.pre_idx = find(fly_ix==f & tc==2); g.post_idx = find(fly_ix==f & tc==4);
    g.r_list = track_r(g.cl_pre_idx); g.r = mean(g.r_list,'omitnan');
    g.mob_pre = R(f,2,1,1); g.mob_post = R(f,4,1,1); g.rho_pre = mean(trial_rho(g.pre_idx)); g.rho_post = mean(trial_rho(g.post_idx));
    G(end+1) = g; %#ok<AGROW>
end
for gi = 1:2
    Gg = G(strcmp({G.geno},geno_keep{gi})); [~,o] = sort([Gg.r],'descend','MissingPlacement','last'); Gg = Gg(o);
    nG = numel(Gg); n_col = 4; n_row = ceil(nG/n_col);
    fig = figure(93+gi); clf; set(fig,'Color','w','Units','pixels','Position',[20 20 1700 320*n_row])
    tl = tiledlayout(n_row,n_col,'TileSpacing','compact','Padding','compact');
    for k = 1:nG
        g = Gg(k); ax = nexttile(tl,k); hold(ax,'on')
        sets = {g.pre_idx, g.post_idx}; cols = {[.35 .35 .35],[.85 .2 .2]}; mobs = [g.mob_pre, g.mob_post]; xmax = 1;
        for c = 1:2
            [bc, bm, bd] = bouts_from_cache(C, sets{c}, sm);
            if isempty(bc), continue, end
            scatter(ax, bc, bm, 10+3*sqrt(bd), cols{c}, 'filled', 'MarkerFaceAlpha', .5); xmax = max([xmax; bc]);
            if ~isnan(mobs(c)), plot(ax, [0 xmax], mobs(c)*[0 xmax], '-', 'Color', cols{c}, 'LineWidth', 1.5); end
        end
        plot(ax, [0 xmax], [0 xmax], 'k:'); xlim(ax,[0 xmax*1.05]); ylim(ax,[0 max(xmax*1.05,1)]); axis(ax,'square'); set(ax,'FontSize',7)
        title(ax, sprintf('#%d %s\nr_{CL}=%.2f %s\nmob %.2f -> %.2f | rho %.2f -> %.2f', k, g.short, g.r, mat2str(round(g.r_list',2)), g.mob_pre, g.mob_post, g.rho_pre, g.rho_post), 'FontSize', 7.5)
        if mod(k-1,n_col)==0, ylabel(ax,'bump path / bout (rad)'); end
        if k > nG-n_col, xlabel(ax,'heading path / bout (rad)'); end
    end
    title(tl, sprintf('%s candidate example flies, sorted by closed-loop tracking r. gray = dark pre bouts, red = dark post bouts, lines = pooled sqrt(dur)-weighted slopes (ungated path length), dotted = 1:1', geno_keep{gi}), 'Interpreter','none','FontSize',10)
    exportgraphics(fig, fullfile(export_dir, sprintf('D_gallery_%s.png', strrep(geno_keep{gi},'>','_'))), 'Resolution', 120)
end

%% 10) example-fly figure (the one that goes side by side with Figure 3)
% Per fly: mean PB image from imgData_reg (the movie the bump was extracted
% from) with mask + centroids; 30 s z-score heatmap snippets for one
% closed-loop trial, one dark pre and one dark post trial, each with a
% separate heading + bump trace panel; and the dark pre/post bout scatter.
% Snippet = among the 30 s windows with >= 90% of the best walking-bout
% coverage, the one with median rho. Override via ex(e).snip_t0.
% CURRENT PICKS (see gallery D): empty 20260728/fly 12 (r_CL .80/.94,
% dark mob .90 -> 1.11), lpsp 20260727/fly 4 (r_CL .88/.88 on trials
% 30/31, dark mob 1.05 -> 3.82; post trial 40 fly barely turned).
clear ex
ex(1).tag = '20260727/fly 4';  ex(1).label = 'lpsp>p2x2';  ex(1).cl = 30; ex(1).dark_pre = 34; ex(1).dark_post = 40;
ex(2).tag = '20260728/fly 12'; ex(2).label = 'empty>p2x2'; ex(2).cl = 98; ex(2).dark_pre = 99; ex(2).dark_post = 103;
snip_len = 30;
for e = 1:2
    f = find(contains(fly_list, ex(e).tag)); assert(isscalar(f), 'fly tag %s not found', ex(e).tag); ex(e).f = f;
    ii = find(fly_ix==f);
    ex(e).i_cl = ii(trial_num(ii)==ex(e).cl); ex(e).i_pre = ii(trial_num(ii)==ex(e).dark_pre); ex(e).i_post = ii(trial_num(ii)==ex(e).dark_post);
    assert(~is_dark(ex(e).i_cl) && is_dark(ex(e).i_pre) && is_dark(ex(e).i_post), 'lighting mismatch for %s', ex(e).tag)
    reg_file = fullfile(all_data(ex(e).i_cl).meta, 'imgData_reg.mat');
    if isfile(reg_file), tmp = load(reg_file,'imgData_reg'); ex(e).mean_img = mean(tmp.imgData_reg,3); clear tmp
    else, warning('%s not reachable (Z: offline?) -- mean PB image panel left blank', reg_file); ex(e).mean_img = []; end
    trial_ii = [ex(e).i_cl, ex(e).i_pre, ex(e).i_post]; ex(e).snip_t0 = nan(1,3);
    for c = 1:3
        i = trial_ii(c); tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf)); w = round(snip_len/dt);
        r_s = smoothdata(tr.ft.r_speed(:),'gaussian',sm); iw = abs(r_s) > turn_thresh; iw = imclose(iw, ones(max_gap+1,1)); iw = bwareaopen(iw, min_walk);
        cover = movmean(double(iw), w, 'Endpoints','discard'); rho_w = movmean(interp1(xb, tr.im.rho, xf, 'linear','extrap'), w, 'Endpoints','discard');
        cand = find(cover >= .9*max(cover)); [~,k] = min(abs(rho_w(cand) - median(rho_w(cand)))); ex(e).snip_t0(c) = xf(cand(k));
    end
end
fig = figure(96); clf; set(fig,'Color','w','Units','pixels','Position',[40 40 1900 900])
tl = tiledlayout(4,5,'TileSpacing','compact','Padding','compact'); heat_clim = [-1.5 3];
for e = 1:2
    r0 = (e-1)*2; i_cl = ex(e).i_cl;
    ax = nexttile(tl, r0*5+1, [2 1]); hold(ax,'on')
    if ~isempty(ex(e).mean_img)
        im = ex(e).mean_img; imagesc(ax, im, [prctile(im(:),1) prctile(im(:),99.5)]); colormap(ax, bone); axis(ax,'image')
        B = bwboundaries(all_data(i_cl).im.mask); for b = 1:numel(B), plot(ax, B{b}(:,2), B{b}(:,1), 'y-', 'LineWidth', .7); end
        cen = all_data(i_cl).im.centroids; plot(ax, cen(:,2), cen(:,1), 'r.', 'MarkerSize', 7)
    end
    set(ax,'XTick',[],'YTick',[]); title(ax, sprintf('%s  (%s)\nmean imgData_reg, trial %d', ex(e).label, strrep(ex(e).tag,'/',' '), ex(e).cl), 'FontSize', 9, 'Interpreter','none')
    trial_ii = [ex(e).i_cl, ex(e).i_pre, ex(e).i_post];
    ttl = {sprintf('closed loop, trial %d', ex(e).cl), sprintf('dark PRE atp, trial %d', ex(e).dark_pre), sprintf('dark POST atp, trial %d', ex(e).dark_post)};
    for c = 1:3
        i = trial_ii(c); tr = all_data(i); t0 = ex(e).snip_t0(c); t1 = t0 + snip_len; xb = tr.ft.xb(:); xf = tr.ft.xf(:); dt = median(diff(xf));
        kb = xb >= t0 & xb <= t1; kf = xf >= t0 & xf <= t1;
        ax = nexttile(tl, r0*5+1+c); imagesc(ax, xb(kb)-t0, 1:32, tr.im.z(:,kb), heat_clim); colormap(ax, parula); hold(ax,'on'); yline(ax, 16.5, 'w:')
        set(ax,'YDir','normal','YTick',[1 16 17 32],'XTickLabel',[],'FontSize',8); xlim(ax,[0 snip_len]); if c == 1, ylabel(ax,'PB glomerulus'); end
        [bc, bm, bd] = bouts_from_cache(C, i, sm); slv = nan; if numel(bc) >= 5, slv = sl(bc, bm, bd); end
        title(ax, sprintf('%s  (mobility %.2f, %d bouts)\n%d s snippet from t = %.0f s, mean rho %.2f', ttl{c}, slv, numel(bc), snip_len, t0, mean(tr.im.rho(kb))), 'FontSize', 8)
        ax = nexttile(tl, (r0+1)*5+1+c); hold(ax,'on')
        mu_u = unwrap(tr.im.mu(kb)); mu_u = mu_u - mu_u(1);
        if is_dark(i), head = cumsum(tr.ft.r_speed(kf))*dt; head_lab = 'heading (\int r\_speed)'; else, head = unwrap(-tr.ft.cue(kf)); head = head - head(1); head_lab = 'heading (-cue)'; end
        plot(ax, xf(kf)-t0, head, '-', 'Color', [0 .6 .8], 'LineWidth', 1.2); plot(ax, xb(kb)-t0, mu_u, 'k-', 'LineWidth', 1.2)
        xlim(ax,[0 snip_len]); xlabel(ax,'time in snippet (s)'); yl = max([abs(head); abs(mu_u); pi]); ylim(ax, [-yl yl]*1.05)
        step = pi*max(1, ceil(yl/(2*pi))); kk = (-floor(yl/step):floor(yl/step)); yticks(ax, kk*step); yticklabels(ax, arrayfun(@(k) pi_label(k*step/pi), kk, 'UniformOutput', false))
        if c == 1, ylabel(ax,'angle (rad, re-zeroed)'); end; legend(ax, {head_lab, 'bump (im.mu)'}, 'Location','best','FontSize',7,'Box','off'); set(ax,'FontSize',8)
    end
    ax = nexttile(tl, r0*5+5, [2 1]); hold(ax,'on'); cols = {[.3 .3 .3], [.85 .2 .2]}; lab = {'dark PRE','dark POST'};
    for c = 1:2
        i = trial_ii(c+1); [bc, bm, bd] = bouts_from_cache(C, i, sm);
        scatter(ax, bc, bm, 20 + 4*sqrt(bd), cols{c}, 'filled', 'MarkerFaceAlpha', .6); slv = sl(bc, bm, bd);
        xl = [0 max([bc; 1])]; plot(ax, xl, slv*xl, '-', 'Color', cols{c}, 'LineWidth', 1.5); lab{c} = sprintf('%s (trial %d): slope %.2f', lab{c}, trial_num(i), slv);
    end
    plot(ax, [0 max(xlim(ax))], [0 max(xlim(ax))], 'k:'); axis(ax,'square'); xlabel(ax,'heading path length per bout (rad)'); ylabel(ax,'bump path length per bout (rad)')
    legend(ax, lab, 'Location','northwest','FontSize',7,'Box','off'); title(ax, sprintf('bump mobility, dark trials\n(fly pooled: pre %.2f, post %.2f)', R(ex(e).f,2,1,1), R(ex(e).f,4,1,1)), 'FontSize', 9)
end
cb = colorbar(nexttile(tl,4)); cb.Label.String = 'z-scored dF/F'; cb.Layout.Tile = 'east';
title(tl, sprintf('E. example flies: PB anatomy, %d s z-scored activity snippets with heading + bump traces below, and dark-trial bump mobility before vs after ATP', snip_len), 'FontSize', 11)
exportgraphics(fig, fullfile(export_dir,'E_example_flies_pb_heatmaps.png'), 'Resolution', 150)
fprintf('exported gallery (D) and example figure (E) to %s\n', export_dir);

%% helpers
function [bc, bm, bd] = bouts_from_cache(C, ii, sm)
    % ungated per-bout heading path, bump path and duration (Figure 3 metric) for trials ii
    bc = []; bm = []; bd = [];
    for i = ii(:)'
        q = C{i}; if isempty(q), continue, end
        ms = smoothdata(q.mu_f, 'gaussian', sm);
        for b = 1:numel(q.b0)
            seg = q.b0(b):q.b1(b);
            bc(end+1,1) = sum(abs(q.r_s(seg)))*q.dt; bm(end+1,1) = sum(abs(diff(ms(seg)))); bd(end+1,1) = numel(seg)*q.dt; %#ok<AGROW>
        end
    end
end

function r = trial_tracking_r(trial, rho_thr, cue_speed_thr)
    xb = trial.ft.xb; xf = trial.ft.xf;
    cue_u = interp1(xf, unwrap(-trial.ft.cue), xb, 'linear', 'extrap');
    cue_sp = abs(gradient(smoothdata(cue_u,'movmean',5), xb));
    good = trial.im.rho > rho_thr & cue_sp > cue_speed_thr;
    r = nan; if sum(good) > 100, r = circ_corrcc(mod(cue_u(good),2*pi), mod(trial.im.mu(good),2*pi)); end
end

function lab = pi_label(m)
    if m == 0, lab = '0'; elseif m == 1, lab = '\pi'; elseif m == -1, lab = '-\pi'; else, lab = [num2str(m) '\pi']; end
end
function [fid, tnum] = walk_meta_fly_trial(meta_path)
    parts = strsplit(meta_path,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    % Figure_3_claude.m / walking script section 3 detector
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
