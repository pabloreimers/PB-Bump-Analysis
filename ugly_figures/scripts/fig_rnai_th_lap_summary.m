%% fig_rnai_th_lap_summary
% One-page summary of what predicts laps in high-gain LPsP > TH-RNAi flies
% (the exploration in fig_rnai_th_lap_predictors.m, which has the full set of
% tests). Same per-turn table: closed loop, included trials, turns >=
% min_turn_deg with a solid bump (mean PVA strength >= rho_min over the
% turn, 10th percentile >= rho_p10_min over the whole scored interval), in the
% high-gain flies that overshoot at least once; outcome = OVERSHOOT, the bump
% rotating >= overshoot_deg past the fly from turn onset to post_s after the
% turn.
%
% Panels: (a) an example lap; (b) the best combination of predictors -- a
% mixed logistic model (fly random intercept) chosen by cross-validated
% subset search over turn speed, duration, size, forward walking before the
% turn, time in trial, bump strength, bar position and bump position --
% odds ratio per SD; (c) how concentrated laps are under that model:
% cumulative share of laps vs share of turns ranked by CROSS-VALIDATED
% predicted risk (10-fold, folds within fly), against fly identity alone and
% chance; (d) lap rate in simple turn categories (size x speed, split by
% where the bar was); (e) where the bar is during lap turns vs other turns
% (ft.cue, 0 = behind the fly); (f) clustering in time vs a within-trial
% shuffle. No interaction term improved cross-validated prediction, so the
% model is additive.
% Needs all_data (loaded if absent), the integrative-gain cache and the turn
% cache. Run from the repo root.

%% parameters
min_turn_deg = 30; overshoot_deg = 90; post_s = 1.0; rho_min = 0.4; gain_split = 0.85;
rho_p10_min = 0.3; % 10th percentile of PVA strength over the whole scored interval
example = {'20250905','5','12'}; % date, fly, trial of the example lap
n_perm = 1000; n_fold = 10; n_bar_bins = 8;
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
addpath(fullfile(pwd,'circ_stats'));
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
TC = load(fullfile(data_dir,'lpsp_rnai_turns_cache.mat'),'turns');
if ~exist('all_data','var')
    tmp = load(fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
g_fly = nan(max(IG.fly_num),1);
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
fly_geno = cell(max(IG.fly_num),1); for f = 1:max(IG.fly_num), fly_geno{f} = IG.genotype{find(IG.fly_num==f,1)}; end
out_dir = fullfile('ugly_figures','rnai_turns','laps_stalls');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
col = [0.60,0.15,0.00]; rng(1)

%% per-turn table (as fig_rnai_th_lap_predictors.m)
hi_fl = find(strcmp(fly_geno,'lpsp>th') & g_fly > gain_split);
tr_use = find(~IG.is_dark & IG.inc & ismember(IG.fly_num,hi_fl));
tt = [TC.turns.trial]';
rows = zeros(0,11); bar_traj = {};
for i = tr_use(:)'
    tr = all_data(i); xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    mu_u = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'));
    fs = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',max(1,round(0.5/median(diff(xf)))));
    a1 = tr.im.alpha(1);
    for k = find(tt==i)'
        T = TC.turns(k); t0 = T.onset_s; t1 = t0 + T.dur_s;
        b0 = find(xb >= t0,1); b1 = find(xb <= t1 + post_s,1,'last'); bt = find(xb <= t1,1,'last');
        f0 = find(xf >= t0,1); f1 = find(xf <= t1,1,'last');
        if isempty(b0) || isempty(b1) || isempty(bt) || bt < b0 || isempty(f1) || f1 <= f0, continue, end
        dfly = rad2deg(hd(f1) - hd(f0)); if abs(dfly) < min_turn_deg, continue, end
        c_in = tr.ft.cue(f0:f1); c_in = c_in(~isnan(c_in)); c_in = c_in(:);
        rows(end+1,:) = [i, IG.fly_num(i), t0, T.dur_s, dfly, rad2deg(mu_u(b1) - mu_u(b0)), ...
            mean(fs(xf >= t0-2 & xf < t0)), mean(rho(b0:bt)), mod(sign(dfly)*(mu(b0) - a1), 2*pi), abs(dfly)/T.dur_s, prctile(rho(b0:b1),10)]; %#ok<AGROW>
        bar_traj{end+1,1} = angle(exp(1i*c_in)); %#ok<AGROW>
    end
end
D = array2table(rows,'VariableNames',{'trial','fly','t0','dur','dfly','dbump','fwd_pre','rho_turn','rel','speed','rho_p10'});
D.bar_traj = bar_traj;
D.os = sign(D.dfly).*D.dbump - abs(D.dfly) >= overshoot_deg;
D = D(D.rho_turn >= rho_min & D.rho_p10 >= rho_p10_min,:);
n_os = accumarray(D.fly,D.os,[max(IG.fly_num),1]); D = D(n_os(D.fly) >= 1,:);
D.bar_behind = cellfun(@(a) mean(cos(a)), D.bar_traj); D.bar_behind(isnan(D.bar_behind)) = 0; % +1 = bar behind all turn, -1 = in front
y = D.os; n = height(D);
fprintf('%d turns in %d flies, %d overshoots (%.1f%%)\n',n,numel(unique(D.fly)),sum(y),100*mean(y));

%% (b) best model: mixed logistic, fly random intercept
Zs = @(v) (v - mean(v,'omitnan'))./std(v,'omitnan');
M = table(y,categorical(D.fly),Zs(log(abs(D.dfly))),Zs(log(D.speed)),Zs(D.rho_turn),Zs(D.bar_behind),Zs(D.t0),cos(D.rel),sin(D.rel), ...
    'VariableNames',{'y','fly','turn_size','turn_speed','bump_strength','bar_behind','time_in_trial','bpos_cos','bpos_sin'});
glme = fitglme(M,'y ~ turn_size + turn_speed + bump_strength + bar_behind + time_in_trial + bpos_cos + bpos_sin + (1|fly)','Distribution','Binomial','FitMethod','Laplace');
disp(glme.Coefficients(:,{'Name','Estimate','SE','pValue'}))

%% (c) cross-validated risk: best model vs fly identity alone
folds = zeros(n,1); for f = unique(D.fly)', idx = find(D.fly==f); folds(idx(randperm(numel(idx)))) = mod(0:numel(idx)-1,n_fold)+1; end
frm_full = 'y ~ fly + turn_size + turn_speed + bump_strength + bar_behind + time_in_trial + bpos_cos + bpos_sin';
p_full = nan(n,1); p_fly = nan(n,1);
warning('off','stats:glmfit:IterationLimit'); warning('off','stats:glmfit:PerfectSeparation')
for k = 1:n_fold
    te = folds==k;
    p_full(te) = predict(fitglm(M(~te,:),frm_full,'Distribution','binomial'),M(te,:));
    p_fly(te) = predict(fitglm(M(~te,:),'y ~ fly','Distribution','binomial'),M(te,:));
end
[~,~,~,auc_full] = perfcurve(y,p_full,true); [~,~,~,auc_fly] = perfcurve(y,p_fly,true);
[~,o_full] = sort(p_full,'descend'); [~,o_fly] = sort(p_fly,'descend');
g_full = cumsum(y(o_full))/sum(y); g_flyonly = cumsum(y(o_fly))/sum(y); xs = (1:n)'/n;
for q = [0.1 0.2], fprintf('top %.0f%% of turns (cross-validated risk): %.0f%% of laps (fly only: %.0f%%)\n',100*q,100*g_full(round(q*n)),100*g_flyonly(round(q*n))); end

%% (d) simple categories
big = abs(D.dfly) > median(abs(D.dfly)); slow = D.speed < median(D.speed); beh = D.bar_behind > 0;
cats = {'large, slow',big & slow; 'small, slow',~big & slow; 'large, fast',big & ~slow; 'small, fast',~big & ~slow};

%% (e) bar occupancy during lap vs other turns, within-fly shuffle
edges_b = linspace(-pi,pi,n_bar_bins+1); occ = zeros(n,n_bar_bins);
for k = 1:n, bb = discretize(D.bar_traj{k},edges_b); bb = bb(~isnan(bb)); if ~isempty(bb), occ(k,:) = accumarray(bb(:),1,[n_bar_bins,1])'/numel(bb); end, end
occ_ratio = @(lab) mean(occ(lab,:),1)./mean(occ(~lab,:),1);
oc_obs = occ_ratio(y); oc_null = zeros(n_perm,n_bar_bins);
for p = 1:n_perm, lab = y; for f = unique(D.fly)', idx = find(D.fly==f); lab(idx) = lab(idx(randperm(numel(idx)))); end, oc_null(p,:) = occ_ratio(lab); end
p_occ = mean(sum(log(oc_null).^2,2) >= sum(log(oc_obs).^2));

%% (f) clustering: rate vs time since the previous overshoot (same trial), within-trial shuffle
G = sortrows(D,{'trial','t0'}); ce = [0,5,15,30,60,120,Inf]; nc = numel(ce);
since = since_last(G); cb = discretize(since,ce); cb(isinf(since)) = nc;
rate_c = accumarray(cb,G.os,[nc,1],@mean); null_c = zeros(n_perm,nc); nul = zeros(n_perm,1);
for p = 1:n_perm
    H = G; for i = unique(H.trial)', idx = find(H.trial==i); H.os(idx) = H.os(idx(randperm(numel(idx)))); end
    s = since_last(H); c = discretize(s,ce); c(isinf(s)) = nc; null_c(p,:) = accumarray(c,H.os,[nc,1],@mean)';
    nul(p) = mean(H.os(s <= 15)) - mean(H.os(s > 15));
end
p_clust = mean(nul >= mean(G.os(since <= 15)) - mean(G.os(since > 15)));

%% figure
figure('Color','w','Position',[40,40,1800,980]); clf
tl = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
% (a) example lap
k_ex = find(contains({all_data.meta}',['\' example{1} '\fly ' example{2}]) & contains({all_data.meta}',[example{1} '-' example{3} '_']),1);
tr = all_data(k_ex); xf = tr.ft.xf(:); nb = size(tr.im.f,2);
if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
[x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); mu = cart2pol(mean(x_,2),mean(y_,2)); mu_u = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/median(diff(xb)))));
cu = -tr.ft.cue(:); ok = ~isnan(cu); cu(ok) = unwrap(cu(ok)); cu = fillmissing(cu,'linear','EndValues','nearest');
w = [85,111]; ib = xb >= w(1) & xb <= w(2); jf = xf >= w(1) & xf <= w(2);
nexttile(tl,1,[1,2]); hold on
imagesc(xb(ib)-w(1),1:size(f_z,1),f_z(:,ib)); clim([-1.5,3]); colormap(gca,flipud(gray)); axis tight; set(gca,'YDir','reverse'); yticks([])
n_half = size(f_z,1)/2; al = tr.im.alpha(1:n_half); dgl = 2*pi/n_half; cb = interp1(xf,cu,xb);
c0 = circ_mean(circ_dist(mu(ib),cb(ib)));                 % cue drawn at the window's mean bump-cue offset
for h = 1:2
    rb = mod(mu(ib) - al(1) + dgl/2,2*pi)/dgl + 0.5; rc = mod(cb(ib) + c0 - al(1) + dgl/2,2*pi)/dgl + 0.5; tb = xb(ib)-w(1);
    rb([false; abs(diff(rb)) > n_half/2]) = NaN; rc([false; abs(diff(rc)) > n_half/2]) = NaN;
    plot(tb,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',2.5); plot(tb,rc+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',2.5)
end
yline(n_half+0.5,'k'); xlim([0,diff(w)])
ylabel('PB glomeruli'); title(sprintf('(a) example lap: %s fly %s trial %s (gain %.2f) -- PB activity, bump (red), bar (blue)',example{:},g_fly(IG.fly_num(k_ex))),'FontWeight','normal')
nexttile(tl,5,[1,2]); hold on
k0 = find(ib,1); j0 = find(jf,1);
plot(xf(jf)-w(1),rad2deg(cu(jf)-cu(j0)),'Color',[0,0.45,0.74],'LineWidth',1.5); plot(xb(ib)-w(1),rad2deg(mu_u(ib)-mu_u(k0)),'Color',[0.85,0.2,0.1],'LineWidth',1.5)
for yv = -1080:360:360, yline(yv,':','Color',[.7,.7,.7],'HandleVisibility','off'), end
xlim([0,diff(w)]); xlabel('time (s)'); ylabel('rotation (deg)'); legend('cue (bar)','bump','Location','southwest')
title('bump laps the bar: ~-900\circ vs ~-90\circ in 20 s (dotted every 360\circ)','FontWeight','normal')
% (b) model
nexttile(tl,3); hold on
C = glme.Coefficients; keep = ~strcmp(C.Name,'(Intercept)') & ~startsWith(C.Name,'bpos');
nm = strrep(C.Name(keep),'_',' '); est = C.Estimate(keep); se = C.SE(keep); pv = C.pValue(keep); yy = numel(nm):-1:1;
errorbar(exp(est),yy,exp(est)-exp(est-1.96*se),exp(est+1.96*se)-exp(est),'horizontal','o','Color','k','MarkerFaceColor',col,'LineWidth',1.2)
for j = 1:numel(nm), text(exp(est(j)),yy(j)+0.32,sprintf('p = %.2g',pv(j)),'FontSize',8,'HorizontalAlignment','center'), end
xline(1,':k'); set(gca,'XScale','log'); yticks(fliplr(yy)); yticklabels(flipud(nm)); ylim([0.5,numel(nm)+0.6])
xlabel('odds ratio per SD (95% CI)'); title({'(b) what raises the odds of a lap','mixed logistic model, fly random effect (+ bump position)'},'FontWeight','normal','FontSize',9)
% (c) concentration
nexttile(tl,4); hold on
plot([0,1],[0,1],':k'); plot(xs,g_flyonly,'Color',[.6,.6,.6],'LineWidth',1.5); plot(xs,g_full,'Color',col,'LineWidth',2)
plot(0.2,g_full(round(0.2*n)),'o','Color',col,'MarkerFaceColor',col)
text(0.22,g_full(round(0.2*n))-0.05,sprintf('top 20%% of turns: %.0f%% of laps',100*g_full(round(0.2*n))),'FontSize',8)
xlim([0,1]); ylim([0,1]); xlabel('fraction of turns, ranked by predicted risk'); ylabel('fraction of laps captured')
legend('chance',sprintf('fly only (AUC %.2f)',auc_fly),sprintf('best model (AUC %.2f)',auc_full),'Location','southeast','FontSize',8)
title({'(c) laps are only partly predictable','cross-validated (10-fold, within fly)'},'FontWeight','normal','FontSize',9)
% (d) categories
nexttile(tl,7); hold on
r = zeros(4,2); nn = r;
for j = 1:4, m = cats{j,2}; r(j,:) = [mean(y(m & beh)), mean(y(m & ~beh))]; nn(j,:) = [sum(m & beh), sum(m & ~beh)]; end
hb = bar(1:4,100*r,0.8); hb(1).FaceColor = col; hb(2).FaceColor = [0.95,0.62,0.45]; hb(1).EdgeColor = 'none'; hb(2).EdgeColor = 'none';
yline(100*mean(y),':k','HandleVisibility','off'); xticks(1:4); xticklabels(cats(:,1)); xtickangle(20)
for j = 1:4, for c = 1:2, text(j+(c-1.5)*0.4,100*r(j,c)+0.3,sprintf('%d',nn(j,c)),'FontSize',7,'HorizontalAlignment','center'), end, end
ylabel('lap rate (%)'); legend('bar mostly behind','bar mostly in front','Location','northeast')
title({'(d) large, slow turns with the bar behind lap most','size / speed split at the median; numbers = turns'},'FontWeight','normal','FontSize',9)
% (e) bar occupancy
nexttile(tl,8); hold on
ctr = rad2deg(edges_b(1:end-1)+pi/n_bar_bins);
patch([ctr, fliplr(ctr)],[prctile(oc_null,2.5), fliplr(prctile(oc_null,97.5))],[.85,.85,.85],'EdgeColor','none')
plot(ctr,oc_obs,'o-','Color',col,'MarkerFaceColor',col,'LineWidth',1.5); yline(1,':k'); set(gca,'YScale','log')
xticks(ctr); xticklabels(arrayfun(@(v) sprintf('%.0f',v),ctr,'UniformOutput',false)); xlim([-180,180])
xlabel({'bar azimuth during the turn (deg)','0 = behind the fly, \pm180 = in front'}); ylabel('time with bar here: lap / other turns')
title({'(e) bar position during lap vs other turns',sprintf('grey = within-fly shuffle 95%%; p = %.2g',p_occ)},'FontWeight','normal','FontSize',9)
% (f) clustering
nexttile(tl,11,[1,2]); hold on
lbl = [arrayfun(@(j) sprintf('%g-%g s',ce(j),ce(j+1)),1:nc-2,'UniformOutput',false), {sprintf('> %g s',ce(nc-1)), 'none earlier'}];
patch([1:nc, nc:-1:1],100*[prctile(null_c,2.5), fliplr(prctile(null_c,97.5))],[.85,.85,.85],'EdgeColor','none')
plot(1:nc,100*rate_c,'o-','Color',col,'MarkerFaceColor',col,'LineWidth',1.5)
xticks(1:nc); xticklabels(lbl); xtickangle(30); xlim([0.5,nc+0.5]); ylabel('lap rate (%)'); xlabel('time since the previous lap (same trial)')
title({'(f) clustering in time (weak)',sprintf('grey = within-trial shuffle 95%%; <= 15 s: p = %.2g',p_clust)},'FontWeight','normal','FontSize',9)
% summary text
nexttile(tl,9,[1,2]); axis off
text(0,1,{'\bfSummary\rm', ...
    sprintf('%d turns >= %d%s (solid bump throughout) in %d high-gain LPsP > TH flies; %d laps (%.1f%%), lap = bump >= %g%s past the fly',n,min_turn_deg,char(176),numel(unique(D.fly)),sum(y),100*mean(y),overshoot_deg,char(176)), ...
    '\bullet No threshold: no turn speed or duration reliably produces a lap; even the riskiest turns lap ~1 in 5 times.', ...
    '\bullet Additive risk factors: larger but slower turns and a somewhat weaker (but still solid) bump are robust;', ...
    '   bar mostly behind the fly and bump position in the PB add a little. No interactions improved prediction.', ...
    '\bullet Weak or not robust: time in trial, where the bar is during the turn (e), clustering in time (f).', ...
    '\bullet Not predictive: absolute bump glomerulus, turn direction, bump amplitude/width before the turn.', ...
    '\bullet Laps also occur in the dark (fig\_rnai\_th\_lapping\_dark): the bar is not required.'}, ...
    'VerticalAlignment','top','FontSize',10)
title(tl,'Laps in high-gain LPsP > TH-RNAi flies: mostly stochastic, with additive risk factors','FontWeight','bold')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_lap_summary.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_lap_summary.pdf'),'ContentType','image','Resolution',300)

%% functions
function s = since_last(G)
    % time since the previous overshoot in the same trial (Inf if none)
    s = inf(height(G),1);
    for i = unique(G.trial)'
        idx = find(G.trial==i); last = -Inf;
        for k = idx', s(k) = G.t0(k) - last; if G.os(k), last = G.t0(k); end, end
    end
end
