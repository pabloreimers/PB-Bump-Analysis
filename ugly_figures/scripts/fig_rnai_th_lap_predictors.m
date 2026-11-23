%% fig_rnai_th_lap_predictors
% What predicts a lap in high-gain LPsP > TH-RNAi flies? Turn by turn, in
% the high-gain flies that lap (fig_rnai_th_lapping_examples.m), closed
% loop, included trials, turns >= min_turn_deg from the turn cache
% (lpsp_rnai_turns_claude.m).
%
%   outcome: OVERSHOOT -- the bump rotates >= overshoot_deg further than the
%     fly, in the turn's direction, from turn onset to post_s after its end
%     (bump lag). Absolute excess rather than the bump/fly ratio: the ratio
%     favours small turns for purely numerical reasons. Only turns with a
%     solid bump (mean PVA strength over the turn >= rho_min, and 10th
%     percentile over the whole scored interval >= rho_p10_min), since
%     weak-bump noise -- often a bump fading right after the turn -- inflates
%     apparent overshoots.
%   predictors: turn speed (mean, deg/s), duration, forward speed in the 2 s
%     before the turn, time in trial, PVA strength during the turn, and the
%     bump's position in the PB at turn onset -- absolute (glomerulus within
%     a half) and relative to the turn direction (counter-clockwise turns
%     mirrored, so "ahead" is the same side).
%   tests: (a) overshoot rate vs speed / duration / forward speed (quintiles,
%     binomial 95% CI); (b) clustering -- rate vs time since the last
%     overshoot in the same trial, against a within-trial label shuffle
%     (keeps each trial's rate); (c) mixed-effects logistic regression,
%     fly random intercept, odds ratio per SD; (d) bump position -- rate per
%     position bin, chi-square across bins against a within-fly label
%     shuffle (keeps each fly's rate and position distribution; also within
%     fly x turn direction), and whether
%     flies' overshoot-onset positions share a preferred direction (Rayleigh
%     on per-fly circular means); (e)-(f) where the bar (cue) went during the
%     turn (ft.cue, 0 = behind the fly): overshoot rate among turns whose bar
%     swept through each azimuth bin, and lap vs non-lap time spent per bin;
%     plus whether the bar passing directly behind the fly matters beyond turn
%     size (likelihood-ratio test in the model).
% Needs all_data (loaded if absent), the integrative-gain cache and the turn
% cache. Run from the repo root.

%% parameters
min_turn_deg = 30; overshoot_deg = 90; post_s = 1.0; rho_min = 0.4; gain_split = 0.85;
rho_p10_min = 0.3; % 10th percentile of PVA strength over the whole scored interval (onset to post_s after the turn): drops overshoots made while the bump fades
n_pos_bins = 8; n_perm = 2000; lap_fly_min_overshoots = 1; % flies need >= this many overshoots to count as lapping
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
rng(1)

%% per-turn table, high-gain LPsP>TH flies
hi_fl = find(strcmp(fly_geno,'lpsp>th') & g_fly > gain_split);
tr_use = find(~IG.is_dark & IG.inc & ismember(IG.fly_num,hi_fl));
tt = [TC.turns.trial]';
rows = zeros(0,16); % trial fly t0 dur dfly dbump fwd_pre rho_turn bump_glom bump_rel t_trial mean_speed cue_az cue_rel crossed_behind rho_p10
bar_traj = {};       % bar azimuth (ft.cue, rad; 0 = behind the fly) at each fictrac sample during the turn
for i = tr_use(:)'
    tr = all_data(i); xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    mu_u = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'));
    fs = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',max(1,round(0.5/median(diff(xf)))));
    n_half = size(tr.im.f,1)/2; a1 = tr.im.alpha(1); dg = 2*pi/n_half;
    for k = find(tt==i)'
        T = TC.turns(k); t0 = T.onset_s; t1 = t0 + T.dur_s;
        b0 = find(xb >= t0,1); b1 = find(xb <= t1 + post_s,1,'last'); bt = find(xb <= t1,1,'last');
        f0 = find(xf >= t0,1); f1 = find(xf <= t1,1,'last');
        if isempty(b0) || isempty(b1) || isempty(bt) || bt < b0 || isempty(f1) || f1 <= f0, continue, end
        dfly = rad2deg(hd(f1) - hd(f0)); if abs(dfly) < min_turn_deg, continue, end
        glom = mod(mu(b0) - a1 + dg/2, 2*pi)/dg + 0.5;            % fractional glomerulus within a half (1..n_half)
        rel = mod(sign(dfly)*(mu(b0) - a1), 2*pi);                % position, mirrored for CCW turns
        cue_on = angle(exp(1i*circ_mean(tr.ft.cue(xf >= t0-0.25 & xf <= t0 & ~isnan(tr.ft.cue(:)))))); % bar azimuth at onset (ft.cue, rad)
        if isnan(cue_on), cue_on = 0; end
        c_in = tr.ft.cue(f0:f1); c_in = c_in(~isnan(c_in)); c_in = c_in(:);
        cu_in = unwrap(c_in); crossed = ~isempty(cu_in) && floor(min(cu_in)/(2*pi)) ~= floor(max(cu_in)/(2*pi)) ... % unwrapped path spans a multiple of 2*pi = passes through 0 (behind)
            || any(abs(angle(exp(1i*c_in))) < deg2rad(5));
        rows(end+1,:) = [i, IG.fly_num(i), t0, T.dur_s, dfly, rad2deg(mu_u(b1) - mu_u(b0)), ...
            mean(fs(xf >= t0-2 & xf < t0)), mean(rho(b0:bt)), glom, rel, t0, abs(dfly)/T.dur_s, ...
            cue_on, mod(sign(dfly)*cue_on + pi, 2*pi) - pi, crossed, prctile(rho(b0:b1),10)]; %#ok<AGROW>
        bar_traj{end+1,1} = angle(exp(1i*c_in)); %#ok<AGROW>
    end
end
D = array2table(rows,'VariableNames',{'trial','fly','t0','dur','dfly','dbump','fwd_pre','rho_turn','glom','rel','t_trial','mean_speed','cue_az','cue_rel','crossed','rho_p10'});
D.bar_traj = bar_traj;
D.exc = sign(D.dfly).*D.dbump - abs(D.dfly);
D.os = D.exc >= overshoot_deg;
D = D(D.rho_turn >= rho_min & D.rho_p10 >= rho_p10_min,:);
n_os = accumarray(D.fly,D.os,[max(IG.fly_num),1]);
D = D(n_os(D.fly) >= lap_fly_min_overshoots,:); % flies that overshoot at all
fprintf('%d turns (>= %d deg, solid bump) in %d high-gain LPsP>TH flies; %d overshoots (%.1f%%)\n',height(D),min_turn_deg,numel(unique(D.fly)),sum(D.os),100*mean(D.os));

%% (c) mixed-effects logistic regression
Zs = @(v) (v - mean(v,'omitnan'))./std(v,'omitnan');
M = table(D.os,categorical(D.fly),Zs(D.mean_speed),Zs(D.dur),Zs(D.fwd_pre),Zs(D.t_trial),Zs(D.rho_turn),cos(2*pi*(D.glom-0.5)/max(D.glom)),sin(2*pi*(D.glom-0.5)/max(D.glom)), ...
    'VariableNames',{'os','fly','turn_speed','duration','fwd_before','time_in_trial','bump_strength','pos_cos','pos_sin'});
M.log_size = Zs(log(abs(D.dfly))); M.crossed = double(D.crossed);
glme0 = fitglme(M,'os ~ turn_speed + duration + fwd_before + time_in_trial + bump_strength + (1|fly)','Distribution','Binomial','FitMethod','Laplace');
glme1 = fitglme(M,'os ~ turn_speed + duration + fwd_before + time_in_trial + bump_strength + pos_cos + pos_sin + (1|fly)','Distribution','Binomial','FitMethod','Laplace');
cmp_pos = compare(glme0,glme1);
% does the bar passing behind the fly during the turn matter, beyond turn size? (bigger turns sweep the bar further)
glme2 = fitglme(M,'os ~ turn_speed + duration + fwd_before + time_in_trial + bump_strength + log_size + (1|fly)','Distribution','Binomial','FitMethod','Laplace');
glme3 = fitglme(M,'os ~ turn_speed + duration + fwd_before + time_in_trial + bump_strength + log_size + crossed + (1|fly)','Distribution','Binomial','FitMethod','Laplace');
cmp_cue = compare(glme2,glme3); or_cross = exp(glme3.Coefficients.Estimate(end));
fprintf('bar passed directly behind during the turn: %.1f%% overshoot (n=%d) vs %.1f%% (n=%d); beyond turn size: odds ratio %.2f, LR p = %.3g\n', ...
    100*mean(D.os(D.crossed==1)),sum(D.crossed==1),100*mean(D.os(D.crossed==0)),sum(D.crossed==0),or_cross,cmp_cue.pValue(2));
disp(glme0.Coefficients(:,{'Name','Estimate','SE','pValue'}))
fprintf('adding bump position (cos + sin of glomerulus) to the model: likelihood-ratio p = %.3g\n',cmp_pos.pValue(2));
pr = predict(glme0,M,'Conditional',false); [~,~,~,auc0] = perfcurve(M.os,pr,true);
fprintf('fixed-effects model AUC = %.2f\n',auc0);

%% (d) bump position: rate per bin, chi-square vs within-fly shuffle; per-fly preferred onset position
edges_g = linspace(0.5,max(D.glom)+1e-9,n_pos_bins+1); edges_r = linspace(0,2*pi,n_pos_bins+1); edges_c = linspace(-pi,pi+1e-9,n_pos_bins+1);
pos = struct('name',{'bump position, absolute','bump position, along travel'},'x',{D.glom,D.rel},'edges',{edges_g,edges_r});
for k = 1:numel(pos)
    b = discretize(pos(k).x,pos(k).edges); pos(k).b = b;
    pos(k).chi = chi_bins(D.os,b,n_pos_bins); null = zeros(n_perm,1);
    for p = 1:n_perm
        lab = D.os; for f = unique(D.fly)', idx = find(D.fly==f); lab(idx) = lab(idx(randperm(numel(idx)))); end
        null(p) = chi_bins(lab,b,n_pos_bins);
    end
    pos(k).p = mean(null >= pos(k).chi);
    % stricter for the direction-relative position: shuffle within fly x turn direction (a fly lapping
    % mostly one way, with its bump at habitual positions, can't fake a position effect)
    grp = findgroups(D.fly,sign(D.dfly)); null2 = zeros(n_perm,1);
    for p = 1:n_perm
        lab = D.os; for g = unique(grp)', idx = find(grp==g); lab(idx) = lab(idx(randperm(numel(idx)))); end
        null2(p) = chi_bins(lab,b,n_pos_bins);
    end
    pos(k).p2 = mean(null2 >= pos(k).chi);
    pos(k).rate = accumarray(b,D.os,[n_pos_bins,1],@mean); pos(k).n = accumarray(b,1,[n_pos_bins,1]); pos(k).nos = accumarray(b,D.os,[n_pos_bins,1]);
    ang = 2*pi*(pos(k).x - pos(k).edges(1))/(pos(k).edges(end)-pos(k).edges(1)); % map to circle
    mf = []; for f = unique(D.fly(D.os))', a = ang(D.os & D.fly==f); if numel(a) >= 3, mf(end+1) = circ_mean(a); end, end %#ok<AGROW>
    pos(k).p_rayleigh = circ_rtest(mf(:)); pos(k).n_fly = numel(mf);
    fprintf('position, %s: rate per bin %s%%, within-fly shuffle p = %.3g (within fly x direction p = %.3g); per-fly preferred onset (%d flies >= 3 overshoots) Rayleigh p = %.3g\n', ...
        pos(k).name,mat2str(round(100*pos(k).rate',1)),pos(k).p,pos(k).p2,pos(k).n_fly,pos(k).p_rayleigh);
end

%% (e) where the bar went during the turn (0 = behind the fly, +-180 = in front)
% pass-through: for each azimuth bin, the overshoot rate among turns whose bar swept through that bin;
% occupancy: time lap turns spent with the bar in each bin, relative to non-lap turns. Significance from
% within-fly label shuffles (keeps each fly's rate and its turns' bar paths).
edges_b = linspace(-pi,pi,n_pos_bins+1); nbb = n_pos_bins;
visit = false(height(D),nbb); occ = zeros(height(D),nbb);
for k = 1:height(D)
    bb = discretize(D.bar_traj{k},edges_b); bb = bb(~isnan(bb));
    if isempty(bb), continue, end
    visit(k,unique(bb)) = true; occ(k,:) = accumarray(bb(:),1,[nbb,1])'/numel(bb);
end
pass_rate = @(y) arrayfun(@(j) mean(y(visit(:,j))),1:nbb);
occ_ratio = @(y) mean(occ(y,:),1)./mean(occ(~y,:),1);
pr_obs = pass_rate(D.os); oc_obs = occ_ratio(D.os);
pass_null = zeros(n_perm,nbb); occ_null = zeros(n_perm,nbb);
for p = 1:n_perm
    lab = D.os; for f = unique(D.fly)', idx = find(D.fly==f); lab(idx) = lab(idx(randperm(numel(idx)))); end
    pass_null(p,:) = pass_rate(lab); occ_null(p,:) = occ_ratio(lab);
end
p_pass = mean(var(pass_null,0,2) >= var(pr_obs));             % overall: do the bins differ more than chance?
p_pass_bin = mean(abs(pass_null - mean(pass_null)) >= abs(pr_obs - mean(pass_null)));
p_occ = mean(sum((log(occ_null)).^2,2) >= sum(log(oc_obs).^2));
fprintf('bar pass-through rate per bin (centers %s deg): %s%%; overall within-fly shuffle p = %.3g\n',mat2str(round(rad2deg(edges_b(1:end-1)+pi/nbb))),mat2str(round(100*pr_obs,1)),p_pass);
fprintf('bar occupancy ratio lap/non-lap: %s; overall p = %.3g\n',mat2str(round(oc_obs,2)),p_occ);
% per fly: lap turns' share of time with the bar behind (|az| < 45) vs non-lap turns
bh = abs(edges_b(1:end-1)+pi/nbb) < pi/4 + 1e-9; fr = [];
for f = unique(D.fly(D.os))', m = D.fly==f; if sum(D.os(m)) >= 2, fr(end+1,:) = [mean(sum(occ(m & D.os,bh),2)), mean(sum(occ(m & ~D.os,bh),2))]; end, end %#ok<AGROW>
p_fly_behind = signrank(fr(:,1),fr(:,2));
fprintf('per fly (>= 2 overshoots, n=%d): time with bar behind (|az| < 45) higher in lap turns for %d flies; signrank p = %.3g\n',size(fr,1),sum(fr(:,1) > fr(:,2)),p_fly_behind);

%% (b) clustering: rate vs time since the previous overshoot in the same trial
G = sortrows(D,{'trial','t0'});
since = since_last(G);
ce = [0,5,15,30,60,120,Inf]; cb = discretize(since,ce); cb(isinf(since)) = numel(ce); % last bin = no earlier overshoot
nc = numel(ce); rate_c = accumarray(cb,G.os,[nc,1],@mean); n_c = accumarray(cb,1,[nc,1]);
null_c = zeros(n_perm,nc);
for p = 1:n_perm
    H = G; for i = unique(H.trial)', idx = find(H.trial==i); H.os(idx) = H.os(idx(randperm(numel(idx)))); end
    s = since_last(H); c = discretize(s,ce); c(isinf(s)) = nc; null_c(p,:) = accumarray(c,H.os,[nc,1],@mean)';
end
hp = since <= 15;
obs = mean(G.os(hp)) - mean(G.os(~hp));
nul = zeros(n_perm,1);
for p = 1:n_perm
    H = G; for i = unique(H.trial)', idx = find(H.trial==i); H.os(idx) = H.os(idx(randperm(numel(idx)))); end
    s = since_last(H); nul(p) = mean(H.os(s <= 15)) - mean(H.os(s > 15));
end
p_clust = mean(nul >= obs);
fprintf('clustering: P(overshoot | one in previous 15 s) %.1f%% vs %.1f%%; within-trial shuffle p = %.3g\n',100*mean(G.os(hp)),100*mean(G.os(~hp)),p_clust);

%% figure
figure('Color','w','Position',[40,40,2100,960]); clf
tl = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
qv = {'mean_speed','turn speed (deg/s)','turn speed'; 'dur','turn duration (s)','turn duration'; 'fwd_pre','forward speed, 2 s before','forward speed before the turn'};
for k = 1:3
    nexttile; hold on
    x = D.(qv{k,1}); q = prctile(x,0:20:100); b = discretize(x,q,'IncludedEdge','right'); b(isnan(b)) = 1;
    [r,lo,hi] = binom_ci(D.os,b,5);
    errorbar(1:5,100*r,100*(r-lo),100*(hi-r),'o-','Color',[0.60,0.15,0.00],'MarkerFaceColor',[0.60,0.15,0.00],'LineWidth',1.5)
    yline(100*mean(D.os),':k'); xlim([0.5,5.5]); ylim([0,max(100*hi)*1.15])
    xticks(1:5); xticklabels(arrayfun(@(j) sprintf('%.3g-%.3g',q(j),q(j+1)),1:5,'UniformOutput',false)); xtickangle(30)
    xlabel(qv{k,2}); if k == 1, ylabel(sprintf('turns with overshoot >= %g%s (%%)',overshoot_deg,char(176))), end
    title(sprintf('(a) by %s (quintiles)',qv{k,3}),'FontWeight','normal')
end
nexttile; hold on
lbl = [arrayfun(@(j) sprintf('%g-%g s',ce(j),ce(j+1)),1:nc-2,'UniformOutput',false), {sprintf('> %g s',ce(nc-1)), 'none earlier'}];
pr_null = mean(null_c,1,'omitnan'); lo_null = prctile(null_c,2.5); hi_null = prctile(null_c,97.5);
patch([1:nc, nc:-1:1],100*[lo_null, fliplr(hi_null)],[.85,.85,.85],'EdgeColor','none')
plot(1:nc,100*pr_null,'-','Color',[.5,.5,.5])
plot(1:nc,100*rate_c,'o-','Color',[0.60,0.15,0.00],'MarkerFaceColor',[0.60,0.15,0.00],'LineWidth',1.5)
xticks(1:nc); xticklabels(lbl); xtickangle(30); xlim([0.5,nc+0.5])
xlabel('time since the previous overshoot (same trial)'); ylabel('overshoot rate (%)')
title({'(b) clustering: grey = within-trial shuffle (95%)',sprintf('<= 15 s: p = %.3g',p_clust)},'FontWeight','normal')
nexttile(6); hold on
C = glme0.Coefficients; nm = C.Name(2:end); est = C.Estimate(2:end); se = C.SE(2:end); pv = C.pValue(2:end);
yy = numel(nm):-1:1;
errorbar(exp(est),yy,exp(est)-exp(est-1.96*se),exp(est+1.96*se)-exp(est),'horizontal','o','Color','k','MarkerFaceColor',[0.60,0.15,0.00],'LineWidth',1.2)
xline(1,':k'); set(gca,'XScale','log'); yticks(fliplr(yy)); yticklabels(flipud(strrep(nm,'_',' '))); ylim([0.5,numel(nm)+0.5])
for j = 1:numel(nm), text(exp(est(j)),yy(j)+0.3,sprintf('p = %.2g',pv(j)),'FontSize',8,'HorizontalAlignment','center'), end
xlabel('odds ratio per SD (95% CI)'); title({sprintf('(c) mixed logistic model (AUC %.2f)',auc0), ...
    sprintf('adding bump position: LR p = %.2g',cmp_pos.pValue(2)), sprintf('bar passed behind (beyond size): OR %.2f, p = %.2g',or_cross,cmp_cue.pValue(2))},'FontWeight','normal','FontSize',9)
nexttile(5); axis off
text(0,1,{sprintf('%d turns >= %d%s with a solid bump (mean \\rho >= %g)',height(D),min_turn_deg,char(176),rho_min), ...
    sprintf('%d high-gain LPsP > TH flies that overshoot',numel(unique(D.fly))), ...
    sprintf('overshoot = bump >= %g%s past the fly,',overshoot_deg,char(176)), sprintf('  onset to %g s after the turn',post_s), ...
    sprintf('%d overshoots (%.1f%% of turns)',sum(D.os),100*mean(D.os)),'', ...
    '(d)-(f): shuffles keep each fly''s overshoot rate;', ...
    'Rayleigh: do flies'' overshoot onsets share a', '  preferred position (per-fly circular means)', '', ...
    'bar azimuth: ft.cue, 0 = behind the fly, \pm180 = in front'},'VerticalAlignment','top','FontSize',9)
for k = 1:numel(pos)
    nexttile(6+k); hold on
    [r,lo,hi] = binom_ci(D.os,pos(k).b,n_pos_bins);
    bar(1:n_pos_bins,100*r,0.7,'FaceColor',[0.60,0.15,0.00],'FaceAlpha',.6,'EdgeColor','none')
    errorbar(1:n_pos_bins,100*r,100*(r-lo),100*(hi-r),'k.','LineWidth',1)
    yline(100*mean(D.os),':k'); xlim([0.5,n_pos_bins+0.5])
    if k == 1
        xticks(1:n_pos_bins); xticklabels(arrayfun(@(j) sprintf('%.0f-%.0f',pos(k).edges(j),pos(k).edges(j+1)),1:n_pos_bins,'UniformOutput',false))
        xlabel('bump glomerulus at onset (each PB half)')
    else
        xticks(1:n_pos_bins); xticklabels(arrayfun(@(j) sprintf('%.0f-%.0f',rad2deg(pos(k).edges(j)),rad2deg(pos(k).edges(j+1))),1:n_pos_bins,'UniformOutput',false))
        if k == 2, xlabel({'bump position along its travel (deg)','0 = PB-half end, 180 = mid-half'})
        end
    end
    xtickangle(30); ylabel('overshoot rate (%)')
    title({sprintf('(d) %s',pos(k).name),sprintf('shuffle p = %.2g / %.2g (fly / fly x dir)',pos(k).p,pos(k).p2), ...
        sprintf('shared onset position: Rayleigh p = %.2g',pos(k).p_rayleigh)},'FontWeight','normal','FontSize',8)
end

ctr = rad2deg(edges_b(1:end-1)+pi/nbb); blbl = arrayfun(@(j) sprintf('%.0f',ctr(j)),1:nbb,'UniformOutput',false);
nexttile(9); hold on
lo_n = prctile(pass_null,2.5); hi_n = prctile(pass_null,97.5);
patch([1:nbb, nbb:-1:1],100*[lo_n, fliplr(hi_n)],[.85,.85,.85],'EdgeColor','none')
bar(1:nbb,100*pr_obs,0.6,'FaceColor',[0.60,0.15,0.00],'FaceAlpha',.7,'EdgeColor','none')
for j = find(p_pass_bin < 0.05), text(j,100*pr_obs(j)+0.4,'*','HorizontalAlignment','center','FontSize',14), end
yline(100*mean(D.os),':k'); xlim([0.5,nbb+0.5]); xticks(1:nbb); xticklabels(blbl)
xlabel({'bar azimuth the turn swept through (deg)','0 = behind the fly, \pm180 = in front'}); ylabel('overshoot rate among turns passing the bin (%)')
title({'(e) where the bar passed during the turn',sprintf('grey = within-fly shuffle 95%%; overall p = %.2g; * bin p < 0.05',p_pass)},'FontWeight','normal','FontSize',8)
nexttile(10); hold on
lo_o = prctile(occ_null,2.5); hi_o = prctile(occ_null,97.5);
patch([1:nbb, nbb:-1:1],[lo_o, fliplr(hi_o)],[.85,.85,.85],'EdgeColor','none')
plot(1:nbb,oc_obs,'o-','Color',[0.60,0.15,0.00],'MarkerFaceColor',[0.60,0.15,0.00],'LineWidth',1.5)
yline(1,':k'); set(gca,'YScale','log'); xlim([0.5,nbb+0.5]); xticks(1:nbb); xticklabels(blbl)
xlabel({'bar azimuth during the turn (deg)','0 = behind the fly, \pm180 = in front'}); ylabel('time with bar here: lap / non-lap turns')
title({'(f) bar occupancy during lap vs other turns',sprintf('overall p = %.2g; per fly, more time behind in laps: p = %.2g (%d flies)',p_occ,p_fly_behind,size(fr,1))},'FontWeight','normal','FontSize',8)
title(tl,'What predicts a lap? High-gain LPsP > TH-RNAi flies, turn by turn (closed loop)')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_lap_predictors.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_lap_predictors.pdf'),'ContentType','vector')

%% functions
function c = chi_bins(y,b,nb)
    % chi-square of event counts across bins against a uniform rate
    ok = ~isnan(b); n = accumarray(b(ok),1,[nb,1]); o = accumarray(b(ok),y(ok),[nb,1]); e = n*mean(y(ok));
    c = sum((o - e).^2./max(e,eps)) + sum(((n-o) - (n-e)).^2./max(n-e,eps));
end

function s = since_last(G)
    % time since the previous overshoot in the same trial (Inf if none)
    s = inf(height(G),1);
    for i = unique(G.trial)'
        idx = find(G.trial==i); last = -Inf;
        for k = idx', s(k) = G.t0(k) - last; if G.os(k), last = G.t0(k); end, end
    end
end

function [r,lo,hi] = binom_ci(y,b,nb)
    % rate and Wilson 95% CI per bin
    ok = ~isnan(b); n = accumarray(b(ok),1,[nb,1]); k = accumarray(b(ok),y(ok),[nb,1]);
    r = k./n; z = 1.96; den = 1 + z^2./n; cen = (r + z^2./(2*n))./den; hw = z*sqrt(r.*(1-r)./n + z^2./(4*n.^2))./den;
    lo = cen - hw; hi = cen + hw;
end
