%% fig_rnai_th_gain_stats
% Is the LPsP > TH-RNAi closed-loop integrative gain really bimodal, and do
% the high- / low-gain groups (fig_rnai_th_gain_bimodality.m) track
% experiment date? Per-fly gain = mean over the fly's closed-loop trials
% that pass the integrative-gain inclusion (fig_rnai_integrative_gain.m
% cache). Every test is also run on the controls.
%
% BIMODALITY, per genotype, on raw and log gain:
%   - 1 vs 2 component Gaussian mixture (unequal variances): dBIC (> 0
%     favors 2) and likelihood-ratio p from a parametric bootstrap under
%     the 1-Gaussian fit. A 2-component mixture also fits plain skew or a
%     lone outlier, so a small p alone is not evidence of two modes --
%     read it with:
%   - Ashman's D = sqrt(2)|mu1-mu2| / sqrt(s1^2+s2^2) (> 2 = cleanly
%     separated components) and the smaller component's weight.
%   - Silverman's critical-bandwidth test for > 1 mode (smoothed
%     bootstrap, uncalibrated, so conservative; sensitive to isolated
%     outliers, which make their own mode).
%   - ICC(1) of single-trial gains within fly (flies with >= 2 closed-loop
%     trials): how much of the spread is stable per fly vs trial noise.
%   - Brown-Forsythe test of spread, LPsP>TH vs each control.
%   No Hartigan dip test (no implementation in the repo).
%
% DATE, LPsP>TH (gain and high/low group), controls for gain only:
%   - Spearman trend of gain over date
%   - clustering by day: one-way F of gain on day, and (LPsP>TH) the
%     fraction of same-day fly pairs in the same group, both against
%     n_perm label permutations; ICC(day)
%   - the same over imaging sessions = runs of dates with no gap longer
%     than session_gaps(k) days (a cross is typically imaged over a few
%     days)
%   - first vs second half of the dataset (Fisher)
%   These are several overlapping tests: p ~ 0.04 in one grouping is weak.
%
% Figure 1: (a) per-fly gain histograms with the 2-component fit, TH set;
% (b) LPsP>TH gain vs date, colored by group, sessions shaded; (c) fraction
% high-gain per session. Figure 2: the same histograms for the VGlut set.
% Run from the repo root.

%% parameters + per-fly gain and date
gain_split   = 0.85;
n_boot       = 1000;   % parametric / smoothed bootstrap draws (bimodality)
n_perm       = 10000;  % label permutations (date tests)
session_gaps = [3,5,7]; % days; a longer gap between imaging dates starts a new session
geno_test    = {'lpsp>th','empty>th','lpsp>mcherry','lpsp>vglut','empty>vglut'};
geno_names   = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry','LPsP > VGlut-RNAi','Empty > VGlut-RNAi'};
geno_colors  = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13; 0.49,0.18,0.56; 0.47,0.67,0.19];
col_lo = [0.95,0.62,0.45]; col_hi = [0.60,0.15,0.00];
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','fly_list','inc','g_trial');
use = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
n_flies = max(IG.fly_num);
g_fly = nan(n_flies,1); fly_geno = cell(n_flies,1); fly_date = nan(n_flies,1); n_tr = zeros(n_flies,1);
for f = 1:n_flies
    fly_geno{f} = IG.genotype{find(IG.fly_num==f,1)};
    tr_i = find(use & IG.fly_num==f); n_tr(f) = numel(tr_i);
    if ~isempty(tr_i), g_fly(f) = mean(vertcat(IG.g_trial{tr_i}),'omitnan'); end
    tok = regexp(IG.fly_list{f},'\\(\d{8})\\','tokens','once');
    if ~isempty(tok), fly_date(f) = datenum(tok{1},'yyyymmdd'); end
end
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
gm_opts = statset('MaxIter',1000);
warning('off','stats:gmdistribution:FailedToConverge'); warning('off','stats:gmdistribution:FailedToConvergeReps');

%% bimodality tests
rng(1)
fprintf('%-13s %3s | %-3s %7s %7s %7s %6s | %7s %7s | %5s\n','genotype','n','', 'dBIC','LRT p','Ashman','w_min','h_crit','Silv p','ICC');
gm_fit = cell(numel(geno_test),1); bm_stats = nan(numel(geno_test),3); % raw scale: LRT p, Silverman p, ICC
for gi = 1:numel(geno_test)
    x_raw = g_fly(strcmp(fly_geno,geno_test{gi}) & ~isnan(g_fly));
    for sc = 1:2
        if sc == 1, x = x_raw; scl = 'raw'; else, x = log(x_raw); scl = 'log'; end
        n = numel(x);
        % 1 vs 2 component mixture, parametric bootstrap LRT under the 1-Gaussian fit
        [L1,L2,m2] = gmm_ll(x,gm_opts);
        dBIC = (-2*L1 + 2*log(n)) - (-2*L2 + 5*log(n));
        lrt = 2*(L2-L1); lrt_b = nan(n_boot,1);
        for b = 1:n_boot
            [l1,l2] = gmm_ll(mean(x) + std(x,1)*randn(n,1),gm_opts); lrt_b(b) = 2*(l2-l1);
        end
        p_lrt = (1+sum(lrt_b >= lrt))/(n_boot+1);
        s = sqrt(squeeze(m2.Sigma)); D = sqrt(2)*abs(diff(m2.mu))/sqrt(sum(s.^2));
        if sc == 1, gm_fit{gi} = m2; end
        % Silverman: smoothed bootstrap at the critical bandwidth, rescaled to the sample variance
        hc = crit_bw(x); sx = std(x); multi = false(n_boot,1);
        for b = 1:n_boot
            xs = x(randi(n,n,1));
            multi(b) = n_modes(mean(xs) + (xs - mean(xs) + hc*randn(n,1))/sqrt(1+hc^2/sx^2),hc) > 1;
        end
        icc_s = '';
        if sc == 1
            y = []; grp = [];
            for f = find(strcmp(fly_geno,geno_test{gi}) & n_tr >= 2)'
                tg = cellfun(@(z) mean(z,'omitnan'),IG.g_trial(use & IG.fly_num==f));
                y = [y; tg(:)]; grp = [grp; f*ones(numel(tg),1)]; %#ok<AGROW>
            end
            [~,~,grp] = unique(grp); icc_v = icc1(y,grp); icc_s = sprintf('%5.2f',icc_v);
            bm_stats(gi,:) = [p_lrt, mean(multi), icc_v];
        end
        fprintf('%-13s %3d | %-3s %7.1f %7.3f %7.2f %6.2f | %7.3f %7.3f | %s\n',geno_test{gi},n,scl,dBIC,p_lrt,D,min(m2.ComponentProportion),hc,mean(multi),icc_s);
    end
end
for c = {'empty>th','lpsp>mcherry'}
    a = g_fly(strcmp(fly_geno,'lpsp>th') & ~isnan(g_fly)); b = g_fly(strcmp(fly_geno,c{1}) & ~isnan(g_fly));
    p = vartestn([a;b],[ones(size(a));2*ones(size(b))],'TestType','BrownForsythe','Display','off');
    fprintf('Brown-Forsythe spread: lpsp>th (SD %.2f) vs %s (SD %.2f), p = %.2g\n',std(a),c{1},std(b),p);
end

%% date tests
rng(3)
for gname = {'lpsp>th','empty>th','lpsp>mcherry'}
    sel = strcmp(fly_geno,gname{1}) & ~isnan(g_fly) & ~isnan(fly_date);
    g = g_fly(sel); d = fly_date(sel); hi = g > gain_split; n = numel(g);
    fprintf('\n== %s: %d flies, %d dates (%s to %s)\n',gname{1},n,numel(unique(d)),datestr(min(d),'yyyy-mm-dd'),datestr(max(d),'yyyy-mm-dd'));
    [r,p] = corr(d,g,'type','Spearman'); fprintf('  trend: Spearman gain ~ date rho = %.2f, p = %.2f\n',r,p);
    [~,~,di] = unique(d);
    [pF,Fobs] = perm_p(@(y) fstat(y,di),g,n_perm);
    fprintf('  gain ~ day (one-way F, permuted): F = %.2f, p = %.3f; ICC(day) = %.2f\n',Fobs,pF,icc1(g,di));
    if ~strcmp(gname{1},'lpsp>th'), continue, end
    [pS,Sobs,Snull] = perm_p(@(y) samefrac(y,di),hi,n_perm);
    fprintf('  same-day pairs in same group: %.2f (chance %.2f), p = %.3f\n',Sobs,mean(Snull),pS);
    ud = unique(d);
    for gap = session_gaps
        si = cumsum([1; diff(ud) > gap]); sess = si(arrayfun(@(x) find(ud==x),d));
        pS = perm_p(@(y) samefrac(y,sess),hi,n_perm); pF = perm_p(@(y) fstat(y,sess),g,n_perm);
        fprintf('  sessions (gap > %d d): %d sessions; same-group pairs p = %.3f; gain ~ session F p = %.3f\n',gap,max(si),pS,pF);
    end
    early = d <= median(d);
    [~,p] = fishertest([sum(hi&early),sum(~hi&early); sum(hi&~early),sum(~hi&~early)]);
    fprintf('  early half %d/%d high, late half %d/%d high, Fisher p = %.2f\n',sum(hi&early),sum(early),sum(hi&~early),sum(~early),p);
end

%% figure
figure('Color','w','Position',[60,60,1500,820]); clf
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for gi = 1:3
    nexttile(tl,gi); hold on
    gain_hist(g_fly(strcmp(fly_geno,geno_test{gi}) & ~isnan(g_fly)),gm_fit{gi},geno_colors(gi,:),geno_names{gi})
    if gi == 1, xline(gain_split,'--k'); ylabel('flies'), end
end

sel = strcmp(fly_geno,'lpsp>th') & ~isnan(g_fly) & ~isnan(fly_date);
g = g_fly(sel); d = fly_date(sel); hi = g > gain_split;
ud = unique(d); si = cumsum([1; diff(ud) > session_gaps(1)]); sess = si(arrayfun(@(x) find(ud==x),d));
ax = nexttile(tl,[1,2]); hold on
for k = 1:max(si)
    dd = ud(si==k);
    patch([dd(1)-1,dd(end)+1,dd(end)+1,dd(1)-1],[0,0,2,2],[.93,.93,.93],'EdgeColor','none')
end
jit = (rand(size(d))-.5)*0.6;
scatter(d(~hi)+jit(~hi),g(~hi),30,col_lo,'filled'); scatter(d(hi)+jit(hi),g(hi),30,col_hi,'filled')
yline(gain_split,'--k'); yline(0.8,':k')
datetick('x','mm-dd','keeplimits'); ylim([0,2]); xlabel('experiment date'); ylabel('integrative gain, closed loop')
title(sprintf('(b) LPsP > TH-RNAi gain vs date; shaded = imaging sessions (no gap > %d d)',session_gaps(1)),'FontWeight','normal')

nexttile(tl); hold on
n_s = max(sess); fh = accumarray(sess,hi,[n_s,1],@mean); ns = accumarray(sess,1,[n_s,1]);
bar(1:n_s,fh,'FaceColor',col_hi,'FaceAlpha',.7,'EdgeColor','none'); yline(mean(hi),'--k')
lbl = arrayfun(@(k) sprintf('%s (n=%d)',datestr(min(ud(si==k)),'mm-dd'),ns(k)),1:n_s,'UniformOutput',false);
xticks(1:n_s); xticklabels(lbl); xtickangle(40); ylim([0,1]); ylabel('fraction high-gain flies')
title('(c) fraction high-gain per session (dashed = overall)','FontWeight','normal')
title(tl,{'LPsP > TH-RNAi integrative gain: bimodality and experiment date (test results in the command window)', ...
    '(a) histograms with 2-component Gaussian mixture fits (dashed = components); dotted = 0.8 cue gain, dashed = group split'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_stats.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_stats.pdf'),'ContentType','vector')

%% figure: VGlut set gain histograms
% Same panels as (a) above for lpsp>vglut, empty>vglut, lpsp>mcherry, with
% the raw-scale test results from the bimodality section in each title.
vg = cellfun(@(g) find(strcmp(geno_test,g)),{'lpsp>vglut','empty>vglut','lpsp>mcherry'});
figure('Color','w','Position',[60,60,1500,460]); clf
tl = tiledlayout(1,numel(vg),'TileSpacing','compact','Padding','compact');
for k = 1:numel(vg)
    gi = vg(k);
    nexttile(tl,k); hold on
    gain_hist(g_fly(strcmp(fly_geno,geno_test{gi}) & ~isnan(g_fly)),gm_fit{gi},geno_colors(gi,:),geno_names{gi})
    if k == 1, ylabel('flies'), end
    t = get(gca,'Title'); t.String = [cellstr(splitlines(string(t.String))); {sprintf('LRT p = %.3f, Silverman p = %.3f, ICC = %.2f',bm_stats(gi,:))}];
end
title(tl,{'VGlut set: per-fly closed-loop integrative gain with 2-component Gaussian mixture fits', ...
    'dashed = components, solid = mixture; dotted = 0.8 cue gain; tests on raw gain (see command window for log scale)'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_vglut_gain_hist.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_vglut_gain_hist.pdf'),'ContentType','vector')

%% functions
function gain_hist(x,m2,c,name)
    % per-fly gain histogram (0.1 bins) with the 2-component mixture fit scaled to counts
    edges = 0:0.1:2.3; xg = linspace(0,2.3,400)';
    histogram(x,edges,'FaceColor',c,'EdgeColor','w','FaceAlpha',.6)
    s = sqrt(squeeze(m2.Sigma)); w = m2.ComponentProportion; sc = numel(x)*0.1;
    for k = 1:2, plot(xg,sc*w(k)*normpdf(xg,m2.mu(k),s(k)),'--','Color',c*.6,'LineWidth',1), end
    plot(xg,sc*pdf(m2,xg),'-','Color',c*.6,'LineWidth',1.5)
    xline(0.8,':k'); xlim([0,2.3]); xlabel('integrative gain, closed loop')
    D = sqrt(2)*abs(diff(m2.mu))/sqrt(sum(s.^2));
    title(sprintf('%s (n=%d)\n2-comp fit: \\mu = %.2f / %.2f, w = %.2f / %.2f, D = %.1f',name,numel(x),m2.mu,w,D),'FontWeight','normal')
end
function [L1,L2,m2] = gmm_ll(x,opts)
    % log-likelihoods of the 1-Gaussian and best-of-5 2-component fits
    L1 = sum(log(normpdf(x,mean(x),std(x,1))));
    best = -inf; m2 = [];
    for r = 1:5 % own replicates so one failed fit doesn't abort
        try
            m = fitgmdist(x,2,'Options',opts,'RegularizationValue',1e-4,'Start','plus');
            if -m.NegativeLogLikelihood > best, best = -m.NegativeLogLikelihood; m2 = m; end
        catch
        end
    end
    L2 = max(best,L1);
end

function k = n_modes(x,h)
    % number of local maxima of a Gaussian KDE with bandwidth h
    g = linspace(min(x)-3*h,max(x)+3*h,512);
    d = sum(exp(-0.5*((g(:)-x(:)')/h).^2),2);
    k = sum(diff(sign(diff(d))) < 0);
end

function hc = crit_bw(x)
    % Silverman's critical bandwidth: smallest h giving a unimodal KDE (bisection)
    lo = 1e-3*std(x); hi = 2*std(x);
    for it = 1:40
        mid = (lo+hi)/2;
        if n_modes(x,mid) > 1, lo = mid; else, hi = mid; end
    end
    hc = hi;
end

function [p,obs,null] = perm_p(stat,y,n_perm)
    % one-sided permutation p of stat(y) against shuffled y
    obs = stat(y); null = zeros(n_perm,1);
    for k = 1:n_perm, null(k) = stat(y(randperm(numel(y)))); end
    p = mean(null >= obs);
end

function F = fstat(y,grp)
    % one-way ANOVA F of y on group index grp
    m = accumarray(grp,y,[],@mean); c = accumarray(grp,1); ok = c > 0;
    ssb = sum(c(ok).*(m(ok)-mean(y)).^2); ssw = sum((y - m(grp)).^2);
    F = (ssb/(sum(ok)-1)) / (ssw/(numel(y)-sum(ok)));
end

function s = samefrac(lab,grp)
    % fraction of within-group pairs that share the binary label
    same = 0; tot = 0;
    for k = unique(grp)'
        l = lab(grp==k); m = numel(l); if m < 2, continue, end
        h = sum(l); same = same + h*(h-1)/2 + (m-h)*(m-h-1)/2; tot = tot + m*(m-1)/2;
    end
    s = same/tot;
end

function r = icc1(y,grp)
    % ICC(1) for unbalanced groups
    [~,tbl] = anova1(y,grp,'off'); MSb = tbl{2,4}; MSw = tbl{3,4};
    c = accumarray(grp,1); n = numel(y); k0 = (n - sum(c.^2)/n)/(numel(c)-1);
    r = (MSb-MSw)/(MSb+(k0-1)*MSw);
end
