%% fig09_figure3_p2x2_lme
% Linear mixed models on the fig08 "ejection trials" data: does the
% LPsP-stimulation effect on EPG build up over successive stims, at
% matched ejected dye, and does it saturate?
%
% Data (same pipeline as fig_figure3_p2x2_traces.m with ej_trials = true):
% lpsp_p2x2_walking_20260728_post0708.mat, 300 s trials only (29 x 0.5 s
% pulses every 10 s), trials with >= 1 confirmed ejection, every commanded
% stim in them. Round 1 / 2 = first / second such trial of a fly's
% same-side run (count restarts when the ejection side switches); later
% trials are left out.
%
% One row per stim:
%   resp - ejection-side EPG GCaMP dF/F, mean over 0-10 s after the stim
%          (baseline -2..0 s), i.e. the fig08 summary measure per stim
%   dye  - ejection-side summed dye dF/F over 0-3 s (how much was
%          ejected; kept out of the response window), centred on its mean
%   lpsp - 0 = empty>P2X2, 1 = lpsp>P2X2;   fly - grouping factor
%   x    - stim count, 0-based, in one of two arrangements:
%          POOLED     = stim # within the trial (0-28), rounds 1 and 2
%                       overlaid (assumes a full reset between trials)
%          CONTINUOUS = count runs on through round 2 (0-57; assumes no
%                       reset across the ~5 min gap)
% Shapes of the build-up, each with a random intercept per fly (per-fly
% random slopes fitted to ~0 in an earlier version):
%   linear      resp ~ dye + x*lpsp + (1|fly)
%   log         resp ~ dye + log2(x+1)*lpsp + (1|fly)
%   saturating  resp ~ dye + sat*lpsp + (1|fly),  sat = 1 - exp(-x/tau),
%               tau (in stims, shared by both genotypes) profiled over a
%               grid by maximum likelihood (counted as one extra parameter
%               in AIC)
% In the saturating model: 'lpsp' = lpsp - empty at the first stim,
% 'sat' = how much empty rises to its plateau, 'sat:lpsp' = how much MORE
% lpsp rises; lpsp - empty at the plateau = lpsp + sat:lpsp (CI from the
% fixed-effect covariance).
%
% Figure: rows = pooled / continuous. A per-fly binned data (5-stim bins)
% with the saturating fit per genotype; B lpsp - empty, data vs model with
% 95% band; C log-likelihood profile over tau; D AIC of all six models.
%
% Run from the repo root: addpath('ugly_figures/ionto_celltypes'); fig09_figure3_p2x2_lme

data_dir = '.data';
out_dir  = fullfile('ugly_figures','ionto_celltypes');
cl = [0.90 0.45 0.05]; ce = [0.10 0.40 0.80];
bin_size = 5; max_stims = 30;
tau_grid = [0.5 1:1:20 22:2:40 45:5:80];

%% per-stim data (fig08 ej_trials pipeline)
R.t = -2:0.1:10; R.base = [-2 0]; R.side_win = 2;
tmp = load(fullfile(data_dir,'lpsp_p2x2_walking_20260728_post0708.mat'),'all_data'); ad = tmp.all_data(:); clear tmp
R = collect(R,ad,max_stims);
R = keep_rows(R,trial_has_ej(R));
R.tord = trial_round(R);
keep = R.tord <= 2;
resp = mean(R.g(:,R.t >= 0 & R.t < 10),2,'omitnan'); dyev = mean(R.r(:,R.t > 0 & R.t <= 3),2,'omitnan');
T = table(resp(keep),dyev(keep),R.snum(keep)-1,double(R.tord(keep)==2),double(R.is_lpsp(keep)),categorical(R.fly(keep)), ...
    'VariableNames',{'resp','dye_raw','stim','round2','lpsp','fly'});
T = T(~isnan(T.resp) & ~isnan(T.dye_raw),:);
T.dye = T.dye_raw - mean(T.dye_raw);
fprintf('%d stims (%d lpsp / %d empty flies), round 1: %d, round 2: %d\n',height(T), ...
    numel(unique(T.fly(T.lpsp==1))),numel(unique(T.fly(T.lpsp==0))),sum(T.round2==0),sum(T.round2==1));

%% fit every arrangement x shape
arr = {'pooled','continuous'};
res = struct();
aic_tab = table('Size',[0 5],'VariableTypes',{'string','string','double','double','double'}, ...
    'VariableNames',{'arrangement','shape','tau','loglik','AIC'});
for a = 1:2
    Ta = T;
    if a == 1; Ta.x = Ta.stim; else; Ta.x = Ta.stim + 29*Ta.round2; end
    % linear and log
    Ta.lx = log2(Ta.x+1);
    mL = fitlme(Ta,'resp ~ dye + x*lpsp + (1|fly)');
    mG = fitlme(Ta,'resp ~ dye + lx*lpsp + (1|fly)');
    % saturating: profile tau
    ll = nan(size(tau_grid));
    for k = 1:numel(tau_grid)
        Ta.sat = 1 - exp(-Ta.x/tau_grid(k));
        ll(k) = fitlme(Ta,'resp ~ dye + sat*lpsp + (1|fly)').LogLikelihood;
    end
    [~,kb] = max(ll); tau = tau_grid(kb);
    Ta.sat = 1 - exp(-Ta.x/tau);
    mS = fitlme(Ta,'resp ~ dye + sat*lpsp + (1|fly)');
    aic_tab = [aic_tab; {arr{a},"linear",nan,mL.LogLikelihood,mL.ModelCriterion.AIC}; ...
                        {arr{a},"log2",nan,mG.LogLikelihood,mG.ModelCriterion.AIC}; ...
                        {arr{a},"saturating",tau,mS.LogLikelihood,mS.ModelCriterion.AIC+2}]; %#ok<AGROW> % +2: tau
    res(a).T = Ta; res(a).mS = mS; res(a).ll = ll; res(a).tau = tau; res(a).mL = mL; res(a).mG = mG;

    cf = mS.Coefficients; ci = coefCI(mS); nm = cf.Name; b = cf.Estimate; V = mS.CoefficientCovariance;
    c = double(strcmp(nm,'lpsp')) + double(is_inter(nm,'sat')); % lpsp - empty at the plateau
    plat = c'*b; plat_se = sqrt(c'*V*c);
    res(a).plat = [plat plat-1.96*plat_se plat+1.96*plat_se 2*(1-normcdf(abs(plat/plat_se)))];
    fprintf('\n=== %s: saturating fit, tau = %g stims (loglik %.1f) ===\n',arr{a},tau,mS.LogLikelihood);
    for k = 1:height(cf)
        fprintf('  %-10s %+9.5f [%+9.5f %+9.5f]  p = %.2g\n',nm{k},b(k),ci(k,1),ci(k,2),cf.pValue(k));
    end
    fprintf('  lpsp - empty at plateau: %+.4f [%+.4f %+.4f], p = %.2g\n',res(a).plat);
    gi = is_inter(mG.Coefficients.Name,'lx'); li = is_inter(mL.Coefficients.Name,'x');
    fprintf('  log2 model: lx:lpsp %+.5f, p = %.2g | linear model: x:lpsp %+.6f, p = %.2g\n', ...
        mG.Coefficients.Estimate(gi),mG.Coefficients.pValue(gi),mL.Coefficients.Estimate(li),mL.Coefficients.pValue(li));
end
aic_tab.dAIC = aic_tab.AIC - min(aic_tab.AIC);
disp(aic_tab)
writetable(aic_tab,fullfile(out_dir,'tables','fig09_model_comparison.csv'));

%% figure
figure(9); clf; set(gcf,'Position',[20 40 1800 820],'Color','w')
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
gc = {ce,cl}; gl = {'empty>P2X2','lpsp>P2X2'};
for a = 1:2
    Ta = res(a).T; mS = res(a).mS; tau = res(a).tau;
    b = mS.Coefficients.Estimate; V = mS.CoefficientCovariance; nm = mS.Coefficients.Name;
    xmax = max(Ta.x); xg = 0:xmax; sg = 1 - exp(-xg/tau);
    % per-fly binned data
    nbin = ceil((xmax+1)/bin_size); bc = ((1:nbin)-0.5)*bin_size - 0.5;
    fl = unique(Ta.fly); is_l = arrayfun(@(f) Ta.lpsp(find(Ta.fly==f,1)),fl);
    M = nan(numel(fl),nbin);
    for k = 1:numel(fl)
        for bb = 1:nbin
            j = Ta.fly==fl(k) & floor(Ta.x/bin_size)+1 == bb;
            if any(j); M(k,bb) = mean(Ta.resp(j)); end
        end
    end
    % A: data + fit
    ax = nexttile(tl,(a-1)*4+1); hold(ax,'on')
    for l = 0:1
        X = design(nm,sg,l); p = X*b; e = sqrt(sum((X*V).*X,2));
        fill(ax,[xg fliplr(xg)]+1,[p+1.96*e; flipud(p-1.96*e)]',gc{l+1},'FaceAlpha',.12,'EdgeColor','none','HandleVisibility','off')
        plot(ax,xg+1,p,'-','Color',gc{l+1},'LineWidth',2,'DisplayName',gl{l+1})
        Mg = M(is_l==l,:); n = sum(~isnan(Mg),1);
        errorbar(ax,bc+1,mean(Mg,1,'omitnan'),std(Mg,0,1,'omitnan')./sqrt(n),'o','Color',gc{l+1},'MarkerFaceColor',gc{l+1},'CapSize',3,'HandleVisibility','off')
    end
    if a == 2; xline(ax,29.5,':k','HandleVisibility','off'); end
    yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[0 xmax+2])
    xlabel(ax,xlab(a)); ylabel(ax,'GCaMP \DeltaF/F, mean 0-10 s')
    legend(ax,'Location','northwest','FontSize',7)
    title(ax,sprintf('%s: data (5-stim bins, \\pm SEM) + saturating fit\n\\tau = %g stims',upper(arr{a}),tau),'FontSize',9)
    % B: lpsp - empty
    ax = nexttile(tl,(a-1)*4+2); hold(ax,'on')
    Xd = design(nm,sg,1) - design(nm,sg,0); p = Xd*b; e = sqrt(sum((Xd*V).*Xd,2));
    fill(ax,[xg fliplr(xg)]+1,[p+1.96*e; flipud(p-1.96*e)]','k','FaceAlpha',.12,'EdgeColor','none')
    plot(ax,xg+1,p,'k-','LineWidth',2)
    D = nan(1,nbin); De = D;
    for bb = 1:nbin
        x1 = M(is_l==1,bb); x0 = M(is_l==0,bb); x1 = x1(~isnan(x1)); x0 = x0(~isnan(x0));
        D(bb) = mean(x1)-mean(x0); De(bb) = sqrt(var(x1)/numel(x1)+var(x0)/numel(x0));
    end
    errorbar(ax,bc+1,D,De,'ko','MarkerFaceColor','w','CapSize',3)
    if a == 2; xline(ax,29.5,':k'); end
    yline(ax,0,':k'); xlim(ax,[0 xmax+2]); xlabel(ax,xlab(a)); ylabel(ax,'lpsp - empty, GCaMP \DeltaF/F')
    cf = mS.Coefficients; i0 = strcmp(nm,'lpsp'); i1 = is_inter(nm,'sat');
    title(ax,sprintf(['lpsp - empty (model 95%% CI)\nat stim 1: %+.4f (p = %.2g); extra rise: %+.4f (p = %.2g)\n' ...
        'at plateau: %+.4f [%+.4f %+.4f], p = %.2g'],cf.Estimate(i0),cf.pValue(i0),cf.Estimate(i1),cf.pValue(i1),res(a).plat),'FontSize',8)
    % C: tau profile
    ax = nexttile(tl,(a-1)*4+3); hold(ax,'on')
    plot(ax,tau_grid,res(a).ll - max(res(a).ll),'k.-')
    yline(ax,-1.92,':r','95% (\Delta logL = -1.92)','FontSize',7,'LabelHorizontalAlignment','left')
    xline(ax,tau,':k'); xlabel(ax,'\tau (stims)'); ylabel(ax,'log-likelihood - max')
    ylim(ax,[max(-12,min(res(a).ll-max(res(a).ll))) 0.5])
    title(ax,sprintf('profile over \\tau (best %g)',tau),'FontSize',9)
    % D: AIC comparison (both arrangements, all shapes) -- same in both rows
    ax = nexttile(tl,(a-1)*4+4); hold(ax,'on')
    lbl = strcat(aic_tab.arrangement,{' / '},aic_tab.shape);
    barh(ax,aic_tab.dAIC,'FaceColor',[.6 .6 .6]); set(ax,'YTick',1:height(aic_tab),'YTickLabel',lbl,'YDir','reverse','FontSize',8)
    for k = 1:height(aic_tab); text(ax,aic_tab.dAIC(k),k,sprintf('  %.1f',aic_tab.dAIC(k)),'FontSize',7,'VerticalAlignment','middle'); end
    xline(ax,2,':k'); xlabel(ax,'\DeltaAIC vs best (lower = better; < 2 ~ equivalent)')
    title(ax,'model comparison','FontSize',9)
end
title(tl,{'LPsP stimulation, walking 0.5 s pulses (fig08 ejection-trials data, rounds 1-2): saturating mixed model', ...
    'resp ~ dye + sat * lpsp + (1 | fly),  sat = 1 - exp(-x / \tau);  resp = ejection-side GCaMP mean 0-10 s; dye = ejection-side dye 0-3 s'})
exportgraphics(gcf,fullfile(out_dir,'fig09_figure3_p2x2_lme.pdf'))
if exist('preview_png','var'); exportgraphics(gcf,preview_png,'Resolution',110); end
fprintf('figure written to %s\n',out_dir);

%% functions
function X = design(nm,sg,l)
    % fixed-effect design rows over the stim grid (sat values sg), dye at its mean
    sg = sg(:); n = numel(sg); X = zeros(n,numel(nm));
    for k = 1:numel(nm)
        switch nm{k}
            case '(Intercept)';           X(:,k) = 1;
            case 'dye';                   X(:,k) = 0;
            case 'sat';                   X(:,k) = sg;
            case 'lpsp';                  X(:,k) = l;
            case {'sat:lpsp','lpsp:sat'}; X(:,k) = sg*l;
        end
    end
end

function tf = is_inter(nm,v)
    % the lpsp x v interaction term, whichever order MATLAB names it in
    tf = strcmp(nm,[v ':lpsp']) | strcmp(nm,['lpsp:' v]);
end

function s = xlab(a)
    if a == 1; s = 'stim # within trial (rounds 1 and 2 pooled)';
    else;      s = 'stim # (round 1: 1-29, round 2: 30-58)'; end
end

function D = collect(D,ad,max_stims)
    % one row per delivered stim, as in fig_figure3_p2x2_traces.m (walking):
    % ejection-side GCaMP / dye traces on D.t (baseline-subtracted), stim #,
    % fly, genotype, trial key, trial side, confirmed-ejection call
    D.g = []; D.r = []; D.snum = []; D.fly = {}; D.is_lpsp = []; D.ej = []; D.tkey = []; D.side = [];
    bidx = D.t >= D.base(1) & D.t < D.base(2);
    for i = 1:numel(ad)
        tr = ad(i);
        if ~isfield(tr,'atp') || ~isstruct(tr.atp) || ~isfield(tr.atp,'d') || isempty(tr.atp.d); continue; end
        if ~isfield(tr.ft,'stims') || isempty(tr.ft.stims) || max(tr.ft.stims) <= 0; continue; end
        st = tr.ft.stims(:) > 0.5*max(tr.ft.stims);
        on = tr.ft.xf(diff([false;st])==1);
        if numel(on) > max_stims; continue; end
        is_lpsp = contains(tr.meta,'_lpsp_');
        p = strsplit(tr.meta,'\'); p(cellfun(@isempty,p)) = [];
        fi = find(startsWith(p,'fly'),1,'last'); n_ = regexp(lower(p{fi}),'^fly\s*(\d+)','tokens','once');
        fid = sprintf('%s|fly %s',p{fi-1},n_{1});
        n  = min(numel(tr.ft.xb),size(tr.im.d,2)); xb = tr.ft.xb(1:n); xb = xb(:);
        A  = tr.atp.d(:,1:n); ev = zeros(size(A,1),1);
        for s = 1:numel(on)
            pw = xb > on(s) & xb <= on(s)+D.side_win; bw = xb >= on(s)+D.base(1) & xb < on(s);
            if any(pw) && any(bw); ev = ev + max(mean(A(:,pw),2) - mean(A(:,bw),2),0); end
        end
        half = size(A,1)/2;
        if sum(ev(1:half)) >= sum(ev(half+1:end)); rows = 1:half; side = 1; else; rows = half+1:2*half; side = 2; end
        tk = regexp(tr.meta,'(\d{8})-(\d+)','tokens','once'); tkey = str2double(tk{1})*1e4 + str2double(tk{2});
        asum = sum(A,1)';
        gs = mean(tr.im.d(rows,1:n),1)'; rs = sum(A(rows,:),1)';
        for s = 1:numel(on)
            g = interp1(xb-on(s),gs,D.t,'linear',nan); r = interp1(xb-on(s),rs,D.t,'linear',nan);
            if all(isnan(g(bidx))); continue; end
            D.g(end+1,:) = g - mean(g(bidx),'omitnan'); D.r(end+1,:) = r - mean(r(bidx),'omitnan'); %#ok<AGROW>
            D.snum(end+1,1) = s; D.fly{end+1,1} = fid; D.is_lpsp(end+1,1) = is_lpsp; %#ok<AGROW>
            bw = xb >= on(s)-5 & xb < on(s); pw = xb > on(s) & xb <= on(s)+3;
            D.ej(end+1,1) = (median(asum(pw)) - median(asum(bw))) > 3*std(asum(bw)); %#ok<AGROW>
            D.tkey(end+1,1) = tkey; D.side(end+1,1) = side; %#ok<AGROW>
        end
    end
    D.is_lpsp = logical(D.is_lpsp);
end

function keep = trial_has_ej(D)
    g = findgroups(D.fly,D.tkey);
    keep = ismember(g,unique(g(logical(D.ej))));
end

function D = keep_rows(D,m)
    m = logical(m);
    for f = {'g','r','snum','fly','is_lpsp','ej','tkey','side'}
        D.(f{1}) = D.(f{1})(m,:);
    end
end

function tord = trial_round(D)
    % which trial (1st, 2nd, ...) of the fly's same-side run each stim is
    % from; the run restarts when the trial's ejection side switches
    tord = nan(numel(D.fly),1);
    for f = unique(D.fly)'
        j = find(strcmp(D.fly,f{1}));
        [~,~,it] = unique(D.tkey(j));
        tside = accumarray(it,D.side(j),[],@(x) x(1));
        seg = cumsum([1; diff(tside) ~= 0]);
        for s = 1:seg(end)
            js = j(seg(it) == s);
            [~,~,tt] = unique(D.tkey(js)); tord(js) = tt;
        end
    end
end
