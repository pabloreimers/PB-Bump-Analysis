%% fig_grab_vs_lpsp_ca_nonlinearity
% Question: can EPG>GRAB(DA2m) be produced from LPsP>syt7f calcium by a
% static (pointwise) nonlinearity, or is DA release/sensing doing something
% more intricate?
%
% The two sensors were imaged in DIFFERENT flies (lpsp_cl_data_20240206.mat),
% so there is no frame-by-frame mapping to fit. Instead we test statistics
% that a static nonlinearity must preserve:
%
%   * Rank (Spearman) correlations between glomeruli are exactly invariant to
%     ANY monotonic pointwise nonlinearity -- even one that differs per
%     glomerulus -- and to a shared linear temporal filter (LN model) only
%     to the extent that inputs are already correlated (filtering ~identical
%     traces keeps them ~identical). So if LPsP glomeruli are ~0.95
%     correlated and GRAB glomeruli are not, no static NL can bridge them.
%   * Sensor shot noise also decorrelates glomeruli, so we use a lagged
%     estimator: im.d was pre-smoothed with a 5-frame gaussian (std ~1
%     frame), so noise is uncorrelated beyond ~3 frames. Cross-covariance at
%     lags 4-6 frames, normalized by autocovariance at the same lags, gives a
%     noise-free estimate of the signal correlation.
%   * Brief (1-3 frame) spatially uniform artifact flashes exist in both
%     sensors and are removed first. They are global, so leaving them in
%     would INFLATE GRAB correlations -- removing them is conservative.

%% load
d = load(fullfile('data','lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); % 'EPG_syt7f' folders are a naming typo: those are GRAB(DA2m) trials
isLpsp = contains(meta,'LPsP_syt7f');
isDark = contains(meta,'dark');
grp    = {find(isLpsp), find(isGrab)};
names  = {'LPsP > syt7f','EPG > GRAB(DA2m)'};
fr     = 9.1;   % imaging frame rate (Hz)
lag    = 10;    % ft samples between |vr| and fluorescence, as in lpsp_cl_script_minimal

%% per-trial: clean artifacts, noise-corrected inter-glomerular rank corr vs distance, speed tuning
spE  = [0 .25 .5 1 1.5 2 3 4 6];
maxd = 15;
res  = struct();
for g = 1:2
    N = numel(grp{g});
    res(g).Cdist = nan(N,maxd+1);   % noise-corrected corr vs bin distance (within side)
    res(g).Craw  = nan(N,maxd+1);
    res(g).TC    = nan(N,numel(spE)-1);
    res(g).D     = cell(N,1);
    for n = 1:N
        i = grp{g}(n);
        D = d(i).im.d;  T = size(D,2);
        D(:,find_artifacts(D)) = nan;
        res(g).D{n} = D;

        [Cs,C0] = lagged_rank_corr(D,4:6);
        for k = 0:maxd
            cs = []; c0 = [];
            for side = 0:1
                idx = side*16 + (1:16);
                cs = [cs; diag(Cs(idx,idx),k)]; c0 = [c0; diag(C0(idx,idx),k)]; %#ok<AGROW>
            end
            res(g).Cdist(n,k+1) = median(cs);
            res(g).Craw(n,k+1)  = median(c0);
        end

        xf = d(i).ft.xf; xb = linspace(xf(1),xf(end),T);
        rs = d(i).ft.r_speed; rs = [rs(lag+1:end); repmat(rs(end),lag,1)];
        sp = interp1(xf,abs(rs),xb);
        m  = mean(D,1,'omitnan'); b = discretize(sp,spE);
        for k = 1:numel(spE)-1
            sel = b==k & ~isnan(m);
            if sum(sel) > 20, res(g).TC(n,k) = median(m(sel)); end
        end
        res(g).TC(n,:) = (res(g).TC(n,:)-res(g).TC(n,1)) / (res(g).TC(n,end-1)-res(g).TC(n,1));
    end
end

%% best-case static NL: quantile-map pooled (per-trial normalized) Ca values onto GRAB values
nz = @(D) (D - median(D,'all','omitnan')) / iqr(D(:));
pool = @(C) cell2mat(cellfun(@(x) reshape(x(~isnan(x)),[],1), cellfun(nz,C,'UniformOutput',false),'UniformOutput',false));
ca = pool(res(1).D); gr = pool(res(2).D);
q  = linspace(.001,.999,999); [qc,ui] = unique(quantile(ca,q)); qg = quantile(gr,q); qg = qg(ui);
fNL = @(x) interp1(qc,qg,x,'linear','extrap');

pk = @(D) (max(D,[],1)-median(D,1)) ./ (mean(abs(D),1)+eps);   % spatial peak contrast per frame
src = {res(1).D, cellfun(@(x) fNL(nz(x)),res(1).D,'UniformOutput',false), res(2).D};
lab = {'LPsP Ca','Ca \rightarrow static NL','EPG GRAB'};
PK  = cell(1,3);
for s = 1:3
    for n = 1:numel(src{s})
        D = smoothdata(nz(src{s}{n}),2,'movmean',5,'omitnan');
        m = mean(D,1); top = m > prctile(m,90);
        PK{s}(n) = median(pk(D(:,top)),'omitnan');
    end
end

%% figure
c = [0.85 0.33 0.10; 0.00 0.45 0.74];
figure(1); clf; set(gcf,'Color','w','Position',[100 100 1300 700]);

subplot(2,3,[1 2]);  % example heatmaps
n1 = 19; n2 = 6;
D = [res(1).D{n1}(:,round(100*fr):round(250*fr)); nan(2,round(150*fr)+1); res(2).D{n2}(:,round(100*fr):round(250*fr))];
imagesc([100 250],[],D./prctile(abs(D),99.5,2)); colormap(gca,parula); caxis([0 1]);
yticks([16 50]); yticklabels(names); xlabel('time (s)'); title('example trials (each row scaled to its 99.5th pct)')

subplot(2,3,3);  % noise-corrected corr vs distance
hold on
for g = 1:2
    plot(0:maxd, median(res(g).Cdist),'-','Color',c(g,:),'LineWidth',2)
    plot(0:maxd, median(res(g).Craw), ':','Color',c(g,:),'LineWidth',1.5)
end
xlabel('distance along PB (bins)'); ylabel('rank corr between glomeruli')
legend({[names{1} ' noise-corr.'],[names{1} ' raw'],[names{2} ' noise-corr.'],[names{2} ' raw']},'Location','southwest','Box','off')
ylim([0 1]); title({'static NL cannot change these'})

subplot(2,3,4);  % per-trial far-pair corr
far = 6:maxd;
x = {mean(res(1).Cdist(:,far+1),2), mean(res(2).Cdist(:,far+1),2)};
hold on
for g = 1:2, swarmchart(g*ones(size(x{g})),x{g},20,c(g,:),'filled'); end
xticks(1:2); xticklabels(names); ylabel('noise-corr. rank corr, \geq6 bins apart'); ylim([0 1.05]);
title(sprintf('ranksum p = %.1g',ranksum(x{1},x{2})))

subplot(2,3,5);  % speed tuning
hold on
spc = (spE(1:end-1)+spE(2:end))/2;
for g = 1:2
    errorbar(spc, median(res(g).TC,'omitnan'), iqr(res(g).TC)/2, '-o','Color',c(g,:),'LineWidth',1.5)
end
xlabel('|rotational speed| (rad/s)'); ylabel('PB-mean dF/F (normalized)'); title('global speed tuning'); legend(names,'Box','off','Location','southeast')

subplot(2,3,6);  % peak contrast
hold on
cc = [c(1,:); 0.6 0.6 0.6; c(2,:)];
for s = 1:3, swarmchart(s*ones(size(PK{s})),PK{s},20,cc(s,:),'filled'); end
xticks(1:3); xticklabels(lab); ylabel('(max - median)/mean |profile|'); title('spatial peakiness, top-10% frames')

fontsize(gcf,12,'pixels')

%% export
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_vs_lpsp_ca_nonlinearity.png'),'Resolution',150);

%% could a sensor floor (thresholding) explain GRAB's low inter-glomerular correlation?
% A floor is monotonic but NOT strictly so: below it only noise is left,
% which does decorrelate glomeruli -- so the invariance argument above
% doesn't cover it. Two tests:
%  (a) data: noise-corrected correlation in frames where DA should be far
%      above any floor (|vr| > 1 or > 2 rad/s). Selecting on behavior, not
%      on the signal itself (selecting on PB level restricts the shared
%      component's range and lowers correlations for both sensors).
%      Floor hypothesis predicts GRAB correlations RISE toward calcium's.
%  (b) simulation: LPsP calcium through a common floor at the p-th
%      percentile, with/without GRAB-like noise (SNR 1, 5-frame gaussian).
spd_thr = [0 1 2];
FLc = cell(1,2);
for g = 1:2
    FLc{g} = nan(numel(grp{g}),numel(spd_thr));
    for n = 1:numel(grp{g})
        i = grp{g}(n); D = res(g).D{n}; T = size(D,2); bad = all(isnan(D),1);
        xf = d(i).ft.xf; xb = linspace(xf(1),xf(end),T);
        rs = d(i).ft.r_speed; rs = [rs(lag+1:end); repmat(rs(end),lag,1)]; sp = interp1(xf,abs(rs),xb);
        for k = 1:numel(spd_thr)
            sel = sp >= spd_thr(k) & ~bad;
            if sum(sel) > 300, FLc{g}(n,k) = masked_lag_corr(D,sel); end
        end
    end
end
fprintf('\nnoise-corrected rank corr (bins 6-10 apart), by |vr| threshold %s rad/s:\n', mat2str(spd_thr));
for g = 1:2, fprintf('  %-18s %s\n', names{g}, mat2str(median(FLc{g},'omitnan'),2)); end

rng(2); thr = [0 25 50 75 90];
SIMc = nan(numel(grp{1}),numel(thr),2);   % (:,:,1) all frames, (:,:,2) |vr|>1
for n = 1:numel(grp{1})
    i = grp{1}(n); D = res(1).D{n}; T = size(D,2); bad = all(isnan(D),1);
    xf = d(i).ft.xf; xb = linspace(xf(1),xf(end),T);
    rs = d(i).ft.r_speed; rs = [rs(lag+1:end); repmat(rs(end),lag,1)]; sp = interp1(xf,abs(rs),xb);
    for a = 1:numel(thr)
        th = prctile(D(:),thr(a)); Df = max(D,th);
        nzz = smoothdata(randn(size(D)),2,'gaussian',5); nzz = nzz./std(nzz,0,2);
        Dn = Df + nzz.*max(std(Df,0,2,'omitnan'),eps); Dn(:,bad) = nan;
        SIMc(n,a,1) = masked_lag_corr(Dn,~bad);
        if sum(sp>1 & ~bad) > 300, SIMc(n,a,2) = masked_lag_corr(Dn,sp>1 & ~bad); end
    end
end
fprintf('LPsP Ca + floor at pct %s + noise: all frames %s, |vr|>1 %s\n', mat2str(thr), ...
    mat2str(median(SIMc(:,:,1),'omitnan'),2), mat2str(median(SIMc(:,:,2),'omitnan'),2));

figure(5); clf; set(gcf,'Color','w','Position',[100 100 900 380]);
subplot(1,2,1); hold on
for g = 1:2, errorbar(spd_thr, median(FLc{g},'omitnan'), iqr(FLc{g})/2, '-o','Color',c(g,:),'MarkerFaceColor',c(g,:),'LineWidth',1.5); end
xlabel('only frames with |vr| \geq (rad/s)'); ylabel('noise-corr. rank corr, 6-10 bins apart'); ylim([-0.2 1]); xlim([-.3 2.3])
legend(names,'Box','off','Location','southwest'); title('data: faster turning (above any floor)')
subplot(1,2,2); hold on
plot(thr, median(SIMc(:,:,1),'omitnan'),'-o','Color',c(1,:),'LineWidth',1.5)
plot(thr, median(SIMc(:,:,2),'omitnan'),'--o','Color',c(1,:),'LineWidth',1.5)
yline(median(FLc{2}(:,1),'omitnan'),'-','Color',c(2,:),'LineWidth',1.5); yline(median(FLc{2}(:,2),'omitnan'),'--','Color',c(2,:),'LineWidth',1.5)
xlabel('floor (percentile of Ca values clipped)'); ylabel('noise-corr. rank corr, 6-10 bins apart'); ylim([-0.2 1])
legend({'Ca + floor + noise, all','Ca + floor + noise, |vr|>1','real GRAB, all','real GRAB, |vr|>1'},'Box','off','Location','southwest')
title('simulation: floor on LPsP calcium')
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_floor_test.png'),'Resolution',150);

%% functions
function bad = find_artifacts(D)
% brief, global upward flashes: fast residual of the PB-mean trace > 6 robust SDs, dilated +/-3 frames
m   = mean(D,1);
r   = m - movmedian(m,21);
bad = r > 6*1.4826*mad(r,1);
bad = conv(double(bad),ones(1,7),'same') > 0;
end

function [Cs,C0] = lagged_rank_corr(D,lags)
% rank-transform each glomerulus, then estimate signal correlation from
% cross-covariance at lags where noise is uncorrelated (see header)
T = size(D,2); Z = nan(size(D));
for j = 1:size(D,1)
    ok = ~isnan(D(j,:)); Z(j,ok) = tiedrank(D(j,ok));
    Z(j,:) = (Z(j,:) - mean(Z(j,:),'omitnan')) / std(Z(j,:),'omitnan');
end
C = zeros(size(D,1));
for L = lags
    A1 = Z(:,1:T-L); A2 = Z(:,1+L:T); ok = ~any(isnan([A1;A2]),1);
    Cl = A1(:,ok)*A2(:,ok)'/sum(ok); C = C + (Cl+Cl')/2;
end
dg = diag(C); Cs = C ./ sqrt(dg*dg');
ok = ~any(isnan(Z),1); C0 = corr(Z(:,ok)');
end

function [far, raw_far] = masked_lag_corr(D, sel, lags)
% lagged_rank_corr restricted to frames in sel (1xT logical): ranks over the
% selected frames, lagged products only over pairs where both frames are
% selected. Returns median over within-side pairs 6-10 bins apart.
if nargin < 3, lags = 4:6; end
T = size(D,2); n = size(D,1); Z = nan(size(D));
for j = 1:n
    ok = sel & ~isnan(D(j,:)); Z(j,ok) = tiedrank(D(j,ok));
    Z(j,:) = (Z(j,:) - mean(Z(j,:),'omitnan')) / std(Z(j,:),'omitnan');
end
C = zeros(n); cnt = 0;
for L = lags
    A1 = Z(:,1:T-L); A2 = Z(:,1+L:T); ok = ~any(isnan([A1;A2]),1);
    if sum(ok) < 50, continue, end
    Cl = A1(:,ok)*A2(:,ok)'/sum(ok); C = C + (Cl+Cl')/2; cnt = cnt+1;
end
C = C/cnt; dg = diag(C); Cs = C ./ sqrt(dg*dg');
ok = ~any(isnan(Z),1); if sum(ok) > 2, C0 = corr(Z(:,ok)'); else, C0 = nan(n); end
farv = []; rawf = [];
for side = 0:1
    idx = side*16 + (1:16);
    for j = 1:16
        for k = j+6:min(j+10,16)
            farv(end+1) = Cs(idx(j),idx(k)); rawf(end+1) = C0(idx(j),idx(k)); %#ok<AGROW>
        end
    end
end
far = median(farv,'omitnan'); raw_far = median(rawf,'omitnan');
end
