%% fig_grab_vs_lpsp_speed_dropouts
% Q: how often is PB fluorescence LOW even though the fly is rotating, for
% LPsP>syt7f vs EPG>GRAB(DA2m)? (Expectation going in: LPsP calcium is high
% every time the fly turns, GRAB has more "missed" turns.)
%
% Method (per trial, lpsp_cl_data_20240206.mat, closed loop + dark):
%   1) F = PB-mean (or PB-max after 3-bin spatial smoothing) dF/F, brief
%      global artifact flashes removed, drift slower than ~60 s removed
%      (running median) so baseline wander/bleaching isn't scored as a miss.
%   2) Speed model: LN prediction of F from |vr| -- ridge-regression linear
%      filter on |vr| at lags -0.5..3 s, then a binned static nonlinearity,
%      each sensor/trial fit separately (so sensor kinetics and gain are
%      absorbed). 2-fold cross-validated over alternating 60 s blocks, so a
%      dip can't be fit away by the model that scores it.
%   3) Residual r = F - prediction, smoothed over 1 s. Sensor noise SD of
%      that smoothed residual estimated per trial from the fast part of the
%      residual (total variance minus autocovariance at lags 4-6 frames --
%      the data were pre-smoothed ~1 frame, so noise is uncorrelated there),
%      scaled for the 1 s smoothing.
%   4) "Moving" = |vr| > 1 rad/s shifted by the sensor's own lag (xcorr peak
%      of PB-mean vs |vr|). A moving frame is flagged LOW when the smoothed
%      residual is both significant (< -3 noise SD) and substantial (below
%      half of the predicted speed-driven increment above rest).
%   Reported: fraction of moving time flagged, and flagged runs >= 0.5 s per
%   minute of moving time. The behavioral "moving" definition gives both
%   sensors the same opportunity (~17-18% of time each).
%
% Result when written: the opposite of the expectation -- LPsP calcium has
% MORE low-while-turning time than GRAB (PB-mean ~12% vs ~8% of moving
% time; ~5.5 vs ~2.8 dropouts per moving minute), in closed loop and dark
% separately, while |vr| predicts both about equally well (noise-corrected
% CV R2 ~0.6). GRAB tracks rotational speed tightly turn by turn; LPsP
% calcium is broader and also has elevations without turning. BUT at the
% trial level GRAB has more trials that barely follow turning at all (7/41
% vs 2/43 with turning-vs-rest d < 0.5, Fisher p ~0.09) -- the per-trial
% model can't flag those (it learns a flat prediction), see the trial-level
% section. So: GRAB misses fewer individual turns, but more whole trials
% look unresponsive (sensor/expression or real -- can't tell here). Caveat: the
% two sensors are in different flies, and LPsP calcium may carry other
% variables (e.g. forward walking) that a |vr|-only model scores as misses.

%% load
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_data_20240206.mat')), data_dir = '.data'; else, data_dir = 'data'; end
d = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isDark = contains(meta,'dark');
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); % 'EPG_syt7f' folders are a naming typo for GRAB(DA2m)
isLpsp = contains(meta,'LPsP_syt7f');
grp    = {find(isLpsp), find(isGrab)};
names  = {'LPsP > syt7f','EPG > GRAB(DA2m)'};
mnames = {'PB-mean','PB-max'};
c      = [0.85 0.33 0.10; 0.00 0.45 0.74];
fr     = 9.1;

opts = struct('lagmin',-0.5,'lagmax',3,'block_s',60,'lambda',50,'nbins',10);
hp_s = 60; smw = round(1*fr); mov_thr = 1; z_thr = 3; frac_thr = 0.5; min_run = round(0.5*fr);
nz = smoothdata(randn(1,1e6),'gaussian',5); kfac = var(movmean(nz,smw))/var(nz); % noise variance kept by the 1 s smoothing

%% per trial: model, residual noise, flags
R = cell(1,2);
for g = 1:2
    for n = 1:numel(grp{g})
        i = grp{g}(n); s = d(i); D = s.im.d; T = size(D,2);
        xf = s.ft.xf; xb = linspace(xf(1),xf(end),T);
        m = mean(D,1); rr = m - movmedian(m,21); bad = rr > 6*1.4826*mad(rr,1); bad = conv(double(bad),ones(1,7),'same') > 0;
        sp = interp1(xf,abs(smoothdata(s.ft.r_speed,'gaussian',30)),xb);
        Fs = {mean(D,1), max(smoothdata(D,1,'movmean',3),[],1)};
        % sensor lag from PB-mean vs |vr|
        lg = -round(1*fr):round(3*fr); x = Fs{1}; x(bad) = mean(x(~bad));
        xc = xcorr(x-mean(x), sp-mean(sp), lg(end), 'coeff'); [~,j] = max(xc(lg+lg(end)+1)); lag = lg(j);
        moving = circshift(sp,lag) > mov_thr & ~bad;
        for mt = 1:2
            F = Fs{mt}; F(bad) = nan;
            F = F - movmedian(F,round(hp_s*fr),'omitnan') + median(F,'omitnan');
            [p,info] = speed_model(F, sp, bad, fr, opts);
            r = F - p; ok = ~isnan(r); x = r(ok) - mean(r(ok));
            vnoise = max(var(x) - mean(arrayfun(@(L) mean(x(1:end-L).*x(1+L:end)),4:6)), eps);
            xF = F(ok) - mean(F(ok)); sigfrac = mean(arrayfun(@(L) mean(xF(1:end-L).*xF(1+L:end)),4:6))/var(xF);
            sdn = sqrt(vnoise*kfac);
            rs = movmean(r,smw,'omitnan');
            inc = p - prctile(p(sp<0.3 & ~bad),50);
            mv = moving & ~isnan(rs);
            flag = mv & rs < -z_thr*sdn & rs < -frac_thr*max(inc,0);
            dd = diff([0 flag 0]); st = find(dd==1); en = find(dd==-1)-1;
            R{g}(n,mt) = struct('R2',info.R2,'R2c',info.R2/sigfrac,'fmov',mean(mv),'fflag',sum(flag)/max(sum(mv),1), ...
                'rate',sum((en-st+1) >= min_run)/max(sum(mv)/fr/60,eps),'dark',isDark(i),'F',F,'p',p,'sp',sp,'flag',flag);
        end
    end
end

%% stats
for mt = 1:2
    fprintf('\n%s:\n', mnames{mt});
    for g = 1:2
        H = R{g}(:,mt);
        fprintf('  %-18s n=%2d  noise-corr CV R2 %.2f [%.2f %.2f]  moving %.2f  flagged low %.3f [%.3f %.3f]  dropouts/moving min %.2f [%.2f %.2f]\n', ...
            names{g}, numel(H), median([H.R2c]), prctile([H.R2c],[25 75]), median([H.fmov]), ...
            median([H.fflag]), prctile([H.fflag],[25 75]), median([H.rate]), prctile([H.rate],[25 75]));
    end
    a = R{1}(:,mt); b = R{2}(:,mt);
    fprintf('  ranksum: flagged fraction p=%.2g, dropout rate p=%.2g, noise-corr R2 p=%.2g\n', ...
        ranksum([a.fflag],[b.fflag]), ranksum([a.rate],[b.rate]), ranksum([a.R2c],[b.R2c]));
    for cond = [false true]
        aa = a([a.dark]==cond); bb = b([b.dark]==cond);
        fprintf('    %-11s flagged low: %.3f (n=%d) vs %.3f (n=%d), p=%.2g\n', char(string(cond)), median([aa.fflag]), numel(aa), median([bb.fflag]), numel(bb), ranksum([aa.fflag],[bb.fflag]));
    end
end

%% trial-level check: trials where the sensor barely follows turning at all
% The per-trial model absorbs a trial that ignores turns entirely (it fits
% a flat prediction, so nothing gets flagged). Score those separately:
% effect size d = (mean F while |vr| > 1 at ~0.4 s lag - mean F at rest) /
% rest SD; d < 0.5 = weak turning response for the whole trial.
for g = 1:2
    H = R{g}(:,1);
    dval = arrayfun(@(h) (mean(h.F(circshift(h.sp,4) > 1),'omitnan') - mean(h.F(h.sp < 0.3),'omitnan')) / std(h.F(h.sp < 0.3),'omitnan'), H);
    R_d{g} = dval; %#ok<SAGROW>
    fprintf('%-18s trials with weak turning response (d < 0.5): %d/%d; median d %.2f\n', names{g}, sum(dval < 0.5), numel(dval), median(dval));
end
[~,p_fish] = fishertest([sum(R_d{1}<0.5) sum(R_d{1}>=0.5); sum(R_d{2}<0.5) sum(R_d{2}>=0.5)]);
fprintf('Fisher p = %.2g\n', p_fish);

%% rises WITHOUT turning: fluorescence well above the speed prediction while the fly is still
% Same residual machinery, flipped: smoothed residual > +3 noise SD AND > half
% the trial's typical speed-driven increment (90th pct of prediction minus
% rest), with |vr| < 0.5 rad/s over the preceding 2 s and the current 0.5 s,
% for >= 1 s. Then the "clean" subset: fly fully still (|vr| < 0.5 and
% |forward| < 20% of the trial's 90th pct, -5..+3 s) AND a genuine rise
% (trial-normalized dF/F up > 0.2 from 2-3 s before the peak -- excludes
% slow-decay-after-walking residuals, ~30% of raw events).
% Caveat: FicTrac can't see grooming or legs off the ball, so "still" means
% the ball isn't moving, not that the fly isn't doing anything.
still_thr = 0.5; rise_min_run = round(1*fr); pre = round(10*fr); post = round(6*fr);
RISE = cell(1,2); CLEAN = cell(1,2);
for g = 1:2
    nstill = zeros(numel(grp{g}),1); nev = nstill;
    for n = 1:numel(grp{g})
        i = grp{g}(n); s = d(i); h = R{g}(n,1); T = numel(h.F);
        fs = interp1(s.ft.xf, smoothdata(s.ft.f_speed,'gaussian',30), linspace(s.ft.xf(1),s.ft.xf(end),T));
        r = h.F - h.p; ok = ~isnan(r); x = r(ok) - mean(r(ok));
        sdn = sqrt(max(var(x) - mean(arrayfun(@(L) mean(x(1:end-L).*x(1+L:end)),4:6)), eps)*kfac);
        rs = movmean(r, smw, 'omitnan');
        base = prctile(h.F(h.sp < 0.3),50); top = prctile(h.F,90);
        inc90 = prctile(h.p,90) - prctile(h.p(h.sp < 0.3),50);
        still = movmax(h.sp,[round(2*fr) round(.5*fr)]) < still_thr;
        flag = still & rs > z_thr*sdn & rs > frac_thr*inc90 & ~isnan(rs);
        dd = diff([0 flag 0]); st = find(dd==1); en = find(dd==-1)-1; keep = (en-st+1) >= rise_min_run;
        st = st(keep); en = en(keep);
        nstill(n) = sum(still)/fr/60; nev(n) = numel(st);
        fthr = 0.2*prctile(abs(fs),90);
        for e = 1:numel(st)
            [~,j] = max(rs(st(e):en(e))); pk = st(e)+j-1;
            if pk-pre < 1 || pk+post > T, continue, end
            dv = median((h.F(pk-3:pk+3)-base)/(top-base),'omitnan') - median((h.F(pk-round(3*fr):pk-round(2*fr))-base)/(top-base),'omitnan');
            w = pk-round(5*fr):pk+round(3*fr);
            isclean = max(h.sp(w)) < still_thr && max(movmean(abs(fs(w)),round(fr))) < fthr && dv > 0.2;
            RISE{g}(end+1,:) = [n pk dv isclean];
        end
    end
    rate = nev./max(nstill,eps);
    CLEAN{g} = RISE{g}(RISE{g}(:,4)==1,:);
    fprintf('%-18s rises without turning: %d events in %d/%d trials, %.2f [%.2f %.2f] per still minute; genuine rises %.0f%%; clean (fully still + genuine) %d in %d trials\n', ...
        names{g}, sum(nev), sum(nev>0), numel(nev), median(rate), prctile(rate,[25 75]), 100*mean(RISE{g}(:,3) > 0.2,'omitnan'), ...
        size(CLEAN{g},1), numel(unique(CLEAN{g}(:,1))));
    RATE_rise{g} = rate; %#ok<SAGROW>
end
fprintf('ranksum rise rate p = %.2g\n', ranksum(RATE_rise{1},RATE_rise{2}));
fprintf('clean LPsP calcium rises (largest first):\n');
cl = sortrows(CLEAN{1},-3);
for k = 1:size(cl,1)
    fprintf('  %s  t = %.1f s  (rise %.2f)\n', regexp(meta{grp{1}(cl(k,1))},'\d{8}-\d+_[^\\]*$','match','once'), cl(k,2)/fr, cl(k,3));
end

figure(3); clf; set(gcf,'Color','w','Position',[50 50 1600 1100]); tl = tiledlayout(4,4,'TileSpacing','compact');
for k = 1:min(8,size(cl,1))
    n = cl(k,1); pk = cl(k,2); i = grp{1}(n); h = R{1}(n,1); D = d(i).im.d; T = numel(h.F);
    fs = interp1(d(i).ft.xf, smoothdata(d(i).ft.f_speed,'gaussian',30), linspace(d(i).ft.xf(1),d(i).ft.xf(end),T));
    a = max(1,pk-round(15*fr)); b = min(T,pk+round(10*fr)); t = (a:b)/fr;
    nexttile(tl,(floor((k-1)/4)*2)*4 + mod(k-1,4)+1);
    imagesc(t,1:32,D(:,a:b)); caxis(prctile(D(:,a:b),[2 99.5],'all')); xline(pk/fr,'w--'); set(gca,'XTick',[]); ylabel('PB bin')
    title(sprintf('%s  t=%.0f s', strrep(regexp(meta{i},'\d{8}-\d+_[^\\]*$','match','once'),'_',' '), pk/fr),'FontSize',8)
    nexttile(tl,(floor((k-1)/4)*2+1)*4 + mod(k-1,4)+1); hold on
    yyaxis left; plot(t,h.F(a:b),'k-'); plot(t,h.p(a:b),'-','Color',[0 .6 0]); ylabel('PB-mean dF/F')
    yyaxis right; plot(t,h.sp(a:b),'m-'); plot(t,abs(fs(a:b))/max(abs(fs(a:b))+eps)*max(h.sp(a:b)+eps),'c-'); xline(pk/fr,'k--'); xlabel('time (s)')
end
title(tl,{'LPsP > syt7f calcium rises while the fly is still (ball not moving)', ...
    'black = PB-mean dF/F, green = speed prediction (CV), magenta = |vr| (rad/s), cyan = |forward| (scaled to |vr| axis)'})
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_lpsp_rises_without_turning.png'),'Resolution',150);

%% figure: summary
figure(1); clf; set(gcf,'Color','w','Position',[100 100 1300 380]);
lab = {'flagged-low fraction of moving time','dropouts per moving minute','noise-corrected CV R^2 (speed model)'};
fld = {'fflag','rate','R2c'};
for k = 1:3
    subplot(1,3,k); hold on
    for g = 1:2
        H = R{g}(:,1); v = [H.(fld{k})]; dk = [H.dark];
        swarmchart(g*ones(1,sum(~dk))-0.12, v(~dk), 25, c(g,:), 'filled');
        swarmchart(g*ones(1,sum(dk))+0.12, v(dk), 25, c(g,:));
        plot(g+[-.3 .3], median(v)*[1 1], 'k-', 'LineWidth', 2)
    end
    xticks(1:2); xticklabels(names); xlim([.5 2.5]); ylabel(lab{k})
    a = R{1}(:,1); b = R{2}(:,1); title(sprintf('PB-mean, ranksum p = %.1g', ranksum([a.(fld{k})],[b.(fld{k})])))
end
sgtitle('filled = closed loop, open = dark; one dot per trial')
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_vs_lpsp_speed_dropouts.png'),'Resolution',150);

%% figure: examples (trials closest to each sensor's median flagged fraction), 200-320 s
figure(2); clf; set(gcf,'Color','w','Position',[100 100 1500 900]); tl = tiledlayout(4,1,'TileSpacing','compact');
for g = 1:2
    H = R{g}(:,1); v = [H.fflag]; [~,o] = sort(abs(v - median(v)));
    for k = 1:2
        h = H(o(k)); t = (1:numel(h.sp))/fr; w = t >= 200 & t <= 320;
        nexttile; hold on
        yyaxis left
        yl = [min(h.F(w)) max(h.F(w))]; yl = yl + [-.05 .05]*range(yl);
        fl = double(h.flag(w)); fl(fl==0) = nan;
        area(t(w), fl*yl(2), yl(1), 'FaceColor', 'r', 'FaceAlpha', .25, 'EdgeColor', 'none')
        plot(t(w), h.F(w), 'k-'); plot(t(w), h.p(w), '-', 'Color', [0 .6 0], 'LineWidth', 1.5); ylim(yl); ylabel('PB-mean dF/F')
        yyaxis right; plot(t(w), h.sp(w), 'm-'); ylabel('|vr| (rad/s)')
        title(sprintf('%s, trial %d (flagged %.2f of moving time)', names{g}, o(k), v(o(k))))
    end
end
legend({'flagged low','measured','speed prediction (CV)','|vr|'},'Location','northeastoutside'); xlabel('time (s)')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_vs_lpsp_speed_dropout_examples.png'),'Resolution',150);

%% functions
function [pred, info] = speed_model(F, sp, bad, fr, opts)
% cross-validated LN prediction of F (1xT) from |vr| sp (1xT): ridge linear
% filter on sp at lags opts.lagmin..lagmax s, then a binned static
% nonlinearity; 2-fold CV over alternating opts.block_s blocks
lags = round(opts.lagmin*fr):round(opts.lagmax*fr);
T = numel(F); X = nan(T,numel(lags));
for j = 1:numel(lags)
    L = lags(j); x = nan(1,T);
    if L >= 0, x(1+L:T) = sp(1:T-L); else, x(1:T+L) = sp(1-L:T); end
    X(:,j) = x;
end
ok = ~bad(:) & ~isnan(F(:)) & all(~isnan(X),2);
blk = mod(floor((0:T-1)'/round(opts.block_s*fr)),2);
pred = nan(1,T); Xc = [X ones(T,1)];
for f = 0:1
    tr = ok & blk~=f; te = blk==f & all(~isnan(X),2);
    A = Xc(tr,:); y = F(tr)';
    Rg = opts.lambda*eye(size(A,2)); Rg(end,end) = 0;
    w = (A'*A + Rg) \ (A'*y);
    lin_tr = A*w; lin_te = Xc(te,:)*w;
    e = unique(quantile(lin_tr, linspace(0,1,opts.nbins+1))); e(1) = -inf; e(end) = inf;
    bt = discretize(lin_tr,e);
    mb = accumarray(bt, y, [numel(e)-1 1], @mean, nan); cb = accumarray(bt, lin_tr, [numel(e)-1 1], @mean, nan);
    okb = ~isnan(mb) & ~isnan(cb);
    pred(te) = interp1(cb(okb), mb(okb), lin_te, 'linear', 'extrap');
end
r = F - pred;
info.R2 = 1 - var(r(ok),'omitnan')/var(F(ok),'omitnan');
end
