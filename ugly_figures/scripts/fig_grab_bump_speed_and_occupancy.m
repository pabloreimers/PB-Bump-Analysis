%% fig_grab_bump_speed_and_occupancy
% Two additions to the EPG GRAB(DA) model GRAB ~ speed gain x (1 + k EPG
% bump) (fig_summary_lpsp_grab_bump_model, panel I), which on its own has no
% term for how the bump moves:
%  A/B  bump speed predicts the GRAB peak beyond fly speed
%  C/D  recent occupancy of the bump's location lowers the GRAB peak
%       (synaptic-depletion-like: release is lower where the bump just was)
% Data: EPG>GRAB(DA2m) trials of lpsp_cl_reproc_nobg_20261117.mat (21 closed
% loop, 21 dark), PMT-offset-corrected per-bin dF/F; PB halves folded after
% aligning by each trial's left-right offset (as in
% fig_grab_peak_fly_vs_bump_speed, sections 4-6 and 15, where the
% derivations and more controls are).
%
% A  per trial, frame-wise joint regression (standardized) of the GRAB peak
%    (PB max) on |fly speed| (0.45 s earlier) and |bump speed| (folded GRAB
%    bump phase, 1 s smoothing, best lag in -0.5..+1 s), frames with a clear
%    bump only (rho > 0.4 over the 1 s velocity window). Closed-loop trials
%    whose bump follows the cue (|r| >= 0.4) and all dark trials. Also with
%    net rotation and forward speed as extra controls (printed).
% B  matched pairs within trials: non-overlapping 1.5 s stretches with the
%    same NET rotation (within 15%, > 1 rad/s, consistent turning >= 0.8),
%    the same forward speed (within 0.15) and bump speed >= 1.5x different;
%    per trial, % of pairs where the faster-bump stretch has the brighter
%    GRAB peak (trials with >= 5 pairs).
% C  closed loop, bump position from the cue (cue-predicted position, time
%    lag taken from EPG GCaMP7f trials of the same file, so it is known when
%    GRAB is dark and is independent of GRAB's own phase). GRAB peak
%    regressed on cubic terms of |fly speed|, net rotation, forward speed and
%    bump speed, mean bump AND fly speed over the past 1/2/4 s (fly: also
%    squared), reversal covariates,
%    and the recent occupancy of the bump's current location (exponentially
%    weighted closeness of past positions, tau 1-8 s); first 30 s excluded.
%    Placebo: occupancy of the location opposite the bump.
% D  GRAB residual after the same controls (no occupancy), bump moving > 0.5
%    rad/s, by time since the bump last overlapped its current location.
%
% Result when written: A bump speed +0.16 (CL) / +0.19 (dark) next to fly
% speed 0.52 / 0.45; B faster bump brighter in 66% / 69% of pairs; C
% occupancy -0.07 (tau 1 s, 12/13 trials, p 0.006), -0.03 (2 s), ~0 beyond;
% opposite location +0.03 (n.s.); D revisit within 1-2 s -0.18 SD, fresh
% (> 8 s) minus revisit (< 2 s) +0.24 SD (9/10, p 0.01). (An earlier version
% without the fly-speed-history controls gave C about twice as large: part
% of it was the gain's integration of recent fly speed.)

%% load + per-trial signals
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_reproc_nobg_20261117.mat')), data_dir = '.data'; else, data_dir = 'data'; end
B = load(fullfile(data_dir,'lpsp_cl_reproc_nobg_20261117.mat'),'all_data'); B = B.all_data(:);
ind = arrayfun(@(s) s.info.indicator, B, 'UniformOutput', false);
gi = find(strcmp(ind,'EPG>GRAB(DA2m)'))'; ei = find(strcmp(ind,'EPG>GCaMP7f'))';
ei = ei(arrayfun(@(k) prctile(double(B(k).flash.meanImg(:)),0.1) < median(prctile(B(k).im.f(:,~logical(B(k).info.badVol(:)')),7,2)), ei));   % drop meanImg/im.f scale mismatch
snr_frac = 0.5; a16 = linspace(-pi, pi-2*pi/16, 16)'; lagsS = -3:0.1:3; rthr = 0.4; edges = [0 0.25 0.5 1 2 inf];
Fd = arrayfun(@(k) grab_folded(B(k), snr_frac, a16), gi); Fe = arrayfun(@(k) grab_folded(B(k), snr_frac, a16), ei);
for j = 1:numel(Fd)
    s = B(gi(j)); Fd(j).fw = interp1(s.ft.xf(:), smoothdata(s.ft.f_speed(:),'gaussian',30), s.ft.xb(:)', 'linear', 'extrap');
    Fd(j).net = abs(smoothdata(Fd(j).vr,'gaussian',max(3,round(Fd(j).fr))));
end
RkF = nan(1,numel(Fd)); for j = find(~[Fd.dark]), RkF(j) = cue_track(Fd(j)); end
QE = []; for j = find(~[Fe.dark]), [r_, l_] = cue_track(Fe(j)); QE(end+1,:) = [r_ l_]; end %#ok<AGROW>
lagE = median(QE(QE(:,1) >= 0.4, 2));
fprintf('GRAB: %d CL (%d follow the cue), %d dark | cue -> bump lag (EPG7f) %.2f s\n', sum(~[Fd.dark]), sum(RkF >= 0.4), sum([Fd.dark]), lagE);

%% (A) bump speed beyond fly speed: joint model per trial
SA = struct();
for d_ = 0:1
    trials = find([Fd.dark] == d_ & (d_ == 1 | RkF >= 0.4));
    [C, lg] = joint_coef(Fd, trials, 'pk', 0.45, rthr, {}, lagsS);
    C4 = joint_coef(Fd, trials, 'pk', 0.45, rthr, {'fly:net','fly:fw'}, lagsS);
    ok = ~isnan(C(:,2)); SA(d_+1).C = C(ok,1:2); SA(d_+1).lag = lg; SA(d_+1).p = signrank(C(ok,2));
    fprintf('A %-11s (%2d trials, bump lag %+.1f s): fly %.2f, bump %.2f (>0 %d/%d, p=%.2g) | + net rotation & forward speed: bump %.2f (>0 %d/%d)\n', ifelse_c(d_), sum(ok), lg, ...
        median(C(ok,1)), median(C(ok,2)), sum(C(ok,2)>0), sum(ok), SA(d_+1).p, median(C4(:,2),'omitnan'), sum(C4(:,2)>0), sum(~isnan(C4(:,2))));
end

%% (B) matched pairs: same net rotation and forward speed, different bump speed
wsP = 1.5; P2 = [];
for j = 1:numel(Fd)
    t = Fd(j); if ~t.dark && ~(RkF(j) >= 0.4), continue, end
    fr = t.fr; n = numel(t.ph); L = round(wsP*fr); st = round(0.5*fr); Lf = round(0.45*fr); Lb = round(0.2*fr);
    bs = abs(gradient(smoothdata(t.ph,'gaussian',max(3,round(fr))))*fr); good = t.rho_raw > rthr & ~t.bad; rows = [];
    for a = 1+max(Lf,Lb):st:n-L
        g_ = a:a+L-1; fi = g_ - Lf; bi = g_ - Lb;
        if mean(good(bi)) < 0.8 || any(t.bad(g_)) || any(isnan(t.pk(g_))), continue, end
        rows(end+1,:) = [a abs(mean(t.vr(fi))) mean(t.fw(fi)) mean(bs(bi(good(bi)))) mean(t.pk(g_)) abs(mean(t.vr(fi)))/mean(abs(t.vr(fi)))]; %#ok<AGROW>
    end
    for p1 = 1:size(rows,1)-1
        for p2 = p1+1:size(rows,1)
            if abs(rows(p1,1)-rows(p2,1)) < L, continue, end
            nA = rows(p1,2); nB = rows(p2,2); if min(nA,nB) < 1 || abs(nA-nB)/mean([nA nB]) > 0.15, continue, end
            if min(rows([p1 p2],6)) < 0.8 || abs(rows(p1,3)-rows(p2,3)) > 0.15, continue, end
            if max(rows(p1,4),rows(p2,4))/max(min(rows(p1,4),rows(p2,4)),eps) < 1.5, continue, end
            if rows(p1,4) > rows(p2,4), q = [p2 p1]; else, q = [p1 p2]; end
            P2(end+1,:) = [j t.dark rows(q(1),5) rows(q(2),5)]; %#ok<AGROW>
        end
    end
end
SB = struct();
for d_ = 0:1
    sel = P2(:,2) == d_; tr = unique(P2(sel,1)); npt = arrayfun(@(k) sum(sel & P2(:,1)==k), tr); keep = npt >= 5;
    pt = arrayfun(@(k) mean(P2(sel & P2(:,1)==k,4) > P2(sel & P2(:,1)==k,3)), tr);
    ptr = arrayfun(@(k) median((P2(sel & P2(:,1)==k,4) - P2(sel & P2(:,1)==k,3))./P2(sel & P2(:,1)==k,3)), tr);
    SB(d_+1).pt = pt(keep); SB(d_+1).p = signrank(pt(keep) - 0.5); SB(d_+1).npairs = sum(npt(keep));
    fprintf('B %-11s %4d pairs, %2d trials: faster bump brighter in per-trial median %.0f%% (p=%.2g), GRAB difference %+.0f%%\n', ifelse_c(d_), SB(d_+1).npairs, sum(keep), ...
        100*median(pt(keep)), SB(d_+1).p, 100*median(ptr(keep)));
end

%% (C, D) recent occupancy of the bump's location (closed loop)
taus = [1 2 4 8]; trials = find(~[Fd.dark] & RkF >= 0.4); cap = 30; tb = [0 1 2 4 8 16 30];
pb = @(z) [z(:) z(:).^2 z(:).^3]; kk = log(2)/(1 - cos(pi/4));
CS = nan(numel(trials), numel(taus)); CPo = CS; RC = nan(numel(trials), numel(tb));
for jj = 1:numel(trials)
    j = trials(jj); t = Fd(j); fr = t.fr; n = numel(t.ph); y = t.pk; Lf = round(0.45*fr); Lb = round(0.3*fr); sh = @(z,L) [nan(1,L) z(1:n-L)];
    [~,~,~,ref,v] = cue_ref_profiles(B(gi(j)), t, snr_frac, lagE, edges, zeros(16,n)); H = bump_history(ref, fr, taus); refi = fillmissing(ref,'nearest');
    sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); fw = smoothdata(t.fw,'gaussian',max(3,round(0.5*fr))); bsp = abs(v);
    vs = smoothdata(v,'gaussian',max(3,round(0.3*fr))); cur = sign(vs); opp = nan(1,n); nrev = nan(1,n); W = round(2*fr);
    for i = W+1:n, seg = vs(i-W:i-1); opp(i) = mean(sign(seg) == -cur(i) & abs(seg) > 0.3); nrev(i) = sum(abs(diff(sign(seg(abs(seg) > 0.3)))) > 0); end
    Cs = [pb(sh(sp,Lf)) pb(sh(t.net,Lf)) pb(sh(fw,Lf)) pb(sh(bsp,Lb)) sh(opp,Lb)' sh(nrev,Lb)'];
    for wsec = [1 2 4], Cs = [Cs movmean(sh(bsp,Lb), [round(wsec*fr) 0])']; end %#ok<AGROW>
    for wsec = [1 2 4], h_ = movmean(sh(sp,Lf), [round(wsec*fr) 0])'; Cs = [Cs h_ h_.^2]; end %#ok<AGROW>   % fly-speed history (and its square)
    early = (1:n)' <= round(30*fr);
    for it = 1:numel(taus)
        K = round(5*taus(it)*fr); w = exp(-(1:K)/(taus(it)*fr)); w = w/sum(w); oO = zeros(1,n);
        for q = 1:K, past = [nan(1,q) refi(1:n-q)]; cO = exp(kk*(cos(past - refi - pi) - 1)); cO(isnan(cO)) = 0; oO = oO + w(q)*cO; end
        oO(1:K) = nan;
        X = [Cs sh(H.occ(it,:),Lb)' sh(oO,Lb)']; m = all(isfinite([X y']),2) & ~t.bad' & ~early;
        Z = zscore([X(m,:) y(m)']); b = [ones(sum(m),1) Z(:,1:end-1)] \ Z(:,end); CS(jj,it) = b(end-1); CPo(jj,it) = b(end);
    end
    m = all(isfinite([Cs y']),2) & ~t.bad' & ~early; b = [ones(sum(m),1) Cs(m,:)] \ y(m)'; r = nan(1,n); r(m) = y(m) - ([ones(sum(m),1) Cs(m,:)]*b)'; r = r/std(r(m));
    tsl = sh(H.tsl,Lb); tsl(isnan(tsl) & ~isnan(r)) = cap + 1; g = ~isnan(r) & sh(bsp,Lb) > 0.5;
    for k = 1:numel(tb)
        if k < numel(tb), mm = g & tsl >= tb(k) & tsl < tb(k+1); else, mm = g & tsl > cap; end
        if sum(mm) >= 20, RC(jj,k) = mean(r(mm)); end
    end
end
fprintf('C occupancy of the bump''s location (tau 1/2/4/8 s): %s, <0 in %s of %d, p %s\n', mat2str(round(median(CS,1),3)), mat2str(sum(CS<0,1)), numel(trials), mat2str(arrayfun(@(k) signrank(CS(:,k)), 1:4),2));
fprintf('  occupancy of the OPPOSITE location:          %s, <0 in %s, p %s\n', mat2str(round(median(CPo,1),3)), mat2str(sum(CPo<0,1)), mat2str(arrayfun(@(k) signrank(CPo(:,k)), 1:4),2));
fprintf('D GRAB residual by time since last visit [1-2 2-4 4-8 8-16 16-30 >30 s]: %s\n', mat2str(round(median(RC(:,2:end),1,'omitnan'),3)));
dFR = mean(RC(:,5:7),2,'omitnan') - mean(RC(:,1:2),2,'omitnan');
fprintf('  fresh (> 8 s) minus revisit (< 2 s): %+.3f SD (>0 %d/%d, p=%.2g)\n', median(dFR,'omitnan'), sum(dFR>0), sum(~isnan(dFR)), signrank(dFR));

%% figure
figure('Color','w','Position',[60 60 1250 820]); tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact'); rng(1);
cF = [0.15 0.15 0.15]; cB = [0 0.45 0.74]; cO = [0.85 0.4 0.1];
% A
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]);
for d_ = 0:1
    C = SA(d_+1).C; xs = d_*3 + [1 2];
    plot(ax, repmat(xs',1,size(C,1)), C', '-', 'Color', [0.88 0.88 0.88]);
    for e = 1:2
        cc = {cF, cB}; scatter(ax, xs(e)*ones(size(C,1),1), C(:,e), 18, cc{e}, 'filled', 'MarkerFaceAlpha', 0.45);
        plot(ax, xs(e) + [-0.3 0.3], median(C(:,e))*[1 1], '-', 'Color', cc{e}, 'LineWidth', 3);
    end
    text(ax, xs(2), 0.95, sprintf('bump %+.2f\n%s', median(C(:,2)), ifelse_p(SA(d_+1).p)), 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color',cB);
end
set(ax,'XTick',[1 2 4 5],'XTickLabel',{'fly speed','bump speed','fly speed','bump speed'},'Box','off','TickDir','out'); xlim(ax,[0.3 5.7]); ylim(ax,[-0.4 1]);
text(ax, 1.5, -0.35, 'closed loop', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]); text(ax, 4.5, -0.35, 'dark', 'HorizontalAlignment','center','Color',[0.4 0.4 0.4]);
ylabel(ax,'joint coefficient on GRAB peak (z)'); title(ax,'A  bump speed adds to fly speed (lines = trials)','FontWeight','normal');
% B
ax = nexttile(tl); hold(ax,'on'); yline(ax, 50, 'Color', [0.6 0.6 0.6]);
for d_ = 0:1
    v_ = 100*SB(d_+1).pt; scatter(ax, d_ + 1 + 0.08*randn(size(v_)), v_, 22, cB, 'filled', 'MarkerFaceAlpha', 0.45);
    plot(ax, d_ + 1 + [-0.25 0.25], median(v_)*[1 1], '-', 'Color', cB, 'LineWidth', 3);
    text(ax, d_ + 1, 114, sprintf('%.0f%%, %s\n%d trials, %d pairs', median(v_), ifelse_p(SB(d_+1).p), numel(v_), SB(d_+1).npairs), 'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8,'Color',cB);
end
set(ax,'XTick',[1 2],'XTickLabel',{'closed loop','dark'},'Box','off','TickDir','out'); xlim(ax,[0.4 2.6]); ylim(ax,[0 115]); yticks(ax, 0:20:100);
ylabel(ax,'% pairs where the faster-bump stretch is brighter');
title(ax,sprintf('B  same net rotation & forward speed, different bump speed\n(1.5 s stretches, bump speed >= 1.5x apart)'),'FontWeight','normal');
% C
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6], 'HandleVisibility','off');
for it = 1:numel(taus)
    for c = 1:2
        V = {CS, CPo}; cc = {cB, cO}; x0 = it + (c - 1.5)*0.3; v_ = V{c}(:,it);
        scatter(ax, x0 + 0.04*randn(size(v_)), v_, 16, cc{c}, 'filled', 'MarkerFaceAlpha', 0.4, 'HandleVisibility','off');
        h_ = plot(ax, x0 + [-0.1 0.1], median(v_)*[1 1], '-', 'Color', cc{c}, 'LineWidth', 3); if it > 1, h_.HandleVisibility = 'off'; end
    end
end
set(ax,'XTick',1:numel(taus),'XTickLabel',arrayfun(@(x) sprintf('%g s',x), taus, 'UniformOutput', false),'Box','off','TickDir','out'); xlim(ax,[0.4 numel(taus)+0.6]);
xlabel(ax,'occupancy window \tau'); ylabel(ax,'partial coefficient on GRAB peak (z)');
legend(ax, {'bump''s own location','opposite location (placebo)'}, 'Box','off','Location','southeast','FontSize',8);
title(ax,sprintf('C  recent occupancy of the bump''s location lowers GRAB\n(closed loop; speeds, speed history, reversals controlled)'),'FontWeight','normal');
% D
ax = nexttile(tl); hold(ax,'on'); yline(ax, 0, 'Color', [0.6 0.6 0.6]); xc = [1.5 3 6 12 23 40]; RCd = RC(:,2:end);
for jj = 1:size(RCd,1), plot(ax, xc, RCd(jj,:), '-', 'Color', [0.85 0.85 0.85]); end
errorbar(ax, xc, median(RCd,1,'omitnan'), std(RCd,0,1,'omitnan')./sqrt(sum(~isnan(RCd),1)), 'o-', 'Color', cB, 'LineWidth', 2, 'MarkerFaceColor', 'w');
set(ax,'XScale','log','XTick',xc,'XTickLabel',{'1-2','2-4','4-8','8-16','16-30','>30'},'Box','off','TickDir','out'); xlim(ax,[1 50]); ylim(ax,[-0.6 0.7]);
xlabel(ax,'time since the bump last visited its current location (s)'); ylabel(ax,'GRAB peak residual (SD)');
text(ax, 0.03, 0.97, sprintf('fresh (> 8 s) - revisit (< 2 s): %+.2f SD, %s', median(dFR,'omitnan'), ifelse_p(signrank(dFR))), 'Units','normalized','VerticalAlignment','top','FontSize',8,'Color',cB);
title(ax,'D  recovery after a visit (moving bump, closed loop)','FontWeight','normal');
title(tl, 'EPG GRAB(DA): beyond the speed x bump product, bump speed adds and recently visited glomeruli release less', 'FontSize', 12);
if ~exist(fullfile('ugly_figures','exports'),'dir'), mkdir(fullfile('ugly_figures','exports')); end
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_bump_speed_and_occupancy.png'),'Resolution',150);

%% functions
function T = grab_folded(s, snr_frac, a16)
% folded GRAB bump: halves aligned by the trial's mean left-right phase
% offset and averaged -> 16-bin profile; phase + rho (|mean(z e^{i alpha})|)
xb = s.ft.xb(:)'; xf = s.ft.xf(:); fr = s.info.fps; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z0 = Z; Z0(isnan(Z0)) = 0;
vL = exp(1i*a16).'*Z0(1:16,:); vR = exp(1i*a16).'*Z0(17:32,:);
dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
Zf = mean(cat(3, Z(1:16,:), ring_shift(Z(17:32,:), dLR)), 3, 'omitnan'); Zf(isnan(Zf)) = 0;
v = exp(1i*a16).'*Zf / 16;
T.fr = fr; T.dark = s.info.dark; T.bad = bad;
T.ph = smoothdata(unwrap(angle(v)), 'gaussian', max(3, round(0.5*fr))); T.rho_raw = abs(v);
T.vr = interp1(xf, smoothdata(s.ft.r_speed(:),'gaussian',30), xb, 'linear', 'extrap');
T.pk = max(smoothdata(D,1,'movmean',3,'omitnan'),[],1);
T.cue = unwrap(angle(interp1(xf, exp(1i*s.ft.cue(:)), xb, 'linear', 'extrap')));
end

function Y = ring_shift(X, d)
% Y(alpha) = X(alpha - d) on the 16-bin ring, linear interpolation
nb = size(X,1); step = 2*pi/nb; a = linspace(-pi, pi-step, nb)';
pos = mod((a - d + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
Y = (1-w).*X(lo+1,:) + w.*X(hi+1,:);
end

function [far, nfar, anti, Df, ftr] = far_field(s, t, snr_frac, Tw, minDeg)
% far-field PB signal: mean folded dF/F over bins that stayed > minDeg from
% the bump phase throughout the last Tw s (so not touched by the bump or
% its decaying trail). nfar = # such bins; anti = mean of bins > minDeg from
% the CURRENT phase only (no trail exclusion). ftr = far-field TROUGH: min
% over those untouched bins of the 3-bin (circular) smoothed profile.
% t = grab_folded output.
a16 = linspace(-pi, pi-2*pi/16, 16)'; bad = logical(s.info.badVol(:)');
off = prctile(double(s.flash.meanImg(:)), 0.1);
bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
V = s.im.f - off - (bgs - median(bgs(~bad)));
x = V(:,~bad); x = x - mean(x,2);
ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
snr = ac ./ max(var(x,0,2) - ac, eps); V(snr < snr_frac*median(snr),:) = nan;
F0 = prctile(V(:,~bad),7,2); D = (V - F0)./F0;
Z = (V - mean(V(:,~bad),2,'omitnan'))./std(V(:,~bad),0,2,'omitnan'); Z(isnan(Z)) = 0;
vL = exp(1i*a16).'*Z(1:16,:); vR = exp(1i*a16).'*Z(17:32,:); dLR = angle(mean(vL(~bad).*conj(vR(~bad))));
Df = mean(cat(3, D(1:16,:), ring_shift(D(17:32,:), dLR)), 3, 'omitnan');
ph = angle(exp(1i*t.ph)); n = numel(ph); K = round(Tw*t.fr); th = minDeg*pi/180;
dmin = inf(16, n);
for k = 0:K
    p = [nan(1,k) ph(1:n-k)];                       % phase k frames ago
    d = abs(angle(exp(1i*(a16 - p))));            % 16 x n distance
    d(isnan(d)) = inf; dmin = min(dmin, d);
end
dmin(:, 1:K) = nan;                                % trail not defined at start
M = dmin > th; Dm = Df; Dm(~M) = nan;
far = mean(Dm, 1, 'omitnan'); nfar = sum(M & ~isnan(Df), 1); far(nfar < 2) = nan;
Ds = (circshift(Df,1,1) + Df + circshift(Df,-1,1))/3; Ds(~M) = nan; ftr = min(Ds, [], 1, 'omitnan'); ftr(nfar < 2) = nan;
dc = abs(angle(exp(1i*(a16 - ph)))); Da = Df; Da(dc <= th) = nan; anti = mean(Da, 1, 'omitnan');
end

function [C, lag, N] = joint_coef(FF, trials, yname, lf, th, xnames, lagsS, bwin)
% per-trial standardized joint regression of FF(j).(yname) on |fly speed|
% (lag lf s) + |folded bump speed| (scanned lag) + optional extra regressors
% xnames (fields of FF, same-time as outcome unless named 'fly:<field>' =
% taken at the fly lag). Only clear-bump frames (rho_raw > th over the 1 s
% velocity window). Bump lag = max |frame-weighted mean coefficient| in bwin.
% C: trials x (2 + numel(xnames)) at that lag, order [fly bump extras].
if nargin < 8, bwin = [-0.5 1]; end
nx = numel(xnames); P = 2 + nx;
Ca = nan(numel(trials), numel(lagsS), P); NF = nan(numel(trials), numel(lagsS));
for jj = 1:numel(trials)
    t = FF(trials(jj)); fr = t.fr; n = numel(t.ph); w = max(3,round(fr)); Lf = round(lf*fr); y = t.(yname);
    sp = smoothdata(abs(t.vr),'gaussian',max(3,round(0.5*fr))); bs = abs(gradient(smoothdata(t.ph,'gaussian',w))*fr);
    good = movmin(t.rho_raw,w) > th & ~movmax(t.bad,w);
    for li = 1:numel(lagsS)
        L = round(lagsS(li)*fr); i1 = max([1,1+L,1+Lf]):min([n,n+L,n+Lf]); i0 = i1 - L; iF = i1 - Lf;
        X = [sp(iF)' bs(i0)'];
        for q = 1:nx
            nm = xnames{q};
            if startsWith(nm,'fly:'), v = t.(nm(5:end)); X(:,end+1) = v(iF)'; else, v = t.(nm); X(:,end+1) = v(i1)'; end %#ok<AGROW>
        end
        m = good(i0) & ~t.bad(i1) & isfinite(y(i1)) & all(isfinite(X),2)'; if sum(m) < 50, continue, end
        Zq = zscore([X(m,:) y(i1(m))']); b = [ones(sum(m),1) Zq(:,1:P)] \ Zq(:,P+1);
        Ca(jj,li,:) = b(2:end); NF(jj,li) = sum(m);
    end
end
W = sum(Ca(:,:,2).*NF,2-1,'omitnan')./sum(NF.*~isnan(Ca(:,:,2)),1,'omitnan');
win = find(lagsS >= bwin(1) & lagsS <= bwin(2)); [~,k] = max(abs(W(win))); jb = win(k);
C = squeeze(Ca(:,jb,:)); if numel(trials) == 1, C = C(:)'; end
lag = lagsS(jb); N = NF(:,jb);
end

function [r, lag, sg, off] = cue_track(t, lagFix)
% How well the bump phase follows the visual cue: bump(t) ~ off - sg*cue(t+lag).
% Scans lag in -2.5..+1 s and sign (or uses lagFix s if given); returns the
% tracking strength r = |mean exp(i(bump + sg*cue(t+lag)))|, the lag, the
% sign and the constant offset off = angle of that mean.
n = numel(t.ph); ph = angle(exp(1i*t.ph)); cue = angle(exp(1i*t.cue)); r = 0; lag = nan; sg = nan; off = nan;
if nargin > 1 && ~isempty(lagFix), tts = round(lagFix*t.fr); else, tts = -round(2.5*t.fr):round(1*t.fr); end
for tt = tts
    idx = 1:n; idx = idx(idx+tt >= 1 & idx+tt <= n & ~t.bad(idx));
    for s_ = [-1 1]
        m = mean(exp(1i*(ph(idx) + s_*cue(idx+tt))));
        if abs(m) > r, r = abs(m); lag = tt/t.fr; sg = s_; off = angle(m); end
    end
end
end

function [P, N, cen, ref, v] = cue_ref_profiles(s, t, snr_frac, lagRef, edges, Df)
% Closed loop: folded dF/F profile aligned on the CUE-predicted bump
% position ref(t) = off - sg*cue(t + lagRef) (lagRef fixed, e.g. from EPG
% 7f; sign and constant offset fitted per trial), mirrored so the bump's
% direction of travel is +, binned by |ref speed| (rad/s). P: 16 x nbins
% (a16 = bins relative to ref, + = ahead), N: frames per bin, cen: circular
% centroid (deg) of P - min(P) per bin. Optional Df (16 x n) replaces the
% measured profile (e.g. a model), aligned with the same reference.
a16 = linspace(-pi, pi-2*pi/16, 16)';
if nargin < 6 || isempty(Df), [~,~,~,Df] = far_field(s, t, snr_frac, 2, 90); end
[~,~,sg,off] = cue_track(t, lagRef);
n = numel(t.ph); tt = round(lagRef*t.fr); cu = unwrap(t.cue(:)');
idx = (1:n) + tt; ok = idx >= 1 & idx <= n; ref = nan(1,n); ref(ok) = off - sg*cu(idx(ok));
v = gradient(smoothdata(ref,'gaussian',max(3,round(0.5*t.fr)),'omitnan'))*t.fr;
Y = shift_cols(Df, -ref);                                   % ref at bin 0
mir = mod(17 - (1:16), 16) + 1; Ym = Y(mir,:);              % theta -> -theta
Y(:, v < 0) = Ym(:, v < 0);                                 % direction of travel = +
[~,~,b] = histcounts(abs(v), edges); nb = numel(edges) - 1; P = nan(16,nb); N = zeros(1,nb); cen = nan(1,nb);
for k = 1:nb
    m = b == k & ~t.bad & ok & ~isnan(v); N(k) = sum(m); if N(k) < 30, continue, end
    P(:,k) = mean(Y(:,m), 2, 'omitnan');
    q = P(:,k) - min(P(:,k)); cen(k) = angle(sum(q.*exp(1i*a16)))*180/pi;
end
end

function Y = shift_cols(X, d)
% per-column fractional circular shift: Y(alpha,t) = X(alpha - d(t), t); NaN d -> NaN column
nb = size(X,1); step = 2*pi/nb; a = linspace(-pi, pi-step, nb)'; dd = d(:)'; bad = isnan(dd); dd(bad) = 0;
pos = mod((a - dd + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
cols = repmat(1:size(X,2), nb, 1);
Y = (1-w).*X(sub2ind(size(X), lo+1, cols)) + w.*X(sub2ind(size(X), hi+1, cols)); Y(:,bad) = nan;
end

function H = bump_history(ref, fr, taus)
% Bump-history ("depletion") variables from a bump-position trace ref (rad,
% e.g. the cue-predicted position in closed loop; may contain NaN).
%  H.occ(i,:)  recent occupancy of the CURRENT location over the past
%              taus(i) s: exponentially weighted mean over past frames of
%              closeness c(dphi) = exp(k (cos dphi - 1)), c(45 deg) = 0.5
%              (excludes the current frame); high = bump has been here.
%  H.tsl       time since the bump last overlapped (|dphi| < 22.5 deg) its
%              current location BEFORE the current stay (the bump must have
%              been > 45 deg away in between); NaN if not within 30 s.
%  H.dwell     time the bump has stayed within 22.5 deg of its current spot.
n = numel(ref); k = log(2)/(1 - cos(pi/4)); refi = fillmissing(ref(:)','nearest');
H.occ = nan(numel(taus), n);
for it = 1:numel(taus)
    K = round(5*taus(it)*fr); w = exp(-(1:K)/(taus(it)*fr)); w = w/sum(w); acc = zeros(1,n); ws = zeros(1,n);
    for q = 1:K
        past = [nan(1,q) refi(1:n-q)]; c = exp(k*(cos(past - refi) - 1)); g = ~isnan(c);
        acc(g) = acc(g) + w(q)*c(g); ws(g) = ws(g) + w(q);
    end
    o = acc./ws; o(ws < 0.9) = nan; H.occ(it,:) = o;
end
H.tsl = nan(1,n); H.dwell = nan(1,n); Lmax = round(30*fr);
for i = 2:n
    if isnan(ref(i)), continue, end
    lo = max(1, i - Lmax); d = abs(angle(exp(1i*(refi(i-1:-1:lo) - refi(i)))));   % back in time
    away = find(d > pi/4, 1);                                                     % end of current stay
    inn = find(d >= pi/8, 1); if isempty(inn), H.dwell(i) = (i-lo)/fr; else, H.dwell(i) = (inn-1)/fr; end
    if isempty(away), continue, end
    prev = find(d(away:end) < pi/8, 1); if ~isempty(prev), H.tsl(i) = (away + prev - 1)/fr; end
end
end

function y = ifelse_c(d), if d, y = 'dark'; else, y = 'closed loop'; end, end

function y = ifelse_p(p), if p < 0.001, y = 'p<0.001'; else, y = sprintf('p=%.2g', p); end, end
