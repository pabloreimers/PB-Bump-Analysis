%% fig_grab_bump_arrival_dynamics
% Does DA release (EPG>GRAB) favor "new" glomeruli, i.e. is local release
% transient when the bump parks on a glomerulus, while LPsP presynaptic
% calcium (LPsP>syt7f) at that glomerulus is sustained? And how does that
% compare with the EPG bump itself (EPG>GCaMP7f)?
%
% Datasets (all 9.1 Hz volumes, 32 PB bins, alpha identical):
%   LPsP>syt7f, EPG>GRAB(DA2m)  data/lpsp_cl_data_20240206.mat
%       NOTE: folders named '*_EPG_syt7f_*' are a naming typo -- those are
%       GRAB(DA2m) trials, and are counted as GRAB here.
%   EPG>GCaMP7f                 data/epg_7f_20261007.mat (Figure_1_claude.m's
%       column-1 'epg7f' source; closed-loop gain 0.7, 300 s trials). 7f
%       roughly matches syt7f kinetics, so it is the reference for how the
%       EPG bump itself behaves at arrivals.
%
% Bump position is predicted from the visual cue in closed-loop trials: for
% each trial we find the sign, offset and lag (-2.5..+1 s) that best lock the
% signal's own PVA phase to the cue, then use cue(t - lag) + offset as the
% bump position. This keeps "where the bump is" independent of the local
% amplitude being measured. Trials that don't lock (R < 0.4) are dropped.
% (lpsp_cl signals all sit ~0.8-1 s behind the cue -- a shared
% imaging/fictrac alignment offset in that dataset's processing -- which
% the per-trial lag fit absorbs.)
%
% An "arrival" = bump-predicted glomerulus changes to j, the bump was >=60deg
% away 1.5 s earlier, and it stays within +/-1.5 bins of j for >= 4 s.
% Signal = (glomerulus dF/F - PB-mean dF/F), L/R homologs averaged,
% z-scored per trial, so the global speed-driven component is removed.
%
% ==CAVEAT -- READ BEFORE INTERPRETING FIGURE 2==
% Arrivals happen at the END of turns (|vr| peaks ~1 s before the bump
% parks, then drops), so the global, speed-driven level falls across every
% arrival. "Local minus PB mean" therefore looks transient for any signal
% whose spatial modulation SCALES with the global level, even if its
% spatial shape never changes. Figure 3 (decomposition) and figure 4
% (bump excess vs. global level over ALL closed-loop frames) show that is
% what happens for GRAB: its bump contrast is constant and its excess is
% proportional to the PB level (multiplicative), whereas LPsP calcium has a
% speed-independent heading offset on top of the global signal (additive).
% So the GRAB "transience" in figure 2 is a gain effect, not evidence that
% release is depressed at glomeruli the bump has parked on.

%% load
fr = 9.1;
d = load(fullfile('data','lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isDark = contains(meta,'dark');
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f');
isLpsp = contains(meta,'LPsP_syt7f');

E = load(fullfile('data','epg_7f_20261007.mat'),'all_data'); E = E.all_data(:);
E = E(~arrayfun(@(x) x.ft.dark, E));

% unify to {D (32 x T dF/F), cue (1 x T, imaging clock), bad (1 x T)}
T = struct('D',{},'cue',{},'bad',{},'sp',{});
U = {T,T,T};
for i = find(isLpsp & ~isDark), U{1}(end+1) = from_lpsp_cl(d(i)); end
for i = find(isGrab & ~isDark), U{2}(end+1) = from_lpsp_cl(d(i)); end
for i = 1:numel(E),             U{3}(end+1) = from_epg7f(E(i));   end
names = {'LPsP > syt7f','EPG > GRAB(DA2m)','EPG > GCaMP7f'};
alpha = d(1).im.alpha;

%% collect arrival-triggered local signals
pre = round(3*fr); post = round(8*fr); minstay = round(4*fr);
ETA = cell(1,3); ETAo = ETA; nTr = zeros(1,3); lagFit = cell(1,3);
for g = 1:3
    for k = 1:numel(U{g})
        [e,eo,ok,lg] = arrival_eta(U{g}(k), alpha, fr, pre, post, minstay);
        if ~ok, continue, end
        nTr(g) = nTr(g)+1; lagFit{g}(end+1) = lg;
        ETA{g} = [ETA{g}; e]; ETAo{g} = [ETAo{g}; eo];
    end
end

%% quantify: sustained (2.5-4 s after arrival) / transient peak, bootstrap CI over arrivals
tt = (-pre:post)/fr;
w_pk = pre + (-round(1.5*fr):round(1*fr)); w_late = pre + (round(2.5*fr):round(4*fr)); w_base = 1:round(1*fr);
ratio = @(E) (mean(E(:,w_late),'all','omitnan') - mean(E(:,w_base),'all','omitnan')) / ...
             (max(mean(E(:,w_pk),'omitnan')) - mean(E(:,w_base),'all','omitnan'));
rng(0); rat = nan(3,3);
for g = 1:3
    n = size(ETA{g},1); b = arrayfun(@(~) ratio(ETA{g}(randi(n,n,1),:)),1:1000);
    rat(g,:) = [ratio(ETA{g}), prctile(b,[2.5 97.5])];
    fprintf('%-18s locked trials %2d/%2d  arrivals %3d  cue lag %+.2f s  sustained/peak %.2f [%.2f %.2f]\n', ...
        names{g}, nTr(g), numel(U{g}), n, median(lagFit{g}), rat(g,:));
end

%% figure
c = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.3 0.3 0.3];
figure(2); clf; set(gcf,'Color','w','Position',[100 100 1300 380]);
for g = 1:3
    subplot(1,4,g); hold on
    mu = mean(ETA{g},'omitnan'); se = std(ETA{g},'omitnan')/sqrt(size(ETA{g},1));
    fill([tt fliplr(tt)],[mu+se fliplr(mu-se)],c(g,:),'FaceAlpha',.2,'EdgeColor','none');
    plot(tt,mu,'Color',c(g,:),'LineWidth',2); plot(tt,mean(ETAo{g},'omitnan'),'k:','LineWidth',1)
    xline(0); xline(minstay/fr,'--');
    title(sprintf('%s (n=%d)',names{g},size(ETA{g},1))); xlabel('time from bump arrival (s)');
    if g == 1, ylabel('glomerulus - PB mean (z)'); end
end
subplot(1,4,4); hold on
for g = 1:3, errorbar(g,rat(g,1),rat(g,1)-rat(g,2),rat(g,3)-rat(g,1),'o','Color',c(g,:),'MarkerFaceColor',c(g,:),'LineWidth',1.5); end
xticks(1:3); xticklabels(names); xlim([.5 3.5]); ylim([0 1.1]); yline(1,':');
ylabel('sustained (2.5-4 s) / peak'); title('1 = sustained, 0 = fully transient')
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_bump_arrival_dynamics.png'),'Resolution',150);

%% decomposition at arrivals: bump glomerulus, PB mean, opposite glomerulus, |vr| (absolute, not mean-subtracted)
DEC = cell(1,3);
for g = 1:3
    X = struct('L',[],'M',[],'O',[],'SP',[]);
    for k = 1:numel(U{g})
        [gl,ok] = cue_bump_glom(U{g}(k), alpha, fr);
        if ~ok, continue, end
        u = U{g}(k); D = u.D; D(:,u.bad) = nan; nT = size(D,2);
        mD = mean(D,1,'omitnan'); D = D/(prctile(mD,95)-prctile(mD,5));      % per-trial global range
        hom = (D(1:16,:)+D(17:32,:))/2; mm = mean(hom,1);
        for t0 = arrival_frames(u, gl, alpha, fr, pre, post, minstay)
            w = t0-pre:t0+post;
            X.L(end+1,:) = hom(gl(t0),w); X.M(end+1,:) = mm(w);
            X.O(end+1,:) = hom(mod(gl(t0)+7,16)+1,w); X.SP(end+1,:) = u.sp(w);
        end
    end
    DEC{g} = X;
end
figure(3); clf; set(gcf,'Color','w','Position',[100 100 1300 650]);
for g = 1:3
    subplot(2,3,g); hold on
    plot(tt,mean(DEC{g}.L,'omitnan'),'Color',c(g,:),'LineWidth',2); plot(tt,mean(DEC{g}.M,'omitnan'),'k','LineWidth',2); plot(tt,mean(DEC{g}.O,'omitnan'),'k:')
    xline(0); xline(minstay/fr,'--'); title(names{g}); if g==1, ylabel('dF/F / trial global range'); end
    if g==1, legend({'bump glom','PB mean','opposite glom'},'Box','off','Location','best'); end
    subplot(2,3,g+3); plot(tt,mean(DEC{g}.SP,'omitnan'),'r','LineWidth',2); xline(0); xline(minstay/fr,'--')
    xlabel('time from bump arrival (s)'); if g==1, ylabel('|rotational speed| (rad/s)'); end
end
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_bump_arrival_decomposition.png'),'Resolution',150);

%% additive vs multiplicative: (bump - opposite glom) vs PB level, all closed-loop frames
% PB level binned into per-trial octiles; both axes divided by the top-bin
% level. Additive heading offset -> flat (slope ~0). Gain-like (fixed
% spatial contrast scaled by global drive) -> rises through ~origin.
% Robust part of this claim: whether the excess changes with level (slope).
% The intercept depends on the dF/F baseline (7th-pct F0) and is softer.
% GCaMP7f is left out: its off-bump dF/F sits near 0, so "PB level" is
% mostly the bump itself and the axes aren't comparable.
nq = 8; AB = cell(1,2); RL = cell(1,2);
for g = 1:2
    for k = 1:numel(U{g})
        [gl,ok] = cue_bump_glom(U{g}(k), alpha, fr);
        if ~ok, continue, end
        u = U{g}(k); D = u.D; nT = size(D,2);
        hom = (D(1:16,:)+D(17:32,:))/2; mm = mean(hom,1);
        ex  = hom(sub2ind(size(hom),gl,1:nT)) - hom(sub2ind(size(hom),mod(gl+7,16)+1,1:nT));
        okf = ~u.bad; qb = discretize(mm,[-inf, quantile(mm(okf),(1:nq-1)/nq), inf]);
        row = nan(2,nq);
        for q = 1:nq, sel = okf & qb==q; row(:,q) = [mean(mm(sel)); mean(ex(sel))]; end
        row = row/row(1,end);
        RL{g} = cat(3,RL{g},row);
        p = polyfit(row(1,:),row(2,:),1); AB{g}(end+1,:) = [p(2) p(1)];
    end
    fprintf('%-18s n=%d  excess vs level: intercept %.3f [%.3f %.3f]  slope %.3f [%.3f %.3f]\n', names{g}, size(AB{g},1), ...
        median(AB{g}(:,1)), prctile(AB{g}(:,1),[25 75]), median(AB{g}(:,2)), prctile(AB{g}(:,2),[25 75]));
end
fprintf('ranksum slope p = %.2g\n', ranksum(AB{1}(:,2),AB{2}(:,2)));

figure(4); clf; set(gcf,'Color','w','Position',[100 100 800 380]);
subplot(1,2,1); hold on
for g = 1:2
    plot(squeeze(RL{g}(1,:,:)),squeeze(RL{g}(2,:,:)),'-','Color',[c(g,:) .25])
    plot(mean(RL{g}(1,:,:),3),mean(RL{g}(2,:,:),3),'-o','Color',c(g,:),'MarkerFaceColor',c(g,:),'LineWidth',2)
end
xlabel('PB-mean dF/F (/ top octile)'); ylabel('bump - opposite glom (/ top octile level)'); title('heading excess vs global level')
subplot(1,2,2); hold on
for g = 1:2, swarmchart(g*ones(size(AB{g},1),1),AB{g}(:,2),25,c(g,:),'filled'); end
yline(0,':'); xticks(1:2); xticklabels(names(1:2)); xlim([.5 2.5]); ylabel('slope (per trial)')
title(sprintf('additive (0) vs gain-like; p = %.1g', ranksum(AB{1}(:,2),AB{2}(:,2))))
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_bump_excess_vs_level.png'),'Resolution',150);

%% functions
function u = from_lpsp_cl(s)
D = s.im.d; n = size(D,2); xf = s.ft.xf; xb = linspace(xf(1),xf(end),n);
u.D   = D;
u.cue = angle(interp1(xf,exp(1i*s.ft.cue),xb));
u.bad = flash_frames(D);
u.sp  = interp1(xf,abs(smoothdata(s.ft.r_speed,'gaussian',30)),xb);
end

function u = from_epg7f(s)
% dF/F from raw im.f the way lpsp_cl was processed: 5-frame gaussian, 7th-pct F0
f  = smoothdata(s.im.f,2,'gaussian',5);
f0 = prctile(f,7,2);
u.D   = (f - f0)./f0;
u.cue = s.ft.cue(:)';                % already on the volume clock (ft.xb)
u.bad = flash_frames(u.D) | logical(s.im.badVol(:)');
u.sp  = interp1(s.ft.xf,abs(smoothdata(s.ft.r_speed_raw,'gaussian',30)),s.ft.xb);
end

function bad = flash_frames(D)
% brief, global upward flashes: fast residual of the PB-mean trace > 6 robust SDs, dilated +/-3 frames
m   = mean(D,1);
r   = m - movmedian(m,21);
bad = r > 6*1.4826*mad(r,1);
bad = conv(double(bad),ones(1,7),'same') > 0;
end

function [ETA,ETAo,ok,lagS] = arrival_eta(u, alpha, fr, pre, post, minstay)
cd_ = @(a,b) abs(angle(exp(1i*(a-b))));
D = u.D; cue = u.cue; bad = u.bad; T = size(D,2);
ETA = []; ETAo = []; lagS = nan;
v = exp(1i*alpha(:)).'*(D - mean(D,2)); ph = angle(v);
okm = find(abs(v) > median(abs(v)) & ~bad);
best = 0;
for sgn = [-1 1]
    for t = -round(2.5*fr):round(1*fr)
        idx = okm(okm+t >= 1 & okm+t <= T); R = mean(exp(1i*(ph(idx) + sgn*cue(idx+t))));
        if abs(R) > best, best = abs(R); bs = sgn; bt = t; c = angle(R); end
    end
end
ok = best >= 0.4;
if ~ok, return, end
lagS = bt/fr;
pred = angle(exp(1i*(c - bs*circshift(cue,-bt))));
[~,gl] = min(cd_(pred(:),alpha(1:16)),[],2); gl = gl';
m  = mean(D,1);
Dn = D - m; Dn(:,bad) = nan; Dn = Dn/std(Dn(:),'omitnan');
hom = (Dn(1:16,:) + Dn(17:32,:))/2;
lb = round(1.5*fr); last = -inf;
for t0 = max(pre,lb)+1:T-post
    j = gl(t0);
    if gl(t0-1) ~= j && cd_(pred(t0-lb),alpha(j)) >= pi/3 && ...
       all(cd_(pred(t0:t0+minstay),alpha(j)) <= 1.5*pi/8) && t0-last > minstay
        ETA(end+1,:)  = hom(j,t0-pre:t0+post);                %#ok<AGROW>
        ETAo(end+1,:) = hom(mod(j+7,16)+1,t0-pre:t0+post);    %#ok<AGROW> opposite glomerulus
        last = t0;
    end
end
end

function [gl,ok,pred] = cue_bump_glom(u, alpha, fr)
% per-trial cue -> bump-glomerulus (1..16) prediction, same fit as arrival_eta
cd_ = @(a,b) abs(angle(exp(1i*(a-b))));
D = u.D; cue = u.cue; T = size(D,2); gl = []; pred = [];
v = exp(1i*alpha(:)).'*(D - mean(D,2)); ph = angle(v);
okm = find(abs(v) > median(abs(v)) & ~u.bad); best = 0;
for sgn = [-1 1]
    for t = -round(2.5*fr):round(1*fr)
        idx = okm(okm+t >= 1 & okm+t <= T); R = mean(exp(1i*(ph(idx) + sgn*cue(idx+t))));
        if abs(R) > best, best = abs(R); bs = sgn; bt = t; c = angle(R); end
    end
end
ok = best >= 0.4;
if ~ok, return, end
pred = angle(exp(1i*(c - bs*circshift(cue,-bt))));
[~,gl] = min(cd_(pred(:),alpha(1:16)),[],2); gl = gl';
end

function t_arr = arrival_frames(u, gl, alpha, fr, pre, post, minstay)
cd_ = @(a,b) abs(angle(exp(1i*(a-b))));
[~,~,pred] = cue_bump_glom(u, alpha, fr);
T = numel(gl); lb = round(1.5*fr); last = -inf; t_arr = [];
for t0 = max(pre,lb)+1:T-post
    j = gl(t0);
    if gl(t0-1) ~= j && cd_(pred(t0-lb),alpha(j)) >= pi/3 && ...
       all(cd_(pred(t0:t0+minstay),alpha(j)) <= 1.5*pi/8) && t0-last > minstay
        t_arr(end+1) = t0; last = t0; %#ok<AGROW>
    end
end
end
