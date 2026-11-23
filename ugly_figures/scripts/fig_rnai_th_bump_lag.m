%% fig_rnai_th_bump_lag
% Phase relationship between the bump and the fly's heading, TH set: does
% the bump turn BEFORE the fly in LPsP > TH-RNAi flies, or is the bump's
% lag behind the fly at least shorter? geno_set = 'th' (default) or 'vglut'
% (LPsP>VGlut, Empty>VGlut, LPsP>mCherry; exports tagged by set). TH groups: LPsP>TH (all), split into
% high / low closed-loop integrative gain, Empty>TH, LPsP>mCherry.
%
%   (a-c) CONTINUOUS cross-correlation: per trial, angular velocity of the
%     bump (unwrapped PVA) and of the fly (unwrapped heading), both on the
%     imaging timebase with the SAME smoothing (xc_smooth_s Gaussian on the
%     angle, then d/dt), correlated at lags of +-max_lag_s using only
%     samples with PVA strength >= rho_min and no raw PVA jump > jump_max_deg
%     at both times; a fly's curve = mean of its trials' curves weighted by
%     samples. Peak lag (parabolic interpolation around the max) = that
%     fly's lag; POSITIVE = the bump lags the fly, negative = it leads.
%     Closed loop (integrative-gain inclusion) and dark (walking + bump;
%     heading corrected for the mounting bias, as in
%     fig_rnai_th_lapping_examples.m).
%   (d-e) TURN-ALIGNED: closed-loop turns >= min_turn_deg from the turn cache
%     (lpsp_rnai_turns_claude.m; fly and bump velocity traces smoothed
%     identically, 0.5 s, bump NOT shifted), per-fly mean traces (>=
%     min_turns turns), each normalized to its own peak so gain doesn't
%     matter; mean +- SEM across flies. Per fly: cross-correlation peak lag of
%     the mean traces, and the time each trace first reaches half its peak
%     (bump minus fly; negative = the bump's turn starts first).
% Note the calcium sensor itself adds a delay (~0.1-0.15 s here), so a
% "true" bump lead would show as a lag shrinking toward or below zero.
% A second figure quantifies PERSISTENCE (the bump rotating on after the fly
% stops): cross-correlation tail length, and turn-offset-aligned coasting.
% p: ranksum on per-fly values, each LPsP>TH group vs Empty>TH and vs
% LPsP>mCherry, and high vs low gain. Needs all_data (loaded if absent), the
% integrative-gain cache and the turn cache. Run from the repo root.

%% parameters
xc_smooth_s = 0.2; max_lag_s = 1.0; rho_min = 0.4; jump_max_deg = 90; min_trial_samples = 300;
min_turn_deg = 30; min_turns = 5; gain_split = 0.85;
geno_set = 'th'; % 'th' (LPsP>TH split by gain), 'th_pooled' or 'vglut'
switch geno_set
    case 'th'
        grp_def = {'LPsP>TH (all)','lpsp>th',[]; 'LPsP>TH high gain','lpsp>th',true; 'LPsP>TH low gain','lpsp>th',false; ...
                   'Empty > TH-RNAi','empty>th',[]; 'LPsP > mCherry','lpsp>mcherry',[]};
        grp_col = [0.85,0.33,0.10; 0.60,0.15,0.00; 0.95,0.62,0.45; 0.00,0.45,0.74; 0.93,0.69,0.13];
        grp_x = [1,2,3,4.4,5.4]; cmp = [1,4; 1,5; 2,4; 2,3];
        geno3 = {'lpsp>th','empty>th','lpsp>mcherry'}; show = [2,3,4,5]; % curves: high, low, Empty, mCherry
        set_tag = 'th'; set_lbl = 'TH set';
    case 'th_pooled' % as 'th', but all LPsP>TH flies together
        grp_def = {'LPsP > TH-RNAi','lpsp>th',[]; 'Empty > TH-RNAi','empty>th',[]; 'LPsP > mCherry','lpsp>mcherry',[]};
        grp_col = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13];
        grp_x = 1:3; cmp = [1,2; 1,3];
        geno3 = {'lpsp>th','empty>th','lpsp>mcherry'}; show = 1:3;
        set_tag = 'th_pooled'; set_lbl = 'TH set';
    case 'vglut'
        grp_def = {'LPsP > VGlut-RNAi','lpsp>vglut',[]; 'Empty > VGlut-RNAi','empty>vglut',[]; 'LPsP > mCherry','lpsp>mcherry',[]};
        grp_col = [0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13];
        grp_x = 1:3; cmp = [1,2; 1,3];
        geno3 = {'lpsp>vglut','empty>vglut','lpsp>mcherry'}; show = 1:3;
        set_tag = 'vglut'; set_lbl = 'VGlut set';
end
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
TC = load(fullfile(data_dir,'lpsp_rnai_turns_cache.mat'),'turns','turn_t');
if ~exist('all_data','var')
    tmp = load(fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
n_flies = max(IG.fly_num); g_fly = nan(n_flies,1);
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = IG.genotype{find(IG.fly_num==f,1)}; end
grp_flies = cell(size(grp_def,1),1);
for gi = 1:size(grp_def,1)
    fl = find(strcmp(fly_geno,grp_def{gi,2}));
    if ~isempty(grp_def{gi,3}), fl = fl(~isnan(g_fly(fl)) & (g_fly(fl) > gain_split) == grp_def{gi,3}); end
    grp_flies{gi} = fl;
end
out_dir = fullfile('ugly_figures','rnai_turns','bump_lag');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% continuous cross-correlation, closed loop and dark
dt_ref = 0.1; lags_s = (-round(max_lag_s/dt_ref):round(max_lag_s/dt_ref))*dt_ref;
XC = struct();
for cond = {'closed','dark'}
    if strcmp(cond{1},'dark')
        cand = find(ismember(IG.genotype,geno3) & IG.is_dark); ok = false(size(cand));
        for j = 1:numel(cand), tr = all_data(cand(j)); ok(j) = sum((tr.ft.f_speed > .1) | (abs(tr.ft.r_speed)*2 > .1)) > numel(tr.ft.f_speed)*.2 && mean(tr.im.rho) > .2; end
        trials = cand(ok);
    else
        trials = find(ismember(IG.genotype,geno3) & ~IG.is_dark & IG.inc);
    end
    R = zeros(n_flies,numel(lags_s)); W = zeros(n_flies,1);
    for i = trials(:)'
        [r,n] = trial_xcorr(all_data(i),lags_s,xc_smooth_s,rho_min,jump_max_deg,strcmp(cond{1},'dark'));
        if n < min_trial_samples || any(isnan(r)), continue, end
        f = IG.fly_num(i); R(f,:) = R(f,:) + n*r; W(f) = W(f) + n;
    end
    R = R./W; R(W==0,:) = NaN;
    pk = nan(n_flies,1); for f = find(W > 0)', pk(f) = peak_lag(R(f,:),lags_s); end
    XC.(cond{1}).R = R; XC.(cond{1}).pk = pk;
    fprintf('\n%s, continuous bump-vs-fly velocity cross-correlation peak lag (s; + = bump lags):\n',cond{1});
    for gi = 1:size(grp_def,1), v = pk(grp_flies{gi}); v = v(~isnan(v)); fprintf('  %-18s median %+.3f s, mean %+.3f +/- %.3f (n=%d flies), %d flies with lead (< 0)\n',grp_def{gi,1},median(v),mean(v),std(v)/sqrt(numel(v)),numel(v),sum(v < 0)); end
    p = stats_p(pk,grp_flies,cmp); XC.(cond{1}).p = p;
    for c = 1:size(cmp,1), fprintf('  %s vs %s: p = %.3g\n',grp_def{cmp(c,1),1},grp_def{cmp(c,2),1},p(c)); end
end

%% turn-aligned, closed loop
tt = [TC.turns.trial]'; t = TC.turn_t(:);
sel = find(~[TC.turns.dark]' & IG.inc(tt) & [TC.turns.size_deg]' >= min_turn_deg & ismember(IG.genotype(tt),geno3));
tf = IG.fly_num(tt(sel));
FV = nan(n_flies,numel(t)); BV = FV; tl_lag = nan(n_flies,1); half_dt = nan(n_flies,1);
lags_t = -round(max_lag_s/median(diff(t))):round(max_lag_s/median(diff(t)));
for f = unique(tf)'
    k = sel(tf==f); if numel(k) < min_turns, continue, end
    fv = mean(cat(2,TC.turns(k).fly_vel_trace),2,'omitnan'); bv = mean(cat(2,TC.turns(k).bump_vel_trace),2,'omitnan');
    FV(f,:) = fv/max(fv); BV(f,:) = bv/max(bv);
    r = arrayfun(@(L) lagcorr(fv,bv,L),lags_t); tl_lag(f) = peak_lag(r,lags_t*median(diff(t)));
    i_f = find(t >= -0.5 & FV(f,:)' >= 0.5,1); i_b = find(t >= -0.5 & BV(f,:)' >= 0.5,1);
    if ~isempty(i_f) && ~isempty(i_b), half_dt(f) = t(i_b) - t(i_f); end
end
p_tl = stats_p(tl_lag,grp_flies,cmp); p_half = stats_p(half_dt,grp_flies,cmp);
fprintf('\nturn-aligned (closed loop, turns >= %d deg, flies with >= %d turns):\n',min_turn_deg,min_turns);
for gi = 1:size(grp_def,1)
    a = tl_lag(grp_flies{gi}); a = a(~isnan(a)); b = half_dt(grp_flies{gi}); b = b(~isnan(b));
    fprintf('  %-18s xcorr lag median %+.3f s (n=%d); half-rise bump - fly median %+.3f s\n',grp_def{gi,1},median(a),numel(a),median(b));
end
for c = 1:size(cmp,1), fprintf('  %s vs %s: xcorr lag p = %.3g, half-rise p = %.3g\n',grp_def{cmp(c,1),1},grp_def{cmp(c,2),1},p_tl(c),p_half(c)); end

%% figure
figure('Color','w','Position',[40,40,1800,900]); clf
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
conds = {'closed','dark'}; ttl = {'closed loop','dark (bias-corrected heading)'};
for ci = 1:2
    nexttile(tl,(ci-1)*4+1); hold on
    for gi = show
        Rg = XC.(conds{ci}).R(grp_flies{gi},:); Rg = Rg(~any(isnan(Rg),2),:);
        m = mean(Rg,1); s = std(Rg,0,1)/sqrt(size(Rg,1));
        patch([lags_s, fliplr(lags_s)],[m+s, fliplr(m-s)],grp_col(gi,:),'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off')
        plot(lags_s,m,'Color',grp_col(gi,:),'LineWidth',1.8,'DisplayName',sprintf('%s (n=%d)',grp_def{gi,1},size(Rg,1)))
    end
    xline(0,':k','HandleVisibility','off'); xlabel('lag (s; + = bump lags the fly)'); ylabel('corr(bump vel, fly vel)')
    legend('Location','southwest','FontSize',8); title({sprintf('(%s) %s: cross-correlation',char('a'+(ci-1)*2),ttl{ci}),'mean \pm SEM across flies'},'FontWeight','normal','FontSize',9)
    nexttile(tl,(ci-1)*4+2); hold on
    dot_panel(XC.(conds{ci}).pk,grp_flies,grp_def,grp_col,grp_x,cmp,XC.(conds{ci}).p)
    ylabel('peak lag (s; + = bump lags)'); title({sprintf('(%s) %s: per-fly peak lag',char('b'+(ci-1)*2),ttl{ci}),'ranksum'},'FontWeight','normal','FontSize',9)
end
nexttile(tl,3); hold on
for gi = show
    Fg = FV(grp_flies{gi},:); Bg = BV(grp_flies{gi},:); ok = ~any(isnan(Fg),2) & ~any(isnan(Bg),2);
    mb = mean(Bg(ok,:),1); sb = std(Bg(ok,:),0,1)/sqrt(sum(ok));
    patch([t', fliplr(t')],[mb+sb, fliplr(mb-sb)],grp_col(gi,:),'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off')
    plot(t,mb,'Color',grp_col(gi,:),'LineWidth',1.8,'DisplayName',sprintf('bump, %s (n=%d)',grp_def{gi,1},sum(ok)))
end
Fall = FV(cell2mat(grp_flies(show)),:); Fall = Fall(~any(isnan(Fall),2),:);
plot(t,mean(Fall,1),'k','LineWidth',2,'DisplayName','fly (all groups)')
xline(0,':k','HandleVisibility','off'); yline(0.5,':','Color',[.6,.6,.6],'HandleVisibility','off'); xlim([-1,2])
xlabel('time from turn onset (s)'); ylabel('angular velocity (normalized to peak)'); legend('Location','northeast','FontSize',7)
title({'(e) turn-aligned velocity, closed loop',sprintf('per-fly mean of turns >= %d%s, normalized; mean \\pm SEM',min_turn_deg,char(176))},'FontWeight','normal','FontSize',9)
nexttile(tl,4); hold on
dot_panel(tl_lag,grp_flies,grp_def,grp_col,grp_x,cmp,p_tl)
ylabel('peak lag of turn-aligned traces (s)'); title({'(f) turn-aligned cross-correlation lag','ranksum'},'FontWeight','normal','FontSize',9)
nexttile(tl,7); hold on
dot_panel(half_dt,grp_flies,grp_def,grp_col,grp_x,cmp,p_half)
ylabel('half-rise time, bump - fly (s)'); title({'(g) turn onset: bump reaches half its peak velocity','relative to the fly (negative = bump first)'},'FontWeight','normal','FontSize',9)
nexttile(tl,8); axis off
text(0,1,{'\bfMethod\rm', sprintf('(a-d) bump and fly angular velocity, %.1f s smoothing, both on the',xc_smooth_s), ...
    sprintf('imaging timebase; samples with \\rho >= %g and no PVA jump; peak lag by',rho_min), 'parabolic interpolation (+ = bump lags the fly).', ...
    '(e-g) closed-loop turns from the turn cache; bump not shifted.', '', ...
    'The GCaMP signal itself delays the bump by ~0.1-0.15 s,', 'so a bump that truly led would show a lag near or below 0.'},'VerticalAlignment','top','FontSize',9)
title(tl,['Bump-heading phase: does the bump turn before the fly? (' set_lbl ')'],'FontSize',12)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_bump_lag.png',set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_bump_lag.pdf',set_tag)),'ContentType','vector')

%% persistence: does the bump keep rotating after the fly stops?
% (1) cross-correlation tail: per fly, from the curves above, the time from the
%     peak until the correlation falls below half its peak on the positive-lag
%     side (bump later) and on the negative side (bump earlier); interpolated
%     on a fine lag grid. A bump that keeps rotating after the fly gives a
%     long positive-side tail; the negative side is a within-fly reference.
% (2) turn-offset aligned (turn cache, turns >= min_turn_deg lasting <= max_dur_s
%     so post_win_s after the turn's end is inside the cached -1..3 s window):
%     the time after offset at which each fly's mean (peak-normalized, so
%     gain-independent) bump / fly velocity first falls below half its peak;
%     bump minus fly. (A post-offset rotation share was tried and dropped:
%     small counter-turns after offset make it unstable.)
max_dur_s = 1.5; post_win_s = 1.5;
lf = (lags_s(1):0.01:lags_s(end));
P = struct();
for cond = {'closed','dark'}
    R = XC.(cond{1}).R; tail_r = nan(n_flies,1); tail_l = nan(n_flies,1);
    for f = find(~any(isnan(R),2))'
        rf = interp1(lags_s,R(f,:),lf,'pchip'); [pk_r,k] = max(rf); h = pk_r/2;
        jr = find(rf(k:end) < h,1); jl = find(fliplr(rf(1:k)) < h,1);
        if ~isempty(jr), tail_r(f) = (jr-1)*0.01; end
        if ~isempty(jl), tail_l(f) = (jl-1)*0.01; end
    end
    P.(cond{1}).tail_r = tail_r; P.(cond{1}).tail_l = tail_l;
    % turn-offset aligned
    is_dk = strcmp(cond{1},'dark');
    if is_dk, okt = IG.is_dark(tt); else, okt = ~IG.is_dark(tt) & IG.inc(tt); end
    sel2 = find(okt & [TC.turns.size_deg]' >= min_turn_deg & [TC.turns.dur_s]' <= max_dur_s & ismember(IG.genotype(tt),geno3));
    dtt = median(diff(t)); t_off = (-0.5:dtt:post_win_s)';
    coast_b = nan(n_flies,1); coast_f = coast_b; dec_b = coast_b; dec_f = coast_b; OB = nan(n_flies,numel(t_off)); OF = OB;
    tf2 = IG.fly_num(tt(sel2));
    for f = unique(tf2)'
        k = sel2(tf2==f); if numel(k) < min_turns, continue, end
        cb = nan(numel(k),1); cf = cb; vb = nan(numel(t_off),numel(k)); vf = vb;
        for j = 1:numel(k)
            T = TC.turns(k(j)); tq = T.dur_s + t_off; iE = find(t <= T.dur_s + post_win_s,1,'last'); iO = find(t >= T.dur_s,1); i0 = find(t >= 0,1);
            btr = double(T.bump_trace); ftr = double(T.fly_trace);
            tot_b = btr(iE) - btr(i0); tot_f = ftr(iE) - ftr(i0);
            if tot_b > 0, cb(j) = (btr(iE) - btr(iO))/tot_b; end
            if tot_f > 0, cf(j) = (ftr(iE) - ftr(iO))/tot_f; end
            vb(:,j) = interp1(t,double(T.bump_vel_trace),tq); vf(:,j) = interp1(t,double(T.fly_vel_trace),tq);
        end
        coast_b(f) = mean(cb,'omitnan'); coast_f(f) = mean(cf,'omitnan');
        mb = mean(vb,2,'omitnan'); mf = mean(vf,2,'omitnan'); OB(f,:) = mb/max(mb); OF(f,:) = mf/max(mf);
        ib = find(t_off >= 0 & OB(f,:)' < 0.5,1); iff = find(t_off >= 0 & OF(f,:)' < 0.5,1);
        if ~isempty(ib), dec_b(f) = t_off(ib); else, dec_b(f) = post_win_s; end
        if ~isempty(iff), dec_f(f) = t_off(iff); else, dec_f(f) = post_win_s; end
    end
    P.(cond{1}).coast_b = coast_b; P.(cond{1}).coast_f = coast_f; P.(cond{1}).dec = dec_b - dec_f; P.(cond{1}).OB = OB; P.(cond{1}).OF = OF; P.(cond{1}).t_off = t_off;
    fprintf('\n%s persistence (per fly; medians):\n',cond{1});
    for gi = 1:size(grp_def,1)
        g = grp_flies{gi}; fprintf('  %-18s xcorr tail: + side %.2f s, - side %.2f s | bump half-decay after offset minus fly %+.2f s\n', grp_def{gi,1}, ...
            median(tail_r(g),'omitnan'), median(tail_l(g),'omitnan'), median(P.(cond{1}).dec(g),'omitnan'));
    end
    P.(cond{1}).p_tail = stats_p(tail_r,grp_flies,cmp); P.(cond{1}).p_asym = stats_p(tail_r - tail_l,grp_flies,cmp);
    P.(cond{1}).p_coast = stats_p(coast_b - coast_f,grp_flies,cmp); P.(cond{1}).p_dec = stats_p(P.(cond{1}).dec,grp_flies,cmp);
    for c = 1:size(cmp,1), fprintf('  %s vs %s: tail p = %.3g, tail asymmetry p = %.3g, decay p = %.3g\n',grp_def{cmp(c,1),1},grp_def{cmp(c,2),1}, ...
        P.(cond{1}).p_tail(c),P.(cond{1}).p_asym(c),P.(cond{1}).p_dec(c)); end
end

%% figure: persistence
figure('Color','w','Position',[40,40,1800,900]); clf
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
for ci = 1:2
    c = conds{ci}; Q = P.(c);
    nexttile(tl,(ci-1)*4+1); hold on
    for gi = show
        Bg = Q.OB(grp_flies{gi},:); Bg = Bg(~any(isnan(Bg),2),:);
        m = mean(Bg,1); s = std(Bg,0,1)/sqrt(size(Bg,1));
        patch([Q.t_off', fliplr(Q.t_off')],[m+s, fliplr(m-s)],grp_col(gi,:),'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off')
        plot(Q.t_off,m,'Color',grp_col(gi,:),'LineWidth',1.8,'DisplayName',sprintf('bump, %s (n=%d)',grp_def{gi,1},size(Bg,1)))
    end
    Fall = Q.OF(cell2mat(grp_flies(show)),:); Fall = Fall(~any(isnan(Fall),2),:);
    plot(Q.t_off,mean(Fall,1),'k','LineWidth',2,'DisplayName','fly (all groups)')
    xline(0,':k','HandleVisibility','off'); yline(0.5,':','Color',[.6,.6,.6],'HandleVisibility','off'); xlim(Q.t_off([1 end]))
    xlabel('time from turn OFFSET (s)'); ylabel('angular velocity (normalized)'); legend('Location','northeast','FontSize',7)
    title({sprintf('(%s) %s: velocity aligned to the end of the turn',char('a'+(ci-1)*4),ttl{ci}),sprintf('turns >= %d%s lasting <= %.1f s; per-fly means, normalized',min_turn_deg,char(176),max_dur_s)},'FontWeight','normal','FontSize',9)
    nexttile(tl,(ci-1)*4+2); hold on
    dot_panel(Q.dec,grp_flies,grp_def,grp_col,grp_x,cmp,Q.p_dec)
    ylabel('half-decay after offset: bump - fly (s)'); title({sprintf('(%s) how long the bump keeps turning',char('b'+(ci-1)*4)),'time to fall to half peak velocity after offset, bump minus fly'},'FontWeight','normal','FontSize',9)
    nexttile(tl,(ci-1)*4+3); hold on
    dot_panel(Q.tail_r,grp_flies,grp_def,grp_col,grp_x,cmp,Q.p_tail)
    ylabel('cross-correlation tail, + side (s)'); title({sprintf('(%s) %s: correlation persistence',char('c'+(ci-1)*4),ttl{ci}),'time from the xcorr peak to half its height, bump-later side'},'FontWeight','normal','FontSize',9)
    nexttile(tl,(ci-1)*4+4); hold on
    dot_panel(Q.tail_r - Q.tail_l,grp_flies,grp_def,grp_col,grp_x,cmp,Q.p_asym)
    ylabel('tail asymmetry: + side - (-) side (s)'); title({sprintf('(%s) tail asymmetry',char('d'+(ci-1)*4)),'> 0: correlation lingers longer with the bump later'},'FontWeight','normal','FontSize',9)
end
title(tl,['Does the bump keep rotating after the fly stops? (' set_lbl '; p = ranksum on per-fly values)'],'FontSize',12)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_bump_persistence.png',set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_bump_persistence.pdf',set_tag)),'ContentType','vector')

%% functions
function [r,n] = trial_xcorr(tr,lags_s,sm_s,rho_min,jump_max_deg,dark)
    % weighted-by-validity Pearson correlation of bump vs fly angular velocity at each lag
    xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    hd = unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'));
    if dark % mounting-bias correction (heading = integral of r_speed)
        wf = max(1,round(0.5/median(diff(xf)))); rr = smoothdata(fillmissing(tr.ft.r_speed(:),'constant',0),'gaussian',wf); ff = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',wf);
        m = ff > 1; if sum(m) > 200, b = [ones(sum(m),1) ff(m)] \ rr(m); hd = hd - cumtrapz(xf,b(2)*ff); end
    end
    dt = median(diff(xb)); w = max(1,round(sm_s/dt));
    bu = smoothdata(unwrap(mu),'gaussian',w); hb = smoothdata(interp1(xf,hd,xb),'gaussian',w);
    bv = gradient(bu)/dt; fv = gradient(hb)/dt;
    good = rho >= rho_min & ~([0; abs(diff(unwrap(mu)))] > deg2rad(jump_max_deg)) & ~isnan(fv);
    good = good & ~movmax(~good,[w w]);  % keep away from excluded stretches
    L = round(lags_s/dt); r = nan(size(L)); n = sum(good);
    for k = 1:numel(L)
        if L(k) >= 0, a = fv(1:end-L(k)); c = bv(1+L(k):end); g = good(1:end-L(k)) & good(1+L(k):end);
        else, a = fv(1-L(k):end); c = bv(1:end+L(k)); g = good(1-L(k):end) & good(1:end+L(k)); end
        if sum(g) > 50, r(k) = corr(a(g),c(g)); end
    end
end

function r = lagcorr(fv,bv,L)
    % corr of bv(t) with fv(t - L); positive L = bump lags
    if L >= 0, a = fv(1:end-L); c = bv(1+L:end); else, a = fv(1-L:end); c = bv(1:end+L); end
    ok = ~isnan(a) & ~isnan(c); r = corr(a(ok),c(ok));
end

function lag = peak_lag(r,lags)
    % lag of the maximum, refined by a parabola through the peak and its neighbours
    [~,k] = max(r); lag = lags(k);
    if k > 1 && k < numel(r)
        y = r(k-1:k+1); den = y(1) - 2*y(2) + y(3);
        if den < 0, lag = lags(k) + 0.5*(y(1) - y(3))/den*(lags(2)-lags(1)); end
    end
end

function p = stats_p(v,grp_flies,cmp)
    p = nan(size(cmp,1),1);
    for c = 1:size(cmp,1)
        a = v(grp_flies{cmp(c,1)}); b = v(grp_flies{cmp(c,2)}); a = a(~isnan(a)); b = b(~isnan(b));
        if numel(a) >= 2 && numel(b) >= 2, p(c) = ranksum(a,b); end
    end
end

function dot_panel(v,grp_flies,grp_def,grp_col,grp_x,cmp,p)
    % per-fly dots + mean +- SEM per group, significance brackets
    for gi = 1:numel(grp_flies)
        y = v(grp_flies{gi}); y = y(~isnan(y));
        scatter(grp_x(gi) + (rand(size(y))-.5)*.35,y,20,grp_col(gi,:),'filled','MarkerFaceAlpha',.55)
        errorbar(grp_x(gi),mean(y),std(y)/sqrt(numel(y)),'o','Color','k','MarkerFaceColor',grp_col(gi,:),'MarkerSize',7,'LineWidth',1.2)
        n(gi) = numel(y); %#ok<AGROW>
    end
    yline(0,':k'); yl = ylim; ylim(yl);
    top = yl(2); h = 0.06*diff(yl);
    for c = 1:size(cmp,1)
        y = top + (c-0.4)*h*1.6;
        plot(grp_x(cmp(c,[1,1,2,2])),[y-h*0.4,y,y,y-h*0.4],'k','LineWidth',0.8)
        text(mean(grp_x(cmp(c,:))),y,sprintf('p = %.2g',p(c)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',7)
    end
    ylim([yl(1), top + (size(cmp,1)+0.6)*h*1.6])
    xticks(grp_x); xticklabels(arrayfun(@(gi) sprintf('%s (n=%d)',grp_def{gi,1},n(gi)),1:numel(grp_flies),'UniformOutput',false)); xtickangle(30); xlim([0.4,max(grp_x)+0.6])
end
