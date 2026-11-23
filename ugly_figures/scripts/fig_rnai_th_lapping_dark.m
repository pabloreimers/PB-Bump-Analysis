%% fig_rnai_th_lapping_dark
% Lapping in the DARK, TH set. With no bar, a lap is the bump over-rotating
% the fly's heading (FicTrac). The catch: a fly mounted at a slight angle to
% the ball reads forward walking as slow rotation, so heading drifts while
% the fly walks straight and the bump (correctly) doesn't follow -- a
% bump-heading mismatch that isn't a lap. So only turns where that mounting
% term is small are used: per trial, the leak of forward speed into measured
% turning is estimated (slope below), and a turn is kept if the heading the
% leak alone would add over the scored interval is < bias_max_deg (well below
% the lap threshold). Trials with too little forward walking to estimate the
% leak fall back to requiring little forward walking around the turn (mean
% |forward| < fwd_max). (An earlier version required little forward walking
% for every turn; that dropped half the high-gain LPsP>TH flies, which turn
% while walking.)
%
%   turns: from the turn cache (lpsp_rnai_turns_claude.m), >= min_turn_deg,
%     mean PVA strength over the turn >= rho_min; trials with walking and a
%     bump (mean raw rho > 0.2; the integrative-gain cue test is skipped in
%     the dark since there is no bar to move). Closed-loop turns from the
%     included closed-loop trials are scored the same way for comparison.
%   overshoot (lap): the bump rotates >= overshoot_deg past the fly, in the
%     turn's direction, from turn onset to post_s after its end -- as in
%     fig_rnai_th_lap_predictors.m. Turns with a raw frame-to-frame PVA jump >
%     jump_max_deg in that interval are dropped (tracking glitches).
%   mounting bias: per trial, slope of smoothed turning speed on forward
%     speed during walking (forward > 1); per turn, the heading that bias
%     alone would add over the scored interval.
%   controls for spontaneous bump drift in the dark (the bump is less
%     anchored without a visual cue): (1) "opposite" -- the bump moving >=
%     |fly| + overshoot_deg the other way, which drift would produce as often
%     as a same-direction overshoot but a turn-driven lap would not; (2) a
%     still-window baseline -- the same measure in a random window of the
%     turn's duration + post_s from the same trial where the fly is still
%     (|turning| < 0.3 rad/s, |forward| < fwd_max) and the bump solid.
%   groups: LPsP>TH (pooled, or split by closed-loop integrative gain with
%     split_gain = true; exports then lose the _pooled tag), Empty>TH, LPsP>mCherry. p: ranksum on per-fly rates (flies with >=
%     min_fly_turns usable dark turns); signrank dark vs closed loop for
%     high-gain flies.
% Needs all_data (loaded if absent), the integrative-gain cache and the turn
% cache. Run from the repo root.

%% parameters
min_turn_deg = 30; overshoot_deg = 150; post_s = 1.0; rho_min = 0.4; fwd_max = 1; gain_split = 0.85;
bias_max_deg = 20; % keep turns whose estimated mounting-bias heading over the scored interval is below this
rho_p10_min = 0.3; % 10th percentile of PVA strength over the whole scored interval (onset to post_s after the turn)
jump_max_deg = 90; % turns with a raw frame-to-frame PVA jump larger than this in the scored interval are dropped (tracking glitches)
min_fly_turns = 5; n_examples = 4; pad_s = 2;
geno_set = {'lpsp>th','empty>th','lpsp>mcherry'};
split_gain = false; % true: LPsP>TH split into high / low closed-loop integrative gain; false: one LPsP>TH group
if split_gain
    grps = {'lpsp>th hi','lpsp>th lo','empty>th','lpsp>mcherry'};
    grp_names = {'LPsP>TH high gain','LPsP>TH low gain','Empty > TH-RNAi','LPsP > mCherry'};
    grp_col = [0.60,0.15,0.00; 0.95,0.62,0.45; 0.00,0.45,0.74; 0.93,0.69,0.13];
    cmp = [1,3; 1,4; 1,2]; main_grp = 'lpsp>th hi'; main_lbl = 'high-gain LPsP>TH flies'; tag = '';
else
    grps = {'lpsp>th','empty>th','lpsp>mcherry'};
    grp_names = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'};
    grp_col = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13];
    cmp = [1,2; 1,3]; main_grp = 'lpsp>th'; main_lbl = 'LPsP>TH flies'; tag = '_pooled';
end
n_g = numel(grps);
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
out_dir = fullfile('ugly_figures','rnai_turns','laps_stalls');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
rng(11)

%% per-turn table, dark and closed loop
ttr = [TC.turns.trial]';
R = zeros(0,14); % trial fly dark dfly dbump rho fwd bias dur exc t0 null_exc jumpy rho_p10
for i = find(ismember(IG.genotype,geno_set))'
    tr = all_data(i);
    walk = sum((tr.ft.f_speed > .1) | (abs(tr.ft.r_speed)*2 > .1)) > numel(tr.ft.f_speed)*.2;
    if ~walk || mean(tr.im.rho) <= .2 || (~IG.is_dark(i) && ~IG.inc(i)), continue, end
    xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    mu_u = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/median(diff(xb)))));
    jmp = [0; abs(diff(unwrap(mu)))] > deg2rad(jump_max_deg);
    hd = unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'));
    w = max(1,round(0.5/median(diff(xf))));
    r = smoothdata(fillmissing(tr.ft.r_speed(:),'constant',0),'gaussian',w);
    f = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',w);
    m = f > 1; slope = NaN; if sum(m) > 200, b = [ones(sum(m),1) f(m)] \ r(m); slope = b(2); end % mounting bias, rad/s per forward unit
    still = interp1(xf,double(abs(r) < 0.3 & abs(f) < fwd_max),xb) > 0.99 & rho >= rho_min;
    cs = cumsum([0; ~still]);
    for k = find(ttr==i)'
        T = TC.turns(k); t0 = T.onset_s; t1 = t0 + T.dur_s;
        b0 = find(xb >= t0,1); b1 = find(xb <= t1 + post_s,1,'last'); bt = find(xb <= t1,1,'last');
        f0 = find(xf >= t0,1); f1 = find(xf <= t1,1,'last');
        if isempty(b0) || isempty(b1) || isempty(bt) || bt < b0 || isempty(f1) || f1 <= f0, continue, end
        dfly = rad2deg(hd(f1) - hd(f0)); if abs(dfly) < min_turn_deg, continue, end
        dbump = rad2deg(mu_u(b1) - mu_u(b0));
        L = b1 - b0; ok = find(cs(1+L:end) - cs(1:end-L) == 0); null_exc = NaN;
        if ~isempty(ok), a = ok(randi(numel(ok))); null_exc = sign(dfly)*rad2deg(mu_u(a+L-1) - mu_u(a)) - abs(dfly); end
        R(end+1,:) = [i, IG.fly_num(i), IG.is_dark(i), dfly, dbump, mean(rho(b0:bt)), mean(abs(f(xf >= t0-0.5 & xf <= t1+0.5))), ...
            rad2deg(abs(slope*mean(f(xf >= t0 & xf <= t1+post_s)))*(T.dur_s+post_s)), T.dur_s, sign(dfly)*dbump - abs(dfly), t0, null_exc, any(jmp(b0:b1)), prctile(rho(b0:b1),10)]; %#ok<AGROW>
    end
end
DT = array2table(R,'VariableNames',{'trial','fly','dark','dfly','dbump','rho','fwd','bias','dur','exc','t0','null_exc','jumpy','rho_p10'});
DT.geno = IG.genotype(DT.trial); DT.gain = g_fly(DT.fly); DT.grp = DT.geno;
if split_gain, m = strcmp(DT.geno,'lpsp>th'); DT.grp(m & DT.gain > gain_split) = {'lpsp>th hi'}; DT.grp(m & ~(DT.gain > gain_split)) = {'lpsp>th lo'}; end
DT.os = DT.exc >= overshoot_deg; DT.anti = -sign(DT.dfly).*DT.dbump - abs(DT.dfly) >= overshoot_deg; DT.null_os = DT.null_exc >= overshoot_deg;
DT.null_os(isnan(DT.null_exc)) = false; has_null = ~isnan(DT.null_exc);
low_bias = DT.bias < bias_max_deg | (isnan(DT.bias) & DT.fwd < fwd_max); % NaN bias = too little walking to estimate the leak
usable = DT.rho >= rho_min & DT.rho_p10 >= rho_p10_min & low_bias & ~DT.jumpy;
DK = usable & DT.dark==1; CL = usable & DT.dark==0;
fprintf('usable turns (mounting bias < %g deg, mean rho >= %g, rho p10 >= %g, no jumps): %d dark, %d closed loop\n',bias_max_deg,rho_min,rho_p10_min,sum(DK),sum(CL));
fprintf('heading added by the mounting bias over the scored interval, usable dark turns: median %.1f deg, 95th pct %.1f deg (all dark turns: %.1f / %.1f)\n', ...
    median(DT.bias(DK),'omitnan'),prctile(DT.bias(DK),95),median(DT.bias(DT.dark==1),'omitnan'),prctile(DT.bias(DT.dark==1),95));

%% rates per group, controls, per-fly stats
rate = nan(numel(grps),2,3); % group x (dark, CL) x (overshoot, opposite, still baseline)
V = cell(numel(grps),1);
for gi = 1:numel(grps)
    for d = 1:2
        mm = strcmp(DT.grp,grps{gi}) & (d==1 & DK | d==2 & CL);
        rate(gi,d,:) = [mean(DT.os(mm)), mean(DT.anti(mm)), mean(DT.null_os(mm & has_null))];
    end
    mm = DK & strcmp(DT.grp,grps{gi}); fl = unique(DT.fly(mm));
    nt = arrayfun(@(f) sum(mm & DT.fly==f),fl); r1 = arrayfun(@(f) mean(DT.os(mm & DT.fly==f)),fl);
    V{gi} = r1(nt >= min_fly_turns);
    fprintf('%-18s dark: overshoot %.1f%%, opposite %.1f%%, still baseline %.1f%% (%d turns); CL %.1f%%; %d flies, %d with >= 1 dark lap\n', ...
        grp_names{gi},100*rate(gi,1,1),100*rate(gi,1,2),100*rate(gi,1,3),sum(mm),100*rate(gi,2,1),numel(V{gi}),sum(V{gi} > 0));
end
p = arrayfun(@(c) ranksum(V{cmp(c,1)},V{cmp(c,2)}),1:size(cmp,1))';
for c = 1:size(cmp,1), fprintf('  %s vs %s: per-fly dark lap rate ranksum p = %.3g\n',grp_names{cmp(c,1)},grp_names{cmp(c,2)},p(c)); end
mh = strcmp(DT.grp,main_grp); fl = intersect(unique(DT.fly(DK & mh)),unique(DT.fly(CL & mh)));
pr = [arrayfun(@(f) mean(DT.os(DK & mh & DT.fly==f)),fl), arrayfun(@(f) mean(DT.os(CL & mh & DT.fly==f)),fl)];
p_dc = signrank(pr(:,1),pr(:,2));
fprintf([main_lbl ' with both dark and CL usable turns (n=%d): dark %.1f%% vs CL %.1f%%, signrank p = %.3g\n'],numel(fl),100*mean(pr(:,1)),100*mean(pr(:,2)),p_dc);

%% figure
figure('Color','w','Position',[40,40,1750,1000]); clf
tl = tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
% (a) mounting bias contamination
nexttile; hold on
edges = 0:5:100;
histogram(min(DT.bias(DT.dark==1),100),edges,'Normalization','probability','FaceColor',[.7,.7,.7],'EdgeColor','none')
histogram(min(DT.bias(DK),100),edges,'Normalization','probability','FaceColor',grp_col(1,:),'FaceAlpha',.6,'EdgeColor','none')
xline(overshoot_deg,'--k'); xlabel('heading added by mounting bias over a turn (deg)'); ylabel('fraction of dark turns')
xline(bias_max_deg,':k','LineWidth',1.2)
legend('all dark turns','usable','lap threshold','bias cutoff','Location','northeast')
title({sprintf('(a) turns kept: mounting bias < %g%s',bias_max_deg,char(176)),sprintf('%d of %d dark turns (last bin = >= 100)',sum(DK),sum(DT.dark==1 & DT.rho >= rho_min & DT.rho_p10 >= rho_p10_min & ~DT.jumpy))},'FontWeight','normal','FontSize',9)
% (b) per-fly dark lap rate
nexttile; hold on
for gi = 1:numel(grps)
    y = 100*V{gi}; scatter(gi + (rand(size(y))-.5)*.35,y,22,grp_col(gi,:),'filled','MarkerFaceAlpha',.55)
    errorbar(gi,mean(y),std(y)/sqrt(numel(y)),'o','Color','k','MarkerFaceColor',grp_col(gi,:),'MarkerSize',8,'LineWidth',1.2)
end
ylim([0,max(cellfun(@(v) max([v;0]),V))*100*1.05]); add_p(gca,1:n_g,cmp,p)
xticks(1:n_g); xticklabels(arrayfun(@(gi) sprintf('%s (n=%d)',grp_names{gi},numel(V{gi})),1:n_g,'UniformOutput',false)); xtickangle(30); xlim([0.4,n_g+0.6])
ylabel('dark laps per usable turn (%, per fly)'); title({['(b) ' main_lbl ' in the dark'],sprintf('flies with >= %d usable dark turns; ranksum',min_fly_turns)},'FontWeight','normal','FontSize',9)
% (c) direction / drift controls
nexttile; hold on
hb = bar(1:n_g,100*squeeze(rate(:,1,:)),0.85,'EdgeColor','none');
hb(1).FaceColor = [0.60,0.15,0.00]; hb(2).FaceColor = [.55,.55,.55]; hb(3).FaceColor = [.85,.85,.85];
for gi = 1:n_g, text(gi+0.27,100*rate(gi,1,3)+0.25,sprintf('%.1f',100*rate(gi,1,3)),'FontSize',7,'HorizontalAlignment','center'), end
xticks(1:n_g); xticklabels(grp_names); xtickangle(30); ylabel('% of usable dark turns')
legend('overshoot in the turn direction','same size, opposite direction','still-window drift baseline','Location','northeast')
title({'(c) dark overshoots follow the turn:','not explained by symmetric or spontaneous bump drift'},'FontWeight','normal','FontSize',9)
% (d) dark vs closed loop, high-gain flies
nexttile; hold on
plot([1,2],100*pr','-','Color',[grp_col(1,:),.35]); scatter(ones(size(pr,1),1),100*pr(:,1),25,grp_col(1,:),'filled'); scatter(2*ones(size(pr,1),1),100*pr(:,2),25,grp_col(1,:),'filled')
errorbar([1,2],100*mean(pr),100*std(pr)/sqrt(size(pr,1)),'ko-','MarkerFaceColor','k','LineWidth',1.5)
xticks([1,2]); xticklabels({'dark','closed loop'}); xlim([0.6,2.4]); ylabel('laps per usable turn (%, per fly)')
title({sprintf(['(d) ' main_lbl ' (n=%d)'],size(pr,1)),sprintf('dark vs closed loop, signrank p = %.2g',p_dc)},'FontWeight','normal','FontSize',9)
% (e) examples: largest dark overshoots in high-gain flies, one per fly
E = DT(DK & mh & DT.os,:); E = sortrows(E,'exc','descend'); [~,first] = unique(E.fly,'stable'); E = E(first(1:min(n_examples,numel(first))),:);
for k = 1:height(E)
    tr = all_data(E.trial(k)); xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); mu = cart2pol(mean(x_,2),mean(y_,2));
    mu_u = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'));
    fw = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',max(1,round(0.5/median(diff(xf)))));
    t0 = E.t0(k); win = [t0-pad_s, t0+E.dur(k)+post_s+pad_s]; ib = xb >= win(1) & xb <= win(2); jf = xf >= win(1) & xf <= win(2);
    nexttile(4+k); hold on
    n_half = size(f_z,1)/2; imagesc(xb(ib)-t0,1:2*n_half,f_z(:,ib)); clim([-1.5,3]); colormap(gca,flipud(gray))
    al = tr.im.alpha(1:n_half); dgl = 2*pi/n_half; b0 = find(xb >= t0,1); f0 = find(xf >= t0,1);
    hh = interp1(xf,hd,xb) - hd(f0) + mu(b0);                 % heading drawn from the bump's onset position
    for h = 1:2
        rb = mod(mu(ib) - al(1) + dgl/2,2*pi)/dgl + 0.5; rh = mod(hh(ib) - al(1) + dgl/2,2*pi)/dgl + 0.5;
        rb([false; abs(diff(rb)) > n_half/2]) = NaN; rh([false; abs(diff(rh)) > n_half/2]) = NaN;
        plot(xb(ib)-t0,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',2); plot(xb(ib)-t0,rh+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',2)
    end
    yline(n_half+0.5,'k'); xline([0,E.dur(k)],'k:'); axis tight; xlim(win-t0); set(gca,'YDir','reverse'); yticks([])
    tok = regexp(tr.meta,'\\(\d{8})\\fly\s*(\d+)[^\\]*\\\d{8}-(\d+)_','tokens','once');
    title(sprintf('(e) dark: %s fly %s trial %s (gain %.2f)',tok{:},E.gain(k)),'FontWeight','normal','FontSize',8)
    if k == 1, ylabel('PB glomeruli'), end
    nexttile(8+k); hold on
    patch([0,E.dur(k),E.dur(k),0],[-1e4,-1e4,1e4,1e4],[.93,.93,.93],'EdgeColor','none','HandleVisibility','off')
    plot(xf(jf)-t0,rad2deg(hd(jf)-hd(f0)),'Color',[0,0.45,0.74],'LineWidth',1.5,'DisplayName','fly heading')
    plot(xb(ib)-t0,rad2deg(mu_u(ib)-mu_u(b0)),'Color',[0.85,0.2,0.1],'LineWidth',1.5,'DisplayName','bump')
    yyaxis right; plot(xf(jf)-t0,fw(jf),'Color',[.6,.6,.6],'HandleVisibility','off'); ylabel('forward speed'); set(gca,'YColor',[.5,.5,.5]); yyaxis left; set(gca,'YColor','k')
    yl = [min([rad2deg(hd(jf)-hd(f0)); rad2deg(mu_u(ib)-mu_u(b0))]), max([rad2deg(hd(jf)-hd(f0)); rad2deg(mu_u(ib)-mu_u(b0))])]; ylim(yl + [-1,1]*0.1*diff(yl))
    xlim(win-t0); xlabel('time from turn onset (s)'); if k == 1, ylabel('rotation (deg)'), legend('Location','best'), end
    title(sprintf('bump %+.0f%s vs fly %+.0f%s (grey = turn)',E.dbump(k),char(176),E.dfly(k),char(176)),'FontWeight','normal','FontSize',8)
end
title(tl,{sprintf('Lapping in the dark: turns with mounting bias < %g%s, lap = bump >= %g%s past the fly within %g s of the turn''s end',bias_max_deg,char(176),overshoot_deg,char(176),post_s), ...
    sprintf('bump strength >= %g throughout (10th pct) and no raw bump jump > %g%s;',rho_p10_min,jump_max_deg,char(176)), '(e) heatmaps: bump (red), fly heading started at the bump (blue); traces: grey right axis = forward speed'},'FontSize',11)
exportgraphics(gcf,fullfile(out_dir,['fig_rnai_th_lapping_dark' tag '.png']),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,['fig_rnai_th_lapping_dark' tag '.pdf']),'ContentType','image','Resolution',300)

%% functions
function add_p(ax,x,cmp,p)
    % significance brackets above the data, stacked
    yl = ylim(ax); top = yl(2); h = 0.06*diff(yl);
    for c = 1:size(cmp,1)
        y = top + (c-0.4)*h*1.6;
        plot(ax,x(cmp(c,[1,1,2,2])),[y-h*0.4,y,y,y-h*0.4],'k','LineWidth',0.8,'HandleVisibility','off')
        if p(c) < 0.001, ps = sprintf('p = %.1e',p(c)); else, ps = sprintf('p = %.3f',p(c)); end
        text(ax,mean(x(cmp(c,:))),y,ps,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    end
    ylim(ax,[yl(1), top + (size(cmp,1)+0.6)*h*1.6])
end
