%% fig_rnai_vglut_stall_examples
% "Stalls": nearly the opposite of lapping (fig_rnai_th_lapping_examples.m)
% -- in closed loop the fly turns the bar (cue) a long way (>= stall_cue_min_deg,
% e.g. 180 deg) but the bump stays stuck, following only a small fraction.
% The reference is the cue, not FicTrac heading: when the bar doesn't move
% (cue not following the fly) a still bump is correct, not a stall. Only
% trials passing the integrative-gain inclusion (walking + working cue +
% bump) are scanned. Searched
% in closed-loop trials of the VGlut set (LPsP > VGlut-RNAi, Empty >
% VGlut-RNAi, LPsP > mCherry); the TH genotypes are scanned too, for
% reference. Same bump / heading estimates and scanning as the lap script.
%
% Detector, per trial, sliding windows of each length in stall_win_s (step
% stall_step_s); a window counts if
%   - the cue turns >= stall_cue_min_deg (unwrapped -ft.cue), and agrees
%     with the fly: same direction, cue / fly rotation within cue_gain_range
%     (closed-loop gain 0.8; rules out a stuck or free-running bar)
%   - the bump moves <= stall_bump_frac x the cue's turn, either direction
%   - the bump stays stuck THROUGHOUT the turn: over 1 s bins where the cue
%     turns >= bin_cue_min_deg (same direction; at least min_turn_frac of the
%     window's bins), median(bump change / cue change) <= stall_ratio_max
%   - the cue-bump gap grows steadily, not in one step (no step_win_s stretch
%     carries > step_frac_max of it; net >= straight_min of its path on a
%     0.5 s grid) -- rejects heading glitches and bump relocations
%   - the bump stays near where it started throughout: never more than
%     stall_follow_max x |cue| away (rejects a bump that tracks the cue for a
%     while, then snaps back -- poor tracking, not a stall)
%   - a solid bump: median PVA strength >= rho_min, 10th percentile >=
%     rho_p10_min, no raw frame-to-frame PVA jump > jump_max_deg -- so a lost /
%     weak bump isn't counted as stuck
%   overlapping passing windows, of any length -> one episode (its window
%   with the largest fly - bump gap).
% Figures: (1) one stall per LPsP > VGlut fly, (2) control flies' stalls,
% (3) summary: an example stall, stalls per 10 min per fly, fraction of
% flies with a stall, p-values (ranksum on per-fly rates, Fisher on flies
% with >= 1). Example panels show the cue (bar) position as the reference
% (closed loop: 0.8 x heading), as does detection.
% Needs all_data (loaded if absent) and the integrative-gain cache (trial
% labels, per-fly gain). Run from the repo root.

%% parameters + labels
stall_win_s = [5,20]; stall_step_s = 1;
stall_cue_min_deg = 180;    % cue (bar) rotation within the window
stall_bump_frac = 0.25;     % |bump| <= this x |cue|
stall_ratio_max = 0.25;     % median bump / cue change over turning 1 s bins
cue_gain_range = [0.5,1.2]; % cue / fly rotation in the window (closed-loop gain 0.8): bar and fly must agree
rho_min = 0.4; jump_max_deg = 90;
rho_p10_min = 0.2;     % 10th percentile of PVA strength in the window: no stretch of lost bump
stall_follow_max = 0.4; % the bump never gets further than this x |cue| from where it started (stuck throughout, not following then snapping back)
step_win_s = 0.3; step_frac_max = 0.35; straight_min = 0.4; % on the cue - bump gap
bin_cue_min_deg = 5; min_turn_frac = 0.2;
pad_s = 3; n_show = 8; gain_split = 0.85;
scan_geno = {'lpsp>vglut','empty>vglut','lpsp>mcherry','lpsp>th','empty>th'};
geno_set = 'vglut'; % 'vglut' or 'th': which set the example galleries and the summary show (all genotypes are always scanned)
switch geno_set
    case 'vglut'
        main_geno = 'lpsp>vglut'; main_name = 'LPsP > VGlut-RNAi'; ctl_geno = {'empty>vglut','lpsp>mcherry'}; set_tag = 'vglut';
        grp_def = {'LPsP > VGlut-RNAi','lpsp>vglut',[]; 'Empty > VGlut-RNAi','empty>vglut',[]; 'LPsP > mCherry','lpsp>mcherry',[]};
        grp_col = [0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13]; grp_x = 1:3; cmp = [1,2; 1,3];
    case 'th'
        main_geno = 'lpsp>th'; main_name = 'LPsP > TH-RNAi'; ctl_geno = {'empty>th','lpsp>mcherry'}; set_tag = 'th';
        grp_def = {'LPsP > TH-RNAi','lpsp>th',[]; 'Empty > TH-RNAi','empty>th',[]; 'LPsP > mCherry','lpsp>mcherry',[]; ...
                   'LPsP>TH high gain','lpsp>th',true; 'LPsP>TH low gain','lpsp>th',false};
        grp_col = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13; 0.60,0.15,0.00; 0.95,0.62,0.45];
        grp_x = [1,2,3,4.6,5.6]; cmp = [1,2; 1,3; 4,5];
end
summary_example = {}; % {date, fly, trial} of the stall shown in the summary; empty = the set's largest main-genotype stall
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
if ~exist('all_data','var')
    tmp = load(fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
g_fly = nan(max(IG.fly_num),1); % closed-loop integrative gain per fly
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
out_dir = fullfile('ugly_figures','rnai_turns','laps_stalls');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% detect stalls in all closed-loop trials
scan_trials = find(ismember(IG.genotype,scan_geno) & ~IG.is_dark & IG.inc); % closed loop, walking + working cue + bump (integrative-gain inclusion)
ev = zeros(0,9); % trial, t0, t1, dbump, dcue, gap (|cue| - |bump| along the turn), median rho, window length, dfly
for i = scan_trials(:)'
    [xb,mu_raw,mu_u,rho,hd] = trial_bump_heading(all_data(i));
    cb = interp1(all_data(i).ft.xf(:),unwrap_cue(all_data(i).ft.cue),xb); % cue on the imaging timebase
    dt = median(diff(xb)); st = max(1,round(stall_step_s/dt)); n1 = max(1,round(1/dt)); n_half_s = max(1,round(0.5/dt)); n_step = max(1,round(step_win_s/dt));
    jump = [0; abs(diff(mu_raw))] > deg2rad(jump_max_deg);
    hits = zeros(0,8);
    for w_s = stall_win_s
    nW = round(w_s/dt);
    for a = 1:st:numel(xb)-nW
        b = a + nW;
        db = rad2deg(mu_u(b)-mu_u(a)); dc = rad2deg(cb(b)-cb(a)); df = rad2deg(hd(b)-hd(a));
        if isnan(dc) || isnan(db) || abs(dc) < stall_cue_min_deg || abs(db) > stall_bump_frac*abs(dc), continue, end
        if isnan(df) || sign(dc) ~= sign(df) || abs(dc) < cue_gain_range(1)*abs(df) || abs(dc) > cue_gain_range(2)*abs(df), continue, end % bar moved with the fly
        if median(rho(a:b)) < rho_min || prctile(rho(a:b),10) < rho_p10_min || any(jump(a:b)), continue, end
        if max(abs(rad2deg(mu_u(a:b)-mu_u(a)))) > stall_follow_max*abs(dc), continue, end % bump stays near its start throughout
        % cue - bump gap, signed so a stall is the gap rising: it must grow steadily
        d = sign(dc)*rad2deg((cb(a:b)-cb(a)) - (mu_u(a:b)-mu_u(a)));
        if any(isnan(d)), continue, end
        gap = d(end);
        step_max = max(abs(d(1+n_step:end) - d(1:end-n_step)));
        if step_max > step_frac_max*gap || gap/sum(abs(diff(d(1:n_half_s:end)))) < straight_min, continue, end
        % bump stuck throughout the turn
        k_b = a:n1:b; dB1 = diff(rad2deg(mu_u(k_b))); dC1 = diff(rad2deg(cb(k_b)));
        turning = abs(dC1) >= bin_cue_min_deg & sign(dC1) == sign(dc);
        if sum(turning) < max(1,ceil(min_turn_frac*numel(dC1))) || median(dB1(turning)./dC1(turning)) > stall_ratio_max, continue, end
        hits(end+1,:) = [xb(a), xb(b), db, dc, gap, median(rho(a:b)), w_s, df]; %#ok<AGROW>
    end
    end
    while ~isempty(hits) % overlapping hits -> one episode, represented by its largest-gap window
        [~,k] = max(hits(:,5)); e = hits(k,:);
        ev(end+1,:) = [i, e]; %#ok<AGROW>
        hits(hits(:,1) < e(2) & hits(:,2) > e(1),:) = [];
    end
end
E = array2table(ev,'VariableNames',{'trial','t0','t1','dbump','dcue','gap','rho','win_s','dfly'});
E.geno = IG.genotype(E.trial); E.fly = IG.fly_num(E.trial); E.gain = g_fly(E.fly);
[E.date,E.fly_n,E.trial_n] = cellfun(@meta_labels,{all_data(E.trial).meta}','UniformOutput',false);
fprintf('\nstalls (cue >= %g deg, |bump| <= %g x cue, median rho >= %g, windows %s s; included trials only):\n',stall_cue_min_deg,stall_bump_frac,rho_min,mat2str(stall_win_s));
fprintf('%-13s %8s %10s %14s %16s\n','genotype','stalls','per 10 min','flies with >=1','stalls/fly');
fly_min = accumarray(IG.fly_num(scan_trials),arrayfun(@(i) all_data(i).ft.xf(end)-all_data(i).ft.xf(1),scan_trials)/60,[max(IG.fly_num),1]);
n_fly_ev = accumarray(E.fly,1,[max(IG.fly_num),1]);
for gi = 1:numel(scan_geno)
    fl = unique(IG.fly_num(scan_trials(strcmp(IG.genotype(scan_trials),scan_geno{gi}))));
    s = strcmp(E.geno,scan_geno{gi});
    fprintf('%-13s %8d %10.2f %8d of %3d %9.2f\n',scan_geno{gi},sum(s),10*sum(s)/sum(fly_min(fl)),sum(n_fly_ev(fl) > 0),numel(fl),mean(n_fly_ev(fl)));
end
E = sortrows(E,'gap','descend');

%% figure: one stall per main-genotype fly
Em = E(strcmp(E.geno,main_geno),:);
[~,first] = unique(Em.fly,'stable'); pick = Em(first(1:min(n_show,numel(first))),:);
if isempty(pick), disp(['no ' main_name ' stalls']), else
    plot_stalls(pick,all_data,gain_split,pad_s, ...
        {sprintf(['stalls in ' main_name ' (closed loop): cue (bar) turns >= %g%s within %s s windows, bump follows <= %g x, median \\rho >= %g; one per fly'],stall_cue_min_deg,char(176),mat2str(stall_win_s),stall_bump_frac,rho_min), ...
        'top: PB activity (z-scored) with bump (red) and cue / bar position (blue) on each half; bottom: unwrapped rotation of bump, cue and fly; shaded = detected window'})
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_examples.png',set_tag)),'Resolution',200)
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_examples.pdf',set_tag)),'ContentType','image','Resolution',300)
end

%% figure: control stalls, same detector, one per fly
Ec = E(ismember(E.geno,ctl_geno),:);
[~,first] = unique(Ec.fly,'stable'); pick_ctl = Ec(first(1:min(n_show,numel(first))),:);
if isempty(pick_ctl), disp('no control stalls'), else
    plot_stalls(pick_ctl,all_data,gain_split,pad_s, ...
        {sprintf('stalls in control flies (%s), same detector, one per fly, largest gap first',strjoin(ctl_geno,', ')), ...
        'top: PB activity (z-scored) with bump (red) and cue / bar position (blue) on each half; bottom: unwrapped rotation of bump, cue and fly; shaded = detected window'})
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_examples_controls.png',set_tag)),'Resolution',200)
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_examples_controls.pdf',set_tag)),'ContentType','image','Resolution',300)
end

%% figure: summary -- an example stall, stalls per fly, fraction of flies with a stall, stats
% groups (grp_def, grp_col, grp_x, cmp) come from geno_set in the parameters
G = struct('rate',{},'has',{},'n',{});
for gi = 1:size(grp_def,1)
    fl = unique(IG.fly_num(scan_trials(strcmp(IG.genotype(scan_trials),grp_def{gi,2}))));
    if ~isempty(grp_def{gi,3}), fl = fl(~isnan(g_fly(fl)) & (g_fly(fl) > gain_split) == grp_def{gi,3}); end
    G(gi).rate = 10*n_fly_ev(fl)./fly_min(fl); G(gi).has = n_fly_ev(fl) > 0; G(gi).n = numel(fl);
end
p_rate = nan(size(cmp,1),1); p_frac = p_rate;
for c = 1:size(cmp,1)
    a = G(cmp(c,1)); b = G(cmp(c,2));
    p_rate(c) = ranksum(a.rate,b.rate);
    [~,p_frac(c)] = fishertest([sum(a.has),sum(~a.has); sum(b.has),sum(~b.has)]);
    fprintf('%s vs %s: stalls / 10 min ranksum p = %.2g; flies with >= 1 stall Fisher p = %.2g\n',grp_def{cmp(c,1),1},grp_def{cmp(c,2),1},p_rate(c),p_frac(c));
end
k = find(strcmp(E.geno,main_geno));
if ~isempty(summary_example), k = find(strcmp(E.date,summary_example{1}) & strcmp(E.fly_n,summary_example{2}) & strcmp(E.trial_n,summary_example{3})); end
figure('Color','w','Position',[40,40,1500,760]); clf
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
if ~isempty(k)
    [~,j] = max(E.gap(k)); ex = E(k(j),:);
    t0 = ex.t0; win = [-pad_s, ex.t1-t0+pad_s]; S = ex_snippet(all_data(ex.trial),t0,win);
    nexttile(tl,1,[1,2]); hold on
    n_half = size(S.f_z,1)/2;
    imagesc(S.xb,1:2*n_half,S.f_z); clim([-1.5,3]); colormap(gca,flipud(gray))
    for h = 1:2
        [tb,rb] = angle_to_row(S.xb,S.bump_abs(:,h),S.alpha); [tc,rc] = angle_to_row(S.xf,S.cue_abs(:,h),S.alpha);
        plot(tb,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',3); plot(tc,rc+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',3)
    end
    yline(n_half+0.5,'k'); xline([0,ex.t1-t0],'k:'); xlim(win); ylim([0.5,2*n_half+0.5]); yticks([]); set(gca,'YDir','reverse')
    ylabel('PB glomeruli'); title(sprintf('(a) example stall: %s, %s fly %s trial %s (gain %.2f)',ex.geno{1},ex.date{1},ex.fly_n{1},ex.trial_n{1},ex.gain),'FontWeight','normal')
    nexttile(tl,5,[1,2]); hold on
    patch([0,ex.t1-t0,ex.t1-t0,0],[-1e4,-1e4,1e4,1e4],[.93,.93,.93],'EdgeColor','none','HandleVisibility','off')
    plot(S.xf,S.fly,'Color',[.6,.6,.6],'LineWidth',1,'DisplayName','fly heading')
    plot(S.xf,S.cue,'Color',[0,0.45,0.74],'LineWidth',1.5,'DisplayName','cue (bar) position')
    plot(S.xb,S.bump,'Color',[0.85,0.2,0.1],'LineWidth',1.5,'DisplayName','bump (PVA)')
    yl = [min([S.fly;S.cue;S.bump]),max([S.fly;S.cue;S.bump])]; ylim(yl + [-1,1]*0.08*diff(yl)); xlim(win)
    xlabel(sprintf('time from %.0f s (s)',t0)); ylabel('rotation (deg)'); legend('Location','best')
    title(sprintf('bump %+.0f%s vs cue %+.0f%s (fly %+.0f%s) in the detected window (shaded)',ex.dbump,char(176),ex.dcue,char(176),ex.dfly,char(176)),'FontWeight','normal')
end
nexttile(tl,3,[2,1]); hold on
for gi = 1:numel(G)
    y = G(gi).rate;
    scatter(grp_x(gi) + (rand(size(y))-.5)*.4,y,22,grp_col(gi,:),'filled','MarkerFaceAlpha',.55)
    errorbar(grp_x(gi),mean(y),std(y)/sqrt(numel(y)),'o','Color','k','MarkerFaceColor',grp_col(gi,:),'MarkerSize',8,'LineWidth',1.2)
end
add_p(gca,grp_x,cmp,p_rate)
xticks(grp_x); xticklabels(arrayfun(@(gi) sprintf('%s (n=%d)',grp_def{gi,1},G(gi).n),1:numel(G),'UniformOutput',false)); xtickangle(35)
xlim([0.4,max(grp_x)+0.6]); ylabel('stalls per 10 min of closed loop (per fly)'); title('(b) stall rate, dots = flies, mean \pm SEM','FontWeight','normal')
nexttile(tl,4,[2,1]); hold on
for gi = 1:numel(G)
    bar(grp_x(gi),mean(G(gi).has),0.7,'FaceColor',grp_col(gi,:),'EdgeColor','none')
    text(grp_x(gi),mean(G(gi).has)+0.02,sprintf('%d/%d',sum(G(gi).has),G(gi).n),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9)
end
ylim([0,1]); add_p(gca,grp_x,cmp,p_frac); yticks(0:0.25:1)
xticks(grp_x); xticklabels(grp_def(:,1)'); xtickangle(35)
xlim([0.4,max(grp_x)+0.6]); ylabel('fraction of flies with >= 1 stall'); title('(c) flies with a stall','FontWeight','normal')
title(tl,{'stalls: the bar turns but the bump stays stuck (closed loop)', ...
    sprintf('stall: in a %s s window the cue (bar) turns >= %g%s and the bump follows <= %g x, throughout the turn, steadily (no jumps), median \\rho >= %g', ...
    strjoin(arrayfun(@num2str,stall_win_s,'UniformOutput',false),' or '),stall_cue_min_deg,char(176),stall_bump_frac,rho_min), ...
    'p: ranksum on per-fly stall rates (b), Fisher exact on flies with >= 1 stall (c)'},'FontSize',11)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_summary.png',set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_stall_summary.pdf',set_tag)),'ContentType','image','Resolution',300)

%% functions
function plot_stalls(pick,all_data,gain_split,pad_s,ttl_lines)
    % one column per episode, 4 per row: heatmap with bump / cue overlay, then unwrapped rotation
    nc = 4; nr = ceil(height(pick)/nc);
    figure('Color','w','Position',[40,40,430*nc,400*nr]); clf
    tl = tiledlayout(2*nr,nc,'TileSpacing','compact','Padding','compact');
    for k = 1:height(pick)
        i = pick.trial(k); t0 = pick.t0(k); win = [-pad_s, pick.t1(k)-t0+pad_s];
        S = ex_snippet(all_data(i),t0,win);
        r = floor((k-1)/nc); c = mod(k-1,nc)+1;
        nexttile(tl,2*r*nc+c); hold on
        n_half = size(S.f_z,1)/2;
        imagesc(S.xb,1:2*n_half,S.f_z); clim([-1.5,3]); colormap(gca,flipud(gray))
        for h = 1:2
            [tb,rb] = angle_to_row(S.xb,S.bump_abs(:,h),S.alpha); [tc,rc] = angle_to_row(S.xf,S.cue_abs(:,h),S.alpha);
            plot(tb,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',2); plot(tc,rc+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',2)
        end
        yline(n_half+0.5,'k'); xline([0,pick.t1(k)-t0],'k:')
        xlim(win); ylim([0.5,2*n_half+0.5]); yticks([]); set(gca,'YDir','reverse')
        if isnan(pick.gain(k)), g_lbl = 'gain n/a'; elseif pick.gain(k) > gain_split, g_lbl = sprintf('gain %.2f, high',pick.gain(k)); else, g_lbl = sprintf('gain %.2f, low',pick.gain(k)); end
        title(sprintf('%s | %s fly %s trial %s (%s)',pick.geno{k},pick.date{k},pick.fly_n{k},pick.trial_n{k},g_lbl),'FontSize',8,'FontWeight','normal','Interpreter','none')
        if c == 1, ylabel('PB glomeruli'), end
        nexttile(tl,(2*r+1)*nc+c); hold on
        patch([0,pick.t1(k)-t0,pick.t1(k)-t0,0],[-1e4,-1e4,1e4,1e4],[.93,.93,.93],'EdgeColor','none','HandleVisibility','off')
        plot(S.xf,S.fly,'Color',[.6,.6,.6],'LineWidth',1,'DisplayName','fly heading')
        plot(S.xf,S.cue,'Color',[0,0.45,0.74],'LineWidth',1.5,'DisplayName','cue (bar)')
        plot(S.xb,S.bump,'Color',[0.85,0.2,0.1],'LineWidth',1.5,'DisplayName','bump (PVA)')
        yl = [min([S.fly;S.cue;S.bump]),max([S.fly;S.cue;S.bump])]; ylim(yl + [-1,1]*0.08*diff(yl)); xlim(win)
        title(sprintf('window: bump %+.0f%s, cue %+.0f%s (fly %+.0f%s)',pick.dbump(k),char(176),pick.dcue(k),char(176),pick.dfly(k),char(176)),'FontSize',9,'FontWeight','normal')
        xlabel(sprintf('time from %.0f s (s)',t0)); if c == 1, ylabel('rotation (deg)'), end
        if k == 1, legend('Location','best'), end
    end
    title(tl,ttl_lines)
end

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

function [xb,mu_raw,mu_u,rho,hd] = trial_bump_heading(tr)
    % PVA (raw unwrapped and 0.2 s smoothed) and unwrapped heading on the imaging timebase
    xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    mu_raw = unwrap(mu);
    mu_u = smoothdata(mu_raw,'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = interp1(xf,unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest')),xb);
end

function [d,f,t] = meta_labels(meta)
    % date, fly number and trial number from the ScanImage path
    tok = regexp(meta,'\\(\d{8})\\fly\s*(\d+)[^\\]*\\\d{8}-(\d+)_','tokens','once'); % fly folder may carry a suffix ('fly 6 _ empty _ vglut')
    if isempty(tok), d = ''; f = ''; t = ''; else, [d,f,t] = deal(tok{:}); end
end

function S = ex_snippet(tr,t0,win_s)
    % as fig_rnai_turns_examples.m: z-scored PB activity, bump and heading in
    % t0+win_s, times relative to t0, angles in deg relative to t0, plus
    % per-PB-half absolute PVA / heading for the heatmap overlay
    S = [];
    xf = tr.ft.xf(:);
    if isfield(tr.ft,'xb') && numel(tr.ft.xb) == size(tr.im.f,2), xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),size(tr.im.f,2))'; end
    win_s = [max(win_s(1),max(xf(1),xb(1))-t0), min(win_s(2),min(xf(end),xb(end))-t0)];
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2)) ./ std(f_s,0,2);
    [x_tmp,y_tmp] = pol2cart(tr.im.alpha(:)',f_z');
    [mu,rho] = cart2pol(mean(x_tmp,2),mean(y_tmp,2)); mu = unwrap(mu);
    mu = smoothdata(mu,'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'); hd = unwrap(hd);
    in_im = xb >= t0+win_s(1) & xb <= t0+win_s(2);
    in_ft = xf >= t0+win_s(1) & xf <= t0+win_s(2);
    S.xb = xb(in_im)-t0; S.f_z = f_z(:,in_im); S.bump = rad2deg(mu(in_im)-interp1(xb,mu,t0)); S.rho = rho(in_im);
    S.xf = xf(in_ft)-t0; S.fly = rad2deg(hd(in_ft)-interp1(xf,hd,t0));
    cu = unwrap_cue(tr.ft.cue); % bar position (closed loop: 0.8 x heading)
    S.cue = rad2deg(cu(in_ft)-interp1(xf,cu,t0));
    n_half = size(tr.im.f,1)/2;
    S.alpha = tr.im.alpha(1:n_half); S.alpha = S.alpha(:)';
    S.bump_abs = nan(sum(in_im),2); S.fly_abs = nan(sum(in_ft),2); S.cue_abs = S.fly_abs;
    for h = 1:2
        rows = (h-1)*n_half + (1:n_half);
        [x_h,y_h] = pol2cart(tr.im.alpha(rows),f_z(rows,:)');
        mu_h = smoothdata(unwrap(cart2pol(mean(x_h,2),mean(y_h,2))),'gaussian',max(1,round(0.2/median(diff(xb)))));
        S.bump_abs(:,h) = mu_h(in_im); S.fly_abs(:,h) = interp1(xb,mu_h,t0) + deg2rad(S.fly);
        S.cue_abs(:,h) = interp1(xb,mu_h,t0) + deg2rad(S.cue);
    end
end

function [t,row] = angle_to_row(t,ang,alpha)
    % PVA angle -> fractional row within one PB half (as fig_rnai_turns_examples.m)
    n = numel(alpha); d = 2*pi/n;
    row = mod(ang(:) - alpha(1) + d/2, 2*pi)/d + 0.5;
    jump = [false; abs(diff(row)) > n/2];
    t = t(:); k = find(jump);
    for j = numel(k):-1:1
        t = [t(1:k(j)-1); NaN; t(k(j):end)]; row = [row(1:k(j)-1); NaN; row(k(j):end)];
    end
end

function cu = unwrap_cue(cue)
    % -ft.cue unwrapped over the valid samples, THEN gaps filled (filling first
    % interpolates across the +-pi wrap the long way round: a fake 360 deg jump)
    cu = -cue(:); ok = ~isnan(cu);
    cu(ok) = unwrap(cu(ok)); cu = fillmissing(cu,'linear','EndValues','nearest');
end
