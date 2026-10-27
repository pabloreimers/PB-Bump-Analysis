%% fig_rnai_th_gain_bimodality
% LPsP > TH-RNAi flies split into two groups by closed-loop integrative
% gain (fig_rnai_integrative_gain.m): ~half sit at 0.3-0.8 like the
% controls, ~half at 0.9-1.7, with a gap near 0.85. This asks what, if
% anything, distinguishes the two groups.
%
%   Figure 1: (a) per-fly closed-loop gain, TH set, split line at
%             gain_split; (b) closed-loop vs. dark gain -- the split
%             persists without the visual cue; (c) per-trial gain of each
%             fly's closed-loop trials -- the group is a property of the
%             fly, not of one trial; (d) turn-aligned fly heading vs. bump
%             (turn cache, closed loop, all turns), low vs. high LPsP>TH
%             flies plus Empty>TH and LPsP>mCherry as reference: for matched fly turns the
%             high group's bump moves ~1x the fly, the low group ~0.5x.
%   Figure 2: per-fly walking / turning / bump / metadata variables vs.
%             closed-loop gain (LPsP>TH colored by group, controls grey),
%             Spearman rho within LPsP>TH and ranksum p low vs. high in
%             each title. These are ~20 uncorrected tests: treat p ~ 0.01-
%             0.05 as leads, not findings.
%
% Walking variables use the same closed-loop trials as the gain (walking +
% working cue + bump). Turn statistics come from the turn cache. Needs the
% integrative-gain cache (run fig_rnai_integrative_gain.m once), the turn
% cache, and all_data (loaded if absent). Run from the repo root.

%% parameters + load caches
gain_split = 0.85;
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
addpath(fullfile(pwd,'circ_stats'));
data_file = fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat');
load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','fly_list','inc','g_trial','xf_trial');
TC = load(fullfile(data_dir,'lpsp_rnai_turns_cache.mat'),'turns','turn_dark','turn_t','min_turns_per_bin');
if ~exist('all_data','var')
    tmp = load(data_file,'all_data'); all_data = tmp.all_data(:); clear tmp
end
set_geno  = {'lpsp>th','empty>th','lpsp>mcherry'};
set_names = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'};
col_lo = [0.95,0.62,0.45]; col_hi = [0.60,0.15,0.00]; col_ctl = [0.6,0.6,0.6];
geno_colors = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13];
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% per-fly gain (closed loop, dark) and groups
use = inc & ~cellfun(@isempty,g_trial);
n_flies = max(fly_num);
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = genotype{find(fly_num==f,1)}; end
g_fly = nan(n_flies,2);
for f = 1:n_flies
    for d = 0:1
        tr_i = find(use & fly_num==f & is_dark==d);
        if ~isempty(tr_i), g_fly(f,d+1) = mean(vertcat(g_trial{tr_i}),'omitnan'); end
    end
end
th  = find(strcmp(fly_geno,'lpsp>th') & ~isnan(g_fly(:,1)));
ctl = find(ismember(fly_geno,{'empty>th','lpsp>mcherry'}) & ~isnan(g_fly(:,1)));
is_hi = false(n_flies,1); is_hi(th) = g_fly(th,1) > gain_split;
lo_f = th(~is_hi(th)); hi_f = th(is_hi(th));
fprintf('LPsP>TH: %d low-gain, %d high-gain flies (split %.2f)\n',numel(lo_f),numel(hi_f),gain_split);

%% per-fly variables
turn_trial = [TC.turns.trial]'; turn_size = [TC.turns.size_deg]'; turn_dur = [TC.turns.dur_s]';
V = nan(n_flies,19);
for f = [th; ctl]'
    tr_i = find(use & fly_num==f & ~is_dark);
    fs = []; rs = []; hv = []; rho = []; fmean = []; width = []; mask_px = []; cue_g = []; t_len = 0;
    for i = tr_i'
        tr = all_data(i); xf = tr.ft.xf(:); h = tr.ft.heading(:);
        fs = [fs; tr.ft.f_speed(:)]; rs = [rs; tr.ft.r_speed(:)]; t_len = t_len + xf(end) - xf(1); %#ok<AGROW>
        w = round(20/median(diff(xf)));
        for s = 1:round(w/2):numel(h)-w, hv(end+1) = 1-abs(mean(exp(1i*h(s:s+w-1)),'omitnan')); end %#ok<AGROW>
        rho = [rho; tr.im.rho(:)]; fmean(end+1) = mean(tr.im.f(:),'omitnan'); mask_px(end+1) = sum(tr.im.mask(:)); %#ok<AGROW>
        n_half = size(tr.im.f,1)/2;
        for hh = 1:2
            fh = tr.im.f((hh-1)*n_half+(1:n_half),:); fh = (fh - min(fh,[],1))./(max(fh,[],1)-min(fh,[],1));
            width = [width, sum(fh > 0.5,1)]; %#ok<AGROW>
        end
        % cue vs heading: slope of 20 s displacements of unwrapped -cue on unwrapped heading
        c = -tr.ft.cue(:); ok = ~isnan(c) & ~isnan(h); cu = unwrap(c(ok)); hu = unwrap(h(ok));
        if numel(cu) > w+10, cue_g(end+1) = (hu(1+w:end)-hu(1:end-w)) \ (cu(1+w:end)-cu(1:end-w)); end %#ok<AGROW>
    end
    walk = fs > 0.1 | 2*abs(rs) > 0.1;
    ti = TC.turn_dark==0 & ismember(turn_trial,tr_i);
    tok = regexp(fly_list{f},'fly\s*(\d+)','tokens','once'); fly_n = str2double(tok{1});
    V(f,:) = [mean(fs,'omitnan'), mean(abs(rs(walk)),'omitnan'), mean(walk), mean(abs(rs)>2), ...
        abs(mean(rs,'omitnan'))/mean(abs(rs),'omitnan'), mean(abs(rs),'omitnan')/mean(abs(fs),'omitnan'), mean(hv,'omitnan'), ...
        sum(ti)/(t_len/60), median(turn_size(ti)), median(turn_size(ti)./turn_dur(ti)), mean(turn_size(ti) > 90), ...
        mean(rho,'omitnan'), mean(fmean), median(mask_px), mean(width,'omitnan'), median(cue_g,'omitnan'), ...
        fly_n, numel(tr_i), g_fly(f,2)];
end
feat_names = {'fwd speed (mean)','|rot speed| when walking (rad/s)','frac time walking','frac time |rot| > 2 rad/s', ...
    'turn bias |mean r| / mean |r|','|rot| / |fwd| speed','heading circ var (20 s windows)', ...
    'turns per min','median turn size (deg)','median turn speed (deg/s)','frac turns > 90 deg', ...
    'mean PVA strength (rho)','raw F (brightness)','PB mask size (px)','bump width (glom > half max)','cue / heading gain (20 s)', ...
    'fly # within day','n closed-loop trials','integrative gain, dark'};

%% figure 1: the split, dark gain, trial consistency, turn traces
figure('Color','w','Position',[60,60,1500,1150]); clf
tl = tiledlayout(4,4,'TileSpacing','compact','Padding','compact');

nexttile(tl,1); hold on
for gi = 1:3
    y = g_fly(strcmp(fly_geno,set_geno{gi}),1); y = y(~isnan(y));
    scatter(gi + (rand(size(y))-.5)*.35,y,22,geno_colors(gi,:),'filled','MarkerFaceAlpha',.5)
end
yline(gain_split,'--k'); yline(0.8,':','Color',[.5,.5,.5])
xticks(1:3); xticklabels(set_names); xtickangle(20); xlim([.4,3.6]); ylabel('integrative gain, closed loop')
title('(a) per-fly gain')

nexttile(tl,2); hold on
scatter(g_fly(ctl,1),g_fly(ctl,2),18,col_ctl,'filled','MarkerFaceAlpha',.5,'DisplayName','controls')
scatter(g_fly(lo_f,1),g_fly(lo_f,2),26,col_lo,'filled','DisplayName','LPsP>TH low')
scatter(g_fly(hi_f,1),g_fly(hi_f,2),26,col_hi,'filled','DisplayName','LPsP>TH high')
xline(gain_split,'--k','HandleVisibility','off'); plot([0,2],[0,2],':k','HandleVisibility','off')
ok = ~isnan(g_fly(th,2)); [r,p] = corr(g_fly(th(ok),1),g_fly(th(ok),2),'type','Spearman');
xlabel('gain, closed loop'); ylabel('gain, dark'); legend('Location','northwest')
title(sprintf('(b) closed loop vs dark (LPsP>TH \\rho=%.2f, p=%.2g)',r,p))

nexttile(tl,3); hold on
for f = th'
    tr_i = find(use & fly_num==f & ~is_dark);
    tg = cellfun(@(x) mean(x,'omitnan'),g_trial(tr_i));
    c = col_lo; if is_hi(f), c = col_hi; end
    plot(g_fly(f,1)*[1,1],[min(tg),max(tg)],'-','Color',c)
    scatter(g_fly(f,1)*ones(size(tg)),tg,16,c,'filled')
end
plot([0,2.2],[0,2.2],':k'); xline(gain_split,'--k'); yline(gain_split,'--k'); xlim([0,2]); ylim([0,2.2])
xlabel('fly gain (mean over trials)'); ylabel('single-trial gain'); title('(c) LPsP>TH: each closed-loop trial')

nexttile(tl,4); axis off
text(0,1,{'groups: LPsP>TH split at gain = 0.85', sprintf('low: n = %d flies, high: n = %d flies',numel(lo_f),numel(hi_f)), '', ...
    'dashed = split, dotted = unity / 0.8 cue gain'},'VerticalAlignment','top')

t = TC.turn_t(:); tf = fly_num(turn_trial);
bins = [30,60,90,135];
ref = find(strcmp(fly_geno,'empty>th') & ~isnan(g_fly(:,1))); mch = find(strcmp(fly_geno,'lpsp>mcherry') & ~isnan(g_fly(:,1)));
grp = {lo_f, hi_f, ref, mch}; grp_c = {col_lo, col_hi, geno_colors(2,:), geno_colors(3,:)}; grp_l = {'LPsP>TH low','LPsP>TH high','Empty>TH','LPsP>mCherry'};
grp_row = [2,2,3,4]; % figure row of each group's panel
ax = gobjects(numel(bins),3);
for b = 1:numel(bins)
    inb = TC.turn_dark==0 & abs(turn_size - bins(b)) < 7.5;
    for g = 1:numel(grp)
        FT = []; BT = [];
        for f = grp{g}'
            idx = find(inb & tf==f);
            if numel(idx) < TC.min_turns_per_bin, continue, end
            FT(end+1,:) = mean(cat(2,TC.turns(idx).fly_trace),2,'omitnan'); %#ok<AGROW>
            BT(end+1,:) = mean(cat(2,TC.turns(idx).bump_trace),2,'omitnan'); %#ok<AGROW>
        end
        r = grp_row(g);
        ax(b,r-1) = nexttile(tl,4*(r-1)+b); hold on
        fly_c = [0,0,0]; if g == 1, fly_c = [.5,.5,.5]; end
        plot_sem(t,FT,fly_c); plot_sem(t,BT,grp_c{g})
        text(0.03,0.97 - 0.08*(g==2),sprintf('%s n=%d',grp_l{g},size(FT,1)),'Units','normalized','VerticalAlignment','top','FontSize',8,'Color',grp_c{g})
        xline(0,':k'); xlim(t([1,end]))
        if r == 2, title(sprintf('(d) %d%s turns: LPsP>TH low vs high',bins(b),char(176))), else, title(sprintf('%d%s turns: %s',bins(b),char(176),grp_l{g})), end
        if r == 4, xlabel('time from turn onset (s)'), end
        if b == 1, ylabel('rotation (deg)'), end
    end
end
for b = 1:numel(bins), linkaxes(ax(b,:),'y'), end
title(tl,{'LPsP > TH-RNAi: high- vs low-gain flies (integrative gain, closed loop)', ...
    '(d) fly heading: black (grey = low group), color = bump; mean \pm SEM across flies'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_bimodality.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_bimodality.pdf'),'ContentType','vector')

%% figure 2: per-fly variables vs gain
n_v = numel(feat_names); nc = 5; nr = ceil(n_v/nc);
figure('Color','w','Position',[60,60,1600,300*nr]); clf
tl2 = tiledlayout(nr,nc,'TileSpacing','compact','Padding','compact');
fprintf('%-36s %8s %8s %9s %7s %8s\n','variable','low med','high med','ranksum p','rho_s','p_s');
for k = 1:n_v
    nexttile(tl2); hold on
    scatter(V(ctl,k),g_fly(ctl,1),14,col_ctl,'filled','MarkerFaceAlpha',.4)
    scatter(V(lo_f,k),g_fly(lo_f,1),22,col_lo,'filled')
    scatter(V(hi_f,k),g_fly(hi_f,1),22,col_hi,'filled')
    yline(gain_split,'--k')
    x = V(th,k); ok = ~isnan(x);
    [rs,ps] = corr(x(ok),g_fly(th(ok),1),'type','Spearman');
    pr = ranksum(V(lo_f(~isnan(V(lo_f,k))),k),V(hi_f(~isnan(V(hi_f,k))),k));
    fprintf('%-36s %8.3g %8.3g %9.3g %7.2f %8.3g\n',feat_names{k},median(V(lo_f,k),'omitnan'),median(V(hi_f,k),'omitnan'),pr,rs,ps);
    xlabel(feat_names{k}); if mod(k,nc) == 1, ylabel('gain, closed loop'), end
    title(sprintf('\\rho=%.2f (p=%.2g); low vs high p=%.2g',rs,ps,pr),'FontWeight','normal','FontSize',9)
end
title(tl2,{'per-fly variables vs closed-loop integrative gain','orange = LPsP>TH low, dark red = LPsP>TH high, grey = Empty>TH + LPsP>mCherry; stats within LPsP>TH only, uncorrected'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_bimodality_variables.png'),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_th_gain_bimodality_variables.pdf'),'ContentType','vector')

%% functions
function plot_sem(t,X,c)
    % mean +/- SEM shading across rows of X, line in color c
    t = t(:)';
    m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m);
    if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end
