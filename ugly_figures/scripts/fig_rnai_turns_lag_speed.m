%% fig_rnai_turns_lag_speed
% Two follow-ups to fig_rnai_turns_quiet_cl.m, same turns (closed loop,
% quiet baseline, from the lpsp_rnai_turns_claude.m turn cache), same three
% genotypes (lpsp>th, empty>th, lpsp>mcherry):
%
%   1) LAG between bump and fly heading, per turn-size bin (15..90 deg).
%      Per fly x bin, the mean rotational-velocity traces of fly and bump
%      (fly_vel_trace / bump_vel_trace: both smoothed identically, 0.5 s
%      Gaussian, so the smoothing itself adds no relative lag) are
%      cross-correlated; the lag of the correlation peak is that fly's
%      lag for the bin (positive = bump lags the fly). Velocity rather than
%      displacement because displacement traces are monotone ramps whose
%      cross-correlation peak is broad and biased toward the window edges.
%      NB the bump trace is NOT shifted by bump_lag_s (0.13 s) in the cache
%      -- the lag measured here is the full behavior->bump lag.
%
%   2) TURNS BINNED BY SPEED instead of size: mean angular speed over the
%      turn (size / duration, deg/s), turns of 7.5..97.5 deg only (same
%      size range as the size bins). (a) fly heading vs. bump overlay per
%      speed bin, rows = speed, columns = genotype (same layout as the
%      size-binned overlay); (b) bump vs. fly displacement over the turn
%      per speed bin (dbump_deg / dfly_deg from the cache: displacement from
%      onset to offset, bump measured bump_lag_s later), one dot per fly
%      (mean over its turns in the bin), genotype means +/- SEM; plus the
%      peak within 1.5 s after turn offset, which lets the lagging bump
%      finish moving.
%      Speed and size are correlated (r ~ 0.6 over these turns), so the
%      fast bins also hold bigger turns on average.
%
% Run from the repo root.

%% load turn cache
cache_file = fullfile('.data','lpsp_rnai_turns_cache.mat');
load(cache_file,'turns','turn_bin','turn_quiet','turn_fly','turn_dark','bin_labels','n_bins', ...
    'turn_t','turn_win_s','turn_bin_edges_deg','min_turns_per_bin','fly_geno','n_flies','geno_keep','geno_colors');

plot_geno  = {'lpsp>th','empty>th','lpsp>mcherry'};
plot_label = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'};
n_g = numel(plot_geno);
gcol = zeros(n_g,3); for gi = 1:n_g, gcol(gi,:) = geno_colors(strcmp(geno_keep,plot_geno{gi}),:); end
size_bins = find(turn_bin_edges_deg(2:end) <= 97.5); % 15..90 deg
sel = turn_dark==0 & turn_quiet;
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% 1) lag per fly x turn-size bin
max_lag_s = 1.0; % search window, +/- s
dt = median(diff(turn_t));
lags = round(-max_lag_s/dt):round(max_lag_s/dt);
lag_fb = nan(n_flies,n_bins); r_fb = nan(n_flies,n_bins);
for f = unique(turn_fly(sel))'
    for b = size_bins
        idx = find(sel & turn_fly==f & turn_bin==b);
        if numel(idx) < min_turns_per_bin, continue, end
        fv = mean(cat(2,turns(idx).fly_vel_trace),2,'omitnan');
        bv = mean(cat(2,turns(idx).bump_vel_trace),2,'omitnan');
        [lag_fb(f,b),r_fb(f,b)] = peak_lag(fv,bv,lags,dt);
    end
end

figure('Color','w','Position',[60,60,1100,420]); clf
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
panels = {lag_fb, 'bump lag behind fly (s)', 'lag of peak velocity cross-correlation';
          r_fb,   'peak correlation',        'correlation at that lag'};
for k = 1:2
    nexttile; hold on
    groupbins(panels{k,1}(:,size_bins),fly_geno,plot_geno,plot_label,gcol)
    xticks(1:numel(size_bins)); xticklabels(bin_labels(size_bins)); xlim([0.5,numel(size_bins)+0.5])
    xlabel('turn size'); ylabel(panels{k,2}); title(panels{k,3})
    if k == 1, yline(0,':k','HandleVisibility','off'); legend('Location','northeast'), end
end
title(tl,{'closed loop, quiet-baseline turns: bump vs. fly rotational velocity, per-fly mean traces cross-correlated', ...
    sprintf('faint = flies, thick = mean \\pm SEM; lag searched over \\pm%.1f s, positive = bump lags',max_lag_s)})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_lag_by_size.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_lag_by_size.png'),'Resolution',200)

for gi = 1:n_g
    g = strcmp(fly_geno,plot_geno{gi});
    fprintf('%-13s lag (s) by size bin: %s\n',plot_geno{gi},strjoin(arrayfun(@(b) sprintf('%s %.2f',bin_labels{b},mean(lag_fb(g,b),'omitnan')),size_bins,'UniformOutput',false),', '));
end

%% 2) bin turns by speed
speed_edges = [20, 35, 50, 65, 80, 100, 150]; % deg/s, mean angular speed over the turn
turn_size  = [turns.size_deg]';
turn_speed = turn_size ./ [turns.dur_s]';
spd_bin = discretize(turn_speed,speed_edges);
spd_bin(turn_size < turn_bin_edges_deg(1) | turn_size > turn_bin_edges_deg(size_bins(end)+1)) = NaN;
n_sb = numel(speed_edges)-1;
spd_labels = arrayfun(@(b) sprintf('%g-%g%s/s',speed_edges(b),speed_edges(b+1),char(176)),1:n_sb,'UniformOutput',false);
for b = 1:n_sb
    fprintf('%s: %d turns, median size %.0f deg\n',spd_labels{b},sum(sel & spd_bin==b),median(turn_size(sel & spd_bin==b)));
end

% per fly x speed bin: mean traces and mean displacements
FMs = nan(n_flies,n_sb,numel(turn_t)); BMs = FMs;
dfly_fs = nan(n_flies,n_sb); dbump_fs = nan(n_flies,n_sb); size_fs = nan(n_flies,n_sb);
for f = unique(turn_fly(sel))'
    for b = 1:n_sb
        idx = find(sel & turn_fly==f & spd_bin==b);
        if numel(idx) < min_turns_per_bin, continue, end
        FMs(f,b,:) = mean(cat(2,turns(idx).fly_trace),2,'omitnan');
        BMs(f,b,:) = mean(cat(2,turns(idx).bump_trace),2,'omitnan');
        dfly_fs(f,b)  = mean([turns(idx).dfly_deg],'omitnan');
        dbump_fs(f,b) = mean([turns(idx).dbump_deg],'omitnan');
        size_fs(f,b)  = mean(turn_size(idx));
    end
end

%% 2a) overlay traces by speed: rows = speed bin, columns = genotype
figure('Color','w','Position',[60,60,240*n_g,170*n_sb]); clf
tl = tiledlayout(n_sb,n_g,'TileSpacing','compact','Padding','compact');
ax = gobjects(n_g,n_sb);
for gi = 1:n_g
    flies_g = find(strcmp(fly_geno,plot_geno{gi}));
    for b = 1:n_sb
        ax(gi,b) = nexttile(tl,(b-1)*n_g+gi); hold on
        FT = reshape(FMs(flies_g,b,:),numel(flies_g),[]);
        BT = reshape(BMs(flies_g,b,:),numel(flies_g),[]);
        n_f = sum(~all(isnan(FT),2));
        if n_f >= 1
            plot_sem(turn_t,FT,[0,0,0])
            plot_sem(turn_t,BT,gcol(gi,:))
            plot(turn_t,0.8*mean(FT,1,'omitnan'),'--','Color',[.5,.5,.5])
        end
        xline(0,':k'); xlim(turn_win_s)
        if b == 1, title(plot_label{gi}), end
        text(0.03,0.97,sprintf('n=%d flies',n_f),'Units','normalized','VerticalAlignment','top','FontSize',8)
        if gi == 1, ylabel(sprintf('%s\nrotation (deg)',spd_labels{b})), end
        if b == n_sb, xlabel('time from turn onset (s)'), end
    end
end
for b = 1:n_sb, linkaxes(ax(:,b),'y'), end
title(tl,{'closed loop, quiet-baseline turns binned by mean speed: fly heading (black) vs. bump (color)','mean \pm SEM across flies; dashed = 0.8 \times fly'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_by_speed_overlay.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_by_speed_overlay.png'),'Resolution',200)

%% 2b) bump vs. fly displacement per speed bin -- two measures
% row 1: displacement over the turn itself (onset -> offset; bump 0.13 s
%        later), as cached. For fast, short turns the bump hasn't caught up
%        by then, so this folds the bump's lag into its apparent gain.
% row 2: peak displacement (from the onset value) within 0..peak_win_s
%        after turn OFFSET, which gives the bump time to finish moving.
%        Turn duration varies within a speed bin, so each turn's traces are
%        re-aligned to its own offset (interpolated onto 0..peak_win_s after
%        offset), averaged per fly, and the peak of that per-fly mean is
%        taken -- same peak-of-mean-trace logic as the other summaries. (A
%        per-turn max was tried first: the max of a single noisy trace over
%        1.5 s is strongly inflated, the noisier bump more so -- ~40 deg
%        bump "peaks" for ~13 deg turns in the slowest bin.) The cached
%        traces end at turn_win_s(2) = 3 s after onset, so turns longer
%        than 1.5 s get a truncated window (NaN beyond it; counted below).
peak_win_s = 1.5;
t_off = (0:median(diff(turn_t)):peak_win_s)';
pkF = nan(n_flies,n_sb); pkB = nan(n_flies,n_sb);
turn_dur = [turns.dur_s]';
for f = unique(turn_fly(sel))'
    for b = 1:n_sb
        idx = find(sel & turn_fly==f & spd_bin==b);
        if numel(idx) < min_turns_per_bin, continue, end
        F = nan(numel(t_off),numel(idx)); B = F;
        for k = 1:numel(idx)
            F(:,k) = interp1(turn_t,double(turns(idx(k)).fly_trace),turn_dur(idx(k))+t_off);
            B(:,k) = interp1(turn_t,double(turns(idx(k)).bump_trace),turn_dur(idx(k))+t_off);
        end
        pkF(f,b) = max(mean(F,2,'omitnan')); pkB(f,b) = max(mean(B,2,'omitnan'));
    end
end
for b = 1:n_sb
    in_b = sel & spd_bin==b;
    fprintf('%s: %.1f%% of turns have a post-offset window truncated by the 3 s trace end\n',spd_labels{b},100*mean(turn_dur(in_b) + peak_win_s > turn_win_s(2)));
end
meas = {dfly_fs, dbump_fs, 'over the turn (onset to offset)';
        pkF,     pkB,      sprintf('peak within %.1f s after offset',peak_win_s)};
figure('Color','w','Position',[60,60,250*n_sb,560]); clf
tl = tiledlayout(2,n_sb,'TileSpacing','compact','Padding','compact');
lim = [0, max([dfly_fs(:);dbump_fs(:);pkF(:);pkB(:)])*1.05];
h = gobjects(1,n_g);
for r = 1:2
    for b = 1:n_sb
        nexttile; hold on
        plot(lim,lim,':k','HandleVisibility','off'); plot(lim,0.8*lim,'--','Color',[.5,.5,.5],'HandleVisibility','off')
        for gi = 1:n_g
            g = strcmp(fly_geno,plot_geno{gi}); x = meas{r,1}(g,b); y = meas{r,2}(g,b); ok = ~isnan(x) & ~isnan(y);
            scatter(x(ok),y(ok),12,gcol(gi,:),'filled','MarkerFaceAlpha',.4,'HandleVisibility','off')
            h(gi) = errorbar(mean(x(ok)),mean(y(ok)),std(y(ok))/sqrt(sum(ok)),std(y(ok))/sqrt(sum(ok)), ...
                std(x(ok))/sqrt(sum(ok)),std(x(ok))/sqrt(sum(ok)),'o','Color','k','MarkerFaceColor',gcol(gi,:),'MarkerSize',7,'LineWidth',1);
        end
        axis square; xlim(lim); ylim([-20,lim(2)])
        if r == 1, title({spd_labels{b},sprintf('(median size %.0f%s)',median(turn_size(sel & spd_bin==b)),char(176))}), end
        xlabel('fly (deg)'); if b == 1, ylabel({meas{r,3},'bump (deg)'}), end
    end
end
legend(h,plot_label,'Location','northwest')
title(tl,{'closed loop, quiet-baseline turns binned by mean speed: bump vs. fly displacement', ...
    'dots = flies (mean over their turns in the bin), outlined = mean \pm SEM; dotted = 1, dashed = 0.8'})
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_by_speed_displacement.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_turns_by_speed_displacement.png'),'Resolution',200)

%% functions
function [lag_s,r_pk] = peak_lag(fv,bv,lags,dt)
    % lag (s) maximizing the Pearson correlation of bv(t) with fv(t - lag);
    % positive = bump lags. Only the overlapping samples count at each lag.
    r = nan(size(lags));
    for k = 1:numel(lags)
        L = lags(k);
        if L >= 0, a = fv(1:end-L); c = bv(1+L:end); else, a = fv(1-L:end); c = bv(1:end+L); end
        ok = ~isnan(a) & ~isnan(c);
        if sum(ok) > 10, r(k) = corr(a(ok),c(ok)); end
    end
    [r_pk,k] = max(r); lag_s = lags(k)*dt;
    if isnan(r_pk), lag_s = NaN; end
end

function groupbins(V,fly_geno,plot_geno,plot_label,gcol)
    % fig9-style: per-fly values (faint) + mean +/- SEM (thick), genotypes offset within each bin
    n_g = numel(plot_geno); n_b = size(V,2);
    for gi = 1:n_g
        Y = V(strcmp(fly_geno,plot_geno{gi}),:);
        xj = (1:n_b) + (gi-(n_g+1)/2)*(0.7/n_g);
        plot(repmat(xj,size(Y,1),1)',Y','.','Color',[gcol(gi,:),.3],'MarkerSize',8,'HandleVisibility','off')
        errorbar(xj,mean(Y,1,'omitnan'),std(Y,0,1,'omitnan')./sqrt(max(1,sum(~isnan(Y),1))),'o-','Color',gcol(gi,:), ...
            'MarkerFaceColor',gcol(gi,:),'LineWidth',2,'DisplayName',sprintf('%s (n=%d flies)',plot_label{gi},sum(any(~isnan(Y),2))))
    end
end

function plot_sem(t,X,c)
    % mean +/- SEM shading across rows of X, line in color c
    t = t(:)';
    m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m);
    if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end
