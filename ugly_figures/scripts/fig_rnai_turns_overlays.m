%% fig_rnai_turns_overlays
% Group-mean fly heading vs. bump traces aligned to turn onset (closed
% loop, quiet-baseline or all turns, from the lpsp_rnai_turns_claude.m turn
% cache), for one genotype set, binned two ways:
%   (a) by TURN SIZE: 15 deg bins from 15 deg, as in fig_rnai_turns_quiet_cl.m,
%       extended past 90 deg while every genotype has > 7 flies per bin
%   (b) by TURN SPEED: mean angular speed over the turn (size / duration,
%       deg/s). No size cap here (fig_rnai_turns_lag_speed.m keeps only
%       7.5..97.5 deg turns), so long, slow turns that persist past the
%       size-binned range are included too.
% Rows = bins, columns = genotypes; mean +/- SEM across flies, each fly =
% mean of its >= min_turns_per_bin turns in the bin; black = fly, color =
% bump, dashed = 0.8 x fly (closed-loop cue gain). Shared y within a row.
% The cached traces only run turn_win_s = -1..3 s around onset.
%
% geno_set picks the genotypes ('vglut', 'th', or 'th_split' = TH set with
% LPsP>TH split into low- / high-gain flies at gain_split, by closed-loop
% integrative gain; needs the fig_rnai_integrative_gain.m cache), turn_set the turns
% ('quiet' = quiet baseline only, 'all' = every turn, counter-corrections
% of a preceding opposite turn included); exports are tagged with both.
% Run from the repo root.

%% parameters + load turn cache
geno_set = 'th_split'; % 'vglut', 'th' or 'th_split'
turn_set = 'all'; % 'quiet' or 'all'
switch geno_set
    case 'vglut'
        plot_geno  = {'lpsp>vglut','empty>vglut','lpsp>mcherry'};
        plot_label = {'LPsP > VGlut-RNAi','Empty > VGlut-RNAi','LPsP > mCherry'};
    case 'th'
        plot_geno  = {'lpsp>th','empty>th','lpsp>mcherry'};
        plot_label = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'};
    case 'th_split'
        % LPsP>TH split by closed-loop integrative gain (fig_rnai_th_gain_bimodality.m)
        plot_geno  = {'lpsp>th lo','lpsp>th hi','empty>th','lpsp>mcherry'};
        plot_label = {'LPsP > TH-RNAi, low gain','LPsP > TH-RNAi, high gain','Empty > TH-RNAi','LPsP > mCherry'};
end
gain_split = 0.85; % th_split only
speed_edges = [20, 35, 50, 65, 80, 100, 150]; % deg/s

cache_file = fullfile('data','lpsp_rnai_turns_cache.mat');
load(cache_file,'turns','turn_bin','turn_quiet','turn_fly','turn_dark','bin_labels','n_bins', ...
    'turn_t','turn_win_s','turn_bin_edges_deg','min_turns_per_bin','fly_geno','n_flies','geno_keep','geno_colors');
if strcmp(geno_set,'th_split')
    % per-fly closed-loop integrative gain from the fig_rnai_integrative_gain.m
    % cache (its own fly numbering, matched to the turn cache's flies via
    % trial index), then relabel LPsP>TH flies 'lpsp>th lo' / 'lpsp>th hi'
    IG = load(fullfile('data','lpsp_rnai_intgain_cache.mat'),'is_dark','fly_num','inc','g_trial');
    ig_use = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
    turn_trial = [turns.trial]';
    for f = find(strcmp(fly_geno,'lpsp>th'))'
        tr_f = unique(turn_trial(turn_fly==f));
        if isempty(tr_f), continue, end
        tr_i = find(ig_use & IG.fly_num==IG.fly_num(tr_f(1)));
        g = mean(vertcat(IG.g_trial{tr_i}),'omitnan');
        if isnan(g), fly_geno{f} = 'lpsp>th no gain'; elseif g > gain_split, fly_geno{f} = 'lpsp>th hi'; else, fly_geno{f} = 'lpsp>th lo'; end
    end
    geno_keep   = [geno_keep(:)', {'lpsp>th lo','lpsp>th hi'}];
    geno_colors = [geno_colors; 0.95,0.62,0.45; 0.60,0.15,0.00];
    fprintf('LPsP>TH: %d low-gain, %d high-gain flies (split %.2f)\n',sum(strcmp(fly_geno,'lpsp>th lo')),sum(strcmp(fly_geno,'lpsp>th hi')),gain_split);
end
n_g = numel(plot_geno);
gcol = zeros(n_g,3); for gi = 1:n_g, gcol(gi,:) = geno_colors(strcmp(geno_keep,plot_geno{gi}),:); end
sel = turn_dark==0;
if strcmp(turn_set,'quiet'), sel = sel & turn_quiet; set_lbl = 'quiet-baseline turns'; set_tag = '';
else, set_lbl = 'all turns'; set_tag = '_allturns'; end
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

turn_size  = [turns.size_deg]';
turn_speed = turn_size ./ [turns.dur_s]';
% size bins: 15 deg wide (centers 15, 30, ...; identical to the cache's
% turn_bin up to 90 deg). 15..90 always shown; beyond 90 deg, bins are
% added one at a time for as long as EVERY genotype in the set has more
% than min_flies_extend flies with >= min_turns_per_bin turns in the bin.
min_flies_extend = 7;
size_edges = 7.5:15:367.5;
size_bin_all = discretize(turn_size,size_edges);
n_keep = find(size_edges(2:end) <= 97.5,1,'last');
for k = n_keep+1:numel(size_edges)-1
    ok_all = true;
    for gi = 1:n_g
        c = accumarray(turn_fly(sel & size_bin_all==k),1,[n_flies,1]);
        ok_all = ok_all && sum(c(strcmp(fly_geno,plot_geno{gi})) >= min_turns_per_bin) > min_flies_extend;
    end
    if ~ok_all, break, end
    n_keep = k;
end
size_labels = arrayfun(@(k) sprintf('%g%s',mean(size_edges(k:k+1)),char(176)),1:n_keep,'UniformOutput',false);
spd_bin    = discretize(turn_speed,speed_edges);
n_sb       = numel(speed_edges)-1;
spd_labels = arrayfun(@(b) sprintf('%g-%g%s/s',speed_edges(b),speed_edges(b+1),char(176)),1:n_sb,'UniformOutput',false);

%% (a) binned by turn size
bin_of = size_bin_all; bin_of(bin_of > n_keep) = NaN;
row_labels = cellfun(@(s) [s ' turns'],size_labels,'UniformOutput',false);
plot_overlay(turns,sel,bin_of,row_labels,turn_fly,fly_geno,plot_geno,plot_label,gcol,turn_t,turn_win_s,min_turns_per_bin)
title(gcf().Children(1),{sprintf('closed loop, %s binned by size: fly heading (black) vs. bump (color)',set_lbl),'mean \pm SEM across flies; dashed = 0.8 \times fly'})
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_size_overlay%s.pdf',geno_set,set_tag)),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_size_overlay%s.png',geno_set,set_tag)),'Resolution',200)

%% (b) binned by turn speed (all sizes)
row_labels = cell(1,n_sb);
for b = 1:n_sb
    in_b = sel & spd_bin==b & ismember(fly_geno(turn_fly),plot_geno);
    row_labels{b} = {spd_labels{b},sprintf('(median %.0f%s, %.1f s)',median(turn_size(in_b)),char(176),median([turns(in_b).dur_s]))};
    fprintf('%s: %d turns, median size %.0f deg, median duration %.2f s, %.1f%% longer than %g s\n',spd_labels{b},sum(in_b), ...
        median(turn_size(in_b)),median([turns(in_b).dur_s]),100*mean([turns(in_b).dur_s] > turn_win_s(2)),turn_win_s(2));
end
plot_overlay(turns,sel,spd_bin,row_labels,turn_fly,fly_geno,plot_geno,plot_label,gcol,turn_t,turn_win_s,min_turns_per_bin)
title(gcf().Children(1),{sprintf('closed loop, %s (all sizes) binned by mean speed',set_lbl),'fly heading (black) vs. bump (color), mean \pm SEM across flies','dashed = 0.8 \times fly; row label = median turn size, duration'})
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_speed_overlay%s.pdf',geno_set,set_tag)),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_speed_overlay%s.png',geno_set,set_tag)),'Resolution',200)

%% (c) summary: peak fly / bump displacement and their ratio, by size (top) and speed (bottom)
% Same measure as the size-binned summary in fig_rnai_turns_quiet_cl.m:
% peak of each fly's mean trace (>= min_turns_per_bin turns in the bin)
% within 0..peak_win_s after onset, displacement from the onset value.
% Faint = flies, thick = genotype mean +/- SEM, genotypes side by side.
peak_win_s = 3;
t_pk = turn_t >= 0 & turn_t <= peak_win_s;
binnings = struct('row_of',{bin_of,spd_bin},'labels',{size_labels,spd_labels},'xlab',{'turn size','mean turn speed (all sizes)'});
figure('Color','w','Position',[60,60,1500,820]); clf
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for r = 1:2
    n_r = numel(binnings(r).labels);
    fly_pk = nan(n_flies,n_r); bump_pk = nan(n_flies,n_r);
    for f = unique(turn_fly(sel & ismember(fly_geno(turn_fly),plot_geno)))'
        for b = 1:n_r
            idx = find(sel & turn_fly==f & binnings(r).row_of==b);
            if numel(idx) < min_turns_per_bin, continue, end
            F = cat(2,turns(idx).fly_trace); B = cat(2,turns(idx).bump_trace);
            fly_pk(f,b)  = max(mean(F(t_pk,:),2,'omitnan'));
            bump_pk(f,b) = max(mean(B(t_pk,:),2,'omitnan'));
        end
    end
    panels = {fly_pk, 'peak fly rotation (deg)'; bump_pk, 'peak bump rotation (deg)'; bump_pk./fly_pk, 'bump peak / fly peak'};
    for k = 1:3
        nexttile; hold on
        groupbins(panels{k,1},fly_geno,plot_geno,plot_label,gcol)
        if k == 2, plot(1:n_r,mean(fly_pk(ismember(fly_geno,plot_geno),:),1,'omitnan'),'k:','LineWidth',1.5,'DisplayName','fly peak, these flies'), end
        if k == 3
            yline(1,':k','HandleVisibility','off'); yline(0.8,'--','Color',[.5,.5,.5],'HandleVisibility','off'); yline(0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
            ylim([-0.5,3])
        end
        if k == 2 && r == 1, legend('Location','northwest'), end
        xticks(1:n_r); xticklabels(binnings(r).labels); xlim([0.5,n_r+0.5])
        xlabel(binnings(r).xlab); ylabel(panels{k,2})
    end
end
title(tl,{sprintf('closed loop, %s: peak of per-fly mean trace, max over 0-%.1f s after onset',set_lbl,peak_win_s), ...
    'top: binned by size; bottom: binned by mean speed (all sizes); faint = flies, thick = mean \pm SEM; dashed = 0.8 cue gain'})
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_peak_summary%s.pdf',geno_set,set_tag)),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_peak_summary%s.png',geno_set,set_tag)),'Resolution',200)

%% (d) bump peak / fly peak for a sub-range of turn sizes, with stats
% The ratio panel of (c), restricted to the size bins whose center lies
% strictly inside ratio_range_deg (e.g. [90,Inf] -> 105..135, [0,90] ->
% 15..75). Ranksum for every pair of genotypes, per bin and pooled (each
% fly's mean ratio over the bins it has in the range).
ratio_range_deg = [0,Inf];
ratio_ymax = 3; % y-axis cap; tiny fly peaks on small turns can give very large ratios
bin_ctr = mean([size_edges(1:n_keep); size_edges(2:n_keep+1)]);
rb = find(bin_ctr > ratio_range_deg(1) & bin_ctr < ratio_range_deg(2)); n_rb = numel(rb);
fly_pk = nan(n_flies,n_rb); bump_pk = nan(n_flies,n_rb);
for f = unique(turn_fly(sel & ismember(fly_geno(turn_fly),plot_geno)))'
    for j = 1:n_rb
        idx = find(sel & turn_fly==f & bin_of==rb(j));
        if numel(idx) < min_turns_per_bin, continue, end
        F = cat(2,turns(idx).fly_trace); B = cat(2,turns(idx).bump_trace);
        fly_pk(f,j) = max(mean(F(t_pk,:),2,'omitnan')); bump_pk(f,j) = max(mean(B(t_pk,:),2,'omitnan'));
    end
end
R = bump_pk./fly_pk;
figure('Color','w','Position',[100,100,140+110*n_rb,560]); clf; hold on
groupbins(R,fly_geno,plot_geno,plot_label,gcol)
yline(1,':k','HandleVisibility','off'); yline(0.8,'--','Color',[.5,.5,.5],'HandleVisibility','off'); yline(0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
xticks(1:n_rb); xticklabels(size_labels(rb)); xlim([0.5,n_rb+0.5]); ylim([-0.2,ratio_ymax])
xlabel('turn size'); ylabel('bump peak / fly peak'); legend('Location','southoutside')
range_lbl = sprintf('%s-%s%s',num2str(bin_ctr(rb(1))),num2str(bin_ctr(rb(end))),char(176));
n_off = sum(R(ismember(fly_geno,plot_geno),:) > ratio_ymax,'all'); % fly values above the axis (still in mean/SEM and stats)
title({sprintf('closed loop, %s, %s: bump peak / fly peak (0-%g s)',set_lbl,range_lbl,peak_win_s), ...
    sprintf('dots = flies, mean \\pm SEM; dotted = 1, dashed = 0.8; %d fly value(s) > %g off-axis',n_off,ratio_ymax)},'FontWeight','normal')
range_tag = sprintf('_%g-%gdeg',bin_ctr(rb(1)),bin_ctr(rb(end)));
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_peak_ratio%s%s.png',geno_set,range_tag,set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_peak_ratio%s%s.pdf',geno_set,range_tag,set_tag)),'ContentType','vector')
for j = 0:n_rb
    if j == 0, V = mean(R,2,'omitnan'); lbl = ['pooled ' range_lbl]; else, V = R(:,j); lbl = size_labels{rb(j)}; end
    Vg = cellfun(@(g) V(strcmp(fly_geno,g) & ~isnan(V)),plot_geno,'UniformOutput',false);
    s = [sprintf('%-14s',lbl), strjoin(cellfun(@(g,v) sprintf('%s %.2f (n=%d)',g,mean(v),numel(v)),plot_geno,Vg,'UniformOutput',false),', ')];
    for g1 = 1:n_g-1
        for g2 = g1+1:n_g
            if numel(Vg{g1}) < 2 || numel(Vg{g2}) < 2, continue, end
            s = [s sprintf(' | %s vs %s p=%.3f',plot_geno{g1},plot_geno{g2},ranksum(Vg{g1},Vg{g2}))]; %#ok<AGROW>
        end
    end
    disp(s)
end

%% functions
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

function plot_overlay(turns,sel,row_of,row_labels,turn_fly,fly_geno,plot_geno,plot_label,gcol,turn_t,turn_win_s,min_turns)
    % rows = bins (row_of: bin index per turn, NaN = unused), columns =
    % genotypes; per-fly mean traces, then mean +/- SEM across flies
    n_g = numel(plot_geno); n_r = numel(row_labels);
    figure('Color','w','Position',[60,60,240*n_g,170*n_r]); clf
    tl = tiledlayout(n_r,n_g,'TileSpacing','compact','Padding','compact'); %#ok<NASGU>
    ax = gobjects(n_g,n_r);
    for gi = 1:n_g
        flies_g = find(strcmp(fly_geno,plot_geno{gi}))';
        for r = 1:n_r
            FT = nan(numel(flies_g),numel(turn_t)); BT = FT;
            for k = 1:numel(flies_g)
                idx = find(sel & turn_fly==flies_g(k) & row_of==r);
                if numel(idx) < min_turns, continue, end
                FT(k,:) = mean(cat(2,turns(idx).fly_trace),2,'omitnan');
                BT(k,:) = mean(cat(2,turns(idx).bump_trace),2,'omitnan');
            end
            ax(gi,r) = nexttile((r-1)*n_g+gi); hold on
            n_f = sum(~all(isnan(FT),2));
            if n_f >= 1
                plot_sem(turn_t,FT,[0,0,0])
                plot_sem(turn_t,BT,gcol(gi,:))
                plot(turn_t,0.8*mean(FT,1,'omitnan'),'--','Color',[.5,.5,.5])
            end
            xline(0,':k'); xlim(turn_win_s)
            if r == 1, title(plot_label{gi}), end
            text(0.03,0.97,sprintf('n=%d flies',n_f),'Units','normalized','VerticalAlignment','top','FontSize',8)
            if gi == 1
                lbl = row_labels{r}; if ~iscell(lbl), lbl = {lbl}; end
                ylabel([lbl(:)', {'rotation (deg)'}])
            end
            if r == n_r, xlabel('time from turn onset (s)'), end
        end
    end
    for r = 1:n_r, linkaxes(ax(:,r),'y'), end
    for a = ax(:)', a.Toolbar.Visible = 'off'; end
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
