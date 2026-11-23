%% fig_rnai_offset_variability
% Bump-cue offset variability per genotype, one figure per set (TH set:
% LPsP>TH, Empty>TH, LPsP>mCherry, plus LPsP>TH split by closed-loop
% integrative gain; VGlut set: LPsP>VGlut, Empty>VGlut, LPsP>mCherry).
%
%   offset variability = circ_var(circ_dist(mu, cue)) per closed-loop trial,
%   over the timepoints where the bump is solid (PVA strength rho >= rho_min)
%   AND the fly is turning (|turning speed| >= turn_min) -- bump drift while
%   fly and bar sit still doesn't count. Per fly = mean over its trials (the
%   offset can legitimately reset between trials, so timepoints are not
%   pooled across trials). 0 = the bump holds a fixed offset to the bar,
%   1 = no consistent relationship.
%   mu / rho: PVA of im.f (movmean 3 frames, z-scored per glomerulus), as in
%   the lap / stall scripts; cue = -ft.cue (same sign convention as mu),
%   unwrapped; turning speed = ft.r_speed, turn_smooth_s Gaussian; both
%   interpolated onto the imaging timebase.
%   Trials: closed loop, passing the integrative-gain inclusion (walking +
%   working cue + bump); a trial is dropped if the fly is still (neither
%   walking nor turning: |turning| < still_turn and |forward| < still_fwd) for
%   more than max_still_frac of it, or has fewer than
%   min_samples timepoints passing both thresholds. Flies with no trials left
%   are excluded.
% p: ranksum on per-fly values, main genotype vs each control (and LPsP>TH
% high vs low gain). Needs all_data (loaded if absent) and the integrative-
% gain cache. Run from the repo root.

%% parameters + labels
rho_min = 0.4;          % PVA strength threshold (0.4 = the lap / stall threshold)
turn_min = 0.5;         % rad/s (~30 deg/s), |turning speed| threshold (integrative-gain sample threshold)
turn_smooth_s = 0.5;    % s, Gaussian smoothing of ft.r_speed before thresholding
max_still_frac = 0.9;   % drop trials where the fly is still for more than this fraction
still_turn = 0.2; still_fwd = 0.1; % still = |turning| < still_turn rad/s and |forward| < still_fwd (smoothed as turning)
min_samples = 100;      % timepoints per trial passing both thresholds
gain_split = 0.85;
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
addpath(fullfile(pwd,'circ_stats'));
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
if ~exist('all_data','var')
    tmp = load(fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
g_fly = nan(max(IG.fly_num),1); % closed-loop integrative gain per fly
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
out_dir = fullfile('ugly_figures','rnai_turns','offset_variability');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% offset variability per trial, then per fly
trials = find(~IG.is_dark & IG.inc & ismember(IG.genotype,{'lpsp>th','empty>th','lpsp>mcherry','lpsp>vglut','empty>vglut'}));
ov_trial = nan(numel(IG.genotype),1); still_frac = nan(numel(IG.genotype),1);
for i = trials(:)'
    tr = all_data(i); xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    cu = -tr.ft.cue(:); ok = ~isnan(cu); cu(ok) = unwrap(cu(ok)); cu = fillmissing(cu,'linear','EndValues','nearest');
    cb = interp1(xf,cu,xb);
    rs = fillmissing(tr.ft.r_speed(:),'constant',0);
    rs = smoothdata(rs,'gaussian',max(1,round(turn_smooth_s/median(diff(xf)))));
    turning = abs(interp1(xf,rs,xb)) >= turn_min;
    fs = smoothdata(fillmissing(tr.ft.f_speed(:),'constant',0),'gaussian',max(1,round(turn_smooth_s/median(diff(xf)))));
    still_frac(i) = mean(abs(rs) < still_turn & abs(fs) < still_fwd); % neither walking nor turning
    idx = rho >= rho_min & turning & ~isnan(cb);
    if still_frac(i) > max_still_frac || sum(idx) < min_samples, continue, end
    ov_trial(i) = circ_var(circ_dist(mu(idx),cb(idx)));
end
n_flies = max(IG.fly_num);
ov_fly = nan(n_flies,1);
for f = unique(IG.fly_num(trials))', ov_fly(f) = mean(ov_trial(IG.fly_num==f),'omitnan'); end
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = IG.genotype{find(IG.fly_num==f,1)}; end

%% figures: one per set
sets = struct('tag',{'th','vglut'}, ...
    'grp',{{'LPsP > TH-RNAi','lpsp>th',[]; 'Empty > TH-RNAi','empty>th',[]; 'LPsP > mCherry','lpsp>mcherry',[]; ...
            'LPsP>TH high gain','lpsp>th',true; 'LPsP>TH low gain','lpsp>th',false}, ...
           {'LPsP > VGlut-RNAi','lpsp>vglut',[]; 'Empty > VGlut-RNAi','empty>vglut',[]; 'LPsP > mCherry','lpsp>mcherry',[]}}, ...
    'col',{[0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13; 0.60,0.15,0.00; 0.95,0.62,0.45], ...
           [0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13]}, ...
    'x',{[1,2,3,4.6,5.6],1:3},'cmp',{[1,2; 1,3; 4,5],[1,2; 1,3]});
for si = 1:numel(sets)
    S = sets(si); n_g = size(S.grp,1); V = cell(n_g,1); n_all = zeros(n_g,1);
    for gi = 1:n_g
        fl = find(strcmp(fly_geno,S.grp{gi,2}) & ismember((1:n_flies)',IG.fly_num(trials)));
        if ~isempty(S.grp{gi,3}), fl = fl(~isnan(g_fly(fl)) & (g_fly(fl) > gain_split) == S.grp{gi,3}); end
        n_all(gi) = numel(fl); V{gi} = ov_fly(fl(~isnan(ov_fly(fl))));
    end
    p = arrayfun(@(c) ranksum(V{S.cmp(c,1)},V{S.cmp(c,2)}),1:size(S.cmp,1))';
    fprintf('\n%s set (rho >= %g, |turning| >= %g rad/s; trials still > %g%% dropped):\n',upper(S.tag),rho_min,turn_min,100*max_still_frac);
    for gi = 1:n_g, fprintf('  %-20s offset variability %.3f +/- %.3f (n=%d of %d flies)\n',S.grp{gi,1},mean(V{gi}),std(V{gi})/sqrt(numel(V{gi})),numel(V{gi}),n_all(gi)); end
    for c = 1:size(S.cmp,1), fprintf('  %s vs %s: ranksum p = %.2g\n',S.grp{S.cmp(c,1),1},S.grp{S.cmp(c,2),1},p(c)); end

    figure('Color','w','Position',[80,80,140+110*n_g,520]); clf; hold on
    for gi = 1:n_g
        y = V{gi};
        scatter(S.x(gi) + (rand(size(y))-.5)*.4,y,22,S.col(gi,:),'filled','MarkerFaceAlpha',.55)
        errorbar(S.x(gi),mean(y),std(y)/sqrt(numel(y)),'o','Color','k','MarkerFaceColor',S.col(gi,:),'MarkerSize',8,'LineWidth',1.2)
    end
    ylim([0,1]); add_p(gca,S.x,S.cmp,p)
    xticks(S.x); xticklabels(arrayfun(@(gi) sprintf('%s (n=%d)',S.grp{gi,1},numel(V{gi})),1:n_g,'UniformOutput',false)); xtickangle(35)
    xlim([0.4,max(S.x)+0.6]); ylabel('offset variability, circ\_var(circ\_dist(\mu, cue))')
    title({'bump-cue offset variability, closed loop', ...
        sprintf('timepoints with \\rho >= %g and |turning| >= %g rad/s; dots = flies, mean \\pm SEM; ranksum',rho_min,turn_min)},'FontWeight','normal','FontSize',10)
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_offset_variability.png',S.tag)),'Resolution',200)
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_offset_variability.pdf',S.tag)),'ContentType','vector')
end

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
