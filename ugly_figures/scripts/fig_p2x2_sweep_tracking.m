%% fig_p2x2_sweep_tracking
% Scratch: how well does the bump follow the SWEEPING bar after an ATP stim,
% LPsP>P2X2 vs Empty>P2X2? Companion to fig_p2x2_displacement_per_stim.m,
% which scores only the period when the bar is still.
%
% Schedule (lpsp_p2x2_reredo_20260518, both trial types): the visual pattern
% is swept in open loop from the trial start until 0.5 s before stim 1, then
% holds still until 18.6 s after each stim, then sweeps (~47 rad of cue
% motion, ~1.3 rad/s) until 0.5 s before the next stim; after stim 4 it
% sweeps until ~274 s. Sweeps are found from ft.cue as in the displacement
% script (|cue velocity| > 0.05 rad/s, 0.5 s smoothing, bouts < 1 s dropped).
%
% Per sweep: bump displacement (im.mu, rho < rho_min masked + interpolated,
% as the displacement script; 1 s means at each end) and cue displacement,
% sign-flipped so that bump = -cue is perfect tracking (the codebase's
% heading convention). Per fly, TRACKING GAIN = least-squares slope of bump
% on -cue displacement through 0, pooled over that fly's sweeps of a kind:
%   pre  -- the sweep at the start of each trial, before stim 1 (only the
%           fly's FIRST trial is truly ATP-naive; later trials follow
%           earlier trials' stims)
%   post -- the sweeps after each confirmed-ejection stim
% Bright-bar trials (0003_4px_brightbar) only in the plot: in background
% trials (0001_background) the bump doesn't follow the sweep at all in either
% genotype (gain ~0.00-0.01; printed).

%% setup
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist('stim_conf','var') || ~exist('trial_is_right','var')
    run(fullfile(repo_root,'ugly_figures','scripts','fig_p2x2_displacement_per_stim.m')); close all
end
rho_cut = 0.1;
export_dir = fullfile(repo_root,'ugly_figures','exports','V2','scratch');
pat = strings(numel(all_data),1);
for i = 1:numel(all_data); p = all_data(i).ft.pattern; if iscell(p); p = p{1}; end; pat(i) = string(p); end
pats = ["0003_4px_brightbar.mat", "0001_background.mat"]; pat_label = {'bright bar','background'};
gn = {'empty','lpsp'}; geno_label = {'Empty>P2X2','LPsP>P2X2'};

%% per fly tracking gain, pre vs post, for each pattern
gainf = @(X) sum(X(:,1).*X(:,2)) / sum(X(:,2).^2);
G = struct('pre',{},'post',{});
for pj = 1:2
    for g = 1:2
        rws = find(strcmp(genotype,gn{g}) & ~isnan(trial_is_right) & pat==pats(pj)); fls = unique(fly_id(rws));
        pre_g = nan(numel(fls),1); post_g = pre_g;
        for f = 1:numel(fls)
            Pre = []; Post = [];
            for i = rws(strcmp(fly_id(rws),fls{f}))'
                tr = all_data(i); xb = tr.ft.xb(:); xf = tr.ft.xf(:); dt = median(diff(xf));
                m = masked_mu(tr.im.mu(:), tr.im.rho(:), xb, rho_cut);
                c = unwrap(tr.ft.cue(:)); v = [0; diff(c)] ./ [1; diff(xf)];
                mv = bwareaopen(abs(smoothdata(v,'movmean',round(.5/dt))) > 0.05, round(1/dt));
                on = xf(find(diff([false; logical(tr.ft.stims(:))])==1));
                t0 = xf(find(mv,1));
                if isempty(on) || isempty(t0); continue; end % no stims or no sweep recorded in this trial
                mb = @(t) mean(m(xb>=t-.5 & xb<t+.5)); mc = @(t) mean(c(xf>=t-.5 & xf<t+.5));
                if on(1) - 0.5 - t0 > 5; Pre(end+1,:) = [mb(on(1)-.5)-mb(t0), -(mc(on(1)-.5)-mc(t0))]; end %#ok<AGROW>
                for s = find(stim_conf{i})'
                    a = xf(find(xf > on(s) & mv, 1)); if isempty(a); continue; end
                    if s < numel(on); b = on(s+1) - 0.5; else; b = xf(find(mv,1,'last')); end
                    if b - a < 5; continue; end
                    Post(end+1,:) = [mb(b)-mb(a), -(mc(b)-mc(a))]; %#ok<AGROW>
                end
            end
            if ~isempty(Pre);  pre_g(f)  = gainf(Pre);  end
            if ~isempty(Post); post_g(f) = gainf(Post); end
        end
        G(g,pj).pre = pre_g; G(g,pj).post = post_g;
    end
    fprintf('\n== %s trials: tracking gain (bump / -cue displacement during sweeps) ==\n', pat_label{pj});
    for g = 1:2
        v = ~isnan(G(g,pj).pre) & ~isnan(G(g,pj).post);
        fprintf('  %-11s n=%2d flies: pre-stim sweep %.2f +/- %.2f | post-stim sweeps %.2f +/- %.2f | signrank pre vs post p=%.3g\n', ...
            geno_label{g}, sum(v), mean(G(g,pj).pre(v)), sem(G(g,pj).pre(v)), mean(G(g,pj).post(v)), sem(G(g,pj).post(v)), signrank(G(g,pj).pre(v), G(g,pj).post(v)));
    end
    d = arrayfun(@(g) G(g,pj).post - G(g,pj).pre, 1:2, 'UniformOutput', false);
    fprintf('  ranksum LPsP vs Empty: pre p=%.3g | post p=%.3g | change (post-pre) p=%.3g\n', ...
        ranksum(rmnan(G(1,pj).pre), rmnan(G(2,pj).pre)), ranksum(rmnan(G(1,pj).post), rmnan(G(2,pj).post)), ranksum(rmnan(d{1}), rmnan(d{2})));
end

%% plot: bright-bar trials, pre vs post per fly, one panel per genotype
gcol = [.5 .5 .5; 0 .2 1];
fig = figure('Color','w','Position',[100 100 700 440]);
tl = tiledlayout(fig,1,2,'TileSpacing','loose');
allv = vertcat(G(1,1).pre, G(1,1).post, G(2,1).pre, G(2,1).post);
yl_all = [min(allv) max(allv)];
for g = 1:2
    ax = nexttile(tl); hold(ax,'on')
    v = ~isnan(G(g,1).pre) & ~isnan(G(g,1).post); a = G(g,1).pre(v); b = G(g,1).post(v);
    plot(ax, [1 2], [a b]', '-', 'Color', [gcol(g,:) .4])
    scatter(ax, ones(size(a)), a, 30, gcol(g,:), 'filled', 'MarkerFaceAlpha', .6)
    scatter(ax, 2*ones(size(b)), b, 30, gcol(g,:), 'filled', 'MarkerFaceAlpha', .6)
    errorbar(ax, [1 2]+[-.15 .15], [mean(a) mean(b)], [sem(a) sem(b)], 'ok-', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 0)
    yline(ax, 0, ':k');
    xlim(ax,[.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,{'before stim 1','after stims'})
    ylim(ax, yl_all + [-.05 .15]*diff(yl_all))
    ylabel(ax,'bar tracking gain (bump / bar displacement)')
    title(ax, sprintf('%s (n=%d flies)\npre vs post p=%.2g', geno_label{g}, sum(v), signrank(a,b)))
end
title(tl, 'bump tracking of the open-loop bright-bar sweeps, one line per fly', 'FontSize', 11)
if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig, fullfile(export_dir,'fig_p2x2_sweep_tracking.png'), 'Resolution', 200)

%% functions
function m = masked_mu(mu, rho, t, rho_min)
    ok = rho >= rho_min & ~isnan(mu);
    m = nan(size(mu)); m(ok) = unwrap(mu(ok));
    m = interp1(t(ok), m(ok), t, 'linear'); m = fillmissing(m, 'nearest');
end
function s = sem(x); s = std(x,'omitnan') / sqrt(sum(~isnan(x))); end
function x = rmnan(x); x = x(~isnan(x)); end
