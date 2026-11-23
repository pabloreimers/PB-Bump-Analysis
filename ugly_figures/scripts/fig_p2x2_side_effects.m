%% fig_p2x2_side_effects
% Scratch: does WHICH PB half got the ATP matter? Two questions, one row each,
% one line per fly (left-side vs right-side trials of the same fly):
%   row 1 -- bump displacement 0-15 s after each stim, RAW sign (not flipped
%            by side), as fig_p2x2_displacement_per_stim.m (rho < 0.1 masked,
%            bar still, confirmed stims). If LPsP pushes the bump toward /
%            away from the stimulated side, left and right trials have
%            opposite signs.
%   row 2 -- change in bump tracking of the open-loop bright-bar sweeps
%            (post-stim sweeps minus the sweep before stim 1 of the SAME
%            trial), as fig_p2x2_sweep_tracking.m, computed separately for
%            the fly's left-side and right-side bright-bar trials. Bright-bar
%            trials only (in background trials the bump doesn't track the
%            sweep at all). Each fly usually has ~1 bright-bar trial per side,
%            so "pre" is a single ~40 s sweep -- noisy per fly.
% Side = each trial's call from its most confident ejection (as everywhere).
% Needs the workspace of fig_p2x2_displacement_per_stim.m; runs it if missing.

%% setup
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist('stim_conf','var') || ~exist('trial_is_right','var') || ~exist('stim_disp','var')
    run(fullfile(repo_root,'ugly_figures','scripts','fig_p2x2_displacement_per_stim.m')); close all
end
rho_cut = 0.1;
export_dir = fullfile(repo_root,'ugly_figures','exports','V2','scratch');
pat = strings(numel(all_data),1);
for i = 1:numel(all_data); p = all_data(i).ft.pattern; if iscell(p); p = p{1}; end; pat(i) = string(p); end
gn = {'empty','lpsp'}; geno_label = {'Empty>P2X2','LPsP>P2X2'}; side_name = {'left','right'};
gainf = @(X) sum(X(:,1).*X(:,2)) / sum(X(:,2).^2);

%% per fly x side: raw displacement (row 1) and tracking change (row 2)
D = struct();
for g = 1:2
    rws = find(strcmp(genotype,gn{g}) & ~isnan(trial_is_right)); fls = unique(fly_id(rws));
    disp_side = nan(numel(fls),2); trk_pre = nan(numel(fls),2); trk_post = nan(numel(fls),2); % columns: left, right
    for f = 1:numel(fls)
        for side = 0:1
            tri = rws(strcmp(fly_id(rws),fls{f}) & trial_is_right(rws)==side)';
            d = [];
            for i = tri; x = stim_disp{i}(stim_conf{i},1); d = [d; x(~isnan(x))]; end %#ok<AGROW>
            if ~isempty(d); disp_side(f,side+1) = mean(d); end
            Pre = []; Post = [];
            for i = tri(pat(tri)=="0003_4px_brightbar.mat")
                tr = all_data(i); xb = tr.ft.xb(:); xf = tr.ft.xf(:); dt = median(diff(xf));
                m = masked_mu(tr.im.mu(:), tr.im.rho(:), xb, rho_cut);
                c = unwrap(tr.ft.cue(:)); v = [0; diff(c)] ./ [1; diff(xf)];
                mv = bwareaopen(abs(smoothdata(v,'movmean',round(.5/dt))) > 0.05, round(1/dt));
                on = xf(find(diff([false; logical(tr.ft.stims(:))])==1)); t0 = xf(find(mv,1));
                if isempty(on) || isempty(t0); continue; end
                mb = @(t) mean(m(xb>=t-.5 & xb<t+.5)); mc = @(t) mean(c(xf>=t-.5 & xf<t+.5));
                if on(1) - 0.5 - t0 > 5; Pre(end+1,:) = [mb(on(1)-.5)-mb(t0), -(mc(on(1)-.5)-mc(t0))]; end %#ok<AGROW>
                for s = find(stim_conf{i})'
                    a = xf(find(xf > on(s) & mv, 1)); if isempty(a); continue; end
                    if s < numel(on); b = on(s+1) - 0.5; else; b = xf(find(mv,1,'last')); end
                    if b - a < 5; continue; end
                    Post(end+1,:) = [mb(b)-mb(a), -(mc(b)-mc(a))]; %#ok<AGROW>
                end
            end
            if ~isempty(Pre) && ~isempty(Post); trk_pre(f,side+1) = gainf(Pre); trk_post(f,side+1) = gainf(Post); end
        end
    end
    D(g).disp = disp_side; D(g).pre = trk_pre; D(g).post = trk_post; D(g).dtrk = trk_post - trk_pre;
end

for g = 1:2
    x = D(g).disp; b = all(~isnan(x),2);
    fprintf('%-11s displacement (raw sign): left %+.2f +/- %.2f, right %+.2f +/- %.2f rad; paired L vs R p=%.3g (n=%d flies)\n', ...
        geno_label{g}, mean(x(:,1),'omitnan'), sem(x(:,1)), mean(x(:,2),'omitnan'), sem(x(:,2)), signrank(x(b,1),x(b,2)), sum(b));
    y = D(g).dtrk; b = all(~isnan(y),2);
    for side = 1:2
        v = ~isnan(D(g).pre(:,side));
        fprintf('            tracking, %s-side bar trials: pre %.2f -> post %.2f (change %+.2f +/- %.2f, signrank p=%.3g, n=%d)\n', ...
            side_name{side}, mean(D(g).pre(v,side)), mean(D(g).post(v,side)), mean(y(v,side)), sem(y(v,side)), signrank(y(v,side)), sum(v));
    end
    fprintf('            change, left vs right paired p=%.3g (n=%d flies with both)\n', signrank(y(b,1),y(b,2)), sum(b));
end

%% plot
gcol = [.5 .5 .5; 0 .2 1];
fig = figure('Color','w','Position',[100 100 820 820]);
tl = tiledlayout(fig,2,2,'TileSpacing','loose');
rows_ = {'disp','dtrk'};
ylabs = {'bump displacement 0-15 s per stim (rad, raw sign)', '\Delta bar tracking gain (post - pre stim)'};
for r = 1:2
    vals = arrayfun(@(g) D(g).(rows_{r}), 1:2, 'UniformOutput', false);
    allv = vertcat(vals{:}); yl = [min(allv(:)) max(allv(:))]; yl = yl + [-.08 .2]*diff(yl);
    for g = 1:2
        ax = nexttile(tl); hold(ax,'on')
        x = vals{g}; b = all(~isnan(x),2);
        plot(ax, [1 2], x(b,:)', '-', 'Color', [gcol(g,:) .4])
        for side = 1:2
            v = ~isnan(x(:,side));
            scatter(ax, side*ones(sum(v),1), x(v,side), 30, gcol(g,:), 'filled', 'MarkerFaceAlpha', .6)
        end
        errorbar(ax, [1 2] + [-.15 .15], mean(x,1,'omitnan'), [sem(x(:,1)) sem(x(:,2))], 'ok-', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 0)
        yline(ax, 0, ':k');
        xlim(ax,[.5 2.5]); xticks(ax,[1 2]); xticklabels(ax,{'left-side stims','right-side stims'}); ylim(ax, yl)
        if g == 1; ylabel(ax, ylabs{r}); end
        title(ax, sprintf('%s: left vs right p=%.2g (n=%d)', geno_label{g}, signrank(x(b,1),x(b,2)), sum(b)))
    end
end
title(tl, {'does the stimulated side matter? one line per fly', ...
    'top: bump push after each stim (0-15 s, bar still) | bottom: change in tracking of the bright-bar sweeps'}, 'FontSize', 11)
if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig, fullfile(export_dir,'fig_p2x2_side_effects.png'), 'Resolution', 200)

%% functions
function m = masked_mu(mu, rho, t, rho_min)
    ok = rho >= rho_min & ~isnan(mu);
    m = nan(size(mu)); m(ok) = unwrap(mu(ok));
    m = interp1(t(ok), m(ok), t, 'linear'); m = fillmissing(m, 'nearest');
end
function s = sem(x); s = std(x,'omitnan') / sqrt(sum(~isnan(x))); end
