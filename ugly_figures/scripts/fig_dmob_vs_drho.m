%% fig_dmob_vs_drho
% Scratch: are the per-fly changes in bump mobility and bump rho (dark
% trials, post - pre ATP) correlated? One dot per fly, the same flies,
% bouts and pre/post pools as Figure 3's two delta panels (walk_delta,
% walk_rho_delta from Figure_3_claude.m). Spearman per genotype and pooled
% (pooled mixes a genotype difference in both axes, so read it with care),
% least-squares line per genotype for the eye.
%
% Context: Figure_3_claude.m's RHO AND MOBILITY ARE COUPLED caveat found
% trial-mean rho vs mobility Spearman -0.40 across dark trials (n=73). This
% asks the per-fly CHANGE version of that question.
%
% REQUIRES the Figure_3_claude.m workspace.

%% check workspace
need = {'walk_delta','walk_rho_delta','walk_fly_group','walk_fly_list','gcamp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');

%% stats
ok = ~isnan(walk_delta) & ~isnan(walk_rho_delta);
x = walk_rho_delta(ok); y = walk_delta(ok); g = walk_fly_group(ok);
geno_label = {'empty>P2X2','lpsp>P2X2'}; gcol = [.5 .5 .5; gcamp_color];
r = nan(1,3); p = nan(1,3); n = nan(1,3);
for gi = 1:2
    [r(gi),p(gi)] = corr(x(g==gi), y(g==gi), 'type','Spearman'); n(gi) = sum(g==gi);
end
[r(3),p(3)] = corr(x, y, 'type','Spearman'); n(3) = numel(x);
for k = 1:3
    lab = [geno_label, {'pooled'}];
    fprintf('%-11s n=%2d  Spearman rho=%+.2f  p=%.3f\n', lab{k}, n(k), r(k), p(k));
end
fly_ok = walk_fly_list(ok);

%% plot
fig = figure('color','w','Visible','off','Position',[100 100 620 640]);
ax = axes(fig); hold(ax,'on')
xline(ax,0,':k','HandleVisibility','off'); yline(ax,0,':k','HandleVisibility','off')
h = gobjects(1,2);
for gi = 1:2
    xi = x(g==gi); yi = y(g==gi);
    b = [ones(numel(xi),1) xi] \ yi; xx = [min(xi) max(xi)];
    plot(ax, xx, b(1)+b(2)*xx, '-', 'Color', gcol(gi,:), 'LineWidth', 1.5, 'HandleVisibility','off')
    h(gi) = scatter(ax, xi, yi, 60, gcol(gi,:), 'filled', 'MarkerFaceAlpha', .7, ...
        'DisplayName', sprintf('%s: Spearman %+.2f, p=%.2g (n=%d)', geno_label{gi}, r(gi), p(gi), n(gi)));
end
xlabel(ax,'\Delta mean bout \rho (post - pre)')
ylabel(ax,'\Delta bump mobility (post - pre)')
legend(ax,h,'Location','southoutside','FontSize',9) % outside: inside, it hid the +2.77 lpsp fly (20260727 fly 4) in the top-right corner
title(ax,{'per-fly change in bump mobility vs. change in bout \rho, dark trials', sprintf('pooled Spearman %+.2f, p=%.2g (n=%d)', r(3), p(3), n(3))},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_dmob_vs_drho.png'), 'Resolution', 200)
fig.Visible = 'on';
