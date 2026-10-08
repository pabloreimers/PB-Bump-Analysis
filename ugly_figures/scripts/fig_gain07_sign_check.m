%% fig_gain07_sign_check
% Is the "glomerulus-order sign" just the direction the skeleton was
% numbered in? For every trial of .data\epg_7f_20261007.mat: the sign the
% pipeline chose (bar_sign), the image-x of bin 1 vs bin 32 (numbering
% direction), and the raw (unflipped) bump-vs-bar relation. Then two example
% bar trials, one of each sign, plotted WITHOUT any sign flipping.

repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(repoRootDir, 'claude'));
S = load(fullfile(repoRootDir, '.data', 'epg_7f_20261007.mat'), 'all_data'); all_data = S.all_data; clear S
exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports', 'epg_7f_gain07');

n = numel(all_data);
T = table('Size', [n 8], 'VariableTypes', {'double','double','string','double','double','string','double','double'}, ...
    'VariableNames', {'fly','trial','cond','bar_sign','slope_mu_vs_heading','numbering','x_bin1','x_bin32'});
for i = 1:n
    d = all_data(i); ci = d.im.clusterIdx;
    [~, x1] = find(ci == 1); [~, x32] = find(ci == max(ci(:)));
    T.fly(i) = d.meta.fly_num; T.trial(i) = d.meta.trial_num;
    if d.ft.dark, T.cond(i) = "dark"; else, T.cond(i) = "bar"; end
    T.bar_sign(i) = d.bump.bar_sign; T.slope_mu_vs_heading(i) = round(d.bump.slope_mu_vs_heading, 2);
    T.x_bin1(i) = round(mean(x1)); T.x_bin32(i) = round(mean(x32));
    if mean(x1) < mean(x32), T.numbering(i) = "left->right"; else, T.numbering(i) = "right->left"; end
end
disp(T)
fprintf('\nsign vs numbering direction:\n'); disp(crosstab(T.numbering, T.bar_sign))
% expected if the sign is purely numbering: all left->right trials share one sign, all right->left the other

%% example trials, raw (no flipping): two good bar trials from different flies,
% and the one trial whose measured bump/heading slope had the other sign
ex = [find(T.fly == 8 & T.cond == "bar", 1), find(T.fly == 16 & T.cond == "bar", 1), find(T.fly == 13 & T.trial == 2, 1)];
ex = ex(~cellfun(@isempty, num2cell(ex)));
fig = figure('Visible', 'off', 'Position', [50 50 2200 1000], 'Color', 'w');
nEx = numel(ex);
for k = 1:nEx
    d = all_data(ex(k)); t = d.ft.xb; ok = d.bump.ok;
    [~, tf] = fileparts(d.meta.trialDir);
    % glomerulus map with numbering
    ax = subplot(4, nEx, k); cmap = [0 0 0; hsv(size(d.im.z, 1))];
    image(ind2rgb(d.im.clusterIdx + 1, cmap)); axis image; hold on
    for c = [1 8 16 17 24 32]
        [cy, cx] = find(d.im.clusterIdx == c); text(mean(cx), mean(cy), num2str(c), 'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
    title(sprintf('fly %d %s: bins numbered %s (pipeline sign %+d)', d.meta.fly_num, tf, T.numbering(ex(k)), d.bump.bar_sign), 'Interpreter', 'none')
    % raw bump vs raw bar, no flipping
    subplot(4, nEx, k + nEx); hold on
    b = rad2deg(d.ft.cue); b(abs(diff([b b(end)])) > 180) = NaN; plot(t, b, '-', 'Color', [.3 .3 .3], 'LineWidth', 1.2)
    m = rad2deg(d.im.mu); m(~ok) = NaN; m(abs(diff([m m(end)])) > 180) = NaN; plot(t, m, 'b.', 'MarkerSize', 4)
    ylim([-180 180]); yticks(-180:90:180); ylabel('deg'); title('raw PVA angle (blue) and raw bar angle (gray), no sign applied')
    % raw bump + bar and raw bump - bar
    subplot(4, nEx, k + 2*nEx); hold on
    o1 = rad2deg(angle(exp(1i * (d.im.mu - d.ft.cue)))); o2 = rad2deg(angle(exp(1i * (d.im.mu + d.ft.cue)))); o1(~ok) = NaN; o2(~ok) = NaN;
    plot(t, o1, '.', 'Color', [0 .5 0], 'MarkerSize', 4); plot(t, o2, '.', 'Color', [.8 0 .5], 'MarkerSize', 4)
    ylim([-180 180]); yticks(-180:90:180); ylabel('deg')
    legend({sprintf('bump - bar (circ std %.0f)', rad2deg(gce_cs(o1))), sprintf('bump + bar (circ std %.0f)', rad2deg(gce_cs(o2)))}, 'Location', 'eastoutside')
    title('which combination is constant tells you the relation')
    % unwrapped raw bump vs heading
    subplot(4, nEx, k + 3*nEx); hold on
    mu_t = d.im.mu; mu_t(~ok) = NaN; iok = ~isnan(mu_t); u = nan(size(mu_t)); u(iok) = unwrap(mu_t(iok));
    h = d.ft.heading - d.ft.heading(1); h = h - median(h(iok) - u(iok) / d.bump.slope_mu_vs_heading, 'omitnan');
    plot(t, rad2deg(u), 'b.', 'MarkerSize', 4); plot(t, rad2deg(h * d.bump.slope_mu_vs_heading), '-', 'Color', [.85 .33 .1])
    ylabel('unwrapped deg'); xlabel('trial time (s)')
    title(sprintf('raw bump (blue) vs %.2f x heading (orange): slope sign is the measured relation', d.bump.slope_mu_vs_heading))
end
pb_export_png(fig, fullfile(exportDir, '_sign_check_two_flies.png'), 120);
writetable(T, fullfile(exportDir, '_sign_check_table.csv'));
fprintf('done\n');

function s = gce_cs(x)
x = deg2rad(x(~isnan(x))); s = sqrt(-2 * log(max(abs(mean(exp(1i * x(:)))), eps)));
end
