%% fig_gain07_heading_clean
% Before/after of pb_clean_heading on the fictrac heading of every trial in
% .data\epg_7f_20261007.mat that has glitches (see fig_gain07_heading_glitches).
% Left: whole trial, raw (black) vs cleaned (red). Right: zoom on the largest
% glitch. Nothing is written back to the dataset here; the cleaning is wired
% into pb_gain_trial for the next batch run.

repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(repoRootDir, 'claude'));
S = load(fullfile(repoRootDir, '.data', 'epg_7f_20261007.mat'), 'all_data'); all_data = S.all_data; clear S
exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports', 'epg_7f_gain07');

n = numel(all_data); rows = {};
for i = 1:n
    d = all_data(i); h = d.ft.heading_raw(:)';
    [hc, info] = pb_clean_heading(h);
    rows(end+1, :) = {d.meta.fly_num, d.meta.trial_num, info.nBad, size(info.segs, 1), round(rad2deg(info.thresh), 1), round(rad2deg(info.netRemoved)), round(rad2deg(max(abs(h - hc))))}; %#ok<SAGROW>
end
T = cell2table(rows, 'VariableNames', {'fly', 'trial', 'frames_replaced', 'n_runs', 'thresh_deg_per_frame', 'net_heading_removed_deg', 'max_abs_change_deg'});
disp(T)
writetable(T, fullfile(exportDir, '_heading_clean_table.csv'));

show = find(T.frames_replaced > 0);
fig = figure('Visible', 'off', 'Position', [50 50 2000 1500], 'Color', 'w');
for k = 1:numel(show)
    d = all_data(show(k)); h = d.ft.heading_raw(:)'; t = d.ft.xf(:)';
    [hc, info] = pb_clean_heading(h);
    subplot(numel(show), 2, 2*k-1); hold on
    plot(t, rad2deg(h), 'k', 'LineWidth', 1); plot(t, rad2deg(hc), 'r', 'LineWidth', 1)
    for s = 1:size(info.segs, 1), xline(t(info.segs(s, 1)), ':', 'Color', [.5 .5 .5]); end
    ylabel('heading (deg)'); legend({'raw fictrac', 'cleaned'}, 'Location', 'best')
    title(sprintf('fly %d trial %d: %d frames replaced in %d runs, net %.0f deg removed', d.meta.fly_num, d.meta.trial_num, info.nBad, size(info.segs, 1), rad2deg(info.netRemoved)))
    subplot(numel(show), 2, 2*k); hold on
    dh = diff(h); [~, j] = max(abs(dh)); w = max(j-150, 1):min(j+150, numel(h));
    plot(t(w), rad2deg(h(w)), 'k.-', 'MarkerSize', 5); plot(t(w), rad2deg(hc(w)), 'r.-', 'MarkerSize', 5)
    title(sprintf('zoom +/-2.5 s around the largest jump (threshold %.1f deg/frame)', rad2deg(info.thresh)))
    if k == numel(show), xlabel('fictrac time (s)'); end
end
pb_export_png(fig, fullfile(exportDir, '_heading_clean_before_after.png'), 110);
fprintf('done\n');
