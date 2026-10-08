%% fig_gain07_heading_glitches
% Where are the "blips" in the fictrac heading traces of
% .data\epg_7f_20261007.mat? Per trial: frame-to-frame heading change at the
% fictrac rate, frames above jump_thresh flagged, classified as SPIKE (trace
% returns within spike_win frames) or STEP (stays offset). Gallery of the
% worst trials with the glitches marked, plus a table.

repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(repoRootDir, 'claude'));
S = load(fullfile(repoRootDir, '.data', 'epg_7f_20261007.mat'), 'all_data'); all_data = S.all_data; clear S
exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports', 'epg_7f_gain07');

jump_thresh = 0.3;   % rad per fictrac frame (~17 deg/frame = ~1000 deg/s at 60 Hz): not a real turn
spike_win   = 30;    % frames (0.5 s): a glitch that returns within this is a spike, else a step

n = numel(all_data);
T = table('Size', [n 8], 'VariableTypes', {'double','double','string','double','double','double','double','double'}, ...
    'VariableNames', {'fly','trial','cond','n_jumps','n_spikes','n_steps','max_jump_deg','p999_jump_deg'});
G = cell(n, 1);
for i = 1:n
    d = all_data(i); h = d.ft.heading_raw(:)'; t = d.ft.xf(:)';
    dh = [0 diff(h)];
    bad = find(abs(dh) > jump_thresh);
    % classify: look for a return jump of opposite sign within spike_win
    isSpike = false(size(bad));
    for k = 1:numel(bad)
        j = bad(k); w = j+1:min(j+spike_win, numel(dh));
        isSpike(k) = any(abs(dh(w)) > jump_thresh & sign(dh(w)) == -sign(dh(j)));
    end
    T.fly(i) = d.meta.fly_num; T.trial(i) = d.meta.trial_num; T.cond(i) = string(char("bar" + "")); if d.ft.dark, T.cond(i) = "dark"; else, T.cond(i) = "bar"; end
    T.n_jumps(i) = numel(bad); T.n_spikes(i) = nnz(isSpike); T.n_steps(i) = nnz(~isSpike);
    T.max_jump_deg(i) = round(rad2deg(max(abs(dh)))); T.p999_jump_deg(i) = round(rad2deg(prctile(abs(dh), 99.9)), 1);
    G{i} = struct('t', t, 'h', h, 'dh', dh, 'bad', bad, 'isSpike', isSpike, 'name', sprintf('fly %d t%d (%s)', d.meta.fly_num, d.meta.trial_num, T.cond(i)));
end
disp(T)
fprintf('%d/%d trials have at least one jump > %.2f rad/frame; total %d jumps (%d spikes, %d steps)\n', nnz(T.n_jumps > 0), n, jump_thresh, sum(T.n_jumps), sum(T.n_spikes), sum(T.n_steps));

%% gallery: the 8 trials with the most jumps
[~, ord] = sort(T.n_jumps, 'descend'); show = ord(1:min(8, n));
fig = figure('Visible', 'off', 'Position', [50 50 2000 1400], 'Color', 'w');
for k = 1:numel(show)
    g = G{show(k)};
    ax = subplot(numel(show), 2, 2*k-1); hold on
    plot(g.t, rad2deg(g.h), 'k', 'LineWidth', 0.8)
    if ~isempty(g.bad)
        plot(g.t(g.bad(g.isSpike)), rad2deg(g.h(g.bad(g.isSpike))), 'ro', 'MarkerSize', 5)
        plot(g.t(g.bad(~g.isSpike)), rad2deg(g.h(g.bad(~g.isSpike))), 'bs', 'MarkerSize', 6, 'LineWidth', 1.2)
    end
    ylabel('heading (deg, unwrapped)'); title(sprintf('%s: %d jumps > %.0f deg/frame (red o = spike, blue square = step)', g.name, numel(g.bad), rad2deg(jump_thresh)))
    if k == numel(show), xlabel('fictrac time (s)'); end
    % zoom on the largest glitch
    subplot(numel(show), 2, 2*k); hold on
    if ~isempty(g.bad)
        [~, im] = max(abs(g.dh(g.bad))); j = g.bad(im); w = max(j-120,1):min(j+120, numel(g.h));
        plot(g.t(w), rad2deg(g.h(w)), 'k.-', 'MarkerSize', 6); plot(g.t(j), rad2deg(g.h(j)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5)
        title(sprintf('zoom +/-2 s around the largest jump (%.0f deg in one frame)', rad2deg(abs(g.dh(j)))))
    else
        title('no jumps')
    end
end
pb_export_png(fig, fullfile(exportDir, '_heading_glitches.png'), 110);
writetable(T, fullfile(exportDir, '_heading_glitches.csv'));
fprintf('done\n');
