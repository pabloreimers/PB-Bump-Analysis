%% overshoot_ttx_occupancy_vs_phase
% Fly 20261005-1. Does the biphasic LED response under TTX (trial_006:
% early 0-2 s excitation in some glomeruli, inhibition in others; see
% overshoot_ttx_early_vs_late.m) line up with where the bump sat in the
% saline trial (trial_005) -- during the LED pulses, or in the LED-off gaps
% between them?
%
% Per glomerulus: TTX early (0-2 s) and late (2-5 s) dF/F (baseline = last
% 1.5 s before onset), and early - late. trial_005: rho-weighted bump
% occupancy of each glomerulus's PVA angle during (a) LED-on, (b) LED-off
% gaps between pulses, (c) the pre-LED part of the hold, (d) the post-LED
% part of the hold; plus per-glomerulus mean z in (a) and (b).
% Correlations with a circular-shift permutation null (profile shifted
% within each hemisphere). Group test: glomeruli split by TTX early
% response (excited / inhibited / neutral) vs occupancy (Kruskal-Wallis).
% Output: ttxOccupancyVsPhase.png/.fig in the trial_006 folder + exports.

%% 0. inputs
flyDir = 'Z:\noah_np123\Data\flyg\overshoot\20261005-1_epg_syt8s_lpsp_cschrimson\';
ttxDir = fullfile(flyDir, 'trial_006\'); salDir = fullfile(flyDir, 'trial_005\');
% Overrides: evp_ttxDir = trial whose per-glomerulus LED response is mapped
% ("ttx" here just means the response-map trial), evp_salDir = trial whose
% bump occupancy is compared against it, evp_tag = export-file prefix.
if exist('evp_ttxDir', 'var') && ~isempty(evp_ttxDir), ttxDir = evp_ttxDir; end
if exist('evp_salDir', 'var') && ~isempty(evp_salDir), salDir = evp_salDir; end
outTag = '20261005-1_trial006_vs_005'; if exist('evp_tag', 'var') && ~isempty(evp_tag), outTag = evp_tag; end
early_win = [0 2]; late_win = [2 5]; base_win = [-1.5 0];
nPerm = 5000; rng(1)
% Which part of the trial_005 LED block to use. The bump changes behaviour
% about halfway through (rotating, then collapsing/parking), so 'second'
% restricts the LED-on / LED-off-gap occupancy to t >= the midpoint of the
% LED span; 'first' the reverse; 'all' uses the whole block.
block_part = 'second';   % 'all' | 'first' | 'second'
if exist('evp_block_part', 'var') && ~isempty(evp_block_part), block_part = evp_block_part; end
thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'ugly_figures', 'exports'); if ~exist(exportDir, 'dir'), mkdir(exportDir); end

%% 1. TTX early / late per glomerulus
S = load(fullfile(ttxDir, 'bump_results.mat')); b6 = S.bump;
t = b6.t_vol(:)'; f = b6.f_cluster; inHold = logical(b6.inHold(:)'); nG = size(f, 1);
pul = find(arrayfun(@(k) all(inHold(t >= b6.led_on(k) & t <= b6.led_off(k))), 1:numel(b6.led_on)));
E = nan(nG, numel(pul)); L = E;
for k = 1:numel(pul)
    t0 = b6.led_on(pul(k));
    F0 = mean(f(:, t >= t0 + base_win(1) & t < t0 + base_win(2)), 2, 'omitnan');
    E(:, k) = (mean(f(:, t >= t0 + early_win(1) & t < t0 + early_win(2)), 2, 'omitnan') - F0) ./ F0;
    L(:, k) = (mean(f(:, t >= t0 + late_win(1)  & t < t0 + late_win(2)),  2, 'omitnan') - F0) ./ F0;
end
ttxE = mean(E, 2); ttxL = mean(L, 2); ttxEL = ttxE - ttxL;
[~, pE0, ~, stE0] = ttest(E');
grp = zeros(nG, 1); grp(pE0 < 0.05 & stE0.tstat > 0) = 1; grp(pE0 < 0.05 & stE0.tstat < 0) = -1;   % +1 early-excited, -1 early-inhibited, 0 neutral
fprintf('TTX groups by early response (uncorrected p<0.05): excited %s; inhibited %s; neutral %d glomeruli\n', mat2str(find(grp == 1)'), mat2str(find(grp == -1)'), nnz(grp == 0));

%% 2. trial_005 occupancy / activity per glomerulus in the different epochs
S = load(fullfile(salDir, 'bump_results.mat')); b5 = S.bump;
t5 = b5.t_vol(:)'; mu = b5.mu(:)'; rho = b5.rho(:)'; z5 = b5.z_cluster; led5 = logical(b5.led(:)'); inH5 = logical(b5.inHold(:)');
alpha = b5.alpha(:)'; aBins = unique(round(alpha, 6)); nB = numel(aBins); binOf = knnsearch(aBins', alpha');
span = [b5.led_on(1) b5.led_off(end)]; tMid = mean(span);
switch block_part
    case 'all',    inPart = true(size(t5));          partLbl = 'whole LED block';
    case 'first',  inPart = t5 <  tMid;              partLbl = sprintf('1st half of LED block (< %.0f s)', tMid);
    case 'second', inPart = t5 >= tMid;              partLbl = sprintf('2nd half of LED block (>= %.0f s)', tMid);
end
ep.on   = led5 & inH5 & inPart;
ep.gap  = ~led5 & inH5 & t5 > span(1) & t5 < span(2) & inPart;
ep.pre  = inH5 & t5 < span(1);
ep.post = inH5 & t5 > span(2);
ep.hold = inH5;
ep.other = (led5 | (t5 > span(1) & t5 < span(2))) & inH5 & ~inPart;   % the complementary part of the LED block, for reference
fprintf('trial_005: LED block %.0f-%.0f s, using %s for the on/gap occupancy\n', span, partLbl);
epNames = fieldnames(ep);
occ = struct(); zm = struct();
for e = 1:numel(epNames)
    sel = ep.(epNames{e}) & rho >= b5.rho_thresh & ~isnan(mu);
    o = accumarray(knnsearch(aBins', angle(exp(1i*mu(sel)))'), rho(sel), [nB 1])' / sum(rho(sel));
    occ.(epNames{e}) = (o(binOf) * nB)';          % relative to uniform (1 = uniform), per glomerulus
    zm.(epNames{e}) = mean(z5(:, ep.(epNames{e})), 2, 'omitnan');
    fprintf('trial_005 epoch %-4s: %5d volumes, %4d with rho>=%.1f, occupancy range %.2f..%.2f x uniform, mean rho %.2f\n', ...
        epNames{e}, nnz(ep.(epNames{e})), nnz(sel), b5.rho_thresh, min(occ.(epNames{e})), max(occ.(epNames{e})), mean(rho(ep.(epNames{e})), 'omitnan'));
end

%% 3. correlations with circular-shift null
hemi = [ones(1, nG/2) 2*ones(1, nG/2)];
permcorr = @(x, y) deal(corr(x(:), y(:)), mean(abs(arrayfun(@(p) corr(cshift(x(:), hemi, randi(nG/2) - 1), y(:)), 1:nPerm)) >= abs(corr(x(:), y(:)))));
ttxVars = {'LED resp. early (0-2 s)', ttxE; 'LED resp. late (2-5 s)', ttxL; 'LED resp. early - late', ttxEL};
salVars = {['occupancy, LED on (' block_part ' half)'], occ.on; ['occupancy, LED-off gaps (' block_part ' half)'], occ.gap; 'occupancy, hold before LED', occ.pre; 'occupancy, hold after LED', occ.post; ...
           ['mean z, LED on (' block_part ')'], zm.on; ['mean z, LED-off gaps (' block_part ')'], zm.gap; 'occupancy, other half of block', occ.other};
if strcmp(block_part, 'all'), salVars = salVars(1:6, :); end
fprintf('\n%-22s', ''); fprintf('%-28s', salVars{:, 1}); fprintf('\n');
Rm = nan(size(ttxVars, 1), size(salVars, 1)); Pm = Rm;
for i = 1:size(ttxVars, 1)
    fprintf('%-22s', ttxVars{i, 1});
    for j = 1:size(salVars, 1)
        [r, p] = permcorr(salVars{j, 2}, ttxVars{i, 2}); Rm(i, j) = r; Pm(i, j) = p;
        fprintf('r=%+.2f p=%.2f%s', r, p, repmat(' ', 1, 28 - 15));
    end
    fprintf('\n');
end
% group comparison
fprintf('\noccupancy (x uniform) by TTX early-response group:  excited(n=%d) / neutral(n=%d) / inhibited(n=%d)\n', nnz(grp == 1), nnz(grp == 0), nnz(grp == -1));
for j = 1:4
    x = salVars{j, 2};
    pKW = kruskalwallis(x, grp, 'off');
    fprintf('  %-28s %.2f / %.2f / %.2f   (Kruskal-Wallis p = %.2f)\n', salVars{j, 1}, mean(x(grp == 1)), mean(x(grp == 0)), mean(x(grp == -1)), pKW);
end

%% 4. figure
figure(32); clf; set(gcf, 'Position', [30 30 1800 950], 'Color', 'w')
gIdx = 1:nG; cE = [0.85 0.4 0.2]; cI = [0.3 0.3 0.8]; cN = [0.6 0.6 0.6];
gcol = cN .* (grp == 0) + cE .* (grp == 1) + cI .* (grp == -1);
% (a) profiles
ax1 = subplot(3, 4, [1 2]); hold on
yyaxis left
bar(gIdx, 100*ttxE, 'FaceColor', 'flat', 'CData', gcol, 'EdgeColor', 'none', 'BarWidth', 0.7, 'HandleVisibility', 'off')
ylabel('LED resp. early dF/F (%)'); yline(0, 'k:', 'HandleVisibility', 'off')
yyaxis right
plot(gIdx, occ.on, '-o', 'Color', [0.1 0.6 0.1], 'MarkerFaceColor', [0.1 0.6 0.1], 'MarkerSize', 4, 'DisplayName', 'occupancy-trial occupancy, LED on')
plot(gIdx, occ.gap, '--s', 'Color', [0 0 0], 'MarkerSize', 4, 'DisplayName', 'occupancy-trial occupancy, LED-off gaps')
plot(gIdx, occ.pre, ':d', 'Color', [0.5 0.2 0.6], 'MarkerSize', 4, 'DisplayName', 'occupancy-trial occupancy, hold before LED')
if ~strcmp(block_part, 'all'), plot(gIdx, occ.other, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'occupancy-trial occupancy, other half of LED block'); end
yline(1, 'k:', 'HandleVisibility', 'off'); ylabel('bump occupancy (x uniform)')
xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5]); xlabel('glomerulus')
legend('Location', 'northwest', 'FontSize', 7)
title(sprintf('LED response (early 0-2 s) (bars: orange = early-excited, blue = early-inhibited, grey = neutral) vs where the bump was in the occupancy trial -- %s', partLbl), 'FontSize', 9)
% (b) group boxes
ax2 = subplot(3, 4, 3); hold on
G = {occ.on, occ.gap}; gn = {'LED on', 'LED-off gaps'};
xpos = 0;
for j = 1:2
    for g = [1 0 -1]
        xpos = xpos + 1; v = G{j}(grp == g);
        scatter(xpos + 0.15*randn(size(v)), v, 25, cN .* (g == 0) + cE .* (g == 1) + cI .* (g == -1), 'filled')
        plot(xpos + [-0.3 0.3], mean(v) * [1 1], 'k-', 'LineWidth', 2)
    end
    xpos = xpos + 1;
end
yline(1, 'k:'); xticks([2 6]); xticklabels(gn); ylabel('occupancy-trial occupancy (x uniform)')
title(sprintf('by LED-response group (exc / neutral / inh); KW p = %.2f (on), %.2f (gaps)', kruskalwallis(occ.on, grp, 'off'), kruskalwallis(occ.gap, grp, 'off')), 'FontSize', 8)
% (c) heatmap of correlations
ax3 = subplot(3, 4, 4);
imagesc(Rm); colormap(ax3, evp_rwb(256)); clim([-1 1]); colorbar
xticks(1:size(salVars, 1)); xticklabels(salVars(:, 1)); xtickangle(30); yticks(1:3); yticklabels(ttxVars(:, 1)); set(gca, 'FontSize', 7)
for i = 1:3, for j = 1:size(salVars, 1), text(j, i, sprintf('%+.2f\np=%.2f', Rm(i, j), Pm(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 7); end, end
title('Pearson r (circ-shift perm p)', 'FontSize', 8)
% (d-i) scatters: early vs on / gap / pre ; late vs on / gap / pre
k = 4;
for i = [1 2]
    for j = 1:3
        k = k + 1; ax = subplot(3, 4, k); hold on
        x = salVars{j, 2}; y = 100*ttxVars{i, 2};
        scatter(x, y, 36, gcol, 'filled'); text(x + 0.02, y, string(gIdx'), 'FontSize', 6)
        pf = polyfit(x, y, 1); xx = linspace(min(x), max(x), 2); plot(xx, polyval(pf, xx), 'k-')
        xline(1, 'k:'); yline(0, 'k:')
        xlabel(['occ.-trial ' salVars{j, 1} ' (x uniform)']); ylabel([ttxVars{i, 1} ' dF/F (%)'])
        title(sprintf('r = %+.2f, perm p = %.2f', Rm(i, j), Pm(i, j)), 'FontSize', 8)
    end
    k = k + 1; % skip column 4 (used by heatmap above / leave empty)
end
% use the two spare panels for mean-z scatters
ax = subplot(3, 4, 8); hold on
x = zm.on; y = 100*ttxE; scatter(x, y, 36, gcol, 'filled'); text(x + 0.005, y, string(gIdx'), 'FontSize', 6)
pf = polyfit(x, y, 1); xx = linspace(min(x), max(x), 2); plot(xx, polyval(pf, xx), 'k-'); yline(0, 'k:')
xlabel('occupancy-trial mean z, LED on'); ylabel('LED resp. early dF/F (%)'); title(sprintf('r = %+.2f, perm p = %.2f', Rm(1, 5), Pm(1, 5)), 'FontSize', 8)
ax = subplot(3, 4, 12); hold on
x = zm.gap; y = 100*ttxE; scatter(x, y, 36, gcol, 'filled'); text(x + 0.005, y, string(gIdx'), 'FontSize', 6)
pf = polyfit(x, y, 1); xx = linspace(min(x), max(x), 2); plot(xx, polyval(pf, xx), 'k-'); yline(0, 'k:')
xlabel('occupancy-trial mean z, LED-off gaps'); ylabel('LED resp. early dF/F (%)'); title(sprintf('r = %+.2f, perm p = %.2f', Rm(1, 6), Pm(1, 6)), 'FontSize', 8)
sgtitle(sprintf('%s: early / late LED response per glomerulus  vs  %s: bump occupancy during LED pulses and in the LED-off gaps -- %s', b6.trialName, b5.trialName, partLbl), 'FontSize', 11, 'Interpreter', 'none')

drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off')
for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow
base = 'ttxOccupancyVsPhase'; if ~strcmp(block_part, 'all'), base = [base '_' block_part 'Half']; end
exportgraphics(gcf, fullfile(ttxDir, [base '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, [outTag '_' base '.png']), 'Resolution', 150);
vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(ttxDir, [base '.fig']), 'compact'); set(gcf, 'Visible', vis);
fprintf('saved %s.png / .fig\n', fullfile(ttxDir, base));

%% local functions
function xs = cshift(x, hemi, s)
xs = x;
for h = 1:2, idx = find(hemi == h); xs(idx) = circshift(x(idx), s); end
end
function cmap = evp_rwb(n)
half = floor(n / 2); up = linspace(0, 1, half)';
cmap = [[up up ones(half, 1)]; [ones(n - half, 1) flipud(linspace(0, 1, n - half)') flipud(linspace(0, 1, n - half)')]];
end
