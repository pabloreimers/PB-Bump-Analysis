%% overshoot_ttx_early_vs_late
% Is the LED response under TTX (fly 20261005-1, trial_006) biphasic -- some
% excitation in the first ~2 s of the pulse, inhibition in the last 3 s?
% Per glomerulus, each of the 30 pulses is a replicate:
%   early = mean dF/F over 0-2 s after LED onset, late = mean over 2-5 s,
%   baseline = mean F over the last 1.5 s before onset (the 5-s window used
%   in overshoot_ttx_inhibition_map.m contains the rebound from the previous
%   pulse's offset; both baselines are reported).
% Tests: paired t (early vs late) per glomerulus; one-sample t (early vs 0)
% per glomerulus; Bonferroni over 32; plus the same for the pooled PB.
% trial_005 (saline, same protocol) is run through the same code for
% comparison. Output: ttxEarlyVsLate.png/.fig in the trial_006 folder and
% ugly_figures\exports\.

%% 0. inputs
flyDir = 'Z:\noah_np123\Data\flyg\overshoot\20261005-1_epg_syt8s_lpsp_cschrimson\';
trials = {fullfile(flyDir, 'trial_006\'), 'TTX'; fullfile(flyDir, 'trial_005\'), 'saline'};
% Override: evl_trials = {dir1, label1; dir2, label2} -- row 1 is the trial
% analysed in detail, row 2 the comparison; evl_tag names the export file.
if exist('evl_trials', 'var') && ~isempty(evl_trials), trials = evl_trials; end
outTag = '20261005-1_trial006'; if exist('evl_tag', 'var') && ~isempty(evl_tag), outTag = evl_tag; end
early_win = [0 2]; late_win = [2 5];
base_short = [-1.5 0]; base_long = [-5 0];
thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'ugly_figures', 'exports'); if ~exist(exportDir, 'dir'), mkdir(exportDir); end

R = struct();
for iT = 1:size(trials, 1)
    S = load(fullfile(trials{iT, 1}, 'bump_results.mat')); b = S.bump;
    t = b.t_vol(:)'; f = b.f_cluster; inHold = logical(b.inHold(:)'); nG = size(f, 1);
    pul = find(arrayfun(@(k) all(inHold(t >= b.led_on(k) & t <= b.led_off(k))), 1:numel(b.led_on)));
    nP = numel(pul);
    dt = median(diff(t)); lags = -5:dt:10;
    trigS = nan(nG, numel(lags), nP); trigL = trigS;   % dF/F with short / long baseline
    early = nan(nG, nP, 2); late = early;              % (:,:,1) short baseline, (:,:,2) long
    for k = 1:nP
        t0 = b.led_on(pul(k));
        F0s = mean(f(:, t >= t0 + base_short(1) & t < t0 + base_short(2)), 2, 'omitnan');
        F0l = mean(f(:, t >= t0 + base_long(1)  & t < t0 + base_long(2)),  2, 'omitnan');
        tk = t0 + lags; okk = tk >= t(1) & tk <= t(end);
        fi = interp1(t, f', tk(okk))';
        trigS(:, okk, k) = (fi - F0s) ./ F0s; trigL(:, okk, k) = (fi - F0l) ./ F0l;
        e = t >= t0 + early_win(1) & t < t0 + early_win(2); l = t >= t0 + late_win(1) & t < t0 + late_win(2);
        early(:, k, 1) = (mean(f(:, e), 2, 'omitnan') - F0s) ./ F0s; late(:, k, 1) = (mean(f(:, l), 2, 'omitnan') - F0s) ./ F0s;
        early(:, k, 2) = (mean(f(:, e), 2, 'omitnan') - F0l) ./ F0l; late(:, k, 2) = (mean(f(:, l), 2, 'omitnan') - F0l) ./ F0l;
    end
    R(iT).label = trials{iT, 2}; R(iT).name = b.trialName; R(iT).nG = nG; R(iT).nP = nP; R(iT).lags = lags;
    R(iT).trigS = trigS; R(iT).trigL = trigL; R(iT).early = early; R(iT).late = late;
    for bl = 1:2
        E = early(:, :, bl); L = late(:, :, bl);
        [~, pEL, ~, stEL] = ttest(E', L');           % paired, per glomerulus
        [~, pE0, ~, stE0] = ttest(E');               % early vs 0
        [~, pL0, ~, stL0] = ttest(L');               % late vs 0
        [~, pPool] = ttest(mean(E, 1, 'omitnan')', mean(L, 1, 'omitnan')');  % pooled PB, paired
        [~, pPoolE] = ttest(mean(E, 1, 'omitnan')');
        R(iT).stats(bl) = struct('pEL', pEL, 'tEL', stEL.tstat, 'pE0', pE0, 'tE0', stE0.tstat, 'pL0', pL0, 'tL0', stL0.tstat, ...
            'mE', mean(E, 2), 'seE', std(E, 0, 2)/sqrt(nP), 'mL', mean(L, 2), 'seL', std(L, 0, 2)/sqrt(nP), 'pPool', pPool, 'pPoolE', pPoolE);
    end
    if iT == 1, blNames = {'baseline = last 1.5 s before onset', 'baseline = 5 s before onset'};
        for bl = 1:2
            st = R(iT).stats(bl); bonf = 0.05 / nG;
            fprintf('\n== %s (%s), %d pulses, %s ==\n', b.trialName, trials{iT, 2}, nP, blNames{bl});
            fprintf('pooled PB: early (0-2 s) %+.1f%%, late (2-5 s) %+.1f%% -> paired p = %.2g; early vs 0: p = %.2g\n', ...
                100*mean(st.mE, 'omitnan'), 100*mean(st.mL, 'omitnan'), st.pPool, st.pPoolE);
            fprintf('per glomerulus, early > late (paired, Bonferroni p<%.4f): %d/%d  -> %s\n', bonf, nnz(st.pEL < bonf & st.tEL > 0), nG, mat2str(find(st.pEL < bonf & st.tEL > 0)));
            fprintf('per glomerulus, early > 0 (Bonferroni): %d/%d -> %s;  early < 0: %d/%d -> %s\n', ...
                nnz(st.pE0 < bonf & st.tE0 > 0), nG, mat2str(find(st.pE0 < bonf & st.tE0 > 0)), nnz(st.pE0 < bonf & st.tE0 < 0), nG, mat2str(find(st.pE0 < bonf & st.tE0 < 0)));
            fprintf('per glomerulus, early > 0 (uncorrected p<0.05): %d/%d -> %s\n', nnz(st.pE0 < 0.05 & st.tE0 > 0), nG, mat2str(find(st.pE0 < 0.05 & st.tE0 > 0)));
            fprintf('per glomerulus, late < 0 (Bonferroni): %d/%d\n', nnz(st.pL0 < bonf & st.tL0 < 0), nG);
            [~, iMax] = sort(st.mE, 'descend');
            fprintf('largest early responses: glomeruli %s = %s %%\n', mat2str(iMax(1:6)'), mat2str(round(100*st.mE(iMax(1:6))', 1)));
        end
    else
        st = R(iT).stats(1);
        fprintf('\n== %s (%s) for comparison, short baseline: pooled early %+.1f%%, late %+.1f%% (paired p = %.2g); early>0 Bonf: %d/%d, late<0 Bonf: %d/%d\n', ...
            b.trialName, trials{iT, 2}, 100*mean(st.mE, 'omitnan'), 100*mean(st.mL, 'omitnan'), st.pPool, nnz(st.pE0 < 0.05/nG & st.tE0 > 0), nG, nnz(st.pL0 < 0.05/nG & st.tL0 < 0), nG);
    end
end

%% figure
iTTX = 1; iSal = 2; lab1 = R(1).label; lab2 = R(2).label;   % row 1 = analysed trial, row 2 = comparison
st = R(iTTX).stats(1); nG = R(iTTX).nG; gIdx = 1:nG; bonf = 0.05 / nG;
figure(31); clf; set(gcf, 'Position', [30 30 1800 1000], 'Color', 'w')
rwb = evl_rwb(256);
% (a) pulse-triggered heatmap, short baseline, with windows
ax1 = subplot(3, 3, 1);
M = 100*mean(R(iTTX).trigS, 3, 'omitnan'); imagesc(R(iTTX).lags, gIdx, M); colormap(ax1, rwb); clim([-15 15]); colorbar; hold on
xline([0 2 5], 'k'); yline(nG/2 + 0.5, 'k:'); xlabel('time from LED onset (s)'); ylabel('glomerulus')
title(sprintf('%s (%s): pulse-triggered dF/F (%%), baseline -1.5..0 s; lines at 0 / 2 / 5 s', R(iTTX).name, lab1), 'Interpreter', 'none', 'FontSize', 8)
% (b) early vs late per glomerulus
ax2 = subplot(3, 3, [2 3]); hold on
bar(gIdx - 0.2, 100*st.mE, 0.4, 'FaceColor', [0.85 0.4 0.2], 'EdgeColor', 'none', 'DisplayName', 'early (0-2 s)')
bar(gIdx + 0.2, 100*st.mL, 0.4, 'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'none', 'DisplayName', 'late (2-5 s)')
errorbar(gIdx - 0.2, 100*st.mE, 100*st.seE, 'k.', 'LineStyle', 'none', 'CapSize', 2, 'HandleVisibility', 'off')
errorbar(gIdx + 0.2, 100*st.mL, 100*st.seL, 'k.', 'LineStyle', 'none', 'CapSize', 2, 'HandleVisibility', 'off')
sigEL = st.pEL < bonf & st.tEL > 0; sigE0 = st.pE0 < bonf & st.tE0 > 0;
plot(gIdx(sigEL), 100*max(st.mE(sigEL), 0) + 100*st.seE(sigEL) + 1, 'kv', 'MarkerFaceColor', 'k', 'MarkerSize', 4, 'DisplayName', 'early > late (Bonf.)')
plot(gIdx(sigE0), 100*st.mE(sigE0) + 100*st.seE(sigE0) + 2.5, 'r*', 'MarkerSize', 6, 'DisplayName', 'early > 0 (Bonf.)')
yline(0, 'k:', 'HandleVisibility', 'off'); xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5])
xlabel('glomerulus'); ylabel('dF/F (%)'); legend('Location', 'southwest', 'FontSize', 7)
title(sprintf('early vs late, mean +- SEM over %d pulses: pooled early %+.1f%% vs late %+.1f%% (paired p = %.1g); early>late in %d/32, early>0 in %d/32 glomeruli (Bonferroni)', ...
    R(iTTX).nP, 100*mean(st.mE, 'omitnan'), 100*mean(st.mL, 'omitnan'), st.pPool, nnz(sigEL), nnz(sigE0)), 'FontSize', 8)
% (c) t-statistics
ax3 = subplot(3, 3, 4); hold on
plot(gIdx, st.tE0, '-o', 'Color', [0.85 0.4 0.2], 'MarkerFaceColor', [0.85 0.4 0.2], 'MarkerSize', 4, 'DisplayName', 'early vs 0')
plot(gIdx, st.tL0, '-o', 'Color', [0.3 0.3 0.8], 'MarkerFaceColor', [0.3 0.3 0.8], 'MarkerSize', 4, 'DisplayName', 'late vs 0')
plot(gIdx, st.tEL, '-s', 'Color', 'k', 'MarkerSize', 4, 'DisplayName', 'early vs late (paired)')
tcrit = tinv(1 - bonf/2, R(iTTX).nP - 1); yline([-tcrit tcrit], 'r--', 'HandleVisibility', 'off'); yline(0, 'k:', 'HandleVisibility', 'off')
xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5]); xlabel('glomerulus'); ylabel('t statistic (29 df)')
legend('Location', 'best', 'FontSize', 7); title('t statistics per glomerulus; red dashed = Bonferroni threshold', 'FontSize', 8)
% (d) time courses: glomeruli ranked by early response
ax4 = subplot(3, 3, 5); hold on
[~, ord] = sort(st.mE, 'descend', 'MissingPlacement', 'last'); ordV = ord(~isnan(st.mE(ord))); grpTop = ordV(1:8); grpBot = ordV(end-7:end);
mTop = 100*squeeze(mean(mean(R(iTTX).trigS(grpTop, :, :), 1, 'omitnan'), 3, 'omitnan'));
mBot = 100*squeeze(mean(mean(R(iTTX).trigS(grpBot, :, :), 1, 'omitnan'), 3, 'omitnan'));
mAll = 100*squeeze(mean(mean(R(iTTX).trigS, 1, 'omitnan'), 3, 'omitnan'));
seTop = 100*squeeze(std(mean(R(iTTX).trigS(grpTop, :, :), 1, 'omitnan'), 0, 3, 'omitnan')) / sqrt(R(iTTX).nP);
seBot = 100*squeeze(std(mean(R(iTTX).trigS(grpBot, :, :), 1, 'omitnan'), 0, 3, 'omitnan')) / sqrt(R(iTTX).nP);
patch([0 5 5 0], [-30 -30 30 30], [1 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off')
lg = R(iTTX).lags;
fill([lg fliplr(lg)], [mTop + seTop, fliplr(mTop - seTop)], [0.85 0.4 0.2], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
fill([lg fliplr(lg)], [mBot + seBot, fliplr(mBot - seBot)], [0.3 0.3 0.8], 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off')
plot(lg, mTop, 'Color', [0.85 0.4 0.2], 'LineWidth', 1.5, 'DisplayName', sprintf('8 glomeruli with largest early response (%s)', mat2str(sort(grpTop)')))
plot(lg, mBot, 'Color', [0.3 0.3 0.8], 'LineWidth', 1.5, 'DisplayName', sprintf('8 with smallest early response (%s)', mat2str(sort(grpBot)')))
plot(lg, mAll, 'k', 'LineWidth', 1, 'DisplayName', 'all 32')
yline(0, 'k:', 'HandleVisibility', 'off'); xlim([lg(1) lg(end)]); ylim([-18 10])
xlabel('time from LED onset (s)'); ylabel('dF/F (%)'); legend('Location', 'southwest', 'FontSize', 7)
title(sprintf('%s: time course, mean +- SEM over pulses', lab1), 'FontSize', 8)
% (e) same for saline trial_005
ax5 = subplot(3, 3, 6); hold on
stS = R(iSal).stats(1);
patch([0 5 5 0], [-30 -30 30 30], [1 0.85 0.85], 'EdgeColor', 'none', 'HandleVisibility', 'off')
mTopS = 100*squeeze(mean(mean(R(iSal).trigS(grpTop, :, :), 1, 'omitnan'), 3, 'omitnan'));
mBotS = 100*squeeze(mean(mean(R(iSal).trigS(grpBot, :, :), 1, 'omitnan'), 3, 'omitnan'));
mAllS = 100*squeeze(mean(mean(R(iSal).trigS, 1, 'omitnan'), 3, 'omitnan'));
plot(R(iSal).lags, mTopS, 'Color', [0.85 0.4 0.2], 'LineWidth', 1.5, 'DisplayName', 'same "early" glomeruli as TTX')
plot(R(iSal).lags, mBotS, 'Color', [0.3 0.3 0.8], 'LineWidth', 1.5, 'DisplayName', 'same "late" glomeruli')
plot(R(iSal).lags, mAllS, 'k', 'LineWidth', 1, 'DisplayName', 'all 32')
yline(0, 'k:', 'HandleVisibility', 'off'); xlim([lg(1) lg(end)]); ylim([-18 25])
xlabel('time from LED onset (s)'); ylabel('dF/F (%)'); legend('Location', 'northeast', 'FontSize', 7)
title(sprintf('%s (%s) for comparison: pooled early %+.1f%%, late %+.1f%%', R(iSal).name, lab2, 100*mean(stS.mE, 'omitnan'), 100*mean(stS.mL, 'omitnan')), 'Interpreter', 'none', 'FontSize', 8)
% (f) per-pulse early response of the top-8 group: is the excitation consistent across pulses, or a few outliers?
ax6 = subplot(3, 3, 7); hold on
Etop = 100*squeeze(mean(R(iTTX).early(grpTop, :, 1), 1, 'omitnan')); Ltop = 100*squeeze(mean(R(iTTX).late(grpTop, :, 1), 1, 'omitnan'));
plot(1:R(iTTX).nP, Etop, '-o', 'Color', [0.85 0.4 0.2], 'MarkerFaceColor', [0.85 0.4 0.2], 'MarkerSize', 4, 'DisplayName', 'early, top-8 glomeruli')
plot(1:R(iTTX).nP, Ltop, '-o', 'Color', [0.3 0.3 0.8], 'MarkerFaceColor', [0.3 0.3 0.8], 'MarkerSize', 4, 'DisplayName', 'late, top-8 glomeruli')
yline(0, 'k:', 'HandleVisibility', 'off'); xlabel('pulse #'); ylabel('dF/F (%)'); legend('Location', 'best', 'FontSize', 7)
title(sprintf('per pulse: early > 0 in %d/%d pulses, early > late in %d/%d', nnz(Etop > 0), R(iTTX).nP, nnz(Etop > Ltop), R(iTTX).nP), 'FontSize', 8)
% (g) early vs late scatter per glomerulus
ax7 = subplot(3, 3, 8); hold on
scatter(100*st.mL, 100*st.mE, 36, gIdx, 'filled'); colormap(ax7, hsv(nG)); text(100*st.mL + 0.2, 100*st.mE, string(gIdx'), 'FontSize', 6)
plot([-16 2], [-16 2], 'k--'); xline(0, 'k:'); yline(0, 'k:'); xlabel('late dF/F (%)'); ylabel('early dF/F (%)')
[rEL, pREL] = corr(st.mL, st.mE);
title(sprintf('early vs late across glomeruli: r = %.2f (p = %.2g)', rEL, pREL), 'FontSize', 8)
% (h) early response in image space (TTX)
ax8 = subplot(3, 3, 9);
B6 = load(fullfile(trials{iTTX, 1}, 'bump_results.mat')); cI = B6.bump.clusterIdx; mapE = nan(size(cI));
for g = 1:nG, mapE(cI == g) = 100*st.mE(g); end
hI = imagesc(mapE); set(hI, 'AlphaData', ~isnan(mapE)); set(ax8, 'Color', [0.92 0.92 0.92]); axis image; colormap(ax8, rwb); clim([-8 8]); colorbar
title('TTX: early (0-2 s) dF/F per glomerulus (%)', 'FontSize', 8)
sgtitle(sprintf('%s (%s): is the LED response biphasic? early (0-2 s) vs late (2-5 s) of each 5-s pulse', R(iTTX).name, lab1), 'FontSize', 11, 'Interpreter', 'none')

drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off')
for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow
base = 'ttxEarlyVsLate';
exportgraphics(gcf, fullfile(trials{iTTX, 1}, [base '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, [outTag '_' base '.png']), 'Resolution', 150);
vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(trials{iTTX, 1}, [base '.fig']), 'compact'); set(gcf, 'Visible', vis);
fprintf('saved %s.png / .fig\n', fullfile(trials{iTTX, 1}, base));

function cmap = evl_rwb(n)
half = floor(n / 2); up = linspace(0, 1, half)';
cmap = [[up up ones(half, 1)]; [ones(n - half, 1) flipud(linspace(0, 1, n - half)') flipud(linspace(0, 1, n - half)')]];
end
