%% overshoot_bar_stim_led_bleed_check
% Is the stim-LED-locked signal in the overshoot bar-stim trials real PB
% activity or light bleeding into the PMT? overshoot_bar_stim_script.m
% subtracts the per-volume mean of the pixels OUTSIDE the PB mask from every
% frame, which removes a spatially uniform additive offset -- but not a
% gradient, a multiplicative (gain) change, or anything else with spatial
% structure. So here we look at the RAW (registered, un-subtracted) data:
%
%   1. per-volume mean fluorescence inside the mask (PB) and outside it
%      (background), over the hold with the LED pulses marked;
%   2. LED-onset-triggered averages of both, and of their difference, in
%      absolute units and as a fraction of the pre-onset baseline;
%   3. the raw LED-on minus LED-off diff image over the WHOLE field, and the
%      same image restricted to outside the mask with its own colour scale,
%      plus row/column profiles of the outside-mask diff. Light artefacts
%      show up as structure outside the PB (gradients, hot corners, the
%      whole field stepping); neural signal does not.
%
% Needs imgData_sum_reg.mat, mask.mat and bump_results.mat in each trial
% folder (all written by overshoot_bar_stim_script.m). Output per trial:
% <trialDir>\ledBleedCheck.png/.fig and a copy in ugly_figures\exports\, and
% a summary table printed to the command window.

%% 0. trials
trialDirs = { ...
    'Z:\noah_np123\Data\flyg\overshoot\20260930-1_epg_syt8s_lpsp_cschrimson\trial_002\', ...
    'Z:\noah_np123\Data\flyg\overshoot\20260930-2_epg_syt8s_lpsp_cschrimson\trial_001\', ...
    'Z:\noah_np123\Data\flyg\overshoot\20260930-2_epg_syt8s_lpsp_cschrimson\trial_002\', ...
    'Z:\noah_np123\Data\flyg\overshoot\20261001-1_epg_syt8s_cyoTM2\trial_001\'};
if exist('olb_trialDirs', 'var') && ~isempty(olb_trialDirs), trialDirs = olb_trialDirs; end
trig_win = [-5 10]; % s around LED onset for the triggered averages

thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'ugly_figures', 'exports');
if ~exist(exportDir, 'dir'), mkdir(exportDir); end

summary = table('Size', [0 9], 'VariableTypes', [{'string'} repmat({'double'}, 1, 8)], ...
    'VariableNames', {'trial', 'bg_off', 'bg_on', 'bg_pct', 'pb_off', 'pb_on', 'pb_pct', 'pbMinusBg_pct', 'outside_diff_p99_over_inside'});

%% 1. loop over trials
for iT = 1:numel(trialDirs)
    trialDir = trialDirs{iT};
    parts = strsplit(strtrim(trialDir), filesep); parts = parts(~cellfun(@isempty, parts));
    trialName = [parts{end-1} filesep parts{end}];
    fprintf('\n[%s]\n', trialName);

    B = load(fullfile(trialDir, 'bump_results.mat'), 'bump'); b = B.bump;
    M = load(fullfile(trialDir, 'mask.mat'), 'mask'); mask = logical(M.mask);
    S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    img = S.imgData_sum_reg; clear S
    nVol = min(size(img, 3), numel(b.t_vol));
    img = img(:, :, 1:nVol);
    t = b.t_vol(1:nVol); led = logical(b.led(1:nVol)); inHold = logical(b.inHold(1:nVol)); bad = logical(b.badVol(1:nVol));
    if isempty(b.led_on), fprintf('  no LED pulses -- skipped\n'); continue; end

    % raw per-volume means inside / outside the mask (no background subtraction)
    flat = reshape(img, [], nVol);
    pb = mean(flat(mask(:), :), 1);
    bg = mean(flat(~mask(:), :), 1);
    clear flat

    % LED-on / LED-off volumes, restricted to the hold when the pulses are in it
    pulseInHold = arrayfun(@(k) any(inHold & t >= b.led_on(k) & t <= b.led_off(k)), 1:numel(b.led_on));
    if nnz(pulseInHold) >= 3
        useP = find(pulseInHold); span = [b.led_on(useP(1)) b.led_off(useP(end))];
        onIdx  = led & inHold & ~bad;
        offIdx = ~led & inHold & ~bad & t >= span(1) & t <= span(2);
    else
        useP = 1:numel(b.led_on); span = [b.led_on(1) b.led_off(end)];
        onIdx  = led & ~bad;
        offIdx = ~led & ~bad & t >= span(1) & t <= span(2);
    end
    bg_on = mean(bg(onIdx)); bg_off = mean(bg(offIdx));
    pb_on = mean(pb(onIdx)); pb_off = mean(pb(offIdx));
    d_on = mean(pb(onIdx) - bg(onIdx)); d_off = mean(pb(offIdx) - bg(offIdx));
    fprintf('  background (outside mask): LED off %.0f, LED on %.0f -> %+.2f%%\n', bg_off, bg_on, 100 * (bg_on - bg_off) / bg_off);
    fprintf('  PB (inside mask):          LED off %.0f, LED on %.0f -> %+.2f%%\n', pb_off, pb_on, 100 * (pb_on - pb_off) / pb_off);
    fprintf('  PB - background:           LED off %.0f, LED on %.0f -> %+.2f%%\n', d_off, d_on, 100 * (d_on - d_off) / d_off);

    % LED-onset-triggered averages
    dt = median(diff(t)); lagVec = trig_win(1):dt:trig_win(2);
    trigBg = nan(numel(useP), numel(lagVec)); trigPb = trigBg;
    for k = 1:numel(useP)
        tk = b.led_on(useP(k)) + lagVec;
        okk = tk >= t(1) & tk <= t(end);
        trigBg(k, okk) = interp1(t, bg, tk(okk)); trigPb(k, okk) = interp1(t, pb, tk(okk));
    end
    base = lagVec < 0;
    trigBgN = 100 * (trigBg ./ mean(trigBg(:, base), 2, 'omitnan') - 1);
    trigPbN = 100 * (trigPb ./ mean(trigPb(:, base), 2, 'omitnan') - 1);
    trigD = trigPb - trigBg; trigDN = 100 * (trigD ./ mean(trigD(:, base), 2, 'omitnan') - 1);

    % raw diff image, whole field
    diffRaw = mean(img(:, :, onIdx), 3) - mean(img(:, :, offIdx), 3);
    meanImg = mean(img, 3);
    insideP99 = prctile(abs(diffRaw(mask)), 99); outsideP99 = prctile(abs(diffRaw(~mask)), 99);
    fprintf('  raw on-off diff: 99th pct |diff| inside mask = %.0f, outside mask = %.0f (ratio %.2f)\n', insideP99, outsideP99, outsideP99 / insideP99);
    summary = [summary; {string(trialName), bg_off, bg_on, 100*(bg_on-bg_off)/bg_off, pb_off, pb_on, 100*(pb_on-pb_off)/pb_off, 100*(d_on-d_off)/d_off, outsideP99/insideP99}]; %#ok<AGROW>

    % ---- figure ----
    figure(20); clf; set(gcf, 'Position', [40 40 1700 950], 'Color', 'w')
    rwb = olb_redwhiteblue(256);

    ax1 = subplot(3, 3, 1:2); hold on
    yl = [min([pb bg]) max([pb bg])]; yl = yl + [-0.05 0.05] * diff(yl);
    for k = 1:numel(b.led_on)
        patch([b.led_on(k) b.led_off(k) b.led_off(k) b.led_on(k)], yl([1 1 2 2]), [1 0.8 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'off');
    end
    plot(t, pb, 'k', 'DisplayName', 'inside mask (PB), raw')
    plot(t, bg, 'Color', [0.2 0.5 0.9], 'DisplayName', 'outside mask (background), raw')
    xlim([span(1) - 30, span(2) + 30]); ylim(yl); ylabel('mean F (raw units)')
    legend('Location', 'northeastoutside'); title(sprintf('%s: raw per-volume mean F over the LED block (red = LED on)', trialName), 'Interpreter', 'none')

    ax2 = subplot(3, 3, 3); hold on
    patch([0 5 5 0], [-100 -100 100 100], [1 0.8 0.8], 'EdgeColor', 'none', 'HandleVisibility', 'off')
    plot(lagVec, mean(trigBgN, 1, 'omitnan'), 'Color', [0.2 0.5 0.9], 'LineWidth', 1.5, 'DisplayName', 'background')
    plot(lagVec, mean(trigPbN, 1, 'omitnan'), 'k', 'LineWidth', 1.5, 'DisplayName', 'PB')
    plot(lagVec, mean(trigDN, 1, 'omitnan'), 'r', 'LineWidth', 1.5, 'DisplayName', 'PB - background')
    yline(0, 'k:', 'HandleVisibility', 'off')
    allN = [mean(trigBgN, 1, 'omitnan') mean(trigPbN, 1, 'omitnan') mean(trigDN, 1, 'omitnan')];
    ylim([min(allN) max(allN)] + [-1 1] * max(1, 0.1 * range(allN)))
    xlabel('time from LED onset (s)'); ylabel('% change from pre-onset'); legend('Location', 'best')
    title(sprintf('LED-onset-triggered average (n = %d pulses)', numel(useP)))

    ax3 = subplot(3, 3, 4);
    imagesc(meanImg); axis image; colormap(ax3, bone); colorbar; hold on
    contour(mask, [0.5 0.5], 'r', 'LineWidth', 1); title('mean raw image + mask')

    ax4 = subplot(3, 3, 5);
    cl = max(insideP99, eps);
    imagesc(diffRaw, [-cl cl]); axis image; colormap(ax4, rwb); colorbar; hold on
    contour(mask, [0.5 0.5], 'k', 'LineWidth', 1)
    title(sprintf('RAW on - off, whole field (clim = 99th pct |diff| inside mask = %.0f)', cl))

    ax5 = subplot(3, 3, 6);
    outsideDiff = diffRaw; outsideDiff(mask) = NaN;
    clo = max(outsideP99, eps);
    hIm = imagesc(outsideDiff, [-clo clo]); set(hIm, 'AlphaData', ~mask); axis image; colormap(ax5, rwb); colorbar; hold on
    set(ax5, 'Color', [0.6 0.6 0.6]); contour(mask, [0.5 0.5], 'k', 'LineWidth', 1)
    title(sprintf('RAW on - off, OUTSIDE mask only (own clim = %.0f; %.0f%% of inside)', clo, 100 * clo / cl))

    ax6 = subplot(3, 3, 7); hold on
    plot(mean(outsideDiff, 2, 'omitnan'), 1:size(diffRaw, 1), 'b', 'LineWidth', 1.5)
    xline(0, 'k:'); set(gca, 'YDir', 'reverse'); ylabel('row'); xlabel('mean outside-mask diff'); title('row profile of outside diff'); grid on

    ax7 = subplot(3, 3, 8); hold on
    plot(1:size(diffRaw, 2), mean(outsideDiff, 1, 'omitnan'), 'b', 'LineWidth', 1.5)
    yline(0, 'k:'); xlabel('column'); ylabel('mean outside-mask diff'); title('column profile of outside diff'); grid on

    ax8 = subplot(3, 3, 9); hold on
    histogram(diffRaw(~mask), 80, 'Normalization', 'pdf', 'FaceColor', [0.2 0.5 0.9], 'EdgeColor', 'none', 'DisplayName', 'outside mask')
    histogram(diffRaw(mask), 80, 'Normalization', 'pdf', 'FaceColor', [0 0 0], 'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'inside mask')
    xline(0, 'k:', 'HandleVisibility', 'off'); xlabel('raw on - off per pixel'); ylabel('pdf'); legend; title('pixel-wise diff distributions')

    sgtitle(sprintf('%s: is the LED-locked signal light bleed-through? bg: %+.2f%%, PB: %+.2f%%, PB-bg: %+.2f%% (LED on vs off)', ...
        trialName, 100*(bg_on-bg_off)/bg_off, 100*(pb_on-pb_off)/pb_off, 100*(d_on-d_off)/d_off), 'Interpreter', 'none')

    drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off');
    for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
    drawnow
    outBase = 'ledBleedCheck';
    exportgraphics(gcf, fullfile(trialDir, [outBase '.png']), 'Resolution', 150);
    exportgraphics(gcf, fullfile(exportDir, [strrep(trialName, filesep, '_') '_' outBase '.png']), 'Resolution', 150);
    vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(trialDir, [outBase '.fig']), 'compact'); set(gcf, 'Visible', vis);
    fprintf('  saved %s\n', fullfile(trialDir, [outBase '.png']));
    clear img
end

fprintf('\n==== summary (LED on vs off, raw units; pct = percent change) ====\n');
disp(summary)

%% local functions
function cmap = olb_redwhiteblue(n)
half = floor(n / 2);
up = linspace(0, 1, half)';
cmap = [[up up ones(half, 1)]; [ones(n - half, 1) flipud(linspace(0, 1, n - half)') flipud(linspace(0, 1, n - half)')]];
end
