%% overshoot_glom_delta_script
% Quantify the per-glomerulus stim response from the "overshoot" experiment
% (Z:\noah_np123\Data\flyg\overshoot\<fly>\trial_NNN_<intensity>_stim_<drug>\)
% as ONE LINE PER TRIAL, then fold the two PB hemispheres onto each other so
% every trial reduces to a single peak + trough across one period of the PB.
%
% Builds directly on overshoot_drug_stim_script.m (same trial discovery,
% same mask.mat, same 40-glomerulus skeleton clustering, same per-glomerulus
% dF/F definition as its sections 7-9). Run that script's sections 1-3 on a
% fly first so mask.mat + imgData_sum_reg.mat exist; this script only reads.
%
% What is computed, per fly:
%   1. delta(g), g = 1..40: per-glomerulus dF/F CHANGE during the LED pulse.
%      dF/F is per glomerulus, baselined to that trial's own out-of-stim mean
%      (identical to overshoot_drug_stim_script section 8). For each of the
%      ~12 pulses in a trial, delta_rep(g) = mean dF/F over the 2 s pulse
%      minus mean dF/F over the pre_win_sec window immediately before onset.
%      delta(g) is the mean of delta_rep over pulses. Low-SNR glomeruli
%      (F0 / out-of-stim std < glom_min_snr, same rule as the heatmaps) are
%      NaN, so the mask's midline seam and arch tips drop out instead of
%      blowing up.
%   2. Hemisphere alignment: the 40-glomerulus line is slid against itself.
%      For each candidate shift P in lag_range, the Pearson correlation
%      between delta(k) and delta(k+P) over all valid pairs is computed per
%      trial; the curves are averaged across the fly's trials (weighted by
%      each trial's response amplitude so the strong-signal trials decide
%      and the flat post-drug trials don't add noise) and the shift with the
%      highest mean correlation is the fly's consensus shift P*. The two
%      halves are then paired as (k, k+P*) restricted to k in the first half
%      (k <= nClusters/2) and k+P* in the second half, so each folded
%      position averages exactly one glomerulus from each hemisphere. Anything
%      that can't be paired that way (arch tips, or midline glomeruli when
%      P* > nClusters/2) is left out and shown gray in the pairing figure.
%      Where one member of a pair is NaN the other is used alone (omitnan).
%   3. Folded line per trial, overlaid: 3 side-by-side axes (low / medium /
%      high intensity), colored by drug stage.
%
% Outputs, written into the fly folder next to the other figures:
%   glomDelta.mat                       -- everything computed (see end of section 4)
%   glomDelta_lines_allTrials.png       -- the raw 40-glomerulus delta line per trial
%   glomDelta_alignment.png             -- the sliding/correlation process + an exemplar fold
%   glomDelta_pairing.png               -- which glomeruli got averaged, drawn on the PB
%   glomDelta_folded_by_intensity.png   -- the requested 1x3 overlay
% and, across flies (section 5), written into the repo rather than the share:
%   ugly_figures/exports/glomDelta_folded_allFlies.png -- same 1x3 overlay, one row per fly
%   data/overshoot_glomDelta_allFlies.mat               -- every fly's glom/align structs
%
% Run one %% section at a time, or top-to-bottom (it's batch-safe).

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort')) && isfolder(repoRootDir)
    addpath(repoRootDir); % graph_sort.m (skeleton ordering) lives at the repo root
end

% One entry per fly to process (all six epg_syt8s_lpsp_cschrimson flies of the
% overshoot experiment; each needs mask.mat + imgData_sum_reg.mat from
% overshoot_drug_stim_script first).
overshoot_root = 'Z:\noah_np123\Data\flyg\overshoot\';
fly_dirs = { ...
    fullfile(overshoot_root, '20260921-1_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260921-2_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260921-3_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-1_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-2_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-3_epg_syt8s_lpsp_cschrimson'), ...
    };

% Cross-fly outputs (section 5) go into the repo, not the data share: the
% overshoot root folder is shared with many unrelated experiments.
summary_fig_dir  = fullfile(fileparts(mfilename('fullpath')), '..', 'ugly_figures', 'exports');
summary_data_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'data');
if ~isfolder(summary_fig_dir);  mkdir(summary_fig_dir);  end
if ~isfolder(summary_data_dir); mkdir(summary_data_dir); end

n_per_hemisphere = 20;   % doubled internally -> 40 glomeruli, same as overshoot_drug_stim_script
glom_min_snr     = 3;    % NaN out glomeruli whose out-of-stim F0/std is below this (same rule as the heatmaps)
stim_win_sec     = [0 2];  % window after onset counted as "during the stimulus" (the LED pulse is 2 s)
pre_win_sec      = [-4 0]; % window before onset used as this pulse's own baseline
lag_range        = [12 28]; % candidate hemisphere shifts (in glomeruli) to search; ~nClusters/2 +/- 8
overwrite_glom_cache = false; % per-trial glomerulus traces are cached in the trial folder (glomTraces_<n>.mat)

% Color follows the drug stage NAME, not its row position, so the same drug
% gets the same color in every fly regardless of application order (the 0921
% flies applied picro before mk801, the 0922 flies the reverse).
drug_colors = containers.Map();
drug_colors('baseline')             = [0.165 0.471 0.839]; % blue
drug_colors('ttx')                  = [0.922 0.408 0.204]; % orange
drug_colors('ttx_mec')              = [0.106 0.686 0.478]; % aqua
drug_colors('ttx_mec_mk801')        = [0.929 0.631 0.000]; % yellow
drug_colors('ttx_mec_picro')        = [0.000 0.514 0.000]; % green
drug_colors('ttx_mec_mk801_picro')  = [0.910 0.482 0.643]; % magenta (full cocktail, either order)
drug_colors('ttx_mec_picro_mk801')  = [0.910 0.482 0.643];
unknown_drug_color = [0.5 0.5 0.5];

intensityOrder = {'low', 'medium', 'high'};
nClusters = 2 * n_per_hemisphere;

%% 1..4 run per fly
clear flies % per-fly results collected for section 5; cleared so a rerun with fewer flies can't keep stale rows
for fi = 1:numel(fly_dirs)
    base_dir = fly_dirs{fi};
    pathParts = strsplit(base_dir, filesep);
    pathParts = pathParts(~cellfun(@isempty, pathParts));
    flyName = pathParts{end};
    fprintf('\n==================== %s ====================\n', flyName);

    %% 1. discover stim trials (same logic as overshoot_drug_stim_script section 1)
    [selected, uDrugs] = ogd_discover_trials(base_dir);

    %% 2. mask -> glomeruli, then per-trial per-glomerulus delta
    maskFile = fullfile(base_dir, 'mask.mat');
    if ~isfile(maskFile)
        error('overshoot_glom_delta_script:noMask', 'No mask.mat in %s -- run overshoot_drug_stim_script sections 1-3 first.', base_dir);
    end
    load(maskFile, 'mask');
    maskInfo = dir(maskFile);
    clusterIdx = ogd_pb_glomeruli(mask, n_per_hemisphere);

    glom = struct('name', {}, 'drug', {}, 'intensity', {}, 'trialDir', {}, ...
        'delta', {}, 'delta_reps', {}, 'delta_sem', {}, 'badGlom', {}, 'nReps', {});
    for i = 1:numel(selected)
        tr = ogd_load_glom_traces(selected(i).trialDir, mask, clusterIdx, nClusters, maskInfo.datenum, overwrite_glom_cache);

        % dF/F per glomerulus, F0 = own out-of-stim mean (overshoot_drug_stim_script section 8)
        f0_cluster  = mean(tr.f_cluster(:, ~tr.stims), 2);
        std_cluster = std(tr.f_cluster(:, ~tr.stims), 0, 2);
        badGlom     = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
        dff_cluster = (tr.f_cluster - f0_cluster) ./ f0_cluster;
        dff_cluster(badGlom, :) = NaN;

        onsetIdx   = find(diff(tr.stims) > 0) + 1;
        onsetTimes = tr.t_volume(onsetIdx);
        delta_reps = nan(numel(onsetTimes), nClusters);
        for k = 1:numel(onsetTimes)
            t0 = onsetTimes(k);
            inIdx  = tr.t_volume >= t0 + stim_win_sec(1) & tr.t_volume < t0 + stim_win_sec(2);
            preIdx = tr.t_volume >= t0 + pre_win_sec(1)  & tr.t_volume < t0 + pre_win_sec(2);
            if t0 + pre_win_sec(1) < tr.t_volume(1) || t0 + stim_win_sec(2) > tr.t_volume(end)
                continue % incomplete window at the trial edge -> leave this rep NaN
            end
            delta_reps(k, :) = mean(dff_cluster(:, inIdx), 2)' - mean(dff_cluster(:, preIdx), 2)';
        end
        nReps = sum(all(~isnan(delta_reps(:, ~badGlom)), 2));

        glom(i).name       = selected(i).name;
        glom(i).drug       = selected(i).drug;
        glom(i).intensity  = selected(i).intensity;
        glom(i).trialDir   = selected(i).trialDir;
        glom(i).delta      = mean(delta_reps, 1, 'omitnan');            % 1 x nClusters  <-- the requested single line
        glom(i).delta_reps = delta_reps;                                  % nReps x nClusters
        glom(i).delta_sem  = std(delta_reps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(delta_reps), 1));
        glom(i).badGlom    = badGlom(:)';
        glom(i).nReps      = nReps;
        fprintf('[%s] delta over %d reps; %d/%d glomeruli excluded (low SNR): %s\n', ...
            glom(i).name, nReps, sum(badGlom), nClusters, mat2str(find(badGlom)'));
    end

    % ---- figure: the raw 40-glomerulus line, one per trial, rows = intensity ----
    figure(11); clf
    set(gcf, 'Position', [50 50 1100 800], 'Color', 'w')
    for c = 1:numel(intensityOrder)
        ax = subplot(numel(intensityOrder), 1, c); hold(ax, 'on')
        yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
        xline(ax, n_per_hemisphere + 0.5, '--', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
        for i = find(strcmp({glom.intensity}, intensityOrder{c}))
            plot(ax, 1:nClusters, glom(i).delta, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
                'Color', ogd_drug_color(glom(i).drug, drug_colors, unknown_drug_color), ...
                'MarkerFaceColor', ogd_drug_color(glom(i).drug, drug_colors, unknown_drug_color), ...
                'DisplayName', glom(i).drug)
        end
        xlim(ax, [0.5 nClusters + 0.5]); xticks(ax, 0:5:nClusters)
        title(ax, [intensityOrder{c} ' intensity'])
        ylabel(ax, '\Delta dF/F (stim - pre)')
        if c == numel(intensityOrder); xlabel(ax, 'glomerulus (along PB arch, dashed = midline)'); end
        legend(ax, 'Interpreter', 'none', 'Location', 'eastoutside', 'Box', 'off')
    end
    linkaxes(findobj(gcf, 'Type', 'Axes'), 'xy')
    sgtitle([strrep(flyName, '_', '\_') ': per-glomerulus \DeltadF/F during the LED pulse, one line per trial (mean over pulses)'])
    exportgraphics(gcf, fullfile(base_dir, 'glomDelta_lines_allTrials.png'), 'Resolution', 200);
    fprintf('saved %s\n', fullfile(base_dir, 'glomDelta_lines_allTrials.png'));

    %% 3. slide the two hemispheres across each other -> consensus shift -> pairing
    lags = lag_range(1):lag_range(2);
    corrByLag = nan(numel(glom), numel(lags));
    for i = 1:numel(glom)
        d = glom(i).delta;
        for li = 1:numel(lags)
            P = lags(li);
            a = d(1:nClusters-P); b = d(1+P:nClusters);
            ok = ~isnan(a) & ~isnan(b);
            if sum(ok) >= 8
                r = corrcoef(a(ok), b(ok));
                corrByLag(i, li) = r(1,2);
            end
        end
    end
    [~, bestLagIdx] = max(corrByLag, [], 2, 'omitnan');
    trialBestLag = lags(bestLagIdx);
    trialAmp = arrayfun(@(g) std(g.delta, 'omitnan'), glom); % response amplitude, used as the weight below
    w = trialAmp(:) / sum(trialAmp);
    meanCorr = sum(corrByLag .* w, 1, 'omitnan') ./ sum(w .* ~isnan(corrByLag), 1);
    [~, bi] = max(meanCorr);
    P = lags(bi);

    fprintf('\nhemisphere shift search (lags %d..%d):\n', lag_range);
    for i = 1:numel(glom)
        fprintf('  %-45s best lag %2d (r = %.2f), weight %.2f\n', glom(i).name, trialBestLag(i), corrByLag(i, bestLagIdx(i)), w(i));
    end
    fprintf('  consensus shift P* = %d glomeruli (amplitude-weighted mean r = %.2f)\n', P, meanCorr(bi));

    % pairing: (k, k+P) with k in the first half and k+P in the second half
    kFirst = max(1, n_per_hemisphere + 1 - P) : min(n_per_hemisphere, nClusters - P);
    pairs  = [kFirst(:), kFirst(:) + P]; % nPairs x 2, one glomerulus from each hemisphere
    nPairs = size(pairs, 1);
    unpaired = setdiff(1:nClusters, pairs(:));
    fprintf('  %d pairs: %s\n', nPairs, strjoin(arrayfun(@(a,b) sprintf('%d+%d', a, b), pairs(:,1), pairs(:,2), 'UniformOutput', false)', ', '));
    if ~isempty(unpaired)
        fprintf('  unpaired (left out of the fold): %s\n', mat2str(unpaired));
    end

    for i = 1:numel(glom)
        glom(i).folded = mean([glom(i).delta(pairs(:,1)); glom(i).delta(pairs(:,2))], 1, 'omitnan'); % 1 x nPairs
        % per-rep fold too, so a rep-level SEM of the folded line is available
        glom(i).folded_reps = (glom(i).delta_reps(:, pairs(:,1)) + glom(i).delta_reps(:, pairs(:,2))) / 2;
    end

    % ---- figure: the process ----
    exemplar = find(trialAmp == max(trialAmp), 1); % strongest-response trial, to illustrate the fold
    figure(12); clf
    set(gcf, 'Position', [50 50 1300 800], 'Color', 'w')

    ax1 = subplot(2,2,1); hold on
    yline(0, 'Color', [0.8 0.8 0.8]); xline(n_per_hemisphere + 0.5, '--', 'Color', [0.6 0.6 0.6]);
    plot(1:nClusters, glom(exemplar).delta, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k')
    xlim([0.5 nClusters + 0.5]); xticks(0:5:nClusters)
    xlabel('glomerulus'); ylabel('\Delta dF/F')
    title(['exemplar trial: ' glom(exemplar).name], 'Interpreter', 'none')

    ax2 = subplot(2,2,2); hold on
    for i = 1:numel(glom)
        plot(lags, corrByLag(i,:), '-', 'Color', [ogd_drug_color(glom(i).drug, drug_colors, unknown_drug_color) 0.5], 'LineWidth', 0.75)
    end
    plot(lags, meanCorr, 'k-', 'LineWidth', 2.5)
    xline(P, 'k--', sprintf('P* = %d', P), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
    yline(0, 'Color', [0.8 0.8 0.8]);
    xlabel('shift P (glomeruli)'); ylabel('corr( \Delta(k), \Delta(k+P) )')
    title('slide the line against itself: thin = each trial (color = drug), thick = amplitude-weighted mean')

    ax3 = subplot(2,2,3); hold on
    yline(0, 'Color', [0.8 0.8 0.8]);
    d = glom(exemplar).delta;
    h1 = plot(1:nPairs, d(pairs(:,1)), '-o', 'Color', [0.2 0.4 0.9], 'LineWidth', 1.25, 'MarkerSize', 3, 'MarkerFaceColor', [0.2 0.4 0.9]);
    h2 = plot(1:nPairs, d(pairs(:,2)), '-o', 'Color', [0.9 0.4 0.2], 'LineWidth', 1.25, 'MarkerSize', 3, 'MarkerFaceColor', [0.9 0.4 0.2]);
    h3 = plot(1:nPairs, glom(exemplar).folded, 'k-', 'LineWidth', 2.5);
    xlim([0.5 nPairs + 0.5]); xticks(1:nPairs)
    xticklabels(arrayfun(@(a,b) sprintf('%d\\newline%d', a, b), pairs(:,1), pairs(:,2), 'UniformOutput', false))
    set(gca, 'FontSize', 8)
    xlabel(sprintf('aligned pair (top = glomerulus from 1st half, bottom = its partner, shift %d)', P))
    ylabel('\Delta dF/F')
    legend([h1 h2 h3], {sprintf('1st half, glomeruli %d..%d', pairs(1,1), pairs(end,1)), ...
        sprintf('2nd half shifted by %d, glomeruli %d..%d', P, pairs(1,2), pairs(end,2)), 'mean of the pair (folded)'}, ...
        'Location', 'best', 'Box', 'off')
    title('exemplar: the two halves after alignment, and their average')

    ax4 = subplot(2,2,4); hold on
    for i = 1:numel(glom)
        plot(1:nPairs, glom(i).folded, '-', 'Color', ogd_drug_color(glom(i).drug, drug_colors, unknown_drug_color), 'LineWidth', 1.25)
    end
    yline(0, 'Color', [0.8 0.8 0.8]);
    xlim([0.5 nPairs + 0.5]); xlabel('aligned pair'); ylabel('\Delta dF/F (folded)')
    title('all trials folded with the consensus shift (color = drug)')

    sgtitle([flyName ': hemisphere alignment'], 'Interpreter', 'none')
    exportgraphics(gcf, fullfile(base_dir, 'glomDelta_alignment.png'), 'Resolution', 200);
    fprintf('saved %s\n', fullfile(base_dir, 'glomDelta_alignment.png'));

    % ---- figure: pairing drawn on the PB itself ----
    S = load(fullfile(glom(exemplar).trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    projMean = double(mean(S.imgData_sum_reg, 3)); clear S % video is stored single; mat2gray wants double limits
    pairOfGlom = zeros(1, nClusters);
    for p = 1:nPairs
        pairOfGlom(pairs(p, :)) = p;
    end
    pairImg = zeros(size(mask));
    pairImg(clusterIdx > 0) = pairOfGlom(clusterIdx(clusterIdx > 0));
    cmapPairs = hsv(nPairs);
    rgb = ones([size(mask) 3]) * 0.75; % unpaired glomeruli stay gray
    for ch = 1:3
        tmp = rgb(:,:,ch);
        tmp(pairImg > 0) = cmapPairs(pairImg(pairImg > 0), ch);
        rgb(:,:,ch) = tmp;
    end
    lum = mat2gray(imgaussfilt(projMean, 1), [prctile(projMean(:), 5), prctile(projMean(:), 99.5)]);
    lum = 0.35 + 0.65 * lum; % keep some brightness in dim glomeruli so their color is still readable
    shaded = rgb .* lum;
    shaded(repmat(~mask, [1 1 3])) = 0.08;

    figure(13); clf
    set(gcf, 'Position', [50 50 1200 520], 'Color', 'w')
    image(shaded); axis image; hold on
    contour(mask, [0.5 0.5], 'w', 'LineWidth', 0.75)
    for c = 1:nClusters
        [cy, cx] = find(clusterIdx == c);
        if pairOfGlom(c) > 0
            lbl = sprintf('%d\n(p%d)', c, pairOfGlom(c));
        else
            lbl = sprintf('%d\n(-)', c);
        end
        text(mean(cx), mean(cy), lbl, 'Color', 'w', 'FontSize', 6, 'FontWeight', 'bold', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
    end
    set(gca, 'XTick', [], 'YTick', [])
    title(sprintf('%s: glomeruli sharing a color are averaged together (shift P* = %d; gray = unpaired). Label = glomerulus (pair #)', flyName, P), 'Interpreter', 'none')
    exportgraphics(gcf, fullfile(base_dir, 'glomDelta_pairing.png'), 'Resolution', 200);
    fprintf('saved %s\n', fullfile(base_dir, 'glomDelta_pairing.png'));

    %% 4. folded lines overlaid: 1 x 3 (low / medium / high), color = drug stage
    figure(14); clf
    set(gcf, 'Position', [50 50 1500 480], 'Color', 'w')
    hLeg = gobjects(0); legNames = {};
    for c = 1:numel(intensityOrder)
        ax = subplot(1, numel(intensityOrder), c); hold(ax, 'on')
        yline(ax, 0, 'Color', [0.8 0.8 0.8]);
        for i = find(strcmp({glom.intensity}, intensityOrder{c}))
            col = ogd_drug_color(glom(i).drug, drug_colors, unknown_drug_color);
            h = plot(ax, 1:nPairs, glom(i).folded, '-o', 'Color', col, 'LineWidth', 2, ...
                'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', glom(i).drug);
            if ~ismember(glom(i).drug, legNames)
                hLeg(end+1) = h; legNames{end+1} = glom(i).drug; %#ok<SAGROW>
            end
        end
        xlim(ax, [0.5 nPairs + 0.5]); xticks(ax, 1:nPairs)
        % single-line "k/k+P" labels, rotated: a two-row tex label (as in the alignment
        % figure) loses its line break once MATLAB auto-rotates crowded tick labels
        xticklabels(ax, arrayfun(@(a,b) sprintf('%d/%d', a, b), pairs(:,1), pairs(:,2), 'UniformOutput', false))
        set(ax, 'FontSize', 8, 'XTickLabelRotation', 90)
        title(ax, [intensityOrder{c} ' intensity'], 'FontSize', 11)
        xlabel(ax, 'aligned glomerulus pair (1st-half glomerulus / its 2nd-half partner)')
        if c == 1; ylabel(ax, '\Delta dF/F during pulse (folded across hemispheres)'); end
    end
    linkaxes(findobj(gcf, 'Type', 'Axes'), 'y')
    % legend ordered by drug application order (uDrugs), placed on the last axes
    [~, ord] = ismember(uDrugs, legNames); ord = ord(ord > 0);
    legend(hLeg(ord), legNames(ord), 'Interpreter', 'none', 'Location', 'best', 'Box', 'off')
    sgtitle(sprintf('%s: hemisphere-folded \\DeltadF/F per trial (shift P* = %d), color = drug stage', strrep(flyName, '_', '\_'), P))
    exportgraphics(gcf, fullfile(base_dir, 'glomDelta_folded_by_intensity.png'), 'Resolution', 200);
    fprintf('saved %s\n', fullfile(base_dir, 'glomDelta_folded_by_intensity.png'));

    % ---- save everything ----
    align = struct('lags', lags, 'corrByLag', corrByLag, 'trialBestLag', trialBestLag, ...
        'trialWeight', w', 'meanCorr', meanCorr, 'P', P, 'pairs', pairs, 'unpaired', unpaired);
    params = struct('n_per_hemisphere', n_per_hemisphere, 'glom_min_snr', glom_min_snr, ...
        'stim_win_sec', stim_win_sec, 'pre_win_sec', pre_win_sec, 'lag_range', lag_range);
    save(fullfile(base_dir, 'glomDelta.mat'), 'glom', 'align', 'params', 'uDrugs', 'clusterIdx', 'flyName');
    fprintf('saved %s\n', fullfile(base_dir, 'glomDelta.mat'));

    % keep this fly's results in memory for the cross-fly figure (section 5)
    flies(fi) = struct('name', flyName, 'base_dir', base_dir, 'glom', glom, 'align', align, 'uDrugs', {uDrugs}); %#ok<SAGROW>
end

%% 5. cross-fly grid: rows = flies, cols = intensity, same folded lines as each fly's own figure
% Identical content/styling to section 4's per-fly figure, stacked so every fly
% is a row. Each row keeps its own pairing (its own P*) on the x-axis, and the
% row label carries the fly's P* so a different shift is visible at a glance.
% y-axes are linked across ALL panels so amplitudes compare across flies.
% The legend is the union of drug stages seen in any fly, in application order
% (the 0921 flies applied picro before mk801, the 0922 flies the reverse, so
% both intermediate stages appear).
if numel(flies) > 1
    nFlies = numel(flies);
    drugCanonicalOrder = {'baseline', 'ttx', 'ttx_mec', 'ttx_mec_mk801', 'ttx_mec_picro', 'ttx_mec_mk801_picro', 'ttx_mec_picro_mk801'};

    figure(15); clf
    set(gcf, 'Position', [50 50 1500 270*nFlies + 80], 'Color', 'w')
    tl = tiledlayout(nFlies, numel(intensityOrder), 'TileSpacing', 'compact', 'Padding', 'compact');
    hLeg = gobjects(0); legNames = {};
    for fi = 1:nFlies
        F = flies(fi);
        pairsF = F.align.pairs; nPairsF = size(pairsF, 1);
        shortName = erase(F.name, '_epg_syt8s_lpsp_cschrimson');
        for c = 1:numel(intensityOrder)
            ax = nexttile(tl, (fi-1)*numel(intensityOrder) + c); hold(ax, 'on')
            yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
            for i = find(strcmp({F.glom.intensity}, intensityOrder{c}))
                col = ogd_drug_color(F.glom(i).drug, drug_colors, unknown_drug_color);
                h = plot(ax, 1:nPairsF, F.glom(i).folded, '-o', 'Color', col, 'LineWidth', 2, ...
                    'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', F.glom(i).drug);
                if ~ismember(F.glom(i).drug, legNames)
                    hLeg(end+1) = h; legNames{end+1} = F.glom(i).drug; %#ok<SAGROW>
                end
            end
            xlim(ax, [0.5 nPairsF + 0.5]); xticks(ax, 1:nPairsF)
            xticklabels(ax, arrayfun(@(a,b) sprintf('%d/%d', a, b), pairsF(:,1), pairsF(:,2), 'UniformOutput', false))
            set(ax, 'FontSize', 7, 'XTickLabelRotation', 90)
            if fi == 1; title(ax, [intensityOrder{c} ' intensity'], 'FontSize', 11); end
            if c == 1
                ylabel(ax, {shortName, sprintf('(P* = %d)', F.align.P), '\Delta dF/F (folded)'}, 'FontSize', 9)
            end
            if fi == nFlies; xlabel(ax, 'aligned glomerulus pair (1st-half glomerulus / its 2nd-half partner)', 'FontSize', 8); end
        end
    end
    linkaxes(findobj(gcf, 'Type', 'Axes'), 'y')

    [~, ord] = ismember(drugCanonicalOrder, legNames); ord = ord(ord > 0);
    ord = [ord, setdiff(1:numel(legNames), ord)]; % any drug name not in the canonical list goes last
    lgd = legend(hLeg(ord), legNames(ord), 'Interpreter', 'none', 'Box', 'off');
    lgd.Layout.Tile = 'east';
    title(tl, 'hemisphere-folded \DeltadF/F per trial: rows = flies, cols = LED intensity, color = drug stage')

    outFile = fullfile(summary_fig_dir, 'glomDelta_folded_allFlies.png');
    exportgraphics(gcf, outFile, 'Resolution', 200);
    fprintf('saved %s\n', outFile);

    outMat = fullfile(summary_data_dir, 'overshoot_glomDelta_allFlies.mat');
    save(outMat, 'flies', 'params', 'intensityOrder');
    fprintf('saved %s\n', outMat);
end

%% local functions

function col = ogd_drug_color(drug, drug_colors, fallback)
if isKey(drug_colors, drug)
    col = drug_colors(drug);
else
    col = fallback;
end
end

function tr = ogd_load_glom_traces(trialDir, mask, clusterIdx, nClusters, maskDatenum, overwrite)
% Per-glomerulus mean fluorescence (background-subtracted per frame, same as
% overshoot_drug_stim_script sections 4/8) + stim timing, cached in the trial
% folder because the registered video is ~450 MB per trial on the network
% share and the reduction to 40 traces is the only expensive part.
cacheFile = fullfile(trialDir, sprintf('glomTraces_%d.mat', nClusters));
if isfile(cacheFile) && ~overwrite
    tr = load(cacheFile);
    if isfield(tr, 'maskDatenum') && tr.maskDatenum == maskDatenum
        return
    end
    fprintf('  (cache %s is from an older mask -- recomputing)\n', cacheFile);
end

S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
img = double(S.imgData_sum_reg); clear S
h5list = dir(fullfile(trialDir, '*.h5'));
sync = ogd_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

nVol = size(img, 3);
if numel(sync.stims) ~= nVol
    n = min(numel(sync.stims), nVol);
    sync.stims = sync.stims(1:n); sync.t_volume = sync.t_volume(1:n);
    img = img(:,:,1:n); nVol = n;
end

img_2d = reshape(img, [], nVol);
bg = mean(img_2d(~mask(:), :), 1);
img_2d = img_2d - bg; % per-frame background (mean outside mask) subtraction

centroidLog = false(nClusters, size(img_2d, 1));
for c = 1:nClusters
    centroidLog(c, clusterIdx(:) == c) = true;
end
f_cluster = double(centroidLog) * img_2d ./ sum(centroidLog, 2); % nClusters x nVol

tr = struct('f_cluster', f_cluster, 't_volume', sync.t_volume(:), 'stims', logical(sync.stims(:)), 'maskDatenum', maskDatenum);
save(cacheFile, '-struct', 'tr');
fprintf('  cached %s\n', cacheFile);
end

function [selected, uDrugs] = ogd_discover_trials(base_dir)
% Verbatim logic of overshoot_drug_stim_script section 1 (minus the manual
% override map): every folder with "stim" in the name, parsed into
% (intensity, drug), a bare "_stim" suffix meaning the pre-drug baseline, and
% duplicate (drug,intensity) cells resolved to the candidate with the most LED
% onsets (ties -> later trial number).
d = dir(base_dir);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
stimTrialDirs = d(contains(lower({d.name}), '_stim_') | endsWith(lower({d.name}), '_stim'));

candidates = struct('name', {}, 'trialNum', {}, 'intensity', {}, 'drug', {}, 'nOnsets', {});
for i = 1:numel(stimTrialDirs)
    name = stimTrialDirs(i).name;
    trialDir = fullfile(stimTrialDirs(i).folder, name);
    tok = regexp(lower(name), '(low|medium|high)_stim_(.+)', 'tokens', 'once');
    if isempty(tok)
        tokBase = regexp(lower(name), '(low|medium|high)_stim$', 'tokens', 'once');
        if isempty(tokBase)
            error('overshoot_glom_delta_script:parseCond', 'Could not parse intensity/drug from "%s".', name);
        end
        tok = {tokBase{1}, 'baseline'};
    end
    tnum = regexp(name, 'trial_(\d+)_', 'tokens', 'once');
    h5list = dir(fullfile(trialDir, '*.h5'));
    if numel(h5list) ~= 1 || ~isfile(fullfile(trialDir, 'imgData_sum_reg.mat'))
        fprintf('  note: %s has no h5 or no imgData_sum_reg.mat -- skipping\n', name);
        continue
    end
    sync = ogd_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));
    candidates(end+1) = struct('name', name, 'trialNum', str2double(tnum{1}), ...
        'intensity', tok{1}, 'drug', tok{2}, 'nOnsets', sum(diff(sync.stims) > 0)); %#ok<AGROW>
end

cellKeys = arrayfun(@(c) [c.drug '|' c.intensity], candidates, 'UniformOutput', false);
uCellKeys = unique(cellKeys, 'stable');
selected = struct('name', {}, 'trialDir', {}, 'intensity', {}, 'drug', {});
for k = 1:numel(uCellKeys)
    inCell = find(strcmp(cellKeys, uCellKeys{k}));
    if numel(inCell) > 1
        [~, bestLocal] = max([candidates(inCell).nOnsets] + 1e-6*[candidates(inCell).trialNum]);
        pick = inCell(bestLocal);
        fprintf('  %s: %d candidates -- using %s\n', uCellKeys{k}, numel(inCell), candidates(pick).name);
    else
        pick = inCell;
    end
    selected(end+1) = struct('name', candidates(pick).name, 'trialDir', fullfile(base_dir, candidates(pick).name), ...
        'intensity', candidates(pick).intensity, 'drug', candidates(pick).drug); %#ok<AGROW>
end
uDrugs = unique({selected.drug}, 'stable');
firstTrialNumByDrug = cellfun(@(dd) min([candidates(strcmp({candidates.drug}, dd)).trialNum]), uDrugs);
[~, order] = sort(firstTrialNumByDrug);
uDrugs = uDrugs(order);
fprintf('%d stim trials; drug stages: %s\n', numel(selected), strjoin(uDrugs, ' -> '));
end

function s = ogd_h5_stim_timing(h5Path)
% Same as overshoot_drug_stim_script's ovs_h5_stim_timing.
info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
fs = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
digi = int32(h5read(h5Path, [sweepNames{1} '/digitalScans']));
volCh = find(diNames == "si_volumeclock", 1);
ledCh = find(diNames == "led_stim_feedback", 1);
volBit = bitget(digi, volCh);
ledBit = bitget(digi, ledCh);
N = numel(digi);
t_fine = (0:N-1)' / fs;
t_volume = t_fine(find(diff(volBit) > 0) + 1);
s.fs = fs;
s.t_volume = t_volume;
s.stims = interp1(t_fine, double(ledBit), t_volume, 'nearest', 'extrap') > 0.5;
end

function clusterIdx = ogd_pb_glomeruli(mask, nPerHemisphere)
% Same as overshoot_drug_stim_script's ovs_pb_glomeruli (skeletonize ->
% order with graph_sort -> resample into centroids -> nearest-centroid assign).
[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));
maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('ogd_pb_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch.');
end
D    = bwdistgeodesic(skel, epIdx);
[~, idxFar] = max(D(:));
D2   = bwdistgeodesic(skel, idxFar);
mid  = (D + D2) == mode(D + D2, 'all');
[y_mid, x_mid] = find(mid);
[x_mid, y_mid] = graph_sort(x_mid, y_mid);
xq    = -min_axis:(length(x_mid) + min_axis);
x_mid = round(interp1(1:length(x_mid), x_mid, xq, 'linear', 'extrap'));
y_mid = round(interp1(1:length(y_mid), y_mid, xq, 'linear', 'extrap'));
keep  = ismember([x_mid', y_mid'], [x_mask, y_mask], 'rows');
x_mid = x_mid(keep); y_mid = y_mid(keep);
nClusters = 2 * nPerHemisphere;
xq2 = linspace(1, length(y_mid), 2*nClusters + 1)';
centroids = [interp1(1:length(y_mid), y_mid, xq2), interp1(1:length(x_mid), x_mid, xq2)];
centroids = centroids(2:2:end-1, :);
[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);
clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end
