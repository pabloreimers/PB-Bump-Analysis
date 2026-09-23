function res = pb_glom_delta_fly(trials, mask, varargin)
% PB_GLOM_DELTA_FLY  Stim-evoked per-glomerulus dF/F change for one fly, folded across hemispheres.
%
%   res = pb_glom_delta_fly(trials, mask, 'Name', value, ...)
%
% One line per trial: the dF/F change in every glomerulus along the PB arch
% during an LED pulse. Then the two PB hemispheres are slid across each other
% to find the shift that best superimposes them, and matched glomeruli are
% averaged, so each trial reduces to a single peak + trough across one
% hemisphere's worth of glomeruli. Dataset-agnostic: the caller discovers
% trials and names conditions; this does everything from the registered video
% onward. Used by data scripts/overshoot_glom_delta_script.m (conditions =
% drug stages) and data scripts/led_stim_berg4_glom_delta_script.m
% (conditions = saline).
%
% ---- inputs ----------------------------------------------------------------
% trials  struct array, one per trial, fields:
%   .name       label for titles/printouts
%   .videoFile  .mat holding the registered, z-summed video [Y x X x nVol]
%               (variable imgData_sum_reg, or imgData_sum as a fallback)
%   .h5Path     WaveSurfer h5 for this acquisition (see pb_h5_stim_timing)
%   .intensity  'low' | 'medium' | 'high'  (column in the figures)
%   .cond       condition string (color in the figures), e.g. 'ttx', 'high_k'
%   .cacheFile  optional path; per-glomerulus traces are cached here because
%               the video is the only expensive thing to load (leave '' to skip)
% mask    logical [Y x X], one per fly (same mask for all its trials)
%
% ---- name/value options ------------------------------------------------------
%   nPerHemisphere  (20)      doubled internally -> nClusters glomeruli along the arch
%   glomMinSnr      (3)       NaN out a glomerulus whose out-of-stim F0 / std is below this
%   stimWinSec      ([0 2])   window after onset counted as "during the pulse"
%   preWinSec       ([-4 0])  window before onset used as that pulse's own baseline
%   lagRange        (auto)    shifts to search; default nClusters/2 +/- nClusters/5
%                             (40 glomeruli -> 12..28, 20 -> 6..14)
%   minConsensusR   (0.6)     if the best amplitude-weighted mean correlation is below
%                             this, the data can't determine the shift and the anatomical
%                             default P* = nPerHemisphere is used instead (align.Psource
%                             = 'default' rather than 'data'); see step 2 below
%   intensityOrder  ({'low','medium','high'})
%   condColors      (Map)     containers.Map cond -> [r g b]; unknown conds get unknownColor
%   unknownColor    ([.5 .5 .5])
%   condOrder       ({})      legend order; conds not listed follow in first-seen order
%   condLabel       ('condition')  what the colors mean, for titles ('drug stage', 'saline')
%   flyName         ('fly')   for titles / filenames
%   outDir          ('')      where to write the four figures; '' = no figures
%   overwriteCache  (false)
%   maskDatenum     (NaN)     stored in the cache; a cache from a different mask is recomputed
%
% ---- what is computed ------------------------------------------------------
% 1. Per glomerulus g and trial: dF/F = (F - F0)/F0 with F the mean
%    background-subtracted fluorescence in the glomerulus per volume (the
%    frame's own mean OUTSIDE the mask is subtracted first) and F0 that
%    glomerulus's mean over the trial's out-of-stim volumes. Baselining to
%    the out-of-stim mean rather than a low percentile keeps inhibition
%    visible (a percentile floor makes dF/F positive almost everywhere by
%    construction). Glomeruli with F0/std below glomMinSnr are NaN -- these
%    are the mask's midline seam and arch tips, where F0 sits at the noise
%    floor and dF/F blows up meaninglessly.
%    For each pulse, delta_rep(g) = mean dF/F over stimWinSec minus mean over
%    preWinSec; delta(g) = mean over pulses.        -> res.glom(i).delta
% 2. Alignment: for each shift P in lagRange, the Pearson correlation between
%    delta(k) and delta(k+P) over valid pairs, per trial. Curves are averaged
%    across trials weighted by each trial's response amplitude (std of its
%    delta), so strong-response trials decide and flat ones don't add noise;
%    argmax = consensus shift P*. Pairs are (k, k+P*) with k in the first
%    half and k+P* in the second, so every folded position averages exactly
%    one glomerulus from each hemisphere. Arch tips / midline glomeruli that
%    can't be paired that way are left out (res.align.unpaired). Where one
%    member of a pair is NaN the other is used alone.
%                                                    -> res.align, res.glom(i).folded
%
% ---- figures (if outDir is set) ---------------------------------------------
%   <outDir>/glomDelta_lines_allTrials.png      raw nClusters-glomerulus line per trial
%   <outDir>/glomDelta_alignment.png            the sliding/correlation process + an exemplar fold
%   <outDir>/glomDelta_pairing.png              which glomeruli got averaged, drawn on the PB
%   <outDir>/glomDelta_folded_by_intensity.png  1x3 overlay (low/medium/high), color = cond

p = inputParser;
p.addParameter('nPerHemisphere', 20);
p.addParameter('glomMinSnr', 3);
p.addParameter('stimWinSec', [0 2]);
p.addParameter('preWinSec', [-4 0]);
p.addParameter('lagRange', []);
p.addParameter('minConsensusR', 0.6);
p.addParameter('intensityOrder', {'low', 'medium', 'high'});
p.addParameter('condColors', containers.Map());
p.addParameter('unknownColor', [0.5 0.5 0.5]);
p.addParameter('condOrder', {});
p.addParameter('condLabel', 'condition');
p.addParameter('flyName', 'fly');
p.addParameter('outDir', '');
p.addParameter('overwriteCache', false);
p.addParameter('maskDatenum', NaN);
p.parse(varargin{:});
o = p.Results;

nClusters = 2 * o.nPerHemisphere;
if isempty(o.lagRange)
    o.lagRange = round([nClusters/2 - nClusters/5, nClusters/2 + nClusters/5]);
end
minOverlap = max(4, round(0.15 * nClusters)); % fewest valid (k, k+P) pairs needed to trust a correlation
colorOf = @(cond) pb_cond_color(cond, o.condColors, o.unknownColor);
flyTex  = strrep(o.flyName, '_', '\_');

clusterIdx = pb_skeleton_glomeruli(mask, o.nPerHemisphere);

%% 1. per-trial per-glomerulus delta
glom = struct('name', {}, 'cond', {}, 'intensity', {}, 'videoFile', {}, ...
    'delta', {}, 'delta_reps', {}, 'delta_sem', {}, 'badGlom', {}, 'nReps', {});
for i = 1:numel(trials)
    tr = load_glom_traces(trials(i), mask, clusterIdx, nClusters, o.maskDatenum, o.overwriteCache);

    f0_cluster  = mean(tr.f_cluster(:, ~tr.stims), 2);
    std_cluster = std(tr.f_cluster(:, ~tr.stims), 0, 2);
    badGlom     = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < o.glomMinSnr;
    dff_cluster = (tr.f_cluster - f0_cluster) ./ f0_cluster;
    dff_cluster(badGlom, :) = NaN;

    onsetIdx   = find(diff(tr.stims) > 0) + 1;
    onsetTimes = tr.t_volume(onsetIdx);
    delta_reps = nan(numel(onsetTimes), nClusters);
    for k = 1:numel(onsetTimes)
        t0 = onsetTimes(k);
        if t0 + o.preWinSec(1) < tr.t_volume(1) || t0 + o.stimWinSec(2) > tr.t_volume(end)
            continue % incomplete window at the trial edge -> this rep stays NaN
        end
        inIdx  = tr.t_volume >= t0 + o.stimWinSec(1) & tr.t_volume < t0 + o.stimWinSec(2);
        preIdx = tr.t_volume >= t0 + o.preWinSec(1)  & tr.t_volume < t0 + o.preWinSec(2);
        delta_reps(k, :) = mean(dff_cluster(:, inIdx), 2)' - mean(dff_cluster(:, preIdx), 2)';
    end

    glom(i).name       = trials(i).name;
    glom(i).cond       = trials(i).cond;
    glom(i).intensity  = trials(i).intensity;
    glom(i).videoFile  = trials(i).videoFile;
    glom(i).delta      = mean(delta_reps, 1, 'omitnan');   % 1 x nClusters  <-- the single line per trial
    glom(i).delta_reps = delta_reps;                        % nReps x nClusters
    glom(i).delta_sem  = std(delta_reps, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(delta_reps), 1));
    glom(i).badGlom    = badGlom(:)';
    glom(i).nReps      = sum(all(~isnan(delta_reps(:, ~badGlom)), 2));
    fprintf('[%s] delta over %d reps; %d/%d glomeruli excluded (low SNR): %s\n', ...
        glom(i).name, glom(i).nReps, sum(badGlom), nClusters, mat2str(find(badGlom)'));
end

%% 2. slide the hemispheres across each other -> consensus shift -> pairing
lags = o.lagRange(1):o.lagRange(2);
corrByLag = nan(numel(glom), numel(lags));
for i = 1:numel(glom)
    d = glom(i).delta;
    for li = 1:numel(lags)
        P = lags(li);
        a = d(1:nClusters-P); b = d(1+P:nClusters);
        ok = ~isnan(a) & ~isnan(b);
        if sum(ok) >= minOverlap
            r = corrcoef(a(ok), b(ok));
            corrByLag(i, li) = r(1,2);
        end
    end
end
[~, bestLagIdx] = max(corrByLag, [], 2, 'omitnan');
trialBestLag = lags(bestLagIdx);
trialBestLag(all(isnan(corrByLag), 2)) = NaN; % a trial with no valid pairs at any lag (e.g. all glomeruli NaN) has no best lag
trialAmp = arrayfun(@(g) std(g.delta, 'omitnan'), glom);
trialAmp(isnan(trialAmp)) = 0; % an all-NaN trial carries no weight rather than poisoning the mean
if sum(trialAmp) > 0
    w = trialAmp(:) / sum(trialAmp);
else
    w = ones(numel(glom), 1) / numel(glom);
end
meanCorr = sum(corrByLag .* w, 1, 'omitnan') ./ sum(w .* ~isnan(corrByLag), 1);
[bestR, bi] = max(meanCorr);
if isnan(bestR) || bestR < o.minConsensusR
    % The data can't pin the shift down (no trial has a clear periodic
    % response), so fall back to the anatomical default: the skeleton
    % clustering splits the arch evenly, so glomerulus k's partner is k +
    % nPerHemisphere. Keeping a noise-driven shift here would fold unrelated
    % glomeruli together and make the flat flies look structured.
    P = o.nPerHemisphere;
    Psource = 'default';
else
    P = lags(bi);
    Psource = 'data';
end

fprintf('\nhemisphere shift search (lags %d..%d):\n', o.lagRange);
for i = 1:numel(glom)
    if isnan(trialBestLag(i))
        fprintf('  %-45s no valid pairs (all glomeruli NaN)\n', glom(i).name);
    else
        fprintf('  %-45s best lag %2d (r = %.2f), weight %.2f\n', glom(i).name, trialBestLag(i), corrByLag(i, bestLagIdx(i)), w(i));
    end
end
if strcmp(Psource, 'data')
    fprintf('  consensus shift P* = %d glomeruli (amplitude-weighted mean r = %.2f)\n', P, bestR);
else
    fprintf('  consensus shift P* = %d glomeruli -- DEFAULT (nPerHemisphere): best amplitude-weighted mean r = %.2f is below minConsensusR = %.2f\n', ...
        P, bestR, o.minConsensusR);
end

kFirst = max(1, o.nPerHemisphere + 1 - P) : min(o.nPerHemisphere, nClusters - P);
pairs  = [kFirst(:), kFirst(:) + P];
nPairs = size(pairs, 1);
unpaired = setdiff(1:nClusters, pairs(:));
fprintf('  %d pairs: %s\n', nPairs, strjoin(arrayfun(@(a,b) sprintf('%d+%d', a, b), pairs(:,1), pairs(:,2), 'UniformOutput', false)', ', '));
if ~isempty(unpaired)
    fprintf('  unpaired (left out of the fold): %s\n', mat2str(unpaired));
end

for i = 1:numel(glom)
    glom(i).folded      = mean([glom(i).delta(pairs(:,1)); glom(i).delta(pairs(:,2))], 1, 'omitnan'); % 1 x nPairs
    glom(i).folded_reps = (glom(i).delta_reps(:, pairs(:,1)) + glom(i).delta_reps(:, pairs(:,2))) / 2;
end

align = struct('lags', lags, 'corrByLag', corrByLag, 'trialBestLag', trialBestLag, ...
    'trialWeight', w', 'meanCorr', meanCorr, 'bestR', bestR, 'P', P, 'Psource', Psource, ...
    'pairs', pairs, 'unpaired', unpaired);

conds = unique({glom.cond}, 'stable');
[~, ord] = ismember(o.condOrder, conds); ord = ord(ord > 0);
conds = conds([ord, setdiff(1:numel(conds), ord, 'stable')]);

res = struct('flyName', o.flyName, 'glom', glom, 'align', align, 'clusterIdx', clusterIdx, ...
    'nClusters', nClusters, 'conds', {conds}, 'params', o);

if isempty(o.outDir)
    return
end
if ~isfolder(o.outDir); mkdir(o.outDir); end

%% figure: raw per-glomerulus line, one per trial, rows = intensity
figure(11); clf
set(gcf, 'Position', [50 50 1100 800], 'Color', 'w')
for c = 1:numel(o.intensityOrder)
    ax = subplot(numel(o.intensityOrder), 1, c); hold(ax, 'on')
    yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    xline(ax, o.nPerHemisphere + 0.5, '--', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    for i = find(strcmp({glom.intensity}, o.intensityOrder{c}))
        plot(ax, 1:nClusters, glom(i).delta, '-o', 'LineWidth', 1.5, 'MarkerSize', 3, ...
            'Color', colorOf(glom(i).cond), 'MarkerFaceColor', colorOf(glom(i).cond), 'DisplayName', glom(i).cond)
    end
    xlim(ax, [0.5 nClusters + 0.5]); xticks(ax, 0:5:nClusters)
    title(ax, [o.intensityOrder{c} ' intensity'])
    ylabel(ax, '\Delta dF/F (stim - pre)')
    if c == numel(o.intensityOrder); xlabel(ax, 'glomerulus (along PB arch, dashed = midline)'); end
    legend(ax, 'Interpreter', 'none', 'Location', 'eastoutside', 'Box', 'off')
end
linkaxes(findobj(gcf, 'Type', 'Axes'), 'xy')
sgtitle([flyTex ': per-glomerulus \DeltadF/F during the LED pulse, one line per trial (mean over pulses)'])
exportgraphics(gcf, fullfile(o.outDir, 'glomDelta_lines_allTrials.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(o.outDir, 'glomDelta_lines_allTrials.png'));

%% figure: the alignment process
exemplar = find(trialAmp == max(trialAmp), 1);
figure(12); clf
set(gcf, 'Position', [50 50 1300 800], 'Color', 'w')

subplot(2,2,1); hold on
yline(0, 'Color', [0.8 0.8 0.8]); xline(o.nPerHemisphere + 0.5, '--', 'Color', [0.6 0.6 0.6]);
plot(1:nClusters, glom(exemplar).delta, 'k-o', 'LineWidth', 1.5, 'MarkerSize', 3, 'MarkerFaceColor', 'k')
xlim([0.5 nClusters + 0.5]); xticks(0:5:nClusters)
xlabel('glomerulus'); ylabel('\Delta dF/F')
title(['exemplar trial: ' glom(exemplar).name], 'Interpreter', 'none')

subplot(2,2,2); hold on
for i = 1:numel(glom)
    plot(lags, corrByLag(i,:), '-', 'Color', [colorOf(glom(i).cond) 0.5], 'LineWidth', 0.75)
end
plot(lags, meanCorr, 'k-', 'LineWidth', 2.5)
xline(P, 'k--', sprintf('P* = %d', P), 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom');
yline(0, 'Color', [0.8 0.8 0.8]);
xlabel('shift P (glomeruli)'); ylabel('corr( \Delta(k), \Delta(k+P) )')
title(sprintf('slide the line against itself: thin = each trial (color = %s), thick = amplitude-weighted mean', o.condLabel))

subplot(2,2,3); hold on
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

subplot(2,2,4); hold on
for i = 1:numel(glom)
    plot(1:nPairs, glom(i).folded, '-', 'Color', colorOf(glom(i).cond), 'LineWidth', 1.25)
end
yline(0, 'Color', [0.8 0.8 0.8]);
xlim([0.5 nPairs + 0.5]); xlabel('aligned pair'); ylabel('\Delta dF/F (folded)')
title(sprintf('all trials folded with the consensus shift (color = %s)', o.condLabel))

sgtitle([o.flyName ': hemisphere alignment'], 'Interpreter', 'none')
exportgraphics(gcf, fullfile(o.outDir, 'glomDelta_alignment.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(o.outDir, 'glomDelta_alignment.png'));

%% figure: pairing drawn on the PB itself
projMean = double(mean(load_video(glom(exemplar).videoFile), 3));
pairOfGlom = zeros(1, nClusters);
for pp = 1:nPairs
    pairOfGlom(pairs(pp, :)) = pp;
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
lum = 0.35 + 0.65 * lum; % keep dim glomeruli readable
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
title(sprintf('%s: glomeruli sharing a color are averaged together (shift P* = %d; gray = unpaired). Label = glomerulus (pair #)', o.flyName, P), 'Interpreter', 'none')
exportgraphics(gcf, fullfile(o.outDir, 'glomDelta_pairing.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(o.outDir, 'glomDelta_pairing.png'));

%% figure: folded lines overlaid, 1 x nIntensities, color = cond
figure(14); clf
set(gcf, 'Position', [50 50 1500 480], 'Color', 'w')
hLeg = gobjects(0); legNames = {};
for c = 1:numel(o.intensityOrder)
    ax = subplot(1, numel(o.intensityOrder), c); hold(ax, 'on')
    yline(ax, 0, 'Color', [0.8 0.8 0.8], 'HandleVisibility', 'off');
    for i = find(strcmp({glom.intensity}, o.intensityOrder{c}))
        col = colorOf(glom(i).cond);
        h = plot(ax, 1:nPairs, glom(i).folded, '-o', 'Color', col, 'LineWidth', 2, ...
            'MarkerSize', 4, 'MarkerFaceColor', col, 'DisplayName', glom(i).cond);
        if ~ismember(glom(i).cond, legNames)
            hLeg(end+1) = h; legNames{end+1} = glom(i).cond; %#ok<AGROW>
        end
    end
    xlim(ax, [0.5 nPairs + 0.5]); xticks(ax, 1:nPairs)
    % single-line "k/k+P" labels, rotated: a two-row tex label (as in the alignment
    % figure) loses its line break once MATLAB auto-rotates crowded tick labels
    xticklabels(ax, arrayfun(@(a,b) sprintf('%d/%d', a, b), pairs(:,1), pairs(:,2), 'UniformOutput', false))
    set(ax, 'FontSize', 8, 'XTickLabelRotation', 90)
    title(ax, [o.intensityOrder{c} ' intensity'], 'FontSize', 11)
    xlabel(ax, 'aligned glomerulus pair (1st-half glomerulus / its 2nd-half partner)')
    if c == 1; ylabel(ax, '\Delta dF/F during pulse (folded across hemispheres)'); end
end
linkaxes(findobj(gcf, 'Type', 'Axes'), 'y')
[~, ord] = ismember(conds, legNames); ord = ord(ord > 0);
legend(hLeg(ord), legNames(ord), 'Interpreter', 'none', 'Location', 'best', 'Box', 'off')
sgtitle(sprintf('%s: hemisphere-folded \\DeltadF/F per trial (shift P* = %d), color = %s', flyTex, P, o.condLabel))
exportgraphics(gcf, fullfile(o.outDir, 'glomDelta_folded_by_intensity.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(o.outDir, 'glomDelta_folded_by_intensity.png'));
end

%% ---- local functions ---------------------------------------------------------

function col = pb_cond_color(cond, condColors, fallback)
if isKey(condColors, cond)
    col = condColors(cond);
else
    col = fallback;
end
end

function img = load_video(videoFile)
vars = who('-file', videoFile);
if ismember('imgData_sum_reg', vars)
    S = load(videoFile, 'imgData_sum_reg'); img = S.imgData_sum_reg;
elseif ismember('imgData_sum', vars)
    S = load(videoFile, 'imgData_sum'); img = S.imgData_sum;
else
    error('pb_glom_delta_fly:video', 'No imgData_sum_reg / imgData_sum variable in %s.', videoFile);
end
end

function tr = load_glom_traces(trial, mask, clusterIdx, nClusters, maskDatenum, overwrite)
% Per-glomerulus mean fluorescence (per-frame background-subtracted) + stim
% timing, cached because the registered video is the only expensive load.
useCache = isfield(trial, 'cacheFile') && ~isempty(trial.cacheFile);
if useCache && isfile(trial.cacheFile) && ~overwrite
    tr = load(trial.cacheFile);
    if isfield(tr, 'maskDatenum') && isequaln(tr.maskDatenum, maskDatenum)
        return
    end
    fprintf('  (cache %s is from a different mask -- recomputing)\n', trial.cacheFile);
end

img  = double(load_video(trial.videoFile));
sync = pb_h5_stim_timing(trial.h5Path);

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
if useCache
    save(trial.cacheFile, '-struct', 'tr');
    fprintf('  cached %s\n', trial.cacheFile);
end
end
