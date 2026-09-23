%% lpsp_cschrimson_reredo_fly_heatmaps
% Per-glomerulus dF/F heatmaps aligned to stim onset, one panel per
% (fly, stim duration) combination, for Z:\pablo\lpsp_cschrimson_reredo\.
% Every commanded pulse is taken at face value as a real stim (no light-
% delivery-confirmation filter this time -- raw traces showed clear,
% visible stim-locked deflections in several flies here, e.g. 20250725 fly 2
% and every trial in 20250807 fly 1, unlike the redo dataset).
%
% Same per-frame background-subtraction + per-glomerulus dF/F convention as
% lpsp_cschrimson_stim_detected_heatmaps.m / overshoot_drug_stim_script.m:
%   1. each frame background-subtracted by that frame's own mean
%      fluorescence OUTSIDE that trial's own mask.mat
%   2. PB divided into 40 glomeruli along the mask skeleton
%   3. per-glomerulus dF/F, F0 = that trial's own out-of-stim mean per
%      glomerulus; low-SNR glomeruli NaN'd out (glom_min_snr=3, the
%      standard threshold used throughout this session -- NOT relaxed the
%      way it had to be for the older GCaMP6m lpsp_cschrimson dataset, since
%      this is the brighter syt8s reporter)
%   4. pulses are pooled PER FLY, but only ever averaged with other pulses
%      of the SAME rounded duration for that fly -- a fly with more than one
%      duration across its trials (e.g. 20250725 fly 2: both ~2s and 10s)
%      gets one panel per duration, not one blended panel.
%
% ---- tricky things worth flagging ------------------------------------------
% 1. Mask is PER-TRIAL (each trial folder has its own mask.mat), not one
%    shared mask per fly -- glomerulus numbering (y-axis) is only guaranteed
%    consistent WITHIN a single trial's own pulses, and only approximately
%    consistent across a fly's different trials/masks. Same caveat as
%    lpsp_cschrimson_rereredo_script.m.
% 2. Peri-stim window is sized PER DURATION: <=2s stim -> [-4,+6]s (matches
%    the 2s convention used elsewhere this session); longer stim -> window
%    scales to [-2*dur, +3*dur]s so the whole pulse stays visible.
% 3. Per the explicit ask this time, every commanded pulse counts as a real
%    stim -- this is a looser assumption than the earlier bleed-through-
%    confirmed analysis for this same dataset (lpsp_cschrimson_stim_detected
%    _heatmaps.m), so trials with a commanded-but-undelivered stim (if any
%    exist here) would still get folded into these averages.
%
% Run one %% section at a time.

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort'))
    addpath(repoRootDir);
end

base_dir = 'Z:\pablo\lpsp_cschrimson_reredo\';
n_per_hemisphere = 20; % -> 40 glomeruli total
glom_min_snr = 3;
apply_snr_exclusion = false; % set true to gray out low-SNR glomeruli again
dt_common = 0.1;

%% 1. discover trials, group by fly
maskFiles = dir(fullfile(base_dir, '**', 'mask.mat'));
trialDirs = {};
for i = 1:numel(maskFiles)
    trialDir = maskFiles(i).folder;
    ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    hasImg = isfile(fullfile(trialDir, 'registration', 'imagingData.mat')) || ...
             isfile(fullfile(trialDir, 'registration_001', 'imagingData.mat'));
    if isempty(ftFile) || ~hasImg; continue; end
    trialDirs{end+1} = trialDir; %#ok<AGROW>
end
fprintf('%d usable trials total\n', numel(trialDirs));

flyIds = cell(numel(trialDirs), 1);
for i = 1:numel(trialDirs)
    tok = regexp(trialDirs{i}, '(\d{8})\\fly\s*(\d+)', 'tokens', 'once');
    if ~isempty(tok)
        flyIds{i} = sprintf('%s_fly%s', tok{1}, tok{2});
    else
        tok2 = regexp(trialDirs{i}, '(\d{8})\\', 'tokens', 'once');
        flyIds{i} = sprintf('%s_fly1', tok2{1});
    end
end
uFlies = unique(flyIds, 'stable');
fprintf('%d flies: %s\n', numel(uFlies), strjoin(uFlies, ', '));

%% 2. per trial: background-subtract, cluster, dF/F, window every pulse by (fly, duration)
groupKeys = {}; % 'flyId|durBucket'
groupPulses = {}; % cell array of [nGlom x nTime] windows
groupTGrid = {};
groupDur = [];
groupFly = {};

for ti = 1:numel(trialDirs)
    trialDir = trialDirs{ti};
    [~, trialName] = fileparts(trialDir);
    fprintf('[%d/%d] %s (%s)\n', ti, numel(trialDirs), trialName, flyIds{ti});

    ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    L = load(fullfile(ftFile(1).folder, ftFile(1).name));
    if ~ismember('stim', L.ftData_DAQ.Properties.VariableNames)
        fprintf('  (no stim column -- skipping)\n');
        continue
    end

    load(fullfile(trialDir, 'mask.mat'), 'mask');
    if isfile(fullfile(trialDir, 'registration', 'imagingData.mat'))
        imgFile = fullfile(trialDir, 'registration', 'imagingData.mat');
    else
        imgFile = fullfile(trialDir, 'registration_001', 'imagingData.mat');
    end
    S = load(imgFile, 'imgData');
    img = double(S.imgData);
    nVol = size(img, 3);

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);
    img_bgsub_2d = reshape(img_bgsub, [], nVol);

    clusterIdx = lcr_pb_glomeruli(mask, n_per_hemisphere);
    nClusters = 2 * n_per_hemisphere;
    centroidLog = false(nClusters, size(img_bgsub_2d, 1));
    for c = 1:nClusters
        centroidLog(c, clusterIdx(:) == c) = true;
    end
    f_cluster = double(centroidLog) * img_bgsub_2d ./ sum(centroidLog, 2);

    trialTime = seconds(L.ftData_DAQ.trialTime{1})';
    volClock  = seconds(L.ftData_DAQ.volClock{1})';
    stim      = L.ftData_DAQ.stim{1}(:)' > 0;
    nb = min(numel(volClock), nVol);
    volClock = volClock(1:nb);
    f_cluster = f_cluster(:, 1:nb);

    stims_vol = logical(interp1(trialTime, double(stim), volClock, 'linear', 'extrap') > 0.5);

    f0_cluster  = mean(f_cluster(:, ~stims_vol), 2);
    std_cluster = std(f_cluster(:, ~stims_vol), 0, 2);
    dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;
    if apply_snr_exclusion
        badGlom = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
        if any(badGlom)
            fprintf('  excluding %d/%d low-SNR glomeruli\n', sum(badGlom), nClusters);
        end
        dff_cluster(badGlom, :) = NaN;
    end

    onsetIdx  = find(diff([0, stims_vol]) == 1);
    offsetIdx = find(diff([stims_vol, 0]) == -1);
    m = min(numel(onsetIdx), numel(offsetIdx));
    fprintf('  %d pulses found\n', m);

    for pi = 1:m
        dur = volClock(offsetIdx(pi)) - volClock(onsetIdx(pi));
        durBucket = round(dur * 2) / 2;
        key = sprintf('%s|%.1f', flyIds{ti}, durBucket);

        gi = find(strcmp(groupKeys, key), 1);
        if isempty(gi)
            if durBucket <= 2
                preSec = 4; postSec = 6;
            else
                preSec = 2 * durBucket; postSec = 3 * durBucket;
            end
            groupKeys{end+1} = key; %#ok<AGROW>
            groupTGrid{end+1} = -preSec:dt_common:postSec; %#ok<AGROW>
            groupPulses{end+1} = []; %#ok<AGROW>
            groupDur(end+1) = durBucket; %#ok<AGROW>
            groupFly{end+1} = flyIds{ti}; %#ok<AGROW>
            gi = numel(groupKeys);
        end

        tGrid = groupTGrid{gi};
        onsetT = volClock(onsetIdx(pi));
        tAbs = volClock - onsetT;
        thisWin = nan(nClusters, numel(tGrid));
        for g = 1:nClusters
            thisWin(g, :) = interp1(tAbs, dff_cluster(g, :), tGrid, 'linear', NaN);
        end
        groupPulses{gi} = cat(3, groupPulses{gi}, thisWin);
    end
end

%% 3. average within each (fly, duration) group, plot one grid figure
rwb = ovs_redwhiteblue(256);

% order groups by fly (in discovery order), then by duration within a fly
[~, flyOrder] = ismember(groupFly, uFlies);
[~, sortIdx] = sortrows([flyOrder(:), groupDur(:)]);

avgHeatmaps = cell(numel(groupKeys), 1);
nPulsesPerGroup = zeros(numel(groupKeys), 1);
for gi = 1:numel(groupKeys)
    avgHeatmaps{gi} = mean(groupPulses{gi}, 3, 'omitnan');
    nPulsesPerGroup(gi) = size(groupPulses{gi}, 3);
    fprintf('%s: %d pulses averaged\n', groupKeys{gi}, nPulsesPerGroup(gi));
end

allVals = cellfun(@(x) x(~isnan(x)), avgHeatmaps, 'UniformOutput', false);
allVals = cat(1, allVals{:});
climVal = prctile(abs(allVals), 99);

nGroups = numel(groupKeys);
nCols = min(5, nGroups);
nRows = ceil(nGroups / nCols);
figure('Color', 'w'); clf
set(gcf, 'Position', [50 50 300*nCols 230*nRows])
for k = 1:nGroups
    gi = sortIdx(k);
    subplot(nRows, nCols, k)
    hm = avgHeatmaps{gi};
    tGrid = groupTGrid{gi};
    imagesc(tGrid, 1:size(hm, 1), hm, 'AlphaData', double(~isnan(hm)))
    set(gca, 'Color', [0.7 0.7 0.7])
    clim([-climVal climVal]); colormap(gca, rwb)
    hold on
    xline(0, 'k-', 'LineWidth', 1); xline(groupDur(gi), 'k--', 'LineWidth', 1);
    axis tight
    title(sprintf('%s, dur=%.1fs (n=%d)', groupFly{gi}, groupDur(gi), nPulsesPerGroup(gi)), ...
        'Interpreter', 'none', 'FontSize', 7)
    if k == 1
        xlabel('time from stim onset (s)'); ylabel('glomerulus')
    end
    set(gca, 'FontSize', 6)
end
snrNote = '';
if ~apply_snr_exclusion
    snrNote = ' -- NO SNR exclusion (every glomerulus shown, including noisy ones)';
end
sgtitle(sprintf('lpsp\\_cschrimson\\_reredo: per-glomerulus dF/F, mean per fly (every command counted as a real stim)%s', snrNote), 'Interpreter', 'tex')

exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports');
if ~isfolder(exportDir); mkdir(exportDir); end
fnameSuffix = '';
if ~apply_snr_exclusion; fnameSuffix = '_noSNRexclusion'; end
outFile = fullfile(exportDir, sprintf('fig_lpsp_cschrimson_reredo_flyHeatmaps%s.png', fnameSuffix));
exportgraphics(gcf, outFile, 'Resolution', 200);
fprintf('saved %s\n', outFile);

%% local functions

function clusterIdx = lcr_pb_glomeruli(mask, nPerHemisphere)
[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('lpsp_cschrimson_reredo_fly_heatmaps:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
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
x_mid = x_mid(keep);
y_mid = y_mid(keep);

nClusters = 2 * nPerHemisphere;
xq2 = linspace(1, length(y_mid), 2*nClusters + 1)';
centroids = [interp1(1:length(y_mid), y_mid, xq2), interp1(1:length(x_mid), x_mid, xq2)];
centroids = centroids(2:2:end-1, :);

[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);

clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end

function cmap = ovs_redwhiteblue(n)
half = n / 2;
blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
whiteToRed  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
cmap = [blueToWhite; whiteToRed];
end
