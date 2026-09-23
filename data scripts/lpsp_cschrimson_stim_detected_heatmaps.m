%% lpsp_cschrimson_stim_detected_heatmaps
% For every trial that lpsp_cschrimson_stim_delivery_check.m flagged
% 'LIGHT DETECTED' (real bleed-through confirming the commanded CsChrimson
% stim actually fired), show a stim-onset-aligned per-glomerulus dF/F
% heatmap, one panel per trial. Same per-frame background-subtraction + per-
% glomerulus dF/F convention as lpsp_cschrimson_rereredo_script.m /
% overshoot_drug_stim_script.m:
%   1. each frame background-subtracted by that frame's own mean
%      fluorescence OUTSIDE that trial's own mask.mat
%   2. PB divided into 40 glomeruli along the mask skeleton
%   3. per-glomerulus dF/F, F0 = that trial's own out-of-stim mean per
%      glomerulus; low-SNR glomeruli (F0 not comfortably above its own out-
%      of-stim noise) NaN'd out and rendered gray, same as elsewhere this
%      session
%   4. one grid figure PER GENOTYPE (control vs experimental) -- see note
%      below on the genotype split -- one panel per TRIAL (not averaged
%      across trials/flies, since these are the individual trials that were
%      specifically flagged as confirmed-light).
%
% ---- genotype split -------------------------------------------------------
% "control" = trial folder name contains 'empty' OR '+' (matches this lab's
% own convention elsewhere, e.g. lpsp_cschrimson_reredo_script.m's
% empty_idx = contains(meta,'empty') | contains(meta,'+')) -- i.e. an empty-
% split-Gal4 negative control, marked either way depending on how that
% session named the fly. "experimental" = contains 'lpsp' and isn't already
% caught by the control check above.
%
% ---- window sizing ---------------------------------------------------------
% Peri-stim window is derived PER TRIAL from that trial's own measured pulse
% duration (median across its pulses), not a fixed number -- this dataset's
% stim durations vary trial-to-trial (unlike the rereredo dataset, which was
% uniformly 2s). For duration <= 2s, window is [-4, +6]s (matches the
% convention used elsewhere this session for a 2s stim); for longer
% durations, window scales to [-2*dur, +3*dur]s so the whole pulse (plus
% context before/after) stays visible.
%
% Run one %% section at a time.

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort'))
    addpath(repoRootDir);
end

resultsFile = fullfile(repoRootDir, '.data', 'lpsp_cschrimson_stim_delivery_results.mat');
load(resultsFile, 'results');

n_per_hemisphere = 20; % -> 40 glomeruli total
% glom_min_snr=3 (used everywhere else this session, tuned on the much
% brighter syt8s/overshoot data) wiped out nearly every glomerulus here --
% checked directly on 20250205-4 (180 pulses, the best-sampled trial in this
% set): background mean=29.6 vs. whole-image mean=33.2 (barely any contrast),
% and even the BRIGHTEST glomeruli only reach SNR~2.4-2.6. This GCaMP6m
% cohort is genuinely dimmer than the newer syt8s data, not a mask bug, so
% the bar is relaxed to 1.5 here -- results should be read with more caution
% than the syt8s-based figures elsewhere this session.
glom_min_snr = 1.5;
apply_snr_exclusion = false; % set true to gray out low-SNR glomeruli again (see header note)
dt_common = 0.1;

detected = results(strcmp({results.verdict}, 'LIGHT DETECTED'));
fprintf('%d trials flagged LIGHT DETECTED -- building heatmaps for all of them\n', numel(detected));

%% 1. per-trial: background-subtract, cluster, dF/F, window every pulse
trialHeatmaps = cell(numel(detected), 1);
trialDur = zeros(numel(detected), 1);
trialGeno = cell(numel(detected), 1);
trialLabel = cell(numel(detected), 1);

for ti = 1:numel(detected)
    trialDir = detected(ti).trialDir;
    [~, trialFolderName] = fileparts(trialDir);
    fprintf('[%d/%d] %s\n', ti, numel(detected), trialFolderName);

    if contains(trialFolderName, 'empty', 'IgnoreCase', true) || contains(trialFolderName, '+')
        geno = 'control';
    else
        geno = 'experimental';
    end
    trialGeno{ti} = geno;
    trialLabel{ti} = trialFolderName;

    load(fullfile(trialDir, 'mask.mat'), 'mask');
    if isfile(fullfile(trialDir, 'registration', 'imgData_reg.mat'))
        regPath = fullfile(trialDir, 'registration', 'imgData_reg.mat');
    else
        regPath = fullfile(trialDir, 'registration_001', 'imgData_reg.mat');
    end
    S = load(regPath);
    if isfield(S, 'imgData')
        img = double(S.imgData);
    else
        img = double(S.imgData_reg);
    end
    nVol = size(img, 3);

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);
    img_bgsub_2d = reshape(img_bgsub, [], nVol);

    clusterIdx = lcs_pb_glomeruli(mask, n_per_hemisphere);
    nClusters = 2 * n_per_hemisphere;
    centroidLog = false(nClusters, size(img_bgsub_2d, 1));
    for c = 1:nClusters
        centroidLog(c, clusterIdx(:) == c) = true;
    end
    f_cluster = double(centroidLog) * img_bgsub_2d ./ sum(centroidLog, 2);

    ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    L = load(fullfile(ftFile(1).folder, ftFile(1).name));
    trialTime = seconds(L.ftData_DAQ.trialTime{1});
    volClock  = seconds(L.ftData_DAQ.volClock{1});
    stim      = L.ftData_DAQ.stim{1}(:)' > 0;
    nb = min(numel(volClock), nVol);
    volClock = volClock(1:nb); volClock = volClock(:)';
    f_cluster = f_cluster(:, 1:nb);

    stims_vol = logical(interp1(trialTime, double(stim), volClock, 'linear', 'extrap') > 0.5);
    stims_vol = stims_vol(:)';

    f0_cluster  = mean(f_cluster(:, ~stims_vol), 2);
    std_cluster = std(f_cluster(:, ~stims_vol), 0, 2);
    dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;
    if apply_snr_exclusion
        badGlom = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
        dff_cluster(badGlom, :) = NaN;
    end

    onsetIdx  = find(diff([0, stims_vol]) == 1);
    offsetIdx = find(diff([stims_vol, 0]) == -1);
    m = min(numel(onsetIdx), numel(offsetIdx));
    onsetIdx = onsetIdx(1:m); offsetIdx = offsetIdx(1:m);
    dur = median(volClock(offsetIdx) - volClock(onsetIdx));
    trialDur(ti) = dur;

    if dur <= 2
        preSec = 4; postSec = 6;
    else
        preSec = 2 * dur; postSec = 3 * dur;
    end
    tGrid = -preSec:dt_common:postSec;

    thisTrialWins = nan(nClusters, numel(tGrid), m);
    for oi = 1:m
        onsetT = volClock(onsetIdx(oi));
        tAbs = volClock - onsetT;
        for g = 1:nClusters
            thisTrialWins(g, :, oi) = interp1(tAbs, dff_cluster(g, :), tGrid, 'linear', NaN);
        end
    end
    trialHeatmaps{ti} = struct('avg', mean(thisTrialWins, 3, 'omitnan'), 'tGrid', tGrid, 'nPulses', m);
end

%% 2. one grid figure per genotype
rwb = ovs_redwhiteblue(256);
exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports');
if ~isfolder(exportDir); mkdir(exportDir); end

genoGroups = {'control', 'experimental'};
for gi = 1:numel(genoGroups)
    idx = find(strcmp(trialGeno, genoGroups{gi}));
    if isempty(idx); continue; end

    allVals = [];
    for k = idx'
        v = trialHeatmaps{k}.avg;
        allVals = [allVals; v(~isnan(v))]; %#ok<AGROW>
    end
    climVal = prctile(abs(allVals), 99);

    nTrials = numel(idx);
    nCols = min(5, nTrials);
    nRows = ceil(nTrials / nCols);
    figure('Color', 'w'); clf
    set(gcf, 'Position', [50 50 280*nCols 230*nRows])
    for kk = 1:nTrials
        k = idx(kk);
        H = trialHeatmaps{k};
        subplot(nRows, nCols, kk)
        imagesc(H.tGrid, 1:size(H.avg, 1), H.avg, 'AlphaData', double(~isnan(H.avg)))
        set(gca, 'Color', [0.7 0.7 0.7])
        clim([-climVal climVal]); colormap(gca, rwb)
        hold on
        xline(0, 'k-', 'LineWidth', 1); xline(trialDur(k), 'k--', 'LineWidth', 1);
        axis tight
        title(sprintf('%s (dur=%.1fs, n=%d)', trialLabel{k}, trialDur(k), H.nPulses), ...
            'Interpreter', 'none', 'FontSize', 7)
        if kk == 1
            xlabel('time from stim onset (s)'); ylabel('glomerulus')
        end
        set(gca, 'FontSize', 6)
    end
    snrNote = '';
    if ~apply_snr_exclusion
        snrNote = ' -- NO SNR exclusion (every glomerulus shown, including noisy ones)';
    end
    sgtitle(sprintf('lpsp\\_cschrimson, %s flies (light confirmed delivered), n=%d trials%s', genoGroups{gi}, nTrials, snrNote), 'Interpreter', 'tex')

    fnameSuffix = '';
    if ~apply_snr_exclusion; fnameSuffix = '_noSNRexclusion'; end
    outFile = fullfile(exportDir, sprintf('fig_lpsp_cschrimson_stimDetected_%s_heatmaps%s.png', genoGroups{gi}, fnameSuffix));
    exportgraphics(gcf, outFile, 'Resolution', 200);
    fprintf('saved %s\n', outFile);
end

%% local functions

function clusterIdx = lcs_pb_glomeruli(mask, nPerHemisphere)
[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('lpsp_cschrimson_stim_detected_heatmaps:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
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
