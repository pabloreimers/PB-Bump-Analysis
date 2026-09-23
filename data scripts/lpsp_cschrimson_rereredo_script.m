%% lpsp_cschrimson_rereredo_script
% Newest LPsP-CsChrimson batch, raw per-trial data at
% Z:\pablo\lpsp_cschrimson_rereredo\<date>\[fly N\]<date>-N_epg_syt8s_lpsp
% _cschrimson\. Unlike the all_data snapshots analyzed earlier this session
% (lpsp_cschrimson_data_full_20241115 / _redo_20250619 / _reredo_20250806,
% which used a global process_im() with a low-percentile-floor dF/F), this
% applies the SAME per-frame background-subtraction + per-glomerulus dF/F
% convention as overshoot_drug_stim_script.m / led_stim_berg4_script.m:
%   1. registered video already exists per trial (registration_001\
%      imgData_reg.mat) -- no motion correction needed here, unlike the raw
%      berg4/overshoot pipeline.
%   2. each frame background-subtracted by that frame's own mean fluorescence
%      OUTSIDE that trial's own mask.mat (per-frame, not dF/F).
%   3. PB divided into 40 glomeruli along the mask skeleton (same
%      skeletonize -> order -> resample method as led_stim_berg4_script.m /
%      overshoot_drug_stim_script.m's *_pb_glomeruli).
%   4. per-glomerulus dF/F, F0 = that trial's own out-of-stim mean per
%      glomerulus; glomeruli whose F0 isn't comfortably above its own
%      out-of-stim noise (glom_min_snr) are NaN'd out (same low-signal-tissue
%      issue as overshoot_drug_stim_script.m section 7 -- mask midline seam /
%      arch endpoints).
%   5. one grid figure: rows/panels = fly (averaged across every stim
%      repetition across ALL of that fly's trials), glomerulus x time,
%      aligned to stim onset.
%
% ---- tricky things worth flagging ------------------------------------------
% 1. Mask is PER-TRIAL here (mask.mat sits inside each trial folder), not one
%    shared per-fly mask like overshoot_drug_stim_script.m built. Per the
%    explicit ask ("you can use the mask and imgData that is in each trial
%    folder"), glomerulus clustering is redone per-trial from that trial's own
%    mask rather than sharing one mask/clustering across a fly's trials. In
%    practice the FOV barely moves within a session so glomerulus numbering
%    should stay close to consistent trial-to-trial for the same fly, but
%    this is NOT guaranteed the way it was when one mask was reused.
% 2. The stim signal (ftData_DAQ.stim{1}) is NOT a clean binary square pulse --
%    it's a smooth analog ramp from 0 up to a plateau value (10.0 here) within
%    each pulse. This script thresholds at >0 (same convention used for the
%    all_data files earlier), so "stim on" spans the ENTIRE ramp-up, not just
%    the plateau. The onset line in each heatmap marks the start of the ramp,
%    not a step to full intensity.
% 3. Stim structure is refreshingly uniform here (unlike the all_data
%    snapshots): every valid trial has exactly 9 pulses of exactly 2.0s each,
%    so there was no need for the duration-bucketing done for
%    fig_lpsp_cschrimson_peristim_heatmaps.m -- only one duration exists.
% 4. Two trials are still mid-upload/processing and were skipped (no
%    registration_001\imgData_reg.mat / mask.mat / ficTracData_DAQ.mat yet):
%    20260420\fly 1\20260420-1 and 20260421\fly 3\20260421-10. Fly 3 on
%    20260421 therefore has ZERO usable trials right now and is excluded
%    entirely -- rerun this script once more trials land.
% 5. 20260415's 5 trials sit directly under the date folder with no "fly N"
%    subfolder (unlike 20260420/20260421) -- treated as a single fly
%    ("20260415_fly1") since nothing in the folder structure suggests
%    otherwise (5 trials in one session is consistent with one fly elsewhere
%    in this dataset).
% 6. Peri-stim window is +/-[4,6]s (pre,post) around each pulse's onset --
%    reused verbatim from overshoot_drug_stim_script.m's own 2s-stim
%    convention (section 9), since the stim here is also exactly 2s.
% 7. Per-trial glomerulus exclusion counts (see printed output) are noisy and
%    sometimes large (a handful of trials excluded 20+/40 glomeruli, others
%    only 1-2) -- directly confirming note 1's worry: per-trial mask quality
%    varies more than the shared per-fly mask in overshoot_drug_stim_script.m
%    did. This mostly washes out at the fly level, since a glomerulus only
%    ends up NaN (gray) in the FINAL averaged heatmap if it failed the SNR
%    check in EVERY one of that fly's trials, not just some -- but it means
%    individual trials' dF/F should be treated with more caution than the
%    fly-level average shown here. Worth revisiting with a shared per-fly
%    mask (like overshoot_drug_stim_script.m) if the per-trial masks turn out
%    to be the culprit.
%
% Run one %% section at a time.

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort'))
    addpath(repoRootDir);
end
claudeDir = fullfile(repoRootDir, 'claude');
if isempty(which('pb_register')) && isfolder(claudeDir)
    addpath(claudeDir);
end

base_dir = 'Z:\pablo\lpsp_cschrimson_rereredo\';
n_per_hemisphere = 20; % -> 40 glomeruli total, matching led_stim_berg4/overshoot convention
glom_min_snr = 3;      % F0 / out-of-stim-noise-SD threshold below which a glomerulus is NaN'd out (see header note)
peristim_pre_sec  = 4;
peristim_post_sec = 6;
dt_common = 0.1; % common resampled time grid, seconds
tGrid = -peristim_pre_sec:dt_common:peristim_post_sec;

%% 1. discover trial folders + group by fly
allTrialDirs = dir(fullfile(base_dir, '**', '*epg_syt8s_lpsp_cschrimson'));
allTrialDirs = allTrialDirs([allTrialDirs.isdir]);
fprintf('found %d candidate trial folders\n', numel(allTrialDirs));

trials = struct('trialDir', {}, 'name', {}, 'flyid', {});
for i = 1:numel(allTrialDirs)
    trialDir = fullfile(allTrialDirs(i).folder, allTrialDirs(i).name);
    hasReg  = isfile(fullfile(trialDir, 'registration_001', 'imgData_reg.mat'));
    hasMask = isfile(fullfile(trialDir, 'mask.mat'));
    ftFile  = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    if ~(hasReg && hasMask && ~isempty(ftFile))
        fprintf('  skip (incomplete -- reg=%d mask=%d ft=%d): %s\n', hasReg, hasMask, ~isempty(ftFile), trialDir);
        continue
    end

    tok = regexp(trialDir, '(\d{8})\\fly\s*(\d+)', 'tokens', 'once');
    if isempty(tok)
        tok2 = regexp(trialDir, '(\d{8})\\', 'tokens', 'once');
        flyid = sprintf('%s_fly1', tok2{1}); % no "fly N" folder -> single fly for that date (see header note 5)
    else
        flyid = sprintf('%s_fly%s', tok{1}, tok{2});
    end

    k = numel(trials) + 1;
    trials(k).trialDir = trialDir;
    trials(k).name = allTrialDirs(i).name;
    trials(k).flyid = flyid;
end
fprintf('%d usable trials across %d flies\n', numel(trials), numel(unique({trials.flyid})));

%% 2. per trial: background-subtract, cluster into glomeruli, compute dF/F, window every pulse
uFlies = unique({trials.flyid}, 'stable');
flyPulseWindows = cell(numel(uFlies), 1); % {fly} -> [nGlom x nTime x nPulses]
flyNPulses = zeros(numel(uFlies), 1);

for ti = 1:numel(trials)
    trialDir = trials(ti).trialDir;
    fIdx = find(strcmp(uFlies, trials(ti).flyid));
    fprintf('[%s] processing (fly %s)...\n', trials(ti).name, trials(ti).flyid);

    load(fullfile(trialDir, 'mask.mat'), 'mask');
    S = load(fullfile(trialDir, 'registration_001', 'imgData_reg.mat'), 'imgData');
    img = double(S.imgData);
    nVol = size(img, 3);

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);
    img_bgsub_2d = reshape(img_bgsub, [], nVol);

    clusterIdx = rrr_pb_glomeruli(mask, n_per_hemisphere);
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
    stim      = L.ftData_DAQ.stim{1}(:)' > 0; % ramp thresholded at >0 -- see header note 2
    nb = min(numel(volClock), nVol);
    volClock = volClock(1:nb);
    volClock = volClock(:)'; % row vector, matches dff_cluster's time dimension below

    stims_vol = logical(interp1(trialTime, double(stim), volClock, 'linear', 'extrap') > 0.5);
    stims_vol = stims_vol(:)'; % row vector -- interp1 shape follows volClock's, which may be a column

    f_cluster = f_cluster(:, 1:nb); % match volClock's length (nb == nVol in practice, guarded defensively)
    f0_cluster  = mean(f_cluster(:, ~stims_vol), 2);
    std_cluster = std(f_cluster(:, ~stims_vol), 0, 2);
    dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;
    badGlom = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
    if any(badGlom)
        fprintf('  excluding %d/%d low-signal glomeruli (F0 too close to noise floor): %s\n', ...
            sum(badGlom), nClusters, mat2str(find(badGlom)'));
    end
    dff_cluster(badGlom, :) = NaN;

    onsetIdx = find(diff([0, stims_vol]) == 1);
    fprintf('  %d stim onsets found\n', numel(onsetIdx));
    for oi = 1:numel(onsetIdx)
        onsetT = volClock(onsetIdx(oi));
        tAbs = volClock - onsetT;
        thisWin = nan(nClusters, numel(tGrid));
        for g = 1:nClusters
            thisWin(g, :) = interp1(tAbs, dff_cluster(g, :), tGrid, 'linear', NaN);
        end
        flyPulseWindows{fIdx} = cat(3, flyPulseWindows{fIdx}, thisWin);
        flyNPulses(fIdx) = flyNPulses(fIdx) + 1;
    end
end

%% 3. average across every pulse for each fly, plot grid figure
rwb = ovs_redwhiteblue(256);

avgHeatmaps = cell(numel(uFlies), 1);
for k = 1:numel(uFlies)
    avgHeatmaps{k} = mean(flyPulseWindows{k}, 3, 'omitnan');
    fprintf('%s: %d pulses averaged\n', uFlies{k}, flyNPulses(k));
end

allVals = cellfun(@(x) x(:), avgHeatmaps, 'UniformOutput', false);
allVals = cat(1, allVals{:});
climVal = prctile(abs(allVals(~isnan(allVals))), 99);

nFlies = numel(uFlies);
nCols = min(4, nFlies);
nRows = ceil(nFlies / nCols);
figure('Color', 'w'); clf
set(gcf, 'Position', [50 50 320*nCols 260*nRows])
for k = 1:nFlies
    subplot(nRows, nCols, k)
    hm = avgHeatmaps{k};
    imagesc(tGrid, 1:size(hm, 1), hm, 'AlphaData', double(~isnan(hm)))
    set(gca, 'Color', [0.7 0.7 0.7])
    clim([-climVal climVal]); colormap(gca, rwb); colorbar
    hold on
    xline(0, 'k-', 'LineWidth', 1); xline(2, 'k--', 'LineWidth', 1); % 2s stim, onset/offset
    axis tight
    title(sprintf('%s (n=%d pulses)', uFlies{k}, flyNPulses(k)), 'Interpreter', 'none')
    xlabel('time from stim onset (s)'); ylabel('glomerulus')
end
sgtitle('lpsp\_cschrimson\_rereredo: per-glomerulus dF/F, mean over all stim repetitions per fly', 'Interpreter', 'tex')

exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports');
if ~isfolder(exportDir); mkdir(exportDir); end
outFile = fullfile(exportDir, 'fig_lpsp_cschrimson_rereredo_periStimHeatmap.png');
exportgraphics(gcf, outFile, 'Resolution', 200);
fprintf('saved %s\n', outFile);

%% local functions

function clusterIdx = rrr_pb_glomeruli(mask, nPerHemisphere)
% Same skeletonize -> order -> resample-into-centroids -> nearest-centroid-
% assign approach as overshoot_drug_stim_script.m's ovs_pb_glomeruli /
% led_stim_berg4_script.m's lsb_pb_glomeruli. Requires graph_sort.m (repo
% root) to order the skeleton into a path.
[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('lpsp_cschrimson_rereredo_script:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
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
