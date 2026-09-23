%% fig_lpsp_cschrimson_peristim_heatmaps
% For each of the three lpsp_cschrimson all_data snapshots (repo's .data/
% folder), show one grid figure per (dataset, stim duration) combination:
% one subplot per fly, each a glomerulus x time dF/F heatmap (using that
% fly's own im.d field -- already-computed per-glomerulus dF/F, see
% process_im in lpsp_cschrimson_reredo_script.m / lpsp_cschrimson_debug.m),
% averaged across every repetition of that duration's stim across all of
% that fly's trials, aligned to stim onset, +/- 15 s window.
%
% ---- dataset / genotype notes (see printed output for exact counts) --------
% All three are EPG-bridge two-photon imaging with a CsChrimson stim, split by
% two effector lines: "lpsp" (experimental -- CsChrimson in LPsP) vs "empty"
% (empty split-Gal4 negative control, otherwise identical light stim). Parsed
% from the folder name (contains 'lpsp' or 'empty'), since all_data carries no
% explicit genotype field.
%   lpsp_cschrimson_data_full_20241115.mat -- oldest, EPG>GCaMP6m reporter
%     ("epg_6m_..." folders). 221 trials / 50 flies (regex-parsed, see below).
%     Stim durations found (rounded to nearest 0.5s): 0.5s (32 pulses, sparse
%     -- likely a different/incidental trial type, not the main protocol),
%     3.0s (640 pulses, the dominant protocol), 10s (100), 15s (49), 30s (202).
%   lpsp_cschrimson_redo_20250619.mat -- newer reporter (EPG>syt8s,
%     "epg_syt8s_..." folders). 125 trials / 21 flies. Durations: 0.5s (572
%     pulses, dominant) and 2.0s (8 pulses, sparse).
%   lpsp_cschrimson_reredo_20250806.mat -- same syt8s reporter, smallest/most
%     recent batch, 27 trials / 7 flies, ALL genotype 'lpsp' in this batch (no
%     empty-Gal4 controls collected here). Durations: 2.0s (16), 10.0s (20),
%     30.0s (8).
%
% ---- tricky things worth flagging ------------------------------------------
% 1. Stim duration is NOT fixed within a dataset -- rounding each pulse's
%    measured on/off time to the nearest 0.5s bin (jitter is a few ms,
%    consistent with DAQ sampling, not real protocol variability) reveals
%    multiple distinct nominal durations coexisting in the SAME all_data file.
%    Per the ask, pulses are only ever averaged together with other pulses of
%    the SAME rounded duration -- never mixed across durations.
% 2. +/-15s window vs. stim duration: several durations here (30s, 15s) are
%    LONGER than the requested 15s post-onset window, so those heatmaps only
%    show the first 15s of the pulse, never its offset -- there is no true
%    "post-stim" data in those panels, only "early-in-stim". Flagged in each
%    such figure's sgtitle.
% 3. im.d's dF/F was computed by the ORIGINAL analysis (process_im, F0 = 7th
%    percentile of that trial's own fluorescence, baked into these all_data
%    files) -- this script does NOT recompute it, so it inherits whatever
%    quirks that baseline choice has (e.g. unlike overshoot_drug_stim_script.m's
%    out-of-stim-mean baseline, a low-percentile floor baseline tends to make
%    dF/F skew positive almost everywhere by construction).
% 4. Glomerulus INDEX (1-32, y-axis) is only consistent WITHIN a fly (same
%    mask/skeleton run for all of that fly's own trials) -- it is NOT
%    anatomically aligned ACROSS flies (which end of the skeleton becomes
%    glomerulus 1 depends on that fly's own mask geometry/orientation). Fine
%    for these per-fly panels, but would need proper across-fly bump
%    alignment (e.g. via im.mu) before ever averaging heatmaps ACROSS flies.
% 5. Fly identity is inferred from the folder path ("<8-digit-date>\fly
%    <N>\...") via regex, since all_data has no explicit fly-ID field --
%    printed counts of unmatched trials (if any) let you sanity-check this.
% 6. Windows that fall too close to a trial's start/end (not enough recorded
%    data for the full +/-15s) are NOT dropped -- they contribute NaN for the
%    missing samples (rendered gray in the heatmap, same convention as
%    overshoot_drug_stim_script.m's low-SNR glomerulus exclusion) rather than
%    silently biasing the average with a truncated window.
% 7. Multiple trials for the same fly can have different imaging frame rates
%    (ft.xb spacing) -- every pulse window is resampled (interp1) onto a
%    common +/-15s grid at dt_common (section 0) before averaging, so this is
%    handled, but it does mean each heatmap's temporal resolution is capped
%    at dt_common even if the raw data was finer.
%
% Run one %% section at a time.

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
if isempty(which('bluewhitered'))
    addpath(repoRootDir);
end

dataDir = fullfile(repoRootDir, '.data');
files = {'lpsp_cschrimson_data_full_20241115.mat', 'lpsp_cschrimson_redo_20250619.mat', 'lpsp_cschrimson_reredo_20250806.mat'};

win_pre_sec  = 15;
win_post_sec = 15;
dt_common    = 0.1; % common resampled time grid, seconds -- see note 7 above
tGrid = -win_pre_sec:dt_common:win_post_sec;

exportDir = fullfile(fileparts(mfilename('fullpath')), '..', 'exports');
if ~isfolder(exportDir); mkdir(exportDir); end

%% 1. process each dataset and make one grid figure per stim duration found
for fi = 1:numel(files)
    fpath = fullfile(dataDir, files{fi});
    [~, fbase] = fileparts(files{fi});
    fprintf('\n==================== %s ====================\n', files{fi});
    S = load(fpath, 'all_data');
    ad = S.all_data;

    % -- parse genotype + fly ID + per-trial pulse list --
    trialInfo = struct('flyid', {}, 'genotype', {}, 'onsets_xb', {}, 'dur', {}, 'xb', {}, 'd', {});
    nUnmatched = 0;
    for i = 1:numel(ad)
        if isempty(ad(i).meta) || isempty(ad(i).ft) || ~isfield(ad(i).ft, 'stims') || isempty(ad(i).ft.stims) || isempty(ad(i).im)
            continue
        end
        m = ad(i).meta;
        if contains(m, 'empty', 'IgnoreCase', true)
            genotype = 'empty';
        elseif contains(m, 'lpsp', 'IgnoreCase', true)
            genotype = 'lpsp';
        else
            genotype = 'unknown';
        end
        tok = regexp(m, '(\d{8})\\fly\s*(\d+)', 'tokens', 'once');
        if isempty(tok)
            nUnmatched = nUnmatched + 1;
            flyid = sprintf('unmatched_trial%d', i);
        else
            flyid = sprintf('%s_fly%s', tok{1}, tok{2});
        end

        xf = ad(i).ft.xf(:)';
        stims = double(ad(i).ft.stims(:)') > 0;
        n = min(numel(xf), numel(stims));
        xf = xf(1:n); stims = stims(1:n);

        xb = ad(i).ft.xb(:)';
        d  = ad(i).im.d; % nGlom x nFrames, already dF/F
        nb = min(numel(xb), size(d, 2));
        xb = xb(1:nb); d = d(:, 1:nb);
        if numel(xb) < 2 || numel(xf) < 2; continue; end

        % pulse on/off measured on the fine xf (behavior) clock -- more
        % precise than the coarser imaging clock, matches the earlier survey
        onsets_xf  = find(diff([0, stims]) == 1);
        offsets_xf = find(diff([stims, 0]) == -1);
        m3 = min(numel(onsets_xf), numel(offsets_xf));
        if m3 == 0; continue; end
        durs = xf(offsets_xf(1:m3)) - xf(onsets_xf(1:m3));
        onsetTimes = xf(onsets_xf(1:m3)); % seconds, same clock as xb (both start at trial start)

        k = numel(trialInfo) + 1;
        trialInfo(k).flyid = flyid;
        trialInfo(k).genotype = genotype;
        trialInfo(k).onsets_xb = onsetTimes;
        trialInfo(k).dur = durs;
        trialInfo(k).xb = xb;
        trialInfo(k).d = d;
    end
    fprintf('parsed %d/%d trials with usable stim+im data (%d unmatched fly IDs)\n', numel(trialInfo), numel(ad), nUnmatched);

    % -- bucket all pulses (across all trials) by rounded duration --
    allDurs = cat(2, trialInfo.dur);
    roundedDurs = round(allDurs * 2) / 2; % nearest 0.5s
    durBuckets = unique(roundedDurs);
    fprintf('duration buckets found (s): %s\n', mat2str(durBuckets));

    for bi = 1:numel(durBuckets)
        thisDur = durBuckets(bi);
        fprintf('\n--- %s, duration %.1fs ---\n', fbase, thisDur);

        flyIDs = {};
        flyGeno = {};
        flyHeatmap = {};
        flyNPulses = [];
        for ti = 1:numel(trialInfo)
            T = trialInfo(ti);
            pulseMask = round(T.dur * 2) / 2 == thisDur;
            if ~any(pulseMask); continue; end
            onsets = T.onsets_xb(pulseMask);

            fIdx = find(strcmp(flyIDs, T.flyid), 1);
            if isempty(fIdx)
                flyIDs{end+1} = T.flyid; %#ok<AGROW>
                flyGeno{end+1} = T.genotype; %#ok<AGROW>
                flyHeatmap{end+1} = []; %#ok<AGROW>
                flyNPulses(end+1) = 0; %#ok<AGROW>
                fIdx = numel(flyIDs);
            end

            nGlom = size(T.d, 1);
            for oi = 1:numel(onsets)
                tAbs = T.xb - onsets(oi); % this pulse's onset-relative time
                thisWin = nan(nGlom, numel(tGrid));
                for g = 1:nGlom
                    thisWin(g, :) = interp1(tAbs, T.d(g, :), tGrid, 'linear', NaN);
                end
                flyHeatmap{fIdx} = cat(3, flyHeatmap{fIdx}, thisWin); %#ok<AGROW>
                flyNPulses(fIdx) = flyNPulses(fIdx) + 1;
            end
        end

        if isempty(flyIDs)
            fprintf('  (no flies with this duration -- skipping)\n');
            continue
        end

        nFlies = numel(flyIDs);
        avgHeatmaps = cell(nFlies, 1);
        for k = 1:nFlies
            avgHeatmaps{k} = mean(flyHeatmap{k}, 3, 'omitnan');
            fprintf('  %-24s [%-7s]: %d pulses averaged\n', flyIDs{k}, flyGeno{k}, flyNPulses(k));
        end

        allVals = cellfun(@(x) x(:), avgHeatmaps, 'UniformOutput', false);
        allVals = cat(1, allVals{:});
        climVal = prctile(abs(allVals(~isnan(allVals))), 99);
        if isempty(climVal) || climVal == 0 || isnan(climVal); climVal = 1; end

        nCols = min(8, nFlies);
        nRows = ceil(nFlies / nCols);
        figure('Color', 'w'); clf
        set(gcf, 'Position', [50 50 220*nCols 170*nRows])
        for k = 1:nFlies
            subplot(nRows, nCols, k)
            hm = avgHeatmaps{k};
            imagesc(tGrid, 1:size(hm, 1), hm, 'AlphaData', double(~isnan(hm)))
            set(gca, 'Color', [0.7 0.7 0.7])
            clim([-climVal climVal]); colormap(gca, bluewhitered(256))
            hold on
            xline(0, 'k-', 'LineWidth', 1);
            xline(thisDur, 'k--', 'LineWidth', 1);
            axis tight
            title(sprintf('%s [%s] (n=%d)', flyIDs{k}, flyGeno{k}, flyNPulses(k)), 'Interpreter', 'none', 'FontSize', 7)
            if k == 1
                xlabel('time from stim onset (s)'); ylabel('glomerulus')
            end
            set(gca, 'FontSize', 6)
        end
        durNote = '';
        if thisDur > win_post_sec
            durNote = sprintf(' -- NOTE: %.1fs stim is LONGER than the %.0fs post-onset window, offset not shown', thisDur, win_post_sec);
        end
        sgtitle(sprintf('%s: stim=%.1fs (n=%d flies)%s', fbase, thisDur, nFlies, durNote), 'Interpreter', 'none')

        outFile = fullfile(exportDir, sprintf('fig_%s_dur%.1fs_periStimHeatmap.png', fbase, thisDur));
        exportgraphics(gcf, outFile, 'Resolution', 150);
        fprintf('  saved %s\n', outFile);
    end
end
