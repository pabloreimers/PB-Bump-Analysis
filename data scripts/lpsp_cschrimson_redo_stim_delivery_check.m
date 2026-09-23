%% lpsp_cschrimson_redo_stim_delivery_check
% Same method as lpsp_cschrimson_stim_delivery_check.m (see that script's
% header for the full rationale), applied to the newer syt8s-reporter batch
% at Z:\pablo\lpsp_cschrimson_redo\ instead of the original GCaMP6m batch.
%
% Some trials have a commanded CsChrimson stim (ftData_DAQ.stim{1} goes
% nonzero) but the LED may not have actually fired. This checks every trial
% for evidence that light was ACTUALLY delivered, using bleed-through into
% the raw (unsubtracted) background -- the mean pixel value OUTSIDE that
% trial's own mask.mat.
%
% ---- the method: a LOCAL, per-pulse, matched-control comparison -----------
% For every individual commanded pulse, compare mean background during that
% exact pulse to mean background in the immediately PRECEDING window of the
% SAME duration (not a trial-wide out-of-stim mean, which is vulnerable to
% un-stim-locked periodic artifacts -- see the original script's header for
% why that matters). A trial is called 'LIGHT DETECTED' only if the per-pulse
% local diff is (a) reliably positive, (b) significant by a one-sample
% t-test across that trial's own pulses (alpha=0.01), and (c) not trivially
% small (>1% of that trial's own baseline).
%
% ---- notes specific to this dataset -----------------------------------
% - This dataset mixes fully-preprocessed trials (mask.mat + a registration
%   folder) with raw-only trials (tif only, no mask/registration yet) --
%   e.g. 20250520-2 has no mask.mat at all. Discovery is by searching for
%   mask.mat, so raw-only trials are silently excluded (not an error).
% - registration_001\ holds TWO different registered-video files:
%   'imagingData_reg_ch1_trial*.mat' (raw multi-plane 'regProduct',
%   [Y X nPlanes nVolumes], plus a 'flyback' plane count -- NOT used here,
%   loading/summing this 4-D int16 array per trial was far slower than
%   needed) and the much simpler 'imagingData.mat' (single pre-summed
%   'imgData', [Y X nVolumes], same shape as the other datasets this
%   session). This script uses the latter.
% - Confirmed on a single fly (20250520 fly 1, 4 trials) before running on
%   the whole directory: background values are well-behaved (mean ~4200-
%   5100, std only ~2.7-5.5), and all 4 trials showed no detectable light
%   (local diff -0.01% to +0.04%, p>0.35 throughout) -- small and tight
%   enough across 4 independent trials that "no light" is a fairly credible
%   call there even at only 5 pulses/trial.
% - Earlier survey (this session) found this dataset's stim durations are
%   0.5s (572 pulses, dominant protocol) and 2.0s (8 pulses, sparse), with a
%   median of only ~5 pulses/trial (much lower per-trial power than
%   lpsp_cschrimson's 60-pulse trials) -- unlike lpsp_cschrimson's GCaMP6m
%   batch, this uses the syt8s reporter (much better SNR elsewhere this
%   session), so hopefully more of these short pulses will still resolve.
%
% Run one %% section at a time.

%% 0. parameters
base_dir = 'Z:\pablo\lpsp_cschrimson_redo\';
min_pct_effect = 1;
alpha = 0.01;

%% 1. discover usable trials
maskFiles = dir(fullfile(base_dir, '**', 'mask.mat'));
fprintf('found %d mask.mat files under %s\n', numel(maskFiles), base_dir);

trials = struct('trialDir', {}, 'regPath', {});
for i = 1:numel(maskFiles)
    trialDir = maskFiles(i).folder;
    ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    rp = '';
    if isfile(fullfile(trialDir, 'registration', 'imagingData.mat'))
        rp = fullfile(trialDir, 'registration', 'imagingData.mat');
    elseif isfile(fullfile(trialDir, 'registration_001', 'imagingData.mat'))
        rp = fullfile(trialDir, 'registration_001', 'imagingData.mat');
    end
    if isempty(ftFile) || isempty(rp); continue; end
    k = numel(trials) + 1;
    trials(k).trialDir = trialDir;
    trials(k).regPath = rp;
end
fprintf('%d usable trials (mask + imgData_reg + ficTracData_DAQ all present)\n', numel(trials));

%% 2. per trial: local pulse-vs-preceding-window background comparison
results = struct('trialDir', {}, 'genotype', {}, 'nPulsesChecked', {}, ...
    'meanPctDiff', {}, 'semPctDiff', {}, 'pVal', {}, 'verdict', {});
nNoStimColumn = 0;
nNoStimVariation = 0;
nErrors = 0;

for ti = 1:numel(trials)
    trialDir = trials(ti).trialDir;
    fprintf('[%d/%d] %s\n', ti, numel(trials), trialDir);
    try
        ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
        L = load(fullfile(ftFile(1).folder, ftFile(1).name));
        if ~ismember('stim', L.ftData_DAQ.Properties.VariableNames)
            nNoStimColumn = nNoStimColumn + 1;
            fprintf('  (no ''stim'' column recorded -- skipping)\n');
            continue
        end
        trialTime = seconds(L.ftData_DAQ.trialTime{1});
        volClock  = seconds(L.ftData_DAQ.volClock{1});
        stim      = L.ftData_DAQ.stim{1}(:)' > 0;

        load(fullfile(trialDir, 'mask.mat'), 'mask');
        S = load(trials(ti).regPath, 'imgData');
        img = double(S.imgData);
        nVol = size(img, 3);
        img_2d = reshape(img, [], nVol);
        bg = mean(img_2d(~mask(:), :), 1);

        nb = min(numel(volClock), nVol);
        volClock = volClock(1:nb); volClock = volClock(:)';
        bg = bg(1:nb);

        stims_vol = logical(interp1(trialTime, double(stim), volClock, 'linear', 'extrap') > 0.5);
        stims_vol = stims_vol(:)';
        if ~any(stims_vol) || all(stims_vol)
            nNoStimVariation = nNoStimVariation + 1;
            fprintf('  (no stim variation in this trial -- skipping)\n');
            continue
        end

        baseline = mean(bg(~stims_vol));

        onsetIdx  = find(diff([0, stims_vol]) == 1);
        offsetIdx = find(diff([stims_vol, 0]) == -1);
        m = min(numel(onsetIdx), numel(offsetIdx));
        onsetIdx = onsetIdx(1:m); offsetIdx = offsetIdx(1:m);

        pulseDiffs = [];
        for pi = 1:m
            inIdx = onsetIdx(pi):(offsetIdx(pi) - 1);
            nPulse = numel(inIdx);
            if nPulse < 1; continue; end
            preStart = onsetIdx(pi) - nPulse;
            if preStart < 1; continue; end
            preIdx = preStart:(onsetIdx(pi) - 1);
            pulseDiffs(end+1) = mean(bg(inIdx)) - mean(bg(preIdx)); %#ok<AGROW>
        end

        if numel(pulseDiffs) < 2
            fprintf('  (fewer than 2 usable pulses -- skipping)\n');
            continue
        end

        pctDiffs = 100 * pulseDiffs / baseline;
        meanPct = mean(pctDiffs);
        semPct  = std(pctDiffs) / sqrt(numel(pctDiffs));
        [~, p] = ttest(pctDiffs);

        if meanPct > min_pct_effect && p < alpha
            verdict = 'LIGHT DETECTED';
        elseif meanPct < -min_pct_effect && p < alpha
            verdict = 'NEGATIVE (unexpected -- check trial)';
        else
            verdict = 'no light detected';
        end

        [~, trialFolderName] = fileparts(trialDir);
        genotype = 'unknown';
        if contains(trialFolderName, 'empty', 'IgnoreCase', true) || contains(trialFolderName, '+'); genotype = 'control'; end
        if contains(trialFolderName, 'lpsp', 'IgnoreCase', true) && ~strcmp(genotype, 'control'); genotype = 'lpsp'; end

        k = numel(results) + 1;
        results(k).trialDir = trialDir;
        results(k).genotype = genotype;
        results(k).nPulsesChecked = numel(pctDiffs);
        results(k).meanPctDiff = meanPct;
        results(k).semPctDiff = semPct;
        results(k).pVal = p;
        results(k).verdict = verdict;
        fprintf('  n=%d pulses, mean local diff=%.2f%% +/- %.2f%%, p=%.4g -> %s\n', ...
            numel(pctDiffs), meanPct, semPct, p, verdict);
    catch ME
        nErrors = nErrors + 1;
        fprintf('  ERROR (skipping): %s\n', ME.message);
    end
end

%% 3. summary table
fprintf('\n==================== SUMMARY ====================\n');
fprintf('%d trials checked, %d skipped (no stim column), %d skipped (no stim variation), %d errored\n', ...
    numel(results), nNoStimColumn, nNoStimVariation, nErrors);

verdicts = {results.verdict};
uV = unique(verdicts);
for i = 1:numel(uV)
    fprintf('  %s: %d trials\n', uV{i}, sum(strcmp(verdicts, uV{i})));
end

fprintf('\n--- trials with NO light detected despite a commanded stim ---\n');
for k = 1:numel(results)
    if strcmp(results(k).verdict, 'no light detected')
        fprintf('  %-100s [%s] mean=%.2f%% p=%.3g (n=%d pulses)\n', ...
            results(k).trialDir, results(k).genotype, results(k).meanPctDiff, results(k).pVal, results(k).nPulsesChecked);
    end
end

fprintf('\n--- trials with light detected (for comparison) ---\n');
for k = 1:numel(results)
    if strcmp(results(k).verdict, 'LIGHT DETECTED')
        fprintf('  %-100s [%s] mean=%.2f%% p=%.3g (n=%d pulses)\n', ...
            results(k).trialDir, results(k).genotype, results(k).meanPctDiff, results(k).pVal, results(k).nPulsesChecked);
    end
end

outFile = fullfile(fileparts(mfilename('fullpath')), '..', '.data', 'lpsp_cschrimson_redo_stim_delivery_results.mat');
save(outFile, 'results');
fprintf('\nsaved results struct to %s\n', outFile);
