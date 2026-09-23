%% lpsp_cschrimson_stim_delivery_check
% Some trials in Z:\pablo\lpsp_cschrimson\ have a commanded CsChrimson stim
% (ftData_DAQ.stim{1} goes nonzero) but the LED may not have actually fired
% (e.g. driver fault, disconnected fiber). This checks every trial for
% evidence that light was ACTUALLY delivered, using bleed-through into the
% raw (unsubtracted) background -- the mean pixel value OUTSIDE that trial's
% own mask.mat -- since real light delivery should show up there regardless
% of whether the fly's own EPG/GCaMP signal responded.
%
% ---- why this is trickier than "mean background during stim vs. outside" --
% A naive whole-trial comparison (mean bg during ALL stim frames vs. ALL
% non-stim frames) is NOT reliable here: background fluorescence in this
% dataset has a strong, continuous ~1s-period oscillation that is present
% BEFORE stim onset and carries on through it completely unchanged (checked
% directly on a pooled peri-stim average across 18 trials/604 pulses --
% amplitude ~3% of baseline, no kink/phase-shift at onset). That oscillation
% is NOT stim-triggered (most likely a residual scan/mechanical artifact --
% this dataset's registration pipeline predates this repo's
% claude/pb_remove_scan_noise.m), but a fixed-duration block average can
% still land on different parts of that cycle for "in-stim" vs "out-of-stim"
% purely by timing coincidence, producing a spurious "significant" whole-
% trial difference that has nothing to do with real light delivery.
%
% ---- the fix: a LOCAL, per-pulse, matched-control comparison ---------------
% For every individual commanded pulse, compare mean background during that
% exact pulse to mean background in the immediately PRECEDING window of the
% SAME duration (not the trial-wide out-of-stim mean). Averaging this local
% diff across every pulse in a trial (typically several to a few dozen)
% cancels the un-stim-locked oscillation, because each pulse samples a
% different, effectively random phase of that ~1s cycle relative to its own
% onset -- so across many pulses the oscillation's contribution to the mean
% averages toward zero, while a genuine light-locked bleed-through (which by
% definition IS aligned to every pulse's own onset) survives and accumulates.
% A trial is called out as "light detected" only if the per-pulse local diff
% is (a) reliably positive (light can only ADD to the background, never
% subtract, so a real effect should never be negative) and (b) significant by
% a one-sample t-test across that trial's own pulses (alpha=0.01, stricter
% than usual since ~100+ trials are being screened at once) and (c) not
% trivially small (>1% of that trial's own baseline -- guards against
% statistically-significant-but-meaningless effects from having many pulses).
%
% ---- other notes ----
% - Many of the earliest trials (Sep-Oct 2024) have no 'stim' column at all in
%   ftData_DAQ (added to the acquisition software later) -- these are skipped
%   and reported separately, since there is no commanded-stim record to check
%   against at all.
% - Mask is per-trial (mask.mat inside each trial folder), reused as-is.
% - Registered video variable name varies across this dataset's history --
%   sometimes 'imgData', sometimes 'imgData_reg' (handled below).
%
% Run one %% section at a time.

%% 0. parameters
base_dir = 'Z:\pablo\lpsp_cschrimson\';
min_pct_effect = 1; % minimum |mean diff| as %% of trial baseline to call a trial (see header note c)
alpha = 0.01;

%% 1. discover usable trials
maskFiles = dir(fullfile(base_dir, '**', 'mask.mat'));
fprintf('found %d mask.mat files under %s\n', numel(maskFiles), base_dir);

trials = struct('trialDir', {}, 'regPath', {});
for i = 1:numel(maskFiles)
    trialDir = maskFiles(i).folder;
    ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
    if isfile(fullfile(trialDir, 'registration', 'imgData_reg.mat'))
        rp = fullfile(trialDir, 'registration', 'imgData_reg.mat');
    elseif isfile(fullfile(trialDir, 'registration_001', 'imgData_reg.mat'))
        rp = fullfile(trialDir, 'registration_001', 'imgData_reg.mat');
    else
        rp = '';
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
        S = load(trials(ti).regPath);
        if isfield(S, 'imgData')
            img = double(S.imgData);
        elseif isfield(S, 'imgData_reg')
            img = double(S.imgData_reg);
        else
            error('no recognized image variable (found: %s)', strjoin(fieldnames(S), ', '));
        end
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
            if preStart < 1; continue; end % not enough preceding data -- skip this pulse
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

        % genotype from the TRIAL folder's own name only -- checking the full
        % trialDir path is wrong here, since "lpsp_cschrimson" is also the
        % top-level dataset folder name and would match every trial
        % regardless of that trial's actual (possibly "empty") genotype.
        [~, trialFolderName] = fileparts(trialDir);
        genotype = 'unknown';
        if contains(trialFolderName, 'empty', 'IgnoreCase', true); genotype = 'empty'; end
        if contains(trialFolderName, 'lpsp', 'IgnoreCase', true); genotype = 'lpsp'; end

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

outFile = fullfile(fileparts(mfilename('fullpath')), '..', '.data', 'lpsp_cschrimson_stim_delivery_results.mat');
save(outFile, 'results');
fprintf('\nsaved results struct to %s\n', outFile);
