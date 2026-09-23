%% lpsp_cschrimson_scan_long_stims
% Scan all four lpsp_cschrimson* raw-data folders under Z:\pablo\ for every
% trial whose commanded stim is longer than 1 second (any single pulse,
% not just the dominant/modal duration for that trial). Pure metadata scan
% -- reads only each trial's *ficTracData_DAQ.mat (small), never touches
% imaging data, so this is fast even across hundreds of trials.
%
% Duration is measured the same way as everywhere else this session:
% threshold ftData_DAQ.stim{1} at >0 (this also catches the analog-ramp-
% shaped stim seen in the rereredo dataset, which counts the whole ramp-up
% as "on"), find onset/offset run boundaries, measure each pulse's duration
% on the fine trialTime clock.
%
% Also flags whether each qualifying trial is imaging-ready (mask.mat +
% either imgData_reg/imagingData.mat present), since that's what actually
% matters if the next step is building a heatmap for it -- but nothing is
% excluded from the list on that basis, since the ask here is just "find the
% trials", not "find the analyzable trials".
%
% Run one %% section at a time.

%% 0. parameters
base_dirs = {
    'Z:\pablo\lpsp_cschrimson\'
    'Z:\pablo\lpsp_cschrimson_redo\'
    'Z:\pablo\lpsp_cschrimson_reredo\'
    'Z:\pablo\lpsp_cschrimson_rereredo\'
    };
min_dur_sec = 1.0;

%% 1. scan
results = struct('baseDir', {}, 'trialDir', {}, 'genotype', {}, 'durations', {}, 'maxDur', {}, 'imagingReady', {});
nNoStimColumn = 0;
nErrors = 0;
nScanned = 0;

for bi = 1:numel(base_dirs)
    base_dir = base_dirs{bi};
    ftFiles = dir(fullfile(base_dir, '**', '*ficTracData_DAQ.mat'));
    fprintf('\n==================== %s ====================\n', base_dir);
    fprintf('found %d ficTracData_DAQ.mat files\n', numel(ftFiles));

    for i = 1:numel(ftFiles)
        trialDir = ftFiles(i).folder;
        nScanned = nScanned + 1;
        try
            L = load(fullfile(ftFiles(i).folder, ftFiles(i).name));
            if ~ismember('stim', L.ftData_DAQ.Properties.VariableNames)
                nNoStimColumn = nNoStimColumn + 1;
                continue
            end
            trialTime = seconds(L.ftData_DAQ.trialTime{1})';
            stim = L.ftData_DAQ.stim{1}(:)' > 0;
            n = min(numel(trialTime), numel(stim));
            trialTime = trialTime(1:n); stim = stim(1:n);

            onsets  = find(diff([0, stim]) == 1);
            offsets = find(diff([stim, 0]) == -1);
            m = min(numel(onsets), numel(offsets));
            if m == 0; continue; end
            durs = trialTime(offsets(1:m)) - trialTime(onsets(1:m));

            if max(durs) > min_dur_sec
                [~, trialFolderName] = fileparts(trialDir);
                genotype = 'unknown';
                if contains(trialFolderName, 'empty', 'IgnoreCase', true) || contains(trialFolderName, '+')
                    genotype = 'control';
                elseif contains(trialFolderName, 'lpsp', 'IgnoreCase', true)
                    genotype = 'lpsp';
                end

                hasMask = isfile(fullfile(trialDir, 'mask.mat'));
                hasImg = isfile(fullfile(trialDir, 'registration', 'imgData_reg.mat')) || ...
                         isfile(fullfile(trialDir, 'registration_001', 'imgData_reg.mat')) || ...
                         isfile(fullfile(trialDir, 'registration', 'imagingData.mat')) || ...
                         isfile(fullfile(trialDir, 'registration_001', 'imagingData.mat'));

                k = numel(results) + 1;
                results(k).baseDir = base_dir;
                results(k).trialDir = trialDir;
                results(k).genotype = genotype;
                results(k).durations = durs;
                results(k).maxDur = max(durs);
                results(k).imagingReady = hasMask && hasImg;
            end
        catch ME
            nErrors = nErrors + 1;
            fprintf('  ERROR on %s: %s\n', trialDir, ME.message);
        end
    end
end

%% 2. report
fprintf('\n==================== SUMMARY ====================\n');
fprintf('%d trials scanned total, %d skipped (no stim column), %d errored\n', nScanned, nNoStimColumn, nErrors);
fprintf('%d trials found with a pulse > %.1fs\n', numel(results), min_dur_sec);

for bi = 1:numel(base_dirs)
    idx = strcmp({results.baseDir}, base_dirs{bi});
    fprintf('\n--- %s (%d qualifying trials) ---\n', base_dirs{bi}, sum(idx));
    thisResults = results(idx);
    for k = 1:numel(thisResults)
        readyStr = 'not imaging-ready';
        if thisResults(k).imagingReady; readyStr = 'imaging-ready'; end
        fprintf('  %-100s [%-7s] max=%.2fs all=%s (%s)\n', ...
            thisResults(k).trialDir, thisResults(k).genotype, thisResults(k).maxDur, ...
            mat2str(round(thisResults(k).durations, 2)), readyStr);
    end
end

outFile = fullfile(fileparts(mfilename('fullpath')), '..', '.data', 'lpsp_cschrimson_long_stim_trials.mat');
save(outFile, 'results');
fprintf('\nsaved results struct to %s\n', outFile);
