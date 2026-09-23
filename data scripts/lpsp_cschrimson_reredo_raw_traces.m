%% lpsp_cschrimson_reredo_raw_traces
% Same simplest-possible sanity check as lpsp_cschrimson_redo_raw_traces.m,
% applied to Z:\pablo\lpsp_cschrimson_reredo\ instead: for every trial, plot
% the whole-frame mean pixel intensity (mean over EVERY pixel, no mask) vs.
% the raw stim command, both on the trial's own time axis. One grid figure
% per fly, one subplot per trial (dual y-axis: intensity left, stim command
% right).
%
% Uses registration_001\imagingData.mat ('imgData', pre-summed [Y X
% nVolumes]) and ficTracData_DAQ's 'stim'/'trialTime'/'volClock' columns --
% confirmed this dataset uses the same simple pre-summed format as
% lpsp_cschrimson_redo (not the raw multi-plane 'regProduct' format).
%
% This dataset is small enough (44 ficTracData_DAQ.mat files total) to just
% show every fly, not a sample. Some trials sit directly under a date folder
% with no "fly N" subfolder (e.g. 20250807, 20250820) -- treated as a single
% fly per date in that case, same convention used in
% lpsp_cschrimson_rereredo_script.m.
%
% Run one %% section at a time.

%% 0. parameters
base_dir = 'Z:\pablo\lpsp_cschrimson_reredo\';

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
        flyIds{i} = sprintf('%s_fly1', tok2{1}); % no "fly N" folder -> single fly for that date
    end
end

uFlies = unique(flyIds, 'stable');
fprintf('%d flies: %s\n', numel(uFlies), strjoin(uFlies, ', '));

exportDir = fullfile(fileparts(mfilename('fullpath')), '..', 'ugly_figures', 'exports');
if ~isfolder(exportDir); mkdir(exportDir); end

%% 2. one grid figure per fly, one subplot per trial
for fi = 1:numel(uFlies)
    flyTrialIdx = find(strcmp(flyIds, uFlies{fi}));
    nTrials = numel(flyTrialIdx);
    fprintf('\n=== %s (%d trials) ===\n', uFlies{fi}, nTrials);

    nCols = min(4, nTrials);
    nRows = ceil(nTrials / nCols);
    figure('Color', 'w'); clf
    set(gcf, 'Position', [50 50 380*nCols 230*nRows])

    for ti = 1:nTrials
        trialDir = trialDirs{flyTrialIdx(ti)};
        [~, trialName] = fileparts(trialDir);
        fprintf('  %s\n', trialName);

        if isfile(fullfile(trialDir, 'registration', 'imagingData.mat'))
            imgFile = fullfile(trialDir, 'registration', 'imagingData.mat');
        else
            imgFile = fullfile(trialDir, 'registration_001', 'imagingData.mat');
        end
        S = load(imgFile, 'imgData');
        img = double(S.imgData);
        nVol = size(img, 3);
        frameMean = squeeze(mean(mean(img, 1), 2))';

        ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
        L = load(fullfile(ftFile(1).folder, ftFile(1).name));
        volClock = seconds(L.ftData_DAQ.volClock{1})';
        nb = min(numel(volClock), nVol);
        volClock = volClock(1:nb); frameMean = frameMean(1:nb);

        hasStim = ismember('stim', L.ftData_DAQ.Properties.VariableNames);
        if hasStim
            trialTime = seconds(L.ftData_DAQ.trialTime{1})';
            stim = L.ftData_DAQ.stim{1}(:)';
        end

        subplot(nRows, nCols, ti)
        yyaxis left
        plot(volClock, frameMean, 'b-')
        ylabel('mean frame intensity')
        if hasStim
            yyaxis right
            plot(trialTime, stim, 'r-')
            ylabel('stim command')
        end
        title(trialName, 'Interpreter', 'none', 'FontSize', 7)
        if ti == 1
            xlabel('time in trial (s)')
        end
        axis tight
        set(gca, 'FontSize', 6)
    end
    sgtitle(sprintf('%s: mean frame intensity (blue, left axis) vs. stim command (red, right axis)', uFlies{fi}), 'Interpreter', 'none')

    outFile = fullfile(exportDir, sprintf('fig_lpsp_cschrimson_reredo_rawTraces_%s.png', uFlies{fi}));
    exportgraphics(gcf, outFile, 'Resolution', 150);
    fprintf('saved %s\n', outFile);
end
