%% overshoot_glom_delta_script
% Quantify the per-glomerulus stim response from the "overshoot" experiment
% (Z:\noah_np123\Data\flyg\overshoot\<fly>\trial_NNN_<intensity>_stim_<drug>\)
% as ONE LINE PER TRIAL, fold the two PB hemispheres onto each other so every
% trial reduces to a single peak + trough, and overlay the folded lines by
% intensity with drug stage as color -- per fly, and as one grid with a row
% per fly.
%
% The analysis itself lives in claude/pb_glom_delta_fly.m (per fly) and
% claude/pb_glom_delta_grid.m (cross-fly grid); read their headers for what is
% computed and why. This script only knows this dataset's layout: it finds the
% stim trials, parses drug/intensity from the folder names, resolves duplicate
% trials, and hands a trial list + mask to pb_glom_delta_fly.
%
% Builds on overshoot_drug_stim_script.m: run that script's sections 1-3 on a
% fly first so mask.mat and each trial's imgData_sum_reg.mat exist.
%
% Outputs, per fly (written into the fly folder next to the heatmaps):
%   glomDelta.mat                       -- res struct from pb_glom_delta_fly
%   glomDelta_lines_allTrials.png       -- raw 40-glomerulus delta line per trial
%   glomDelta_alignment.png             -- the sliding/correlation process + an exemplar fold
%   glomDelta_pairing.png               -- which glomeruli got averaged, drawn on the PB
%   glomDelta_folded_by_intensity.png   -- 1x3 overlay, color = drug stage
% and across flies (written into the repo, not the share -- the overshoot
% root folder is shared with many unrelated experiments):
%   ugly_figures/exports/glomDelta_folded_allFlies_overshoot.png
%   data/overshoot_glomDelta_allFlies.mat
%
% Run one %% section at a time, or top-to-bottom (it's batch-safe).

%% 0. parameters
claudeDir = fullfile(fileparts(mfilename('fullpath')), '..', 'claude');
if isempty(which('pb_glom_delta_fly')) && isfolder(claudeDir)
    addpath(claudeDir);
end
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');

overshoot_root = 'Z:\noah_np123\Data\flyg\overshoot\';
fly_dirs = { ...
    fullfile(overshoot_root, '20260921-1_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260921-2_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260921-3_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-1_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-2_epg_syt8s_lpsp_cschrimson'), ...
    fullfile(overshoot_root, '20260922-3_epg_syt8s_lpsp_cschrimson'), ...
    };

normalization    = 'zscore'; % 'dff' = (F-F0)/F0, 'zscore' = (F-F0)/std(F out-of-stim); outputs are prefixed glomDelta_ / glomDeltaZ_ respectively
n_per_hemisphere = 20;    % doubled internally -> 40 glomeruli, same as overshoot_drug_stim_script
glom_min_snr     = 3;     % NaN out glomeruli whose out-of-stim F0/std is below this (same rule as the heatmaps)
stim_win_sec     = [0 2]; % "during the stimulus" window after onset (the LED pulse is 2 s)
pre_win_sec      = [-4 0];% each pulse's own baseline window before onset
overwrite_glom_cache = false; % per-trial glomerulus traces are cached in the trial folder (glomTraces_40.mat)

% Color follows the drug stage NAME, not its row position, so the same drug
% gets the same color in every fly regardless of application order (the 0921
% flies applied picro before mk801, the 0922 flies the reverse). Both
% four-drug cocktails share one color: same stage, just named in a different
% order on different days.
drug_colors = containers.Map();
drug_colors('baseline')            = [0.165 0.471 0.839]; % blue
drug_colors('ttx')                 = [0.922 0.408 0.204]; % orange
drug_colors('ttx_mec')             = [0.106 0.686 0.478]; % aqua
drug_colors('ttx_mec_mk801')       = [0.929 0.631 0.000]; % yellow
drug_colors('ttx_mec_picro')       = [0.000 0.514 0.000]; % green
drug_colors('ttx_mec_mk801_picro') = [0.910 0.482 0.643]; % magenta
drug_colors('ttx_mec_picro_mk801') = [0.910 0.482 0.643];
drug_order = {'baseline', 'ttx', 'ttx_mec', 'ttx_mec_mk801', 'ttx_mec_picro', 'ttx_mec_mk801_picro', 'ttx_mec_picro_mk801'};

grid_ylink = 'none'; % cross-fly grid y-axes: 'all' (shared, amplitudes comparable across flies) | 'row' | 'none' (each panel at its own dynamic range)

if strcmpi(normalization, 'zscore'); out_tag = 'glomDeltaZ'; else; out_tag = 'glomDelta'; end
if strcmpi(grid_ylink, 'all'); y_tag = ''; else; y_tag = ['_' grid_ylink 'Y']; end % e.g. _noneY, so the shared-axis png isn't overwritten
summary_fig  = fullfile(repoRootDir, 'ugly_figures', 'exports', [out_tag '_folded_allFlies_overshoot' y_tag '.png']);
summary_data = fullfile(repoRootDir, 'data', ['overshoot_' out_tag '_allFlies.mat']);

%% 1. per fly: discover trials -> pb_glom_delta_fly
clear flies
for fi = 1:numel(fly_dirs)
    base_dir = fly_dirs{fi};
    [~, flyName] = fileparts(base_dir);
    fprintf('\n==================== %s ====================\n', flyName);

    selected = ogd_discover_trials(base_dir);

    maskFile = fullfile(base_dir, 'mask.mat');
    if ~isfile(maskFile)
        error('overshoot_glom_delta_script:noMask', 'No mask.mat in %s -- run overshoot_drug_stim_script sections 1-3 first.', base_dir);
    end
    load(maskFile, 'mask');
    maskInfo = dir(maskFile);

    trials = struct('name', {}, 'videoFile', {}, 'h5Path', {}, 'intensity', {}, 'cond', {}, 'cacheFile', {});
    for i = 1:numel(selected)
        h5list = dir(fullfile(selected(i).trialDir, '*.h5'));
        trials(i).name      = selected(i).name;
        trials(i).videoFile = fullfile(selected(i).trialDir, 'imgData_sum_reg.mat');
        trials(i).h5Path    = fullfile(h5list(1).folder, h5list(1).name);
        trials(i).intensity = selected(i).intensity;
        trials(i).cond      = selected(i).drug;
        trials(i).cacheFile = fullfile(selected(i).trialDir, sprintf('glomTraces_%d.mat', 2*n_per_hemisphere));
    end

    res = pb_glom_delta_fly(trials, mask, 'normalization', normalization, ...
        'nPerHemisphere', n_per_hemisphere, 'glomMinSnr', glom_min_snr, ...
        'stimWinSec', stim_win_sec, 'preWinSec', pre_win_sec, ...
        'condColors', drug_colors, 'condOrder', drug_order, 'condLabel', 'drug stage', ...
        'flyName', flyName, 'outDir', base_dir, ...
        'overwriteCache', overwrite_glom_cache, 'maskDatenum', maskInfo.datenum); %#ok<NASGU>

    save(fullfile(base_dir, [res.tag '.mat']), 'res');
    fprintf('saved %s\n', fullfile(base_dir, [res.tag '.mat']));

    flies(fi) = struct('label', erase(flyName, '_epg_syt8s_lpsp_cschrimson'), 'res', res); %#ok<SAGROW>
end

%% 2. cross-fly grid: rows = flies, cols = intensity, color = drug stage
if numel(flies) > 1
    pb_glom_delta_grid(flies, 'condColors', drug_colors, 'condOrder', drug_order, ...
        'condLabel', 'drug stage', 'yLink', grid_ylink, 'outFile', summary_fig);
    % same grid without the hemisphere fold: the raw line across all glomeruli
    pb_glom_delta_grid(flies, 'condColors', drug_colors, 'condOrder', drug_order, ...
        'condLabel', 'drug stage', 'yLink', grid_ylink, 'fold', false, ...
        'outFile', strrep(summary_fig, '_folded_', '_unfolded_'));
    % in-stim minus out-of-stim mean images, rows = (fly, drug stage), cols = intensity
    pb_stim_diff_grid(flies, 'condOrder', drug_order, ...
        'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', 'stimDiffImages_allFlies_overshoot.png'));
    % and the same per fly (a copy of each fly folder's stimDiffImages.png, kept here so
    % the whole set can be flipped through in one place)
    for k = 1:numel(flies)
        pb_stim_diff_grid(flies(k), 'condOrder', drug_order, 'figNum', 16, ...
            'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', sprintf('stimDiffImages_overshoot_%s.png', flies(k).label)));
        % same, but every panel on its own color scale (weak stages become visible)
        pb_stim_diff_grid(flies(k), 'condOrder', drug_order, 'figNum', 16, 'climMode', 'panel', ...
            'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', sprintf('stimDiffImages_overshoot_%s_panelClim.png', flies(k).label)));
    end
    if ~isfolder(fileparts(summary_data)); mkdir(fileparts(summary_data)); end
    save(summary_data, 'flies');
    fprintf('saved %s\n', summary_data);
end

%% local functions

function selected = ogd_discover_trials(base_dir)
% Same logic as overshoot_drug_stim_script section 1 (minus the manual
% override map): every folder with "stim" in the name, parsed into
% (intensity, drug), a bare "_stim" suffix meaning the pre-drug baseline, and
% duplicate (drug,intensity) cells resolved to the candidate with the most LED
% onsets (ties -> later trial number). Trials without a registered video yet
% are skipped with a note.
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
    sync = pb_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));
    candidates(end+1) = struct('name', name, 'trialNum', str2double(tnum{1}), ...
        'intensity', tok{1}, 'drug', tok{2}, 'nOnsets', sum(diff(sync.stims) > 0)); %#ok<AGROW>
end

cellKeys  = arrayfun(@(c) [c.drug '|' c.intensity], candidates, 'UniformOutput', false);
uCellKeys = unique(cellKeys, 'stable');
selected  = struct('name', {}, 'trialDir', {}, 'intensity', {}, 'drug', {});
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
fprintf('%d stim trials; drug stages: %s\n', numel(selected), strjoin(uDrugs(order), ' -> '));
end
