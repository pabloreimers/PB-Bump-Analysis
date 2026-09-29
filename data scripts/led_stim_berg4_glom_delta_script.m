%% led_stim_berg4_glom_delta_script
% Per-glomerulus stim response as ONE LINE PER TRIAL, folded across the two PB
% hemispheres, for the led_stim_berg4 dataset
% (Z:\pablo\lpsp_cschrimson_epg_syt8s\led_stim_berg4\<fly>\<trial>\): 6 flies x
% 3 LED intensities x {normal saline, high-[K+] saline}, all in TTX. Same
% analysis as overshoot_glom_delta_script.m, so the two datasets' figures read
% the same way -- here the color dimension is saline condition instead of drug
% stage.
%
% The analysis lives in claude/pb_glom_delta_fly.m (per fly) and
% claude/pb_glom_delta_grid.m (cross-fly grid); read their headers for what is
% computed and why. This script only knows this dataset's layout.
%
% Builds on led_stim_berg4_script.m: run its sections 1-2 first so each fly
% has mask.mat and each trial has imgData_sum_reg<label>.mat. A trial folder
% normally holds one acquisition (imgData_sum_reg.mat + wvsrfr_0001.h5); fly
% 4's medium_stim has two complete takes kept as separate trials
% (imgData_sum_reg_0001.mat/_0002.mat with wvsrfr_0001.h5/_0002.h5).
%
% ---- 20 glomeruli, not 40 --------------------------------------------------
% This FOV is only 32x64 px (vs. 200x512 in the overshoot data), so the arch
% is ~80 skeleton pixels long. led_stim_berg4_script.m already divides it into
% 20 (n_per_hemisphere = 10); 40 would be ~2 px per glomerulus, so this keeps
% 20 and pb_glom_delta_fly scales its shift search accordingly (6..14).
%
% Outputs, per fly (in the fly folder): glomDelta.mat + the same four pngs as
% the overshoot script (glomDelta_lines_allTrials / _alignment / _pairing /
% _folded_by_intensity). Across flies, written into the repo:
%   ugly_figures/exports/glomDelta_folded_allFlies_led_stim_berg4.png
%   data/led_stim_berg4_glomDelta_allFlies.mat
%
% Run one %% section at a time, or top-to-bottom (it's batch-safe).

%% 0. parameters
claudeDir = fullfile(fileparts(mfilename('fullpath')), '..', 'claude');
if isempty(which('pb_glom_delta_fly')) && isfolder(claudeDir)
    addpath(claudeDir);
end
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');

base_dir = 'Z:\pablo\lpsp_cschrimson_epg_syt8s\led_stim_berg4\';

normalization    = 'zscore'; % 'dff' = (F-F0)/F0, 'zscore' = (F-F0)/std(F out-of-stim); outputs are prefixed glomDelta_ / glomDeltaZ_ respectively
n_per_hemisphere = 10;    % doubled internally -> 20 glomeruli (see header)
glom_min_snr     = 3;
stim_win_sec     = [0 2];
pre_win_sec      = [-4 0];
overwrite_glom_cache = false; % per-trial traces cached as glomTraces_20<label>.mat in the trial folder

% Every trial here is in TTX; the two conditions differ only in the saline.
cond_colors = containers.Map();
cond_colors('ttx')        = [0.165 0.471 0.839]; % blue   -- normal saline
cond_colors('ttx_high_k') = [0.922 0.408 0.204]; % orange -- high-[K+] saline
cond_order = {'ttx', 'ttx_high_k'};

grid_ylink = 'none'; % cross-fly grid y-axes: 'all' (shared, amplitudes comparable across flies) | 'row' | 'none' (each panel at its own dynamic range)

if strcmpi(normalization, 'zscore'); out_tag = 'glomDeltaZ'; else; out_tag = 'glomDelta'; end
if strcmpi(grid_ylink, 'all'); y_tag = ''; else; y_tag = ['_' grid_ylink 'Y']; end % e.g. _noneY, so the shared-axis png isn't overwritten
summary_fig  = fullfile(repoRootDir, 'ugly_figures', 'exports', [out_tag '_folded_allFlies_led_stim_berg4' y_tag '.png']);
summary_data = fullfile(repoRootDir, 'data', ['led_stim_berg4_' out_tag '_allFlies.mat']);

d = dir(base_dir);
fly_dirs = d([d.isdir] & ~startsWith({d.name}, '.'));
fprintf('found %d fly folders in %s\n', numel(fly_dirs), base_dir);

%% 1. per fly: list acquisitions -> pb_glom_delta_fly
clear flies
for fi = 1:numel(fly_dirs)
    flyDir  = fullfile(fly_dirs(fi).folder, fly_dirs(fi).name);
    flyName = fly_dirs(fi).name;
    fprintf('\n==================== %s ====================\n', flyName);

    maskFile = fullfile(flyDir, 'mask.mat');
    if ~isfile(maskFile)
        error('led_stim_berg4_glom_delta_script:noMask', 'No mask.mat in %s -- run led_stim_berg4_script sections 1-2 first.', flyDir);
    end
    load(maskFile, 'mask');
    maskInfo = dir(maskFile);

    trials = lsbg_list_trials(flyDir, 2*n_per_hemisphere);
    fprintf('%d acquisitions\n', numel(trials));

    res = pb_glom_delta_fly(trials, mask, 'normalization', normalization, ...
        'nPerHemisphere', n_per_hemisphere, 'glomMinSnr', glom_min_snr, ...
        'stimWinSec', stim_win_sec, 'preWinSec', pre_win_sec, ...
        'condColors', cond_colors, 'condOrder', cond_order, 'condLabel', 'saline', ...
        'flyName', flyName, 'outDir', flyDir, ...
        'overwriteCache', overwrite_glom_cache, 'maskDatenum', maskInfo.datenum); %#ok<NASGU>

    save(fullfile(flyDir, [res.tag '.mat']), 'res');
    fprintf('saved %s\n', fullfile(flyDir, [res.tag '.mat']));

    flies(fi) = struct('label', erase(flyName, '_epg_syt8s_lpsp_cschrimson'), 'res', res); %#ok<SAGROW>
end

%% 2. cross-fly grid: rows = flies, cols = intensity, color = saline
if numel(flies) > 1
    pb_glom_delta_grid(flies, 'condColors', cond_colors, 'condOrder', cond_order, ...
        'condLabel', 'saline (all in TTX)', 'yLink', grid_ylink, 'outFile', summary_fig);
    % same grid without the hemisphere fold: the raw line across all glomeruli
    pb_glom_delta_grid(flies, 'condColors', cond_colors, 'condOrder', cond_order, ...
        'condLabel', 'saline (all in TTX)', 'yLink', grid_ylink, 'fold', false, ...
        'outFile', strrep(summary_fig, '_folded_', '_unfolded_'));
    % in-stim minus out-of-stim mean images, rows = (fly, saline), cols = intensity
    pb_stim_diff_grid(flies, 'condOrder', cond_order, ...
        'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', 'stimDiffImages_allFlies_led_stim_berg4.png'));
    % and the same per fly (a copy of each fly folder's stimDiffImages.png, kept here so
    % the whole set can be flipped through in one place)
    for k = 1:numel(flies)
        pb_stim_diff_grid(flies(k), 'condOrder', cond_order, 'figNum', 16, ...
            'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', sprintf('stimDiffImages_led_stim_berg4_%s.png', flies(k).label)));
        % same, but every panel on its own color scale (weak conditions become visible)
        pb_stim_diff_grid(flies(k), 'condOrder', cond_order, 'figNum', 16, 'climMode', 'panel', ...
            'outFile', fullfile(repoRootDir, 'ugly_figures', 'exports', sprintf('stimDiffImages_led_stim_berg4_%s_panelClim.png', flies(k).label)));
    end
    if ~isfolder(fileparts(summary_data)); mkdir(fileparts(summary_data)); end
    save(summary_data, 'flies');
    fprintf('saved %s\n', summary_data);
end

%% local functions

function trials = lsbg_list_trials(flyDir, nClusters)
% One entry per registered acquisition under flyDir/<trial>/. Mirrors
% led_stim_berg4_script's lsb_list_acquisitions + lsb_parse_condition: the
% video is imgData_sum_reg<label>.mat where label is '' for the usual single
% take or '_000N' when a folder holds several takes, and the matching h5 is
% wvsrfr_000N.h5 (wvsrfr_0001.h5 for the unlabeled case). Folder names are
% not perfectly consistent across flies ('low_stim_3', 'med_stim_high_k'), so
% intensity is matched on the leading substring and high-K on contains().
d = dir(flyDir);
trialDirs = d([d.isdir] & ~startsWith({d.name}, '.'));

trials = struct('name', {}, 'videoFile', {}, 'h5Path', {}, 'intensity', {}, 'cond', {}, 'cacheFile', {});
for t = 1:numel(trialDirs)
    trialDir = fullfile(trialDirs(t).folder, trialDirs(t).name);
    n = lower(trialDirs(t).name);
    if startsWith(n, 'low');      intensity = 'low';
    elseif startsWith(n, 'med');  intensity = 'medium';
    elseif startsWith(n, 'high'); intensity = 'high';
    else
        fprintf('  note: skipping %s (not a stim trial folder)\n', trialDirs(t).name);
        continue
    end
    if contains(n, 'high_k'); cond = 'ttx_high_k'; else; cond = 'ttx'; end

    vids = dir(fullfile(trialDir, 'imgData_sum_reg*.mat'));
    if isempty(vids)
        fprintf('  note: %s has no imgData_sum_reg*.mat -- skipping\n', trialDirs(t).name);
        continue
    end
    for v = 1:numel(vids)
        label = erase(erase(vids(v).name, 'imgData_sum_reg'), '.mat'); % '' or '_0001'
        if isempty(label)
            h5Path = fullfile(trialDir, 'wvsrfr_0001.h5');
        else
            h5Path = fullfile(trialDir, ['wvsrfr' label '.h5']);
        end
        if ~isfile(h5Path)
            error('led_stim_berg4_glom_delta_script:h5', 'No %s to go with %s.', h5Path, vids(v).name);
        end
        trials(end+1) = struct('name', [trialDirs(t).name label], ...
            'videoFile', fullfile(trialDir, vids(v).name), 'h5Path', h5Path, ...
            'intensity', intensity, 'cond', cond, ...
            'cacheFile', fullfile(trialDir, sprintf('glomTraces_%d%s.mat', nClusters, label))); %#ok<AGROW>
    end
end
end
