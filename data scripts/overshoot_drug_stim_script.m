%% overshoot_drug_stim_script
% One fly per run from the "overshoot" experiment, in
% Z:\noah_np123\Data\flyg\overshoot\<flyname>\<trial_NNN_...>\. LED (CsChrimson)
% stim at 3 intensities (low/medium/high), repeated after each drug in a
% sequentially-applied cocktail (TTX -> +mecamylamine -> +picrotoxin ->
% +MK801, i.e. each stage includes every drug applied so far; some flies also
% have a pre-drug "baseline" stage -- see below). Only trial folders with
% "stim" in the name matter -- other trial_NNN folders are a different trial
% type entirely (e.g. closed-loop/VR walking trials, no tif at all) and are
% ignored here.
%
% This is the same raw-data format as led_stim_berg4_script.m (raw ScanImage
% fastZ tif + one WaveSurfer h5, both handled the same way), so this script
% reuses that one's approach where it applies -- but several things about this
% dataset are different enough that this is its own script rather than a
% change to that one.
%
% ---- FLY-SPECIFIC PARAMETERS: change these for each new fly -----------------
% base_dir and ovs_fallback_scan_params (section 0) are the two things that
% actually change between flies; everything else below is shared logic.
% scan_params.json has been missing for large stretches of every fly checked
% so far, so ovs_fallback_scan_params can't just be read from one -- see the
% per-fly notes below for how each fly's config was actually determined.
%
% 20260921-1: scan_params.json exists for trial_009-016 (nplanes=13, 3
%   flyback, 10 valid, 100x256 px) and is trusted there; missing from
%   trial_017 on (ttx_mec_picro / ttx_mec_picro_mk801 stages) -- those trials'
%   tif dims/file size match trial_012_high_stim_ttx almost exactly, so this
%   script assumes they share its scan config.
% 20260921-2: scan_params.json is missing for EVERY trial (checked directly --
%   none of them have it). si_frameclock/si_volumeclock edge counts give
%   framesPerVolume=8 exactly (nplanes=8), and reading a raw sample of volumes
%   showed planes 1-6 have real signal (mean ~30-75, std ~400-670) while
%   planes 7-8 are blank (mean ~0, std ~95) -- see fly2_planes.png from
%   workshopping this -- so nplanes_valid=6, flyback_planes=[6,7] (0-indexed).
%   tif is 200x512 px (exactly 2x each dimension vs. 20260921-1's 100x256);
%   um_per_px is a GUESS (halved from 20260921-1's 0.253125, assuming the same
%   physical FOV at 2x the raster resolution, not a 2x-smaller FOV at the same
%   resolution) -- only used to scale the mask smoothing/morphology
%   parameters, and checked against mask_qc.png before trusting it.
% This fly also has a bare "..._high_stim" trial (trial_005, no drug suffix at
% all) -- a pre-drug BASELINE condition, not a parsing mistake. Section 1
% treats a trial ending in "_stim" with nothing after it as drug='baseline'.
% 20260921-3: same missing-json situation for most trials; checked directly
% (200x512 px, planesPerVolume=8, planes 1-6 real / 7-8 blank -- fly3_planes.png)
% and it's identical to 20260921-2's config. One trial (trial_008_low_stim_ttx
% _mec) unexpectedly DID have scan_params.json -- confirmed nplanes/flyback/px
% dims exactly, and gave a real um_per_px (0.140625) to replace the earlier
% 0.126563 guess.
% 20260922-1, 20260922-2: scan_params.json present for every stim trial (new
%   day, new cfg: pb_10plane_10p5hz.cfg, zoom=5.0) -- nplanes=8, valid=6,
%   flyback=[6,7], 200x512 px, um_per_px=0.1265625 (real, not guessed). Note
%   this drug cocktail's order changed from the 921 flies: now
%   ttx -> ttx_mec -> ttx_mec_mk801 -> ttx_mec_mk801_picro (mk801 before
%   picro). Not hardcoded anywhere -- section 1 derives row order from each
%   drug string's first trial number, so this just works either way.
%   20260922-1 also has a duplicate trial_004/005_low_stim_ttx_mec (picked
%   005, 12 onsets, over 004's 3 onsets -- an aborted/restarted trial).
% 20260922-3: scan_params.json missing for every trial; checked directly
%   (same method as 20260921-2/3) and confirmed identical config to
%   20260922-1/2 above (200x512, 8 planes/volume, planes 1-6 real / 7-8
%   blank), so reuses their fallback. Also has a bare "..._high_stim" baseline
%   trial (trial_004), same as 20260921-2's pattern.
%
% ---- duplicate-labeled trials: picks the more complete one, not both --------
% If a (drug,intensity) cell has more than one candidate trial folder, this
% script picks whichever has the most LED stim onsets (ties broken by the
% later trial number) and prints what it picked/excluded. See led_stim_berg4
% _script.m's fly-4 case for the alternative (both equally complete -> keep
% both) -- override this heuristic per-cell via ovs_trial_overrides if needed.
%
% ---- what's actually computed here: -----------------------------------------
%   1. raw tif -> motion-corrected, z-summed video (pb_register, same as
%      led_stim_berg4_script.m's section 1/1b)
%   2. one mask for this fly; each frame background-subtracted by that frame's
%      own mean fluorescence OUTSIDE the mask (a per-frame subtraction, not a
%      dF/F normalization -- sections 4-6 use this directly as a raw F diff)
%   3. mean(in-stim frames) - mean(out-of-stim frames), per trial
%   4. one grid figure (drugStim_diffGrid.png): rows = drug stage
%      (chronological), cols = intensity, one diff image per trial averaged
%      over all its repetitions
%   5. one per-repetition breakdown figure (perRepDiffGrid.png): rows =
%      (drug,intensity) condition, cols = repetition number, each cell = that
%      one 2 s pulse vs. its own preceding 8 s (not the trial-wide average)
%   6. (sections 7-9, NOT a per-pixel dF/F -- see section 7's comment for why)
%      the PB divided into 40 glomeruli along the mask skeleton (same method
%      as led_stim_berg4_script.m's section 7); a per-glomerulus dF/F
%      (background-subtracted F, baselined to that trial's own out-of-stim
%      mean per glomerulus, NaN'd out below glom_min_snr -- see section 0);
%      one full-trial heatmap per trial (glomHeatmap_fullTrial_*.png) and one
%      stim-averaged peri-onset heatmap per trial (glomHeatmap_periStim_*.png)
%
% Run one %% section at a time.

%% 0. parameters
claudeDir = fullfile(fileparts(mfilename('fullpath')), '..', 'claude');
if isempty(which('pb_register')) && isfolder(claudeDir)
    addpath(claudeDir);
end
% graph_sort.m (used by section 7's glomerulus clustering) lives at the repo root.
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort')) && isfolder(repoRootDir)
    addpath(repoRootDir);
end

% ---- change these two for a new fly (see header) ----
% 20260922-3: missing scan_params.json for every trial; verified directly
% (same method as 20260921-2/3) that it shares the identical config as its
% same-day siblings 20260922-1/2 -- 200x512 px, 8 planes/volume, planes 1-6
% real signal / 7-8 blank flyback -- so it reuses their fallback below.
base_dir = 'Z:\noah_np123\Data\flyg\overshoot\20260922-3_epg_syt8s_lpsp_cschrimson\';
ovs_fallback_scan_params = struct( ...
    'nplanes', 8, 'nplanes_valid', 6, 'flyback_planes', [6,7], ...
    'nchannels', 1, 'channels_saved', 1, 'px_height', 200, 'px_width', 512, ...
    'um_per_px', [0.1265625, 0.1265625, 5]);
% um_per_px above is MEASURED (from 20260922-1/2's real scan_params.json,
% zoom=5.0), not guessed -- this is a different rig config/day than
% 20260921-3's 0.140625 (zoom=4.5), so don't reuse that value here.
% ---- everything below is shared across flies ----

do_motion_correction = true;  % rigid NoRMCorre registration (claude/pb_register.m)
overwrite_video       = false; % re-convert tif -> imgData_sum.mat even if it already exists
overwrite_reg         = false; % redo motion correction even if imgData_sum_reg.mat already exists
overwrite_mask        = false; % rebuild mask.mat even if it already exists

n_per_hemisphere = 20; % sections 7-9's glomerulus clustering, doubled internally -> 40 total, as requested
% Some glomeruli end up sitting on near-zero real signal -- the mask's
% midline "bridge" (imclose seam joining the two PB hemispheres, section 3)
% and the arch's very endpoints both average out to F0 close to the noise
% floor. Dividing by that near-zero, noisy F0 makes dF/F blow up there with
% no real meaning, so sections 8-9 NaN out any glomerulus whose out-of-stim
% F0 isn't comfortably above its own out-of-stim noise (F0 / std(F, out-of-
% stim) < glom_min_snr), rather than plotting the resulting garbage.
glom_min_snr = 3;

% mask smoothing/morphology tuned in PHYSICAL units (um), carried over from
% workshopping led_stim_berg4_script.m's mask on that dataset's 1.0125 um/px
% scale, then converted to pixels below using THIS dataset's own um_per_px --
% so the same effective smoothing applies despite the differing pixel scale.
mask_smooth_sigma_um = 2.5 * 1.0125; % ~2.53 um
mask_open_radius_um  = 1   * 1.0125; % ~1.01 um
mask_join_line_um    = 15  * 1.0125; % ~15.2 um

% Manual override for the "most complete duplicate wins" heuristic (section 1)
% -- map from '<drug>|<intensity>' to the trial folder name to force, if the
% automatic pick (most LED onsets, ties to later trial number) picks wrong.
ovs_trial_overrides = containers.Map('KeyType', 'char', 'ValueType', 'char');

pathParts = strsplit(base_dir, filesep);
pathParts = pathParts(~cellfun(@isempty, pathParts));
flyName = pathParts{end}; % e.g. '20260921-2_epg_syt8s_lpsp_cschrimson', for figure titles

%% 1. discover stim trials, parse drug/intensity, pick the most complete of any duplicates
allTrialDirs = lsb_list_subdirs_ovs(base_dir);
stimTrialDirs = allTrialDirs(contains(lower({allTrialDirs.name}), '_stim_') | endsWith(lower({allTrialDirs.name}), '_stim'));
fprintf('found %d trial folders, %d look like stim trials\n', numel(allTrialDirs), numel(stimTrialDirs));

candidates = struct('name', {}, 'trialNum', {}, 'intensity', {}, 'drug', {}, 'nOnsets', {}, 'durationSec', {});
for i = 1:numel(stimTrialDirs)
    name = stimTrialDirs(i).name;
    trialDir = fullfile(stimTrialDirs(i).folder, name);

    tok = regexp(lower(name), '(low|medium|high)_stim_(.+)', 'tokens', 'once');
    if isempty(tok)
        % no drug suffix at all (e.g. "..._high_stim") -- a pre-drug baseline trial
        tokBase = regexp(lower(name), '(low|medium|high)_stim$', 'tokens', 'once');
        if isempty(tokBase)
            error('overshoot_drug_stim_script:parseCond', 'Could not parse intensity/drug from "%s".', name);
        end
        tok = {tokBase{1}, 'baseline'};
    end
    tnum = regexp(name, 'trial_(\d+)_', 'tokens', 'once');

    % A trial folder can exist with its tif fully written but h5/metadata
    % still landing (this data arrives incrementally) -- skip it for now
    % rather than crash the whole discovery pass; it'll show up next run.
    h5list = dir(fullfile(trialDir, '*.h5'));
    tiflist = dir(fullfile(trialDir, '*.tif'));
    if numel(h5list) ~= 1 || numel(tiflist) ~= 1
        fprintf('  note: %s looks incomplete (found %d tif, %d h5) -- skipping for now, still uploading?\n', ...
            name, numel(tiflist), numel(h5list));
        continue
    end
    sync = ovs_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

    candidates(end+1) = struct('name', name, 'trialNum', str2double(tnum{1}), ... %#ok<SAGROW>
        'intensity', tok{1}, 'drug', tok{2}, 'nOnsets', sum(diff(sync.stims)>0), ...
        'durationSec', sync.durationSec);
end

fprintf('\n%-45s %8s %10s %10s\n', 'trial', 'nOnsets', 'duration(s)', 'trialNum');
for i = 1:numel(candidates)
    fprintf('%-45s %8d %10.1f %10d\n', candidates(i).name, candidates(i).nOnsets, candidates(i).durationSec, candidates(i).trialNum);
end

% one row per unique (drug,intensity), keeping the most-complete candidate
cellKeys = arrayfun(@(c) [c.drug '|' c.intensity], candidates, 'UniformOutput', false);
uCellKeys = unique(cellKeys, 'stable');
selected = struct('name', {}, 'trialDir', {}, 'intensity', {}, 'drug', {});
fprintf('\n--- duplicate resolution ---\n');
for k = 1:numel(uCellKeys)
    inCell = find(strcmp(cellKeys, uCellKeys{k}));
    if isKey(ovs_trial_overrides, uCellKeys{k})
        pickName = ovs_trial_overrides(uCellKeys{k});
        pick = inCell(strcmp({candidates(inCell).name}, pickName));
        fprintf('%s: %d candidate(s), OVERRIDDEN to %s\n', uCellKeys{k}, numel(inCell), pickName);
    elseif numel(inCell) > 1
        [~, bestLocal] = max([candidates(inCell).nOnsets] + 1e-6*[candidates(inCell).trialNum]);
        pick = inCell(bestLocal);
        others = inCell(inCell ~= pick);
        fprintf('%s: %d candidates -- using %s (%d onsets), excluding %s\n', uCellKeys{k}, numel(inCell), ...
            candidates(pick).name, candidates(pick).nOnsets, strjoin({candidates(others).name}, ', '));
    else
        pick = inCell;
    end
    selected(end+1) = struct('name', candidates(pick).name, ...
        'trialDir', fullfile(base_dir, candidates(pick).name), ...
        'intensity', candidates(pick).intensity, 'drug', candidates(pick).drug); %#ok<SAGROW>
end

% chronological drug order (drugs are applied sequentially, so first-seen
% trial number order = application order) and fixed intensity column order
drugOrderKey = arrayfun(@(c) [c.drug '|' num2str(min([candidates(strcmp({candidates.drug},c.drug)).trialNum]))], selected, 'UniformOutput', false); %#ok<NASGU>
[uDrugs, ~, ~] = unique({selected.drug}, 'stable');
firstTrialNumByDrug = cellfun(@(d) min([candidates(strcmp({candidates.drug}, d)).trialNum]), uDrugs);
[~, order] = sort(firstTrialNumByDrug);
uDrugs = uDrugs(order);
intensityOrder = {'low', 'medium', 'high'};

fprintf('\ndrug stages (row order): %s\n', strjoin(uDrugs, ' -> '));
fprintf('%d x %d grid (%d cells)\n', numel(uDrugs), numel(intensityOrder), numel(selected));

%% 2. convert each selected trial's raw tif into a z-summed video, saved in place
for i = 1:numel(selected)
    trialDir = selected(i).trialDir;
    outFile = fullfile(trialDir, 'imgData_sum.mat');
    if isfile(outFile) && ~overwrite_video
        fprintf('skip (exists): %s\n', outFile);
        continue
    end
    fprintf('[%s] reading raw tif...\n', selected(i).name);
    [imgData_sum, sp] = ovs_read_raw_tif_summed(trialDir, ovs_fallback_scan_params); %#ok<ASGLU>
    save(outFile, 'imgData_sum', 'sp', '-v7.3');
    fprintf('  saved %s  (%d x %d x %d volumes)\n', outFile, size(imgData_sum,1), size(imgData_sum,2), size(imgData_sum,3));
end

%% 2b. motion-correct each trial's summed-z video (rigid NoRMCorre via pb_register)
if do_motion_correction
    for i = 1:numel(selected)
        trialDir = selected(i).trialDir;
        outFile = fullfile(trialDir, 'imgData_sum_reg.mat');
        if isfile(outFile) && ~overwrite_reg
            fprintf('skip (exists): %s\n', outFile);
            continue
        end
        S = load(fullfile(trialDir, 'imgData_sum.mat'), 'imgData_sum');
        fprintf('[%s] motion-correcting (%d volumes)...\n', selected(i).name, size(S.imgData_sum,3));
        [imgData_sum_reg, regShifts, regDiag] = pb_register(S.imgData_sum, ...
            'smoothWindow', 5, 'maxShift', [30 45], 'verbose', false, 'figureVisible', false); %#ok<ASGLU>
        close(regDiag.figHandles);
        save(outFile, 'imgData_sum_reg', 'regShifts', '-v7.3');
        fprintf('  saved %s\n', outFile);
    end
end

%% 3. build one PB mask for this fly (mean projection across all selected trials)
maskFile = fullfile(base_dir, 'mask.mat');
if ~isfile(maskFile) || overwrite_mask
    projMean = []; nProj = 0;
    % um_per_px comes from ovs_fallback_scan_params (section 0), NOT from a
    % cached trial's saved 'sp' -- that 'sp' is whatever scan params were true
    % the LAST time that trial's imgData_sum.mat was written, which silently
    % goes stale here if you fix/refine ovs_fallback_scan_params.um_per_px
    % later without re-converting every trial's raw tif (a big re-read this
    % script deliberately avoids by caching). um_per_px is a fly-wide constant
    % anyway, so reading it from section 0's params instead is both cheap and
    % actually correct.
    umPerPx = ovs_fallback_scan_params.um_per_px(1);
    for i = 1:numel(selected)
        S = load(fullfile(selected(i).trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
        proj = mean(S.imgData_sum_reg, 3);
        if isempty(projMean); projMean = proj; else; projMean = projMean + proj; end
        nProj = nProj + 1;
    end
    projMean = projMean / nProj;

    sigma_px = max(mask_smooth_sigma_um / umPerPx, 0.5);
    open_radius_px = max(round(mask_open_radius_um / umPerPx), 1);
    join_line_px = max(round(mask_join_line_um / umPerPx), 3);
    fprintf('mask params converted to px (um_per_px=%.4f): sigma=%.2f, open_radius=%d, join_line=%d\n', ...
        umPerPx, sigma_px, open_radius_px, join_line_px);

    projSmooth = imgaussfilt(projMean, sigma_px);
    top_pct = prctile(projSmooth, 98, 'all');
    bot_pct = prctile(projSmooth, 5, 'all');
    projClipped = projSmooth;
    projClipped(projClipped > top_pct) = top_pct;
    projClipped(projClipped < bot_pct) = bot_pct;

    mask = projClipped > prctile(projClipped(:), 60);
    mask = imopen(mask, strel('disk', open_radius_px));
    mask = imfill(mask, 'holes');

    areas = sort([regionprops(mask).Area], 'descend');
    if numel(areas) < 2 || (areas(1) / areas(2)) > 1.5
        mask = bwareafilt(mask, 1);
    else
        mask = bwareafilt(mask, 2);
        mask = imclose(mask, strel('line', join_line_px, 0));
    end
    mask = imfill(mask, 'holes');

    figure(1); clf
    imagesc(projMean); axis equal tight; colormap(bone); hold on
    contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
    title([flyName ' mask QC'], 'Interpreter', 'none');
    exportgraphics(figure(1), fullfile(base_dir, 'mask_qc.png'), 'Resolution', 150);

    save(maskFile, 'mask');
    fprintf('saved %s\n', maskFile);
else
    fprintf('skip (exists): %s\n', maskFile);
end
load(maskFile, 'mask');

%% 4. per-frame background subtraction (mean outside mask, per volume) + in-stim/out-stim diff image
for i = 1:numel(selected)
    trialDir = selected(i).trialDir;
    S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    img = double(S.imgData_sum_reg);

    h5list = dir(fullfile(trialDir, '*.h5'));
    sync = ovs_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

    nVol = size(img, 3);
    if numel(sync.stims) ~= nVol
        warning('overshoot_drug_stim_script:volMismatch', ...
            '%s: h5 stim logical has %d volumes but video has %d; truncating to the shorter.', ...
            selected(i).name, numel(sync.stims), nVol);
        n = min(numel(sync.stims), nVol);
        sync.stims = sync.stims(1:n);
        img = img(:,:,1:n);
        nVol = n;
    end

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1); % 1 x nVolumes, mean OUTSIDE the mask, per frame
    img_bgsub = img - reshape(bg, 1, 1, nVol); % subtract that frame's own background from every pixel

    selected(i).in_stim  = mean(img_bgsub(:,:, sync.stims), 3);
    selected(i).out_stim = mean(img_bgsub(:,:,~sync.stims), 3);
    selected(i).diff     = selected(i).in_stim - selected(i).out_stim;
    fprintf('[%s] diff image done\n', selected(i).name);
end

%% 5. grid figure: rows = drug stage, cols = intensity
rwb = ovs_redwhiteblue(256);
diffImgs = {selected.diff};
diff_pix = prctile(abs(cat(1, diffImgs{:})), 99, 'all');

% Classic subplot + 'axis equal tight' (NOT 'axis off'/'axis image off') --
% tiledlayout+nexttile+'axis off' turned into a losing fight: 'axis off' hides
% an axes' XLabel/YLabel even set afterward, tiledlayout ignores extra figure
% width unless its OuterPosition is pinned explicitly, and pinning that then
% clipped row 1 and the sgtitle. 'axis equal tight' never hides labels in the
% first place (ticks/box stay, stripped below via XTick/YTick instead), and
% subplot's own margins are generous enough that row labels just work.
figure(2); clf
nRows = numel(uDrugs); nCols = numel(intensityOrder);
set(gcf, 'Position', [50 50 300*nCols 260*nRows], 'Color', 'w')

for r = 1:nRows
    for c = 1:nCols
        i = find(strcmp({selected.drug}, uDrugs{r}) & strcmp({selected.intensity}, intensityOrder{c}), 1);
        subplot(nRows, nCols, (r-1)*nCols + c)
        if isempty(i)
            axis equal tight
            text(0.5, 0.5, 'missing', 'Units', 'normalized', 'HorizontalAlignment', 'center')
            continue
        end
        imagesc(selected(i).diff); clim([-diff_pix diff_pix]); axis equal tight
        colormap(gca, rwb)
        if r == 1
            title(intensityOrder{c}, 'Interpreter', 'none')
        end
        if c == 1
            ylabel(uDrugs{r}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'Interpreter', 'none', 'FontWeight', 'bold')
        end
    end
end
set(findobj(gcf, 'Type', 'Axes'), 'XTick', [], 'YTick', [])
sgtitle([flyName ': in-out diff (background-subtracted), rows=drug stage, cols=intensity (red=brighter in stim, blue=dimmer)'], 'Interpreter', 'none')

exportgraphics(figure(2), fullfile(base_dir, 'drugStim_diffGrid.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(base_dir, 'drugStim_diffGrid.png'));

%% 6. per-repetition diff images: rows=condition, cols=repetition, one figure per fly
% Unlike section 5 (which lumps ALL in-stim/out-of-stim volumes across the
% whole trial into one flat average), each cell here is exactly ONE stim
% repetition: mean(its own 2 s pulse) - mean(the 8 s immediately PRECEDING
% that pulse's onset) -- not the trial-wide out-of-stim mean. Since the
% protocol is a clean 10 s period with a 2 s pulse (checked directly against
% every trial's h5 so far, always 12 pulses/trial), that preceding 8 s window
% is exactly the inter-pulse gap right before this repetition, matching "the
% preceding outside-of-stim time" as asked. Rows are (drug,intensity)
% conditions in the same order as section 5's grid, flattened into one list;
% columns are repetition number.
pre_stim_sec = 8;
stim_dur_sec = 2;

rowOrder = [];
rowLabels = {};
for r = 1:numel(uDrugs)
    for c = 1:numel(intensityOrder)
        i = find(strcmp({selected.drug}, uDrugs{r}) & strcmp({selected.intensity}, intensityOrder{c}), 1);
        if ~isempty(i)
            rowOrder(end+1) = i; %#ok<AGROW>
            rowLabels{end+1} = sprintf('%s %s', uDrugs{r}, intensityOrder{c}); %#ok<AGROW>
        end
    end
end

repDiffs = cell(numel(rowOrder), 1); % repDiffs{row} = [Y x X x nReps for that trial]
nRepsPerRow = zeros(numel(rowOrder), 1);
for ri = 1:numel(rowOrder)
    i = rowOrder(ri);
    trialDir = selected(i).trialDir;
    S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    img = double(S.imgData_sum_reg);

    h5list = dir(fullfile(trialDir, '*.h5'));
    sync = ovs_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

    nVol = size(img, 3);
    if numel(sync.stims) ~= nVol
        n = min(numel(sync.stims), nVol);
        sync.stims    = sync.stims(1:n);
        sync.t_volume = sync.t_volume(1:n);
        img = img(:,:,1:n);
        nVol = n;
    end

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);

    onsetIdx   = find(diff(sync.stims) > 0) + 1;
    onsetTimes = sync.t_volume(onsetIdx);

    thisReps = nan(size(img,1), size(img,2), numel(onsetTimes));
    for k = 1:numel(onsetTimes)
        t0 = onsetTimes(k);
        inIdx  = find(sync.t_volume >= t0 & sync.t_volume < t0 + stim_dur_sec);
        preIdx = find(sync.t_volume >= t0 - pre_stim_sec & sync.t_volume < t0);
        if isempty(inIdx) || isempty(preIdx)
            continue % pulse too close to the trial edge for a full window -- leave as NaN
        end
        thisReps(:,:,k) = mean(img_bgsub(:,:,inIdx), 3) - mean(img_bgsub(:,:,preIdx), 3);
    end
    repDiffs{ri} = thisReps;
    nRepsPerRow(ri) = numel(onsetTimes);
    fprintf('[%s] %d repetitions windowed\n', selected(i).name, numel(onsetTimes));
end

nRepsMax = max(nRepsPerRow);
allVals  = cellfun(@(x) x(:), repDiffs, 'UniformOutput', false);
allVals  = cat(1, allVals{:});
diff_pix2 = prctile(abs(allVals(~isnan(allVals))), 99, 'all');

figure(3); clf
nRowsRep = numel(rowOrder);
set(gcf, 'Position', [50 50 90*nRepsMax 140*nRowsRep], 'Color', 'w')
for ri = 1:nRowsRep
    for k = 1:nRepsPerRow(ri)
        subplot(nRowsRep, nRepsMax, (ri-1)*nRepsMax + k)
        d = repDiffs{ri}(:,:,k);
        if all(isnan(d(:)))
            axis equal tight
            continue
        end
        imagesc(d); clim([-diff_pix2 diff_pix2]); axis equal tight
        colormap(gca, rwb)
        if ri == 1
            title(sprintf('rep %d', k))
        end
        if k == 1
            ylabel(rowLabels{ri}, 'Rotation', 0, 'HorizontalAlignment', 'right', 'Interpreter', 'none', 'FontSize', 8, 'FontWeight', 'bold')
        end
    end
end
set(findobj(gcf, 'Type', 'Axes'), 'XTick', [], 'YTick', [])
sgtitle([flyName ': per-repetition in-out diff (2s stim vs its own preceding 8s), rows=condition, cols=repetition'], 'Interpreter', 'none')

exportgraphics(figure(3), fullfile(base_dir, 'perRepDiffGrid.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(base_dir, 'perRepDiffGrid.png'));

%% 7. divide the PB into glomeruli (dF/F normalization, troubleshooting pass)
% Sections 4-6 above show raw background-subtracted F differences, which is
% fine pixel-by-pixel as long as motion correction is perfect, but any residual
% motion turns a real dF/F into a spurious edge artifact when done per-pixel.
% Averaging within a glomerulus (a chunk of the PB arch, not a single pixel)
% is much more robust to a few pixels of residual jitter, so from here on
% normalization happens per-glomerulus instead of per-pixel.
%
% Same clustering approach as led_stim_berg4_script.m's section 7 (see that
% script's lsb_pb_glomeruli for the fuller explanation): skeletonize the mask,
% order the skeleton into one continuous path (graph_sort.m), resample it into
% evenly spaced centroids, then assign every mask pixel to its nearest
% centroid. n_per_hemisphere (section 0) is doubled internally -> 40 total.
clusterIdx = ovs_pb_glomeruli(mask, n_per_hemisphere);
nClusters  = 2 * n_per_hemisphere;

% QC figure: same numbered-path-over-mean-image overlay as led_stim_berg4's
% mask QC, so a bad skeleton/ordering (crossed lines, doubled-back numbering)
% is visible here instead of only showing up as a scrambled heatmap later.
qcCentroids = zeros(nClusters, 2);
for c = 1:nClusters
    [cy, cx] = find(clusterIdx == c);
    qcCentroids(c,:) = [mean(cx), mean(cy)];
end
S = load(fullfile(selected(1).trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
projMean = mean(S.imgData_sum_reg, 3);

figure(4); clf
imagesc(projMean); axis equal tight; colormap(bone); hold on
contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
plot(qcCentroids(:,1), qcCentroids(:,2), 'y-', 'LineWidth', 1)
plot(qcCentroids(:,1), qcCentroids(:,2), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 4)
for c = 1:nClusters
    text(qcCentroids(c,1), qcCentroids(c,2), num2str(c), 'Color', 'g', ...
        'FontSize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end
title([flyName ': glomerulus numbering (' num2str(nClusters) ' total)'], 'Interpreter', 'none')
exportgraphics(figure(4), fullfile(base_dir, 'glomeruli_qc.png'), 'Resolution', 150);
fprintf('saved %s\n', fullfile(base_dir, 'glomeruli_qc.png'));

%% 8. per-glomerulus dF/F, full trial and peri-stim -- every trial for this fly
% Troubleshot on the baseline trial only at first (see git history / prior
% run); now expanded to every trial now that the low-signal-glomerulus NaN
% exclusion (glom_min_snr, section 0) checked out.
trialsToShow = 1:numel(selected);

for i = trialsToShow
    trialDir = selected(i).trialDir;
    S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    img = double(S.imgData_sum_reg);

    h5list = dir(fullfile(trialDir, '*.h5'));
    sync = ovs_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

    nVol = size(img, 3);
    if numel(sync.stims) ~= nVol
        n = min(numel(sync.stims), nVol);
        sync.stims    = sync.stims(1:n);
        sync.t_volume = sync.t_volume(1:n);
        img = img(:,:,1:n);
        nVol = n;
    end

    % step 1 (per the user's ask): PMT/background offset subtraction, same
    % per-frame convention as section 4 -- NOT a per-pixel dF/F, just removing
    % the frame's own mean-outside-mask offset before any averaging.
    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);
    img_bgsub_2d = reshape(img_bgsub, [], nVol);

    % step 2: mean fluorescence per glomerulus per volume (nClusters x nVolumes)
    centroidLog = false(nClusters, size(img_bgsub_2d,1));
    for c = 1:nClusters
        centroidLog(c, clusterIdx(:) == c) = true;
    end
    f_cluster = double(centroidLog) * img_bgsub_2d ./ sum(centroidLog, 2);

    % F0 = this glomerulus's own mean OUT-OF-STIM fluorescence for this trial
    % (same convention as led_stim_berg4_script.m's section 7 -- baselining to
    % a low percentile would make dFF positive almost everywhere by
    % construction and hide real inhibition).
    f0_cluster  = mean(f_cluster(:, ~sync.stims), 2);
    dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;

    std_cluster = std(f_cluster(:, ~sync.stims), 0, 2);
    badGlom = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
    if any(badGlom)
        fprintf('  [%s] excluding %d/%d low-signal glomeruli from dF/F (F0 too close to noise floor): %s\n', ...
            selected(i).name, sum(badGlom), nClusters, mat2str(find(badGlom)'));
    end
    dff_cluster(badGlom, :) = NaN;

    onsetIdx   = find(diff(sync.stims) > 0) + 1;
    offsetIdx  = find(diff(sync.stims) < 0) + 1;
    onsetTimes = sync.t_volume(onsetIdx);
    offsetTimes = sync.t_volume(offsetIdx);

    figure(5); clf
    set(gcf, 'Position', [50 50 1400 500], 'Color', 'w')
    % imagesc maps NaN to the colormap's first entry (solid blue here) instead
    % of leaving it blank -- AlphaData + a gray axes color renders excluded
    % (badGlom) glomeruli as gray instead of a misleadingly extreme color.
    imagesc(sync.t_volume, 1:nClusters, dff_cluster, 'AlphaData', double(~isnan(dff_cluster)))
    set(gca, 'Color', [0.7 0.7 0.7])
    dff_pix = prctile(abs(dff_cluster(:)), 99);
    clim([-dff_pix dff_pix]); colormap(gca, rwb); colorbar
    hold on
    for k = 1:numel(onsetTimes)
        xline(onsetTimes(k), 'k-', 'LineWidth', 0.75);
        if k <= numel(offsetTimes)
            xline(offsetTimes(k), 'k--', 'LineWidth', 0.75);
        end
    end
    axis tight
    xlabel('time in trial (s)'); ylabel('glomerulus')
    title([flyName ' - ' selected(i).name ': per-glomerulus dF/F, full trial (solid=stim on, dashed=stim off)'], 'Interpreter', 'none')

    outFile = fullfile(base_dir, sprintf('glomHeatmap_fullTrial_%s.png', selected(i).name));
    exportgraphics(figure(5), outFile, 'Resolution', 200);
    fprintf('saved %s\n', outFile);
end

%% 9. per-glomerulus dF/F, averaged across stim repetitions, aligned to onset
% Same window logic as section 6 (clean 10s cycle: 8s gap + 2s pulse) but here
% averaging dF/F traces across all repetitions instead of taking single-pulse
% diff images -- shows the full time-course, not just a before/after snapshot.
peristim_pre_sec  = 4;
peristim_post_sec = 6;

for i = trialsToShow
    trialDir = selected(i).trialDir;
    S = load(fullfile(trialDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg');
    img = double(S.imgData_sum_reg);

    h5list = dir(fullfile(trialDir, '*.h5'));
    sync = ovs_h5_stim_timing(fullfile(h5list(1).folder, h5list(1).name));

    nVol = size(img, 3);
    if numel(sync.stims) ~= nVol
        n = min(numel(sync.stims), nVol);
        sync.stims    = sync.stims(1:n);
        sync.t_volume = sync.t_volume(1:n);
        img = img(:,:,1:n);
        nVol = n;
    end

    img_2d = reshape(img, [], nVol);
    bg = mean(img_2d(~mask(:), :), 1);
    img_bgsub = img - reshape(bg, 1, 1, nVol);
    img_bgsub_2d = reshape(img_bgsub, [], nVol);

    centroidLog = false(nClusters, size(img_bgsub_2d,1));
    for c = 1:nClusters
        centroidLog(c, clusterIdx(:) == c) = true;
    end
    f_cluster = double(centroidLog) * img_bgsub_2d ./ sum(centroidLog, 2);
    f0_cluster  = mean(f_cluster(:, ~sync.stims), 2);
    dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;

    std_cluster = std(f_cluster(:, ~sync.stims), 0, 2);
    badGlom = f0_cluster <= 0 | (f0_cluster ./ std_cluster) < glom_min_snr;
    dff_cluster(badGlom, :) = NaN;

    fs_vol  = 1 / mean(diff(sync.t_volume));
    winSamp = round(-peristim_pre_sec*fs_vol) : round(peristim_post_sec*fs_vol);
    onsetIdx = find(diff(sync.stims) > 0) + 1;
    valid    = onsetIdx + winSamp(1) >= 1 & onsetIdx + winSamp(end) <= nVol;
    onsetIdx = onsetIdx(valid);

    slices = nan(numel(onsetIdx), nClusters, numel(winSamp));
    for o = 1:numel(onsetIdx)
        slices(o,:,:) = dff_cluster(:, onsetIdx(o) + winSamp);
    end
    avgHeatmap = squeeze(mean(slices, 1, 'omitnan'));

    figure(6); clf
    set(gcf, 'Position', [50 50 900 500], 'Color', 'w')
    t_win = winSamp / fs_vol;
    imagesc(t_win, 1:nClusters, avgHeatmap, 'AlphaData', double(~isnan(avgHeatmap)))
    set(gca, 'Color', [0.7 0.7 0.7])
    dff_pix2 = prctile(abs(avgHeatmap(:)), 99);
    clim([-dff_pix2 dff_pix2]); colormap(gca, rwb); colorbar
    hold on
    xline(0, 'k-', 'LineWidth', 1); xline(2, 'k--', 'LineWidth', 1); % stim on/off, 2s pulse
    axis tight
    xlabel('time from stim onset (s)'); ylabel('glomerulus')
    title([flyName ' - ' selected(i).name ': per-glomerulus dF/F, mean over ' num2str(numel(onsetIdx)) ' repetitions'], 'Interpreter', 'none')

    outFile = fullfile(base_dir, sprintf('glomHeatmap_periStim_%s.png', selected(i).name));
    exportgraphics(figure(6), outFile, 'Resolution', 200);
    fprintf('saved %s\n', outFile);
end

%% local functions

function d = lsb_list_subdirs_ovs(parentDir)
d = dir(parentDir);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
end

function n = ovs_count_tif_frames(tifPath)
t = Tiff(tifPath, 'r');
cleanupTiff = onCleanup(@() close(t)); %#ok<NASGU>
n = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    n = n + 1;
end
end

function sp = ovs_scan_params(trialDir, fallback)
jsonFile = fullfile(trialDir, 'scan_params.json');
if isfile(jsonFile)
    sp = jsondecode(fileread(jsonFile));
    sp.provenance = 'scan_params.json';
else
    sp = fallback;
    sp.provenance = 'FALLBACK (scan_params.json missing -- copied from trial_012)';
    fprintf('  note: %s has no scan_params.json; using fallback scan params (see script header).\n', trialDir);
end
end

function [imgSum, sp] = ovs_read_raw_tif_summed(trialDir, fallbackSp)
sp = ovs_scan_params(trialDir, fallbackSp);
if sp.nchannels ~= 1
    error('overshoot_drug_stim_script:multiChannel', '%s has %d saved channels; this script assumes 1.', trialDir, sp.nchannels);
end
validPlanes = setdiff(1:sp.nplanes, sp.flyback_planes(:)' + 1); % json flyback is 0-indexed
if numel(validPlanes) ~= sp.nplanes_valid
    error('overshoot_drug_stim_script:flybackMismatch', ...
        '%s: nplanes_valid=%d but nplanes/flyback_planes leaves %d.', trialDir, sp.nplanes_valid, numel(validPlanes));
end

tifList = dir(fullfile(trialDir, '*.tif'));
if numel(tifList) ~= 1
    error('overshoot_drug_stim_script:tif', 'Expected exactly one .tif in %s, found %d.', trialDir, numel(tifList));
end
tifPath = fullfile(tifList(1).folder, tifList(1).name);

nFramesOnDisk = ovs_count_tif_frames(tifPath);
framesPerVol = sp.nplanes * sp.nchannels;
nVolumes = floor(nFramesOnDisk / framesPerVol);
leftover = nFramesOnDisk - nVolumes * framesPerVol;
if leftover > 0
    fprintf('  note: %d leftover raw frame(s) after %d complete volumes, dropping the trailing partial volume\n', leftover, nVolumes);
end

t = Tiff(tifPath, 'r');
cleanupTiff = onCleanup(@() close(t)); %#ok<NASGU>
t.setDirectory(1);
nRawToRead = nVolumes * framesPerVol;
raw = zeros(sp.px_height, sp.px_width, nRawToRead, 'int16');
tic
for k = 1:nRawToRead
    raw(:,:,k) = t.read();
    if k < nRawToRead; t.nextDirectory(); end
    if mod(k, 20000) == 0 || k == nRawToRead
        fprintf('  read %d/%d raw frames (%.0f frames/s)\n', k, nRawToRead, k/toc);
    end
end

raw = reshape(raw, sp.px_height, sp.px_width, sp.nplanes, nVolumes);
imgSum = squeeze(sum(raw(:,:,validPlanes,:), 3));
sp.nVolumes = nVolumes;
end

function s = ovs_h5_stim_timing(h5Path)
% Same idea as led_stim_berg4_script.m's lsb_h5_stim_timing, but returns the
% stim logical already downsampled onto the volume clock (si_volumeclock),
% since that's the only resolution this script needs (no fine-DAQ-clock use).
info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
if numel(sweepNames) ~= 1
    error('overshoot_drug_stim_script:h5sweep', 'Expected exactly one sweep group in %s, found %d.', h5Path, numel(sweepNames));
end

fs = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
digi = int32(h5read(h5Path, [sweepNames{1} '/digitalScans']));

volCh = find(diNames == "si_volumeclock", 1);
ledCh = find(diNames == "led_stim_feedback", 1);
if isempty(volCh) || isempty(ledCh)
    error('overshoot_drug_stim_script:h5chan', ...
        'Could not find si_volumeclock/led_stim_feedback in %s DIChannelNames (found: %s).', h5Path, strjoin(diNames, ', '));
end

volBit = bitget(digi, volCh);
ledBit = bitget(digi, ledCh);
N = numel(digi);
t_fine = (0:N-1)' / fs;
t_volume = t_fine(find(diff(volBit) > 0) + 1);

s.fs = fs;
s.durationSec = N / fs; % real trial duration -- NOT numel(s.stims)/fs, since s.stims is per-volume, not per-DAQ-sample
s.t_volume = t_volume; % real per-volume timestamps (DAQ seconds) -- needed to window individual repetitions precisely
s.stims = interp1(t_fine, double(ledBit), t_volume, 'nearest', 'extrap') > 0.5; % per-volume logical
end

function cmap = ovs_redwhiteblue(n)
half = n / 2;
blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
whiteToRed  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
cmap = [blueToWhite; whiteToRed];
end

function clusterIdx = ovs_pb_glomeruli(mask, nPerHemisphere)
% Same skeletonize -> order -> resample-into-centroids -> nearest-centroid-
% assign approach as led_stim_berg4_script.m's lsb_pb_glomeruli (itself
% trimmed from this lab's process_im, e.g. epg_dlight_script.m,
% dopamine_ionto_script.m). Requires graph_sort.m (repo root) to order the
% skeleton into a path -- see its own comments for why this assumes mask is a
% single open arch, not a closed loop.
%
% clusterIdx: same [Y,X] size as mask; 0 outside the mask, 1..2*nPerHemisphere
% inside, numbered along the skeleton from one end to the other.

[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('ovs_pb_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
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
centroids = centroids(2:2:end-1, :); % drop the edge samples so every cluster is the same size

[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);

clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end
