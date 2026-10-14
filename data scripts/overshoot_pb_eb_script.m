%% overshoot_pb_eb_script
% Two-structure version of overshoot_bar_stim_script.m for volumes that
% span both the protocerebral bridge (PB) and, deeper, the ellipsoid body
% (EB): 20261012-1_epg_syt8s_lpsp_cschrimson onwards (24 frames/volume =
% 20 slices at 25-120 um + 4 flyback, 64 x 128 px, zoom 3, 8.25 Hz).
%
% PB (planes 1-7, z 25-55 um): exactly the single-structure pipeline --
%   z-sum, rigid registration, arch mask, 2N glomeruli along the skeleton
%   (N per hemisphere, alpha repeated over the two halves), z-score, PVA.
% EB (planes 10-20, z 70-120 um): the ring is imaged EDGE-ON in this prep --
%   its two walls run vertically through z (planes x columns projection is a
%   ring; the x-y projection is a 28 x 10 px blob), so summing the planes
%   would collapse one ring diameter. The EB planes are therefore kept
%   SEPARATE (cropped to rows 1-40, cols 25-100), registered with the PB's
%   shifts (same volumes; the EB's own rigid registration ran away on 14 % of
%   volumes because the structure is small), masked in 3-D, and split into N
%   angular wedges by the angle of each voxel in the x-z plane around the
%   ring centre (claude\pb_eb_wedges3d). The PVA is then the PB recipe on
%   the N wedges (alpha once around the ring).
% Both: per-volume background subtraction (mean outside the structure's own
%   mask), z-scored activity heatmap, bump vs bar with sign + lag fitted on
%   the sweeps, bump_results_PB.mat / bump_results_EB.mat (same fields as
%   bump_results.mat), PB-vs-EB comparison, stim on-off diff images.
% Sync: claude\pb_h5_sync (LED from the analog vr_led_stim_copy line; this
%   rig configuration has no led_stim_feedback DI line).
%
% ---- trial structure, 20261012-1 trial_002 (dangling_burn_in_8pos_ledstim
% _pre5s_rev_endburnin_826s) ----
%   0-185 s VR: burn-in sweeps at 36 deg/s with 10-s holds; then 16 holds of
%   25 s at 45-deg steps (8 positions forward, then the same 8 in reverse),
%   each = 5 s pre + 10 s LED ON + 10 s post, separated by 405-deg sweeps;
%   then two more burn-in sweeps. LED: 16 pulses x 10 s (DAQ 232.7-794.0 s).
%   DAQ: 6900 volume clocks = 165600 frames / 24.
%
% ---- 20261012-2 trial_002 (190-s burn-in, 1585 volumes) lessons ----
%   * The right lateral PB tip is ~1/3 as bright as the rest and came out as a
%     separate blob at the mask threshold; the old "keep the largest blob"
%     dropped it (glomeruli 29-32 missing). Section 3 now builds a hysteresis
%     mask: all blobs >= 5 % of the largest at the tight threshold are seeds,
%     gaps between seeds are filled at a threshold 15 pct-points lower.
%   * NoRMCorre locked onto spurious -15..-20 px peaks for ~7 volumes every
%     ~45 volumes (16 % of frames; true motion < 2 px), painting a shifted
%     ghost of the PB (and, via the shared shifts, the EB) into the mean.
%     Fix: PB maxShift [6 8] here + pb_register's outlierJumpPx cleanup
%     (31-frame median reference). Check the printed shift range (should be
%     ~ +/-2 px) and the mean image for a faint duplicate arch.
%
% ---- 20261012-3 trial_001 lesson ----
%   The left hemisphere is ~1/3 as bright as the right, and a bright non-PB
%   blob sits just below the bridge midline. The percentile + disk-closing
%   mask bridged the blob into the arch (glomeruli 19-20 ended up on it and
%   the numbering wrapped). PB masking moved to claude\pb_arch_mask.m:
%   column-normalised threshold for the dim arm, components kept only if
%   they are a second hemisphere or within 1.5 half-widths of the main arch
%   skeleton (and >= 20 % as bright), gap fill restricted to the arch.
%   Verified on flies -1, -2, -3 of 20261012.
%
% Run one %% section at a time, or headless:
%   matlab -batch "obs2_trialDir='Z:\...\trial_002\'; run('overshoot_pb_eb_script.m')"

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'claude'));
if isempty(which('graph_sort')), addpath(fullfile(thisDir, '..')); end
if isempty(which('circ_mean')), addpath(fullfile(thisDir, '..', 'circ_stats')); end

trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20261012-1_epg_syt8s_lpsp_cschrimson\trial_001\';
if exist('obs2_trialDir', 'var') && ~isempty(obs2_trialDir), trialDir = obs2_trialDir; end
if trialDir(end) ~= '\' && trialDir(end) ~= '/', trialDir(end+1) = filesep; end
parts = strsplit(trialDir, filesep); parts = parts(~cellfun(@isempty, parts));
trialName = [parts{end-1} filesep parts{end}];

n_per_hemisphere = 16;              % PB: 2N glomeruli; EB: N wedges
S = struct('name', {'PB', 'EB'}, ...
    'planes', {1:7, 10:20}, ...     % 1-based frame indices within the volume
    'crop', {[], [1 40 25 100]}, ...% EB sub-volume [row1 row2 col1 col2] (keeps the 4-D stack < 1 GB)
    'mask_smooth_sigma_um', {1.5, 1.2}, ...
    'mask_thresh_pct', {88, 92}, ...% percentile of the smoothed mean image/volume kept as mask
    'mask_open_radius_um', {1.0, 0}, ...
    'mask_join_line_um', {15, 0}, ...
    'maxShift', {[6 8], []}, ...    % PB: real motion in these preps is < 2 px; a wide search let NoRMCorre jump to spurious +/-15-20 px peaks (ghost PB in the mean). EB: not estimated, PB shifts applied
    'nClusters', {2*n_per_hemisphere, n_per_hemisphere});
if exist('obs2_planes', 'var') && ~isempty(obs2_planes), S(1).planes = obs2_planes{1}; S(2).planes = obs2_planes{2}; end

min_hold_sec = 4; hold_flat_tol_V = 0.1;
lag_range_sec = [-2 5]; rho_thresh = 0.1; max_ok_shift_px = 8;
% Cache control. Sections 1-3 skip work whose output file already exists in
% the trial folder; set these true (or the obs2_* base-workspace overrides
% when running headless) to recompute after a pipeline change:
%   overwrite_mask  -> masks + segmentation (after a pb_arch_mask / wedge change; cheap)
%   overwrite_reg   -> registration (also set overwrite_mask; ~10 min)
%   overwrite_video -> re-read the tif (plane assignment / crop change; also set the two above)
% Sections 4-7 (sync, heatmaps, PVA, figures) are always recomputed.
overwrite_video = true; overwrite_reg = true; overwrite_mask = true;
if exist('obs2_overwrite_video', 'var') && ~isempty(obs2_overwrite_video), overwrite_video = obs2_overwrite_video; end
if exist('obs2_overwrite_reg', 'var') && ~isempty(obs2_overwrite_reg), overwrite_reg = obs2_overwrite_reg; end
if exist('obs2_overwrite_mask', 'var') && ~isempty(obs2_overwrite_mask), overwrite_mask = obs2_overwrite_mask; end
if overwrite_video, overwrite_reg = true; end   % downstream caches are invalid if an upstream one is redone
if overwrite_reg, overwrite_mask = true; end

% scan geometry from the tif header
tifList = dir(fullfile(trialDir, '*.tif')); assert(numel(tifList) == 1, 'expected one tif in %s', trialDir);
tifPath = fullfile(tifList(1).folder, tifList(1).name);
fid = fopen(tifPath, 'r'); hdr = fread(fid, 80000, 'uint8=>char')'; fclose(fid);
getNum = @(key) str2double(regexp(hdr, ['SI\.' regexptranslate('escape', key) '\s*=\s*([-0-9.eE]+)'], 'tokens', 'once'));
framesPerVol = getNum('hStackManager.numFramesPerVolumeWithFlyback');
zoom = getNum('hRoiManager.scanZoomFactor'); W = getNum('hRoiManager.pixelsPerLine'); H = getNum('hRoiManager.linesPerFrame');
zsStr = regexp(hdr, 'SI\.hStackManager\.zs = \[([^\]]*)\]', 'tokens', 'once'); zs = str2num(zsStr{1}); %#ok<ST2NM>
zStep = median(diff(zs));
umPerPx = 0.3164 * (4 / zoom) * (256 / W);  % nominal; 0.3164 um/px at zoom 4 with 256 px across the FOV
fprintf('[%s] %dx%d px, %d frames/volume, zoom %g -> %.3f um/px (nominal); slices z = %s um\n', trialName, H, W, framesPerVol, zoom, umPerPx, mat2str(zs));
for k = 1:2, fprintf('  %s planes %s (z %d-%d um)\n', S(k).name, mat2str(S(k).planes), zs(S(k).planes(1)), zs(S(k).planes(end))); end

%% 1. raw tif -> PB z-sum movie + EB per-plane sub-volume (cached)
vidFiles = arrayfun(@(s) fullfile(trialDir, sprintf('imgData_%s.mat', s.name)), S, 'UniformOutput', false);
if overwrite_video || ~all(cellfun(@isfile, vidFiles))
    fprintf('[%s] reading raw tif (PB z-sum + EB planes)...\n', trialName);
    o = pb_read_raw_tif_planes(tifPath, framesPerVol, {S.planes}, 'separate', [false true], 'crop', {S.crop}, 'planeMeans', true);
    for k = 1:2
        imgData = o.stacks{k}; planes = S(k).planes; nVolumes = o.nVolumes; crop = S(k).crop; rows = o.rows{k}; cols = o.cols{k}; %#ok<NASGU>
        save(vidFiles{k}, 'imgData', 'planes', 'nVolumes', 'crop', 'rows', 'cols', '-v7.3');
        fprintf('  saved %s %s\n', vidFiles{k}, mat2str(size(imgData)));
    end
    planeMean = o.planeMean; save(fullfile(trialDir, 'planeMeans.mat'), 'planeMean', 'zs'); %#ok<NASGU>
    figure(9); clf; set(gcf, 'Position', [50 50 1800 900], 'Color', 'w')
    for k = 1:framesPerVol
        subplot(4, 6, k); imagesc(o.planeMean(:,:,k)); axis image off; colormap gray
        if k <= numel(zs), ttl = sprintf('plane %d, z = %d um', k, zs(k)); else, ttl = sprintf('plane %d (flyback)', k); end
        if ismember(k, S(1).planes), ttl = [ttl ' [PB]']; elseif ismember(k, S(2).planes), ttl = [ttl ' [EB]']; end
        title(ttl, 'FontSize', 8)
    end
    sgtitle(sprintf('%s: per-plane mean image; planes used for PB / EB', trialName), 'Interpreter', 'none')
    exportgraphics(figure(9), fullfile(trialDir, 'planeMeans.png'), 'Resolution', 110);
    clear o imgData
else
    fprintf('skip (exists): %s\n', strjoin(vidFiles, ', '));
end

%% 2. motion correction: PB estimated, EB = PB shifts applied plane by plane (cached)
regFiles = arrayfun(@(s) fullfile(trialDir, sprintf('imgData_%s_reg.mat', s.name)), S, 'UniformOutput', false);
if overwrite_reg || ~isfile(regFiles{1})
    V = load(vidFiles{1}, 'imgData');
    fprintf('[%s] motion-correcting PB (%d volumes)...\n', trialName, size(V.imgData, 3));
    [imgData_reg, regShifts] = pb_register(V.imgData, 'smoothWindow', 5, 'maxShift', S(1).maxShift, 'verbose', false, 'figureVisible', false);
    sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), regShifts(:), 'UniformOutput', false));
    fprintf('  PB shift range row = [%.1f %.1f] px, col = [%.1f %.1f] px\n', min(sh(:,1)), max(sh(:,1)), min(sh(:,2)), max(sh(:,2)));
    save(regFiles{1}, 'imgData_reg', 'regShifts', '-v7.3'); clear V imgData_reg
else
    fprintf('skip (exists): %s\n', regFiles{1});
end
if overwrite_reg || ~isfile(regFiles{2})
    R = load(regFiles{1}, 'regShifts'); regShifts = R.regShifts; clear R
    V = load(vidFiles{2}, 'imgData'); [nR, nC, nP, nV] = size(V.imgData);
    fprintf('[%s] applying the PB shifts to the %d EB planes (%d volumes)...\n', trialName, nP, nV);
    imgData_reg = zeros(nR, nC, nP, nV, 'single');
    for pIdx = 1:nP
        imgData_reg(:, :, pIdx, :) = reshape(pb_register(squeeze(V.imgData(:, :, pIdx, :)), 'applyShifts', regShifts, 'smoothWindow', 0, ...
            'maxShift', S(1).maxShift, 'verbose', false, 'figureVisible', false), nR, nC, 1, nV);
    end
    save(regFiles{2}, 'imgData_reg', 'regShifts', '-v7.3'); clear V imgData_reg regShifts
    fprintf('  saved %s\n', regFiles{2});
else
    fprintf('skip (exists): %s\n', regFiles{2});
end

%% 3. masks + segmentation (cached mask), QC figure
maskFiles = arrayfun(@(s) fullfile(trialDir, sprintf('mask_%s.mat', s.name)), S, 'UniformOutput', false);
% ---- PB: arch mask + skeleton glomeruli ----
V = load(regFiles{1}, 'imgData_reg'); S(1).projMean = mean(V.imgData_reg, 3); clear V
if overwrite_mask || ~isfile(maskFiles{1})
    sigma_px = max(S(1).mask_smooth_sigma_um / umPerPx, 0.5); open_px = max(round(S(1).mask_open_radius_um / umPerPx), 1);
    ps = imgaussfilt(S(1).projMean, sigma_px); pc = min(max(ps, prctile(ps, 5, 'all')), prctile(ps, 99, 'all'));
    % claude\pb_arch_mask: global + column-normalised threshold (keeps a dim
    % hemisphere), components far from the main arch's skeleton rejected
    % (20261012-3: a bright blob below the bridge used to get bridged in and
    % wreck the glomerulus order), gaps filled only along the arch.
    [mask, maskCore, pbInfo] = pb_arch_mask(S(1).projMean, umPerPx, 'sigma_um', S(1).mask_smooth_sigma_um, ...
        'thresh_pct', S(1).mask_thresh_pct, 'open_um', S(1).mask_open_radius_um, 'join_um', S(1).mask_join_line_um);
    fprintf('PB mask: %s\n', pbInfo.summary);
    save(maskFiles{1}, 'mask', 'maskCore', 'pbInfo'); fprintf('saved %s (%d px)\n', maskFiles{1}, nnz(mask));
else, fprintf('skip (exists): %s\n', maskFiles{1}); end
M = load(maskFiles{1}); S(1).mask = M.mask; S(1).maskCore = M.maskCore;
S(1).clusterIdx = obs_pb_glomeruli(S(1).mask, n_per_hemisphere);
S(1).alpha = repmat(linspace(-pi, pi - 2*pi/n_per_hemisphere, n_per_hemisphere), 1, 2);
% ---- EB: 3-D mask on the mean sub-volume, wedges in the x-z plane ----
V = load(regFiles{2}, 'imgData_reg'); S(2).meanVol = mean(V.imgData_reg, 4); clear V      % rows x cols x planes
if overwrite_mask || ~isfile(maskFiles{2})
    sig = max(S(2).mask_smooth_sigma_um / umPerPx, 0.5);
    vs = imgaussfilt3(S(2).meanVol, [sig sig 0.6]);
    vc = min(max(vs, prctile(vs, 5, 'all')), prctile(vs, 99.5, 'all'));
    mask3 = vc > prctile(vc(:), S(2).mask_thresh_pct);
    cc = bwconncomp(mask3, 26); [~, iBig] = max(cellfun(@numel, cc.PixelIdxList));
    mask3 = false(size(mask3)); mask3(cc.PixelIdxList{iBig}) = true;
    mask3 = imclose(mask3, strel('sphere', 1));
    save(maskFiles{2}, 'mask3'); fprintf('saved %s (%d voxels)\n', maskFiles{2}, nnz(mask3));
else, fprintf('skip (exists): %s\n', maskFiles{2}); end
M = load(maskFiles{2}); S(2).mask3 = M.mask3;
[S(2).clusterIdx, S(2).alpha, ebInfo] = pb_eb_wedges3d(S(2).mask3, n_per_hemisphere, 'zScale', zStep / umPerPx);
S(2).mask = squeeze(any(S(2).mask3, 3)); S(2).maskCore = S(2).mask3;   % x-y footprint for display; core = the 3-D mask
S(2).projMean = sum(S(2).meanVol, 3);
fprintf('EB: ring centre col %.1f (image col %.1f), plane %.1f (z %.0f um); semi-axes %.1f px x %.1f planes (%.0f x %.0f um nominal); %d-%d voxels per wedge\n', ...
    ebInfo.center(1), ebInfo.center(1) + S(2).crop(3) - 1, ebInfo.center(2), interp1(1:numel(S(2).planes), zs(S(2).planes), ebInfo.center(2)), ...
    ebInfo.semiAxes, ebInfo.semiAxes(1) * umPerPx, ebInfo.semiAxes(2) * zStep, min(ebInfo.nVox), max(ebInfo.nVox));
% ---- QC figure ----
figure(1); clf; set(gcf, 'Position', [50 50 1800 800], 'Color', 'w')
subplot(2, 3, 1); imagesc(S(1).projMean); axis image; colormap(gca, bone); hold on
contour(S(1).mask, [0.5 0.5], 'r', 'LineWidth', 1.5); if any(S(1).mask(:) & ~S(1).maskCore(:)), contour(S(1).maskCore, [0.5 0.5], 'c', 'LineWidth', 1); end
title(sprintf('PB: mean registered z-sum, planes %s + mask', mat2str(S(1).planes)), 'FontSize', 9)
subplot(2, 3, 4); cI = S(1).clusterIdx; imagesc(cI); axis image; colormap(gca, [0 0 0; hsv(S(1).nClusters)]); hold on
for c = 1:S(1).nClusters, [cy, cx] = find(cI == c); if ~isempty(cx), text(mean(cx), mean(cy), num2str(c), 'Color', 'w', 'FontSize', 7, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end, end
title('PB segmentation: 32 glomeruli along the skeleton (16 per hemisphere)', 'FontSize', 9)
subplot(2, 3, 2); imagesc(S(2).projMean); axis image; colormap(gca, bone); hold on; contour(S(2).mask, [0.5 0.5], 'r', 'LineWidth', 1.5)
title(sprintf('EB: x-y projection of the sub-volume (rows %d-%d, cols %d-%d, planes %s) + mask footprint', S(2).crop(1), S(2).crop(2), S(2).crop(3), S(2).crop(4), mat2str(S(2).planes)), 'FontSize', 8)
subplot(2, 3, 3); xz = squeeze(sum(S(2).meanVol, 1))'; imagesc(1:size(xz, 2), zs(S(2).planes), xz); hold on
set(gca, 'YDir', 'reverse'); daspect([1 zStep / umPerPx 1]); colormap(gca, hot)
contour(1:size(xz, 2), zs(S(2).planes), double(ebInfo.footprint), [0.5 0.5], 'c', 'LineWidth', 1)
plot(ebInfo.center(1), interp1(1:numel(S(2).planes), zs(S(2).planes), ebInfo.center(2)), 'c+', 'MarkerSize', 12, 'LineWidth', 2)
xlabel('column (px, cropped)'); ylabel('z (um)'); title('EB: x-z projection (sum over rows) -- the ring is edge-on; cyan = 3-D mask footprint + centre', 'FontSize', 8)
subplot(2, 3, 5); lab = ebInfo.xzLabel; imagesc(1:size(lab, 2), zs(S(2).planes), lab); set(gca, 'YDir', 'reverse'); daspect([1 zStep / umPerPx 1])
colormap(gca, [0 0 0; hsv(S(2).nClusters)]); hold on
for c = 1:S(2).nClusters, [pz, cx] = find(lab == c); if ~isempty(cx), text(mean(cx), interp1(1:numel(S(2).planes), zs(S(2).planes), mean(pz)), num2str(c), 'Color', 'w', 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end, end
xlabel('column (px, cropped)'); ylabel('z (um)'); title('EB segmentation: 16 wedges by angle in the x-z plane (alpha = -pi at #1)', 'FontSize', 9)
subplot(2, 3, 6); bar(1:S(2).nClusters, ebInfo.nVox); xlabel('EB wedge'); ylabel('voxels'); title('voxels per wedge', 'FontSize', 9)
sgtitle([trialName ': masks and segmentation'], 'Interpreter', 'none')
obs2_export(figure(1), fullfile(trialDir, 'mask_qc_PB_EB.png'));

%% 4. sync; per-structure background subtraction; bad volumes
h5list = dir(fullfile(trialDir, '*.h5')); assert(numel(h5list) == 1, 'expected one h5 in %s', trialDir);
sync = pb_h5_sync(fullfile(h5list(1).folder, h5list(1).name));
fprintf('LED source: %s\n', sync.ledSource);
R = load(regFiles{1}, 'regShifts'); sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), R.regShifts(:), 'UniformOutput', false)); clear R
badVol = any(abs(sh) > max_ok_shift_px, 2)';
V = load(regFiles{1}, 'imgData_reg'); img = single(V.imgData_reg); clear V
V = load(regFiles{2}, 'imgData_reg'); img4 = single(V.imgData_reg); clear V
nVol = min([size(img, 3), size(img4, 4), numel(sync.t_volume)]);
img = img(:, :, 1:nVol); img4 = img4(:, :, :, 1:nVol); badVol = badVol(1:nVol);
t_vol = sync.t_volume(1:nVol)';
fprintf('%d volumes; %d with |PB shift| > %d px excluded (applies to both structures)\n', nVol, nnz(badVol), max_ok_shift_px);
% PB: 2-D pixels x time; EB: 3-D voxels x time. Background = mean outside the structure's own mask, per volume.
S(1).vox = reshape(img, [], nVol); bg = mean(S(1).vox(~S(1).mask(:), :), 1); S(1).vox = S(1).vox - bg; S(1).img = reshape(S(1).vox, size(img)); clear img
S(2).vox = reshape(img4, [], nVol); bg = mean(S(2).vox(~S(2).mask3(:), :), 1); S(2).vox = S(2).vox - bg; S(2).dims = size(img4); clear img4
for k = 1:2, S(k).badVol = badVol; end
led_vol = interp1(sync.t_fine, double(sync.led), t_vol, 'nearest', 'extrap') > 0.5;
bar_rad_vol = mod(interp1(sync.t_fine, sync.bar_rad, t_vol, 'nearest', 'extrap'), 2*pi);
bar_deg_vol = rad2deg(bar_rad_vol);
holdSegs = obs2_flat_segments(sync.t_fine, sync.bar_rad, sync.vr_start, sync.fs, hold_flat_tol_V, min_hold_sec);
holdSegs = holdSegs(holdSegs(:,1) < sync.vr_end - 0.5 & holdSegs(:,1) < t_vol(end), :);
inHold_vol = false(1, nVol); pulseInHold = false(size(sync.led_on));
for k = 1:size(holdSegs, 1)
    inHold_vol = inHold_vol | (t_vol >= holdSegs(k,1) & t_vol <= holdSegs(k,2));
    pulseInHold = pulseInHold | (sync.led_on >= holdSegs(k,1) & sync.led_off <= holdSegs(k,2));
end
fprintf('%d held-bar segment(s) >= %g s; LED: %d pulses of %.1f s (%.1f-%.1f s DAQ), %d fully inside holds\n', size(holdSegs, 1), min_hold_sec, ...
    numel(sync.led_on), median(sync.led_off - sync.led_on), min([sync.led_on NaN]), max([sync.led_off NaN]), nnz(pulseInHold));
zoomWin = [];
if ~isempty(sync.led_on)
    zoomWin = [max(sync.led_on(1) - 30, t_vol(1)), min(sync.led_off(end) + 30, t_vol(end))];
    if diff(zoomWin) > 0.9 * (t_vol(end) - t_vol(1)), zoomWin = []; end
end
rwb = obs2_redwhiteblue(256);
bar_deg_plot = bar_deg_vol; bar_deg_plot(abs(diff([bar_deg_plot bar_deg_plot(end)])) > 180) = NaN;

%% 5. z-scored activity per glomerulus / wedge -> heatmaps
for k = 1:2
    nC = S(k).nClusters; cI = S(k).clusterIdx(:); core = S(k).maskCore(:);
    f = nan(nC, nVol); nPix = zeros(nC, 1);
    for c = 1:nC
        pix = cI == c & core; nPix(c) = nnz(pix);
        if nPix(c) > 0, f(c, :) = mean(S(k).vox(pix, :), 1); end
    end
    empty = nPix < max(5, 0.2 * median(nPix(nPix > 0)));
    if any(empty), f(empty, :) = NaN; fprintf('%s: %d/%d units with (almost) no pixels -> NaN: %s\n', S(k).name, nnz(empty), nC, mat2str(find(empty)')); end
    f(:, S(k).badVol) = NaN;
    S(k).f_cluster = f; S(k).nCorePx = nPix;
    S(k).z_cluster = (f - mean(f, 2, 'omitnan')) ./ std(f, 0, 2, 'omitnan');
    for zoomPass = 1:(1 + ~isempty(zoomWin))
        if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
        figure(2); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
        ax1 = subplot(6, 1, 1:3);
        imagesc(t_vol, 1:nC, S(k).z_cluster, 'AlphaData', double(~isnan(S(k).z_cluster))); colormap(ax1, rwb); clim([-3 3]); hold on
        set(ax1, 'Color', [0.7 0.7 0.7]); obs2_mark_holds(holdSegs)
        if strcmp(S(k).name, 'PB'), ylabel('PB glomerulus (1-16 left, 17-32 right)'); else, ylabel('EB wedge (angle around the ring, x-z plane)'); end
        title(sprintf('%s -- %s: z-scored activity (background-subtracted, planes %s, %d units); dashed = held-bar segments (%d)', ...
            trialName, S(k).name, mat2str(S(k).planes), nC, size(holdSegs, 1)), 'Interpreter', 'none')
        cb = colorbar(ax1); cb.Label.String = 'z-score';
        ax2 = subplot(6, 1, 4:5);
        obs2_shade_pulses(sync.led_on, sync.led_off, [0 360]); hold on
        plot(t_vol, bar_deg_plot, 'k', 'LineWidth', 1); ylim([0 360]); yticks(0:90:360); ylabel('bar position (deg)')
        title('bar position on the LED display (red shading = stim LED on)')
        ax3 = subplot(6, 1, 6);
        plot(t_vol, double(led_vol), 'r', 'LineWidth', 1); ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'off', 'on'}); ylabel('stim LED'); xlabel('DAQ time (s)')
        linkaxes([ax1 ax2 ax3], 'x'); xlim(xl); ax1.Position([1 3]) = ax2.Position([1 3]);
        obs2_export(figure(2), fullfile(trialDir, sprintf('glomHeatmap_zscore_%s_%s.png', S(k).name, tag)));
    end
end

%% 6. PVA per structure -> bump vs bar (sign + lag fitted during the sweeps)
bar_pm0 = angle(exp(1i * bar_rad_vol));
for k = 1:2
    [xt, yt] = pol2cart(S(k).alpha, S(k).z_cluster');
    [mu, rho] = cart2pol(mean(xt, 2, 'omitnan'), mean(yt, 2, 'omitnan')); mu = mu(:)'; rho = rho(:)';
    sweepOK = ~inHold_vol & rho >= rho_thresh & t_vol > sync.vr_start;
    cs_same = obs2_circ_std(mu(sweepOK) - bar_pm0(sweepOK)); cs_flip = obs2_circ_std(mu(sweepOK) + bar_pm0(sweepOK));
    if cs_flip < cs_same, bar_sign = -1; signNote = sprintf('%s order is MIRRORED relative to the display (bump ~ -bar)', S(k).name);
    else, bar_sign = 1; signNote = sprintf('%s order runs WITH the display (bump ~ +bar)', S(k).name); end
    fprintf('%s bump-vs-bar: circ std of sweep offset %.2f rad (same) vs %.2f rad (flipped) -> %s\n', S(k).name, cs_same, cs_flip, signNote);
    bar_c = exp(1i * bar_rad_vol); dt_vol = median(diff(t_vol)); lagScan = lag_range_sec(1):dt_vol:lag_range_sec(2); lagScore = nan(size(lagScan));
    for L = 1:numel(lagScan)
        barL = angle(interp1(t_vol, bar_c, t_vol - lagScan(L), 'linear', 'extrap'));
        lagScore(L) = obs2_circ_std(mu(sweepOK) - bar_sign * barL(sweepOK));
    end
    [~, iBest] = min(lagScore); lag_sec = lagScan(iBest); if isnan(lag_sec) || lag_sec < 0, lag_sec = 0; end
    bar_rad_lag = angle(interp1(t_vol, bar_c, t_vol - lag_sec, 'linear', 'extrap')); bar_pm = angle(exp(1i * bar_rad_lag));
    offset = angle(exp(1i * (mu - bar_sign * bar_pm)));
    iUsed = find(abs(lagScan - lag_sec) == min(abs(lagScan - lag_sec)), 1); iZero = find(abs(lagScan) == min(abs(lagScan)), 1);
    fprintf('  lag %.2f s (sweep offset circ std %.1f deg; %.1f at lag 0); mean sweep offset %.1f deg; hold LED off %.1f deg (%.1f), LED on %.1f deg (%.1f); mean rho sweep %.2f, hold off %.2f, hold on %.2f\n', ...
        lag_sec, rad2deg(lagScore(iUsed)), rad2deg(lagScore(iZero)), rad2deg(obs2_circ_mean(offset(sweepOK))), ...
        rad2deg(obs2_circ_mean(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), rad2deg(obs2_circ_std(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), ...
        rad2deg(obs2_circ_mean(offset(inHold_vol & led_vol & rho >= rho_thresh))), rad2deg(obs2_circ_std(offset(inHold_vol & led_vol & rho >= rho_thresh))), ...
        mean(rho(sweepOK), 'omitnan'), mean(rho(inHold_vol & ~led_vol), 'omitnan'), mean(rho(inHold_vol & led_vol), 'omitnan'));
    S(k).mu = mu; S(k).rho = rho; S(k).bar_sign = bar_sign; S(k).lag_sec = lag_sec; S(k).offset = offset; S(k).signNote = signNote; S(k).bar_pm = bar_pm;
    bump = struct('t_vol', t_vol, 'mu', mu, 'rho', rho, 'bar_rad', bar_rad_vol, 'bar_sign', bar_sign, 'bar_rad_lag', bar_rad_lag, 'lag_sec', lag_sec, ...
        'offset', offset, 'led', led_vol, 'inHold', inHold_vol, 'badVol', S(k).badVol, 'holdSegs', holdSegs, 'led_on', sync.led_on, 'led_off', sync.led_off, ...
        'vr_start', sync.vr_start, 'vr_end', sync.vr_end, 'rho_thresh', rho_thresh, 'z_cluster', S(k).z_cluster, 'f_cluster', S(k).f_cluster, ...
        'clusterIdx', S(k).clusterIdx, 'nCorePx', S(k).nCorePx, 'alpha', S(k).alpha, 'trialName', trialName, 'signNote', signNote, ...
        'structure', S(k).name, 'planes', S(k).planes); %#ok<NASGU>
    save(fullfile(trialDir, sprintf('bump_results_%s.mat', S(k).name)), 'bump');
    mu_plot = rad2deg(mu); mu_plot(rho < rho_thresh) = NaN; mu_plot(abs(diff([mu_plot mu_plot(end)])) > 180) = NaN;
    barS_plot = rad2deg(angle(exp(1i * bar_sign * bar_pm))); barS_plot(abs(diff([barS_plot barS_plot(end)])) > 180) = NaN;
    off_plot = rad2deg(offset); off_plot(rho < rho_thresh) = NaN;
    for zoomPass = 1:(1 + ~isempty(zoomWin))
        if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
        figure(3); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
        ax1 = subplot(3, 1, 1); obs2_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
        plot(t_vol, barS_plot, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); plot(t_vol, mu_plot, 'b.', 'MarkerSize', 4); obs2_mark_holds(holdSegs)
        ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg)')
        legend({sprintf('bar (signed, lagged %.2f s)', lag_sec), sprintf('%s bump (PVA), rho >= %.1f', S(k).name, rho_thresh)}, 'Location', 'northeastoutside')
        title(sprintf('%s -- %s bump vs bar: %s; lag %.2f s; red = stim LED on', trialName, S(k).name, signNote, lag_sec), 'Interpreter', 'none')
        ax2 = subplot(3, 1, 2); obs2_shade_pulses(sync.led_on, sync.led_off, [0 max(rho) * 1.05]); hold on
        plot(t_vol, rho, 'k'); yline(rho_thresh, 'r:'); ylabel(sprintf('%s bump strength (rho)', S(k).name))
        ax3 = subplot(3, 1, 3); obs2_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
        plot(t_vol, off_plot, 'b.', 'MarkerSize', 4); yline(0, 'k:'); obs2_mark_holds(holdSegs)
        ylim([-180 180]); yticks(-180:90:180); ylabel(sprintf('bump - bar(t - %.2f s) (deg)', lag_sec)); xlabel('DAQ time (s)')
        linkaxes([ax1 ax2 ax3], 'x'); xlim(xl); for ax = [ax2 ax3], ax.Position([1 3]) = ax1.Position([1 3]); end
        obs2_export(figure(3), fullfile(trialDir, sprintf('bumpVsBar_%s_%s.png', S(k).name, tag)));
    end
end

% PB vs EB directly: both bumps in the display frame (each with its own fitted sign), and their difference
figure(5); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
muPB = rad2deg(angle(exp(1i * S(1).bar_sign * S(1).mu))); muEB = rad2deg(angle(exp(1i * S(2).bar_sign * S(2).mu)));
muPB(S(1).rho < rho_thresh) = NaN; muEB(S(2).rho < rho_thresh) = NaN;
muPB(abs(diff([muPB muPB(end)])) > 180) = NaN; muEB(abs(diff([muEB muEB(end)])) > 180) = NaN;
dPE = rad2deg(angle(exp(1i * (S(1).bar_sign * S(1).mu - S(2).bar_sign * S(2).mu)))); dPE(S(1).rho < rho_thresh | S(2).rho < rho_thresh) = NaN;
barD = rad2deg(angle(exp(1i * bar_rad_vol))); barD(abs(diff([barD barD(end)])) > 180) = NaN;
ax1 = subplot(3, 1, 1); obs2_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
plot(t_vol, barD, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); plot(t_vol, muPB, 'b.', 'MarkerSize', 4); plot(t_vol, muEB, 'r.', 'MarkerSize', 4); obs2_mark_holds(holdSegs)
ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg, display frame)'); legend({'bar', 'PB bump', 'EB bump'}, 'Location', 'northeastoutside')
title(sprintf('%s: PB and EB bumps in display coordinates (each with its own fitted sign); red shading = stim LED on', trialName), 'Interpreter', 'none')
ax2 = subplot(3, 1, 2); obs2_shade_pulses(sync.led_on, sync.led_off, [0 1.05 * max([S(1).rho S(2).rho])]); hold on
plot(t_vol, S(1).rho, 'b'); plot(t_vol, S(2).rho, 'r'); ylabel('bump strength (rho)'); legend({'PB', 'EB'}, 'Location', 'northeastoutside')
ax3 = subplot(3, 1, 3); obs2_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
plot(t_vol, dPE, 'k.', 'MarkerSize', 4); yline(0, 'k:'); obs2_mark_holds(holdSegs)
ylim([-180 180]); yticks(-180:90:180); ylabel('PB bump - EB bump (deg)'); xlabel('DAQ time (s)')
okB = ~isnan(dPE); fprintf('PB - EB bump offset: circ mean %.1f deg, circ std %.1f deg over %d volumes with both rho >= %.1f\n', ...
    rad2deg(obs2_circ_mean(deg2rad(dPE(okB)))), rad2deg(obs2_circ_std(deg2rad(dPE(okB)))), nnz(okB), rho_thresh);
linkaxes([ax1 ax2 ax3], 'x'); xlim([t_vol(1) t_vol(end)]); for ax = [ax2 ax3], ax.Position([1 3]) = ax1.Position([1 3]); end
obs2_export(figure(5), fullfile(trialDir, 'bumpPBvsEB_fullTrial.png'));

%% 7. stim on - off diff images (pulses inside holds): PB x-y; EB x-y and x-z projections
if ~isempty(sync.led_on)
    if nnz(pulseInHold) >= 3, pulseIdx = find(pulseInHold); offMask = inHold_vol; note = sprintf('%d pulses inside holds', numel(pulseIdx));
    else, pulseIdx = 1:numel(sync.led_on); offMask = true(1, nVol); note = sprintf('all %d pulses', numel(pulseIdx)); end
    onIdx = false(1, nVol); for k = pulseIdx, onIdx = onIdx | (t_vol >= sync.led_on(k) & t_vol < sync.led_off(k)); end
    offIdx = offMask & ~led_vol & t_vol >= sync.led_on(pulseIdx(1)) - 5 & t_vol <= sync.led_off(pulseIdx(end)) + 10;
    on = onIdx & ~badVol; off = offIdx & ~badVol;
    figure(4); clf; set(gcf, 'Position', [50 50 1800 800], 'Color', 'w')
    dPB = reshape(mean(S(1).vox(:, on), 2) - mean(S(1).vox(:, off), 2), size(S(1).mask));
    subplot(2, 3, 1); imagesc(S(1).projMean); axis image; colormap(gca, bone); hold on; contour(S(1).mask, [0.5 0.5], 'r', 'LineWidth', 1); colorbar; title('PB mean image + mask')
    subplot(2, 3, 4); dl = max(prctile(abs(dPB(S(1).mask)), 99), eps); imagesc(dPB, [-dl dl]); axis image; colormap(gca, rwb); hold on; contour(S(1).mask, [0.5 0.5], 'k', 'LineWidth', 0.75); colorbar
    title(sprintf('PB: mean(LED on) - mean(LED off), %s; %d on / %d off', note, nnz(on), nnz(off)), 'FontSize', 9)
    dEB = reshape(mean(S(2).vox(:, on), 2) - mean(S(2).vox(:, off), 2), S(2).dims(1:3));
    dEB(~S(2).mask3) = 0;
    subplot(2, 3, 2); imagesc(S(2).projMean); axis image; colormap(gca, bone); hold on; contour(S(2).mask, [0.5 0.5], 'r', 'LineWidth', 1); colorbar; title('EB x-y projection + footprint')
    dxy = sum(dEB, 3); subplot(2, 3, 5); dl = max(prctile(abs(dxy(S(2).mask)), 99), eps); imagesc(dxy, [-dl dl]); axis image; colormap(gca, rwb); hold on; contour(S(2).mask, [0.5 0.5], 'k', 'LineWidth', 0.75); colorbar
    title('EB: on - off, x-y projection (sum over planes, mask only)', 'FontSize', 9)
    dxz = squeeze(sum(dEB, 1))'; subplot(2, 3, 3); xzM = squeeze(sum(S(2).meanVol, 1))'; imagesc(1:size(xzM, 2), zs(S(2).planes), xzM); set(gca, 'YDir', 'reverse'); daspect([1 zStep / umPerPx 1]); colormap(gca, hot); colorbar
    hold on; contour(1:size(xzM, 2), zs(S(2).planes), double(ebInfo.footprint), [0.5 0.5], 'c'); title('EB x-z projection (mean) + footprint'); ylabel('z (um)')
    subplot(2, 3, 6); dl = max(prctile(abs(dxz(ebInfo.footprint)), 99), eps); imagesc(1:size(dxz, 2), zs(S(2).planes), dxz, [-dl dl]); set(gca, 'YDir', 'reverse'); daspect([1 zStep / umPerPx 1]); colormap(gca, rwb); colorbar
    hold on; contour(1:size(dxz, 2), zs(S(2).planes), double(ebInfo.footprint), [0.5 0.5], 'k'); title('EB: on - off, x-z projection (sum over rows, mask only)', 'FontSize', 9); ylabel('z (um)')
    sgtitle([trialName ': stim on - off diff images (red = brighter with LED on)'], 'Interpreter', 'none')
    obs2_export(figure(4), fullfile(trialDir, 'holdDiff_stimOn_minus_off_PB_EB.png'));
end
fprintf('done.\n');

%% local functions
function segs = obs2_flat_segments(t, bar, tStart, fs, tolV, minSec)
barF = movmedian(bar(:), round(0.25 * fs)); win = round(2 * fs);
flat = (movmax(barF, win) - movmin(barF, win)) < tolV & t(:) > tStart;
d = diff([0; flat; 0]); s = find(d == 1); e = find(d == -1) - 1; keep = (e - s + 1) / fs >= minSec;
segs = [t(s(keep)), t(e(keep))]; if isempty(segs), segs = zeros(0, 2); end
end
function obs2_export(fig, outFile)
drawnow; set(findall(fig, 'Type', 'axestoolbar'), 'Visible', 'off');
for ax = findall(fig, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow; exportgraphics(fig, outFile, 'Resolution', 150); fprintf('saved %s\n', outFile);
[p, n] = fileparts(outFile); vis = get(fig, 'Visible'); set(fig, 'Visible', 'on'); savefig(fig, fullfile(p, [n '.fig'])); set(fig, 'Visible', vis);
end
function obs2_mark_holds(segs)
for k = 1:size(segs, 1), xline(segs(k,1), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); xline(segs(k,2), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); end
end
function obs2_shade_pulses(onsets, offsets, yl)
for k = 1:numel(onsets), patch([onsets(k) offsets(k) offsets(k) onsets(k)], [yl(1) yl(1) yl(2) yl(2)], [1 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off'); end
ylim(yl)
end
function m = obs2_circ_mean(x), m = angle(mean(exp(1i * x(:)), 'omitnan')); end
function s = obs2_circ_std(x)
x = x(~isnan(x)); if isempty(x), s = NaN; return; end
R = abs(mean(exp(1i * x(:)))); s = sqrt(-2 * log(max(R, eps)));
end
function cmap = obs2_redwhiteblue(n)
half = n / 2; cmap = [[linspace(0,1,half)', linspace(0,1,half)', ones(half,1)]; [ones(half,1), linspace(1,0,half)', linspace(1,0,half)']];
end
function clusterIdx = obs_pb_glomeruli(mask, nPerHemisphere)
% skeleton -> graph_sort -> 2N evenly spaced centroids -> nearest centroid (process_im convention)
[y_mask, x_mask] = find(mask); min_axis = min(range(x_mask), range(y_mask));
skel = bwskel(bwmorph(mask, 'majority')); epIdx = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx), error('obs_pb_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch.'); end
D = bwdistgeodesic(skel, epIdx); [~, idxFar] = max(D(:)); D2 = bwdistgeodesic(skel, idxFar); mid = (D + D2) == mode(D + D2, 'all');
[y_mid, x_mid] = find(mid); [x_mid, y_mid] = graph_sort(x_mid, y_mid);
xq = -min_axis:(length(x_mid) + min_axis);
x_mid = round(interp1(1:length(x_mid), x_mid, xq, 'linear', 'extrap')); y_mid = round(interp1(1:length(y_mid), y_mid, xq, 'linear', 'extrap'));
keep = ismember([x_mid', y_mid'], [x_mask, y_mask], 'rows'); x_mid = x_mid(keep); y_mid = y_mid(keep);
nClusters = 2 * nPerHemisphere; xq2 = linspace(1, length(y_mid), 2*nClusters + 1)';
centroids = [interp1(1:length(y_mid), y_mid, xq2), interp1(1:length(x_mid), x_mid, xq2)]; centroids = centroids(2:2:end-1, :);
[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);
clusterIdx = zeros(size(mask)); clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end
