%% overshoot_eb_led_script
% EB-only version of overshoot_pb_led_script.m: bar + stim-LED trials (same
% rig configuration, same WaveSurfer / ViRMEn logging) where the volume is
% centred on the ellipsoid body instead of the protocerebral bridge.
%
% Same analysis as the PB script, but the EB is segmented in 3-D:
%   * the planes are kept SEPARATE (summing them collapses the ring when it
%     is tilted or edge-on -- see overshoot_pb_eb_script.m), cropped to the
%     EB. The crop is found automatically from the first oel_crop_probe_vols
%     volumes (largest bright 3-D blob + margin) unless oel_crop is given.
%   * motion is estimated on the full-FOV z-sum (more texture than the small
%     EB crop, which can make rigid registration run away) with pb_register,
%     and the shifts are applied to every EB plane.
%   * 3-D mask on the mean registered sub-volume: anisotropic Gaussian
%     smoothing, Otsu threshold (or a percentile, mask_thresh_pct), largest
%     26-connected component, 1-voxel closing. The ring hole is NOT filled.
%   * claude\pb_eb_ring3d fits the ring plane by PCA of the mask voxels in
%     um (works face-on, edge-on or tilted; force it with ring_plane) and
%     splits the ring into N wedges by angle around the ring centre.
% Then exactly the PB recipe on the N wedges (alpha once around the ring):
% per-volume background subtraction (mean of the crop outside the mask),
% z-score, PVA, sign and lag vs the bar fitted during the sweeps.
%
% Outputs in the trial folder (cached stages skip if the file exists; EB
% suffix so they never collide with PB outputs):
%   imgData_EB.mat (cropped per-plane movie + z-sum for registration),
%   imgData_EB_reg.mat, mask_EB.mat, planeMeans_EB.png, mask_qc_EB.png
%   wedgeHeatmap_zscore_EB_{fullTrial,stim}.png, bumpVsBar_EB_{fullTrial,stim}.png
%   holdDiff_stimOn_minus_off_EB.png, bump_results_EB.mat (same fields as
%   the PB bump_results.mat; z_cluster/f_cluster are wedges; plus bump.ring)
%
% Check mask_qc_EB.png before trusting anything downstream: the ring-plane
% panel should show a closed ring (or disc) with wedges 1..N running once
% around it, and the planarity in the printout should be >~0.4. If the
% plane fit picked the wrong orientation, set ring_plane = 'xy' / 'xz' /
% 'yz' (or oel_ring_plane) and rerun section 3 (cheap).
%
% Headless:
%   matlab -batch "oel_trialDir='Z:\...\trial_002\'; run('overshoot_eb_led_script.m')"
% Overrides: oel_overwrite_mask / oel_overwrite_reg / oel_overwrite_video,
% oel_min_hold_sec, oel_crop ([row1 row2 col1 col2], [] = auto),
% oel_planes (plane indices within the volume, default all valid planes),
% oel_ring_plane ('auto' | 'xy' | 'xz' | 'yz').

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'claude'));
if isempty(which('circ_mean')), addpath(fullfile(thisDir, '..', 'circ_stats')); end

trialDir = '';                       % '' -> pick the trial folder in a dialog
if exist('oel_trialDir', 'var') && ~isempty(oel_trialDir), trialDir = oel_trialDir; end
if isempty(trialDir), trialDir = uigetdir('Z:\noah_np123\Data\flyg\overshoot\', 'EB trial folder'); assert(ischar(trialDir), 'no trial folder chosen'); end
if trialDir(end) ~= '\' && trialDir(end) ~= '/', trialDir(end+1) = filesep; end
parts = strsplit(trialDir, filesep); parts = parts(~cellfun(@isempty, parts));
trialName = [parts{end-1} filesep parts{end}];

n_wedges = 16;                       % matches the PB's 16 glomeruli per hemisphere
crop = [];                           % [row1 row2 col1 col2] of the EB sub-volume; [] = automatic
crop_margin_um = 8; crop_probe_vols = 300;
planes = [];                         % plane indices within the volume; [] = all valid (non-flyback) planes
mask_smooth_sigma_um = 1.2; mask_smooth_sigma_planes = 0.6;
mask_thresh_pct = [];                % [] = Otsu on the smoothed mean volume; or a percentile (e.g. 92 as in overshoot_pb_eb_script)
ring_plane = 'auto';                 % pb_eb_ring3d 'plane': 'auto' (PCA) | 'xy' | 'xz' | 'yz'
maxShift_px = [6 8]; max_ok_shift_px = 8;
min_hold_sec = 4; hold_flat_tol_V = 0.1;
lag_range_sec = [-2 5]; rho_thresh = 0.1;
pattern_names = containers.Map({2}, {'multigauss'});    % pattern_idx -> label; unknown idx shows as 'pattern N'
overwrite_video = false; overwrite_reg = false; overwrite_mask = false;
if exist('oel_min_hold_sec', 'var') && ~isempty(oel_min_hold_sec), min_hold_sec = oel_min_hold_sec; end
if exist('oel_crop', 'var') && ~isempty(oel_crop), crop = oel_crop; end
if exist('oel_planes', 'var') && ~isempty(oel_planes), planes = oel_planes; end
if exist('oel_ring_plane', 'var') && ~isempty(oel_ring_plane), ring_plane = oel_ring_plane; end
if exist('oel_overwrite_video', 'var') && ~isempty(oel_overwrite_video), overwrite_video = oel_overwrite_video; end
if exist('oel_overwrite_reg', 'var') && ~isempty(oel_overwrite_reg), overwrite_reg = oel_overwrite_reg; end
if exist('oel_overwrite_mask', 'var') && ~isempty(oel_overwrite_mask), overwrite_mask = oel_overwrite_mask; end
if overwrite_video, overwrite_reg = true; end
if overwrite_reg, overwrite_mask = true; end

% scan geometry from the tif header
tifList = dir(fullfile(trialDir, '*.tif')); assert(numel(tifList) == 1, 'expected exactly one .tif in %s', trialDir);
tifPath = fullfile(tifList(1).folder, tifList(1).name);
fid = fopen(tifPath, 'r'); hdr = fread(fid, 80000, 'uint8=>char')'; fclose(fid);
getNum = @(key) str2double(regexp(hdr, ['SI\.' regexptranslate('escape', key) '\s*=\s*([-0-9.eE]+)'], 'tokens', 'once'));
framesPerVol = getNum('hStackManager.numFramesPerVolumeWithFlyback'); nValid = getNum('hStackManager.numFramesPerVolume');
zoom = getNum('hRoiManager.scanZoomFactor'); W = getNum('hRoiManager.pixelsPerLine'); H = getNum('hRoiManager.linesPerFrame');
assert(all(~isnan([framesPerVol nValid zoom W H])), 'could not parse the ScanImage header of %s', tifPath);
umPerPx = 0.3164 * (4 / zoom) * (256 / W);               % 0.3164 um/px at zoom 4, 256 px
zsTok = regexp(hdr, 'SI\.hStackManager\.zs\s*=\s*\[([^\]]*)\]', 'tokens', 'once');
if ~isempty(zsTok), zs = str2num(zsTok{1}); else, zs = []; end %#ok<ST2NM>
if numel(zs) < nValid
    zStepHdr = getNum('hStackManager.stackZStepSize'); if isnan(zStepHdr), zStepHdr = 4; end
    warning('no SI.hStackManager.zs list in the header; assuming %g um between planes.', zStepHdr); zs = (0:nValid-1) * zStepHdr;
end
if isempty(planes), planes = 1:nValid; end                % flyback frames are the last ones in the volume
zPlanes = zs(planes); zStep = median(diff(zPlanes));
fprintf('[%s] %dx%d px, %d frames/volume (%d valid), zoom %g -> %.4f um/px; EB planes %s (z %g-%g um, %g um steps)\n', ...
    trialName, H, W, framesPerVol, nValid, zoom, umPerPx, mat2str(planes), zPlanes(1), zPlanes(end), zStep);

%% 1. raw tif -> EB per-plane sub-volume movie (auto crop) + full-FOV z-sum for registration (cached)
vidFile = fullfile(trialDir, 'imgData_EB.mat');
if overwrite_video || ~isfile(vidFile)
    if isempty(crop)
        fprintf('[%s] probing the first %d volumes for the EB crop...\n', trialName, crop_probe_vols);
        o = pb_read_raw_tif_planes(tifPath, framesPerVol, {}, 'maxVolumes', crop_probe_vols, 'planeMeans', true, 'verbose', false);
        crop = oel_auto_crop(o.planeMean(:, :, planes), umPerPx, mask_smooth_sigma_um, crop_margin_um);
        fprintf('  auto crop rows %d-%d, cols %d-%d (of %d x %d)\n', crop, H, W); clear o
    end
    fprintf('[%s] reading raw tif (EB planes, separate + full-FOV z-sum)...\n', trialName);
    o = pb_read_raw_tif_planes(tifPath, framesPerVol, {planes, planes}, 'separate', [true false], 'crop', {crop, []}, 'planeMeans', true);
    imgData = o.stacks{1}; imgData_sum = o.stacks{2}; nVolumes = o.nVolumes; planeMean = o.planeMean; %#ok<NASGU>
    save(vidFile, 'imgData', 'imgData_sum', 'planes', 'zPlanes', 'crop', 'nVolumes', 'planeMean', '-v7.3');
    fprintf('  saved %s %s\n', vidFile, mat2str(size(imgData)));
    clear o imgData imgData_sum
else
    fprintf('skip (exists): %s\n', vidFile);
end
V = load(vidFile, 'crop', 'planeMean', 'planes'); crop = V.crop; planeMean = V.planeMean;
assert(isequal(V.planes(:), planes(:)), 'imgData_EB.mat was read with planes %s, now %s -- set overwrite_video', mat2str(V.planes), mat2str(planes)); clear V
figure(9); clf; set(gcf, 'Position', [50 50 1600 800], 'Color', 'w')
nCol = ceil(sqrt(numel(planes) * 1.6)); nRow = ceil(numel(planes) / nCol);
for k = 1:numel(planes)
    subplot(nRow, nCol, k); imagesc(planeMean(:, :, planes(k))); axis image off; colormap gray; hold on
    rectangle('Position', [crop(3) - 0.5, crop(1) - 0.5, crop(4) - crop(3) + 1, crop(2) - crop(1) + 1], 'EdgeColor', 'y', 'LineWidth', 1)
    title(sprintf('plane %d, z = %g um', planes(k), zs(planes(k))), 'FontSize', 8)
end
sgtitle(sprintf('%s: per-plane mean image; yellow = EB crop', trialName), 'Interpreter', 'none')
oel_export(figure(9), fullfile(trialDir, 'planeMeans_EB.png'));

%% 2. motion correction: shifts from the full-FOV z-sum, applied to every EB plane (cached)
regFile = fullfile(trialDir, 'imgData_EB_reg.mat');
if overwrite_reg || ~isfile(regFile)
    V = load(vidFile, 'imgData_sum');
    fprintf('[%s] motion-correcting the z-sum (%d volumes)...\n', trialName, size(V.imgData_sum, 3));
    [~, regShifts] = pb_register(V.imgData_sum, 'smoothWindow', 5, 'maxShift', maxShift_px, 'verbose', false, 'figureVisible', false); clear V
    sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), regShifts(:), 'UniformOutput', false));
    fprintf('  shift range row = [%.1f %.1f] px, col = [%.1f %.1f] px\n', min(sh(:,1)), max(sh(:,1)), min(sh(:,2)), max(sh(:,2)));
    V = load(vidFile, 'imgData'); [nR, nC, nP, nV] = size(V.imgData);
    fprintf('[%s] applying the shifts to the %d EB planes...\n', trialName, nP);
    imgData_reg = zeros(nR, nC, nP, nV, 'single');
    for pIdx = 1:nP
        imgData_reg(:, :, pIdx, :) = reshape(pb_register(squeeze(V.imgData(:, :, pIdx, :)), 'applyShifts', regShifts, 'smoothWindow', 0, ...
            'maxShift', maxShift_px, 'verbose', false, 'figureVisible', false), nR, nC, 1, nV);
    end
    save(regFile, 'imgData_reg', 'regShifts', '-v7.3'); fprintf('  saved %s\n', regFile); clear V imgData_reg regShifts
else
    fprintf('skip (exists): %s\n', regFile);
end

%% 3. 3-D mask + ring wedges (cached mask), QC figure
maskFile = fullfile(trialDir, 'mask_EB.mat');
V = load(regFile, 'imgData_reg'); meanVol = mean(V.imgData_reg, 4); clear V   % rows x cols x planes
if overwrite_mask || ~isfile(maskFile)
    sig = max(mask_smooth_sigma_um / umPerPx, 0.5);
    vs = imgaussfilt3(meanVol, [sig sig mask_smooth_sigma_planes]);
    vc = min(max(vs, prctile(vs, 5, 'all')), prctile(vs, 99.5, 'all'));
    if isempty(mask_thresh_pct), thr = multithresh(vc(:), 1); thrNote = 'Otsu';
    else, thr = prctile(vc(:), mask_thresh_pct); thrNote = sprintf('%g th percentile', mask_thresh_pct); end
    mask3 = vc > thr;
    cc = bwconncomp(mask3, 26); nComp = cc.NumObjects; [~, iBig] = max(cellfun(@numel, cc.PixelIdxList));
    mask3 = false(size(mask3)); mask3(cc.PixelIdxList{iBig}) = true;
    mask3 = imclose(mask3, strel('sphere', 1));
    save(maskFile, 'mask3', 'thr', 'thrNote'); fprintf('EB mask: %s threshold, largest of %d components, %d voxels; saved %s\n', thrNote, nComp, nnz(mask3), maskFile);
else
    fprintf('skip (exists): %s\n', maskFile);
end
M = load(maskFile); mask3 = M.mask3; clear M
[clusterIdx, alpha, ring] = pb_eb_ring3d(mask3, n_wedges, umPerPx, zPlanes, 'plane', ring_plane, 'weights', meanVol);
nClusters = n_wedges;
fprintf('EB ring (%s plane): normal [x y z] = %s, planarity %.2f, centre col %.1f row %.1f (image), z %.1f um; semi-axes %.1f x %.1f um; hole %d; %d-%d voxels per wedge\n', ...
    ring.plane, mat2str(round(ring.axes(:, 3)', 2)), ring.planarity, ring.center_vox(1) + crop(3) - 1, ring.center_vox(2) + crop(1) - 1, ring.center_um(3), ...
    ring.semiAxes_um, ring.hasHole, min(ring.nVox), max(ring.nVox));
mask = any(mask3, 3);                % x-y footprint, for display
cmapW = [0 0 0; hsv(nClusters)];
figure(1); clf; set(gcf, 'Position', [50 50 1800 850], 'Color', 'w')
subplot(2, 4, 1); imagesc(max(meanVol, [], 3)); axis image; colormap(gca, bone); hold on; contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.25)
title(sprintf('x-y max projection (crop rows %d-%d, cols %d-%d) + mask', crop), 'FontSize', 8)
subplot(2, 4, 2); xz = squeeze(max(meanVol, [], 1))'; imagesc(1:size(xz, 2), zPlanes, xz); hold on
contour(1:size(xz, 2), zPlanes, double(squeeze(any(mask3, 1))'), [0.5 0.5], 'r', 'LineWidth', 1.25)
set(gca, 'YDir', 'reverse'); daspect([1 umPerPx 1]); colormap(gca, bone); xlabel('column (px)'); ylabel('z (um)'); title('x-z max projection + mask', 'FontSize', 8)
subplot(2, 4, 3); yz = squeeze(max(meanVol, [], 2))'; imagesc(1:size(yz, 2), zPlanes, yz); hold on
contour(1:size(yz, 2), zPlanes, double(squeeze(any(mask3, 2))'), [0.5 0.5], 'r', 'LineWidth', 1.25)
set(gca, 'YDir', 'reverse'); daspect([1 umPerPx 1]); colormap(gca, bone); xlabel('row (px)'); ylabel('z (um)'); title('y-z max projection + mask', 'FontSize', 8)
subplot(2, 4, 4); [xg, yg, zg] = meshgrid((1:size(mask3, 2)) * umPerPx, (1:size(mask3, 1)) * umPerPx, zPlanes); hold on
for c = 1:nClusters
    m = smooth3(double(clusterIdx == c), 'box', 3); if ~any(m(:) > 0.3), continue, end
    patch(isosurface(xg, yg, zg, m, 0.3), 'FaceColor', cmapW(c + 1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.9);
end
set(gca, 'ZDir', 'reverse', 'YDir', 'reverse'); axis equal tight; view(3); camlight; lighting gouraud; grid on
xlabel('x (um)'); ylabel('y (um)'); zlabel('z (um)'); title('3-D wedges', 'FontSize', 8)
subplot(2, 4, 5); imagesc(ring.uGrid, ring.vGrid, ring.planeImg, 'AlphaData', double(~isnan(ring.planeImg))); axis image; colormap(gca, bone); set(gca, 'Color', 'k')
hold on; plot(0, 0, 'c+', 'MarkerSize', 10, 'LineWidth', 1.5); xlabel('e1 (um)'); ylabel('e2 (um)')
title(sprintf('mean F in the ring plane (%s; planarity %.2f)', ring.plane, ring.planarity), 'FontSize', 8)
subplot(2, 4, 6); imagesc(ring.uGrid, ring.vGrid, ring.planeLabel); axis image; colormap(gca, cmapW); clim([-0.5 nClusters + 0.5]); hold on
for c = 1:nClusters, [vy, ux] = find(ring.planeLabel == c); if ~isempty(ux), text(ring.uGrid(round(mean(ux))), ring.vGrid(round(mean(vy))), num2str(c), 'Color', 'w', 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end, end
xlabel('e1 (um)'); ylabel('e2 (um)'); title(sprintf('%d wedges in the ring plane (alpha = -pi at #1)', nClusters), 'FontSize', 8)
subplot(2, 4, 7); lab = zeros(size(mask)); for c = 1:nClusters, lab(any(clusterIdx == c, 3) & lab == 0) = c; end
imagesc(lab); axis image; colormap(gca, cmapW); clim([-0.5 nClusters + 0.5]); title('wedges, x-y view (top-most label per pixel)', 'FontSize', 8)
subplot(2, 4, 8); bar(1:nClusters, ring.nVox); xlabel('EB wedge'); ylabel('voxels'); title('voxels per wedge', 'FontSize', 8)
sgtitle(sprintf('%s: EB 3-D mask (%d voxels) and ring segmentation', trialName, nnz(mask3)), 'Interpreter', 'none')
oel_export(figure(1), fullfile(trialDir, 'mask_qc_EB.png'));

%% 4. sync (volume times, LED, bar, holds), bad volumes, background subtraction
h5list = dir(fullfile(trialDir, '*.h5')); assert(numel(h5list) == 1, 'expected exactly one .h5 in %s', trialDir);
sync = pb_h5_sync(fullfile(h5list(1).folder, h5list(1).name));
fprintf('LED source: %s\n', sync.ledSource);
if isempty(sync.vr_start)
    % vr_ao0_copy only carries the DOT angle: with the dot off it sits at 0 V
    % and vr_start/vr_end come back empty -> cue angle from the iteration log
    sync = oel_cue_from_iter_log(sync, trialDir);
    fprintf('bar channel flat -> cue angle from the iteration log (%s); VR %.1f-%.1f s DAQ\n', sync.barSource, sync.vr_start, sync.vr_end);
end
S = load(regFile, 'imgData_reg', 'regShifts'); img4 = single(S.imgData_reg);
sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), S.regShifts(:), 'UniformOutput', false)); clear S
badVol = any(abs(sh) > max_ok_shift_px, 2)';
nVol = min(size(img4, 4), numel(sync.t_volume));
if size(img4, 4) ~= numel(sync.t_volume), fprintf('note: %d volumes in the video vs %d volume clocks -> using %d\n', size(img4, 4), numel(sync.t_volume), nVol); end
dims = size(img4); dims = dims(1:3);
vox = reshape(img4(:, :, :, 1:nVol), [], nVol); clear img4
badVol = badVol(1:nVol); t_vol = sync.t_volume(1:nVol)';
fprintf('%d volumes (DAQ %.1f-%.1f s); %d with |shift| > %d px excluded\n', nVol, t_vol(1), t_vol(end), nnz(badVol), max_ok_shift_px);
led_vol = interp1(sync.t_fine, double(sync.led), t_vol, 'nearest', 'extrap') > 0.5;
bar_rad_vol = mod(interp1(sync.t_fine, sync.bar_rad, t_vol, 'nearest', 'extrap'), 2*pi); bar_deg_vol = rad2deg(bar_rad_vol);
holdSegs = oel_flat_segments(sync.t_fine, sync.bar_rad, sync.vr_start, sync.fs, hold_flat_tol_V, min_hold_sec);
holdSegs = holdSegs(holdSegs(:,1) < sync.vr_end - 0.5 & holdSegs(:,1) < t_vol(end), :);
inHold_vol = false(1, nVol); pulseInHold = false(size(sync.led_on));
for k = 1:size(holdSegs, 1)
    inHold_vol = inHold_vol | (t_vol >= holdSegs(k,1) & t_vol <= holdSegs(k,2));
    pulseInHold = pulseInHold | (sync.led_on >= holdSegs(k,1) & sync.led_off <= holdSegs(k,2));
end
fprintf('%d held-bar segment(s) >= %g s\n', size(holdSegs, 1), min_hold_sec);
for k = 1:size(holdSegs, 1)
    segAng = rad2deg(median(sync.bar_rad(sync.t_fine >= holdSegs(k,1) & sync.t_fine <= holdSegs(k,2))));
    iL = find(sync.led_on >= holdSegs(k,1) & sync.led_off <= holdSegs(k,2), 1);
    if isempty(iL), tag = ''; else, tag = sprintf('   LED pulse %d', iL); end
    fprintf('  %7.1f - %7.1f s (%5.1f s) at %6.1f deg%s\n', holdSegs(k,1), holdSegs(k,2), diff(holdSegs(k,:)), segAng, tag);
end
if isempty(sync.led_on), fprintf('LED: no pulses.\n');
else, fprintf('LED: %d pulses of %.1f s (%.1f-%.1f s DAQ), %d fully inside holds\n', numel(sync.led_on), median(sync.led_off - sync.led_on), sync.led_on(1), sync.led_off(end), nnz(pulseInHold)); end
visEp = oel_visual_epochs(trialDir, sync.vr_start, sync.vr_end, pattern_names);
if isempty(visEp), fprintf('visual epochs: no pattern log found (sync_info.mat iter / protocol csv).\n');
else
    fprintf('visual epochs (from %s):\n', visEp(1).src);
    for k = 1:numel(visEp), fprintf('  %7.1f - %7.1f s  %s\n', visEp(k).t0, visEp(k).t1, visEp(k).label); end
end
zoomWin = [];
if ~isempty(sync.led_on)
    zoomWin = [max(sync.led_on(1) - 30, t_vol(1)), min(sync.led_off(end) + 30, t_vol(end))];
    if diff(zoomWin) > 0.9 * (t_vol(end) - t_vol(1)), zoomWin = []; end
end
bg = mean(vox(~mask3(:), :), 1); vox = vox - bg;           % background = mean of the crop outside the 3-D mask, per volume
rwb = oel_redwhiteblue(256);
bar_deg_plot = bar_deg_vol; bar_deg_plot(abs(diff([bar_deg_plot bar_deg_plot(end)])) > 180) = NaN;

%% 5. wedge traces -> z-score -> heatmap
f_cluster = nan(nClusters, nVol); nCorePx = zeros(nClusters, 1);
for c = 1:nClusters
    pix = clusterIdx(:) == c; nCorePx(c) = nnz(pix);
    if nCorePx(c) > 0, f_cluster(c, :) = mean(vox(pix, :), 1); end
end
emptyGlom = nCorePx < max(10, 0.2 * median(nCorePx(nCorePx > 0)));
if any(emptyGlom), f_cluster(emptyGlom, :) = NaN; fprintf('%d/%d wedges with (almost) no voxels -> NaN: %s\n', nnz(emptyGlom), nClusters, mat2str(find(emptyGlom)')); end
f_cluster(:, badVol) = NaN;
z_cluster = (f_cluster - mean(f_cluster, 2, 'omitnan')) ./ std(f_cluster, 0, 2, 'omitnan');
for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
    figure(2); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(7, 1, 1:3);
    imagesc(t_vol, 1:nClusters, z_cluster, 'AlphaData', double(~isnan(z_cluster))); colormap(ax1, rwb); clim([-3 3]); hold on
    set(ax1, 'Color', [0.7 0.7 0.7]); oel_mark_holds(holdSegs); oel_mark_visual(visEp); ylabel('EB wedge (angle around the ring)')
    title(sprintf('%s: z-scored EB wedge activity (background-subtracted, %d wedges, 3-D); dashed = held-bar segments (%d)', trialName, nClusters, size(holdSegs, 1)), 'Interpreter', 'none')
    cb = colorbar(ax1); cb.Label.String = 'z-score';
    ax2 = subplot(7, 1, 4:5); oel_shade_pulses(sync.led_on, sync.led_off, [0 360]); hold on
    plot(t_vol, bar_deg_plot, 'k', 'LineWidth', 1); ylim([0 360]); yticks(0:90:360); ylabel('bar position (deg)'); oel_mark_visual(visEp)
    title('bar position on the LED display (red shading = stim LED on)')
    ax3 = subplot(7, 1, 6); plot(t_vol, double(led_vol), 'r', 'LineWidth', 1); ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'off', 'on'}); ylabel('stim LED')
    ax4 = subplot(7, 1, 7); oel_draw_visual(visEp); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim(xl); for ax = [ax1 ax4], ax.Position([1 3]) = ax2.Position([1 3]); end
    oel_export(figure(2), fullfile(trialDir, sprintf('wedgeHeatmap_zscore_EB_%s.png', tag)));
end

%% 6. PVA -> bump vs bar (sign + lag fitted during the sweeps), bump_results_EB.mat
[xt, yt] = pol2cart(alpha, z_cluster');
[mu, rho] = cart2pol(mean(xt, 2, 'omitnan'), mean(yt, 2, 'omitnan')); mu = mu(:)'; rho = rho(:)';
bar_pm0 = angle(exp(1i * bar_rad_vol));
sweepOK = ~inHold_vol & rho >= rho_thresh & t_vol > sync.vr_start;
cs_same = oel_circ_std(mu(sweepOK) - bar_pm0(sweepOK)); cs_flip = oel_circ_std(mu(sweepOK) + bar_pm0(sweepOK));
if cs_flip < cs_same, bar_sign = -1; signNote = 'wedge order is MIRRORED relative to the display (bump ~ -bar)';
else, bar_sign = 1; signNote = 'wedge order runs WITH the display (bump ~ +bar)'; end
fprintf('bump-vs-bar: circ std of sweep offset %.2f rad (same) vs %.2f rad (flipped) -> %s\n', cs_same, cs_flip, signNote);
bar_c = exp(1i * bar_rad_vol); dt_vol = median(diff(t_vol)); lagScan = lag_range_sec(1):dt_vol:lag_range_sec(2); lagScore = nan(size(lagScan));
for L = 1:numel(lagScan)
    barL = angle(interp1(t_vol, bar_c, t_vol - lagScan(L), 'linear', 'extrap'));
    lagScore(L) = oel_circ_std(mu(sweepOK) - bar_sign * barL(sweepOK));
end
[~, iBest] = min(lagScore); lag_sec = lagScan(iBest); if isnan(lag_sec) || lag_sec < 0, lag_sec = 0; end
bar_rad_lag = angle(interp1(t_vol, bar_c, t_vol - lag_sec, 'linear', 'extrap')); bar_pm = angle(exp(1i * bar_rad_lag));
offset = angle(exp(1i * (mu - bar_sign * bar_pm)));
iUsed = find(abs(lagScan - lag_sec) == min(abs(lagScan - lag_sec)), 1); iZero = find(abs(lagScan) == min(abs(lagScan)), 1);
fprintf('  lag %.2f s (sweep offset circ std %.1f deg; %.1f at lag 0); mean sweep offset %.1f deg; hold LED off %.1f deg (%.1f), LED on %.1f deg (%.1f); mean rho sweep %.2f, hold off %.2f, hold on %.2f\n', ...
    lag_sec, rad2deg(lagScore(iUsed)), rad2deg(lagScore(iZero)), rad2deg(oel_circ_mean(offset(sweepOK))), ...
    rad2deg(oel_circ_mean(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), rad2deg(oel_circ_std(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), ...
    rad2deg(oel_circ_mean(offset(inHold_vol & led_vol & rho >= rho_thresh))), rad2deg(oel_circ_std(offset(inHold_vol & led_vol & rho >= rho_thresh))), ...
    mean(rho(sweepOK), 'omitnan'), mean(rho(inHold_vol & ~led_vol), 'omitnan'), mean(rho(inHold_vol & led_vol), 'omitnan'));
bump = struct('t_vol', t_vol, 'mu', mu, 'rho', rho, 'bar_rad', bar_rad_vol, 'bar_sign', bar_sign, 'bar_rad_lag', bar_rad_lag, 'lag_sec', lag_sec, ...
    'offset', offset, 'led', led_vol, 'inHold', inHold_vol, 'badVol', badVol, 'holdSegs', holdSegs, 'led_on', sync.led_on, 'led_off', sync.led_off, ...
    'vr_start', sync.vr_start, 'vr_end', sync.vr_end, 'rho_thresh', rho_thresh, 'z_cluster', z_cluster, 'f_cluster', f_cluster, ...
    'clusterIdx', clusterIdx, 'nCorePx', nCorePx, 'alpha', alpha, 'trialName', trialName, 'signNote', signNote, 'structure', 'EB'); %#ok<NASGU>
bump.visEpochs = visEp; bump.crop = crop; bump.planes = planes; bump.zPlanes = zPlanes;
bump.ring = rmfield(ring, 'angleMap');
save(fullfile(trialDir, 'bump_results_EB.mat'), 'bump'); fprintf('saved %s\n', fullfile(trialDir, 'bump_results_EB.mat'));
mu_plot = rad2deg(mu); mu_plot(rho < rho_thresh) = NaN; mu_plot(abs(diff([mu_plot mu_plot(end)])) > 180) = NaN;
barS_plot = rad2deg(angle(exp(1i * bar_sign * bar_pm))); barS_plot(abs(diff([barS_plot barS_plot(end)])) > 180) = NaN;
off_plot = rad2deg(offset); off_plot(rho < rho_thresh) = NaN;
for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
    figure(3); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(7, 1, 1:2); oel_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, barS_plot, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); plot(t_vol, mu_plot, 'b.', 'MarkerSize', 4); oel_mark_holds(holdSegs); oel_mark_visual(visEp)
    ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg)')
    legend({sprintf('bar (signed, lagged %.2f s)', lag_sec), sprintf('EB bump (PVA), rho >= %.1f', rho_thresh)}, 'Location', 'northeastoutside')
    title(sprintf('%s: EB bump vs bar -- %s; lag %.2f s; red = stim LED on', trialName, signNote, lag_sec), 'Interpreter', 'none')
    ax2 = subplot(7, 1, 3:4); oel_shade_pulses(sync.led_on, sync.led_off, [0 max(rho) * 1.05]); hold on
    plot(t_vol, rho, 'k'); yline(rho_thresh, 'r:'); ylabel('bump strength (rho)'); oel_mark_visual(visEp)
    ax3 = subplot(7, 1, 5:6); oel_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, off_plot, 'b.', 'MarkerSize', 4); yline(0, 'k:'); oel_mark_holds(holdSegs); oel_mark_visual(visEp)
    ylim([-180 180]); yticks(-180:90:180); ylabel(sprintf('bump - bar(t - %.2f s) (deg)', lag_sec))
    ax4 = subplot(7, 1, 7); oel_draw_visual(visEp); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim(xl); for ax = [ax2 ax3 ax4], ax.Position([1 3]) = ax1.Position([1 3]); end
    oel_export(figure(3), fullfile(trialDir, sprintf('bumpVsBar_EB_%s.png', tag)));
end

%% 7. stim on - off diff (pulses inside holds; LED-off = the rest of those holds): x-y projection + ring plane
if isempty(sync.led_on), fprintf('no LED pulses -- no diff image.\ndone.\n'); return; end
if nnz(pulseInHold) >= 3, pulseIdx = find(pulseInHold); offMask = inHold_vol; note = sprintf('%d pulses inside holds', numel(pulseIdx));
else, pulseIdx = 1:numel(sync.led_on); offMask = true(1, nVol); note = sprintf('all %d pulses', numel(pulseIdx)); end
onIdx = false(1, nVol); for k = pulseIdx, onIdx = onIdx | (t_vol >= sync.led_on(k) & t_vol < sync.led_off(k)); end
offIdx = offMask & ~led_vol & t_vol >= sync.led_on(pulseIdx(1)) - 15 & t_vol <= sync.led_off(pulseIdx(end)) + 15;
onIdx = onIdx & ~badVol; offIdx = offIdx & ~badVol;
if nnz(onIdx) < 5 || nnz(offIdx) < 5, fprintf('diff image skipped: %d on / %d off volumes.\ndone.\n', nnz(onIdx), nnz(offIdx)); return; end
dVol = reshape(mean(vox(:, onIdx), 2) - mean(vox(:, offIdx), 2), dims); dVol(~mask3) = NaN;
fprintf('diff (%s): %d LED-on vs %d LED-off volumes\n', note, nnz(onIdx), nnz(offIdx));
dWedge = accumarray(clusterIdx(mask3), dVol(mask3), [nClusters 1], @mean)';
dxy = mean(dVol, 3, 'omitnan');
figure(4); clf; set(gcf, 'Position', [50 50 1600 450], 'Color', 'w')
subplot(1, 3, 1); imagesc(max(meanVol, [], 3)); axis image; colormap(gca, bone); hold on; contour(mask, [0.5 0.5], 'r', 'LineWidth', 1); colorbar; title('x-y max projection + mask footprint')
subplot(1, 3, 2); dl = max(prctile(abs(dxy(mask)), 99), eps); imagesc(dxy, 'AlphaData', double(mask)); clim([-dl dl]); axis image; colormap(gca, rwb); set(gca, 'Color', [0.7 0.7 0.7]); colorbar
title(sprintf('mean(LED on) - mean(LED off), averaged over z inside the mask; %s', note), 'FontSize', 9)
subplot(1, 3, 3); bar(1:nClusters, dWedge, 'FaceColor', [0.8 0.3 0.3]); xlabel('EB wedge'); ylabel('\DeltaF (on - off)'); yline(0, 'k:')
title(sprintf('per-wedge on - off; %d on / %d off volumes', nnz(onIdx), nnz(offIdx)), 'FontSize', 9)
sgtitle([trialName ': EB stim on - off (red = brighter with LED on)'], 'Interpreter', 'none')
oel_export(figure(4), fullfile(trialDir, 'holdDiff_stimOn_minus_off_EB.png'));
fprintf('done.\n');

%% local functions
function crop = oel_auto_crop(vol, umPerPx, sigma_um, margin_um)
% EB bounding box from a short mean volume: smooth, Otsu, largest 26-connected
% blob, padded by margin_um. Falls back to the full frame if the blob covers
% most of the FOV (nothing distinctly brighter to crop to).
[nR, nC, ~] = size(vol); sig = max(sigma_um / umPerPx, 0.5);
vs = imgaussfilt3(double(vol), [sig sig 0.6]); vs = min(max(vs, prctile(vs, 5, 'all')), prctile(vs, 99.5, 'all'));
bw = vs > multithresh(vs(:), 1); cc = bwconncomp(bw, 26);
[~, iBig] = max(cellfun(@numel, cc.PixelIdxList)); [r, c, ~] = ind2sub(size(bw), cc.PixelIdxList{iBig});
m = ceil(margin_um / umPerPx);
crop = [max(min(r) - m, 1), min(max(r) + m, nR), max(min(c) - m, 1), min(max(c) + m, nC)];
if (crop(2) - crop(1) + 1) * (crop(4) - crop(3) + 1) > 0.8 * nR * nC
    warning('auto crop covers > 80 percent of the frame -- using the full frame; set oel_crop if the EB is only part of it.'); crop = [1 nR 1 nC];
end
end
function segs = oel_flat_segments(t, bar, tStart, fs, tolV, minSec)
barF = movmedian(bar(:), round(0.25 * fs)); win = round(2 * fs);
flat = (movmax(barF, win) - movmin(barF, win)) < tolV & t(:) > tStart;
d = diff([0; flat; 0]); s = find(d == 1); e = find(d == -1) - 1; keep = (e - s + 1) / fs >= minSec;
segs = [t(s(keep)), t(e(keep))]; if isempty(segs), segs = zeros(0, 2); end
end
function oel_export(fig, outFile)
drawnow; set(findall(fig, 'Type', 'axestoolbar'), 'Visible', 'off');
for ax = findall(fig, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow; exportgraphics(fig, outFile, 'Resolution', 150); fprintf('saved %s\n', outFile);
[p, n] = fileparts(outFile); vis = get(fig, 'Visible'); set(fig, 'Visible', 'on'); savefig(fig, fullfile(p, [n '.fig'])); set(fig, 'Visible', vis);
end
function oel_mark_holds(segs)
for k = 1:size(segs, 1), xline(segs(k,1), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); xline(segs(k,2), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); end
end
function oel_shade_pulses(onsets, offsets, yl)
for k = 1:numel(onsets), patch([onsets(k) offsets(k) offsets(k) onsets(k)], [yl(1) yl(1) yl(2) yl(2)], [1 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off'); end
ylim(yl)
end
function ep = oel_visual_epochs(trialDir, vrStart, vrEnd, patternNames)
% piecewise-constant visual scene: prefer the per-iteration log (actual
% values, DAQ clock), else the protocol csv (VR clock -> + vrStart)
ep = struct('t0', {}, 't1', {}, 'pattern_idx', {}, 'pattern_brightness', {}, 'dot_brightness', {}, 'is_cl', {}, 'label', {}, 'src', {});
t = []; src = ''; cols = {'pattern_idx', 'pattern_brightness', 'dot_brightness', 'pattern_is_cl', 'dot_is_cl'};
f = fullfile(trialDir, 'sync_info.mat');
vars = {}; if isfile(f), w = whos('-file', f); vars = {w.name}; end
if all(ismember({'iter', 'iter_times_on_daq'}, vars))
    S = load(f, 'iter', 'iter_times_on_daq');
    if all(isfield(S.iter, cols(1:3)))
        t = S.iter_times_on_daq(:); X = zeros(numel(t), numel(cols));
        for c = 1:numel(cols), if isfield(S.iter, cols{c}), v = double(S.iter.(cols{c})(:)); X(:, c) = v(1:numel(t)); end, end
        src = 'sync_info.mat iteration log';
    end
end
if isempty(src)
    csvs = dir(fullfile(trialDir, '*.csv')); csvs = csvs(~endsWith({csvs.name}, '_metadata.csv'));
    for k = 1:numel(csvs)
        T = readtable(fullfile(csvs(k).folder, csvs(k).name));
        if ~all(ismember([{'time'} cols(1:3)], T.Properties.VariableNames)), continue, end
        t = vrStart + T.time(:); X = zeros(numel(t), numel(cols));
        for c = 1:numel(cols), if ismember(cols{c}, T.Properties.VariableNames), X(:, c) = T.(cols{c})(:); end, end
        src = csvs(k).name; break
    end
end
if isempty(src), return, end
ok = ~isnan(t); t = t(ok); X = X(ok, :);
ch = [1; find(any(diff(X(:, 1:3)) ~= 0, 2)) + 1]; tEnd = [t(ch(2:end)); max(vrEnd, t(end))];
for k = 1:numel(ch)
    x = X(ch(k), :); parts = {};
    if x(1) > 0
        if isKey(patternNames, x(1)), nm = patternNames(x(1)); else, nm = sprintf('pattern %d', x(1)); end
        parts{end+1} = sprintf('%s @ %g', nm, x(2)); %#ok<AGROW>
    end
    if x(3) > 0, parts{end+1} = sprintf('dot @ %g', x(3)); end %#ok<AGROW>
    if isempty(parts), lbl = 'dark'; else, lbl = strjoin(parts, ' + '); end
    isCl = (x(1) > 0 && x(4) > 0) || (x(3) > 0 && x(5) > 0);
    if ~isempty(parts) && isCl, lbl = [lbl ' (CL)']; end
    ep(end+1) = struct('t0', t(ch(k)), 't1', tEnd(k), 'pattern_idx', x(1), 'pattern_brightness', x(2), 'dot_brightness', x(3), ...
        'is_cl', isCl, 'label', lbl, 'src', src); %#ok<AGROW>
end
end
function sync = oel_cue_from_iter_log(sync, trialDir)
f = fullfile(trialDir, 'sync_info.mat'); vars = {};
if isfile(f), w = whos('-file', f); vars = {w.name}; end
assert(all(ismember({'iter', 'iter_times_on_daq'}, vars)), ...
    'bar channel (vr_ao0_copy) is flat and %s has no iteration log -- cannot recover the cue angle.', f);
S = load(f, 'iter', 'iter_times_on_daq'); it = S.iter; t = S.iter_times_on_daq(:); n = numel(t);
deg = double(it.pattern_deg_ego(1:n)); deg = deg(:); src = 'pattern_deg_ego';
if isfield(it, 'dot_brightness') && isfield(it, 'dot_deg_ego')
    dotOn = double(it.dot_brightness(1:n)) > 0; dotDeg = double(it.dot_deg_ego(1:n));
    deg(dotOn) = dotDeg(dotOn); if any(dotOn), src = 'pattern_deg_ego / dot_deg_ego'; end
end
ok = ~isnan(t) & ~isnan(deg); t = t(ok); deg = deg(ok);
ang = interp1(t, unwrap(deg2rad(deg)), sync.t_fine, 'linear', NaN);
ang(sync.t_fine < t(1)) = ang(find(sync.t_fine >= t(1), 1)); ang(sync.t_fine > t(end)) = ang(find(sync.t_fine <= t(end), 1, 'last'));
sync.bar_rad = mod(ang, 2*pi); sync.vr_start = t(1); sync.vr_end = t(end); sync.barSource = ['iteration log ' src];
end
function oel_draw_visual(ep)
% one patch per epoch, hue = what is shown, lightness = brightness (0-255)
hold on
if isempty(ep), text(0.5, 0.5, 'no visual-pattern log found', 'Units', 'normalized', 'HorizontalAlignment', 'center'); end
pal = lines(7); dotCol = [1 0.75 0.1];
for k = 1:numel(ep)
    e = ep(k);
    if e.pattern_idx > 0, base = pal(mod(e.pattern_idx - 1, 7) + 1, :); b = e.pattern_brightness;
    elseif e.dot_brightness > 0, base = dotCol; b = e.dot_brightness;
    else, base = [0 0 0]; b = 0; end
    w = 0.25 + 0.75 * min(b / 255, 1); col = w * base + (1 - w) * [1 1 1];
    if b == 0, col = [0.15 0.15 0.15]; end
    patch([e.t0 e.t1 e.t1 e.t0], [0 0 1 1], col, 'EdgeColor', 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');
    if b == 0, tc = 'w'; else, tc = 'k'; end
    text((e.t0 + e.t1) / 2, 0.5, e.label, 'HorizontalAlignment', 'center', 'Color', tc, 'FontSize', 8, 'Clipping', 'on', 'Interpreter', 'none');
end
ylim([0 1]); yticks([]); ylabel('visual')
end
function oel_mark_visual(ep)
for k = 2:numel(ep), xline(ep(k).t0, ':', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.25, 'HandleVisibility', 'off'); end
end
function m = oel_circ_mean(x), m = angle(mean(exp(1i * x(:)), 'omitnan')); end
function s = oel_circ_std(x)
x = x(~isnan(x)); if isempty(x), s = NaN; return; end
R = abs(mean(exp(1i * x(:)))); s = sqrt(-2 * log(max(R, eps)));
end
function cmap = oel_redwhiteblue(n)
half = n / 2; cmap = [[linspace(0,1,half)', linspace(0,1,half)', ones(half,1)]; [ones(half,1), linspace(1,0,half)', linspace(1,0,half)']];
end
