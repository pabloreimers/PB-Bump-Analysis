%% overshoot_pb_led_script
% PB-only bar + stim-LED trials recorded with the Oct-2026 rig configuration
% (20261013-* flies onward): 13 frames/volume (10 slices at 25-61 um + 3
% flyback), 100 x 256 px, zoom 4, 10.5 Hz; WaveSurfer logs the stim LED only
% as the ANALOG line vr_led_stim_copy (no led_stim_feedback DI channel), and
% there is usually no scan_params.json (geometry comes from the tif header).
%
% Same analysis as overshoot_bar_stim_script.m (the Sep/Oct-5 PB trials),
% rebuilt on the shared helpers in claude\ that the PB+EB pipeline uses:
%   pb_read_raw_tif_planes  BigTIFF IFD walker (libtiff dies at 65535 frames),
%                           sums the valid planes per volume
%   pb_register             NoRMCorre rigid, maxShift [6 8], with the
%                           outlier-shift cleanup (spurious +/-15-20 px locks
%                           painted a ghost PB into the mean on 20261012-2)
%   pb_arch_mask            column-normalised threshold (dim hemisphere ok),
%                           non-PB blobs rejected by distance from the arch
%   pb_h5_sync              volume clock, bar (vr_ao0_copy, 1 V/rad), LED
%                           from DI if present else the analog copy
% Glomeruli: 2N along the skeleton (obs_pb_glomeruli, N per hemisphere,
% alpha repeated), per-volume background subtraction (mean outside the
% mask), z-score, PVA; sign and lag vs the bar fitted during the sweeps.
%
% Outputs in the trial folder (all cached stages skip if the file exists):
%   imgData_sum.mat, imgData_sum_reg.mat, mask.mat, mask_qc.png
%   glomHeatmap_zscore_{fullTrial,stim}.png, bumpVsBar_{fullTrial,stim}.png
%   holdDiff_stimOn_minus_off.png, bump_results.mat (same fields as before,
%   so fig_bump_vs_bar_holds.m, fig_bump_vs_cue_scatter_blocks.m,
%   overshoot_ttx_early_vs_late.m, overshoot_crosstrial_glomerulus_alignment.m
%   all work unchanged)
%
% ---- protocols seen so far ----
%   dangling_burn_in_sweep36_holds10_305s: 36 deg/s sweeps alternating
%     direction with 10-s holds, no LED.
%   dangling_burn_in_4pos_ledstim_pre10s_2rev_endburnin_652s: the same
%     burn-in to VR 185 s, then 8 holds of 30 s (bar at 0/90/180/270 deg
%     going forward, then 180/90/0/270 coming back; each 10 s pre, 10 s LED
%     ON, 10 s post) separated by 810-deg sweeps, then two burn-in sweeps.
%     LED = 8 pulses x 10 s; e.g. 20261013-2 trial_002: DAQ 239-621 s,
%     6959 volumes.
%   lpsp_cl_multigauss150_dark_cl_dot150_360s (20261021-2 onward): no bar
%     sweeps; 3 x 120-s closed-loop blocks with a different visual scene
%     each: multigauss pattern (pattern_idx 2) at brightness 150, dark, then
%     the dot at brightness 150.
%
% Visual epochs (which pattern / dot is shown, and at what brightness) come
% from the per-iteration ViRMEn log in sync_info.mat (pattern_idx,
% pattern_brightness, dot_brightness on the DAQ clock) or, if that is
% missing, from the protocol .csv in the trial folder (time column is VR
% time, so offset by vr_start). They are drawn as a strip under figures 2
% and 3, dotted lines mark the epoch boundaries, and they go into
% bump_results.mat as bump.visEpochs. pattern_idx -> name is not logged
% anywhere, so pattern_names below is filled in by hand.
%
% Headless:
%   matlab -batch "opl_trialDir='Z:\...\trial_002\'; run('overshoot_pb_led_script.m')"
% Overrides: opl_overwrite_mask / opl_overwrite_reg / opl_overwrite_video,
% opl_min_hold_sec (default 4; the 10-s protocol holds erode to ~8 s, the
% 30-s LED holds to ~28 s, so 4 catches both).

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'claude'));
if isempty(which('graph_sort')), addpath(fullfile(thisDir, '..')); end
if isempty(which('circ_mean')), addpath(fullfile(thisDir, '..', 'circ_stats')); end

trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20261023-3_epg_syt8s_lpsp_cschrimson\trial_004\';
if exist('opl_trialDir', 'var') && ~isempty(opl_trialDir), trialDir = opl_trialDir; end
if trialDir(end) ~= '\' && trialDir(end) ~= '/', trialDir(end+1) = filesep; end
parts = strsplit(trialDir, filesep); parts = parts(~cellfun(@isempty, parts));
trialName = [parts{end-1} filesep parts{end}];

n_per_hemisphere = 16;
mask_smooth_sigma_um = 1.5; mask_thresh_pct = 88; mask_open_radius_um = 1.0; mask_join_line_um = 15;
maxShift_px = [6 8]; max_ok_shift_px = 8;
min_hold_sec = 4; hold_flat_tol_V = 0.1;
if exist('opl_min_hold_sec', 'var') && ~isempty(opl_min_hold_sec), min_hold_sec = opl_min_hold_sec; end
lag_range_sec = [-2 5]; rho_thresh = 0.1;
pattern_names = containers.Map({2}, {'multigauss'});    % pattern_idx -> label; unknown idx shows as 'pattern N'
overwrite_video = false; overwrite_reg = false; overwrite_mask = false;
if exist('opl_overwrite_video', 'var') && ~isempty(opl_overwrite_video), overwrite_video = opl_overwrite_video; end
if exist('opl_overwrite_reg', 'var') && ~isempty(opl_overwrite_reg), overwrite_reg = opl_overwrite_reg; end
if exist('opl_overwrite_mask', 'var') && ~isempty(opl_overwrite_mask), overwrite_mask = opl_overwrite_mask; end
if overwrite_video, overwrite_reg = true; end
if overwrite_reg, overwrite_mask = true; end

% scan geometry from the tif header (scan_params.json, if present, is only
% used as a cross-check for the plane layout)
tifList = dir(fullfile(trialDir, '*.tif')); assert(numel(tifList) == 1, 'expected exactly one .tif in %s', trialDir);
tifPath = fullfile(tifList(1).folder, tifList(1).name);
fid = fopen(tifPath, 'r'); hdr = fread(fid, 80000, 'uint8=>char')'; fclose(fid);
getNum = @(key) str2double(regexp(hdr, ['SI\.' regexptranslate('escape', key) '\s*=\s*([-0-9.eE]+)'], 'tokens', 'once'));
framesPerVol = getNum('hStackManager.numFramesPerVolumeWithFlyback'); nValid = getNum('hStackManager.numFramesPerVolume');
zoom = getNum('hRoiManager.scanZoomFactor'); W = getNum('hRoiManager.pixelsPerLine'); H = getNum('hRoiManager.linesPerFrame');
assert(all(~isnan([framesPerVol nValid zoom W H])), 'could not parse the ScanImage header of %s', tifPath);
validPlanes = 1:nValid;                                   % flyback frames are the last ones in the volume
umPerPx = 0.3164 * (4 / zoom) * (256 / W);               % 0.3164 um/px at zoom 4, 256 px
jsonFile = fullfile(trialDir, 'scan_params.json');
if isfile(jsonFile)
    sp = jsondecode(fileread(jsonFile));
    if sp.nplanes ~= framesPerVol || sp.nplanes_valid ~= nValid
        warning('scan_params.json (%d planes, %d valid) disagrees with the tif header (%d, %d); using the header.', sp.nplanes, sp.nplanes_valid, framesPerVol, nValid);
    end
end
fprintf('[%s] %dx%d px, %d frames/volume (%d valid), zoom %g -> %.4f um/px\n', trialName, H, W, framesPerVol, nValid, zoom, umPerPx);

%% 1. raw tif -> z-summed video (cached)
sumFile = fullfile(trialDir, 'imgData_sum.mat');
if overwrite_video || ~isfile(sumFile)
    fprintf('[%s] reading raw tif...\n', trialName);
    o = pb_read_raw_tif_planes(tifPath, framesPerVol, {validPlanes});
    imgData_sum = o.stacks{1}; nVolumes = o.nVolumes; %#ok<NASGU>
    sp = struct('nplanes', framesPerVol, 'nplanes_valid', nValid, 'px_height', H, 'px_width', W, 'zoom', zoom, 'um_per_px', [umPerPx umPerPx NaN], ...
        'nVolumes', o.nVolumes, 'nFramesRead', o.nFramesRead, 'tifChainOk', o.chainOk, 'provenance', 'tif header'); %#ok<NASGU>
    save(sumFile, 'imgData_sum', 'sp', '-v7.3'); fprintf('  saved %s (%d x %d x %d)\n', sumFile, H, W, o.nVolumes);
    clear o imgData_sum
else
    fprintf('skip (exists): %s\n', sumFile);
end

%% 2. motion correction (cached)
regFile = fullfile(trialDir, 'imgData_sum_reg.mat');
if overwrite_reg || ~isfile(regFile)
    S = load(sumFile, 'imgData_sum');
    fprintf('[%s] motion-correcting (%d volumes)...\n', trialName, size(S.imgData_sum, 3));
    [imgData_sum_reg, regShifts] = pb_register(S.imgData_sum, 'smoothWindow', 5, 'maxShift', maxShift_px, 'verbose', false, 'figureVisible', false);
    sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), regShifts(:), 'UniformOutput', false));
    fprintf('  shift range row = [%.1f %.1f] px, col = [%.1f %.1f] px\n', min(sh(:,1)), max(sh(:,1)), min(sh(:,2)), max(sh(:,2)));
    save(regFile, 'imgData_sum_reg', 'regShifts', '-v7.3'); clear S imgData_sum_reg regShifts
else
    fprintf('skip (exists): %s\n', regFile);
end

%% 3. mask + glomeruli (cached mask), QC figure
maskFile = fullfile(trialDir, 'mask.mat');
S = load(regFile, 'imgData_sum_reg'); projMean = mean(S.imgData_sum_reg, 3); clear S
if overwrite_mask || ~isfile(maskFile)
    [mask, maskCore, pbInfo] = pb_arch_mask(projMean, umPerPx, 'sigma_um', mask_smooth_sigma_um, 'thresh_pct', mask_thresh_pct, ...
        'open_um', mask_open_radius_um, 'join_um', mask_join_line_um);
    fprintf('PB mask: %s\n', pbInfo.summary);
    save(maskFile, 'mask', 'maskCore', 'pbInfo'); fprintf('saved %s\n', maskFile);
else
    fprintf('skip (exists): %s\n', maskFile);
end
M = load(maskFile); mask = M.mask; maskCore = M.maskCore; clear M
clusterIdx = obs_pb_glomeruli(mask, n_per_hemisphere);
nClusters = 2 * n_per_hemisphere;
alpha = repmat(linspace(-pi, pi - 2*pi/n_per_hemisphere, n_per_hemisphere), 1, 2);
figure(1); clf; set(gcf, 'Position', [50 50 1400 600], 'Color', 'w')
subplot(2, 1, 1); imagesc(projMean); axis image; colormap(gca, bone); hold on
contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5); if any(mask(:) & ~maskCore(:)), contour(maskCore, [0.5 0.5], 'c', 'LineWidth', 1); end
title(sprintf('%s: mean registered z-sum + mask (%d px)', trialName, nnz(mask)), 'Interpreter', 'none')
subplot(2, 1, 2); imagesc(clusterIdx); axis image; colormap(gca, [0 0 0; hsv(nClusters)]); hold on
for c = 1:nClusters, [cy, cx] = find(clusterIdx == c); if ~isempty(cx), text(mean(cx), mean(cy), num2str(c), 'Color', 'w', 'FontSize', 7, 'HorizontalAlignment', 'center', 'FontWeight', 'bold'); end, end
title(sprintf('%d glomeruli along the skeleton (16 per hemisphere), alpha = -pi at #1', nClusters))
opl_export(figure(1), fullfile(trialDir, 'mask_qc.png'));

%% 4. sync (volume times, LED, bar, holds), bad volumes, background subtraction
h5list = dir(fullfile(trialDir, '*.h5')); assert(numel(h5list) == 1, 'expected exactly one .h5 in %s', trialDir);
sync = pb_h5_sync(fullfile(h5list(1).folder, h5list(1).name));
fprintf('LED source: %s\n', sync.ledSource);
if isempty(sync.vr_start)
    % vr_ao0_copy only carries the DOT angle: with the dot off (e.g. the
    % open-loop dangling_pattern_* sweeps) it sits at 0 V and vr_start/vr_end
    % come back empty. Take the cue angle from the ViRMEn iteration log
    % instead (pattern_deg_ego, or dot_deg_ego while the dot is shown).
    sync = opl_cue_from_iter_log(sync, trialDir);
    fprintf('bar channel flat -> cue angle from the iteration log (%s); VR %.1f-%.1f s DAQ\n', sync.barSource, sync.vr_start, sync.vr_end);
end
S = load(regFile, 'imgData_sum_reg', 'regShifts'); img = single(S.imgData_sum_reg);
sh = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), S.regShifts(:), 'UniformOutput', false)); clear S
badVol = any(abs(sh) > max_ok_shift_px, 2)';
nVol = min(size(img, 3), numel(sync.t_volume));
if size(img, 3) ~= numel(sync.t_volume), fprintf('note: %d volumes in the video vs %d volume clocks -> using %d\n', size(img, 3), numel(sync.t_volume), nVol); end
img = img(:, :, 1:nVol); badVol = badVol(1:nVol); t_vol = sync.t_volume(1:nVol)';
fprintf('%d volumes (DAQ %.1f-%.1f s); %d with |shift| > %d px excluded\n', nVol, t_vol(1), t_vol(end), nnz(badVol), max_ok_shift_px);
led_vol = interp1(sync.t_fine, double(sync.led), t_vol, 'nearest', 'extrap') > 0.5;
bar_rad_vol = mod(interp1(sync.t_fine, sync.bar_rad, t_vol, 'nearest', 'extrap'), 2*pi); bar_deg_vol = rad2deg(bar_rad_vol);
holdSegs = opl_flat_segments(sync.t_fine, sync.bar_rad, sync.vr_start, sync.fs, hold_flat_tol_V, min_hold_sec);
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
visEp = opl_visual_epochs(trialDir, sync.vr_start, sync.vr_end, pattern_names);
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
img_2d = reshape(img, [], nVol); bg = mean(img_2d(~mask(:), :), 1); img = img - reshape(single(bg), 1, 1, nVol); clear img_2d
rwb = opl_redwhiteblue(256);
bar_deg_plot = bar_deg_vol; bar_deg_plot(abs(diff([bar_deg_plot bar_deg_plot(end)])) > 180) = NaN;

%% 5. glomerulus traces -> z-score -> heatmap
img_2d = reshape(img, [], nVol); f_cluster = nan(nClusters, nVol); nCorePx = zeros(nClusters, 1);
for c = 1:nClusters
    pix = clusterIdx(:) == c & maskCore(:); nCorePx(c) = nnz(pix);
    if nCorePx(c) > 0, f_cluster(c, :) = mean(img_2d(pix, :), 1); end
end
clear img_2d
emptyGlom = nCorePx < max(10, 0.2 * median(nCorePx(nCorePx > 0)));
if any(emptyGlom), f_cluster(emptyGlom, :) = NaN; fprintf('%d/%d glomeruli with (almost) no signal pixels -> NaN: %s\n', nnz(emptyGlom), nClusters, mat2str(find(emptyGlom)')); end
f_cluster(:, badVol) = NaN;
z_cluster = (f_cluster - mean(f_cluster, 2, 'omitnan')) ./ std(f_cluster, 0, 2, 'omitnan');
for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
    figure(2); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(7, 1, 1:3);
    imagesc(t_vol, 1:nClusters, z_cluster, 'AlphaData', double(~isnan(z_cluster))); colormap(ax1, rwb); clim([-3 3]); hold on
    set(ax1, 'Color', [0.7 0.7 0.7]); opl_mark_holds(holdSegs); opl_mark_visual(visEp); ylabel('glomerulus (1-16 left, 17-32 right)')
    title(sprintf('%s: z-scored glomerulus activity (background-subtracted, %d glomeruli); dashed = held-bar segments (%d)', trialName, nClusters, size(holdSegs, 1)), 'Interpreter', 'none')
    cb = colorbar(ax1); cb.Label.String = 'z-score';
    ax2 = subplot(7, 1, 4:5); opl_shade_pulses(sync.led_on, sync.led_off, [0 360]); hold on
    plot(t_vol, bar_deg_plot, 'k', 'LineWidth', 1); ylim([0 360]); yticks(0:90:360); ylabel('bar position (deg)'); opl_mark_visual(visEp)
    title('bar position on the LED display (red shading = stim LED on)')
    ax3 = subplot(7, 1, 6); plot(t_vol, double(led_vol), 'r', 'LineWidth', 1); ylim([-0.1 1.1]); yticks([0 1]); yticklabels({'off', 'on'}); ylabel('stim LED')
    ax4 = subplot(7, 1, 7); opl_draw_visual(visEp); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim(xl); for ax = [ax1 ax4], ax.Position([1 3]) = ax2.Position([1 3]); end
    opl_export(figure(2), fullfile(trialDir, sprintf('glomHeatmap_zscore_%s.png', tag)));
end

%% 6. PVA -> bump vs bar (sign + lag fitted during the sweeps), bump_results.mat
[xt, yt] = pol2cart(alpha, z_cluster');
[mu, rho] = cart2pol(mean(xt, 2, 'omitnan'), mean(yt, 2, 'omitnan')); mu = mu(:)'; rho = rho(:)';
bar_pm0 = angle(exp(1i * bar_rad_vol));
sweepOK = ~inHold_vol & rho >= rho_thresh & t_vol > sync.vr_start;
cs_same = opl_circ_std(mu(sweepOK) - bar_pm0(sweepOK)); cs_flip = opl_circ_std(mu(sweepOK) + bar_pm0(sweepOK));
if cs_flip < cs_same, bar_sign = -1; signNote = 'glomerulus order is MIRRORED relative to the display (bump ~ -bar)';
else, bar_sign = 1; signNote = 'glomerulus order runs WITH the display (bump ~ +bar)'; end
fprintf('bump-vs-bar: circ std of sweep offset %.2f rad (same) vs %.2f rad (flipped) -> %s\n', cs_same, cs_flip, signNote);
bar_c = exp(1i * bar_rad_vol); dt_vol = median(diff(t_vol)); lagScan = lag_range_sec(1):dt_vol:lag_range_sec(2); lagScore = nan(size(lagScan));
for L = 1:numel(lagScan)
    barL = angle(interp1(t_vol, bar_c, t_vol - lagScan(L), 'linear', 'extrap'));
    lagScore(L) = opl_circ_std(mu(sweepOK) - bar_sign * barL(sweepOK));
end
[~, iBest] = min(lagScore); lag_sec = lagScan(iBest); if isnan(lag_sec) || lag_sec < 0, lag_sec = 0; end
bar_rad_lag = angle(interp1(t_vol, bar_c, t_vol - lag_sec, 'linear', 'extrap')); bar_pm = angle(exp(1i * bar_rad_lag));
offset = angle(exp(1i * (mu - bar_sign * bar_pm)));
iUsed = find(abs(lagScan - lag_sec) == min(abs(lagScan - lag_sec)), 1); iZero = find(abs(lagScan) == min(abs(lagScan)), 1);
fprintf('  lag %.2f s (sweep offset circ std %.1f deg; %.1f at lag 0); mean sweep offset %.1f deg; hold LED off %.1f deg (%.1f), LED on %.1f deg (%.1f); mean rho sweep %.2f, hold off %.2f, hold on %.2f\n', ...
    lag_sec, rad2deg(lagScore(iUsed)), rad2deg(lagScore(iZero)), rad2deg(opl_circ_mean(offset(sweepOK))), ...
    rad2deg(opl_circ_mean(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), rad2deg(opl_circ_std(offset(inHold_vol & ~led_vol & rho >= rho_thresh))), ...
    rad2deg(opl_circ_mean(offset(inHold_vol & led_vol & rho >= rho_thresh))), rad2deg(opl_circ_std(offset(inHold_vol & led_vol & rho >= rho_thresh))), ...
    mean(rho(sweepOK), 'omitnan'), mean(rho(inHold_vol & ~led_vol), 'omitnan'), mean(rho(inHold_vol & led_vol), 'omitnan'));
bump = struct('t_vol', t_vol, 'mu', mu, 'rho', rho, 'bar_rad', bar_rad_vol, 'bar_sign', bar_sign, 'bar_rad_lag', bar_rad_lag, 'lag_sec', lag_sec, ...
    'offset', offset, 'led', led_vol, 'inHold', inHold_vol, 'badVol', badVol, 'holdSegs', holdSegs, 'led_on', sync.led_on, 'led_off', sync.led_off, ...
    'vr_start', sync.vr_start, 'vr_end', sync.vr_end, 'rho_thresh', rho_thresh, 'z_cluster', z_cluster, 'f_cluster', f_cluster, ...
    'clusterIdx', clusterIdx, 'nCorePx', nCorePx, 'alpha', alpha, 'trialName', trialName, 'signNote', signNote, 'structure', 'PB'); %#ok<NASGU>
bump.visEpochs = visEp;
save(fullfile(trialDir, 'bump_results.mat'), 'bump'); fprintf('saved %s\n', fullfile(trialDir, 'bump_results.mat'));
mu_plot = rad2deg(mu); mu_plot(rho < rho_thresh) = NaN; mu_plot(abs(diff([mu_plot mu_plot(end)])) > 180) = NaN;
barS_plot = rad2deg(angle(exp(1i * bar_sign * bar_pm))); barS_plot(abs(diff([barS_plot barS_plot(end)])) > 180) = NaN;
off_plot = rad2deg(offset); off_plot(rho < rho_thresh) = NaN;
for zoomPass = 1:(1 + ~isempty(zoomWin))
    if zoomPass == 1, xl = [t_vol(1) t_vol(end)]; tag = 'fullTrial'; else, xl = zoomWin; tag = 'stim'; end
    figure(3); clf; set(gcf, 'Position', [50 50 1600 900], 'Color', 'w')
    ax1 = subplot(7, 1, 1:2); opl_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, barS_plot, '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5); plot(t_vol, mu_plot, 'b.', 'MarkerSize', 4); opl_mark_holds(holdSegs); opl_mark_visual(visEp)
    ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg)')
    legend({sprintf('bar (signed, lagged %.2f s)', lag_sec), sprintf('bump (PVA), rho >= %.1f', rho_thresh)}, 'Location', 'northeastoutside')
    title(sprintf('%s: bump vs bar -- %s; lag %.2f s; red = stim LED on', trialName, signNote, lag_sec), 'Interpreter', 'none')
    ax2 = subplot(7, 1, 3:4); opl_shade_pulses(sync.led_on, sync.led_off, [0 max(rho) * 1.05]); hold on
    plot(t_vol, rho, 'k'); yline(rho_thresh, 'r:'); ylabel('bump strength (rho)'); opl_mark_visual(visEp)
    ax3 = subplot(7, 1, 5:6); opl_shade_pulses(sync.led_on, sync.led_off, [-180 180]); hold on
    plot(t_vol, off_plot, 'b.', 'MarkerSize', 4); yline(0, 'k:'); opl_mark_holds(holdSegs); opl_mark_visual(visEp)
    ylim([-180 180]); yticks(-180:90:180); ylabel(sprintf('bump - bar(t - %.2f s) (deg)', lag_sec))
    ax4 = subplot(7, 1, 7); opl_draw_visual(visEp); xlabel('DAQ time (s)')
    linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim(xl); for ax = [ax2 ax3 ax4], ax.Position([1 3]) = ax1.Position([1 3]); end
    opl_export(figure(3), fullfile(trialDir, sprintf('bumpVsBar_%s.png', tag)));
end

%% 7. stim on - off diff image (pulses inside holds; LED-off = the rest of those holds)
if isempty(sync.led_on), fprintf('no LED pulses -- no diff image.\ndone.\n'); return; end
if nnz(pulseInHold) >= 3, pulseIdx = find(pulseInHold); offMask = inHold_vol; note = sprintf('%d pulses inside holds', numel(pulseIdx));
else, pulseIdx = 1:numel(sync.led_on); offMask = true(1, nVol); note = sprintf('all %d pulses', numel(pulseIdx)); end
onIdx = false(1, nVol); for k = pulseIdx, onIdx = onIdx | (t_vol >= sync.led_on(k) & t_vol < sync.led_off(k)); end
offIdx = offMask & ~led_vol & t_vol >= sync.led_on(pulseIdx(1)) - 15 & t_vol <= sync.led_off(pulseIdx(end)) + 15;
onIdx = onIdx & ~badVol; offIdx = offIdx & ~badVol;
if nnz(onIdx) < 5 || nnz(offIdx) < 5, fprintf('diff image skipped: %d on / %d off volumes.\ndone.\n', nnz(onIdx), nnz(offIdx)); return; end
diffImg = mean(img(:, :, onIdx), 3) - mean(img(:, :, offIdx), 3);
fprintf('diff image (%s): %d LED-on vs %d LED-off volumes\n', note, nnz(onIdx), nnz(offIdx));
figure(4); clf; set(gcf, 'Position', [50 50 1400 450], 'Color', 'w')
subplot(1, 2, 1); imagesc(projMean); axis image; colormap(gca, bone); hold on; contour(mask, [0.5 0.5], 'r', 'LineWidth', 1); colorbar; title('mean image + mask')
subplot(1, 2, 2); dl = max(prctile(abs(diffImg(mask)), 99), eps); imagesc(diffImg, [-dl dl]); axis image; colormap(gca, rwb); hold on
contour(mask, [0.5 0.5], 'k', 'LineWidth', 0.75); colorbar
title(sprintf('mean(LED on) - mean(LED off): %s; %d on / %d off volumes', note, nnz(onIdx), nnz(offIdx)), 'FontSize', 9)
sgtitle([trialName ': stim on - off diff image (red = brighter with LED on)'], 'Interpreter', 'none')
opl_export(figure(4), fullfile(trialDir, 'holdDiff_stimOn_minus_off.png'));
fprintf('done.\n');

%% local functions
function segs = opl_flat_segments(t, bar, tStart, fs, tolV, minSec)
barF = movmedian(bar(:), round(0.25 * fs)); win = round(2 * fs);
flat = (movmax(barF, win) - movmin(barF, win)) < tolV & t(:) > tStart;
d = diff([0; flat; 0]); s = find(d == 1); e = find(d == -1) - 1; keep = (e - s + 1) / fs >= minSec;
segs = [t(s(keep)), t(e(keep))]; if isempty(segs), segs = zeros(0, 2); end
end
function opl_export(fig, outFile)
drawnow; set(findall(fig, 'Type', 'axestoolbar'), 'Visible', 'off');
for ax = findall(fig, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow; exportgraphics(fig, outFile, 'Resolution', 150); fprintf('saved %s\n', outFile);
[p, n] = fileparts(outFile); vis = get(fig, 'Visible'); set(fig, 'Visible', 'on'); savefig(fig, fullfile(p, [n '.fig'])); set(fig, 'Visible', vis);
end
function opl_mark_holds(segs)
for k = 1:size(segs, 1), xline(segs(k,1), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); xline(segs(k,2), 'k--', 'LineWidth', 0.75, 'HandleVisibility', 'off'); end
end
function opl_shade_pulses(onsets, offsets, yl)
for k = 1:numel(onsets), patch([onsets(k) offsets(k) offsets(k) onsets(k)], [yl(1) yl(1) yl(2) yl(2)], [1 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.35, 'HandleVisibility', 'off'); end
ylim(yl)
end
function ep = opl_visual_epochs(trialDir, vrStart, vrEnd, patternNames)
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
function sync = opl_cue_from_iter_log(sync, trialDir)
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
function opl_draw_visual(ep)
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
function opl_mark_visual(ep)
for k = 2:numel(ep), xline(ep(k).t0, ':', 'Color', [0.2 0.6 0.2], 'LineWidth', 1.25, 'HandleVisibility', 'off'); end
end
function m = opl_circ_mean(x), m = angle(mean(exp(1i * x(:)), 'omitnan')); end
function s = opl_circ_std(x)
x = x(~isnan(x)); if isempty(x), s = NaN; return; end
R = abs(mean(exp(1i * x(:)))); s = sqrt(-2 * log(max(R, eps)));
end
function cmap = opl_redwhiteblue(n)
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
