%% gain_change_example_fly
% One closed-loop gain-change trial, end to end, section by section:
%   raw ScanImage tif -> z-summed movie -> motion correction -> automatic PB
%   mask -> glomeruli along the PB skeleton -> z-scored activity per
%   glomerulus -> bump position (PVA) -> compared with the bar and with the
%   fly's heading from fictrac.
%
% Built from data scripts/overshoot_bar_stim_script.m (sections 1-3, 5, 6)
% with the WaveSurfer h5 sync replaced by the trial's ficTracData_DAQ.mat
% (berg1 gain_change rig). Everything this script writes goes to a
% claude_pipeline\ subfolder of the trial, so the hand-drawn mask.mat and the
% old registration_001\ in the trial folder are left untouched.
%
% Default trial: 20230110 fly 1, trial 1 = closed loop, bar, gain 0.7
% (first trial of the fly; see data scripts/gain_change_fly_inventory.csv).
%
% Run one %% section at a time, or in batch:
%   matlab -batch "cd('<repo>'); addpath('data scripts'); gain_change_example_fly"

%% 0. parameters
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(repoRootDir, 'claude'));
addpath(fullfile(repoRootDir, 'circ_stats'));
addpath(repoRootDir); % graph_sort.m

trialDir = 'Z:\pablo\gain_change\20230110\fly 1\20230110-1_EPG_7f_.7\';
if exist('gce_trialDir', 'var') && ~isempty(gce_trialDir), trialDir = gce_trialDir; end
outDir = fullfile(trialDir, 'claude_pipeline');
if ~isfolder(outDir), mkdir(outDir); end
exportDir = fullfile(repoRootDir, 'ugly_figures', 'exports', 'gain_change_example');
if ~isfolder(exportDir), mkdir(exportDir); end

pathParts = strsplit(trialDir, filesep); pathParts = pathParts(~cellfun(@isempty, pathParts));
trialName = [pathParts{end-2} ' ' pathParts{end-1} ' ' pathParts{end}];

overwrite_video = false;
overwrite_reg   = false;
overwrite_mask  = false;

n_per_hemisphere = 16; % process_im convention (n_centroid = 16): 32 bins across the PB, alpha repeated over the halves
max_ok_shift_px  = 10; % volumes whose rigid shift exceeds this are NaN'd out of the traces

% mask parameters in physical units (same scheme as overshoot_bar_stim_script;
% this rig images at ~0.49 um/px, zoom 5)
mask_smooth_sigma_um = 1.2; % ~2.5 px here; tested best on the CV image (see section 3)
mask_open_radius_um  = 1;
mask_join_line_um    = 15;
mask_thresh_pct      = 75; % 80 missed the upper tip of one hemisphere on 20230110 fly 1; 75 + 1 px dilation recovers it
mask_dilate_px       = 1;

% closed-loop sanity exclusion: volumes where the bar has not moved for
% bar_frozen_win_s although the fly is turning faster than frozen_yaw_thresh
% (bar parked at the end of 20230110-1 while the bump kept moving). Those
% volumes are dropped from every bump-vs-bar / bump-vs-heading statistic.
bar_frozen_win_s  = 3;    % window over which the bar must be flat
bar_frozen_tol    = deg2rad(2.5); % flat = total bar excursion below this (bar is quantized at 1.875 deg)
frozen_head_tol   = deg2rad(8); % fly "turned" = heading excursion over the window above this (bar should then have moved > 5 deg at gain 0.7)

rho_thresh    = 0.1;      % PVA length below which the bump position is not plotted
lag_range_sec = [-2 5];   % cue-bump lag scan (positive = bump follows the cue)
ft_rate       = 60;       % fictrac/DAQ table rate (Hz)

%% 1. raw tif -> z-summed video (cached)
sumFile = fullfile(outDir, 'imgData_sum.mat');
if ~isfile(sumFile) || overwrite_video
    fprintf('[%s] reading raw tif...\n', trialName);
    [imgData_sum, sp] = gce_read_raw_tif_summed(trialDir);
    save(sumFile, 'imgData_sum', 'sp', '-v7.3');
    fprintf('  saved %s  (%d x %d x %d volumes)\n', sumFile, size(imgData_sum,1), size(imgData_sum,2), size(imgData_sum,3));
else
    fprintf('skip (exists): %s\n', sumFile);
end
S = load(sumFile, 'sp'); sp = S.sp; clear S

%% 2. motion correction (cached)
regFile = fullfile(outDir, 'imgData_sum_reg.mat');
if ~isfile(regFile) || overwrite_reg
    S = load(sumFile, 'imgData_sum');
    fprintf('[%s] motion-correcting (%d volumes)...\n', trialName, size(S.imgData_sum,3));
    [imgData_sum_reg, regShifts, regDiag] = pb_register(S.imgData_sum, ...
        'smoothWindow', 5, 'maxShift', [30 45], 'verbose', false, 'figureVisible', false); %#ok<ASGLU>
    close(regDiag.figHandles);
    save(regFile, 'imgData_sum_reg', 'regShifts', '-v7.3');
    fprintf('  saved %s\n', regFile);
    clear S
else
    fprintf('skip (exists): %s\n', regFile);
end
S = load(regFile); img = single(S.imgData_sum_reg); regShifts = S.regShifts; clear S
nVol = size(img, 3);
projMean = mean(img, 3);
regShiftXY = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), regShifts(:), 'UniformOutput', false));
badVol = any(abs(regShiftXY) > max_ok_shift_px, 2)';
fprintf('%d volumes; %d (%.1f%%) have a registration shift > %d px and are excluded\n', nVol, sum(badVol), 100*mean(badVol), max_ok_shift_px);

%% 3. automatic PB mask (cached), compared with the hand-drawn mask.mat if one exists
maskFile = fullfile(outDir, 'mask_auto.mat');
umPerPx = sp.um_per_px(1);
if ~isfile(maskFile) || overwrite_mask
    sigma_px       = max(mask_smooth_sigma_um / umPerPx, 0.5);
    open_radius_px = max(round(mask_open_radius_um / umPerPx), 1);
    join_line_px   = max(round(mask_join_line_um / umPerPx), 3);
    fprintf('mask params in px (um_per_px=%.4f): sigma=%.2f, open_radius=%d, join_line=%d\n', umPerPx, sigma_px, open_radius_px, join_line_px);

    % Summary image = coefficient of variation of the temporally smoothed
    % movie (std over time / mean), NOT the mean image. On this rig the mean
    % image is dominated by bright static puncta that outshine the bridge
    % (an 80th-percentile threshold on the mean grabbed one hemisphere plus
    % the puncta; Jaccard 0.26 vs the hand mask). The PB is where the bump
    % moves, so pixels with large SLOW fluctuations relative to their
    % brightness are PB: 2-s movmean kills shot noise, std/mean removes the
    % brightness bias. Tested against the hand-drawn mask for this trial:
    % mean 0.26-0.43, std 0.5, CV 0.71 (ugly_figures/exports/gain_change_example/mask_test*.png).
    projCV = std(movmean(img, round(2 * sp.fps), 3), 0, 3) ./ max(projMean, prctile(projMean, 20, 'all'));
    projSmooth  = imgaussfilt(projCV, sigma_px);
    projClipped = min(max(projSmooth, prctile(projSmooth, 5, 'all')), prctile(projSmooth, 98, 'all'));
    mask = projClipped > prctile(projClipped(:), mask_thresh_pct);
    mask = imopen(mask, strel('disk', open_radius_px));
    mask = imfill(mask, 'holes');
    areas = sort([regionprops(mask).Area], 'descend');
    if numel(areas) < 2 || (areas(1) / areas(2)) > 1.5
        mask = bwareafilt(mask, 1); maskCore = mask;
    else
        mask = bwareafilt(mask, 2); maskCore = mask;
        mask = imclose(mask, strel('line', join_line_px, 0)); % bridge the two hemispheres so the skeleton is one path
    end
    mask = imfill(mask, 'holes');
    if mask_dilate_px > 0
        mask = imdilate(mask, strel('disk', mask_dilate_px));
        maskCore = imdilate(maskCore, strel('disk', mask_dilate_px));
    end
    maskCore = maskCore & mask;
    save(maskFile, 'mask', 'maskCore', 'projCV');
    fprintf('saved %s\n', maskFile);
else
    M = load(maskFile); mask = M.mask; maskCore = M.maskCore; projCV = M.projCV; clear M
    fprintf('skip (exists): %s\n', maskFile);
end
handMask = [];
if isfile(fullfile(trialDir, 'mask.mat'))
    M = load(fullfile(trialDir, 'mask.mat'));
    if isfield(M, 'mask') && isequal(size(M.mask), size(mask)), handMask = logical(M.mask); end
    clear M
end

%% 4. glomeruli along the PB skeleton, and the QC figure (mean image | mask | glomeruli)
clusterIdx = pb_skeleton_glomeruli(mask, n_per_hemisphere);
nClusters  = 2 * n_per_hemisphere;
cen = zeros(nClusters, 2);
for c = 1:nClusters
    [cy, cx] = find(clusterIdx == c); cen(c,:) = [mean(cx), mean(cy)];
end

fig1 = figure(1); clf; set(fig1, 'Position', [50 50 1500 900], 'Color', 'w')
dispImg = min(max(projMean, prctile(projMean, 1, 'all')), prctile(projMean, 99.5, 'all'));
subplot(3,1,1); imagesc(dispImg); axis image; colormap(gca, bone); colorbar
title(sprintf('%s: mean of %d motion-corrected, z-summed volumes', trialName, nVol), 'Interpreter', 'none')
subplot(3,1,2); imagesc(imgaussfilt(projCV, 1)); axis image; colormap(gca, bone); hold on
contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
if any(mask(:) & ~maskCore(:)), contour(maskCore, [0.5 0.5], 'c', 'LineWidth', 1); end
leg = {'auto mask (red)'};
if ~isempty(handMask)
    contour(handMask, [0.5 0.5], 'y--', 'LineWidth', 1); leg{end+1} = 'hand-drawn mask.mat (yellow dashed)';
    fprintf('auto mask: %d px; hand mask: %d px; overlap (Jaccard) = %.2f\n', nnz(mask), nnz(handMask), nnz(mask & handMask) / nnz(mask | handMask));
end
title(['PB mask on the CV image (std of 2-s-smoothed movie / mean): ' strjoin(leg, ', ') sprintf(' -- %d px, threshold = %dth pct', nnz(mask), mask_thresh_pct)])
subplot(3,1,3)
cmap = [0 0 0; hsv(nClusters)];
image(ind2rgb(clusterIdx + 1, cmap)); axis image; hold on
plot(cen(:,1), cen(:,2), 'w-', 'LineWidth', 1); plot(cen(:,1), cen(:,2), 'wo', 'MarkerFaceColor', 'w', 'MarkerSize', 3)
for c = 1:nClusters
    text(cen(c,1), cen(c,2), num2str(c), 'Color', 'w', 'FontSize', 7, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
end
title(sprintf('%d glomeruli (%d per hemisphere) assigned by nearest point on the resampled skeleton', nClusters, n_per_hemisphere))
pb_export_png(fig1, fullfile(exportDir, 'A_mean_mask_glomeruli.png'), 150);

%% 5. fictrac: bar (cue) angle and fly heading on the imaging timebase
ftFile = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
if numel(ftFile) ~= 1, error('gain_change_example_fly:ft', 'expected one ficTracData_DAQ.mat in %s', trialDir); end
S = load(fullfile(ftFile(1).folder, ftFile(1).name)); ft = S.ftData_DAQ; clear S
t_ft    = seconds(ft.trialTime{1}(:))';
cue_ft  = ft.cuePos{1}(:)' / 192 * 2*pi;        % bar position on the arena, 0..2pi (process_ft convention)
cue_ft(abs(gradient(unwrap(cue_ft))) > 2) = NaN; % drop the frame-step artefacts process_ft also removes
head_ft = unwrap(ft.intHD{1}(:)');              % fictrac integrated heading (rad, unwrapped)
rspd_ft = ft.velYaw{1}(:)';
fspd_ft = ft.velFor{1}(:)';
% The DAQ starts fictrac and triggers ScanImage together, so imaging volume
% k is taken at (k-1)/fps s of trial time. There is no volume clock in this
% rig's DAQ table, so that is the alignment used; the cue-bump lag scan in
% section 7 would show a gross offset if it were wrong.
t_vol = (0:nVol-1) / sp.fps;
fprintf('imaging spans 0-%.1f s (%d volumes @ %.3f Hz), fictrac 0-%.1f s\n', t_vol(end), nVol, sp.fps, t_ft(end));
cue_vol  = angle(interp1(t_ft, exp(1i*fillmissing(cue_ft, 'linear')), t_vol, 'linear', 'extrap')); % (-pi, pi]
head_vol = interp1(t_ft, head_ft, t_vol, 'linear', 'extrap');                                     % unwrapped
rspd_vol = interp1(t_ft, smoothdata(rspd_ft, 'movmean', 10), t_vol, 'linear', 'extrap');
fspd_vol = interp1(t_ft, smoothdata(fspd_ft, 'movmean', 10), t_vol, 'linear', 'extrap');
% empirical closed-loop gain: slope of bar angle vs heading over 1-s bins where the fly turned
ib = 1:ft_rate:numel(t_ft); dc = diff(unwrap(fillmissing(cue_ft(ib), 'linear'))); dh = diff(head_ft(ib)); mv = abs(dh) > 0.15;
gain_emp = median(dc(mv) ./ dh(mv));
fprintf('empirical closed-loop gain (bar/heading) = %.2f (%d moving 1-s bins; negative = bar moves against the fly)\n', gain_emp, nnz(mv));

% bar frozen while the fly turns -> closed loop is not running, exclude
cue_u_vol = unwrap(cue_vol);
winV = max(3, round(bar_frozen_win_s * sp.fps));
barExcursion = movmax(cue_u_vol, winV) - movmin(cue_u_vol, winV);
barFrozen = barExcursion < bar_frozen_tol;
flyTurning = (movmax(head_vol, winV) - movmin(head_vol, winV)) > frozen_head_tol;
excl_vol = barFrozen & flyTurning;
% grow each excluded stretch by the window so its edges are not half-included
excl_vol = movmax(double(excl_vol), winV) > 0;
exclSegs = gce_segments(excl_vol, t_vol);
fprintf('closed-loop exclusion: %d volumes (%.1f%%) where the bar was frozen (< %.1f deg over %g s) while the heading moved > %.0f deg', ...
    nnz(excl_vol), 100*mean(excl_vol), rad2deg(bar_frozen_tol), bar_frozen_win_s, rad2deg(frozen_head_tol));
if isempty(exclSegs), fprintf('.\n'); else, fprintf(':\n'); fprintf('  %6.1f - %6.1f s\n', exclSegs'); end

%% 6. glomerulus traces -> z-score -> bump position (PVA)
img_2d = reshape(img, [], nVol);
bg  = mean(img_2d(~mask(:), :), 1);                 % per-volume background = mean outside the mask
f_cluster = nan(nClusters, nVol); nCorePx = zeros(nClusters, 1);
for c = 1:nClusters
    pix = clusterIdx(:) == c & maskCore(:);
    nCorePx(c) = nnz(pix);
    if nCorePx(c) > 0, f_cluster(c,:) = mean(img_2d(pix, :), 1) - bg; end
end
clear img_2d
emptyGlom = nCorePx < max(10, 0.2 * median(nCorePx(nCorePx > 0)));
f_cluster(emptyGlom, :) = NaN;
f_cluster(:, badVol) = NaN;
if any(emptyGlom), fprintf('%d glomeruli with (almost) no real-signal pixels -> NaN: %s\n', nnz(emptyGlom), mat2str(find(emptyGlom)')); end
z_cluster = (f_cluster - mean(f_cluster, 2, 'omitnan')) ./ std(f_cluster, 0, 2, 'omitnan');
alpha = repmat(linspace(-pi, pi - 2*pi/n_per_hemisphere, n_per_hemisphere), 1, 2); % same angle in both hemispheres

[x_tmp, y_tmp] = pol2cart(alpha, z_cluster');
[mu, rho] = cart2pol(mean(x_tmp, 2, 'omitnan'), mean(y_tmp, 2, 'omitnan'));
mu = mu(:)'; rho = rho(:)';

%% 7. bump vs bar and vs heading: sign, lag, offset
% the glomerulus numbering starts at whichever skeleton end bwmorph finds, so
% the bump may track the bar directly or mirrored -- pick the tighter one
ok = rho >= rho_thresh & ~isnan(mu) & ~excl_vol;
cs_same = gce_circ_std(mu(ok) - cue_vol(ok)); cs_flip = gce_circ_std(mu(ok) + cue_vol(ok));
if cs_flip < cs_same, bar_sign = -1; signNote = 'glomerulus order MIRRORED vs display (bump ~ -bar)';
else,                 bar_sign =  1; signNote = 'glomerulus order runs WITH the display (bump ~ +bar)'; end
cue_c = exp(1i * cue_vol); dt_vol = median(diff(t_vol));
lagScan = lag_range_sec(1):dt_vol:lag_range_sec(2); lagScore = nan(size(lagScan));
for L = 1:numel(lagScan)
    cueL = angle(interp1(t_vol, cue_c, t_vol - lagScan(L), 'linear', 'extrap'));
    lagScore(L) = gce_circ_std(mu(ok) - bar_sign * cueL(ok));
end
[~, iBest] = min(lagScore); lag_sec = max(lagScan(iBest), 0);
cue_lag = angle(interp1(t_vol, cue_c, t_vol - lag_sec, 'linear', 'extrap'));
offset  = angle(exp(1i * (mu - bar_sign * cue_lag)));
fprintf('%s; cue-bump lag = %.2f s; bump-bar offset circ std = %.1f deg (%.1f deg at lag 0); %.0f%% of volumes used (rho >= %.2f and closed loop running)\n', ...
    signNote, lag_sec, rad2deg(min(lagScore)), rad2deg(lagScore(find(lagScan >= 0, 1))), 100*mean(ok), rho_thresh);

% heading comparison: in closed loop bar = gain_emp*heading + const (gain_emp
% is signed, ~ -0.7: the bar moves against the fly). bar_sign*mu ~ bar, so
% unwrap bar_sign*mu where the bump is reliable and compare with gain_emp*heading.
mu_u = mu; mu_u(~ok) = NaN;
mu_unwrap = gce_unwrap_nan(bar_sign * mu_u);
head_scaled = gain_emp * (head_vol - head_vol(1));
% remove the arbitrary offset between the two by matching their medians
head_scaled = head_scaled - median(head_scaled(ok) - mu_unwrap(ok), 'omitnan');
bump_gain = gce_slope_bins(mu_unwrap, head_vol, ok, round(1/dt_vol));
fprintf('bump vs heading slope (1-s bins, fly turning) = %.2f; bar vs heading = %.2f\n', bump_gain, gain_emp);

bump = struct('trialDir', trialDir, 't_vol', t_vol, 'mu', mu, 'rho', rho, 'bar_sign', bar_sign, 'lag_sec', lag_sec, ...
    'cue_rad', cue_vol, 'cue_rad_lag', cue_lag, 'offset', offset, 'heading_rad_unwrapped', head_vol, 'r_speed', rspd_vol, 'f_speed', fspd_vol, ...
    'gain_empirical', gain_emp, 'z_cluster', z_cluster, 'f_cluster', f_cluster, 'clusterIdx', clusterIdx, 'alpha', alpha, ...
    'badVol', badVol, 'excl_vol', excl_vol, 'exclSegs', exclSegs, 'ok', ok, 'rho_thresh', rho_thresh, 'signNote', signNote); %#ok<NASGU>
save(fullfile(outDir, 'bump_results.mat'), 'bump');

%% 8. figure: heatmap, bump vs bar, bump vs heading
rwb = gce_redwhiteblue(256);
brk = @(x) gce_break_wraps(x);
fig2 = figure(2); clf; set(fig2, 'Position', [50 50 1600 1000], 'Color', 'w')
ax1 = subplot(5,1,1:2);
imagesc(t_vol, 1:nClusters, z_cluster, 'AlphaData', double(~isnan(z_cluster))); colormap(ax1, rwb); clim([-3 3]); set(ax1, 'Color', [0.7 0.7 0.7])
ylabel('glomerulus'); cb = colorbar(ax1); cb.Label.String = 'z-score';
title(sprintf('%s: z-scored glomerulus activity (closed loop, bar, empirical gain %.2f)', trialName, abs(gain_emp)), 'Interpreter', 'none')
ax2 = subplot(5,1,3); hold on
gce_shade(exclSegs, [-180 180])
plot(t_vol, brk(rad2deg(angle(exp(1i*bar_sign*cue_lag)))), '-', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5)
mu_plot = rad2deg(mu); mu_plot(~ok) = NaN;
plot(t_vol, brk(mu_plot), 'b.', 'MarkerSize', 5)
ylim([-180 180]); yticks(-180:90:180); ylabel('angle (deg)')
legend({sprintf('bar (signed, lagged %.2f s)', lag_sec), sprintf('bump PVA (rho >= %.1f)', rho_thresh)}, 'Location', 'northeastoutside')
title(sprintf('bump vs bar -- %s; offset circ std %.0f deg; gray = bar frozen while fly turning (excluded)', signNote, rad2deg(min(lagScore))))
ax3 = subplot(5,1,4); hold on
gce_shade(exclSegs, [-1 1] * 1e4)
plot(t_vol, rad2deg(head_scaled), '-', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.5)
plot(t_vol, rad2deg(mu_unwrap), 'b.', 'MarkerSize', 5)
ylabel('unwrapped angle (deg)'); ylim([min([mu_unwrap head_scaled], [], 'omitnan') max([mu_unwrap head_scaled], [], 'omitnan')] * 180/pi + [-30 30])
legend({sprintf('%.2f x fly heading (fictrac, shifted)', gain_emp), 'bump PVA (x bar\_sign), unwrapped'}, 'Location', 'northeastoutside')
title(sprintf('bump vs fly heading: bump/heading slope = %.2f (bar/heading = %.2f)', bump_gain, gain_emp))
ax4 = subplot(5,1,5); hold on
yyaxis left;  plot(t_vol, rho, 'k'); yline(rho_thresh, 'r:'); ylabel('bump strength (rho)'); set(gca, 'YColor', 'k')
yyaxis right; plot(t_vol, abs(rspd_vol), '-', 'Color', [0 0.6 0]); ylabel('|yaw speed| (rad/s)'); set(gca, 'YColor', [0 0.6 0])
xlabel('trial time (s)')
linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim([t_vol(1) t_vol(end)])
for ax = [ax1 ax4]; ax.Position([1 3]) = ax2.Position([1 3]); end
pb_export_png(fig2, fullfile(exportDir, 'B_heatmap_bump_vs_bar_heading.png'), 150);

% scatter: bump vs bar, and offset histogram
fig3 = figure(3); clf; set(fig3, 'Position', [50 50 1100 450], 'Color', 'w')
subplot(1,2,1)
scatter(rad2deg(angle(exp(1i*bar_sign*cue_lag(ok)))), rad2deg(mu(ok)), 6, t_vol(ok), 'filled'); axis equal; xlim([-180 180]); ylim([-180 180])
hold on; plot([-180 180], [-180 180], 'k:'); xlabel('bar (signed, lagged, deg)'); ylabel('bump PVA (deg)'); cb = colorbar; cb.Label.String = 'time (s)';
title(sprintf('bump vs bar, %d volumes with rho >= %.1f', nnz(ok), rho_thresh))
subplot(1,2,2)
histogram(rad2deg(offset(ok)), -180:10:180, 'FaceColor', [0.2 0.4 0.8]); xlim([-180 180]); xlabel('bump - bar (deg)'); ylabel('volumes')
title(sprintf('offset: circ mean %.0f deg, circ std %.0f deg', rad2deg(gce_circ_mean(offset(ok))), rad2deg(gce_circ_std(offset(ok)))))
pb_export_png(fig3, fullfile(exportDir, 'C_bump_vs_bar_scatter.png'), 150);
fprintf('done: figures in %s\n', exportDir);

%% local functions
function [imgSum, sp] = gce_read_raw_tif_summed(trialDir)
% BigTIFF IFD walker (copied from overshoot_bar_stim_script.m): MATLAB's Tiff
% class truncates ScanImage stacks at 65535 frames, so read frames by offset
% and sum the valid z-planes of each volume on the fly.
sp = jsondecode(fileread(fullfile(trialDir, 'scan_params.json')));
if sp.nchannels ~= 1, error('gce:multiChannel', '%s has %d channels; expected 1', trialDir, sp.nchannels); end
validPlanes = setdiff(1:sp.nplanes, sp.flyback_planes(:)' + 1);
tifList = dir(fullfile(trialDir, '*.tif'));
if numel(tifList) ~= 1, error('gce:tif', 'expected one .tif in %s, found %d', trialDir, numel(tifList)); end
tifPath = fullfile(tifList(1).folder, tifList(1).name);
framesPerVol = sp.nplanes; H = sp.px_height; W = sp.px_width;
d = dir(tifPath); fileSize = d.bytes;
nVolAlloc = ceil(fileSize / (H*W*2) / framesPerVol) + 10;
fid = fopen(tifPath, 'r', 'ieee-le'); cleanupFid = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid, 2, 'uint8=>char')'; version = fread(fid, 1, 'uint16');
if ~strcmp(magic, 'II') || version ~= 43, error('gce:notBigTiff', '%s is not a little-endian BigTIFF', tifPath); end
fread(fid, 1, 'uint16'); fread(fid, 1, 'uint16');
p = fread(fid, 1, 'uint64');
imgSum = zeros(H, W, nVolAlloc, 'single'); vol = zeros(H, W, framesPerVol, 'single');
k = 0; v = 0; f = 0; tic
while p > 0 && p + 8 <= fileSize
    fseek(fid, p, 'bof');
    nEntries = fread(fid, 1, 'uint64');
    if isempty(nEntries) || nEntries < 1 || nEntries > 500, break; end
    ent = fread(fid, [10 nEntries], 'uint16=>double');
    tags = ent(1,:); vals = ent(7,:) + ent(8,:)*2^16 + ent(9,:)*2^32 + ent(10,:)*2^48;
    nextIFD = fread(fid, 1, 'uint64');
    fseek(fid, vals(tags == 273), 'bof');
    frame = fread(fid, [W H], 'int16=>single');
    if numel(frame) ~= W*H, break; end
    k = k + 1; f = f + 1; vol(:,:,f) = frame';
    if f == framesPerVol
        v = v + 1; imgSum(:,:,v) = sum(vol(:,:,validPlanes), 3); f = 0;
        if mod(v, 500) == 0, fprintf('  read %d volumes (%.0f frames/s)\n', v, k/toc); end
    end
    p = nextIFD;
end
imgSum = imgSum(:,:,1:v);
fprintf('  %d frames -> %d complete volumes (%d leftover frames dropped)\n', k, v, k - v*framesPerVol);
sp.nVolumes = v;
end

function segs = gce_segments(flag, t)
% [start end] times of each run of true in a logical vector
d = diff([false flag(:)' false]); s = find(d == 1); e = find(d == -1) - 1;
segs = [t(s)' t(e)'];
end
function gce_shade(segs, yl)
for k = 1:size(segs, 1)
    patch([segs(k,1) segs(k,2) segs(k,2) segs(k,1)], [yl(1) yl(1) yl(2) yl(2)], [0.75 0.75 0.75], 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'HandleVisibility', 'off');
end
end
function m = gce_circ_mean(x), m = angle(mean(exp(1i * x(:)), 'omitnan')); end
function s = gce_circ_std(x)
x = x(~isnan(x)); if isempty(x), s = NaN; return; end
s = sqrt(-2 * log(max(abs(mean(exp(1i * x(:)))), eps)));
end
function y = gce_break_wraps(y), y(abs(diff([y y(end)])) > 180) = NaN; end
function u = gce_unwrap_nan(x)
% unwrap across NaN gaps by unwrapping the valid samples only
u = nan(size(x)); i = ~isnan(x); u(i) = unwrap(x(i));
end
function g = gce_slope_bins(y, x, ok, step)
% median ratio of 1-s changes in y and x, over bins where the fly turned
idx = 1:step:numel(x); dy = diff(y(idx)); dx = diff(x(idx)); good = ok(idx(1:end-1)) & ok(idx(2:end)) & abs(dx) > 0.15 & ~isnan(dy);
if nnz(good) < 5, g = NaN; else, g = median(dy(good) ./ dx(good)); end
end
function cmap = gce_redwhiteblue(n)
half = n/2; cmap = [[linspace(0,1,half)', linspace(0,1,half)', ones(half,1)]; [ones(half,1), linspace(1,0,half)', linspace(1,0,half)']];
end
