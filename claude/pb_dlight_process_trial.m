%% pb_dlight_process_trial
% Interactive, section-by-section pipeline for very low-SNR dLight trials
% (Z:\pablo\epg_dlight): raw ScanImage tif -> photon counts (kills line
% noise) -> z-sum -> motion correction on binned photons -> PB mask from
% turning correlation -> cross-validated low-rank denoising -> clean z-summed
% movie at full volume rate.
%
% Why photon counting instead of pb_remove_scan_noise: at these light levels
% the no-photon baseline is confined to about +/-50 counts and the scan/line
% noise (~13 cyc/line, ~11 counts) lives entirely inside it, while every
% photon is a single-pixel pulse of 50-1500 counts. Thresholding each frame
% at ~50 therefore removes the line noise in ONE frame, with no temporal
% averaging (see ugly_figures/scripts/fig_dlight_noise_anatomy.m). What is
% left is Poisson photon noise, which only pooling can beat -- that is the
% job of pb_lowrank_denoise (see its header for the cross-validation).
%
% Run one %% section at a time (Ctrl+Enter). Sections 2-3 take ~1 min each
% on a 600 s / 85k-frame trial over the network.

%% 0. parameters
expDir      = 'Z:\pablo\epg_dlight\20260216\fly 1\20260216-1_epg_dLight';
photonThr   = 50;     % counts; check in section 1 that the baseline never reaches it
outSubdir   = 'dlight_photon_001';
exportDir   = fullfile(fileparts(mfilename('fullpath')), '..', 'ugly_figures', 'exports');
addpath(fileparts(mfilename('fullpath')));

md = pb_scan_metadata(expDir);
[~, trialName] = fileparts(expDir);

%% 1. check the photon threshold on a short raw chunk
% The negative side of the histogram is pure baseline (photons are only
% positive), so mirroring it shows how far the baseline reaches. The fraction
% of baseline pixels above photonThr should be ~1e-5 or less.
raw = pb_load_raw_tif(expDir, 'volumes', [1 100]);
v = double(raw(:)); clear raw
fracBaselineAbove = mean(v < -photonThr) / mean(v < 0);
fprintf('baseline pixels beyond -%d: %.2g of negative pixels (want << 1e-3)\n', photonThr, fracBaselineAbove);
fprintf('photon pixels (> %d): %.2f%%, median pulse %.0f counts\n', photonThr, 100*mean(v > photonThr), median(v(v > photonThr)));
figure(Name = 'photon threshold check'); clf
histogram(v, -150:5:1500, 'EdgeColor', 'none'); set(gca, 'YScale', 'log'); hold on
histogram(-v(v < 0), -150:5:1500, 'DisplayStyle', 'stairs'); xline(photonThr, 'r--');
legend('raw', 'mirrored negative side (= baseline)', 'threshold'); xlabel('counts')
clear v

%% 2. stream the whole trial into photon-count z-sums (one pass, ~50 s)
% Reads by IFD offset, so >65535-frame BigTIFFs are not truncated.
o = pb_read_raw_tif_planes(md.tifPath, md.nplanes, {md.valid_planes}, 'photonThreshold', photonThr);
counts = o.stacks{1};   % Y x X x T, photons per pixel per volume (0..nplanes)
nV = size(counts, 3);
fprintf('%d volumes, %.0f photons/volume\n', nV, mean(sum(counts, [1 2])));

%% 3. behaviour: yaw velocity at each volume
f = dir(fullfile(expDir, '*ficTracData_DAQ.mat'));
ft = load(fullfile(f(1).folder, f(1).name)).ftData_DAQ;
tFt = seconds(ft.trialTime{1}); tVol = seconds(ft.volClock{1});
assert(numel(tVol) == nV, 'volClock has %d edges but the tif has %d volumes', numel(tVol), nV);
yaw = interp1(tFt, ft.velYaw{1}, tVol, 'linear', 'extrap');           % rad/s
heading = interp1(tFt, ft.intHD{1}, tVol, 'linear', 'extrap');        % unwrapped rad
fps = 1 / median(diff(tVol));

%% 4. motion correction on binned photons
[reg, shifts, regInfo] = pb_photon_register(counts);
figure(Name = 'photon registration shifts'); clf
tB = regInfo.binCentres / fps;
plot(tB, regInfo.binShifts, '.', 'Color', [0.7 0.7 0.7]); hold on
plot(tVol, shifts, 'LineWidth', 1.5); plot(tVol, abs(yaw)/5 - 4, 'k')
legend('dy (bin)', 'dx (bin)', 'dy', 'dx', '|yaw|/5 - 4'); xlabel('s'); ylabel('px')

%% 5. PB mask from turning correlation
rm = pb_rotation_map(reg, yaw, fps);
figure(Name = 'rotation map'); clf; tiledlayout(1, 3)
nexttile; imagesc(imgaussfilt(mean(reg, 3), 1.5)); axis image off; colorbar; title('registered photon mean')
nexttile; imagesc(rm.map); axis image off; colorbar; hold on; contour(rm.mask, [.5 .5], 'w'); title('corr with turning regressor')
nexttile; plot(rm.kernelLagsSec, rm.kernel, 'k'); xline(0, ':'); xlabel('lag (s)'); ylabel('photons per rad/s'); title('|yaw| \rightarrow PB kernel')
fprintf('PB mask: %d px; total PB photons vs regressor r = %.2f\n', nnz(rm.mask), ...
    corr(squeeze(sum(reg .* rm.mask, [1 2])), rm.regressor));

%% 6. low-rank denoising (thinning cross-validation picks rank and smoothing)
dn = pb_lowrank_denoise(counts, shifts, rm.mask);
clean = pb_lowrank_movie(dn);     % Y x X x T, denoised z-summed photon rate

%% 7. save (factors + raw counts; the full movie is rebuilt with pb_lowrank_movie)
outDir = fullfile(expDir, outSubdir);
if ~isfolder(outDir), mkdir(outDir); end
counts_u8 = uint8(counts); %#ok<NASGU>
save(fullfile(outDir, 'dlight_denoised.mat'), 'dn', 'rm', 'shifts', 'regInfo', 'counts_u8', ...
    'yaw', 'heading', 'tVol', 'photonThr', 'md', '-v7.3');
fprintf('saved %s\n', fullfile(outDir, 'dlight_denoised.mat'));

%% 8. validation figure: raw vs denoised at single-volume resolution
% PB columns: mask pixels binned by x (8 px) for a kymograph of the bridge.
xb = 1:8:size(reg, 2);
[~, xpix] = find(rm.mask);
lin = find(rm.mask);
colOf = discretize(xpix, [xb size(reg,2)+1]);
useCols = find(accumarray(colOf, 1, [numel(xb) 1]) >= 15);
kym = @(M) cell2mat(arrayfun(@(c) mean(M(lin(colOf == c), :), 1), useCols, 'UniformOutput', false));
Rm = reshape(imgaussfilt(reg, dn.sigma), [], nV); Cm = reshape(clean, [], nV);
kRaw = kym(Rm); kClean = kym(Cm);
pbRaw = sum(Rm(rm.mask(:), :), 1)'; pbClean = sum(Cm(rm.mask(:), :), 1)';
clear Rm Cm

% example window: the 60 s around the biggest turning bout
[~, iPk] = max(movmean(abs(yaw), round(5*fps)));
win = max(1, iPk - round(20*fps)):min(nV, iPk + round(40*fps));
[~, tEx] = max(abs(yaw(win))); tEx = win(tEx) + round(1*fps);   % 1 s after peak yaw

fig = figure(Name = 'dLight denoising', Color = 'w', Position = [40 40 1500 950]);
tl = tiledlayout(fig, 5, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('%s: photon counting + low-rank (k = %d, w = %d vols = %.2f s)', trialName, dn.rank, dn.smoothWin, dn.smoothWin/fps), 'Interpreter', 'none')
cl = [0 prctile(dn.meanImage(rm.mask), 99.5) * 1.6];
nexttile; imagesc(sum(reg(:,:,tEx), 3)); axis image off; title(sprintf('raw photons, 1 vol (t = %.1f s)', tVol(tEx))); colorbar
nexttile; imagesc(mean(imgaussfilt(reg(:,:,tEx-5:tEx+4), dn.sigma), 3)); axis image off; clim(cl); title('raw, 10-vol mean (~1 s)'); colorbar
nexttile; imagesc(clean(:,:,tEx)); axis image off; clim(cl); title('denoised, 1 vol'); colorbar
nexttile; imagesc(rm.map); axis image off; hold on; contour(rm.mask, [.5 .5], 'w'); title('turning map + mask'); colorbar
nexttile([1 4]); plot(tVol(win), yaw(win), 'k'); ylabel('yaw (rad/s)'); xlim(tVol(win([1 end]))); box off
nexttile([1 4]); imagesc(tVol(win), 1:numel(useCols), movmean(kRaw(:, win), 3, 2) ./ mean(kRaw, 2) - 1); clim([-0.5 0.8]); ylabel('PB column'); title('raw (3-vol mean), \DeltaF/F_{mean}'); colorbar
nexttile([1 4]); imagesc(tVol(win), 1:numel(useCols), kClean(:, win) ./ mean(kClean, 2) - 1); clim([-0.5 0.8]); ylabel('PB column'); title('denoised, \DeltaF/F_{mean}'); colorbar
nexttile([1 3]); plot(tVol, pbRaw / mean(pbRaw), 'Color', [0.75 0.75 0.75]); hold on
plot(tVol, pbClean / mean(pbClean), 'b'); plot(tVol, 0.6 + abs(yaw)/max(abs(yaw)) * 0.3, 'k'); xlim(tVol([1 end])); box off
legend('raw PB total', 'denoised PB total', '|yaw| (scaled)', 'Location', 'northwest', 'Box', 'off'); xlabel('s')
nexttile; G = reshape(dn.cv.heldoutLLgain, numel(unique(dn.cv.rank)), []);
imagesc(unique(dn.cv.smoothWin), unique(dn.cv.rank), G); colorbar; xlabel('w (vols)'); ylabel('rank'); title('held-out LL gain (nats/vol)')

if ~isfolder(exportDir), mkdir(exportDir); end
exportgraphics(fig, fullfile(exportDir, sprintf('dlight_denoise_%s.png', trialName)), 'Resolution', 120);
