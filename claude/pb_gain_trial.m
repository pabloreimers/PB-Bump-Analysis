function out = pb_gain_trial(trialDir, tifIdx, opts)
% PB_GAIN_TRIAL  Raw ScanImage tif -> bump (PVA) + fictrac for one gain-change trial.
%
%   out = pb_gain_trial(trialDir)             one tif in the folder
%   out = pb_gain_trial(trialDir, k)          k-th tif (trial_00k) in a folder that holds several
%   out = pb_gain_trial(trialDir, k, opts)    override parameters (see defaults below)
%
% Function form of data scripts/gain_change_example_fly.m (developed on
% 20230110 fly 1 trial 1), so a batch script can run it over many flies:
%   1. raw BigTIFF read by walking the IFD chain (MATLAB's Tiff class dies at
%      65535 frames), summing the valid z-planes of each volume
%   2. rigid motion correction (pb_register)
%   3. automatic PB mask from the COEFFICIENT OF VARIATION image (std of the
%      2-s-smoothed movie / mean). The mean image is dominated by static
%      bright puncta on this rig; the PB is where slow bump fluctuations
%      live. 75th-percentile threshold + 1 px dilation (Jaccard 0.70 vs the
%      hand mask on the development trial; 0.26 for a mean-image threshold).
%   4. 2*nPerHemisphere glomerulus bins along the mask skeleton
%      (pb_skeleton_glomeruli), alpha repeated over the two hemispheres
%   5. per-bin mean F minus the per-volume mean outside the mask, z-scored
%      over the trial, bump position/strength as the population vector (mu, rho)
%   6. fictrac bar (cue) angle, heading and speeds on the imaging timebase.
%      Imaging volume k is placed at (k-1)/fps s of trial time (the DAQ starts
%      fictrac and triggers ScanImage together; this rig logs no volume clock).
%      Source: <trial>\*ficTracData_DAQ.mat row tifIdx; if that file is
%      missing (never preprocessed) the raw daqData_*trial_00k.mat is used:
%      g4panels and ficTracYaw are 0-10 V per turn, ficTracIntForward is the
%      integrated forward displacement on the same 0-10 V scale.
%   7. closed-loop sanity flag: volumes where the bar did not move (< 2.5 deg
%      over 3 s) although the heading moved > 8 deg (bar parked / closed loop
%      not running) -> out.excl_vol. Bump-vs-bar statistics (sign of the
%      glomerulus order, cue-bump lag, offset) use only volumes with
%      rho >= rho_thresh and ~excl_vol.
%
% Cached in <trialDir>\claude_pipeline\: imgData_sum_reg_trial_00k.mat (the
% registered movie + shifts), mask_auto_trial_00k.mat. The hand-drawn mask.mat
% and any registration_00k\ folder in the trial are never touched. Set
% opts.overwrite to redo.
%
% out fields (lab all_data-style):
%   .ft   xf (fictrac time), xb (volume time), cue (-pi..pi, bar on arena),
%         heading (unwrapped rad), r_speed, f_speed, cue_src ('ficTracData_DAQ'|'daqData'),
%         gain_empirical (signed slope bar/heading), dark (pattern is background)
%   .im   f (nGlom x T, bg-subtracted F), z, mu, rho, alpha, mask, maskCore,
%         clusterIdx, projMean, projCV, nCorePx, badVol
%   .bump bar_sign, signNote, lag_sec, cue_lag, offset, ok, excl_vol, exclSegs,
%         offset_circstd_deg, offset_circmean_deg
%   .meta trialDir, tifPath, tifIdx, fps, nVol, um_per_px, mask_jaccard_vs_hand (NaN if no mask.mat)

arguments
    trialDir (1,:) char
    tifIdx (1,1) double = 1
    opts.nPerHemisphere (1,1) double = 16
    opts.maxOkShiftPx (1,1) double = 10
    opts.maskSigmaUm (1,1) double = 1.2
    opts.maskOpenRadiusUm (1,1) double = 1
    opts.maskJoinLineUm (1,1) double = 15
    opts.maskThreshPct (1,1) double = 75
    opts.maskDilatePx (1,1) double = 1
    opts.rhoThresh (1,1) double = 0.1
    opts.lagRangeSec (1,2) double = [-2 5]
    opts.barFrozenWinS (1,1) double = 3
    opts.barFrozenTolRad (1,1) double = deg2rad(2.5)
    opts.frozenHeadTolRad (1,1) double = deg2rad(8)
    opts.overwrite (1,1) logical = false
    opts.verbose (1,1) logical = true
    opts.refMask = []   % struct('mask',..,'maskCore',..,'projMean',..) from another trial of the SAME fly:
                        % that mask is shifted onto this trial (rigid xcorr of the mean images) instead
                        % of being recomputed. Needed because the CV image has no contrast when the bump
                        % is weak (dark trials): e.g. 20221221 fly 1 trial 2 gave a random blob (Jaccard 0.13).
    opts.refMaxShiftPx (1,1) double = 15
    opts.barSign (1,1) double = -1 % FIXED convention, not estimated per trial. With bins numbered left->right in
                                    % the image, the PVA angle runs opposite to the bar: bump ~ -bar + const. Measured
                                    % on all 32 trials of the 2026-10-07 epg_7f batch (31/32 agree; the one exception
                                    % is a dark trial with 2% walking). A bump moving against the bar is physiology,
                                    % not a display artefact, so it is NOT compensated; the per-trial slope of bump
                                    % vs heading is still saved as a diagnostic (bump.slope_mu_vs_heading).
    opts.preferHandMask (1,1) logical = true % use <trial>\mask.mat (hand-drawn) when it exists; the auto CV
                                             % mask failed on 2 of the first 5 reference flies of the gain_change
                                             % batch (one hemisphere only), so the hand mask is the safer default
end

say = @(varargin) fprintf([varargin{1} '\n'], varargin{2:end});
if ~opts.verbose, say = @(varargin) []; end

outDir = fullfile(trialDir, 'claude_pipeline');
if ~isfolder(outDir), mkdir(outDir); end
sp = jsondecode(fileread(fullfile(trialDir, 'scan_params.json')));
tifList = dir(fullfile(trialDir, '*.tif'));
tifList = tifList(~cellfun(@isempty, regexp({tifList.name}, '_trial_\d+_\d+\.tif$', 'once')));
if isempty(tifList), error('pb_gain_trial:noTif', 'no ScanImage tif in %s', trialDir); end
tifTrial = cellfun(@(n) str2double(regexp(n, '_trial_(\d+)_', 'tokens', 'once')), {tifList.name});
[tifTrial, ord] = sort(tifTrial); tifList = tifList(ord);
k = find(tifTrial == tifIdx, 1);
if isempty(k), error('pb_gain_trial:tifIdx', '%s has no tif for trial %d (has %s)', trialDir, tifIdx, mat2str(tifTrial)); end
tifPath = fullfile(tifList(k).folder, tifList(k).name);
tag = sprintf('trial_%03d', tifIdx);

%% 1-2. read + register (cached)
regFile = fullfile(outDir, ['imgData_sum_reg_' tag '.mat']);
if ~isfile(regFile) || opts.overwrite
    say('  reading %s', tifList(k).name);
    imgData_sum = local_read_tif_summed(tifPath, sp);
    say('  registering %d volumes', size(imgData_sum, 3));
    [imgData_sum_reg, regShifts, regDiag] = pb_register(imgData_sum, 'smoothWindow', 5, 'maxShift', [30 45], 'verbose', false, 'figureVisible', false); %#ok<ASGLU>
    close(regDiag.figHandles);
    save(regFile, 'imgData_sum_reg', 'regShifts', '-v7.3');
    clear imgData_sum
else
    S = load(regFile); imgData_sum_reg = S.imgData_sum_reg; regShifts = S.regShifts; clear S
end
img = single(imgData_sum_reg); clear imgData_sum_reg
nVol = size(img, 3);
regShiftXY = cell2mat(arrayfun(@(s) reshape(s.shifts, 1, []), regShifts(:), 'UniformOutput', false));
% bad volumes: registration ran away, or the frame went dark (laser/shutter
% dropout -- 20221221 fly 3 trial 1 has ~100 volumes near zero intensity at
% 45-75 s, during which the shifts also pinned to the +/-45 px bound)
frameMean = squeeze(mean(img, [1 2]))';
badVol = any(abs(regShiftXY) > opts.maxOkShiftPx, 2)' | frameMean < 0.5 * median(frameMean);
goodVol = ~badVol;
projMean = mean(img(:,:,goodVol), 3);
% rows/cols that registration ever shifted into the frame are zero-filled on
% some volumes -> huge coefficient of variation at the border. Mask them out
% of the CV image (margin = 99th pct of the good-volume shifts + 2 px).
shGood = abs(regShiftXY(goodVol, :));
border = min(ceil(prctile(shGood(:), 99)) + 2, 20);
borderMask = false(size(projMean)); borderMask([1:border, end-border+1:end], :) = true; borderMask(:, [1:border, end-border+1:end]) = true;

%% 3. mask (cached) -- own CV mask, or the fly's reference mask shifted onto this trial
umPerPx = sp.um_per_px(1);
refShift = [0 0];
handFile = fullfile(trialDir, 'mask.mat');
useHand = false;
if opts.preferHandMask && isfile(handFile)
    M = load(handFile);
    if isfield(M, 'mask') && isequal(size(M.mask), size(projMean)) && local_mask_ok(logical(M.mask))
        useHand = true; mask = logical(M.mask); maskCore = mask;
    end
    clear M
end
if useHand
    maskFile = fullfile(outDir, ['mask_hand_' tag '.mat']);
    projCV = local_cv_image(img, goodVol, sp.fps, projMean, borderMask);
    mask_source = 'hand mask.mat';
    save(maskFile, 'mask', 'maskCore', 'projCV', 'mask_source');
elseif ~isempty(opts.refMask)
    maskFile = fullfile(outDir, ['mask_ref_' tag '.mat']);
    R = opts.refMask;
    % integer rigid shift that best aligns this trial's mean image to the
    % reference trial's (smoothed, mean-removed, fft cross-correlation)
    a = imgaussfilt(projMean, 2); a = a - mean(a(:));
    b = imgaussfilt(R.projMean, 2); b = b - mean(b(:));
    xc = fftshift(real(ifft2(fft2(a) .* conj(fft2(b)))));
    [H, W] = size(a); cy = floor(H/2) + 1; cx = floor(W/2) + 1;
    win = false(H, W); win(max(cy-opts.refMaxShiftPx,1):min(cy+opts.refMaxShiftPx,H), max(cx-opts.refMaxShiftPx,1):min(cx+opts.refMaxShiftPx,W)) = true;
    xc(~win) = -Inf; [~, iMax] = max(xc(:)); [py, px] = ind2sub([H W], iMax);
    refShift = [py - cy, px - cx]; % [rows cols] to move the reference onto this trial
    mask = imtranslate(R.mask, fliplr(refShift)) > 0.5;
    maskCore = imtranslate(R.maskCore, fliplr(refShift)) > 0.5;
    projCV = local_cv_image(img, goodVol, sp.fps, projMean, borderMask);
    mask_source = 'reference trial mask, shifted';
    save(maskFile, 'mask', 'maskCore', 'projCV', 'refShift');
    say('  reference mask shifted by [%d %d] px (row, col)', refShift(1), refShift(2));
else
maskFile = fullfile(outDir, ['mask_auto2_' tag '.mat']); % "2": good-volume-only CV image with border mask (v1 files are stale)
if ~isfile(maskFile) || opts.overwrite
    sigma_px = max(opts.maskSigmaUm / umPerPx, 0.5);
    open_px  = max(round(opts.maskOpenRadiusUm / umPerPx), 1);
    join_px  = max(round(opts.maskJoinLineUm / umPerPx), 3);
    projCV = local_cv_image(img, goodVol, sp.fps, projMean, borderMask);
    projS = imgaussfilt(projCV, sigma_px);
    inner = projS(~borderMask); % percentiles over the non-border pixels, else the zeroed border drags the threshold down
    projC = min(max(projS, prctile(inner, 5)), prctile(inner, 98));
    mask = projC > prctile(projC(~borderMask), opts.maskThreshPct);
    mask = imopen(mask, strel('disk', open_px)); mask = imfill(mask, 'holes');
    areas = sort([regionprops(mask).Area], 'descend');
    if numel(areas) < 2 || areas(1) / areas(2) > 1.5
        mask = bwareafilt(mask, 1); maskCore = mask;
    else
        mask = bwareafilt(mask, 2); maskCore = mask;
        mask = imclose(mask, strel('line', join_px, 0));
    end
    mask = imfill(mask, 'holes');
    if opts.maskDilatePx > 0
        mask = imdilate(mask, strel('disk', opts.maskDilatePx)); maskCore = imdilate(maskCore, strel('disk', opts.maskDilatePx));
    end
    maskCore = maskCore & mask;
    mask_source = 'own CV mask';
    % sanity: the PB arch spans most of the frame width AND is not a thin
    % strip (height >= 15% of the frame). If it fails and a hand-drawn
    % mask.mat exists in the trial folder, use that instead.
    if ~local_mask_ok(mask)
        if isfile(fullfile(trialDir, 'mask.mat'))
            M = load(fullfile(trialDir, 'mask.mat'));
            if isfield(M, 'mask') && isequal(size(M.mask), size(mask))
                mask = logical(M.mask); maskCore = mask; mask_source = 'hand mask.mat (auto mask failed sanity)';
                say('  auto mask failed sanity -> using the hand-drawn mask.mat');
            end
            clear M
        else
            say('  WARNING: auto mask failed sanity and there is no hand mask.mat');
        end
    end
    save(maskFile, 'mask', 'maskCore', 'projCV', 'mask_source');
else
    M = load(maskFile); mask = M.mask; maskCore = M.maskCore; projCV = M.projCV; mask_source = M.mask_source; clear M
end
end
[ym, xm] = find(mask); mask_width_frac = range(xm) / size(mask, 2); mask_height_frac = range(ym) / size(mask, 1); %#ok<ASGLU>
jac = NaN; % overlap of the mask in use with the hand-drawn one (1 when the hand mask itself is used)
if isfile(handFile)
    M = load(handFile);
    if isfield(M, 'mask') && isequal(size(M.mask), size(mask)), h = logical(M.mask); jac = nnz(mask & h) / nnz(mask | h); end
    clear M
end

%% 4-5. glomeruli, traces, PVA
clusterIdx = pb_skeleton_glomeruli(mask, opts.nPerHemisphere);
nClusters = 2 * opts.nPerHemisphere;
% number the bins left -> right in the image, always. pb_skeleton_glomeruli
% starts at whichever skeleton endpoint bwmorph finds first, which is not
% guaranteed; a consistent numbering is what makes the bump/bar relation a
% property of the fly rather than of the mask.
[~, x1] = find(clusterIdx == 1); [~, xN] = find(clusterIdx == nClusters);
if mean(x1) > mean(xN)
    tmp = clusterIdx; tmp(clusterIdx > 0) = nClusters + 1 - clusterIdx(clusterIdx > 0); clusterIdx = tmp; clear tmp
end
img_2d = reshape(img, [], nVol); clear img
bg = mean(img_2d(~mask(:), :), 1);
f_cluster = nan(nClusters, nVol); nCorePx = zeros(nClusters, 1);
for c = 1:nClusters
    pix = clusterIdx(:) == c & maskCore(:); nCorePx(c) = nnz(pix);
    if nCorePx(c) > 0, f_cluster(c,:) = mean(img_2d(pix, :), 1) - bg; end
end
clear img_2d
emptyGlom = nCorePx < max(10, 0.2 * median(nCorePx(nCorePx > 0)));
f_cluster(emptyGlom, :) = NaN; f_cluster(:, badVol) = NaN;
z_cluster = (f_cluster - mean(f_cluster, 2, 'omitnan')) ./ std(f_cluster, 0, 2, 'omitnan');
alpha = repmat(linspace(-pi, pi - 2*pi/opts.nPerHemisphere, opts.nPerHemisphere), 1, 2);
[x_tmp, y_tmp] = pol2cart(alpha, z_cluster');
[mu, rho] = cart2pol(mean(x_tmp, 2, 'omitnan'), mean(y_tmp, 2, 'omitnan'));
mu = mu(:)'; rho = rho(:)';

%% 6. fictrac
t_vol = (0:nVol-1) / sp.fps;
ft = local_load_ft(trialDir, tifIdx);
cue_vol  = angle(interp1(ft.t, exp(1i*ft.cue), t_vol, 'linear', 'extrap'));
head_vol = interp1(ft.t, ft.heading, t_vol, 'linear', 'extrap');
rspd_vol = interp1(ft.t, smoothdata(ft.r_speed, 'movmean', 10), t_vol, 'linear', 'extrap');
fspd_vol = interp1(ft.t, smoothdata(ft.f_speed, 'movmean', 10), t_vol, 'linear', 'extrap');
ib = 1:round(ft.rate):numel(ft.t); dc = diff(unwrap(ft.cue(ib))); dh = diff(ft.heading(ib)); mv = abs(dh) > 0.15;
if nnz(mv) >= 5, gain_emp = median(dc(mv) ./ dh(mv)); else, gain_emp = NaN; end

%% 7. closed-loop sanity flag + bump vs bar
winV = max(3, round(opts.barFrozenWinS * sp.fps));
cue_u = unwrap(cue_vol);
barFrozen  = (movmax(cue_u, winV) - movmin(cue_u, winV)) < opts.barFrozenTolRad;
flyTurning = (movmax(head_vol, winV) - movmin(head_vol, winV)) > opts.frozenHeadTolRad;
excl_vol = movmax(double(barFrozen & flyTurning), winV) > 0;
exclSegs = local_segments(excl_vol, t_vol);

ok = rho >= opts.rhoThresh & ~isnan(mu) & ~excl_vol;
% Diagnostics of the bump/bar relation. The sign itself is NOT estimated per
% trial any more (see opts.barSign): with left->right numbering the bump
% runs as -bar in 31/32 trials of this dataset, and a bump moving against
% the bar is physiology, not something to compensate. Saved per trial:
% slope of unwrapped bump vs heading (1-s bins where the fly turned), and
% the circular spread of bump-bar and bump+bar.
mu_tmp = mu; mu_tmp(~ok) = NaN; iOk = ~isnan(mu_tmp); mu_u = nan(size(mu)); mu_u(iOk) = unwrap(mu_tmp(iOk));
stepV = max(1, round(1 / median(diff(t_vol))));
ib2 = 1:stepV:nVol; dmu = diff(mu_u(ib2)); dhd = diff(head_vol(ib2));
good = ok(ib2(1:end-1)) & ok(ib2(2:end)) & abs(dhd) > 0.15 & ~isnan(dmu);
if nnz(good) >= 5, slope_mu_head = median(dmu(good) ./ dhd(good)); else, slope_mu_head = NaN; end
cs_same = local_circ_std(mu(ok) - cue_vol(ok)); cs_flip = local_circ_std(mu(ok) + cue_vol(ok));
bar_sign = opts.barSign;
if bar_sign < 0, signNote = 'fixed convention: bump ~ -bar'; else, signNote = 'fixed convention: bump ~ +bar'; end
measuredSign = sign(slope_mu_head) * sign(gain_emp);
if isfinite(measuredSign) && measuredSign ~= bar_sign
    signNote = [signNote sprintf(' (NB measured bump/heading slope %.2f disagrees, %d turning bins)', slope_mu_head, nnz(good))];
end
cue_c = exp(1i * cue_vol); dt = median(diff(t_vol));
lagScan = opts.lagRangeSec(1):dt:opts.lagRangeSec(2); lagScore = nan(size(lagScan));
for L = 1:numel(lagScan)
    cueL = angle(interp1(t_vol, cue_c, t_vol - lagScan(L), 'linear', 'extrap'));
    lagScore(L) = local_circ_std(mu(ok) - bar_sign * cueL(ok));
end
[~, iBest] = min(lagScore); lag_sec = max(lagScan(iBest), 0);
if isnan(lag_sec), lag_sec = 0; end
cue_lag = angle(interp1(t_vol, cue_c, t_vol - lag_sec, 'linear', 'extrap'));
offset = angle(exp(1i * (mu - bar_sign * cue_lag)));

say('  %d vol @ %.2f Hz (%d bad) | %s: %d px, width %.0f%% (Jaccard vs hand %.2f) | gain %.2f | bump/heading slope %.2f | %s | lag %.2f s | offset std %.0f deg | %.0f%% vol used | excl %.0f%%', ...
    nVol, sp.fps, nnz(badVol), mask_source, nnz(mask), 100*mask_width_frac, jac, gain_emp, slope_mu_head, signNote, lag_sec, rad2deg(local_circ_std(offset(ok))), 100*mean(ok), 100*mean(excl_vol));

%% pack
out.ft = struct('xf', ft.t, 'xb', t_vol, 'cue', cue_vol, 'heading', head_vol, 'r_speed', rspd_vol, 'f_speed', fspd_vol, ...
    'cue_src', ft.src, 'gain_empirical', gain_emp, 'dark', ft.dark, 'pattern', ft.pattern, ...
    'cue_raw', ft.cue, 'heading_raw', ft.heading, 'heading_raw_uncleaned', ft.heading_uncleaned, 'heading_clean_info', ft.heading_clean_info, ...
    'r_speed_raw', ft.r_speed, 'f_speed_raw', ft.f_speed);
out.im = struct('f', f_cluster, 'z', z_cluster, 'mu', mu, 'rho', rho, 'alpha', alpha, 'mask', mask, 'maskCore', maskCore, ...
    'clusterIdx', clusterIdx, 'projMean', projMean, 'projCV', projCV, 'nCorePx', nCorePx, 'badVol', badVol);
out.bump = struct('bar_sign', bar_sign, 'signNote', signNote, 'slope_mu_vs_heading', slope_mu_head, 'n_turning_bins', nnz(good), ...
    'bar_offset_circstd_same', cs_same, 'bar_offset_circstd_flip', cs_flip, 'lag_sec', lag_sec, 'cue_lag', cue_lag, 'offset', offset, ...
    'ok', ok, 'excl_vol', excl_vol, 'exclSegs', exclSegs, 'offset_circstd_deg', rad2deg(local_circ_std(offset(ok))), ...
    'offset_circmean_deg', rad2deg(angle(mean(exp(1i*offset(ok)), 'omitnan'))), 'rho_thresh', opts.rhoThresh);
out.meta = struct('trialDir', trialDir, 'tifPath', tifPath, 'tifIdx', tifIdx, 'fps', sp.fps, 'nVol', nVol, ...
    'um_per_px', umPerPx, 'mask_jaccard_vs_hand', jac, 'mask_width_frac', mask_width_frac, ...
    'mask_source', mask_source, 'refShift', refShift, 'border_px', border, 'n_badVol', nnz(badVol), ...
    'opts', rmfield(opts, 'refMask'), 'processed', datestr(now, 'yyyy-mm-dd HH:MM'));
end

%% ------------------------------------------------------------------------
function imgSum = local_read_tif_summed(tifPath, sp)
validPlanes = setdiff(1:sp.nplanes, sp.flyback_planes(:)' + 1);
framesPerVol = sp.nplanes * sp.nchannels; H = sp.px_height; W = sp.px_width;
if sp.nchannels ~= 1, error('pb_gain_trial:channels', '%s: %d channels, expected 1', tifPath, sp.nchannels); end
d = dir(tifPath); fileSize = d.bytes;
nVolAlloc = ceil(fileSize / (H*W*2) / framesPerVol) + 10;
fid = fopen(tifPath, 'r', 'ieee-le'); c = onCleanup(@() fclose(fid)); %#ok<NASGU>
magic = fread(fid, 2, 'uint8=>char')'; version = fread(fid, 1, 'uint16');
if ~strcmp(magic, 'II') || version ~= 43, error('pb_gain_trial:notBigTiff', '%s is not a little-endian BigTIFF', tifPath); end
fread(fid, 2, 'uint16'); p = fread(fid, 1, 'uint64');
imgSum = zeros(H, W, nVolAlloc, 'single'); vol = zeros(H, W, framesPerVol, 'single');
k = 0; v = 0; f = 0;
while p > 0 && p + 8 <= fileSize
    fseek(fid, p, 'bof'); nEntries = fread(fid, 1, 'uint64');
    if isempty(nEntries) || nEntries < 1 || nEntries > 500, break; end
    ent = fread(fid, [10 nEntries], 'uint16=>double');
    tags = ent(1,:); vals = ent(7,:) + ent(8,:)*2^16 + ent(9,:)*2^32 + ent(10,:)*2^48;
    nextIFD = fread(fid, 1, 'uint64');
    fseek(fid, vals(tags == 273), 'bof'); frame = fread(fid, [W H], 'int16=>single');
    if numel(frame) ~= W*H, break; end
    k = k + 1; f = f + 1; vol(:,:,f) = frame';
    if f == framesPerVol, v = v + 1; imgSum(:,:,v) = sum(vol(:,:,validPlanes), 3); f = 0; end
    p = nextIFD;
end
if v == 0, error('pb_gain_trial:noFrames', 'no complete volume read from %s', tifPath); end
imgSum = imgSum(:,:,1:v);
end

function ft = local_load_ft(trialDir, tifIdx)
% bar/heading/speeds at the fictrac/DAQ-table rate; from ficTracData_DAQ.mat if present, else raw daqData
ft.pattern = ''; ft.dark = false;
ts = fullfile(trialDir, 'csv', 'trialSettings.csv');
if isfile(ts)
    try, T = readtable(ts, 'Delimiter', ','); ft.pattern = char(string(T.patternPath(end))); catch, end
end
ft.dark = contains(ft.pattern, 'background');
f = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
% ficTracData_DAQ.cuePos is wrong for pattern 0006_4px_brightbar_with_starfield
% (apparent gain 3-5; the DAQ g4panels channel gives the true 0.7), so that
% pattern always goes through the raw daqData path.
if contains(ft.pattern, 'with_starfield'), f = []; end
if numel(f) == 1
    S = load(fullfile(f(1).folder, f(1).name)); T = S.ftData_DAQ; clear S
    row = find(T.trialNum == tifIdx, 1); if isempty(row), row = min(tifIdx, height(T)); end
    ft.t = seconds(T.trialTime{row}(:))';
    cue = T.cuePos{row}(:)' / 192 * 2*pi;
    cue(abs(gradient(unwrap(cue))) > 2) = NaN;
    ft.cue = angle(exp(1i * fillmissing(cue, 'linear', 'EndValues', 'nearest')));
    ft.heading_uncleaned = unwrap(T.intHD{row}(:)');
    ft.r_speed = T.velYaw{row}(:)';
    ft.f_speed = T.velFor{row}(:)';
    ft.rate = 1 / median(diff(ft.t));
    ft.src = 'ficTracData_DAQ';
else
    d = dir(fullfile(trialDir, sprintf('*daqData*trial_%03d.mat', tifIdx)));
    if numel(d) ~= 1, error('pb_gain_trial:noFt', '%s: no ficTracData_DAQ.mat and no daqData for trial %d', trialDir, tifIdx); end
    S = load(fullfile(d(1).folder, d(1).name), 'trialData'); D = S.trialData; clear S
    fs = D.Properties.SampleRate; if isempty(fs) || isnan(fs), fs = 1/seconds(median(diff(D.Properties.RowTimes))); end
    vn = D.Properties.VariableNames;
    pan = D.(vn{find(contains(lower(vn), 'panel'), 1)}); yaw = D.(vn{find(contains(lower(vn), 'yaw'), 1)}); fwd = D.(vn{find(contains(lower(vn), 'forward'), 1)});
    ds = round(fs / 60); % to ~60 Hz like the fictrac table
    % unwrap the 0-10 V (= one turn) channels BEFORE any filtering: a median
    % filter across the 10 V -> 0 V wrap manufactures intermediate values
    % that unwrap into spurious heading blips
    yawU = unwrap(yaw(:)' / 10 * 2*pi); fwdU = unwrap(fwd(:)' / 10 * 2*pi); panU = unwrap(pan(:)' / 10 * 2*pi);
    yawU = movmedian(yawU, ds); fwdU = movmedian(fwdU, ds); panU = movmedian(panU, ds);
    idx = 1:ds:numel(yawU);
    ft.t = (idx - 1) / fs;
    ft.cue = angle(exp(1i * panU(idx)));
    ft.heading_uncleaned = yawU(idx);
    ft.rate = fs / ds;
    ft.f_speed = gradient(fwdU(idx)) * ft.rate; % integrated-forward channel; units = (ball rad)/s on the 0-10 V = 2 pi convention
    ft.src = 'daqData';
end
% single-frame fictrac glitches (spikes while locking on, steps after losing
% the ball): remove them from the integrated heading, see pb_clean_heading
[ft.heading, ft.heading_clean_info] = pb_clean_heading(ft.heading_uncleaned);
if strcmp(ft.src, 'daqData'), ft.r_speed = gradient(ft.heading) * ft.rate; end
end

function cv = local_cv_image(img, goodVol, fps, projMean, borderMask)
% coefficient of variation of the 2-s-smoothed movie over the good volumes,
% zeroed on the registration border
g = img(:,:,goodVol);
cv = std(movmean(g, round(2 * fps), 3), 0, 3) ./ max(projMean, prctile(projMean, 20, 'all'));
cv(borderMask) = 0;
end

function ok = local_mask_ok(mask)
[ym, xm] = find(mask);
% width >= 65% of the frame: a one-hemisphere blob (20230106 fly 2, 53%) must fail
ok = ~isempty(xm) && range(xm) / size(mask, 2) >= 0.65 && range(ym) / size(mask, 1) >= 0.15 && nnz(mask) >= 0.04 * numel(mask);
end

function segs = local_segments(flag, t)
d = diff([false flag(:)' false]); s = find(d == 1); e = find(d == -1) - 1; segs = [t(s)' t(e)'];
end
function s = local_circ_std(x)
x = x(~isnan(x)); if isempty(x), s = NaN; return; end
s = sqrt(-2 * log(max(abs(mean(exp(1i * x(:)))), eps)));
end
