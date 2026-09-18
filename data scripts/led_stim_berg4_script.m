%% led_stim_berg4_script
% PB imaging + LED (CsChrimson) stim, 6 flies x 6 trials in
% Z:\pablo\lpsp_cschrimson_epg_syt8s\led_stim_berg4\<fly>\<trial>\, where trial in
% {high_stim, high_stim_high_k, low_stim(_3), low_stim_high_k, medium_stim,
% med(ium)_stim_high_k} -- 3 LED intensities x normal/high-[K+] saline, one raw
% ScanImage tif + one WaveSurfer h5 per trial.
%
% Modeled on data scripts/d7_pen_script.m (same 3-section shape: convert raw
% movies -> make one mask per fly -> extract F/dF-F traces into all_data),
% but the raw data here is stored very differently, which this header explains
% before you read the sections below.
%
% ---- reading the raw tif -------------------------------------------------
% Each trial's tif is a raw (unregistered) ScanImage fastZ stack: scan_params.json
% says nplanes=10 per volume (5 real z-planes + 5 flyback planes to ignore).
% IMPORTANT: scan_params.json's nframes_total field is STALE for this dataset --
% it says 67889 for every single trial checked, but every actual tif on disk
% has exactly 42910 raw frames (4291 volumes; confirmed against both the
% ScanImage per-frame timestamps in the tif's own ImageDescription tags, and
% independently against the si_frameclock/si_volumeclock line counts in the h5,
% see below -- both agree exactly with 42910/4291 and disagree with 67889).
% So this script never trusts nframes_total; it counts the real number of IFDs
% in the tif itself (lsb_count_tif_frames, a fast header-only pass) and derives
% volume count from that.
%
% ---- stimulus timing: WaveSurfer h5, not ficTrac -------------------------
% There's no ficTracData_DAQ.mat here. Instead each trial has a WaveSurfer
% wvsrfr_0001.h5 (Janelia WaveSurfer's HDF5 layout: /header + one /sweep_000N),
% sampled at 2 kHz, with digital lines (see /header/DIChannelNames) including
% si_frameclock (one pulse per raw ScanImage frame), si_volumeclock (one pulse
% per completed volume), and led_stim_feedback (the actual LED TTL). Verified
% directly against the data: si_frameclock has exactly 42910 rising edges and
% si_volumeclock exactly 4291, i.e. exactly one edge per raw frame/volume in
% every trial checked -- so si_volumeclock edge times ARE the volume timestamps
% on the DAQ clock, no separate alignment step needed. led_stim_feedback has
% exactly 12 onsets per trial, each a clean 2 s pulse on a 10 s period, matching
% the described stim protocol.
%
% ---- how this differs from d7_pen_script ---------------------------------
% d7_pen_script interpolates the WHOLE image stack up onto ficTrac's fine time
% base before masking. Here the "fine" clock is the 2 kHz DAQ clock, and the
% images are only 32x64 px but ~4300 volumes/trial -- upsampling every pixel to
% 2 kHz (~144000 samples) would mean holding a ~[32x64x144000] array per trial
% for no benefit. Instead only the extracted 1-D traces (F, dF/F, and their
% background-subtracted versions) are interpolated up to the fine DAQ clock;
% the still images (in_stim/out_stim) are built by downsampling the stim
% logical onto the (coarser, volume-rate) image time base instead. Everything
% is otherwise named the same way as d7_pen_script (all_data(i).ft.xb/xf/stims,
% all_data(i).im.*) so it should feel familiar.
%
% ---- F / dF-F definitions --------------------------------------------------
% Per trial, inside the fly's mask:
%   im.F_raw     = mean fluorescence inside the mask (PMT dark-offset
%                  subtracted first, same 1st-percentile trick as d7_pen_script)
%   im.F_bgsub   = im.F_raw minus the mean fluorescence OUTSIDE the mask
%                  (this is exactly d7_pen_script's `trace`)
%   im.dFF_raw   = (F_raw   - F0_raw)   / F0_raw,   F0 = prctile(trace, f0_pct)
%   im.dFF_bgsub = (F_bgsub - F0_bgsub) / F0_bgsub
% F0 for each version is taken from its OWN trace's low percentile (this
% repo's usual convention, see f0_pct usage in epg_dlight_script.m etc.).
% Caveat worth knowing before trusting dFF_bgsub: F_bgsub is already a
% difference of two means and can sit close to zero (or go negative), so its
% own baseline percentile is a much less stable denominator than F_raw's --
% dFF_raw is the more robust of the two; dFF_bgsub is offered mainly as a
% diagnostic for how much the surrounding-tissue background is contributing.
%
% ---- mask: one per fly, not one per trial ---------------------------------
% Like d7_pen_script, ONE mask is made per fly (not per trial) and reused
% across that fly's 6 trials, since the FOV shouldn't move within a session.
% It's built from the across-trial-averaged time-projection of that fly's 6
% videos, then thresholded/joined exactly like d7_pen_script (threshold at the
% 60th percentile, keep the 2 largest components and dilate them together if
% comparably sized -- this is what joins the PB's two brightest halves into one
% mask -- or just the single largest component if one dominates). If you'd
% rather have a mask per trial instead (e.g. if you expect the FOV to drift
% across a session), that's a one-line change in section 2 -- ask and I'll
% adjust it.
%
% Run one %% section at a time.

%% 0. parameters
% pb_register.m (used by section 1b) lives in claude/, a sibling of this
% data scripts/ folder -- not necessarily on the path already, so add it here
% rather than relying on whatever folder MATLAB happened to start in.
claudeDir = fullfile(fileparts(mfilename('fullpath')), '..', 'claude');
if isempty(which('pb_register')) && isfolder(claudeDir)
    addpath(claudeDir);
end
% graph_sort.m (used by section 7's glomerulus clustering) lives at the repo root.
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
if isempty(which('graph_sort')) && isfolder(repoRootDir)
    addpath(repoRootDir);
end

base_dir = 'Z:\pablo\lpsp_cschrimson_epg_syt8s\led_stim_berg4\';

f0_pct              = 7;      % baseline percentile for dF/F (this repo's usual convention)
win_edges_sec       = [-3,9]; % seconds around each stim onset to extract for the aligned-trace plots
overwrite_video     = false;  % re-convert tif -> imgData_sum.mat even if it already exists
overwrite_mask      = false;  % rebuild mask.mat even if it already exists
do_motion_correction = true;  % rigid NoRMCorre registration (via claude/pb_register.m) on the summed-z video
overwrite_reg       = false;  % redo motion correction even if imgData_sum_reg*.mat already exists
n_per_hemisphere    = 10;     % section 7's glomerulus clustering, doubled internally -> 20 total; also used by section 2's mask QC overlay

fly_dirs = lsb_list_subdirs(base_dir);
fprintf('found %d fly folders in %s\n', numel(fly_dirs), base_dir);

%% 1. convert each trial's raw tif into a z-summed, single-plane video, saved in place
% Normally one tif + one h5 per trial folder. A folder can hold more than one
% complete (tif, h5) pair if an acquisition was restarted and both takes were
% kept (e.g. fly 4's medium_stim: two independent, complete 130 s recordings,
% wvsrfr_0001.h5/sweep_0001 and wvsrfr_0002.h5/sweep_0002) -- confirmed valid,
% not a corrupted retake, and kept as two separate trials per your call.
for f = 1:numel(fly_dirs)
    flyDir = fullfile(fly_dirs(f).folder, fly_dirs(f).name);
    trial_dirs = lsb_list_subdirs(flyDir);

    for t = 1:numel(trial_dirs)
        trialDir = fullfile(trial_dirs(t).folder, trial_dirs(t).name);
        acqs = lsb_list_acquisitions(trialDir);

        for a = 1:numel(acqs)
            outFile = fullfile(trialDir, ['imgData_sum' acqs(a).label '.mat']);
            if isfile(outFile) && ~overwrite_video
                fprintf('skip (exists): %s\n', outFile);
                continue
            end

            fprintf('[%s / %s%s] reading raw tif...\n', fly_dirs(f).name, trial_dirs(t).name, acqs(a).label);
            [imgData_sum, sp] = lsb_read_raw_tif_summed(trialDir, acqs(a).tifPath); %#ok<ASGLU>
            save(outFile, 'imgData_sum', 'sp', '-v7.3');
            fprintf('  saved %s  (%d x %d x %d volumes)\n', outFile, size(imgData_sum,1), size(imgData_sum,2), size(imgData_sum,3));
        end
    end
end

%% 1b. motion-correct each trial's summed-z video (rigid NoRMCorre via pb_register)
% The single-plane video can still drift within a 130 s trial. pb_register
% (claude/pb_register.m) runs the same rigid, translation-only NoRMCorre
% registration already used elsewhere in this repo. Saved as a separate
% imgData_sum_reg*.mat next to the raw imgData_sum*.mat (raw is kept for
% comparison, e.g. to judge whether registration actually helped).
if do_motion_correction
    for f = 1:numel(fly_dirs)
        flyDir = fullfile(fly_dirs(f).folder, fly_dirs(f).name);
        trial_dirs = lsb_list_subdirs(flyDir);

        for t = 1:numel(trial_dirs)
            trialDir = fullfile(trial_dirs(t).folder, trial_dirs(t).name);
            acqs = lsb_list_acquisitions(trialDir);

            for a = 1:numel(acqs)
                outFile = fullfile(trialDir, ['imgData_sum_reg' acqs(a).label '.mat']);
                if isfile(outFile) && ~overwrite_reg
                    fprintf('skip (exists): %s\n', outFile);
                    continue
                end

                S = load(fullfile(trialDir, ['imgData_sum' acqs(a).label '.mat']), 'imgData_sum');
                fprintf('[%s / %s%s] motion-correcting (%d volumes)...\n', ...
                    fly_dirs(f).name, trial_dirs(t).name, acqs(a).label, size(S.imgData_sum,3));

                [imgData_sum_reg, regShifts, regDiag] = pb_register(S.imgData_sum, ...
                    'smoothWindow', 5, 'maxShift', [10 15], 'verbose', false, 'figureVisible', false); %#ok<ASGLU>
                close(regDiag.figHandles);

                save(outFile, 'imgData_sum_reg', 'regShifts', '-v7.3');
                fprintf('  saved %s\n', outFile);
            end
        end
    end
end

%% 2. build one PB mask per fly (mean projection averaged across that fly's trials)
for f = 1:numel(fly_dirs)
    flyDir   = fullfile(fly_dirs(f).folder, fly_dirs(f).name);
    maskFile = fullfile(flyDir, 'mask.mat');

    if isfile(maskFile) && ~overwrite_mask
        fprintf('skip (exists): %s\n', maskFile);
        continue
    end

    trial_dirs = lsb_list_subdirs(flyDir);
    projMean = [];
    nProj = 0;
    for t = 1:numel(trial_dirs)
        trialDir = fullfile(trial_dirs(t).folder, trial_dirs(t).name);
        acqs = lsb_list_acquisitions(trialDir);
        for a = 1:numel(acqs)
            img  = lsb_load_video(trialDir, acqs(a).label, do_motion_correction);
            proj = mean(img, 3);
            if isempty(projMean); projMean = proj; else; projMean = projMean + proj; end
            nProj = nProj + 1;
        end
    end
    projMean = projMean / nProj;

    % Originally this was d7_pen_script's threshold-and-join logic verbatim
    % (threshold at the 60th percentile, keep/dilate the 2 biggest components),
    % but on this data it produced jagged, spiky masks -- fine for a whole-mask
    % F/dF-F trace, but bad for section 7's skeletonize-based glomerulus
    % clustering (spiky edges -> spurious skeleton branches -> scrambled
    % glomerulus ordering). Added: a Gaussian smooth before thresholding, an
    % opening to strip thin noise spurs, hole-filling, and a proper imclose
    % (dilate+erode, not dilate-only) to bridge the two hemispheres without
    % permanently ballooning the boundary.
    sigma = 2.5; % px, tuned for this 32x64 FOV by workshopping it against all 6 flies' masks. The dim
    % strip between the two PB hemispheres is right at the threshold cutoff, so with weaker smoothing,
    % pixel noise there produced an unstable, spiky boundary that scrambled the glomerulus ordering
    % (a bright spurious "peak" wandering to a different column per fly, not a fixed scanner artifact --
    % checked directly across all 6 flies before landing on this fix). No single sigma is perfect for
    % all 6 flies: at 2.5, 5/6 are clean and fly 4 has one small detour around glomerulus #12; pushing to
    % 3.0 fixes fly 4 but introduces a worse zigzag in fly 3's ordering, so 2.5 is the better trade-off.
    projSmooth = imgaussfilt(projMean, sigma);
    top_pct = prctile(projSmooth, 98, 'all');
    bot_pct = prctile(projSmooth, 5, 'all');
    projClipped = projSmooth;
    projClipped(projClipped > top_pct) = top_pct;
    projClipped(projClipped < bot_pct) = bot_pct;

    mask = projClipped > prctile(projClipped(:), 60);
    mask = imopen(mask, strel('disk', 1)); % strip thin spurs before we start counting/keeping components
    mask = imfill(mask, 'holes');

    areas = sort([regionprops(mask).Area], 'descend');
    if numel(areas) < 2 || (areas(1) / areas(2)) > 1.5
        mask = bwareafilt(mask, 1);
    else
        mask = bwareafilt(mask, 2);
        mask = imclose(mask, strel('line', 15, 0)); % join the PB's two hemispheres horizontally --
        % a disk-shaped close bridges wherever is geometrically closest, which on this data built a
        % spurious bridge well above the true (much lower, more horizontal) gap between hemispheres
    end
    mask = imfill(mask, 'holes');

    % QC figure: mask boundary + the section 7 glomerulus numbering it implies,
    % so a bad skeleton/ordering is visible here instead of only showing up as
    % a scrambled heatmap several sections later.
    figure(1); clf
    imagesc(projMean); axis equal tight; colormap(bone); hold on
    contour(mask, [0.5 0.5], 'r', 'LineWidth', 1.5);
    try
        qcClusters = lsb_pb_glomeruli(mask, n_per_hemisphere);
        qcCentroids = zeros(2*n_per_hemisphere, 2);
        for c = 1:2*n_per_hemisphere
            [cy, cx] = find(qcClusters == c);
            qcCentroids(c,:) = [mean(cx), mean(cy)];
        end
        % numbered path, not a second colormap sharing this axes with imagesc
        % (a scatter 'CData' index would just get mapped through the same
        % 'bone' colormap as the image, which isn't what we want here)
        plot(qcCentroids(:,1), qcCentroids(:,2), 'y-', 'LineWidth', 1)
        plot(qcCentroids(:,1), qcCentroids(:,2), 'yo', 'MarkerFaceColor', 'y', 'MarkerSize', 4)
        for c = 1:2*n_per_hemisphere
            text(qcCentroids(c,1), qcCentroids(c,2), num2str(c), 'Color', 'g', ...
                'FontSize', 6, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom')
        end
    catch ME
        fprintf('  note: glomerulus QC overlay failed (%s) -- mask itself still saved.\n', ME.message);
    end
    title(fly_dirs(f).name, 'Interpreter', 'none');
    exportgraphics(figure(1), fullfile(flyDir, 'mask_qc.png'), 'Resolution', 150);

    save(maskFile, 'mask');
    fprintf('saved %s\n', maskFile);
end

%% 3. extract F / dF-F traces into all_data, synced to the LED stim via each trial's h5
% Count the total number of (tif,h5) acquisitions up front (usually 6/fly, but
% a fly can have more -- see section 1) just so the ETR estimate below is accurate.
nTotalAcq = 0;
for f = 1:numel(fly_dirs)
    trial_dirs = lsb_list_subdirs(fullfile(fly_dirs(f).folder, fly_dirs(f).name));
    for t = 1:numel(trial_dirs)
        nTotalAcq = nTotalAcq + numel(lsb_list_acquisitions(fullfile(trial_dirs(t).folder, trial_dirs(t).name)));
    end
end

all_data = struct();
i = 0;
tic
for f = 1:numel(fly_dirs)
    flyDir = fullfile(fly_dirs(f).folder, fly_dirs(f).name);
    load(fullfile(flyDir, 'mask.mat'), 'mask');
    trial_dirs = lsb_list_subdirs(flyDir);

    for t = 1:numel(trial_dirs)
        trialDir = fullfile(trial_dirs(t).folder, trial_dirs(t).name);
        acqs = lsb_list_acquisitions(trialDir);
        [intensity, isHighK] = lsb_parse_condition(trial_dirs(t).name);

        for a = 1:numel(acqs)
            i = i + 1;
            fprintf('[%d] processing %s / %s%s ... ', i, fly_dirs(f).name, trial_dirs(t).name, acqs(a).label);

            imgData_sum = lsb_load_video(trialDir, acqs(a).label, do_motion_correction);
            sync = lsb_h5_stim_timing(acqs(a).h5Path);

            nVol = size(imgData_sum, 3);
            if numel(sync.t_volume) ~= nVol
                warning('led_stim_berg4_script:volMismatch', ...
                    '%s%s: %d volume-clock edges in the h5 but %d volumes in the video; truncating to the shorter of the two.', ...
                    trialDir, acqs(a).label, numel(sync.t_volume), nVol);
                n = min(numel(sync.t_volume), nVol);
                sync.t_volume = sync.t_volume(1:n);
                imgData_sum   = imgData_sum(:,:,1:n);
                nVol = n;
            end

            % remove PMT dark offset (same trick as d7_pen_script)
            imgData_sum = imgData_sum - prctile(imgData_sum, 1, 'all');

            imgData_2D  = reshape(imgData_sum, [], nVol);
            F_raw_vol   = mean(imgData_2D(mask(:), :), 1);
            F_bg_vol    = mean(imgData_2D(~mask(:), :), 1);
            F_bgsub_vol = F_raw_vol - F_bg_vol;

            F0_raw   = prctile(F_raw_vol, f0_pct);
            F0_bgsub = prctile(F_bgsub_vol, f0_pct);
            dFF_raw_vol   = (F_raw_vol   - F0_raw)   / F0_raw;
            dFF_bgsub_vol = (F_bgsub_vol - F0_bgsub) / F0_bgsub;

            % still images: downsample the stim logical onto the (coarser) volume clock
            stims_vol = interp1(sync.t_fine, double(sync.stims), sync.t_volume, 'nearest', 'extrap') > 0.5;

            all_data(i).meta.trialDir  = trialDir;
            all_data(i).meta.fly       = fly_dirs(f).name;
            all_data(i).meta.trial     = [trial_dirs(t).name acqs(a).label]; % label distinguishes repeated acquisitions
            all_data(i).meta.acqLabel  = acqs(a).label; % '' or '_0001'/'_0002' -- needed to reload this acquisition's video (section 7)
            all_data(i).meta.intensity = intensity; % 1=low, 2=medium, 3=high
            all_data(i).meta.highK     = isHighK;

            all_data(i).ft.xb        = sync.t_volume; % volume-clock timestamps (DAQ seconds)
            all_data(i).ft.xf        = sync.t_fine;   % native 2 kHz DAQ clock (DAQ seconds)
            all_data(i).ft.stims     = sync.stims;    % LED stim logical, native 2 kHz resolution
            all_data(i).ft.stims_vol = stims_vol;     % LED stim logical, downsampled onto ft.xb

            % native per-volume resolution (~33 Hz) -- use these for plotting/
            % analysis; they're the real sampling rate of the calcium signal
            all_data(i).im.F_raw_vol     = F_raw_vol;
            all_data(i).im.F_bgsub_vol   = F_bgsub_vol;
            all_data(i).im.dFF_raw_vol   = dFF_raw_vol;
            all_data(i).im.dFF_bgsub_vol = dFF_bgsub_vol;

            % upsampled onto the native 2 kHz DAQ clock (ft.xf) -- only useful if you
            % need the trace on the exact same clock as some other fine (DAQ-rate)
            % signal; for plotting the calcium signal alone, prefer the _vol fields
            % above (upsampling is pure linear interpolation, it adds no information,
            % and turns dense per-pulse plots into a noisy-looking smear)
            all_data(i).im.F_raw     = interp1(sync.t_volume, F_raw_vol,     sync.t_fine, 'linear', 'extrap');
            all_data(i).im.F_bgsub   = interp1(sync.t_volume, F_bgsub_vol,   sync.t_fine, 'linear', 'extrap');
            all_data(i).im.dFF_raw   = interp1(sync.t_volume, dFF_raw_vol,   sync.t_fine, 'linear', 'extrap');
            all_data(i).im.dFF_bgsub = interp1(sync.t_volume, dFF_bgsub_vol, sync.t_fine, 'linear', 'extrap');

            all_data(i).im.in_stim  = mean(imgData_sum(:,:, stims_vol), 3);
            all_data(i).im.out_stim = mean(imgData_sum(:,:,~stims_vol), 3);

            fprintf('ETR: %.2f min\n', toc/i * (nTotalAcq - i) / 60);
        end
    end
end

save(fullfile(base_dir, 'all_data.mat'), 'all_data', '-v7.3');
fprintf('saved %s\n', fullfile(base_dir, 'all_data.mat'));

%% 4. condition indices
intensity_ind  = arrayfun(@(a) a.meta.intensity, all_data); % 1=low, 2=medium, 3=high
highK_ind      = arrayfun(@(a) a.meta.highK, all_data);
intensity_name = {'low','medium','high'};

%% 5. still images: in-stim vs. out-of-stim vs. difference, one row per trial
outStims = arrayfun(@(a) a.im.out_stim, all_data, 'UniformOutput', false);
max_pix  = prctile(cat(1, outStims{:}), 99, 'all'); % shared color scale, set from baseline images
rwb      = lsb_redwhiteblue(256); % diverging colormap for the difference column only

% Diff color scale set from the difference images themselves (not max_pix,
% which is a raw-intensity scale and made the diff column washed out/pale --
% the actual in-out signal is much smaller than the raw fluorescence range).
diffImgs  = arrayfun(@(a) a.im.in_stim - a.im.out_stim, all_data, 'UniformOutput', false);
diff_pix  = prctile(abs(cat(1, diffImgs{:})), 99, 'all');

figure(2); clf
set(gcf, 'Position', [50 50 900 max(400, 90*numel(all_data))], 'Color', 'w')
n = numel(all_data);
for i = 1:n
    subplot(n,3,3*(i-1)+1)
    imagesc(all_data(i).im.out_stim); clim([0, max_pix]); axis equal tight
    ylabel({all_data(i).meta.fly, all_data(i).meta.trial}, 'Rotation', 0, 'Interpreter', 'none')

    subplot(n,3,3*(i-1)+2)
    imagesc(all_data(i).im.in_stim); clim([0, max_pix]); axis equal tight

    subplot(n,3,3*(i-1)+3)
    imagesc(diffImgs{i}); clim([-diff_pix diff_pix]); axis equal tight
    colormap(gca, rwb); % red = brighter during stim, blue = dimmer during stim, white = no change
end
tmp = get(gcf,'Children');
set(tmp, 'XTick', [], 'YTick', [])
subplot(n,3,1); title('Out Stim')
subplot(n,3,2); title('In Stim')
subplot(n,3,3); title('In - Out')

exportgraphics(figure(2), fullfile(base_dir, 'stillImages_allTrials.png'), 'Resolution', 150);
fprintf('saved %s\n', fullfile(base_dir, 'stillImages_allTrials.png'));

%% 5b. difference images only, tiled tightly so each one is bigger
% Same diff images/color scale as section 5's third column, just laid out as a
% grid (rows = fly, cols = that fly's trials in processing order) instead of a
% single narrow column, so each tile is much larger and easier to read.
flyOfTrial = arrayfun(@(a) string(a.meta.fly), all_data);
uFlies  = unique(flyOfTrial, 'stable');
nRows   = numel(uFlies);
nCols   = max(arrayfun(@(f) sum(flyOfTrial == f), uFlies));

figure(4); clf
set(gcf, 'Position', [50 50 220*nCols 220*nRows], 'Color', 'w')
tiledlayout(nRows, nCols, 'TileSpacing', 'compact', 'Padding', 'compact');

curFly = ""; r = 0; c = 0;
for i = 1:numel(all_data)
    if flyOfTrial(i) ~= curFly
        curFly = flyOfTrial(i);
        r = r + 1;
        c = 0;
    end
    c = c + 1;

    nexttile((r-1)*nCols + c)
    imagesc(diffImgs{i}); clim([-diff_pix diff_pix]); axis image off
    colormap(gca, rwb)
    title(all_data(i).meta.trial, 'Interpreter', 'none', 'FontSize', 8)
    if c == 1
        ylabel(erase(curFly, '_epg_syt8s_lpsp_cschrimson'), 'Rotation', 0, ...
            'HorizontalAlignment', 'right', 'Interpreter', 'none', 'FontWeight', 'bold')
    end
end
sgtitle('In - Out (red = brighter in stim, blue = dimmer)')

exportgraphics(figure(4), fullfile(base_dir, 'diffImages_tiled.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(base_dir, 'diffImages_tiled.png'));

%% 6. overlapped traces per condition, aligned to stim onset
% Plotted at the native per-volume (~33 Hz) sampling rate, not the upsampled
% 2 kHz DAQ clock -- upsampling is pure linear interpolation (adds no real
% information) and made individual pulses look like a noisy smear rather than
% a legible line.
trace_field = 'dFF_bgsub_vol'; % switch to 'dFF_raw_vol' / 'F_raw_vol' / 'F_bgsub_vol' to compare

figure(3); clf
set(gcf, 'Position', [50 50 1400 800], 'Color', 'w')
for row = 0:1                    % 0 = normal [K+], 1 = high [K+]
    for col = 1:3                 % low / medium / high intensity
        ax = subplot(2,3,col + row*3); hold(ax,'on')
        trial_idx = find(intensity_ind == col & highK_ind == row);

        pulses = {};
        for i = trial_idx
            trace   = all_data(i).im.(trace_field)(:)'; % force row vector
            fs      = 1 / mean(diff(all_data(i).ft.xb));
            winSamp = round(win_edges_sec(1)*fs) : round(win_edges_sec(2)*fs); % 1 x M
            onsets  = find(diff(all_data(i).ft.stims_vol) > 0) + 1;
            onsets  = onsets(:); % K x 1

            stim_wins = onsets + winSamp; % K x M index matrix (broadcast), one row per onset
            valid = stim_wins(:,1) >= 1 & stim_wins(:,end) <= numel(trace); % drop pulses too close to trial edges
            stim_wins = stim_wins(valid, :);

            if ~isempty(stim_wins)
                pulses{end+1} = trace(stim_wins); %#ok<SAGROW> % (valid onsets) x M
            end
        end

        if ~isempty(pulses)
            pulses = cell2mat(pulses(:));
            t_win  = winSamp / fs;
            plot(ax, t_win, pulses', 'Color', [0,0,0,.08])
            h = lsb_plot_sem(ax, t_win, pulses); h.FaceColor = 'r'; h.EdgeColor = 'none';
            tmp_y = ylim(ax);
            patch(ax, [0 2 2 0], tmp_y([1 1 2 2]), 'r', 'FaceAlpha', .1, 'EdgeColor', 'none')
            uistack(findobj(ax,'Type','patch'), 'bottom')
        end

        title(ax, sprintf('%s%s', intensity_name{col}, repmat(' + high K^+', 1, row)))
        if col == 1; ylabel(ax, trace_field, 'Interpreter', 'none'); end
        if row == 1; xlabel(ax, 'time from stim onset (s)'); end
    end
end
tmp = get(gcf,'Children');
linkaxes(tmp, 'xy')
set(tmp, 'Color', 'none')
sgtitle(sprintf('peri-stim %s, mean \\pm SEM (rows: normal / high K^+; cols: intensity)', trace_field), 'Interpreter', 'tex')

exportgraphics(figure(3), fullfile(base_dir, 'peristim_traces_by_condition.png'), 'Resolution', 200);
fprintf('saved %s\n', fullfile(base_dir, 'peristim_traces_by_condition.png'));

%% 7. divide the PB into glomeruli and plot per-fly peri-stim heatmaps
% Same clustering approach as this lab's other scripts (see process_im in e.g.
% epg_dlight_script.m / dopamine_ionto_script.m, and graph_sort.m): skeletonize
% the mask, order the skeleton into one continuous path, resample it into
% evenly spaced centroids, then assign every mask pixel to its nearest
% centroid. n_per_hemisphere (section 0) is doubled internally (matches
% process_im's own convention) -- 10 gives 20 total glomeruli, as requested.
glom_trace_field = 'dFF'; % per-glomerulus dF/F, same f0_pct baseline convention as F_raw_vol/dFF_raw_vol

for fi = 1:numel(uFlies)
    flyName = uFlies(fi);
    load(fullfile(base_dir, flyName, 'mask.mat'), 'mask');
    clusterIdx = lsb_pb_glomeruli(mask, n_per_hemisphere);
    nClusters  = 2 * n_per_hemisphere;

    trialsThisFly = find(flyOfTrial == flyName);

    figure; clf
    set(gcf, 'Position', [50 50 260*numel(trialsThisFly) 450], 'Color', 'w')
    tiledlayout(1, numel(trialsThisFly), 'TileSpacing', 'compact', 'Padding', 'loose');

    for k = 1:numel(trialsThisFly)
        idx = trialsThisFly(k);
        trialDir = all_data(idx).meta.trialDir;

        img = lsb_load_video(trialDir, all_data(idx).meta.acqLabel, do_motion_correction);
        img = img(:, :, 1:size(all_data(idx).ft.xb,1)); % match any truncation applied in section 3
        img = img - prctile(img, 1, 'all'); % PMT dark offset, same as section 3

        img_2d      = reshape(img, [], size(img,3));
        centroidLog = false(nClusters, size(img_2d,1));
        for c = 1:nClusters
            centroidLog(c, clusterIdx(:) == c) = true;
        end
        f_cluster   = double(centroidLog) * double(img_2d) ./ sum(centroidLog, 2); % nClusters x nVolumes
        f0_cluster  = prctile(f_cluster, f0_pct, 2);
        dff_cluster = (f_cluster - f0_cluster) ./ f0_cluster;

        fs_vol  = 1 / mean(diff(all_data(idx).ft.xb));
        winSamp = round(win_edges_sec(1)*fs_vol) : round(win_edges_sec(2)*fs_vol);
        onsets  = find(diff(all_data(idx).ft.stims_vol) > 0) + 1;
        valid   = onsets + winSamp(1) >= 1 & onsets + winSamp(end) <= size(dff_cluster,2);
        onsets  = onsets(valid);

        slices = nan(numel(onsets), nClusters, numel(winSamp));
        for o = 1:numel(onsets)
            slices(o,:,:) = dff_cluster(:, onsets(o) + winSamp);
        end
        avgHeatmap = squeeze(mean(slices, 1, 'omitnan')); % nClusters x nWinSamples, averaged over repetitions

        nexttile(k)
        t_win = winSamp / fs_vol;
        imagesc(t_win, 1:nClusters, avgHeatmap)
        hold on
        plot([0 0], [0.5 nClusters+0.5], 'k--', [2 2], [0.5 nClusters+0.5], 'k--')
        colormap(gca, rwb)
        clim([-1 1] * max(abs(avgHeatmap(:))))
        axis tight
        title(all_data(idx).meta.trial, 'Interpreter', 'none', 'FontSize', 9)
        xlabel('time from stim onset (s)')
        if k == 1; ylabel('glomerulus'); end
        colorbar
    end
    sgtitle(sprintf('%s -- per-glomerulus %s, mean over stim repetitions', flyName, glom_trace_field), 'Interpreter', 'none')

    outFile = fullfile(base_dir, sprintf('glomHeatmap_%s.png', flyName));
    exportgraphics(gcf, outFile, 'Resolution', 200);
    fprintf('saved %s\n', outFile);
end

%% local functions

function clusterIdx = lsb_pb_glomeruli(mask, nPerHemisphere)
% Divide a single-connected-arch PB mask into 2*nPerHemisphere equally sized
% glomeruli. This is the same skeletonize -> order -> resample-into-centroids
% -> nearest-centroid-assign approach as process_im in this lab's other
% scripts (e.g. epg_dlight_script.m, dopamine_ionto_script.m), trimmed down to
% just the clustering step (no bump-fitting) since that's all this needs.
% Requires graph_sort.m (repo root) to order the skeleton into a path -- see
% its own comments for why this assumes mask is a single open arch, not a loop.
%
% clusterIdx: same [Y,X] size as mask; 0 outside the mask, 1..2*nPerHemisphere
% inside, numbered along the skeleton from one end to the other.

[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('lsb_pb_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
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

function d = lsb_list_subdirs(parentDir)
% Non-hidden subfolders of parentDir, alphabetical (dir()'s default order).
d = dir(parentDir);
d = d([d.isdir] & ~startsWith({d.name}, '.'));
end

function n = lsb_count_tif_frames(tifPath)
% Fast frame count: walks IFD headers only, never decodes pixel data.
t = Tiff(tifPath, 'r');
cleanupTiff = onCleanup(@() close(t)); %#ok<NASGU>
n = 1;
while ~t.lastDirectory()
    t.nextDirectory();
    n = n + 1;
end
end

function img = lsb_load_video(trialDir, label, useReg)
% Load either the raw or motion-corrected summed-z video for one acquisition,
% abstracting away the two different file/variable names (see section 1b).
if useReg
    S = load(fullfile(trialDir, ['imgData_sum_reg' label '.mat']), 'imgData_sum_reg');
    img = S.imgData_sum_reg;
else
    S = load(fullfile(trialDir, ['imgData_sum' label '.mat']), 'imgData_sum');
    img = S.imgData_sum;
end
end

function acqs = lsb_list_acquisitions(trialDir)
% One entry per (tif, matching h5) pair found in trialDir. Normally exactly
% one; a folder can hold more if an acquisition was restarted and both takes
% were kept (see this script's header / section 1) -- each is then treated as
% its own trial. Tif files are named '..._<seq>.tif' and their matching h5 is
% 'wvsrfr_<seq>.h5' (same sequence number, zero-padded to 4 digits).
tifList = dir(fullfile(trialDir, '*.tif'));
if isempty(tifList)
    error('led_stim_berg4_script:tif', 'No .tif found in %s.', trialDir);
end

acqs = struct('tifPath', {}, 'h5Path', {}, 'seq', {}, 'label', {});
for k = 1:numel(tifList)
    tok = regexp(tifList(k).name, '_(\d+)\.tif$', 'tokens', 'once');
    if isempty(tok)
        error('led_stim_berg4_script:tifName', 'Could not parse a trailing sequence number from tif name "%s".', tifList(k).name);
    end
    seq = str2double(tok{1});
    h5Path = fullfile(trialDir, sprintf('wvsrfr_%04d.h5', seq));
    if ~isfile(h5Path)
        error('led_stim_berg4_script:h5match', 'No matching %s for %s in %s.', sprintf('wvsrfr_%04d.h5', seq), tifList(k).name, trialDir);
    end
    acqs(end+1).tifPath = fullfile(tifList(k).folder, tifList(k).name); %#ok<AGROW>
    acqs(end).h5Path    = h5Path;
    acqs(end).seq       = seq;
end
[~, order] = sort([acqs.seq]);
acqs = acqs(order);

if numel(acqs) == 1
    acqs.label = ''; % keep the plain imgData_sum.mat name for the (overwhelmingly common) single-acquisition case
else
    fprintf('  note: %s has %d separate (tif,h5) acquisitions; keeping all as separate trials.\n', trialDir, numel(acqs));
    for k = 1:numel(acqs)
        acqs(k).label = sprintf('_%04d', acqs(k).seq);
    end
end
end

function [imgSum, sp] = lsb_read_raw_tif_summed(trialDir, tifPath)
% Read every raw ScanImage frame in tifPath and collapse it to a single-plane
% video [Y,X,nVolumes], summed over the valid (non-flyback) z-planes only.
% scan_params.json (shared by every acquisition in trialDir) gives the plane
% layout; see this script's header for why its nframes_total is NOT used to
% size the read.
sp = jsondecode(fileread(fullfile(trialDir, 'scan_params.json')));
if sp.nchannels ~= 1
    error('led_stim_berg4_script:multiChannel', '%s has %d saved channels; this script assumes 1.', trialDir, sp.nchannels);
end
validPlanes = setdiff(1:sp.nplanes, sp.flyback_planes(:)' + 1); % json flyback is 0-indexed
if numel(validPlanes) ~= sp.nplanes_valid
    error('led_stim_berg4_script:flybackMismatch', ...
        '%s: scan_params.json says nplanes_valid=%d but nplanes/flyback_planes leaves %d.', ...
        trialDir, sp.nplanes_valid, numel(validPlanes));
end

nFramesOnDisk = lsb_count_tif_frames(tifPath);
framesPerVol  = sp.nplanes * sp.nchannels;
nVolumes      = floor(nFramesOnDisk / framesPerVol);
leftover      = nFramesOnDisk - nVolumes * framesPerVol;
if leftover > 0
    fprintf('  note: %d leftover raw frame(s) after %d complete volumes, dropping the trailing partial volume\n', leftover, nVolumes);
end
if nFramesOnDisk ~= sp.nframes_total
    fprintf('  note: scan_params.json nframes_total=%d does not match the %d raw frames actually on disk; using the on-disk count.\n', ...
        sp.nframes_total, nFramesOnDisk);
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
    if mod(k, 10000) == 0 || k == nRawToRead
        fprintf('  read %d/%d raw frames (%.0f frames/s)\n', k, nRawToRead, k/toc);
    end
end

raw = reshape(raw, sp.px_height, sp.px_width, sp.nplanes, nVolumes);
imgSum = squeeze(sum(raw(:,:,validPlanes,:), 3)); % Y x X x volume
sp.nVolumes = nVolumes;
end

function s = lsb_h5_stim_timing(h5Path)
% Pull the LED stim TTL and per-volume timestamps out of a trial's WaveSurfer
% h5, all on the h5's own (native 2 kHz) DAQ clock. See this script's header
% for why si_volumeclock edges can be used directly as volume timestamps.
% NB: the sweep group name isn't always '/sweep_0001' -- it's found dynamically
% below, since a restarted acquisition's second h5 can be '/sweep_0002' etc.
info = h5info(h5Path);
sweepNames = {info.Groups.Name};
sweepNames = sweepNames(~strcmp(sweepNames, '/header'));
if numel(sweepNames) ~= 1
    error('led_stim_berg4_script:h5sweep', 'Expected exactly one sweep group in %s, found %d.', h5Path, numel(sweepNames));
end

fs      = double(h5read(h5Path, '/header/AcquisitionSampleRate'));
diNames = strtrim(string(h5read(h5Path, '/header/DIChannelNames')));
digi    = int32(h5read(h5Path, [sweepNames{1} '/digitalScans']));

volCh = find(diNames == "si_volumeclock", 1);
ledCh = find(diNames == "led_stim_feedback", 1);
if isempty(volCh) || isempty(ledCh)
    error('led_stim_berg4_script:h5chan', ...
        'Could not find si_volumeclock/led_stim_feedback in %s DIChannelNames (found: %s).', h5Path, strjoin(diNames, ', '));
end

volBit = bitget(digi, volCh);
ledBit = bitget(digi, ledCh);

N = numel(digi);
s.fs       = fs;
s.t_fine   = (0:N-1)' / fs;
s.stims    = logical(ledBit);
s.t_volume = s.t_fine(find(diff(volBit) > 0) + 1);
end

function [intensity, isHighK] = lsb_parse_condition(trialName)
% Trial folder names aren't perfectly consistent across flies (e.g.
% 'low_stim_3' vs 'low_stim', 'med_stim_high_k' vs 'medium_stim_high_k'), so
% match on leading substrings/contains rather than exact names.
n = lower(trialName);
if startsWith(n, 'low')
    intensity = 1;
elseif startsWith(n, 'med')
    intensity = 2;
elseif startsWith(n, 'high')
    intensity = 3;
else
    error('led_stim_berg4_script:cond', 'Could not parse stim intensity from trial folder name "%s".', trialName);
end
isHighK = contains(n, 'high_k');
end

function cmap = lsb_redwhiteblue(n)
% Diverging colormap: blue (low) -> white (zero) -> red (high). n should be even.
half = n / 2;
blueToWhite = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
whiteToRed  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];
cmap = [blueToWhite; whiteToRed];
end

function h = lsb_plot_sem(ax, t, x)
t  = reshape(t, 1, []);
m1 = mean(x, 1, 'omitnan');
s1 = std(x, 1, 'omitnan') ./ sqrt(sum(~isnan(x), 1));

idx = ~isnan(m1);
m1 = m1(idx); s1 = s1(idx); t = t(idx);

h = patch(ax, [t, fliplr(t)], [m1+s1, fliplr(m1-s1)], 'r', 'FaceAlpha', .5);
end
