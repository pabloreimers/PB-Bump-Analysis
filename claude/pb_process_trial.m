%% pb_process_trial
% Interactive, section-by-section pipeline: raw ScanImage tif -> scan-noise
% removed -> (optionally) shot-noise reduced -> motion-corrected (PB-stable)
% movie, for one trial folder.
%
% Scan noise (the periodic cross-hatching) and shot noise (sparse bright
% single-pixel puncta from low photon counts) are different noise sources.
% Scan noise is removed with a short TEMPORAL mean (pb_remove_scan_noise) --
% NOT a spatial filter. An earlier spatial-frequency notch approach
% (pb_remove_line_noise, now deprecated) was found to smear real glomerulus
% structure, because the noise and the real signal are not spectrally
% separable along the fast-scan axis on this data. See pb_remove_scan_noise.m
% and pb_remove_line_noise.m for the full account.
%
% Run this one %% section at a time (Ctrl+Enter). Every section prints what it
% did and/or pops a figure, so if something looks wrong you can stop and poke at
% the workspace right there instead of re-running a 10-minute batch job to find
% out. Nothing here silently swallows errors -- if a step fails, fix it before
% moving to the next section.
%
% Expects a trial folder like:
%   Z:\pablo\lpsp_p2x2_walking\<date>\fly <n>\<trial>\
% containing scan_params.json and one raw ScanImage .tif (see pb_scan_metadata).

%% 0. parameters you might actually want to change
channel_to_process = 1;      % 1 = GCaMP (matches show_single_PB_trial.m convention)
test_n_volumes     = 300;    % how many volumes to use for the quick look/tuning pass
scan_noise_window  = 10;     % frames, temporal mean window for pb_remove_scan_noise;
                              % re-check with pb_diagnose_line_noise -- this is NOT the
                              % same as the 30-frame smoothing already used downstream in
                              % show_single_PB_trial.m; the point is to get away with less
do_shot_noise_removal = false; % 3x3 median filter -- has a real cost (blurs small sharp
                                % features along with noise); inspect section 3b before enabling
smooth_window_for_reg = 5;   % frames, temporal smoothing before shift estimation (0 = off)
max_shift_px = [25 25];      % [row col] max rigid shift NoRMCorre is allowed to find

%% 1. pick a trial folder and read its metadata (fails loudly if json/tif are inconsistent)
expDir = uigetdir('Z:\pablo\lpsp_p2x2_walking\', 'Select a trial folder (contains scan_params.json + raw .tif)');
if isequal(expDir, 0); error('pb_process_trial:cancelled', 'No folder selected.'); end

md = pb_scan_metadata(expDir);
fprintf('%s\n  %d valid planes (of %d), %d channel(s), %d complete volumes @ %.2f Hz\n', ...
    md.tifPath, md.nplanes_valid, md.nplanes, md.nchannels, md.n_volumes, md.fps);

%% 2. quick look: load a small chunk and diagnose the line noise before committing to a full read
[testMovie, testMd] = pb_load_raw_tif(expDir, 'channel', channel_to_process, ...
    'volumes', [1, min(test_n_volumes, md.n_volumes)]);
testMovieSum = squeeze(sum(testMovie, 3)); % Y x X x T, summed over valid z-planes

figure(Name = 'raw test chunk, single volume vs time-average'); clf
subplot(1,2,1); imagesc(testMovieSum(:,:,1)); axis image; colorbar; title('single volume, raw')
subplot(1,2,2); imagesc(mean(testMovieSum,3)); axis image; colorbar; title(sprintf('mean of %d volumes, raw', size(testMovieSum,3)))

statsBefore = pb_diagnose_line_noise(testMovieSum, 'title', 'raw test chunk'); %#ok<NASGU>

%% 3. remove scan (line) noise on the test chunk and confirm the peak is gone
% Judge this on a SINGLE frame (not time-averaged) zoomed on the brightest part
% of the bridge -- that's where the earlier spatial-notch approach failed and a
% time-average would hide both the noise and any artifact either way.
testMovieClean = pb_remove_scan_noise(testMovieSum, scan_noise_window);
statsAfter = pb_diagnose_line_noise(testMovieClean, 'title', 'after scan-noise removal'); %#ok<NASGU>

midFrame = round(size(testMovieSum,3)/2);
figure(Name = 'scan noise removal, before/after (single frame, no time-averaging)'); clf
subplot(1,2,1); imagesc(testMovieSum(:,:,midFrame));   axis image; colorbar; title('before (raw)')
subplot(1,2,2); imagesc(testMovieClean(:,:,midFrame)); axis image; colorbar; title('after pb\_remove\_scan\_noise')

% ---> if statsAfter still shows a strong peak at the same frequency as
%      statsBefore, or the bright glomerulus clusters in the "after" panel
%      look smeared/blobby compared to "before", STOP: increase
%      scan_noise_window, or something is wrong -- do not proceed assuming
%      more temporal smoothing will fix a spatial artifact (it won't).

%% 3b. optional: reduce shot noise (sparse bright single-pixel puncta) on the test chunk
% Off by default -- compare the two panels below and decide if it's worth the
% tradeoff (see pb_remove_shot_noise.m) before turning it on for the full trial.
if do_shot_noise_removal
    testMovieDenoised = pb_remove_shot_noise(testMovieClean);
else
    testMovieDenoised = testMovieClean;
end

figure(Name = 'shot noise removal, before/after (test chunk, single frame)'); clf
subplot(1,2,1); imagesc(testMovieClean(:,:,midFrame));    axis image; colorbar; title('before (scan-noise-removed only)')
subplot(1,2,2); imagesc(testMovieDenoised(:,:,midFrame)); axis image; colorbar; title(sprintf('after (shot-noise removal %s)', mat2str(do_shot_noise_removal)))

%% 4. quick look: register the test chunk to see if the PB stabilizes
[testMovieReg, ~, regDiagTest] = pb_register(testMovieDenoised, ... %#ok<ASGLU>
    'smoothWindow', smooth_window_for_reg, 'maxShift', max_shift_px);

% ---> compare the two kymograph panels in the figure this just made: the PB
%      glomeruli columns (vertical stripes) should be flat/stationary in the
%      "after registration" panel. If they still drift, increase max_shift_px
%      or reconsider smooth_window_for_reg before running on the full trial.

%% 5. happy with the test chunk? now do the full trial (this is the slow part)
[fullMovie, fullMd] = pb_load_raw_tif(expDir, 'channel', channel_to_process); %#ok<ASGLU>
fullMovieSum = squeeze(sum(fullMovie, 3)); % Y x X x T
clear fullMovie % the per-plane data is no longer needed and is large

fullMovieClean = pb_remove_scan_noise(fullMovieSum, scan_noise_window);
clear fullMovieSum

if do_shot_noise_removal
    fullMovieDenoised = pb_remove_shot_noise(fullMovieClean);
else
    fullMovieDenoised = fullMovieClean;
end

[fullMovieReg, regShifts, regDiagFull] = pb_register(fullMovieDenoised, ... %#ok<ASGLU>
    'smoothWindow', smooth_window_for_reg, 'maxShift', max_shift_px);

%% 6. save, following this repo's existing imgData_reg.mat / imgData_smooth_reg.mat convention
outDir = fullfile(expDir, 'registration');
if ~isfolder(outDir); mkdir(outDir); end

imgData_denoised = fullMovieDenoised; %#ok<NASGU>
save(fullfile(outDir, 'imgData_denoised.mat'), 'imgData_denoised', 'scan_noise_window', ...
    'do_shot_noise_removal', '-v7.3');

imgData_denoised_reg = fullMovieReg; %#ok<NASGU>
save(fullfile(outDir, 'imgData_denoised_reg.mat'), 'imgData_denoised_reg', 'regShifts', ...
    'scan_noise_window', 'smooth_window_for_reg', 'max_shift_px', '-v7.3');

fprintf('saved %s\n', fullfile(outDir, 'imgData_denoised.mat'));
fprintf('saved %s\n', fullfile(outDir, 'imgData_denoised_reg.mat'));

%% 7. side-by-side playback: raw vs. denoised+registered (nothing written to disk)
% Self-contained -- only needs expDir/md from section 1 and the parameters from
% section 0. Loads its own chunk so it doesn't care whether you've run sections
% 2-6. Judge on THIS (minimally smoothed / effectively raw playback), not on the
% heavily time-averaged figures from sections 2-4 -- those only show that a big
% temporal average hides the noise either way, which is exactly the crutch we're
% trying to get rid of.
video_n_volumes = 1000; % frames to compare -- enough to judge real registration drift
video_playback_pause = 0.02; % seconds between displayed frames; 0 for as-fast-as-possible

[videoRaw, ~] = pb_load_raw_tif(expDir, 'channel', channel_to_process, ...
    'volumes', [1, min(video_n_volumes, md.n_volumes)]);
videoRawSum = squeeze(sum(videoRaw, 3)); % Y x X x T
clear videoRaw

videoClean = pb_remove_scan_noise(videoRawSum, scan_noise_window);
if do_shot_noise_removal
    videoClean = pb_remove_shot_noise(videoClean);
end
videoReg = pb_register(videoClean, 'smoothWindow', smooth_window_for_reg, ...
    'maxShift', max_shift_px, 'figureVisible', false);

clim = [0, prctile(videoRawSum(:), 99.5)]; % shared color scale: raw sets it so brightness is comparable

figure(Name = 'raw vs. denoised+registered playback'); clf
subplot(1,2,1)
hRaw = imagesc(videoRawSum(:,:,1), clim); axis image; colormap gray; colorbar
title('raw (no denoise, no registration)')
subplot(1,2,2)
hReg = imagesc(videoReg(:,:,1), clim); axis image; colormap gray; colorbar
title('denoised + registered')
supT = sgtitle('');

for f = 1:size(videoRawSum,3)
    hRaw.CData = videoRawSum(:,:,f);
    hReg.CData = videoReg(:,:,f);
    supT.String = sprintf('frame %d / %d', f, size(videoRawSum,3));
    drawnow
    if video_playback_pause > 0; pause(video_playback_pause); end
end
