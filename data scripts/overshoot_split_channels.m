%% overshoot_split_channels
% Two-channel (ScanImage channelSave = [1 2]) overshoot trials: the tif holds
% the channels interleaved frame by frame (plane 1 ch1, plane 1 ch2, plane 2
% ch1, ...; 13 planes incl. 3 flyback -> 26 frames/volume), which
% overshoot_bar_stim_script.m cannot read (it assumes one channel). This
% reads the tif once, z-sums the 10 valid planes of each channel, and writes
% one working folder per channel inside the trial folder:
%   <trial>\ch1\ and <trial>\ch2\ , each with imgData_sum.mat (imgData_sum, sp)
%   -- the cache file overshoot_bar_stim_script.m looks for -- plus a copy of
%   the WaveSurfer h5, so the script can then be run with obs_trialDir pointing
%   at each folder (motion correction, mask, glomeruli, bump, figures all land
%   in that folder). Also saves planeMeans_ch1/2.png (per-plane mean images) so
%   it is easy to see which channel is which.
% um/px follows the nominal 0.3164 um/px at zoom 4 (256 px across).
%   matlab -batch "osc_trialDir='Z:\...\trial_002\'; run('overshoot_split_channels.m')"

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
addpath(fullfile(thisDir, '..', 'claude'));
trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20261019-1_epg_syt8s_empty_cschrimson\trial_002\';
if exist('osc_trialDir', 'var') && ~isempty(osc_trialDir), trialDir = osc_trialDir; end
if trialDir(end) ~= filesep, trialDir(end+1) = filesep; end

tifList = dir(fullfile(trialDir, '*.tif')); assert(numel(tifList) == 1, 'expected one tif in %s', trialDir);
h5List = dir(fullfile(trialDir, '*.h5')); assert(numel(h5List) == 1, 'expected one h5 in %s', trialDir);
tifPath = fullfile(tifList(1).folder, tifList(1).name);
fid = fopen(tifPath, 'r'); hdr = fread(fid, 80000, 'uint8=>char')'; fclose(fid);
getNum = @(key) str2double(regexp(hdr, ['SI\.' regexptranslate('escape', key) '\s*=\s*([-0-9.eE]+)'], 'tokens', 'once'));
framesPerVol = getNum('hStackManager.numFramesPerVolumeWithFlyback');
nPlanesValid = getNum('hStackManager.actualNumSlices');
zoom = getNum('hRoiManager.scanZoomFactor'); W = getNum('hRoiManager.pixelsPerLine'); H = getNum('hRoiManager.linesPerFrame');
chTok = regexp(hdr, 'SI\.hChannels\.channelSave = \[([^\]]*)\]', 'tokens', 'once'); chSave = str2num(chTok{1}); %#ok<ST2NM>
nCh = numel(chSave); assert(nCh == 2, 'expected 2 saved channels, found %d', nCh);
umPerPx = 0.3164 * (4 / zoom) * (256 / W);
fprintf('%dx%d px, zoom %g (%.4f um/px), %d planes (%d valid) x %d channels, channels saved %s\n', H, W, zoom, umPerPx, framesPerVol, nPlanesValid, nCh, mat2str(chSave));

%% 1. read tif once: z-sum of each channel (valid planes only), per-plane means
fpv = framesPerVol * nCh;
groups = arrayfun(@(c) c:nCh:(nPlanesValid * nCh), 1:nCh, 'UniformOutput', false);   % frame indices within a volume, channel c
o = pb_read_raw_tif_planes(tifPath, fpv, groups, 'planeMeans', true);
fprintf('%d volumes, IFD chain ok = %d\n', o.nVolumes, o.chainOk);

%% 2. one working folder per channel
for c = 1:nCh
    chDir = fullfile(trialDir, sprintf('ch%d', c)); if ~isfolder(chDir), mkdir(chDir); end
    imgData_sum = o.stacks{c};
    sp = struct('nplanes', framesPerVol, 'nplanes_valid', nPlanesValid, 'flyback_planes', nPlanesValid:(framesPerVol - 1), 'nchannels', 1, ...
        'channels_saved', chSave(c), 'px_height', H, 'px_width', W, 'um_per_px', [umPerPx umPerPx 4], 'zoom', zoom, 'nVolumes', o.nVolumes, ...
        'provenance', sprintf('overshoot_split_channels (ScanImage channel %d)', chSave(c)));
    save(fullfile(chDir, 'imgData_sum.mat'), 'imgData_sum', 'sp', '-v7.3');
    copyfile(fullfile(h5List(1).folder, h5List(1).name), chDir);
    fig = figure('Visible', 'off', 'Position', [50 50 1500 500], 'Color', 'w');
    pm = o.planeMean(:, :, c:nCh:(nPlanesValid * nCh));
    for k = 1:nPlanesValid, subplot(2, 5, k); imagesc(pm(:, :, k)); axis image off; colormap gray; title(sprintf('plane %d', k), 'FontSize', 8); end
    sgtitle(sprintf('ScanImage channel %d: per-plane mean image', chSave(c)))
    exportgraphics(fig, fullfile(trialDir, sprintf('planeMeans_ch%d.png', c)), 'Resolution', 110); close(fig);
    fprintf('saved %s (%s)\n', chDir, mat2str(size(imgData_sum)));
end
