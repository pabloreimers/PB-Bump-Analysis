function [movie, md] = pb_load_raw_tif(expDir, options)
%PB_LOAD_RAW_TIF Load raw ScanImage frames for a PB imaging trial into [Y,X,plane,T].
%
%   [movie, md] = pb_load_raw_tif(expDir) loads every volume of channel 1 for the
%   trial in expDir (a folder containing scan_params.json and one raw .tif, e.g.
%   Z:\pablo\lpsp_p2x2_walking\<date>\fly <n>\<trial>), keeping only the valid
%   (non-flyback) planes, and returns movie sized [Y, X, nplanes_valid, T].
%
%   [movie, md] = pb_load_raw_tif(expDir, 'channel', 2, 'volumes', [1 50]) loads
%   only the first 50 volumes of channel 2 -- use a small 'volumes' range while
%   developing/debugging so you are not waiting on a multi-GB network read every
%   time you tweak a filter parameter.
%
%   Name-value options:
%     channel   (1,1) double = 1        which saved channel to return (not a list;
%                                        call twice if you need both channels)
%     volumes   (1,2) double = [1 Inf]  inclusive [first last] volume to load;
%                                        Inf means "through the last complete volume"
%     dropFlyback (1,1) logical = true  if false, all md.nplanes planes are kept
%                                        (in raw plane order) instead of just
%                                        md.valid_planes
%
%   Frames are read sequentially (Tiff.nextDirectory), which is dramatically
%   faster than ScanImageTiffReader for a partial read of a large file (a few
%   hundred frames/sec vs. ~2 minutes just to open an 8+ GB file). Both channels
%   are necessarily read off disk for the requested volumes (channel is the
%   fastest-varying index in the raw frame sequence and frames can only be read
%   in order), but only the requested channel is kept in memory.
%
%   See also PB_SCAN_METADATA.

arguments
    expDir (1,:) char
    options.channel (1,1) double {mustBeMember(options.channel,[1 2 3 4])} = 1
    options.volumes (1,2) double = [1 Inf]
    options.dropFlyback (1,1) logical = true
end

md = pb_scan_metadata(expDir);

if options.channel > md.nchannels
    error('pb_load_raw_tif:badChannel', ...
        '%s only has %d channel(s) saved, cannot load channel %d.', expDir, md.nchannels, options.channel);
end

firstVol = options.volumes(1);
lastVol  = min(options.volumes(2), md.n_volumes);
if ~isfinite(lastVol); lastVol = md.n_volumes; end
if firstVol < 1 || firstVol > lastVol
    error('pb_load_raw_tif:badVolumeRange', ...
        'Requested volumes [%g %g] is not a valid, non-empty range within 1:%d.', ...
        options.volumes(1), options.volumes(2), md.n_volumes);
end
nVolLoad = lastVol - firstVol + 1;

framesPerVolume = md.nplanes * md.nchannels;
firstIFD = (firstVol - 1) * framesPerVolume + 1;
nIFDtoRead = nVolLoad * framesPerVolume;

fprintf('pb_load_raw_tif: %s\n', md.tifPath);
fprintf('  loading volumes %d:%d of %d (channel %d of %d, %d raw frames)\n', ...
    firstVol, lastVol, md.n_volumes, options.channel, md.nchannels, nIFDtoRead);

raw = zeros(md.px_height, md.px_width, nIFDtoRead, 'int16');
t = Tiff(md.tifPath, 'r');
cleanupTiff = onCleanup(@() close(t));
t.setDirectory(firstIFD);
tic
for i = 1:nIFDtoRead
    raw(:,:,i) = t.read();
    if i < nIFDtoRead
        t.nextDirectory();
    end
    if mod(i, 20000) == 0 || i == nIFDtoRead
        elapsed = toc;
        fprintf('  read %d/%d raw frames (%.1f s, %.0f frames/s)\n', i, nIFDtoRead, elapsed, i/elapsed);
    end
end

raw5 = reshape(raw, md.px_height, md.px_width, md.nchannels, md.nplanes, nVolLoad);
chanMovie = squeeze(raw5(:, :, options.channel, :, :)); % Y x X x plane x volume

if options.dropFlyback
    chanMovie = chanMovie(:, :, md.valid_planes, :);
end

movie = chanMovie;
md.volumesLoaded = [firstVol, lastVol];
md.channelLoaded = options.channel;

end
