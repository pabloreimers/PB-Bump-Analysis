function md = pb_scan_metadata(expDir)
%PB_SCAN_METADATA Read and sanity-check scan_params.json for a PB imaging trial.
%
%   md = pb_scan_metadata(expDir) reads <expDir>/scan_params.json (written by the
%   scopa acquisition tools) and the single raw ScanImage .tif in that folder, and
%   returns a struct describing how to carve the raw frame sequence into
%   [Y, X, plane, channel, volume]. Everything is cross-checked against the actual
%   .tif on disk (frame size, file size / frame count) so that a stale or
%   hand-edited json fails loudly here instead of silently corrupting a reshape
%   three functions downstream.
%
%   Fields of md:
%     expDir             the folder passed in
%     tifPath             full path to the raw .tif
%     px_height, px_width frame dimensions [Y, X]
%     nplanes             planes per volume, INCLUDING flyback
%     nchannels           channels acquired (channel is the fastest-varying index
%                          in the raw frame sequence, per ScanImage convention)
%     flyback_planes       1-indexed plane numbers that are flyback (junk)
%     valid_planes         1-indexed plane numbers that contain real data
%     nframes_total         raw IFD count actually used (derived from file size,
%                          cross-checked against scan_params.json)
%     n_volumes            complete volumes available (frames_total is truncated
%                          to a whole number of volumes; leftover frames from an
%                          interrupted acquisition are dropped, loudly)
%
%   See also PB_LOAD_RAW_TIF.

arguments
    expDir (1,:) char
end

if ~isfolder(expDir)
    error('pb_scan_metadata:missingDir', 'Experiment folder does not exist: %s', expDir);
end

jsonFile = fullfile(expDir, 'scan_params.json');
if ~isfile(jsonFile)
    error('pb_scan_metadata:missingJson', 'No scan_params.json found in %s', expDir);
end
sp = jsondecode(fileread(jsonFile));

tifList = dir(fullfile(expDir, '*.tif'));
tifList = tifList(~[tifList.isdir]);
if isempty(tifList)
    error('pb_scan_metadata:missingTif', 'No raw .tif found in %s', expDir);
elseif numel(tifList) > 1
    error('pb_scan_metadata:multipleTifs', ...
        'Expected exactly one raw .tif in %s, found %d:\n%s', ...
        expDir, numel(tifList), strjoin({tifList.name}, newline));
end

md = struct();
md.expDir          = expDir;
md.tifPath         = fullfile(tifList(1).folder, tifList(1).name);
md.px_height       = sp.px_height;
md.px_width        = sp.px_width;
md.nplanes         = sp.nplanes;                 % total planes/volume, incl. flyback
md.nchannels       = sp.nchannels;
md.channels_saved  = sp.channels_saved(:)';
md.fps             = sp.fps;
md.flyback_planes  = sort(sp.flyback_planes(:)') + 1;  % scopa writes 0-indexed
md.valid_planes    = setdiff(1:md.nplanes, md.flyback_planes);
md.nplanes_valid   = sp.nplanes_valid;

if numel(md.valid_planes) ~= md.nplanes_valid
    error('pb_scan_metadata:flybackMismatch', ...
        ['scan_params.json is internally inconsistent in %s: nplanes_valid=%d but ' ...
         'nplanes(%d) minus flyback_planes(%s, 1-indexed) leaves %d valid planes.'], ...
        expDir, md.nplanes_valid, md.nplanes, mat2str(md.flyback_planes), numel(md.valid_planes));
end

% cross-check frame size against the tif itself
t = Tiff(md.tifPath, 'r');
tifHeight = t.getTag('ImageLength');
tifWidth  = t.getTag('ImageWidth');
close(t);
if tifHeight ~= md.px_height || tifWidth ~= md.px_width
    error('pb_scan_metadata:dimMismatch', ...
        'scan_params.json says %d x %d px but %s frames are %d x %d.', ...
        md.px_height, md.px_width, md.tifPath, tifHeight, tifWidth);
end

% loose cross-check on frame count. Each IFD carries its own tags + an
% ImageDescription string on top of the raw pixel bytes, so file size is only
% an approximate lower bound on frame count, not an exact multiple -- trust
% scan_params.json's nframes_total (written by the acquisition software) and
% just flag it if the file size implies something wildly different (e.g. a
% truncated/corrupt file, or stale metadata copied from a different trial).
bytesPerFrame     = md.px_height * md.px_width * 2;
nFramesLowerBound = tifList(1).bytes / bytesPerFrame; % ignores per-IFD tag/description overhead
nFramesJson       = sp.nframes_total;
if nFramesLowerBound > nFramesJson * 1.02 || nFramesLowerBound < nFramesJson * 0.90
    warning('pb_scan_metadata:frameCountMismatch', ...
        ['%s: file size implies at least %.0f raw frames (ignoring per-frame TIFF overhead) but ' ...
         'scan_params.json says nframes_total=%d. These should be close (json count slightly below ' ...
         'the size-based bound); double check this scan_params.json actually belongs to this .tif.'], ...
        md.tifPath, nFramesLowerBound, nFramesJson);
end
md.nframes_total = nFramesJson;

framesPerVolume = md.nplanes * md.nchannels;
md.n_volumes    = floor(md.nframes_total / framesPerVolume);
leftover        = md.nframes_total - md.n_volumes * framesPerVolume;
if leftover > 0
    fprintf(['pb_scan_metadata: %s has %d leftover raw frame(s) (%.1f%% of one volume) after ' ...
             '%d complete volumes -- dropping the incomplete trailing volume.\n'], ...
             expDir, leftover, 100*leftover/framesPerVolume, md.n_volumes);
end
if md.n_volumes < 1
    error('pb_scan_metadata:noCompleteVolumes', ...
        '%s does not contain even one complete volume (%d raw frames, %d needed per volume).', ...
        expDir, md.nframes_total, framesPerVolume);
end

end
