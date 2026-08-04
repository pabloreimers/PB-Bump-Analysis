function [movie, info] = pb_read_tiff_movie(tifPath, options)
%PB_READ_TIFF_MOVIE Read a multi-page TIFF movie into a [Y,X,T] array.
%
%   movie = pb_read_tiff_movie(tifPath) reads every page of the TIFF at
%   tifPath (a single-channel movie, one frame per page) and returns it as a
%   3-D array sized [Y, X, T] in the file's native class (e.g. single/uint16).
%
%   movie = pb_read_tiff_movie(tifPath, 'frames', [1 200]) reads only the
%   inclusive [first last] frame range -- Inf as the last frame means
%   "through the final page".
%
%   [movie, info] = pb_read_tiff_movie(...) also returns the imfinfo struct
%   array for the file.
%
%   If tifPath is a folder, the function looks for a single .tif/.tiff inside
%   it; if there is exactly one it is used, otherwise an error lists the
%   candidates. This lets you pass the denoised_ch0 folder directly.
%
%   Example:
%     d = ['Z:\pablo\epg_gain_change\20250624\fly 1\' ...
%          '20250624-1_epg_syt8m_gain_change_8\zproj_ds_srdt_asinh_s001\denoised_ch0'];
%     movie = pb_read_tiff_movie(fullfile(d, 'denoised_srdtrans.tif'));
%
%   See also PB_PLAY_MOVIE, PB_LOAD_RAW_TIF.

arguments
    tifPath (1,:) char
    options.frames (1,2) double = [1 Inf]
end

% Allow passing a folder that contains a single tiff.
if isfolder(tifPath)
    hits = [dir(fullfile(tifPath, '*.tif')); dir(fullfile(tifPath, '*.tiff'))];
    if numel(hits) == 1
        tifPath = fullfile(hits.folder, hits.name);
    elseif isempty(hits)
        error('pb_read_tiff_movie:noTiff', 'No .tif/.tiff found in folder %s.', tifPath);
    else
        error('pb_read_tiff_movie:manyTiff', ...
            '%d tiffs in %s; pass the full file path. Found: %s', ...
            numel(hits), tifPath, strjoin({hits.name}, ', '));
    end
end
if ~isfile(tifPath)
    error('pb_read_tiff_movie:notFound', 'File not found: %s', tifPath);
end

info    = imfinfo(tifPath);
nPages  = numel(info);
Y       = info(1).Height;
X       = info(1).Width;

firstFr = max(1, options.frames(1));
lastFr  = min(options.frames(2), nPages);
if ~isfinite(lastFr); lastFr = nPages; end
if firstFr > lastFr
    error('pb_read_tiff_movie:badFrameRange', ...
        'Requested frames [%g %g] is not a valid range within 1:%d.', ...
        options.frames(1), options.frames(2), nPages);
end
nFrames = lastFr - firstFr + 1;

fprintf('pb_read_tiff_movie: %s\n', tifPath);
fprintf('  %d x %d, reading frames %d:%d of %d\n', Y, X, firstFr, lastFr, nPages);

t = Tiff(tifPath, 'r');
cleanupTiff = onCleanup(@() close(t));

t.setDirectory(firstFr);
firstImg = t.read();
movie = zeros(Y, X, nFrames, 'like', firstImg);
movie(:, :, 1) = firstImg;

tic
for i = 2:nFrames
    t.nextDirectory();
    movie(:, :, i) = t.read();
    if mod(i, 500) == 0 || i == nFrames
        elapsed = toc;
        fprintf('  read %d/%d frames (%.1f s, %.0f frames/s)\n', i, nFrames, elapsed, i/elapsed);
    end
end

end
