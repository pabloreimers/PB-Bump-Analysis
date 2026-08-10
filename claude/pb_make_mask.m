function [mask, summaryImg, info] = pb_make_mask(movie, options)
%PB_MAKE_MASK Automatically segment the W-shaped protocerebral bridge (PB).
%
%   mask = pb_make_mask(movie) collapses a [Y,X,T] movie to a temporal
%   projection (the PB is wherever fluorescence lives), then thresholds and
%   cleans it up to return a logical [Y,X] mask that captures just the PB.
%
%   The pipeline (see also PB_READ_TIFF_MOVIE):
%     1. temporal projection over time    (default 'mean'; 'std' or 'max' too)
%     2. Gaussian smooth (sigma px)        bridge glomerular gaps so the W is
%                                          one connected piece
%     3. Otsu threshold (graythresh)       automatic, per-image -- no hand-tuned
%                                          constant that breaks across flies
%     4. fill holes
%     5. keep connected component(s) >= minAreaFrac of the frame (or, if none
%        qualify, the single largest) -- drops isolated bright noise specks
%     6. morphological close + fill        smooth the boundary
%
%   [mask, summaryImg, info] = pb_make_mask(...) also returns the projection
%   image used and an info struct (threshold, area, projection type).
%
%   You can also pass a precomputed 2-D summary image instead of a movie; it
%   is used directly as the projection.
%
%   Name-value options:
%     projection (1,:) char = 'mean'   'mean' | 'std' | 'max' temporal projection
%     sigma      (1,1) double = 1      Gaussian smoothing sigma in pixels
%     threshScale(1,1) double = 1      multiply the Otsu level (<1 grows mask,
%                                       >1 shrinks it) for manual fine-tuning
%     minAreaFrac(1,1) double = 0.02   min component area as a fraction of the
%                                       frame to be kept
%     closeRadius(1,1) double = 1      radius (px) of the disk used for closing
%
%   Requires the Image Processing Toolbox (graythresh, imbinarize, imfill,
%   bwconncomp, imclose, imgaussfilt).
%
%   Example:
%     load('...\denoised_srdtrans.mat','movie');   % [Y,X,T]
%     [mask, proj] = pb_make_mask(movie);
%     figure; imshowpair(mat2gray(proj), mask);    % overlay
%
%   See also PB_READ_TIFF_MOVIE, PB_MASK_ALL.

arguments
    movie
    options.projection (1,:) char {mustBeMember(options.projection,{'mean','std','max'})} = 'mean'
    options.sigma (1,1) double {mustBeNonnegative} = 1
    options.threshScale (1,1) double {mustBePositive} = 1
    options.minAreaFrac (1,1) double {mustBeInRange(options.minAreaFrac,0,1)} = 0.02
    options.closeRadius (1,1) double {mustBeNonnegative} = 1
end

% ---- 1. temporal projection ---------------------------------------------
if ndims(movie) == 3
    m = single(movie);
    switch options.projection
        case 'mean'; summaryImg = mean(m, 3);
        case 'std';  summaryImg = std(m, 0, 3);
        case 'max';  summaryImg = max(m, [], 3);
    end
elseif ismatrix(movie)
    summaryImg = single(movie);   % caller passed a precomputed projection
else
    error('pb_make_mask:badInput', 'movie must be [Y,X,T] or a 2-D image.');
end

% ---- 2. smooth ----------------------------------------------------------
if options.sigma > 0
    sm = imgaussfilt(summaryImg, options.sigma);
else
    sm = summaryImg;
end

% ---- 3. Otsu threshold (on a [0 1]-normalized copy) ---------------------
smN   = mat2gray(sm);
level = graythresh(smN) * options.threshScale;
bw    = imbinarize(smN, min(max(level, 0), 1));

% ---- 4. fill holes ------------------------------------------------------
bw = imfill(bw, 'holes');

% ---- 5. keep the large connected component(s) ---------------------------
cc = bwconncomp(bw);
if cc.NumObjects == 0
    warning('pb_make_mask:empty', 'Threshold produced an empty mask.');
    mask = false(size(bw));
else
    areas   = cellfun(@numel, cc.PixelIdxList);
    minArea = options.minAreaFrac * numel(bw);
    keep    = find(areas >= minArea);
    if isempty(keep)
        [~, keep] = max(areas);   % fall back to the single largest
    end
    mask = false(size(bw));
    mask(vertcat(cc.PixelIdxList{keep})) = true;
end

% ---- 6. close + fill ----------------------------------------------------
if options.closeRadius > 0
    mask = imclose(mask, strel('disk', round(options.closeRadius)));
end
mask = imfill(mask, 'holes');

% ---- info ---------------------------------------------------------------
info = struct('projection', options.projection, ...
              'otsuLevelNorm', level, ...
              'sigma', options.sigma, ...
              'threshScale', options.threshScale, ...
              'areaPx', nnz(mask), ...
              'areaFrac', nnz(mask)/numel(mask));

end
