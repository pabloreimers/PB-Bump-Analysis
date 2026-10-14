function [mask, maskCore, info] = pb_arch_mask(projMean, umPerPx, varargin)
% PB_ARCH_MASK  Mask the protocerebral bridge (a single thin arch) in a mean
% image, robust to a brightness gradient along the bridge and to bright
% non-PB structures elsewhere in the frame.
%
%   [mask, maskCore, info] = pb_arch_mask(projMean, umPerPx)
%   [...] = pb_arch_mask(projMean, umPerPx, 'thresh_pct', 88, 'sigma_um', 1.5, ...)
%
% Steps (all sizes in um, converted with umPerPx):
%   1. smooth (sigma_um), clip to the 5th-99th percentile.
%   2. candidate pixels = above the global percentile threshold (thresh_pct)
%      OR above local_frac of the local column maximum (window local_win_um)
%      while above the local_floor_pct percentile. The column-normalised
%      criterion keeps a hemisphere that is several times dimmer than the
%      other (20261012-3: left arm 1/3 the brightness of the right).
%   3. open (open_um), fill; connected components. The largest is the main
%      arch; its skeleton gives a half-width (area / 2 length). Other
%      components are kept only if (a) they are >= hemi_frac of the largest
%      and lie in the same row band (a separate hemisphere), or (b) they
%      are >= min_seed_frac of the largest and within dist_factor half-
%      widths of the main skeleton (a split-off tip glomerulus). Anything
%      else -- e.g. the bright blob below the bridge in 20261012-3 -- is
%      dropped. Previous versions bridged such blobs into the mask with a
%      disk closing, which then wrecked the skeleton-based glomerulus order.
%   4. gaps between kept components are filled from a looser threshold
%      (thresh_pct - loose_drop) but only inside a HORIZONTAL closing hull of
%      the seeds (join_um), i.e. along the arch, never towards structures
%      above or below it. If two hemispheres are still disconnected they are
%      bridged with the same horizontal line closing (the bridge pixels are
%      not in maskCore).
%
% mask     : logical, the (bridged) arch; maskCore: pixels with real signal
% info     : struct with the component table (area, kept, reason), the
%            half-width, thresholds and the skeleton of the main component

p = inputParser;
p.addParameter('sigma_um', 1.5); p.addParameter('thresh_pct', 88); p.addParameter('loose_drop', 15);
p.addParameter('open_um', 1.0); p.addParameter('join_um', 15);
p.addParameter('local_frac', 0.45); p.addParameter('local_floor_pct', 75); p.addParameter('local_win_um', 12);
p.addParameter('min_seed_frac', 0.05); p.addParameter('hemi_frac', 0.3); p.addParameter('dist_factor', 1.5);
p.addParameter('bright_frac', 0.2);   % a kept component / local-criterion pixel must be at least this fraction of the main arch's median brightness
p.parse(varargin{:}); o = p.Results;

sigma_px = max(o.sigma_um / umPerPx, 0.5); open_px = max(round(o.open_um / umPerPx), 1);
join_px = max(round(o.join_um / umPerPx), 3); win_px = max(round(o.local_win_um / umPerPx), 3);
ps = imgaussfilt(double(projMean), sigma_px);
pc = min(max(ps, prctile(ps, 5, 'all')), prctile(ps, 99, 'all'));
thrHi = prctile(pc(:), o.thresh_pct); thrLo = prctile(pc(:), o.thresh_pct - o.loose_drop); thrFloor = prctile(pc(:), o.local_floor_pct);
% reference brightness = the main arch at the global threshold, so the local
% criterion and the component checks scale with this fly's signal
hiAbs = imopen(pc > thrHi, strel('disk', open_px));
ccA = bwconncomp(hiAbs); [~, iA] = max(cellfun(@numel, ccA.PixelIdxList));
refBright = median(pc(ccA.PixelIdxList{iA}));
thrFloor = max(thrFloor, o.bright_frac * refBright);
colmax = movmax(max(pc, [], 1), win_px);
loc = pc ./ max(colmax, eps);
hi = (pc > thrHi) | (loc > o.local_frac & pc > thrFloor);
hi = imfill(imopen(hi, strel('disk', open_px)), 'holes');

cc = bwconncomp(hi); areas = cellfun(@numel, cc.PixelIdxList);
assert(~isempty(areas), 'pb_arch_mask: nothing above threshold');
[~, iMain] = max(areas);
main = false(size(hi)); main(cc.PixelIdxList{iMain}) = true;
skel = bwskel(bwmorph(main, 'majority'));
if nnz(skel) < 3, skel = bwskel(main); end
halfWidth = max(nnz(main) / (2 * max(nnz(skel), 1)), 1.5);
D = bwdist(skel);
[rMain, ~] = find(main); band = [min(rMain) - 2 * halfWidth, max(rMain) + 2 * halfWidth];
keep = false(1, cc.NumObjects); reason = repmat({'dropped'}, 1, cc.NumObjects); keep(iMain) = true; reason{iMain} = 'main arch';
for c = setdiff(1:cc.NumObjects, iMain)
    if areas(c) < o.min_seed_frac * areas(iMain), reason{c} = 'too small'; continue; end
    [rc, ~] = ind2sub(size(hi), cc.PixelIdxList{c}); dmin = min(D(cc.PixelIdxList{c}));
    bright = median(pc(cc.PixelIdxList{c})) / refBright;
    if areas(c) >= o.hemi_frac * areas(iMain) && mean(rc) >= band(1) && mean(rc) <= band(2)
        keep(c) = true; reason{c} = 'second hemisphere';
    elseif dmin <= o.dist_factor * halfWidth && bright >= o.bright_frac
        keep(c) = true; reason{c} = sprintf('near arch (%.1f px, %.0f%% bright)', dmin, 100 * bright);
    elseif dmin <= o.dist_factor * halfWidth
        reason{c} = sprintf('near but dim (%.1f px, %.0f%% bright)', dmin, 100 * bright);
    else
        reason{c} = sprintf('far from arch (%.1f px)', dmin);
    end
end
seeds = false(size(hi)); for c = find(keep), seeds(cc.PixelIdxList{c}) = true; end

lo = pc > thrLo;
hull = (imclose(seeds, strel('line', join_px, 0)) | imclose(seeds, strel('disk', 2))) & (bwdist(seeds) <= halfWidth);
gap = imdilate(lo & hull & ~seeds, strel('disk', 1)) & hull;
mask = imfill(seeds | gap, 'holes');
maskCore = mask;
cc2 = bwconncomp(mask);
if cc2.NumObjects > 1
    a2 = cellfun(@numel, cc2.PixelIdxList); big = find(a2 >= o.hemi_frac * max(a2));
    two = false(size(mask)); for c = big, two(cc2.PixelIdxList{c}) = true; end
    maskCore = two;
    mask = imfill(imclose(two, strel('line', join_px, 0)), 'holes');
end
mask = bwareafilt(mask, 1); maskCore = maskCore & mask;

info.areas = areas; info.keep = keep; info.reason = reason; info.halfWidth_px = halfWidth; info.skelMain = skel;
info.thresholds = [thrHi thrLo thrFloor]; info.nPix = nnz(mask); info.nCorePix = nnz(maskCore);
info.summary = sprintf('%d component(s): %s; half-width %.1f px; %d px (core %d)', cc.NumObjects, ...
    strjoin(arrayfun(@(c) sprintf('%d px %s', areas(c), reason{c}), 1:cc.NumObjects, 'UniformOutput', false), ', '), halfWidth, nnz(mask), nnz(maskCore));
end
