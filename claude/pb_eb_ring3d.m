function [clusterIdx, alpha, info] = pb_eb_ring3d(mask3, nWedges, umPerPx, zPos_um, varargin)
% PB_EB_RING3D  Split a 3-D ellipsoid-body (EB) mask into N angular wedges
% around the ring, for ANY orientation of the ring in the imaging volume.
%
%   [clusterIdx, alpha, info] = pb_eb_ring3d(mask3, nWedges, umPerPx, zPos_um)
%   [...] = pb_eb_ring3d(..., 'plane', 'xy', 'weights', meanVol, 'direction', -1)
%
% mask3     : logical rows x cols x planes EB mask
% nWedges   : number of wedges (16 to match the PB's 16 glomeruli/hemisphere)
% umPerPx   : x-y pixel size (um)
% zPos_um   : 1 x planes, z position of each plane (um); uneven spacing ok
% 'plane'   : 'auto' (default) = fit the ring plane by PCA of the mask voxels
%             in physical (um) coordinates; or force it: 'xy' (ring seen
%             face-on), 'xz' (edge-on, ring spans columns x depth -- the
%             PB+EB prep, same geometry as pb_eb_wedges3d), 'yz'
% 'weights' : rows x cols x planes intensity (e.g. the mean volume) to weight
%             the PCA / centre and to build the ring-plane QC image (default
%             [] = unweighted, QC image = voxel counts)
% 'direction': +1 (default) or -1, sense of wedge numbering around the ring
%             (resolved downstream by the bump-vs-bar sign test, as for the PB)
% 'theta0'  : in-plane angle (rad, 0 = +e1, increasing towards +e2) of the
%             centre of wedge 1; default -pi
%
% Why a plane fit: pb_eb_wedges (2-D) assumes the ring lies in x-y and
% pb_eb_wedges3d assumes it lies in x-z. Which one holds depends on the
% prep and head tilt, and a tilted ring under either assumption gives
% wedges of very unequal arc length. Here every mask voxel is placed in um
% (x = col, y = row, z = plane depth), the two largest principal axes span
% the ring plane (e1, e2) and the smallest is its normal (e3). Coarse z
% sampling (a few um per plane vs ~0.3-0.6 um/px) is handled because the
% PCA runs in um, not voxel units. A filled-in or thick EB still has a
% well-defined plane as long as its diameter clearly exceeds its thickness;
% info.planarity (1 - lambda3/lambda2) reports this, and a warning is
% issued below 0.4 -- then force 'plane' from the QC projections.
%
% In the plane, voxels are rasterised at max(umPerPx, plane spacing / 2)
% into a footprint image;
% the ring centre is the centroid of the footprint's outer boundary (so an
% uneven ring or a bright filled centre does not bias it, as in
% pb_eb_wedges), the semi-axes are the 98th percentiles of |u|, |v| about
% that centre, and each voxel gets angle atan2(v/b, u/a) (ellipse-
% normalised, so an oval ring gets equal-arc wedges) and the nearest wedge.
% Axis signs are made deterministic (largest component of e1 and e3
% positive, e2 = e3 x e1), so the numbering is reproducible across runs.
%
% clusterIdx: rows x cols x planes, 0 outside the mask, 1..N inside
% alpha     : 1 x N, alpha(k) = -pi + (k-1) * 2pi/N (PVA angle of wedge k)
% info      : plane ('auto'/'xy'/...), axes [e1 e2 e3] (3x3, columns, in
%             [x y z]), eigvals, planarity, center_um ([x y z] um), center_vox
%             ([col row plane]), semiAxes_um [a b], hasHole, nVox (1 x N),
%             angleMap (NaN outside), and ring-plane QC images: footprint,
%             planeImg (mean weight per bin), planeLabel (modal wedge per
%             bin), uGrid / vGrid (bin centres, um, relative to the centre)
%
% The PVA then follows the PB recipe exactly: unit vector at alpha(k)
% weighted by the wedge's z-scored dF, summed over the N wedges.

p = inputParser;
p.addParameter('plane', 'auto', @(s) any(strcmpi(s, {'auto', 'xy', 'xz', 'yz'})));
p.addParameter('weights', [], @(v) isempty(v) || isnumeric(v));
p.addParameter('direction', 1, @(v) any(v == [-1 1]));
p.addParameter('theta0', -pi, @isnumeric);
p.parse(varargin{:});
mask3 = logical(mask3); [nR, nC, nP] = size(mask3);
zPos_um = double(zPos_um(:)); assert(numel(zPos_um) == nP, 'pb_eb_ring3d: zPos_um has %d entries, mask has %d planes', numel(zPos_um), nP);
assert(nnz(mask3) >= nWedges, 'pb_eb_ring3d: mask has only %d voxels', nnz(mask3));

idx = find(mask3); [rr, cc, pp] = ind2sub(size(mask3), idx);
P = [cc * umPerPx, rr * umPerPx, zPos_um(pp)];                  % [x y z] um
wts = p.Results.weights;
if isempty(wts), w = ones(numel(idx), 1); else, w = double(wts(idx)); w = max(w - min(w), 0) + eps; end
w = w / sum(w);
mu = sum(P .* w, 1); X = P - mu;

% ---- ring plane ----
lam = nan(1, 3);
switch lower(p.Results.plane)
    case 'auto'
        C = X' * (X .* w); [V, D] = eig((C + C') / 2); [lam, o] = sort(diag(D)', 'descend'); E = V(:, o);
    case 'xy', E = [1 0 0; 0 1 0; 0 0 1]';
    case 'xz', E = [1 0 0; 0 0 1; 0 1 0]';
    case 'yz', E = [0 1 0; 0 0 1; 1 0 0]';
end
e1 = E(:, 1); e3 = E(:, 3);
[~, m] = max(abs(e1)); e1 = e1 * sign(e1(m));
[~, m] = max(abs(e3)); e3 = e3 * sign(e3(m));
e2 = cross(e3, e1);
if strcmpi(p.Results.plane, 'auto')
    planarity = 1 - lam(3) / lam(2);
    if planarity < 0.4
        warning('pb_eb_ring3d:flat', 'ring plane poorly defined (planarity %.2f, eigenvalues %s um^2) -- check the QC projections and force ''plane'' if needed.', planarity, mat2str(round(lam, 1)));
    end
else
    planarity = NaN;
end
u = X * e1; v = X * e2;

% ---- centre + semi-axes in the ring plane (footprint outer boundary) ----
% bin size: half the plane spacing (or the pixel size if larger) -- planes
% project into the ring plane as stripes up to dz apart (edge-on ring), and
% a finer raster leaves gaps that split the footprint and wreck the centre
dz = 0; if nP > 1, dz = median(abs(diff(zPos_um))); end
res = max(umPerPx, dz / 2); u0 = min(u); v0 = min(v);
iu = round((u - u0) / res) + 1; iv = round((v - v0) / res) + 1; sz = [max(iv) max(iu)];
F = false(sz); F(sub2ind(sz, iv, iu)) = true;
F = imclose(F, strel('disk', 1));
B = bwboundaries(F, 'noholes');
if ~isempty(B)
    [~, iB] = max(cellfun(@numel, B)); b = B{iB};
    uc = (mean(b(:, 2)) - 1) * res + u0; vc = (mean(b(:, 1)) - 1) * res + v0;
else
    uc = 0; vc = 0;
end
Ff = imfill(F, 'holes'); hasHole = nnz(Ff) > 1.05 * nnz(F);
a = max(prctile(abs(u - uc), 98), res); bb = max(prctile(abs(v - vc), 98), res);

% ---- wedges ----
th = atan2((v - vc) / bb, (u - uc) / a);
rel = p.Results.direction * (th - p.Results.theta0);
k = mod(round(rel / (2*pi/nWedges)), nWedges) + 1;
clusterIdx = zeros(nR, nC, nP); clusterIdx(idx) = k;
alpha = linspace(-pi, pi - 2*pi/nWedges, nWedges);
angMap = nan(nR, nC, nP); angMap(idx) = th;

% ---- ring-plane QC images ----
lin = sub2ind(sz, iv, iu);
if isempty(wts), wv = ones(numel(idx), 1); else, wv = double(wts(idx)); end
planeImg = accumarray(lin, wv, [prod(sz) 1], @mean, NaN); planeImg = reshape(planeImg, sz);
planeLabel = accumarray(lin, k, [prod(sz) 1], @mode, 0); planeLabel = reshape(planeLabel, sz);

c3 = mu + uc * e1' + vc * e2';
info.plane = lower(p.Results.plane); info.axes = [e1 e2 e3]; info.eigvals = lam; info.planarity = planarity;
info.center_um = c3; info.center_vox = [c3(1) / umPerPx, c3(2) / umPerPx, interp1(zPos_um, 1:nP, c3(3), 'linear', 'extrap')];
info.semiAxes_um = [a bb]; info.hasHole = hasHole; info.nVox = accumarray(k, 1, [nWedges 1])';
info.angleMap = angMap; info.footprint = F; info.planeImg = planeImg; info.planeLabel = planeLabel;
info.uGrid = ((1:sz(2)) - 1) * res + u0 - uc; info.vGrid = ((1:sz(1)) - 1) * res + v0 - vc;
info.direction = p.Results.direction; info.theta0 = p.Results.theta0;
end
