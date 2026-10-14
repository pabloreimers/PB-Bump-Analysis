function [clusterIdx, alpha, info] = pb_eb_wedges3d(mask3, nWedges, varargin)
% PB_EB_WEDGES3D  Split a 3-D ellipsoid-body (EB) mask into angular wedges
% when the ring is imaged EDGE-ON, i.e. the ring plane is spanned by the
% image x axis and the z (plane) axis, not by x and y.
%
%   [clusterIdx, alpha, info] = pb_eb_wedges3d(mask3, nWedges)
%   [...] = pb_eb_wedges3d(mask3, nWedges, 'zScale', 5/0.84, 'direction', -1)
%
% mask3     : logical rows x cols x planes EB mask
% nWedges   : number of wedges (16 to match the PB's 16 glomeruli/hemisphere)
% 'zScale'  : plane spacing in units of xy pixels (um/plane / um/px); used
%             only for the reported radii -- the angle itself is computed in
%             ellipse-normalised coordinates, so it is scale-free
% 'direction': +1 (default) or -1, sense of wedge numbering around the ring
%             (resolved downstream by the bump-vs-bar sign test, as for the PB)
% 'theta0'  : angle (rad, in the normalised x-z plane: 0 = +x, increasing
%             towards +z / deeper planes) of the centre of wedge 1; default -pi
%
% Geometry: the mask is projected onto the x-z plane (any row), the ring
% centre is the centroid of that footprint's outer boundary, and the semi-
% axes are the footprint's half-extents in x and z. Every voxel (r, c, p) in
% the mask gets angle atan2((p - zc)/bz, (c - xc)/ax) and is assigned to the
% nearest wedge centre. Rows (the ring's thickness, ~10 px here) are pooled.
%
% clusterIdx: rows x cols x planes, 0 outside the mask, 1..N inside
% alpha     : 1 x N, alpha(k) = -pi + (k-1) * 2pi/N (PVA angle of wedge k)
% info      : center [xc zc] (col, plane), semiAxes [ax bz] (px, planes),
%             nVox per wedge, footprint (planes x cols logical), the per-voxel
%             angle map (NaN outside), and the x-z wedge label map used for QC
%
% The PVA then follows the PB recipe exactly: unit vector at alpha(k)
% weighted by the wedge's z-scored dF, summed over the N wedges.

p = inputParser;
p.addParameter('zScale', 1, @isnumeric);
p.addParameter('direction', 1, @(v) any(v == [-1 1]));
p.addParameter('theta0', -pi, @isnumeric);
p.parse(varargin{:});
mask3 = logical(mask3);
[nR, nC, nP] = size(mask3);
foot = squeeze(any(mask3, 1));                     % cols x planes
if size(foot, 1) ~= nC, foot = foot'; end          % make it cols x planes robustly
footPZ = foot';                                     % planes x cols (image convention: row = plane)
B = bwboundaries(footPZ, 'noholes');
if ~isempty(B)
    [~, iB] = max(cellfun(@numel, B)); b = B{iB};  % [plane col]
    center = [mean(b(:,2)) mean(b(:,1))];           % [xc zc]
else
    [pz, cx] = find(footPZ); center = [mean(cx) mean(pz)];
end
[pz, cx] = find(footPZ);
semi = [max(cx) - min(cx), max(pz) - min(pz)] / 2; semi = max(semi, 0.5);
[rr, cc, pp] = ind2sub(size(mask3), find(mask3)); %#ok<ASGLU>
th = atan2((pp - center(2)) / semi(2), (cc - center(1)) / semi(1));
rel = p.Results.direction * (th - p.Results.theta0);
k = mod(round(rel / (2*pi/nWedges)), nWedges) + 1;
clusterIdx = zeros(nR, nC, nP);
clusterIdx(mask3) = k;
alpha = linspace(-pi, pi - 2*pi/nWedges, nWedges);
angMap = nan(nR, nC, nP); angMap(mask3) = th;
% x-z label map for QC (mode over rows)
xzLabel = zeros(nP, nC);
for pi_ = 1:nP
    for ci = 1:nC
        v = clusterIdx(:, ci, pi_); v = v(v > 0);
        if ~isempty(v), xzLabel(pi_, ci) = mode(v); end
    end
end
info.center = center; info.semiAxes = semi; info.nVox = accumarray(k, 1, [nWedges 1])';
info.footprint = footPZ; info.angleMap = angMap; info.xzLabel = xzLabel;
info.radius_px = [semi(1), semi(2) * p.Results.zScale];
info.direction = p.Results.direction; info.theta0 = p.Results.theta0;
end
