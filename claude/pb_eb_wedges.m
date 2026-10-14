function [clusterIdx, alpha, info] = pb_eb_wedges(mask, nWedges, varargin)
% PB_EB_WEDGES  Split an ellipsoid-body (EB) ring mask into angular wedges.
%
%   [clusterIdx, alpha, info] = pb_eb_wedges(mask, nWedges)
%   [...] = pb_eb_wedges(mask, nWedges, 'center', [x y], 'theta0', th, 'direction', -1)
%
% mask      : logical EB mask (ring or filled disc)
% nWedges   : number of wedges (16 to match the PB's 16 glomeruli/hemisphere)
% 'center'  : [x y] ring centre in pixels; default = centroid of the mask,
%             refined by one pass of the mean of the mask's boundary points
%             (so a filled centre or an uneven ring does not bias it)
% 'theta0'  : angle (rad, image coordinates: 0 = +x, increasing towards +y,
%             i.e. clockwise on screen) of the centre of wedge 1; default -pi
%             so wedge k has alpha = -pi + (k-1) * 2pi/N, matching the PB
%             convention alpha = linspace(-pi, pi - 2pi/N, N)
% 'direction': +1 (default) wedges numbered with increasing image angle,
%             -1 for the opposite sense (use to match the PB's glomerulus
%             order; the sign is otherwise resolved downstream by the
%             bump-vs-bar sign test, as for the PB)
%
% clusterIdx: same size as mask, 0 outside the mask, 1..N wedge id inside
% alpha     : 1 x N wedge angles, alpha(k) = -pi + (k-1) * 2pi/N  (the PVA
%             angle each wedge votes for, regardless of 'theta0'/'direction')
% info      : struct with center, theta0, direction, nPix per wedge, and
%             the mean radius of the mask
%
% The PVA then follows exactly the PB recipe: each wedge is a unit vector at
% alpha(k) weighted by its z-scored dF, summed over the N wedges. Unlike the
% PB there is only one copy of the ring, so no hemisphere averaging applies.

p = inputParser;
p.addParameter('center', [], @(v) isempty(v) || numel(v) == 2);
p.addParameter('theta0', -pi, @isnumeric);
p.addParameter('direction', 1, @(v) any(v == [-1 1]));
p.parse(varargin{:});
mask = logical(mask);
[yy, xx] = find(mask);
if isempty(p.Results.center)
    c0 = [mean(xx) mean(yy)];
    B = bwboundaries(mask, 'noholes');
    if ~isempty(B)
        [~, iB] = max(cellfun(@numel, B)); b = B{iB};           % outer boundary (largest)
        c0 = [mean(b(:,2)) mean(b(:,1))];
    end
    center = c0;
else
    center = p.Results.center(:)';
end
th = atan2(yy - center(2), xx - center(1));                   % image-coordinate angle of every mask pixel
rel = p.Results.direction * (th - p.Results.theta0);          % angle relative to wedge-1 centre, signed by direction
k = mod(round(rel / (2*pi/nWedges)), nWedges) + 1;            % nearest wedge centre
clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), yy, xx)) = k;
alpha = linspace(-pi, pi - 2*pi/nWedges, nWedges);
info.center = center; info.theta0 = p.Results.theta0; info.direction = p.Results.direction;
info.nPix = accumarray(k, 1, [nWedges 1])';
info.meanRadius = mean(hypot(xx - center(1), yy - center(2)));
end
