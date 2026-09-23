function clusterIdx = pb_skeleton_glomeruli(mask, nPerHemisphere)
% PB_SKELETON_GLOMERULI  Divide a PB arch mask into 2*nPerHemisphere equal "glomeruli".
%
%   clusterIdx = pb_skeleton_glomeruli(mask, nPerHemisphere)
%
% Same skeletonize -> order into a path (graph_sort.m, repo root) -> resample
% into evenly spaced centroids -> nearest-centroid-assign approach as this
% lab's process_im (e.g. epg_dlight_script.m, dopamine_ionto_script.m),
% trimmed to just the clustering step. Assumes mask is a single open arch,
% not a closed loop (graph_sort needs an endpoint to start from).
%
% clusterIdx: same [Y,X] size as mask; 0 outside the mask, 1..2*nPerHemisphere
% inside, numbered along the skeleton from one end of the arch to the other.
% Which end is #1 depends on which skeleton endpoint bwmorph finds first, so
% the numbering is consistent within a fly (one mask) but not guaranteed to
% start on the same side across flies.

if isempty(which('graph_sort'))
    repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
    if isfile(fullfile(repoRoot, 'graph_sort.m'))
        addpath(repoRoot);
    else
        error('pb_skeleton_glomeruli:graphSort', 'graph_sort.m (repo root) is not on the path.');
    end
end

[y_mask, x_mask] = find(mask);
min_axis = min(range(x_mask), range(y_mask));

maskClean = bwmorph(mask, 'majority');
skel      = bwskel(maskClean);
epIdx     = find(bwmorph(skel, 'endpoints'), 1);
if isempty(epIdx)
    error('pb_skeleton_glomeruli:noEndpoint', 'PB mask skeleton has no endpoint -- expected a single open arch, not a closed loop.');
end
D    = bwdistgeodesic(skel, epIdx);
[~, idxFar] = max(D(:));
D2   = bwdistgeodesic(skel, idxFar);
mid  = (D + D2) == mode(D + D2, 'all');

[y_mid, x_mid] = find(mid);
[x_mid, y_mid] = graph_sort(x_mid, y_mid);

xq    = -min_axis:(length(x_mid) + min_axis);
x_mid = round(interp1(1:length(x_mid), x_mid, xq, 'linear', 'extrap'));
y_mid = round(interp1(1:length(y_mid), y_mid, xq, 'linear', 'extrap'));

keep  = ismember([x_mid', y_mid'], [x_mask, y_mask], 'rows');
x_mid = x_mid(keep);
y_mid = y_mid(keep);

nClusters = 2 * nPerHemisphere;
xq2 = linspace(1, length(y_mid), 2*nClusters + 1)';
centroids = [interp1(1:length(y_mid), y_mid, xq2), interp1(1:length(x_mid), x_mid, xq2)];
centroids = centroids(2:2:end-1, :); % drop the edge samples so every cluster is the same size

[~, idx] = pdist2(centroids, [y_mask, x_mask], 'euclidean', 'smallest', 1);

clusterIdx = zeros(size(mask));
clusterIdx(sub2ind(size(mask), y_mask, x_mask)) = idx;
end
