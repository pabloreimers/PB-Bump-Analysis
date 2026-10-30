function rm = pb_rotation_map(movie, yaw, fps, options)
%PB_ROTATION_MAP Find the PB as the pixels whose photon rate follows turning.
%
%   rm = pb_rotation_map(movie, yaw, fps) correlates every pixel of a
%   registered photon movie (Y x X x T) with a turning regressor built from
%   yaw (T x 1, rad/s, sampled at the volumes) and thresholds that map into a
%   PB mask. Use this when the PB is too dim to see in the mean image -- on
%   20260216-1 (epg>dLight) the PB is only ~1.1-1.4x its surround in the
%   whole-trial mean, but in this map it is the only structure in the frame.
%
%   The regressor is |yaw| convolved with the measured |yaw| -> PB photon
%   response, not raw |yaw|: dLight rises ~0.3 s after a turn starts and stays
%   up for seconds, so the lag-0 |yaw| map is much weaker (peak corr ~0.15 vs
%   ~0.3 at sigma 1.5 px). The kernel is estimated by ridge regression of the
%   summed photons in a provisional mask on lagged |yaw| (lags 'lagsSec'),
%   where the provisional mask comes from a first pass with a generic
%   regressor (|yaw| lagged 0.5 s, 2 s moving mean) -- or from 'refMask'
%   (e.g. a hand-drawn mask.mat) if you pass one. Only the causal part of the
%   kernel (lag >= 0) is used to build the regressor.
%
%   Caveat: because the mask is chosen for turning-related signal, anything
%   in the PB that does NOT co-vary with turning can be under-represented at
%   the mask border. pb_lowrank_denoise itself is unsupervised (does not see
%   yaw) so it is not biased toward turning inside the mask.
%
%   Name-value options:
%     sigma      (1,1) double = 1.5   px, spatial pooling before correlating
%     lagsSec    (1,2) double = [-1 4] kernel lag range, s (+ = dLight after yaw)
%     ridge      (1,1) double = 10    ridge penalty, x mean regressor energy
%     refMask    (:,:) logical = []   provisional mask (skip the first pass)
%     maskFrac   (1,1) double = 0.3   mask = map > maskFrac * max(map)
%     minBlob    (1,1) double = 150   px, smaller connected blobs are dropped
%     dilatePx   (1,1) double = 2     final dilation of the mask
%     edgePx     (1,1) double = 4     ignore this border (interp. zero-fill)
%
%   rm fields: map (Y x X correlation), mask, regressor (T x 1), kernel,
%   kernelLagsSec, provisionalMask.
%
%   See also PB_PHOTON_REGISTER, PB_LOWRANK_DENOISE.

arguments
    movie (:,:,:) single
    yaw (:,1) double
    fps (1,1) double
    options.sigma (1,1) double = 1.5
    options.lagsSec (1,2) double = [-1 4]
    options.ridge (1,1) double = 10
    options.refMask logical = false(0)
    options.maskFrac (1,1) double = 0.3
    options.minBlob (1,1) double = 150
    options.dilatePx (1,1) double = 2
    options.edgePx (1,1) double = 4
end

[H, W, T] = size(movie);
assert(numel(yaw) == T, 'pb_rotation_map: %d yaw samples for %d volumes', numel(yaw), T);
ay = abs(yaw);
zs = @(x) (x - mean(x)) / std(x);

Ms = imgaussfilt(movie, options.sigma);
Mz = (Ms - mean(Ms, 3)) ./ max(std(Ms, 0, 3), eps);
clear Ms
corrMap = @(r) mean(Mz .* reshape(single(zs(r)), 1, 1, []), 3);
edgeOk = false(H, W); e = options.edgePx; edgeOk(e+1:H-e, e+1:W-e) = true;

% provisional mask
if isempty(options.refMask)
    r0 = lagged(movmean(ay, max(1, round(2*fps))), round(0.5*fps));
    provMask = toMask(corrMap(r0), edgeOk, options);
else
    provMask = logical(options.refMask);
end

% |yaw| -> photons kernel
L = round(options.lagsSec(1)*fps):round(options.lagsSec(2)*fps);
X = zeros(T, numel(L));
for i = 1:numel(L), X(:,i) = lagged(ay, L(i)); end
X = X - mean(X);
y = squeeze(sum(movie .* provMask, [1 2])); y = double(y - mean(y));
lam = options.ridge * mean(sum(X.^2));
kernel = (X'*X + lam*eye(numel(L))) \ (X'*y);
kc = kernel; kc(L < 0) = 0;
regressor = X * kc;

rm.map = corrMap(regressor);
rm.mask = toMask(rm.map, edgeOk, options);
rm.regressor = regressor;
rm.kernel = kernel;
rm.kernelLagsSec = L(:) / fps;
rm.provisionalMask = provMask;
end

function y = lagged(x, L)
% y(t) = x(t - L), zero padded (L > 0 means y lags x)
y = zeros(size(x));
if L >= 0, y(L+1:end) = x(1:end-L); else, y(1:end+L) = x(1-L:end); end
end

function m = toMask(map, edgeOk, o)
s = imgaussfilt(map, 1.5); s(~edgeOk) = 0;
m = bwareaopen(s > o.maskFrac * max(s(:)), o.minBlob);
if o.dilatePx > 0, m = imdilate(m, strel('disk', o.dilatePx)); end
end
