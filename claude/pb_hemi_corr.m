function out = pb_hemi_corr(x, y, varargin)
% PB_HEMI_CORR  Correlate two PB glomerulus profiles as hemisphere-averaged
% 16-angle compass maps, with an exact circular-rotation null.
%
%   out = pb_hemi_corr(x, y)
%   out = pb_hemi_corr(x, y, 'alignLR', true)      % also estimate / remove a
%                                                  % left-right numbering offset
%   out = pb_hemi_corr(x, y, 'alignLR', true, 'minGain', 0.1)
%
% x, y : 2N x 1 profiles over glomeruli ordered left hemisphere (1..N) then
%        right hemisphere (N+1..2N), the two halves sharing the same angle
%        order (process_im convention: alpha = repmat(linspace(-pi, pi-2pi/N,
%        N), 1, 2); N = 16 in this repo). NaNs allowed (omitted pairwise).
%
% Why not a plain 32-point correlation: the two hemispheres are two copies of
% the same 360-degree map, so 32 glomeruli are ~16 independent angles and a
% 32-point r double-counts; a one-ring circular-shift null would also let the
% left copy slide onto the right, inflating the null. Here each profile is
% averaged across hemispheres into a 16-angle profile, r is computed on those
% 16 points, and the null is the 16 rotations of the compass (s = 0..15,
% enumerated exactly). Because the observed configuration (s = 0) is one of
% the 16, the smallest achievable two-sided p is 1/16 = 0.0625: a single
% profile pair can only ever say "the true alignment beats all 15 rotations";
% significance has to come from replication across trials / holds / flies
% (e.g. (1/16)^k for k independent pairs all best at s = 0).
%
% Left/right alignment ('alignLR'): overshoot_bar_stim_script.m splits the
% skeleton into 16 glomeruli per hemisphere independently, so the right
% half's numbering can be off by a glomerulus relative to the left (seen in
% fly 20261005-1: L/R corr of the mean-z profile 0.49 -> 0.84 after a -1
% shift). The offset is estimated from x ONLY (the reference profile, never
% the response being tested), as the circular shift of the right half that
% maximises corr(left, right), and applied identically to x and y if it
% raises that correlation by at least minGain (default 0.1).
%
% out fields:
%   r, p_perm   16-point Pearson r and exact rotation-null p (>= 1/16)
%   rNull       1 x N: r at each rotation s = 0..N-1 (rNull(1) = r)
%   xH, yH      the hemisphere-averaged N x 1 profiles actually correlated
%   lrOffset    glomeruli the right half was rotated by (0 if not aligned)
%   rLR_x, rLR_y  left-right correlation of x and y, [before after] alignment
%   r32, p32_perm  for reference: the 2N-point r with the same shared-shift
%               (within-hemisphere wrap) null used in earlier analyses
%
% See overshoot_crosstrial_glomerulus_alignment.m for use.

p = inputParser;
p.addParameter('alignLR', false, @(v) islogical(v) || isnumeric(v));
p.addParameter('minGain', 0.1, @isnumeric);
p.parse(varargin{:});
x = x(:); y = y(:); n = numel(x);
assert(numel(y) == n && mod(n, 2) == 0, 'pb_hemi_corr: x and y must be the same even length (left then right hemisphere)');
h = n / 2;
xL = x(1:h); xR = x(h+1:end); yL = y(1:h); yR = y(h+1:end);
c = @(a, b) corr(a, b, 'rows', 'complete');

lrOffset = 0; rLRx0 = c(xL, xR); rLRy0 = c(yL, yR);
if p.Results.alignLR
    rs = arrayfun(@(s) c(xL, circshift(xR, s)), 0:h-1);
    [rBest, iBest] = max(rs);
    if rBest - rLRx0 >= p.Results.minGain
        lrOffset = iBest - 1;
        xR = circshift(xR, lrOffset); yR = circshift(yR, lrOffset);
    end
end
out.lrOffset = lrOffset - h * (lrOffset > h/2);   % report signed, in (-h/2, h/2]: 15 -> -1
out.rLR_x = [rLRx0, c(xL, xR)];
out.rLR_y = [rLRy0, c(yL, yR)];

xH = mean([xL xR], 2, 'omitnan'); yH = mean([yL yR], 2, 'omitnan');
out.xH = xH; out.yH = yH;

rNull = arrayfun(@(s) c(circshift(xH, s), yH), 0:h-1);
out.r = rNull(1);
out.rNull = rNull;
out.p_perm = mean(abs(rNull) >= abs(out.r) - 1e-12);

r32 = arrayfun(@(s) c([circshift(x(1:h), s); circshift(x(h+1:end), s)], y), 0:h-1);
out.r32 = r32(1);
out.p32_perm = mean(abs(r32) >= abs(out.r32) - 1e-12);
end
