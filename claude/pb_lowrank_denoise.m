function dn = pb_lowrank_denoise(counts, shifts, mask, options)
%PB_LOWRANK_DENOISE Denoise a photon-count PB movie with a cross-validated low-rank model.
%
%   dn = pb_lowrank_denoise(counts, shifts, mask) takes the UNREGISTERED
%   photon counts (Y x X x T integers, pb_read_raw_tif_planes with
%   'photonThreshold'), the per-volume shifts from pb_photon_register, and a
%   PB mask (pb_rotation_map), and returns a rank-k model of the photon rate
%   inside the mask:
%
%       rate(px,t) = mu(px) + sqrt(mu(px)) * U(px,1:k) * S(1:k,1:k) * V(t,1:k)'
%
%   where V's columns are lightly smoothed over 'w' volumes. Use
%   pb_lowrank_movie(dn) to get the denoised z-summed movie back.
%
%   WHY LOW RANK: after thresholding away the line noise (see
%   fig_dlight_noise_anatomy), what is left is Poisson photon noise -- ~40
%   photons per glomerulus per volume, ~15% noise per volume. A per-pixel
%   temporal filter can only beat that by averaging over seconds. A low-rank
%   model instead pools photons ACROSS the ~3000 PB pixels that share a time
%   course, so each volume keeps its own estimate. The data are whitened by
%   sqrt(mean) first (Poisson variance = mean) so bright pixels do not
%   dominate the SVD.
%
%   HOW k AND w ARE CHOSEN -- split-half cross-validation: the model is fit
%   on half A of the photons and scored by the Poisson log-likelihood of half
%   B's raw pixel counts, for every k in 'ranks' and w in 'smoothWins',
%   relative to a static (mean-only) model. Nothing about B (or yaw) is used
%   to fit, so this measures how well the denoised movie predicts photons it
%   never saw -- PROVIDED A and B are independent.
%
%   USE 'halves' = {oddPlanes, evenPlanes} (pb_read_raw_tif_planes with
%   groups {1:2:9, 2:2:10}). Do NOT rely on the fallback, binomial thinning of
%   the z-summed counts: ~6% of photon pulses spill into the next pixel along
%   the fast axis (event autocorrelation 0.065 at x-lag 1), and thinning pixel
%   by pixel sends the two halves of one photon to A and B. That makes the
%   halves share white noise -- on 20260216-1 trial 1 it showed up as ~40% of
%   the per-cluster "shared signal" sitting above 1.5 Hz, which vanished
%   (7.8e-3 -> -1.8e-4) with the odd/even plane split, since different planes
%   are different laser scans. Odd and even planes see slightly different
%   PB brightness, so B is scored against A's prediction rescaled by the
%   ratio of the two halves' per-pixel means (dynamics, not baseline, are
%   what is being compared; the static reference gets the same rescaling).
%   The halves are registered with the same shifts (hence counts + shifts in,
%   rather than an already registered movie).
%
%   On 20260216-1 trial 1 with the (biased) binomial split, the low-rank
%   model beat every per-pixel moving mean on held-out photons; see
%   ugly_figures/scripts/fig_dlight_pva_snr.m for the plane-split SNR
%   comparison, which is the number to trust. k/w are chosen on half the
%   photons, so the full-data optimum may be slightly higher k.
%
%   STATUS (2026-10-28, 20260216-1 trial 1) -- DOES NOT HELP THE BUMP YET.
%   Scored at the 32-cluster level with odd/even planes (fit odd, predict
%   even), rank 2-3 recovers ~50% of the cluster PATTERN signal vs <0 for
%   the gauss30 smoothing in epg_dlight_script, but the resulting PVA does
%   not follow the cue: corr of 1-s PVA changes with 1-s cue changes during
%   movement is ~0 for k = 1-5 (vs -0.14 old pipeline, -0.17 photon counts +
%   gauss30, null 95% ~0.1). Its higher static R(mu, cue) is an occupancy
%   artefact (low-rank mu sits at the fly's usual heading; shuffled-cue null
%   R is 0.3-0.4). Pixel-level CV with the odd/even split is also unusable
%   (planes 4 um apart differ in fine structure, so every k scores below
%   static). Treat the rank choice below as unvalidated.
%
%   What this can NOT recover: structure that is local in both space and
%   time and never repeats (it would need its own component). On this trial
%   the heading-binned PB profile changed mostly in overall amplitude
%   (moving vs still), not as a moving bump, so that is not the regime here --
%   but check on any fly where you expect a bump.
%
%   Name-value options:
%     sigma      (1,1) double = 1      px, spatial pooling before the SVD
%     ranks      (1,:) double = 1:8    k candidates for CV
%     smoothWins (1,:) double = [1 3 5 9]  w candidates (volumes, moving mean on V)
%     rank       (1,1) double = NaN    fix k (skips CV for k)
%     smoothWin  (1,1) double = NaN    fix w (skips CV for w)
%     halves     (1,2) cell = {}       {countsA, countsB}: independent photon
%                                      halves (odd/even planes) for the CV
%     seed       (1,1) double = 1      rng seed for the fallback thinning
%     verbose    (1,1) logical = true
%
%   dn fields: U (p x k), S (k x 1), V (T x k, smoothed), mu (p x 1), mask,
%   rank, smoothWin, cv (table, if run), meanImage (registered full-frame
%   mean, used outside the mask by pb_lowrank_movie), sigma.
%
%   See also PB_LOWRANK_MOVIE, PB_PHOTON_REGISTER, PB_ROTATION_MAP.

arguments
    counts (:,:,:) {mustBeNumeric}
    shifts (:,2) double
    mask (:,:) logical
    options.sigma (1,1) double = 1
    options.ranks (1,:) double = 1:8
    options.smoothWins (1,:) double = [1 3 5 9]
    options.rank (1,1) double = NaN
    options.smoothWin (1,1) double = NaN
    options.halves cell = {}
    options.seed (1,1) double = 1
    options.verbose (1,1) logical = true
end

T = size(counts, 3);
idx = find(mask(:));
ranks = options.ranks; wins = options.smoothWins;
if ~isnan(options.rank), ranks = options.rank; end
if ~isnan(options.smoothWin), wins = options.smoothWin; end
cv = table();

if numel(ranks) > 1 || numel(wins) > 1
    if isempty(options.halves)
        warning('pb_lowrank_denoise:thinning', ['no ''halves'' given: falling back to binomial thinning, ' ...
            'which shares noise between halves (photon pulses span adjacent pixels) -- pass odd/even planes']);
        rng(options.seed);
        A = single(binornd(double(counts), 0.5)); B = counts - A;
    else
        A = options.halves{1}; B = options.halves{2};
    end
    Br = pb_photon_register(B, 'shifts', shifts, 'verbose', false); clear B
    Ar = pb_photon_register(A, 'shifts', shifts, 'verbose', false); clear A
    Yb = reshape(Br, [], T); Yb = double(Yb(idx, :)); clear Br
    [Ua, Sa, Va, mua] = fitSvd(Ar, idx, options.sigma, max(ranks)); clear Ar
    scale = max(mean(Yb, 2), eps) ./ max(mua, eps);   % B/A brightness per pixel
    ll = @(lam) mean(sum(Yb .* log(max(lam .* scale, 1e-6)) - lam .* scale, 1));
    base = ll(repmat(mua, 1, T));
    [kk, ww] = ndgrid(ranks, wins);
    gain = zeros(numel(kk), 1);
    for i = 1:numel(kk)
        gain(i) = ll(rateFrom(mua, Ua, Sa, Va, kk(i), ww(i))) - base;
    end
    cv = table(kk(:), ww(:), gain, 'VariableNames', {'rank', 'smoothWin', 'heldoutLLgain'});
    [~, best] = max(gain); ranks = kk(best); wins = ww(best);
    if options.verbose
        fprintf('pb_lowrank_denoise: thinning CV, held-out log-lik gain over static (nats/vol):\n');
        G = reshape(gain, size(kk));
        fprintf('   w\\k'); fprintf('%7d', kk(:,1)); fprintf('\n');
        for j = 1:size(G,2), fprintf('%5d ', ww(1,j)); fprintf('%+7.2f', G(:,j)); fprintf('\n'); end
        fprintf('  -> k = %d, w = %d volumes\n', ranks, wins);
    end
end

reg = pb_photon_register(counts, 'shifts', shifts, 'verbose', false);
[U, S, V, mu] = fitSvd(reg, idx, options.sigma, ranks);
dn.U = U(:, 1:ranks);
dn.S = S(1:ranks);
dn.V = movmean(V(:, 1:ranks), wins, 1);
dn.mu = mu;
dn.mask = mask;
dn.rank = ranks;
dn.smoothWin = wins;
dn.cv = cv;
dn.meanImage = mean(imgaussfilt(reg, options.sigma), 3);
dn.sigma = options.sigma;
end

function [U, S, V, mu] = fitSvd(movie, idx, sigma, k)
T = size(movie, 3);
Y = reshape(imgaussfilt(movie, sigma), [], T);
Y = double(Y(idx, :));
mu = mean(Y, 2);
Wt = (Y - mu) ./ sqrt(max(mu, eps));
[U, S, V] = svd(Wt, 'econ');   % p ~ 3000 px: full econ SVD is a few s, more robust than svds
S = diag(S);
U = U(:, 1:k); S = S(1:k); V = V(:, 1:k);
end

function lam = rateFrom(mu, U, S, V, k, w)
lam = mu + sqrt(mu) .* (U(:,1:k) * diag(S(1:k)) * movmean(V(:,1:k), w, 1)');
end
