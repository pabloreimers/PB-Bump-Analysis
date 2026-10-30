function hm = pb_bump_hmm(y, alpha, options)
%PB_BUMP_HMM Bump phase from photon counts with a hidden-Markov (phase random walk) model.
%
%   hm = pb_bump_hmm(y, alpha) takes photon counts per PB cluster per volume
%   (y: nCluster x T, e.g. 32 glomerulus clusters from process_im's
%   centroids) and each cluster's angle alpha (1 x nCluster, process_im's
%   im.alpha: 16 angles repeated over the two hemispheres), and returns the
%   posterior over bump phase at every volume.
%
%   WHY: at dLight photon rates a cluster gets ~40 photons per volume and the
%   bump modulates that by only a few %, so a per-volume population vector is
%   mostly noise, and smoothing it in time (epg_dlight_script's 30-frame
%   Gaussian) smears turns. A generic low-rank denoiser (pb_lowrank_denoise)
%   kept the slow non-bump structure and lost the bump entirely. This model
%   instead says exactly what a bump is and lets every photon vote:
%
%     y_c(t) ~ Poisson( b_c * g(t) * (1 + a * (f_kappa(alpha_c - phi(t)) - m_kappa)) )
%     phi(t) = phi(t-1) + N(0, sigmaPhi^2)  (wrapped), plus a small jump prob.
%
%   f_kappa is a von Mises bump (peak 1), m_kappa its mean over the cluster
%   angles (so the bump redistributes photons rather than adding them),
%   b_c each cluster's baseline, g(t) the PB-wide gain (Gaussian-smoothed
%   total counts, 'gainWin' volumes -- that is the non-bump dopamine signal
%   and is estimated separately). phi lives on an nPhi-point grid and is
%   inferred by forward-backward, so each volume's estimate uses the
%   photons of that volume plus the continuity of the phase, not a fixed
%   smoothing window: the walk can move fast when the photons demand it.
%   The cue/heading is never used.
%
%   TUNING (defaults below), on 20260305-10, the best-tracking closed-loop
%   trial in .data/epg_dlight_20260417.mat: run the model separately on odd
%   and even planes (independent photons, same bump) and maximise the
%   agreement of the two phase estimates, |mean(exp(i(mu_odd - mu_even)))|.
%   This never looks at the cue. Held-out photon likelihood (hm.predict on
%   the other half) is nearly flat in sigmaPhi and slightly prefers jumps
%   (pJump 1e-2: 2.20 vs pJump 0: 2.14 nats/vol), but jumps let the posterior
%   mean hop ~pi/2-pi between modes while the fly is still, which wrecks
%   tracking -- so pJump = 0 and the agreement criterion are used.
%     odd-even agreement: 0.92 (vs 0.77 for epg_dlight_script's gauss30 PVA)
%     corr(1-s change in mu, 1-s change in cue) while moving: -0.48
%       (vs -0.345 for the stored process_im PVA; null 5% -0.09; on HALF the
%       photons the HMM still gives -0.49 / -0.44)
%     same at 0.3 s / 0.5 s steps: -0.27 / -0.37 (vs -0.18 / -0.24)
%     R(mu, -cue) 0.60 (vs 0.58), median error after offset 37 deg (vs 39)
%   The cue-vs-bump sign is negative throughout this dataset (90% of
%   closed-loop trials), i.e. the bump phase moves as -cue in im.alpha's
%   coordinates.
%
%   Name-value options:
%     kappa    (1,1) double = 2      bump concentration
%     amp      (1,1) double = 0.1    fractional bump amplitude a
%     sigmaPhi (1,1) double = 0.06   rad per volume, phase random-walk SD
%     pJump    (1,1) double = 0      prob. per volume of a uniform phase jump
%     gainWin  (1,1) double = 10     volumes, Gaussian window for g(t)
%     nPhi     (1,1) double = 64     phase grid points
%
%   hm fields: mu (T x 1 posterior circular mean phase), rho (posterior
%   resultant length, 0..1 = confidence), post (nPhi x T), phiGrid,
%   b, g, logLik (data log-likelihood), rate (nCluster x T posterior mean
%   rate), predict(@(bOther) rate for another half with baselines bOther).
%
%   See also PB_LOWRANK_DENOISE.

arguments
    y (:,:) double
    alpha (1,:) double
    options.kappa (1,1) double = 2
    options.amp (1,1) double = 0.1
    options.sigmaPhi (1,1) double = 0.06
    options.pJump (1,1) double = 0
    options.gainWin (1,1) double = 10
    options.nPhi (1,1) double = 64
end

[nC, T] = size(y);
assert(numel(alpha) == nC, 'pb_bump_hmm: %d angles for %d clusters', numel(alpha), nC);
tot = sum(y, 1);
g = smoothdata(tot, 2, 'gaussian', options.gainWin); g = g / mean(g);
b = mean(y, 2) ./ mean(g);                       % nC x 1
phi = (0:options.nPhi-1) * 2*pi / options.nPhi;  % 1 x nPhi

% bump shape per (cluster, phase), mean-zero over clusters of one hemisphere set
f = exp(options.kappa * (cos(alpha(:) - phi) - 1));      % nC x nPhi
f = f - mean(f, 1);
shape = max(1 + options.amp * f, 1e-6);                  % nC x nPhi

% emission log-lik (up to constants): sum_c y log(lam) - lam
%   lam = b_c g_t shape_c,phi -> y'*log(shape) + const(t) - g_t * (b'*shape)
logE = log(shape)' * y - (b' * shape)' * g;              % nPhi x T
logE = logE - max(logE, [], 1);
E = exp(logE);

% transition: wrapped Gaussian kernel + uniform jumps (circulant)
dphi = angle(exp(1i * phi));
k = exp(-0.5 * (dphi / options.sigmaPhi).^2); k = k / sum(k);
k = (1 - options.pJump) * k + options.pJump / options.nPhi;
K = fft(k(:));
step = @(p) real(ifft(fft(p) .* K));                     % circular convolution

fw = zeros(options.nPhi, T); c = zeros(1, T);
p = ones(options.nPhi, 1) / options.nPhi;
for t = 1:T
    p = (t > 1) * step(p) + (t == 1) * p;
    p = max(p, 0) .* E(:, t); c(t) = sum(p); p = p / c(t); fw(:, t) = p;
end
bw = ones(options.nPhi, 1); post = zeros(options.nPhi, T); post(:, T) = fw(:, T);
for t = T-1:-1:1
    bw = step(max(bw, 0) .* E(:, t+1)); bw = bw / sum(bw);
    q = fw(:, t) .* bw; post(:, t) = q / sum(q);
end

z = exp(1i * phi) * post;
hm.mu = angle(z(:));
hm.rho = abs(z(:));
hm.post = post;
hm.phiGrid = phi;
hm.b = b; hm.g = g;
hm.logLik = sum(log(c)) + sum(max(log(shape)' * y - (b' * shape)' * g, [], 1));
hm.rate = (b .* (shape * post)) .* g;                    % posterior mean rate
hm.predict = @(bOther) (bOther .* (shape * post)) .* g;
hm.options = options;
end
