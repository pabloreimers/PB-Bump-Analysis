%% fig_dlight_bump_hmm
% Bump phase from photon counts with pb_bump_hmm vs the stored process_im PVA
% (data scripts\epg_dlight_script.m), on one epg>dLight trial. Phase SNR is
% measured without the cue: the model is run separately on odd and even
% planes (independent photons) and the two phase estimates are compared.
% ~2 min (one streaming read of the raw tif).

%% 0. trial + data
repo = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(fullfile(repo, 'claude')); addpath(repo);
dataFile = fullfile(repo, '.data', 'epg_dlight_20260417.mat');
trialIdx = 95;            % 20260305-10: best cue tracking of the closed-loop trials
photonThr = 50;
fps = 10.51; lagVol = 10; % 1-s steps for the change-vs-change tracking metric

L = load(dataFile, 'all_data'); a = L.all_data(trialIdx); clear L
expDir = fileparts(a.meta);
md = pb_scan_metadata(expDir);
mask = load(fullfile(expDir, 'mask.mat')).mask;

%% 1. photons: odd and even planes in one pass, register on their sum
vp = md.valid_planes;
o = pb_read_raw_tif_planes(md.tifPath, md.nplanes, {vp(1:2:end), vp(2:2:end)}, 'photonThreshold', photonThr);
cO = o.stacks{1}; cE = o.stacks{2}; clear o
[~, shifts] = pb_photon_register(cO + cE);

%% 2. cluster counts with the same 32 centroids/angles as process_im
[ym, xm] = find(mask);
[~, cid] = pdist2(a.im.centroids, [ym xm], 'euclidean', 'smallest', 1);
M = sparse(cid, sub2ind(size(mask), ym, xm), 1, 32, numel(mask));
clus = @(c) full(M * double(reshape(pb_photon_register(c, 'shifts', shifts, 'verbose', false), [], size(c,3))));
yO = clus(cO); yE = clus(cE); yAll = yO + yE; clear cO cE
alpha = a.im.alpha;

%% 3. behaviour on the imaging clock
xb = a.ft.xb(:); nV = size(yAll, 2);
cue = angle(interp1(a.ft.xf, exp(1i*a.ft.cue(:)), xb(1:nV)));
yawAbs = abs(interp1(a.ft.xf, a.ft.r_speed(:), xb(1:nV)));
moving = movmax(yawAbs, 21) > 0.3;

%% 4. bump HMM (defaults = tuned on this trial by odd/even agreement)
hA = pb_bump_hmm(yAll, alpha); hO = pb_bump_hmm(yO, alpha); hE = pb_bump_hmm(yE, alpha);
muOld = a.im.mu(1:nV);
pvaG30 = @(y) pva_gauss(y, alpha, 30);   % process_im's PVA on cluster counts
agree = @(m1, m2) abs(mean(exp(1i*(m1 - m2))));
fprintf('odd/even phase agreement: HMM %.3f, gauss30 PVA %.3f\n', agree(hO.mu, hE.mu), agree(pvaG30(yO), pvaG30(yE)));
steps = round([0.3 0.5 1 2 4] * fps);
dcorr = @(mu, L) stepcorr(mu, cue, moving, L);
fprintf('corr(d mu, d cue) while moving, by step:\n');
for L = steps, fprintf('  %.1f s: HMM %.3f, stored PVA %.3f\n', L/fps, dcorr(hA.mu, L), dcorr(muOld, L)); end

%% 5. figure
t = (1:nV)'/fps;
[~, pk] = max(movsum(moving, round(60*fps))); win = max(1, pk-round(30*fps)):min(nV, pk+round(30*fps));
off = angle(mean(exp(1i*(hA.mu(win) + cue(win))), 'omitnan'));
wrapNaN = @(x) setjumps(mod(x, 2*pi));
fig = figure('Color', 'w', 'Position', [40 40 1500 950]);
tl = tiledlayout(fig, 4, 3, 'TileSpacing', 'compact');
[~, nm] = fileparts(expDir);
title(tl, sprintf('%s: bump phase from photon counts (pb_bump_hmm) vs stored PVA', nm), 'Interpreter', 'none')
nexttile([1 3]); imagesc(t(win), hA.phiGrid, hA.post(:, win)); axis xy; colormap(gca, flipud(gray)); hold on
plot(t(win), wrapNaN(-cue(win) + off), 'r', 'LineWidth', 1.2); ylabel('phase (rad)')
title('HMM posterior over bump phase (all photons); red = -cue + offset')
nexttile([1 3]); plot(t(win), wrapNaN(muOld(win)), 'Color', [0.2 0.6 0.2]); hold on
plot(t(win), wrapNaN(hA.mu(win)), 'b'); plot(t(win), wrapNaN(-cue(win) + off), 'r'); ylim([0 2*pi]); xlim(t(win([1 end])))
legend('stored PVA (gauss30)', 'HMM', '-cue + offset', 'Box', 'off'); ylabel('rad')
nexttile([1 3]); ok = ~isnan(cue);
plot(t, unwrap(hA.mu) - hA.mu(1), 'b'); hold on; plot(t, unwrap(muOld) - muOld(1), 'Color', [0.2 0.6 0.2])
cu = -unwrap(cue(ok)); plot(t(ok), cu - cu(1), 'r'); xlim(t([1 end])); ylabel('unwrapped rad'); xlabel('s')
legend('HMM', 'stored PVA', '-cue', 'Box', 'off', 'Location', 'northwest')
nexttile; b = [arrayfun(@(L) dcorr(hA.mu, L), steps); arrayfun(@(L) dcorr(muOld, L), steps)]';
bar(steps/fps, -b); legend('HMM', 'stored PVA', 'Box', 'off', 'Location', 'northwest'); xlabel('step (s)'); ylabel('-corr(\Delta\mu, \Delta cue)'); title('tracking while moving')
nexttile; bar([agree(hO.mu, hE.mu) agree(pvaG30(yO), pvaG30(yE))]); xticklabels({'HMM', 'gauss30 PVA'}); ylim([0 1])
ylabel('|mean e^{i(\mu_{odd}-\mu_{even})}|'); title('phase agreement, independent halves')
nexttile; d1 = angle(exp(1i*(hO.mu - hE.mu))); d2 = angle(exp(1i*(pvaG30(yO) - pvaG30(yE))));
histogram(d2*180/pi, -180:10:180, 'Normalization', 'probability', 'EdgeColor', 'none'); hold on
histogram(d1*180/pi, -180:10:180, 'Normalization', 'probability', 'EdgeColor', 'none'); xlabel('\mu_{odd} - \mu_{even} (deg)'); legend('gauss30 PVA', 'HMM', 'Box', 'off'); title('half-vs-half phase difference')
exportgraphics(fig, fullfile(repo, 'ugly_figures', 'exports', sprintf('dlight_bump_hmm_%s.png', nm)), 'Resolution', 110);

%% functions
function mu = pva_gauss(y, alpha, w)
% process_im's PVA steps on cluster counts: gaussian smoothing, dF/F (7th pct F0), z-score, population vector
y = smoothdata(y, 2, 'gaussian', w); f0 = prctile(y, 7, 2); z = zscore((y - f0) ./ f0, [], 2);
mu = angle(exp(1i*alpha(:)).' * z).';
end

function r = stepcorr(mu, cue, moving, L)
t = find(moving); t = t(t > L);
d = @(x) angle(exp(1i*(x(t) - x(t-L))));
dm = d(mu(:)); dc = d(cue(:)); ok = ~isnan(dc);
r = corr(dm(ok), dc(ok));
end

function x = setjumps(x)
x(abs(diff([x; x(end)])) > pi) = NaN;
end
