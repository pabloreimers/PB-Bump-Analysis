%% overshoot_crosstrial_glomerulus_alignment
% Each trial in overshoot_bar_stim_script.m gets its own mask, skeleton and
% 32-glomerulus split, so glomerulus k in one trial need not be the same
% tissue as glomerulus k in another. This script quantifies that for a pair
% of trials of the same fly and redoes the cross-trial comparison with ONE
% segmentation:
%   1. rigid xy shift between the two trials' mean images (brute-force
%      correlation search inside the union of the two masks);
%   2. the reference trial's clusterIdx shifted into the test trial's frame;
%      per test-trial glomerulus: which reference glomerulus it mostly
%      overlaps, the fraction of its pixels in that glomerulus, Dice, index
%      offset, centroid distance;
%   3. the test trial's per-glomerulus LED response (early 0-2 s / late 2-5 s
%      dF/F, baseline -1.5..0 s, pulses as replicates) recomputed with the
%      REFERENCE segmentation, and correlated with the reference trial's bump
%      occupancy (LED on / LED-off gaps / pre / post), with circular-shift
%      permutation p -- the same analysis as overshoot_ttx_occupancy_vs_phase
%      but with indices that line up by construction. Index-matched (own
%      segmentation) numbers are printed alongside for comparison.
% Output: crossTrialAlignment.png/.fig in the test trial folder + exports.

%% 0. inputs
refDir  = 'Z:\noah_np123\Data\flyg\overshoot\20261005-3_epg_syt8_lpsp_cschrimson\trial_003\';  % segmentation + bump occupancy
testDir = 'Z:\noah_np123\Data\flyg\overshoot\20261005-3_epg_syt8_lpsp_cschrimson\trial_004\';  % LED response map
if exist('cta_refDir', 'var') && ~isempty(cta_refDir), refDir = cta_refDir; end
if exist('cta_testDir', 'var') && ~isempty(cta_testDir), testDir = cta_testDir; end
maxShift = 15; early_win = [0 2]; late_win = [2 5]; base_win = [-1.5 0]; nPerm = 5000; rng(1)
thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'ugly_figures', 'exports'); if ~exist(exportDir, 'dir'), mkdir(exportDir); end
addpath(fullfile(thisDir, '..', 'claude'));   % pb_hemi_corr

%% 1. load both trials
R = load(fullfile(refDir, 'bump_results.mat')); bR = R.bump; MR = load(fullfile(refDir, 'mask.mat'));
T = load(fullfile(testDir, 'bump_results.mat')); bT = T.bump; MT = load(fullfile(testDir, 'mask.mat'));
coreR = MR.mask; if isfield(MR, 'maskCore'), coreR = MR.maskCore; end
coreT = MT.mask; if isfield(MT, 'maskCore'), coreT = MT.maskCore; end
V = load(fullfile(refDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg'); meanR = mean(V.imgData_sum_reg, 3); clear V
V = load(fullfile(testDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg'); imgT = V.imgData_sum_reg; clear V
nT = min(size(imgT, 3), numel(bT.t_vol)); imgT = imgT(:, :, 1:nT); meanT = mean(imgT, 3);
cR = bR.clusterIdx; cT = bT.clusterIdx; nG = max(cR(:));
fprintf('reference %s, test %s\n', bR.trialName, bT.trialName);

%% 2. rigid shift between the mean images
roi = imdilate(MR.mask | MT.mask, strel('disk', 5));
a = meanR; a = (a - mean(a(roi))) / std(a(roi)); bimg = meanT; bimg = (bimg - mean(bimg(roi))) / std(bimg(roi));
best = -inf; shift = [0 0];
for dy = -maxShift:maxShift
    for dx = -maxShift:maxShift
        bs = imtranslate(bimg, [dx dy], 'FillValues', NaN);          % move test image by (dx,dy) to land on the reference
        m = roi & ~isnan(bs);
        c = corr(a(m), bs(m));
        if c > best, best = c; shift = [dy dx]; end
    end
end
c0 = corr(a(roi), bimg(roi));
fprintf('test -> reference shift: dy = %+d px (rows), dx = %+d px (cols); image corr %.3f (unshifted %.3f); %.2f um/px -> %.1f um\n', ...
    shift(1), shift(2), best, c0, 0.2531, hypot(shift(1), shift(2)) * 0.2531);
% reference clusters in the test frame = shift by (-dy, -dx)
cR_inT = imtranslate(cR, [-shift(2) -shift(1)], 'nearest', 'FillValues', 0);
maskR_inT = imtranslate(double(MR.mask), [-shift(2) -shift(1)], 'nearest', 'FillValues', 0) > 0;
dice_mask = 2 * nnz(maskR_inT & MT.mask) / (nnz(maskR_inT) + nnz(MT.mask));
fprintf('whole-mask Dice after shift: %.2f (before: %.2f)\n', dice_mask, 2 * nnz(MR.mask & MT.mask) / (nnz(MR.mask) + nnz(MT.mask)));

%% 3. per-glomerulus correspondence (test's own glomeruli vs reference glomeruli in the test frame)
map = zeros(nG, 1); frac = zeros(nG, 1); dice = zeros(nG, 1); dcen = zeros(nG, 1);
cenT = zeros(nG, 2); cenR = zeros(nG, 2);
for g = 1:nG
    pixT = cT == g; [yy, xx] = find(pixT); cenT(g, :) = [mean(xx) mean(yy)];
    [yy, xx] = find(cR_inT == g); if ~isempty(xx), cenR(g, :) = [mean(xx) mean(yy)]; else, cenR(g, :) = NaN; end
    lab = cR_inT(pixT); lab = lab(lab > 0);
    if isempty(lab), map(g) = 0; continue; end
    map(g) = mode(lab); frac(g) = mean(lab == map(g));
    dice(g) = 2 * nnz(pixT & cR_inT == map(g)) / (nnz(pixT) + nnz(cR_inT == map(g)));
end
for g = 1:nG, if map(g) > 0, dcen(g) = hypot(cenT(g,1) - cenR(map(g),1), cenT(g,2) - cenR(map(g),2)); else, dcen(g) = NaN; end, end
offs = map - (1:nG)'; offs(map == 0) = NaN;
fprintf('\nper test glomerulus: best-matching reference glomerulus (index offset), fraction of pixels in it, Dice:\n');
for g = 1:nG
    fprintf('  t%02d -> r%02d (%+d)  frac %.2f  Dice %.2f  centroid dist %.1f px\n', g, map(g), offs(g), frac(g), dice(g), dcen(g));
end
fprintf('summary: same index in %d/%d glomeruli; |offset| = 1 in %d; >= 2 in %d; median Dice %.2f; median centroid distance %.1f px (%.1f um)\n', ...
    nnz(offs == 0), nG, nnz(abs(offs) == 1), nnz(abs(offs) >= 2), median(dice(map > 0)), median(dcen, 'omitnan'), median(dcen, 'omitnan') * 0.2531);

%% 4. test-trial LED response with the REFERENCE segmentation
t = bT.t_vol(1:nT); inHold = logical(bT.inHold(1:nT));
% background-subtract per volume (mean outside the test mask), as the main script does
flat = reshape(imgT, [], nT);
bg = mean(flat(~MT.mask(:), :), 1);
fA = nan(nG, nT);   % aligned: reference clusters (shifted) restricted to real-signal pixels of the test trial
nPix = zeros(nG, 1);
for g = 1:nG
    pix = cR_inT(:) == g & coreT(:);
    nPix(g) = nnz(pix);
    if nPix(g) >= 10, fA(g, :) = mean(flat(pix, :), 1) - bg; end
end
clear flat imgT
fprintf('\naligned segmentation: %d/%d reference glomeruli have >= 10 signal pixels in the test trial (dropped: %s)\n', nnz(nPix >= 10), nG, mat2str(find(nPix < 10)'));
fO = bT.f_cluster(:, 1:nT);                                 % own segmentation (index-matched comparison)
pul = find(arrayfun(@(k) all(inHold(t >= bT.led_on(k) & t <= bT.led_off(k))), 1:numel(bT.led_on)));
resp = @(f) deal(nan(nG, numel(pul)), nan(nG, numel(pul)));
[EA, LA] = resp(fA); [EO, LO] = resp(fO);
for k = 1:numel(pul)
    t0 = bT.led_on(pul(k));
    bsel = t >= t0 + base_win(1) & t < t0 + base_win(2); esel = t >= t0 + early_win(1) & t < t0 + early_win(2); lsel = t >= t0 + late_win(1) & t < t0 + late_win(2);
    for which = 1:2
        if which == 1, f = fA; else, f = fO; end
        F0 = mean(f(:, bsel), 2, 'omitnan');
        e = (mean(f(:, esel), 2, 'omitnan') - F0) ./ F0; l = (mean(f(:, lsel), 2, 'omitnan') - F0) ./ F0;
        if which == 1, EA(:, k) = e; LA(:, k) = l; else, EO(:, k) = e; LO(:, k) = l; end
    end
end
mEA = mean(EA, 2); mLA = mean(LA, 2); mEO = mean(EO, 2); mLO = mean(LO, 2);
[~, pEA] = ttest(EA'); [~, sA] = deal([], mean(EA, 2) ./ (std(EA, 0, 2) / sqrt(numel(pul))));
excA = find(pEA' < 0.05 & mEA > 0); inhA = find(pEA' < 0.05 & mEA < 0);
fprintf('aligned-segmentation LED response: early-excited (p<0.05) %s, early-inhibited %s; early vs own-segmentation profile r = %.2f\n', ...
    mat2str(excA'), mat2str(inhA'), corr(mEA(~isnan(mEA) & ~isnan(mEO)), mEO(~isnan(mEA) & ~isnan(mEO))));

%% 5. reference-trial mean z per (reference) glomerulus by epoch, and correlations with the aligned early / late response
tR = bR.t_vol(:)'; zR = bR.z_cluster; ledR = logical(bR.led(:)'); inHR = logical(bR.inHold(:)');
span = [bR.led_on(1) bR.led_off(end)];
ep.on   = ledR & inHR;
ep.gap  = ~ledR & inHR & tR > span(1) & tR < span(2);
ep.pre  = inHR & tR < span(1);
ep.post = inHR & tR > span(2);
en = fieldnames(ep); zm = struct();
for e = 1:numel(en), zm.(en{e}) = mean(zR(:, ep.(en{e})), 2, 'omitnan'); end
zm.onMinusGap = zm.on - zm.gap;                     % the reference trial's own LED-locked relocation
en2 = [en; {'onMinusGap'}];
% Correlations are done on hemisphere-averaged 16-angle profiles with the
% exact 16-rotation null (claude\pb_hemi_corr.m): the two hemispheres are two
% copies of one compass map, so 32 points double-count and the smallest
% achievable p is 1/16. 'alignLR' first removes any left/right numbering
% offset, estimated from the reference profile only.
hemi = [ones(1, nG/2) 2*ones(1, nG/2)]; %#ok<NASGU>
respVars = {'aligned early (0-2 s)', mEA; 'aligned late (2-5 s)', mLA; 'aligned early - late', mEA - mLA; 'index-matched early', mEO; 'index-matched late', mLO};
Rz = nan(size(respVars, 1), numel(en2)); Pz = Rz; R32 = Rz; P32 = Rz; LRoff = Rz;
H = cell(size(respVars, 1), numel(en2));   % hemisphere-averaged (xH, yH) per pair, for the scatters
fprintf('\nreference-trial mean z vs test-trial LED response, hemisphere-averaged 16-angle profiles -- Pearson r (exact 16-rotation p; floor 1/16 = 0.0625):\n%-24s', '');
fprintf('%-20s', en2{:}); fprintf('\n');
for i = 1:size(respVars, 1)
    fprintf('%-24s', respVars{i, 1});
    for j = 1:numel(en2)
        o = pb_hemi_corr(zm.(en2{j}), respVars{i, 2}, 'alignLR', true);
        Rz(i, j) = o.r; Pz(i, j) = o.p_perm; R32(i, j) = o.r32; P32(i, j) = o.p32_perm; LRoff(i, j) = o.lrOffset; H{i, j} = o;
        fprintf('%+.2f (p=%.3f)     ', o.r, o.p_perm);
    end
    fprintf('\n');
end
% groups: aligned early-excited / early-inhibited (uncorrected p<0.05) vs the rest
[~, pLA] = ttest(LA'); inhLate = find(pLA' < 0.05 & mLA < 0); excLate = find(pLA' < 0.05 & mLA > 0);
fprintf('\nmean z in the reference trial at aligned early-excited glomeruli %s vs the rest (rank-sum p):\n', mat2str(excA'));
for j = 1:numel(en2)
    x = zm.(en2{j}); g1 = ismember((1:nG)', excA) & ~isnan(mEA); g0 = ~ismember((1:nG)', excA) & ~isnan(mEA);
    if nnz(g1) >= 2, p = ranksum(x(g1), x(g0)); else, p = NaN; end
    fprintf('  %-11s %+.3f vs %+.3f (p = %.2f)\n', en2{j}, mean(x(g1)), mean(x(g0)), p);
end
if ~isempty(inhA)
    fprintf('... and at aligned early-inhibited glomeruli %s vs the rest:\n', mat2str(inhA'));
    for j = 1:numel(en2)
        x = zm.(en2{j}); g1 = ismember((1:nG)', inhA); g0 = ~g1 & ~isnan(mEA);
        fprintf('  %-11s %+.3f vs %+.3f (p = %.2f)\n', en2{j}, mean(x(g1)), mean(x(g0)), ranksum(x(g1), x(g0)));
    end
end

%% 6. figure
figure(40); clf; set(gcf, 'Position', [30 30 1900 1050], 'Color', 'w')
cE = [0.85 0.4 0.2]; cL = [0.3 0.3 0.8];
ax1 = subplot(3, 4, 1); imagesc(meanR); axis image; colormap(ax1, gray); hold on
B = bwboundaries(cR > 0); for k = 1:numel(B), plot(B{k}(:,2), B{k}(:,1), 'r', 'LineWidth', 0.5); end
for g = 1:nG, [yy, xx] = find(cR == g); if ~isempty(xx), text(mean(xx), mean(yy), num2str(g), 'Color', 'y', 'FontSize', 6, 'HorizontalAlignment', 'center'); end, end
title(sprintf('%s (reference): mean image + its glomeruli', regexprep(bR.trialName, '.*\\', '')), 'Interpreter', 'none', 'FontSize', 8)
ax2 = subplot(3, 4, 2); imagesc(meanT); axis image; colormap(ax2, gray); hold on
B = bwboundaries(cT > 0); for k = 1:numel(B), plot(B{k}(:,2), B{k}(:,1), 'c', 'LineWidth', 0.5); end
for g = 1:nG, [yy, xx] = find(cT == g); if ~isempty(xx), text(mean(xx), mean(yy), num2str(g), 'Color', 'c', 'FontSize', 6, 'HorizontalAlignment', 'center'); end, end
for g = 1:nG, if all(~isnan(cenR(g, :))), text(cenR(g,1), cenR(g,2) + 6, num2str(g), 'Color', 'r', 'FontSize', 6, 'HorizontalAlignment', 'center'); end, end
contour(maskR_inT, [0.5 0.5], 'r', 'LineWidth', 0.75)
title(sprintf('%s (test): own glomeruli (cyan) + reference shifted dy=%+d, dx=%+d px (red)', regexprep(bT.trialName, '.*\\', ''), shift(1), shift(2)), 'Interpreter', 'none', 'FontSize', 8)
ax3 = subplot(3, 4, 3); hold on
bar(1:nG, offs, 'FaceColor', [0.4 0.4 0.4], 'EdgeColor', 'none'); ylabel('index offset (ref - test)')
yyaxis right; plot(1:nG, dice, 'o-', 'Color', cE, 'MarkerFaceColor', cE, 'MarkerSize', 4); ylabel('Dice'); ylim([0 1])
xline(nG/2 + 0.5, 'k:'); xlim([0.5 nG + 0.5]); xlabel('test-trial glomerulus')
title(sprintf('correspondence: same index %d/%d, |offset|=1: %d, >=2: %d; median Dice %.2f', nnz(offs == 0), nG, nnz(abs(offs) == 1), nnz(abs(offs) >= 2), median(dice(map > 0))), 'FontSize', 8)
ax4 = subplot(3, 4, 4); hold on
bar((1:nG) - 0.2, 100*mEO, 0.4, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none', 'DisplayName', 'early, own segmentation')
bar((1:nG) + 0.2, 100*mEA, 0.4, 'FaceColor', cE, 'EdgeColor', 'none', 'DisplayName', 'early, reference segmentation (aligned)')
yline(0, 'k:', 'HandleVisibility', 'off'); xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5])
xlabel('glomerulus'); ylabel('early dF/F (%)'); legend('Location', 'best', 'FontSize', 7)
title(sprintf('test-trial early response, two segmentations (r = %.2f)', corr(mEA, mEO, 'rows', 'complete')), 'FontSize', 8)
% row 2: aligned early/late response with the reference mean z overlaid, and the r heatmap
ax5 = subplot(3, 4, [5 6]); hold on
yyaxis left
bar((1:nG) - 0.2, 100*mEA, 0.4, 'FaceColor', cE, 'EdgeColor', 'none', 'DisplayName', 'test: early (0-2 s) dF/F, aligned')
bar((1:nG) + 0.2, 100*mLA, 0.4, 'FaceColor', cL, 'EdgeColor', 'none', 'DisplayName', 'test: late (2-5 s) dF/F, aligned')
plot(excA, 100*max(mEA(excA), mLA(excA)) + 0.3, 'k*', 'MarkerSize', 5, 'DisplayName', 'early > 0 (p<0.05)')
if ~isempty(inhA), plot(inhA, 100*min(mEA(inhA), mLA(inhA)) - 0.3, 'kv', 'MarkerSize', 4, 'MarkerFaceColor', 'k', 'DisplayName', 'early < 0 (p<0.05)'); end
ylabel('test-trial LED response, dF/F (%)'); yline(0, 'k:', 'HandleVisibility', 'off')
yyaxis right
plot(1:nG, zm.on, '-o', 'Color', [0.1 0.6 0.1], 'MarkerFaceColor', [0.1 0.6 0.1], 'MarkerSize', 4, 'LineWidth', 1.2, 'DisplayName', 'reference: mean z, LED on')
plot(1:nG, zm.gap, '--s', 'Color', 'k', 'MarkerSize', 4, 'LineWidth', 1, 'DisplayName', 'reference: mean z, LED-off gaps')
yline(0, 'k:', 'HandleVisibility', 'off'); ylabel('reference-trial mean z-score')
xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5]); xlabel('reference glomerulus')
legend('Location', 'best', 'FontSize', 7); title('aligned LED response (test trial) vs mean glomerulus activity during the LED block (reference trial)', 'FontSize', 8)
ax6 = subplot(3, 4, 7);
rwb = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
imagesc(Rz); colormap(ax6, rwb); clim([-1 1]); colorbar
xticks(1:numel(en2)); xticklabels(strrep(en2, 'onMinusGap', 'on - gap')); yticks(1:size(respVars, 1)); yticklabels(respVars(:, 1)); set(ax6, 'FontSize', 7)
for i = 1:size(respVars, 1), for j = 1:numel(en2), text(j, i, sprintf('%+.2f\np=%.2f', Rz(i, j), Pz(i, j)), 'HorizontalAlignment', 'center', 'FontSize', 7); end, end
title('r(test LED response, reference mean z by epoch); hemisphere-averaged 16 angles, exact rotation p (floor 0.0625)', 'FontSize', 8)
ax7 = subplot(3, 4, 8); hold on
scatter(100*mLA, 100*mEA, 36, 1:nG, 'filled'); colormap(ax7, hsv(nG)); text(100*mLA + 0.1, 100*mEA, string(1:nG), 'FontSize', 6)
lim = [min([100*mLA; 100*mEA]) max([100*mLA; 100*mEA])] + [-1 1]; plot(lim, lim, 'k--'); xline(0, 'k:'); yline(0, 'k:')
xlabel('aligned late dF/F (%)'); ylabel('aligned early dF/F (%)'); title(sprintf('test trial: early vs late across glomeruli, r = %.2f', corr(mLA, mEA, 'rows', 'complete')), 'FontSize', 8)
% row 3: the four scatters the question is about
pairs = {'early', mEA, 'on'; 'early', mEA, 'gap'; 'late', mLA, 'on'; 'late', mLA, 'gap'};
for k = 1:4
    ax = subplot(3, 4, 8 + k); hold on
    i = strcmp(pairs{k, 1}, 'early') * 1 + strcmp(pairs{k, 1}, 'late') * 2; j = find(strcmp(en2, pairs{k, 3}));
    o = H{i, j}; x = o.xH; y = 100*o.yH; ok = ~isnan(x) & ~isnan(y);        % 16 hemisphere-averaged angles
    if strcmp(pairs{k, 1}, 'early'), col = cE; else, col = cL; end
    scatter(x(ok), y(ok), 48, col, 'filled')
    text(x(ok) + 0.01*range(x(ok)), y(ok), string(find(ok)), 'FontSize', 7)
    pf = polyfit(x(ok), y(ok), 1); xx = linspace(min(x(ok)), max(x(ok)), 2); plot(xx, polyval(pf, xx), 'k-')
    xline(0, 'k:'); yline(0, 'k:')
    xlabel(sprintf('reference mean z, LED %s (L/R averaged)', strrep(pairs{k, 3}, 'gap', 'off (gaps)'))); ylabel(sprintf('test aligned %s dF/F (%%), L/R averaged', pairs{k, 1}))
    title(sprintf('%s vs z(%s): r = %+.2f, exact rotation p = %.3f (16 angles; 32-pt r = %+.2f)', pairs{k, 1}, pairs{k, 3}, Rz(i, j), Pz(i, j), R32(i, j)), 'FontSize', 8)
end
sgtitle(sprintf('%s (reference segmentation + bump activity)  vs  %s (LED response): mean z during LED on / off-gaps vs early / late dF/F -- hemisphere-averaged, L/R offset %d glomerulus', ...
    bR.trialName, bT.trialName, mode(LRoff(:))), 'Interpreter', 'none', 'FontSize', 11)
drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off')
for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow
base = 'crossTrialAlignment_meanZ';
exportgraphics(gcf, fullfile(testDir, [base '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, [strrep(bT.trialName, '\', '_') '_vs_' regexprep(bR.trialName, '.*\\', '') '_' base '.png']), 'Resolution', 150);
vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(testDir, [base '.fig']), 'compact'); set(gcf, 'Visible', vis);
fprintf('saved %s.png / .fig\n', fullfile(testDir, base));
save(fullfile(testDir, 'crossTrialAlignment_results.mat'), 'shift', 'map', 'frac', 'dice', 'offs', 'dcen', 'mEA', 'mLA', 'mEO', 'mLO', 'zm', 'Rz', 'Pz', 'excA', 'inhA', 'cR_inT');
