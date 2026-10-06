%% overshoot_ttx_inhibition_map
% Fly 20261005-1 (EPG>syt-GCaMP8s, LPsP>CsChrimson): trial_005 = LED stim
% during a 7-min held bar in normal saline; trial_006 = identical protocol
% with TTX in the bath. With spiking blocked, any LED-locked change in EPG
% fluorescence in trial_006 should be the direct (graded) LPsP -> EPG input,
% free of circuit dynamics and of the fly's visual response to the LED.
%
% Questions: (1) is the LED-evoked change in trial_006 non-uniform across
% the PB, and is that non-uniformity real (vs pulse-to-pulse noise)?
% (2) do the glomeruli with the strongest change in trial_006 coincide with
% where the bump sat during the stimulation in trial_005?
%
% Per glomerulus (32, same skeleton ordering in both trials; alpha = the PVA
% angle assigned to each glomerulus, repeated over the two hemispheres):
%   - trial_006: dF/F = (F_on - F_off) / F_off within the hold, with every
%     LED pulse as a replicate (on = the 5-s pulse, off = the 5 s before it)
%     -> mean, SEM, one-way ANOVA across glomeruli, left-right correlation,
%     dependence on baseline F0 (a purely multiplicative artefact would give
%     constant dF/F).
%   - trial_005: bump occupancy during the LED block = rho-weighted fraction
%     of LED-on volumes whose PVA falls in each glomerulus's angle bin; also
%     the per-glomerulus mean z during LED-on, and the same dF/F as above.
%   - correlation of trial_006 dF/F with each trial_005 quantity, with a
%     circular-shift permutation null (shifting the 16-bin profile within
%     each hemisphere), since neighbouring glomeruli are not independent.
%
% Needs bump_results.mat (with f_cluster, clusterIdx) + mask.mat +
% imgData_sum_reg.mat in each trial folder (overshoot_bar_stim_script.m).
% Output: ttxInhibitionMap.png/.fig in the trial_006 folder and
% ugly_figures\exports\, plus numbers in the command window.

%% 0. inputs
flyDir  = 'Z:\noah_np123\Data\flyg\overshoot\20261005-1_epg_syt8s_lpsp_cschrimson\';
ctlDir  = fullfile(flyDir, 'trial_005\');   % normal saline
ttxDir  = fullfile(flyDir, 'trial_006\');   % TTX
pre_sec = 5;                                % LED-off baseline window before each pulse
nPerm   = 2000;
thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'ugly_figures', 'exports'); if ~exist(exportDir, 'dir'), mkdir(exportDir); end
rng(1)

%% 1. per-glomerulus LED responses in both trials
T = struct();
for which = {'ctl', 'ttx'}
    w = which{1};
    if strcmp(w, 'ctl'), d = ctlDir; else, d = ttxDir; end
    S = load(fullfile(d, 'bump_results.mat')); b = S.bump;
    t = b.t_vol(:)'; f = b.f_cluster; z = b.z_cluster; led = logical(b.led(:)'); inHold = logical(b.inHold(:)');
    nG = size(f, 1);
    % pulses fully inside the hold
    pul = find(arrayfun(@(k) all(inHold(t >= b.led_on(k) & t <= b.led_off(k))), 1:numel(b.led_on)));
    dFF = nan(nG, numel(pul)); dF = dFF; F0 = dFF;
    for k = 1:numel(pul)
        on  = t >= b.led_on(pul(k)) & t < b.led_off(pul(k));
        off = t >= b.led_on(pul(k)) - pre_sec & t < b.led_on(pul(k));
        F0(:, k) = mean(f(:, off), 2, 'omitnan');
        dF(:, k) = mean(f(:, on), 2, 'omitnan') - F0(:, k);
        dFF(:, k) = dF(:, k) ./ F0(:, k);
    end
    % block-wise (all on vs all off within the LED span) for the image-space map
    span = [b.led_on(pul(1)) b.led_off(pul(end))];
    onAll = led & inHold; offAll = ~led & inHold & t >= span(1) & t <= span(2);
    % pulse-triggered average of dF/F0 per glomerulus (-5 .. 10 s)
    dt = median(diff(t)); lags = -pre_sec:dt:10;
    trig = nan(nG, numel(lags), numel(pul));
    for k = 1:numel(pul)
        tk = b.led_on(pul(k)) + lags; okk = tk >= t(1) & tk <= t(end);
        fi = interp1(t, f', tk(okk))';                   % nG x nLags
        trig(:, okk, k) = (fi - F0(:, k)) ./ F0(:, k);
    end
    % bump occupancy during LED-on (and LED-off) in the hold, per alpha bin
    alpha = b.alpha(:)'; aBins = unique(round(alpha, 6)); nB = numel(aBins);   % 16 angles
    mu = b.mu(:)'; rho = b.rho(:)';
    occ = @(sel) accumarray(knnsearch(aBins', angle(exp(1i*mu(sel)))'), rho(sel), [nB 1])' / sum(rho(sel));
    occOn = occ(onAll & rho >= b.rho_thresh); occOff = occ(offAll & rho >= b.rho_thresh);
    T.(w) = struct('name', b.trialName, 'dFF', dFF, 'dF', dF, 'F0', F0, 'pul', pul, 'trig', trig, 'lags', lags, ...
        'zOn', mean(z(:, onAll), 2, 'omitnan'), 'zOff', mean(z(:, offAll), 2, 'omitnan'), ...
        'occOn', occOn, 'occOff', occOff, 'aBins', aBins, 'alpha', alpha, 'nG', nG, ...
        'onAll', onAll, 'offAll', offAll, 'clusterIdx', b.clusterIdx, 'dir', d, 'rhoOn', mean(rho(onAll), 'omitnan'), 'rhoOff', mean(rho(offAll), 'omitnan'));
    fprintf('%s: %d pulses inside the hold; mean dF/F over all glomeruli = %+.1f%% (per-glomerulus range %+.1f .. %+.1f%%); rho on/off %.2f/%.2f\n', ...
        b.trialName, numel(pul), 100*mean(dFF(:), 'omitnan'), 100*min(mean(dFF, 2)), 100*max(mean(dFF, 2)), T.(w).rhoOn, T.(w).rhoOff);
end
nG = T.ttx.nG; gIdx = 1:nG; hemi = [ones(1, nG/2) 2*ones(1, nG/2)];
% each glomerulus's alpha bin index (same in both trials by construction)
binOf = knnsearch(T.ttx.aBins', T.ttx.alpha');

%% 2. is the TTX-trial response non-uniform across the PB?
D = T.ttx.dFF;                                 % nG x nPulses
m = mean(D, 2); se = std(D, 0, 2) / sqrt(size(D, 2));
% one-way ANOVA: glomerulus as factor, pulses as replicates
[pA, tbl] = anova1(D', [], 'off'); Fstat = tbl{2, 5};
% glomeruli individually different from the grand mean (t-test across pulses, Bonferroni)
gm = mean(D(:));
[~, pG] = ttest(D' - gm); sigG = pG < 0.05 / nG;
% left-right symmetry of the profile and dependence on baseline
mL = m(hemi == 1); mR = m(hemi == 2);
[rLR, pLR] = corr(mL, mR);
F0m = mean(T.ttx.F0, 2);
[rF0, pF0] = corr(F0m, m, 'type', 'Spearman');
fprintf('\n== trial_006 (TTX) non-uniformity ==\n');
fprintf('dF/F across glomeruli: mean %+.1f%%, SD %.1f%%, range %+.1f .. %+.1f%%; per-pulse SEM per glomerulus ~%.1f%%\n', 100*mean(m), 100*std(m), 100*min(m), 100*max(m), 100*mean(se));
fprintf('one-way ANOVA (glomerulus, %d pulses as replicates): F = %.1f, p = %.2g\n', size(D, 2), Fstat, pA);
fprintf('%d/%d glomeruli differ from the grand mean at Bonferroni p<0.05: %s\n', nnz(sigG), nG, mat2str(find(sigG)));
fprintf('left vs right hemisphere profile (16 vs 16): r = %.2f (p = %.3f)\n', rLR, pLR);
fprintf('dF/F vs baseline F0 across glomeruli: Spearman r = %.2f (p = %.3f)  [a multiplicative/light artefact predicts ~0]\n', rF0, pF0);
[~, iMost] = sort(m);
fprintf('most inhibited glomeruli: %s (%.1f .. %.1f%%); least: %s (%.1f .. %.1f%%)\n', mat2str(iMost(1:5)'), 100*m(iMost(1)), 100*m(iMost(5)), mat2str(iMost(end-4:end)'), 100*m(iMost(end-4)), 100*m(iMost(end)));

%% 3. does the TTX map line up with the bump in trial_005 during stimulation?
% trial_005 quantities per glomerulus
occ5  = T.ctl.occOn(binOf)';            % bump occupancy of this glomerulus's angle during LED-on (both hemispheres share a bin)
z5on  = T.ctl.zOn;                      % mean z during LED-on
dz5   = T.ctl.zOn - T.ctl.zOff;         % LED-on minus LED-off z in trial_005
dff5  = mean(T.ctl.dFF, 2);             % trial_005 dF/F (same measure as the TTX map)
cands = {'bump occupancy (LED on, trial 005)', occ5; 'mean z during LED on (trial 005)', z5on; ...
         'z(LED on) - z(LED off) (trial 005)', dz5; 'dF/F LED on vs off (trial 005)', dff5};
fprintf('\n== trial_006 dF/F map vs trial_005 bump measures (32 glomeruli) ==\n');
res = struct('name', {}, 'r', {}, 'p_perm', {}, 'rho_s', {});
for c = 1:size(cands, 1)
    x = cands{c, 2}(:); y = m(:);
    r = corr(x, y); rs = corr(x, y, 'type', 'Spearman');
    % null: circularly shift x within each hemisphere by the same random amount (and optionally mirror)
    rNull = nan(nPerm, 1);
    for p = 1:nPerm
        s = randi(nG/2) - 1; xs = x;
        for h = 1:2
            idx = find(hemi == h); xs(idx) = circshift(x(idx), s);
        end
        rNull(p) = corr(xs, y);
    end
    pPerm = mean(abs(rNull) >= abs(r));
    res(c) = struct('name', cands{c, 1}, 'r', r, 'p_perm', pPerm, 'rho_s', rs);
    fprintf('  %-40s Pearson r = %+.2f (circ-shift perm p = %.3f), Spearman r = %+.2f\n', cands{c, 1}, r, pPerm, rs);
end
% also: where was the bump in trial_005 during LED-on vs LED-off, and does the
% TTX map differ between glomeruli the bump visited vs did not
[~, iOcc] = sort(occ5, 'descend'); top = iOcc(1:8); bot = iOcc(end-7:end);
fprintf('TTX dF/F at the 8 glomeruli the trial_005 bump occupied most: %+.1f%% vs least: %+.1f%%\n', 100*mean(m(top)), 100*mean(m(bot)));

%% 4. figure
figure(30); clf; set(gcf, 'Position', [30 30 1800 1000], 'Color', 'w')
cols = lines(2);
% (a) TTX dF/F per glomerulus with SEM, significance
ax1 = subplot(3, 3, [1 2]); hold on
bar(gIdx, 100*m, 'FaceColor', [0.3 0.3 0.8], 'EdgeColor', 'none', 'BarWidth', 0.8)
errorbar(gIdx, 100*m, 100*se, 'k.', 'LineStyle', 'none', 'CapSize', 3)
plot(gIdx(sigG), 100*m(sigG) - sign(m(sigG)').*100*se(sigG) - 1.5*sign(m(sigG)'), 'r*', 'MarkerSize', 5)
yline(100*gm, 'k--', sprintf('mean %+.1f%%', 100*gm), 'LabelHorizontalAlignment', 'left')
xline(nG/2 + 0.5, 'k:'); xlim([0.5 nG + 0.5]); xlabel('glomerulus (1-16 left, 17-32 right)'); ylabel('dF/F, LED on vs 5 s before (%)')
title(sprintf('%s (TTX): per-glomerulus LED response, mean +- SEM over %d pulses; * = differs from mean (Bonf. p<0.05); ANOVA F=%.1f p=%.1g; L/R r=%.2f', ...
    T.ttx.name, size(D, 2), Fstat, pA, rLR), 'Interpreter', 'none', 'FontSize', 9)
% (b) pulse-triggered dF/F heatmap (TTX)
ax2 = subplot(3, 3, 3);
imagesc(T.ttx.lags, gIdx, 100*mean(T.ttx.trig, 3, 'omitnan')); colormap(ax2, ttx_rwb(256)); cl = prctile(abs(100*mean(T.ttx.trig, 3, 'omitnan')), 98, 'all'); clim([-cl cl]); colorbar
hold on; xline(0, 'k'); xline(5, 'k'); yline(nG/2 + 0.5, 'k:')
xlabel('time from LED onset (s)'); ylabel('glomerulus'); title('TTX: pulse-triggered dF/F (%)', 'FontSize', 9)
% (c) image-space dF/F map for TTX
ax3 = subplot(3, 3, 4);
S6 = load(fullfile(ttxDir, 'imgData_sum_reg.mat'), 'imgData_sum_reg'); M6 = load(fullfile(ttxDir, 'mask.mat'));
img6 = S6.imgData_sum_reg; n6 = min(size(img6, 3), numel(T.ttx.onAll)); clear S6
bgv = squeeze(mean(reshape(img6(:, :, 1:n6), [], n6) .* ~M6.mask(:), 1)) * numel(M6.mask) / nnz(~M6.mask);
img6 = img6(:, :, 1:n6) - reshape(single(bgv), 1, 1, []);
Fon = mean(img6(:, :, T.ttx.onAll(1:n6)), 3); Foff = mean(img6(:, :, T.ttx.offAll(1:n6)), 3); clear img6
dffImg = (Fon - Foff) ./ max(Foff, prctile(Foff(M6.mask), 10)); dffImg(~M6.mask) = NaN;
hI = imagesc(100*dffImg); set(hI, 'AlphaData', M6.mask); set(ax3, 'Color', [0.92 0.92 0.92]); axis image; colormap(ax3, ttx_rwb(256));
clim([-25 25]); colorbar; hold on
contour(M6.mask, [0.5 0.5], 'k', 'LineWidth', 0.5)
cI = T.ttx.clusterIdx; for g = 1:nG, [yy, xx] = find(cI == g); if ~isempty(xx), text(mean(xx), mean(yy), num2str(g), 'FontSize', 6, 'HorizontalAlignment', 'center'); end; end
title('TTX: dF/F image (LED on - off) / off, %', 'FontSize', 9)
% (d) trial_005: bump occupancy during LED on/off + z profile
ax4 = subplot(3, 3, 5); hold on
plot(gIdx, occ5 / mean(occ5), '-o', 'Color', cols(1, :), 'MarkerFaceColor', cols(1, :), 'MarkerSize', 4, 'DisplayName', 'bump occupancy, LED on (rel. to uniform)')
plot(gIdx, T.ctl.occOff(binOf) / mean(T.ctl.occOff), '--o', 'Color', [0.5 0.5 0.5], 'MarkerSize', 3, 'DisplayName', 'bump occupancy, LED off')
yyaxis right; plot(gIdx, z5on, '-s', 'Color', cols(2, :), 'MarkerFaceColor', cols(2, :), 'MarkerSize', 4, 'DisplayName', 'mean z, LED on'); ylabel('mean z-score (LED on)')
yyaxis left; ylabel('occupancy / uniform'); xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5]); xlabel('glomerulus')
legend('Location', 'best', 'FontSize', 7); title(sprintf('%s (saline): where the bump was during the LED block', T.ctl.name), 'Interpreter', 'none', 'FontSize', 9)
% (e) overlay of the two profiles
ax5 = subplot(3, 3, 6); hold on
plot(gIdx, zscore(100*m), '-o', 'Color', [0.3 0.3 0.8], 'MarkerFaceColor', [0.3 0.3 0.8], 'MarkerSize', 4, 'DisplayName', 'TTX dF/F (z-scored across glomeruli)')
plot(gIdx, zscore(occ5), '-o', 'Color', cols(1, :), 'MarkerSize', 4, 'DisplayName', 'trial 005 bump occupancy (z-scored)')
plot(gIdx, zscore(z5on), '-s', 'Color', cols(2, :), 'MarkerSize', 4, 'DisplayName', 'trial 005 mean z LED on (z-scored)')
xline(nG/2 + 0.5, 'k:', 'HandleVisibility', 'off'); yline(0, 'k:', 'HandleVisibility', 'off'); xlim([0.5 nG + 0.5]); xlabel('glomerulus'); ylabel('z across glomeruli')
legend('Location', 'best', 'FontSize', 7); title('profiles overlaid', 'FontSize', 9)
% (f-h) scatters
for c = 1:3
    ax = subplot(3, 3, 6 + c); hold on
    x = cands{c, 2}(:);
    scatter(x(hemi == 1), 100*m(hemi == 1), 30, cols(1, :), 'filled', 'DisplayName', 'left')
    scatter(x(hemi == 2), 100*m(hemi == 2), 30, cols(2, :), 'filled', 'DisplayName', 'right')
    text(x + 0.01*range(x), 100*m, string(gIdx'), 'FontSize', 6)
    pf = polyfit(x, 100*m, 1); xx = linspace(min(x), max(x), 2); plot(xx, polyval(pf, xx), 'k-', 'HandleVisibility', 'off')
    xlabel(cands{c, 1}); ylabel('TTX dF/F (%)')
    title(sprintf('r = %+.2f, circ-shift perm p = %.2f (Spearman %+.2f)', res(c).r, res(c).p_perm, res(c).rho_s), 'FontSize', 9)
    if c == 1, legend('Location', 'best', 'FontSize', 7); end
end
sgtitle('LPsP>CsChrimson LED response map under TTX (trial_006) vs bump position during stimulation in saline (trial_005) -- fly 20261005-1', 'FontSize', 11, 'Interpreter', 'none')

drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off')
for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow
base = 'ttxInhibitionMap';
exportgraphics(gcf, fullfile(ttxDir, [base '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, ['20261005-1_trial006_vs_005_' base '.png']), 'Resolution', 150);
vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(ttxDir, [base '.fig']), 'compact'); set(gcf, 'Visible', vis);
fprintf('saved %s.png / .fig\n', fullfile(ttxDir, base));
save(fullfile(ttxDir, 'ttxInhibitionMap_results.mat'), 'T', 'm', 'se', 'sigG', 'pA', 'Fstat', 'rLR', 'rF0', 'res', 'occ5', 'z5on', 'dz5', 'dff5', 'binOf', 'hemi');

%% local functions
function cmap = ttx_rwb(n)
half = floor(n / 2); up = linspace(0, 1, half)';
cmap = [[up up ones(half, 1)]; [ones(n - half, 1) flipud(linspace(0, 1, n - half)') flipud(linspace(0, 1, n - half)')]];
end
