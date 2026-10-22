%% fig_channel_compare
% Is channel 2 a (noisier) copy of channel 1, or a different signal / static
% fluorescence? Both channels' registered z-sum videos (ch1\ and ch2\ folders
% from data scripts\overshoot_split_channels.m + overshoot_bar_stim_script.m)
% are read out with the SAME glomerulus segmentation (ch1's bump.clusterIdx),
% background-subtracted per volume (mean outside the mask), and compared:
%  * mean-image spatial correlation inside the mask
%  * slow trend (bleaching) of each channel
%  * per-glomerulus temporal variability (std / mean)
%  * per-glomerulus correlation between channels vs. its noise ceiling:
%    odd/even-volume split-half reliability of each channel (Spearman-Brown),
%    ceiling = sqrt(rel1 * rel2). If ch2 were ch1 + noise, observed r ~ ceiling;
%    r << ceiling with reliable ch2 = a different dynamic signal; low rel2 with
%    low std/mean = mostly static + noise.
%  * PVA bump position of each channel (same glomeruli, same alpha) vs each other.
% Signals are high-passed (> ~0.05 Hz, 20-s moving mean subtracted) first so
% slow bleaching does not dominate the correlations.

%% 0. parameters
thisDir = fileparts(mfilename('fullpath')); outDir = fullfile(thisDir, '..', 'exports');
trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20261019-1_epg_syt8s_empty_cschrimson\trial_002';
tag = '20261019-1_trial_002';
rho_thresh = 0.1;

%% 1. load
b1 = load(fullfile(trialDir, 'ch1', 'bump_results.mat')).bump; b2 = load(fullfile(trialDir, 'ch2', 'bump_results.mat')).bump;
nVol = numel(b1.t_vol); fs = 1 / median(diff(b1.t_vol));
cI = b1.clusterIdx; nC = max(cI(:)); inMask = cI > 0;
F = cell(1, 2); M = cell(1, 2);
for c = 1:2
    V = load(fullfile(trialDir, sprintf('ch%d', c), 'imgData_sum_reg.mat'), 'imgData_sum_reg'); img = single(V.imgData_sum_reg(:, :, 1:nVol)); clear V
    M{c} = mean(img, 3); X = reshape(img, [], nVol); clear img
    X = X - mean(X(~inMask(:), :), 1);                       % per-volume background
    f = nan(nC, nVol); for g = 1:nC, f(g, :) = mean(X(cI(:) == g, :), 1); end
    F{c} = f; clear X
end
bad = b1.badVol(:)'; for c = 1:2, F{c}(:, bad) = NaN; end

%% 2. statistics
bp = @(x) x - movmean(x, round(20 * fs));      % high-pass (> ~0.05 Hz): removes bleaching; no low-pass so odd/even volume noise stays independent
bp = @(x) bp(fillmissing(x, 'linear', 'EndValues', 'nearest'));
Z = cell(1, 2); R = nan(nC, 2); CV = nan(nC, 2);
for c = 1:2
    Z{c} = nan(nC, nVol);
    for g = 1:nC
        x = F{c}(g, :); CV(g, c) = std(x, 'omitnan') / mean(x, 'omitnan');
        y = bp(x')'; Z{c}(g, :) = (y - mean(y)) / std(y);
        r = corr(y(1:2:end)', y(2:2:end)'); R(g, c) = 2 * r / (1 + r);           % split-half reliability (odd vs even volumes), Spearman-Brown
    end
end
rObs = arrayfun(@(g) corr(Z{1}(g, :)', Z{2}(g, :)'), 1:nC)'; ceil_ = sqrt(max(R(:, 1), 0) .* max(R(:, 2), 0));
rMean = corr(M{1}(inMask), M{2}(inMask));
m1 = mean(F{1}, 1, 'omitnan'); m2 = mean(F{2}, 1, 'omitnan');
fprintf('mean-image correlation inside mask: r = %.2f\n', rMean);
fprintf('mean F (bg-subtracted, a.u.): ch1 %.1f, ch2 %.1f; last/first 10%% of trial: ch1 %.2f, ch2 %.2f\n', mean(m1), mean(m2), mean(m1(end-round(nVol/10):end)) / mean(m1(1:round(nVol/10))), mean(m2(end-round(nVol/10):end)) / mean(m2(1:round(nVol/10))));
fprintf('std/mean over time (median across %d glomeruli): ch1 %.2f, ch2 %.2f\n', nC, median(CV(:, 1), 'omitnan'), median(CV(:, 2), 'omitnan'));
fprintf('split-half reliability (> 0.05 Hz): ch1 %.2f, ch2 %.2f (medians)\n', median(R(:, 1)), median(R(:, 2)));
fprintf('ch1-ch2 r per glomerulus: median %.2f (range %.2f-%.2f); noise ceiling median %.2f; observed/ceiling median %.2f\n', median(rObs), min(rObs), max(rObs), median(ceil_), median(rObs ./ ceil_));
% bump from each channel on the same glomeruli
mu = nan(2, nVol); rho = nan(2, nVol);
for c = 1:2
    [xt, yt] = pol2cart(b1.alpha, Z{c}'); [mu(c, :), rho(c, :)] = cart2pol(mean(xt, 2, 'omitnan'), mean(yt, 2, 'omitnan'));
end
both = rho(1, :) >= rho_thresh & rho(2, :) >= rho_thresh & ~bad; d = angle(exp(1i * (mu(2, :) - mu(1, :))));
Rv = abs(mean(exp(1i * d(both)))); fprintf('bump (same glomeruli): rho ch1 %.2f, ch2 %.2f (means); on %d volumes with both rho>=%.1f: ch2-ch1 circ mean %.0f deg, mean resultant length %.2f (1 = identical)\n', ...
    mean(rho(1, :), 'omitnan'), mean(rho(2, :), 'omitnan'), nnz(both), rho_thresh, rad2deg(angle(mean(exp(1i * d(both))))), Rv);

%% 3. figure
fig = figure('Position', [50 50 1600 900], 'Color', 'w');
subplot(2, 3, 1); imagesc(M{1}); axis image off; colormap(gca, gray); hold on; contour(inMask, [0.5 0.5], 'r'); title('ch1 mean image + ch1 mask')
subplot(2, 3, 2); imagesc(M{2}); axis image off; colormap(gca, gray); hold on; contour(inMask, [0.5 0.5], 'r'); title(sprintf('ch2 mean image (same mask); spatial r = %.2f', rMean))
subplot(2, 3, 3); plot(b1.t_vol, m1 / mean(m1), 'g', b1.t_vol, m2 / mean(m2), 'r'); xlabel('DAQ time (s)'); ylabel('mask-mean F / its mean'); legend('ch1', 'ch2'); title('slow trends (bleaching)')
subplot(2, 3, 4); gsel = find(ceil_ == max(ceil_), 1); tw = b1.t_vol >= 100 & b1.t_vol <= 160;
plot(b1.t_vol(tw), Z{1}(gsel, tw), 'g', b1.t_vol(tw), Z{2}(gsel, tw), 'r'); xlabel('DAQ time (s)'); ylabel('z (> 0.05 Hz)'); title(sprintf('glomerulus %d (highest ceiling): r = %.2f', gsel, rObs(gsel))); legend('ch1', 'ch2')
subplot(2, 3, 5); plot([0 1], [0 1], 'k:'); hold on; plot(ceil_, rObs, 'o', 'MarkerFaceColor', [0.4 0.4 0.4]); xlabel('noise ceiling sqrt(rel1 rel2)'); ylabel('observed ch1-ch2 r'); axis([0 1 -0.2 1]); axis square
title(sprintf('per glomerulus; rel ch1 %.2f, ch2 %.2f', median(R(:, 1)), median(R(:, 2))))
subplot(2, 3, 6); plot(rad2deg(mu(1, both)), rad2deg(mu(2, both)), '.', 'MarkerSize', 4); axis([-180 180 -180 180]); axis square; xlabel('ch1 bump (deg)'); ylabel('ch2 bump (deg)'); title(sprintf('PVA, same glomeruli; %d volumes both rho>=%.1f; R = %.2f', nnz(both), rho_thresh, Rv))
sgtitle([tag ': channel 2 vs channel 1'], 'Interpreter', 'none')
exportgraphics(fig, fullfile(outDir, [tag '_channelCompare.png']), 'Resolution', 130);
