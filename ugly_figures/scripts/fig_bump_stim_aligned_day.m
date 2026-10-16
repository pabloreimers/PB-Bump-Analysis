%% fig_bump_stim_aligned_day
% One figure, one subpanel per imaged trial of a single day: PB bump position
% aligned to each LED onset (window `win`, minus the circular mean bump angle over
% the 5 s before onset in 'baseline' mode). Light = single pulses, thick black = circular mean across
% pulses. Only pulses that fall during a bar hold are used (a held-bar segment
% from claude/overshoot scripts' holdSegs covers [onset-1, onset+5] s); samples
% with rho < rho_thresh or in bad volumes are NaN and pulses whose onset sample
% is NaN are dropped. Trials without a bump_results.mat (no tif) are not shown.
% Bump angle is in the PB's own order (not sign-flipped to the display).
% Three versions (modes): 'baseline' = bump - circular mean of the bump over the 5 s before onset (pulses with < 3 valid volumes there are dropped); 'raw' = bump angle as
% is (does the bump land in the same place every time?); 'offset' = bump -
% bar_sign*bar(t - lag) (bump.offset from the bump_results; the bump-bar offset
% through the hold, not baseline subtracted).
% Section 2 adds three stimulus-aligned glomerulus heatmaps (bump.z_cluster, mean
% over the in-hold pulses of each trial): raw, rotated so the mean bump of the 5 s
% pre-onset sits at the centre glomerulus, and rotated so the bar position at hold
% onset (bar_sign * bar, PB order) sits at the centre. The rotation is done per
% hemisphere (glomeruli 1-16 left, 17-32 right; each half spans the full circle).

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, '..', 'exports');
root = 'Z:\noah_np123\Data\flyg\overshoot';
dateStr = '20261013';
if exist('figDate', 'var') && ~isempty(figDate), dateStr = figDate; end   % batch override
win = [-10 20];

%% 1. collect trials and plot
flies = dir(fullfile(root, [dateStr '-*'])); flies = flies([flies.isdir] & ~contains({flies.name}, [dateStr '-99_']));
tdirs = {}; tags = {};
for fi = 1:numel(flies)
    tr_ = dir(fullfile(flies(fi).folder, flies(fi).name, 'trial_*'));
    for ti = 1:numel(tr_)
        f = fullfile(tr_(ti).folder, tr_(ti).name, 'bump_results.mat');
        if isfile(f), tdirs{end+1} = f; tags{end+1} = [extractBefore(flies(fi).name, '_') ' ' tr_(ti).name]; end %#ok<SAGROW>
    end
end
nT = numel(tdirs); nCol = min(4, nT); nRow = ceil(nT / nCol);
modes = {'baseline', 'raw', 'offset'};
ylabs = {'bump - mean(pre 5 s) (deg)', 'bump position (deg)', 'bump - bar (deg)'};
ttls = {'PB bump position minus its mean over the 5 s before stim onset', 'PB bump position (raw, PB order)', 'PB bump minus bar position (offset; sign + lag fitted per trial)'};
for mi = 1:numel(modes)
fig = figure('Position', [50 50 420*nCol 330*nRow], 'Color', 'w');
for i = 1:nT
    b = load(tdirs{i}).bump;
    dt = median(diff(b.t_vol)); tq = win(1):dt:win(2);
    if strcmp(modes{mi}, 'offset'), ang = b.offset(:)'; else, ang = b.mu(:)'; end
    z = exp(1i * ang); z(b.rho < b.rho_thresh | b.badVol) = NaN;
    nP = numel(b.led_on); tr = nan(nP, numel(tq));
    for p = 1:nP
        if ~any(b.holdSegs(:,1) <= b.led_on(p) - 1 & b.holdSegs(:,2) >= b.led_on(p) + 5), continue; end
        if strcmp(modes{mi}, 'baseline')
            zp = z(b.t_vol >= b.led_on(p) - 5 & b.t_vol < b.led_on(p));
            if nnz(~isnan(zp)) < 3, continue; end
            z0 = exp(1i * angle(mean(zp, 'omitnan')));
        else, z0 = 1; end
        tr(p, :) = rad2deg(angle(interp1(b.t_vol, z, b.led_on(p) + tq, 'linear', NaN) .* conj(z0)));
    end
    ok = any(~isnan(tr), 2);
    subplot(nRow, nCol, i); hold on
    patch([0 median(b.led_off - b.led_on) median(b.led_off - b.led_on) 0], [-180 -180 180 180], [1 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    for p = find(ok)'
        y = tr(p, :); y(abs(diff([y y(end)])) > 180) = NaN;
        plot(tq, y, '-', 'Color', [0.2 0.45 0.9 0.35], 'LineWidth', 0.75);   % wraps broken so +/-180 jumps aren't drawn across
    end
    if any(ok)
        m = rad2deg(angle(mean(exp(1i * deg2rad(tr(ok, :))), 1, 'omitnan'))); m(abs(diff([m m(end)])) > 180) = NaN;
        plot(tq, m, 'k', 'LineWidth', 2.5);
    end
    xline(0, 'k:'); yline(0, 'k:'); ylim([-180 180]); yticks(-180:90:180); xlim(win)
    title(sprintf('%s: %d/%d pulses in holds', tags{i}, nnz(ok), nP), 'FontSize', 9, 'Interpreter', 'none')
    if mod(i - 1, nCol) == 0, ylabel(ylabs{mi}); end
    if i > nT - nCol, xlabel('time from LED onset (s)'); end
end
sgtitle(sprintf('%s: %s, aligned to stim onset (pulses during bar holds only; red = LED on)', dateStr, ttls{mi}))
fn = {'baseline_allTrials', 'raw_allTrials', 'barOffset_allTrials'};
exportgraphics(fig, fullfile(outDir, sprintf('%s_bumpStimAligned_%s.png', dateStr, fn{mi})), 'Resolution', 130);
end

%% 2. stimulus-aligned glomerulus heatmaps (mean over in-hold pulses per trial)
hmModes = {'raw', 'shiftBump', 'shiftBar'};
hmTtl = {'raw glomerulus z-score', 'rotated to mean bump of the 5 s pre-onset (centre = bump)', 'rotated to bar position at hold onset (centre = bar, PB order)'};
cm = [[linspace(0,1,128)', linspace(0,1,128)', ones(128,1)]; [ones(128,1), linspace(1,0,128)', linspace(1,0,128)']];
for mi = 1:numel(hmModes)
fig = figure('Position', [50 50 420*nCol 330*nRow], 'Color', 'w');
for i = 1:nT
    b = load(tdirs{i}).bump;
    dt = median(diff(b.t_vol)); tq = win(1):dt:win(2);
    nC = size(b.z_cluster, 1); nBlk = nC / 2; cen = nBlk / 2 + 1;       % centre glomerulus of each hemisphere (alpha = 0)
    z = exp(1i * b.mu(:)'); z(b.rho < b.rho_thresh | b.badVol) = NaN;
    nP = numel(b.led_on); H = nan(nC, numel(tq), nP);
    for p = 1:nP
        h = find(b.holdSegs(:,1) <= b.led_on(p) - 1 & b.holdSegs(:,2) >= b.led_on(p) + 5, 1);
        if isempty(h), continue; end
        switch hmModes{mi}
            case 'raw', ctr = 0;
            case 'shiftBump'
                zp = z(b.t_vol >= b.led_on(p) - 5 & b.t_vol < b.led_on(p));
                if nnz(~isnan(zp)) < 3, continue; end
                ctr = angle(mean(zp, 'omitnan'));
            case 'shiftBar'
                ctr = angle(exp(1i * b.bar_sign * interp1(b.t_vol, b.bar_rad(:)', b.holdSegs(h,1), 'nearest', 'extrap')));
        end
        Hp = interp1(b.t_vol, b.z_cluster', b.led_on(p) + tq, 'linear', NaN)';     % nC x time
        if ~strcmp(hmModes{mi}, 'raw')
            sh = round(cen - (ctr / (2*pi/nBlk) + cen));                          % glomeruli to move the target angle to the centre
            Hp = [circshift(Hp(1:nBlk, :), sh, 1); circshift(Hp(nBlk+1:end, :), sh, 1)];
        end
        H(:, :, p) = Hp;
    end
    ok = squeeze(any(any(~isnan(H), 1), 2));
    subplot(nRow, nCol, i); hold on
    imagesc(tq, 1:nC, mean(H(:, :, ok), 3, 'omitnan')); axis tight; set(gca, 'YDir', 'reverse'); colormap(gca, cm); clim([-2 2])
    xline(0, 'k'); xline(median(b.led_off - b.led_on), 'k--'); yline(nBlk + 0.5, 'k:')
    if ~strcmp(hmModes{mi}, 'raw'), yline(cen, 'k:'); yline(nBlk + cen, 'k:'); end
    title(sprintf('%s: %d/%d pulses in holds', tags{i}, nnz(ok), nP), 'FontSize', 9, 'Interpreter', 'none')
    if mod(i - 1, nCol) == 0, ylabel('glomerulus (1-16 L, 17-32 R)'); end
    if i > nT - nCol, xlabel('time from LED onset (s)'); end
end
cb = colorbar('Position', [0.93 0.3 0.01 0.4]); cb.Label.String = 'z-score';
sgtitle(sprintf('%s: %s, aligned to stim onset (mean over in-hold pulses; solid = LED on, dashed = LED off)', dateStr, hmTtl{mi}))
exportgraphics(fig, fullfile(outDir, sprintf('%s_heatmapStimAligned_%s.png', dateStr, hmModes{mi})), 'Resolution', 130);
end
