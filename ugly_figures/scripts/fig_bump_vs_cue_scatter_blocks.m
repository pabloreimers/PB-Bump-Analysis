%% fig_bump_vs_cue_scatter_blocks
% Scatter of bump phase (PVA of z-scored glomerulus activity) against cue
% (bar) position, chopped into 1-minute blocks, for the sweeping epochs
% BEFORE and AFTER the held-bar + stim-LED period of an overshoot bar-stim
% trial. One row of scatters per sweep epoch (pre / post), one column per
% minute.
%
% Input: <trialDir>\bump_results.mat written by section 6 of
% data scripts\overshoot_bar_stim_script.m (per-volume t, mu, rho, bar angle,
% LED, hold flags, registration-failure flags). Run that script first.
%
% Conventions (same as overshoot_bar_stim_script.m): bar is plotted "signed"
% (multiplied by bar_sign, so the bump ~ +bar diagonal is the identity line);
% volumes with rho < rho_thresh or a failed registration are dropped; the
% circular mean / std of (bump - bar) per block is printed in each title and
% drawn as the red diagonal.
%
% Output: <trialDir>\bumpVsCue_scatter_1minBlocks.png/.fig and a copy in
% ugly_figures\exports\.

%% 0. parameters
trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20260930-2_epg_syt8s_lpsp_cschrimson\trial_002\';
if exist('fbs_trialDir', 'var') && ~isempty(fbs_trialDir), trialDir = fbs_trialDir; end
block_sec = 60;      % block length
min_block_frac = 0.5; % keep a trailing partial block if it is at least this fraction of block_sec
% Time lag between cue and bump. The calcium indicator (and the circuit)
% delay the bump relative to the cue, which smears the scatter along the
% diagonal in a sweep-direction-dependent way. 'optimal' scans lag_range_sec
% and picks the lag that minimises the circular std of (bump - cue) over
% both sweep epochs; 'none' uses the raw traces. Positive lag = bump follows
% the cue (cue is compared at t - lag).
lag_mode = 'optimal';   % 'optimal' | 'none'
lag_range_sec = [-2 5];

thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'exports');
if ~exist(exportDir, 'dir'), mkdir(exportDir); end

%% 1. load and define the two sweep epochs
S = load(fullfile(trialDir, 'bump_results.mat'), 'bump');
b = S.bump;
t = b.t_vol(:)';
ok = b.rho(:)' >= b.rho_thresh & ~b.badVol(:)';

mu_deg = rad2deg(angle(exp(1i * b.mu(:)')));

if isempty(b.holdSegs)
    error('fig_bump_vs_cue_scatter_blocks:noHold', '%s has no held-bar segment -> no pre/post sweep split.', trialDir);
end
% The pre/post split is the LONGEST held-bar segment (the stim hold); short
% protocol holds (e.g. the 5-s hold that ends sweep_12_24_36_hold270_sweep)
% occasionally survive the detector and must not define the split.
[~, iMain] = max(diff(b.holdSegs(:, 1:2), [], 2)); mainHold = b.holdSegs(iMain, 1:2);
% pre-sweep: from the start of the VR protocol (or imaging, whichever is later)
% to the first hold; post-sweep: from the end of the last hold to the end of
% the protocol (or imaging, whichever is earlier). The LED stim sits inside
% the hold, so "post" = after stimulation.
epochs = struct('name', {'pre-stim sweeps', 'post-stim sweeps'}, ...
    'win', {[max(t(1), b.vr_start), mainHold(1)], [mainHold(2), min(t(end), b.vr_end)]});

for e = 1:numel(epochs)
    w = epochs(e).win;
    nFull = floor(diff(w) / block_sec);
    rem_sec = diff(w) - nFull * block_sec;
    edges = w(1) + (0:nFull) * block_sec;
    if rem_sec >= min_block_frac * block_sec
        edges(end+1) = w(2); %#ok<SAGROW>
    end
    epochs(e).edges = edges;
    fprintf('%s: DAQ %.1f-%.1f s (%.1f s) -> %d blocks of %d s%s\n', epochs(e).name, w, diff(w), ...
        numel(edges) - 1, block_sec, ternary(rem_sec >= min_block_frac * block_sec, ...
        sprintf(' (last one %.0f s)', rem_sec), sprintf(' (%.0f s remainder dropped)', rem_sec)));
end

%% 1b. cue-bump lag: scan lags, pick the one that tightens the sweep offset most
% Cue at time t - lag is compared with the bump at t. The cue is interpolated
% on the unit circle (complex linear interp of exp(1i*bar)) so wrap-around at
% +/-180 deg is handled. The lag is chosen once for the whole trial (both
% sweep epochs pooled), then also reported per epoch for reference.
bar_c = exp(1i * b.bar_rad(:)');
inSweep = false(size(t));
for e = 1:numel(epochs), inSweep = inSweep | (t >= epochs(e).win(1) & t < epochs(e).win(2)); end
dt = median(diff(t));
lags = lag_range_sec(1):dt:lag_range_sec(2);
lagScore = nan(numel(epochs) + 1, numel(lags)); % rows: pooled, then each epoch
for L = 1:numel(lags)
    barL = angle(interp1(t, bar_c, t - lags(L), 'linear', 'extrap'));
    offL = angle(exp(1i * (b.mu(:)' - b.bar_sign * barL)));
    lagScore(1, L) = circ_std_local(offL(ok & inSweep));
    for e = 1:numel(epochs)
        selE = ok & t >= epochs(e).win(1) & t < epochs(e).win(2);
        lagScore(1 + e, L) = circ_std_local(offL(selE));
    end
end
[~, iBest] = min(lagScore, [], 2);
switch lag_mode
    case 'optimal'
        lag_used = lags(iBest(1));
        if lag_used < 0 % a bump LEADING the cue is not physical; only happens when tracking is too poor to fit
            fprintf('  note: scan minimum at %.2f s is negative -> using lag 0 instead\n', lag_used);
            lag_used = 0;
        end
    case 'none',    lag_used = 0;
    otherwise, error('unknown lag_mode %s', lag_mode);
end
fprintf('cue-bump lag scan (%.1f to %.1f s): pooled optimum %.2f s (circ std %.1f deg; %.1f deg at lag 0)', ...
    lag_range_sec, lags(iBest(1)), rad2deg(lagScore(1, iBest(1))), rad2deg(lagScore(1, abs(lags) == min(abs(lags)))));
for e = 1:numel(epochs)
    fprintf('; %s optimum %.2f s (%.1f deg)', epochs(e).name, lags(iBest(1 + e)), rad2deg(lagScore(1 + e, iBest(1 + e))));
end
fprintf('\n  -> using lag = %.2f s (lag_mode = %s)\n', lag_used, lag_mode);

bar_lag_rad = angle(interp1(t, bar_c, t - lag_used, 'linear', 'extrap'));
% signed, lagged cue (so bump ~ +cue) and the per-volume offset, (-180, 180]
bar_deg = rad2deg(angle(exp(1i * b.bar_sign * bar_lag_rad)));
off_rad = angle(exp(1i * (b.mu(:)' - b.bar_sign * bar_lag_rad)));

% lag-scan figure
figure(11); clf
set(gcf, 'Position', [100 100 640 420], 'Color', 'w'); hold on
cols = [0 0 0; 0 0.45 0.75; 0.85 0.33 0.1];
lbl = [{'both sweep epochs'}, {epochs.name}];
for r = 1:size(lagScore, 1)
    plot(lags, rad2deg(lagScore(r, :)), '-', 'Color', cols(r, :), 'LineWidth', 1 + (r == 1))
    plot(lags(iBest(r)), rad2deg(lagScore(r, iBest(r))), 'o', 'Color', cols(r, :), 'MarkerFaceColor', cols(r, :), 'HandleVisibility', 'off')
end
xline(lag_used, 'k--', sprintf('used: %.2f s', lag_used), 'HandleVisibility', 'off')
xline(0, 'k:', 'HandleVisibility', 'off')
xlabel('lag (s): cue at t - lag vs bump at t'); ylabel('circ std of bump - cue (deg)')
legend(lbl, 'Location', 'best'); grid on; box on
title(sprintf('%s: cue-bump lag scan', b.trialName), 'Interpreter', 'none')
lagBase = sprintf('bumpVsCue_lagScan');
hide_toolbars(gcf)
exportgraphics(gcf, fullfile(trialDir, [lagBase '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, [strrep(b.trialName, '\', '_') '_' lagBase '.png']), 'Resolution', 150);
fprintf('saved %s\n', fullfile(trialDir, [lagBase '.png']));

%% 2. figure: rows = epochs, columns = 1-min blocks
nCols = max(arrayfun(@(e) numel(e.edges) - 1, epochs));
figure(10); clf
set(gcf, 'Position', [50 50 320 * nCols + 100, 760], 'Color', 'w')
for e = 1:numel(epochs)
    edges = epochs(e).edges;
    for k = 1:numel(edges) - 1
        sel = ok & t >= edges(k) & t < edges(k+1);
        ax = subplot(numel(epochs), nCols, (e-1) * nCols + k);
        hold on
        plot([-180 180], [-180 180], 'k--', 'LineWidth', 0.75)      % bump = bar
        if any(sel)
            moff = rad2deg(circ_mean_local(off_rad(sel)));
            soff = rad2deg(circ_std_local(off_rad(sel)));
            % identity shifted by the block's mean offset (drawn in two
            % pieces so it wraps)
            xx = linspace(-180, 180, 361);
            yy = rad2deg(angle(exp(1i * deg2rad(xx + moff))));
            yy(abs(diff([yy yy(end)])) > 180) = NaN;
            plot(xx, yy, 'r-', 'LineWidth', 1)
            scatter(bar_deg(sel), mu_deg(sel), 6, [0 0.2 0.8], 'filled', ...
                'MarkerFaceAlpha', 0.35, 'MarkerEdgeAlpha', 0)
            ttl = sprintf('%s, min %d (%.0f-%.0f s)\noffset %.0f deg, circ std %.0f deg, n = %d', ...
                epochs(e).name, k, edges(k), edges(k+1), moff, soff, nnz(sel));
        else
            ttl = sprintf('%s, min %d (%.0f-%.0f s)\nno data', epochs(e).name, k, edges(k), edges(k+1));
        end
        axis square; xlim([-180 180]); ylim([-180 180])
        xticks(-180:90:180); yticks(-180:90:180); grid on; box on
        title(ttl, 'FontSize', 8, 'FontWeight', 'normal')
        if k == 1, ylabel('bump phase (deg)'); end
        if e == numel(epochs), xlabel(sprintf('cue position at t - %.2f s (deg)', lag_used)); end
    end
end
sgtitle(sprintf('%s: bump phase vs cue position in %d-s blocks, cue lagged by %.2f s (%s) -- %s; red = identity + block mean offset; rho >= %.1f', ...
    b.trialName, block_sec, lag_used, lag_mode, b.signNote, b.rho_thresh), 'Interpreter', 'none', 'FontSize', 10)

%% 3. export (.png + .fig, in the trial folder and ugly_figures/exports)
if lag_used == 0
    baseName = sprintf('bumpVsCue_scatter_%dsBlocks', block_sec);
else
    baseName = sprintf('bumpVsCue_scatter_%dsBlocks_lag%.2fs', block_sec, lag_used);
end
hide_toolbars(gcf)
for outDir = {trialDir, exportDir}
    if strcmp(outDir{1}, exportDir)
        fn = fullfile(outDir{1}, [strrep(b.trialName, '\', '_') '_' baseName]);
    else
        fn = fullfile(outDir{1}, baseName);
    end
    exportgraphics(gcf, [fn '.png'], 'Resolution', 150);
    vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); % saved .fig should open visibly even from a headless run
    savefig(gcf, [fn '.fig'], 'compact'); set(gcf, 'Visible', vis);
    fprintf('saved %s.png / .fig\n', fn);
end

%% local functions
function hide_toolbars(fig)
% exportgraphics bakes the interactive axes toolbar into the image when run
% headless. The toolbars are created lazily on the first draw, so draw first,
% then hide them (both via findall and via each axes' Toolbar property).
drawnow
set(findall(fig, 'Type', 'axestoolbar'), 'Visible', 'off');
for ax = findall(fig, 'Type', 'axes')'
    try, ax.Toolbar.Visible = 'off'; catch, end
end
drawnow
end

function m = circ_mean_local(a)
m = angle(mean(exp(1i * a(:))));
end

function s = circ_std_local(a)
% circular standard deviation, sqrt(-2 ln R) (same as circ_stats' circ_std s0)
R = abs(mean(exp(1i * a(:))));
s = sqrt(-2 * log(max(R, eps)));
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
