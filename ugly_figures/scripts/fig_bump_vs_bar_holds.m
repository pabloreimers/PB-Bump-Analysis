%% fig_bump_vs_bar_holds
% Bump position vs bar position during the short (5-s) held-bar periods of
% the dangling_sweep_12_24_36_prehold270 protocol (or any protocol whose csv
% has holds). Hold windows come from the protocol csv, aligned to the
% measured bar trace by a constant offset, so the full hold is used rather
% than the hold detector's eroded window.
%
% Panels: (1) per-volume bump angle vs held-bar angle, with the circular
% mean +- circ std of each hold; identity and identity+mean offset lines.
% (2) bump angle vs time since hold onset, one line per hold, coloured by
% the held-bar angle (so a bump that stops / drifts / rotates is visible).
% (3) bump - bar offset vs time since hold onset, same colouring. (4) rho.
%
% Input: <trialDir>\bump_results.mat (overshoot_bar_stim_script.m section 6)
% and the protocol csv in the trial folder. Output: bumpVsBar_holds.png/.fig
% in the trial folder and a copy in ugly_figures\exports\.

%% 0. parameters
trialDir = 'Z:\noah_np123\Data\flyg\overshoot\20261005-1_epg_syt8s_lpsp_cschrimson\trial_003\';
if exist('fbh_trialDir', 'var') && ~isempty(fbh_trialDir), trialDir = fbh_trialDir; end
min_hold_sec = 2;         % csv segments with constant bar at least this long count as holds
offset_scan_sec = 20:0.02:80;
thisDir = fileparts(mfilename('fullpath'));
exportDir = fullfile(thisDir, '..', 'exports'); if ~exist(exportDir, 'dir'), mkdir(exportDir); end

%% 1. load, align protocol to the measured bar, find holds
S = load(fullfile(trialDir, 'bump_results.mat')); b = S.bump;
t = b.t_vol(:)'; mu = b.mu(:)'; rho = b.rho(:)'; bar = b.bar_rad(:)';
ok = rho >= b.rho_thresh & ~b.badVol(:)';
csvList = dir(fullfile(trialDir, '*.csv')); csvList = csvList(~contains({csvList.name}, '__'));
assert(numel(csvList) == 1, 'expected one protocol csv in %s', trialDir);
P = readtable(fullfile(trialDir, csvList(1).name));
protoBar = @(tv) deg2rad(interp1(P.time, P.dot_deg_ego, tv, 'linear', NaN));
best = [-inf 0];
for off = offset_scan_sec
    pb = protoBar(t - off); m = ~isnan(pb);
    r = abs(mean(exp(1i*(pb(m) - bar(m)))));
    if r > best(1), best = [r off]; end
end
t_off = best(2);
fprintf('protocol -> DAQ offset %.2f s (R = %.3f)\n', t_off, best(1));
holds = [];  % [t0 t1 bar_deg]
for k = 1:height(P)-1
    if abs(P.dot_deg_ego(k+1) - P.dot_deg_ego(k)) < 1e-6 && P.time(k+1) - P.time(k) >= min_hold_sec
        holds(end+1,:) = [P.time(k) + t_off, P.time(k+1) + t_off, mod(P.dot_deg_ego(k), 360)]; %#ok<SAGROW>
    end
end
holds = holds(holds(:,1) >= t(1) & holds(:,2) <= t(end), :);
fprintf('%d holds inside the imaging window: bar at %s deg, durations %.1f-%.1f s\n', size(holds,1), mat2str(holds(:,3)'), min(diff(holds(:,1:2),[],2)), max(diff(holds(:,1:2),[],2)));

sgn = b.bar_sign;
bar_deg = rad2deg(angle(exp(1i*sgn*bar)));         % signed bar, (-180,180]
mu_deg  = rad2deg(angle(exp(1i*mu)));
cmap = hsv(360);

%% 2. figure
figure(12); clf; set(gcf, 'Position', [40 40 1700 900], 'Color', 'w')
ax1 = subplot(2,2,[1 3]); hold on
plot([-180 180], [-180 180], 'k--', 'LineWidth', 0.75)
allOff = [];
for k = 1:size(holds,1)
    sel = ok & t >= holds(k,1) & t <= holds(k,2);
    hb = rad2deg(angle(exp(1i*deg2rad(sgn*holds(k,3)))));
    c = cmap(mod(round(holds(k,3)),360)+1, :);
    jit = (rand(1, nnz(sel)) - 0.5) * 6;                     % spread the column a little
    scatter(hb + jit, mu_deg(sel), 8, c, 'filled', 'MarkerFaceAlpha', 0.35)
    off_k = angle(exp(1i*(mu(sel) - sgn*deg2rad(holds(k,3)))));
    allOff = [allOff off_k]; %#ok<AGROW>
    m = angle(mean(exp(1i*mu(sel)))); R = abs(mean(exp(1i*mu(sel)))); s = sqrt(-2*log(max(R,eps)));
    errorbar(hb, rad2deg(m), rad2deg(min(s, pi)), 'o', 'Color', c*0.7, 'MarkerFaceColor', c*0.7, 'MarkerSize', 7, 'LineWidth', 1.5, 'CapSize', 8)
    text(hb + 8, rad2deg(m), sprintf('h%d R=%.2f', k, R), 'FontSize', 7, 'Color', c*0.6)
end
mOff = angle(mean(exp(1i*allOff))); sOff = sqrt(-2*log(max(abs(mean(exp(1i*allOff))),eps)));
xx = linspace(-180,180,361); yy = rad2deg(angle(exp(1i*deg2rad(xx + rad2deg(mOff))))); yy(abs(diff([yy yy(end)])) > 180) = NaN;
plot(xx, yy, 'r-', 'LineWidth', 1)
axis square; xlim([-180 180]); ylim([-180 180]); xticks(-180:60:180); yticks(-180:90:180); grid on; box on
xlabel(sprintf('held bar position (deg, signed x%+d)', sgn)); ylabel('bump position (deg)')
title(sprintf('bump vs bar during %d holds (dots = volumes, markers = circ mean +- circ std per hold)\nall-hold mean offset %.0f deg, circ std %.0f deg; red = identity + mean offset', ...
    size(holds,1), rad2deg(mOff), rad2deg(sOff)), 'FontSize', 9)

ax2 = subplot(2,2,2); hold on
for k = 1:size(holds,1)
    sel = t >= holds(k,1) & t <= holds(k,2);
    c = cmap(mod(round(holds(k,3)),360)+1, :);
    y = mu_deg(sel); y(~ok(sel)) = NaN; y(abs(diff([y y(end)])) > 180) = NaN;
    plot(t(sel) - holds(k,1), y, '-', 'Color', c, 'LineWidth', 1.2)
    plot(t(sel) - holds(k,1), rad2deg(angle(exp(1i*deg2rad(sgn*holds(k,3))))) * ones(1, nnz(sel)), ':', 'Color', c, 'LineWidth', 0.8)
end
ylim([-180 180]); yticks(-180:90:180); xlabel('time since hold onset (s)'); ylabel('bump position (deg)'); grid on; box on
title('bump trajectory within each hold (solid; colour = held-bar angle, dotted = that bar angle)', 'FontSize', 9)

ax3 = subplot(2,2,4); hold on
for k = 1:size(holds,1)
    sel = t >= holds(k,1) & t <= holds(k,2);
    c = cmap(mod(round(holds(k,3)),360)+1, :);
    y = rad2deg(angle(exp(1i*(mu(sel) - sgn*deg2rad(holds(k,3)))))); y(~ok(sel)) = NaN; y(abs(diff([y y(end)])) > 180) = NaN;
    plot(t(sel) - holds(k,1), y, '-', 'Color', c, 'LineWidth', 1.2)
end
yline(0, 'k:'); ylim([-180 180]); yticks(-180:90:180); xlabel('time since hold onset (s)'); ylabel('bump - bar (deg)'); grid on; box on
title('offset within each hold', 'FontSize', 9)
colormap(ax1, cmap); cb = colorbar(ax1, 'Location', 'eastoutside'); clim(ax1, [0 360]); cb.Label.String = 'held bar angle (deg, display)'; cb.Ticks = 0:60:360;
sgtitle(sprintf('%s: bump vs bar position during the held-bar periods (%s; rho >= %.1f)', b.trialName, b.signNote, b.rho_thresh), 'Interpreter', 'none', 'FontSize', 10)

%% 3. export
drawnow; set(findall(gcf, 'Type', 'axestoolbar'), 'Visible', 'off');
for ax = findall(gcf, 'Type', 'axes')', try, ax.Toolbar.Visible = 'off'; catch, end, end
drawnow
base = 'bumpVsBar_holds';
exportgraphics(gcf, fullfile(trialDir, [base '.png']), 'Resolution', 150);
exportgraphics(gcf, fullfile(exportDir, [strrep(b.trialName, '\', '_') '_' base '.png']), 'Resolution', 150);
vis = get(gcf, 'Visible'); set(gcf, 'Visible', 'on'); savefig(gcf, fullfile(trialDir, [base '.fig']), 'compact'); set(gcf, 'Visible', vis);
fprintf('saved %s.png / .fig\n', fullfile(trialDir, base));
