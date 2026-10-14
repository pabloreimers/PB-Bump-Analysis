%% fig_bump_stim_aligned
% Stimulus-aligned bump position for every trial of the overshoot flies since
% firstDate (PB+EB panels where bump_results_PB/EB.mat exist, otherwise one
% panel from the older single-structure bump_results.mat). For each LED onset a [-5 +5] s window of the PVA bump
% angle (bump_results_PB/EB.mat, 'mu') is cut and the bump angle at t = 0 is
% subtracted (circular difference, wrapped to +/-180 deg). Samples with
% rho < rho_thresh are NaN; pulses whose t = 0 sample is NaN are dropped.
% Thin = single pulses (green = bar held from 1 s before onset to 5 s after, grey = not),
% thick = circular mean across pulses. Bump angle is in the structure's own
% frame (not sign-flipped to the display), so the sign of a shift is
% PB/EB-order dependent.
% One png per trial in ugly_figures/exports.

%% 0. parameters
thisDir = fileparts(mfilename('fullpath'));
outDir = fullfile(thisDir, '..', 'exports');
root = 'Z:\noah_np123\Data\flyg\overshoot';
firstDate = '20260930';   % all flies with folder date >= this (-99 test folders skipped)
win = [-5 5];

%% 1. per trial
flies = dir(root); flies = flies([flies.isdir] & cellfun(@(n) ~isempty(regexp(n, '^\d{8}-\d', 'once')) && str2double(n(1:8)) >= str2double(firstDate) && isempty(regexp(n, '^\d{8}-99_', 'once')), {flies.name}));
for fi = 1:numel(flies)
    trials = dir(fullfile(flies(fi).folder, flies(fi).name, 'trial_*'));
    for ti = 1:numel(trials)
        tdir = fullfile(trials(ti).folder, trials(ti).name);
        if isfile(fullfile(tdir, 'bump_results_PB.mat')) && isfile(fullfile(tdir, 'bump_results_EB.mat'))
            files = {fullfile(tdir, 'bump_results_PB.mat'), fullfile(tdir, 'bump_results_EB.mat')}; names = {'PB', 'EB'};
        elseif isfile(fullfile(tdir, 'bump_results.mat'))
            files = {fullfile(tdir, 'bump_results.mat')}; names = {'PB'};
        else, continue; end
        nS = numel(files);
        tag = [flies(fi).name '_' trials(ti).name];
        fig = figure('Visible', 'off', 'Position', [50 50 700*nS 520], 'Color', 'w');
        for k = 1:nS
            B = load(files{k}); b = B.bump;
            if isempty(b.led_on), continue; end
            dt = median(diff(b.t_vol)); tq = win(1):dt:win(2);
            z = exp(1i * b.mu(:)'); z(b.rho < b.rho_thresh | b.badVol) = NaN;
            nP = numel(b.led_on); tr = nan(nP, numel(tq)); inHold = false(nP, 1);
            for p = 1:nP
                zz = interp1(b.t_vol, z, b.led_on(p) + tq, 'linear', NaN);
                z0 = interp1(b.t_vol, z, b.led_on(p), 'linear', NaN);
                if isnan(z0), continue; end
                tr(p, :) = rad2deg(angle(zz .* conj(z0)));
                inHold(p) = any(b.holdSegs(:,1) <= b.led_on(p) - 1 & b.holdSegs(:,2) >= b.led_on(p) + 5);
            end
            ok = ~isnan(tr(:, find(tq >= 0, 1)));
            subplot(1, nS, k); hold on
            patch([0 median(b.led_off - b.led_on) median(b.led_off - b.led_on) 0], [-180 -180 180 180], [1 0.55 0.55], 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            for p = find(ok)'
                % break lines at wraps so +/-180 jumps aren't drawn across the panel
                y = tr(p, :); y(abs(diff([y y(end)])) > 180) = NaN;
                plot(tq, y, '-', 'Color', [0.55 0.55 0.55 0.5] + inHold(p) * [-0.55 0.2 -0.35 0], 'LineWidth', 0.75);
            end
            m = rad2deg(angle(mean(exp(1i * deg2rad(tr(ok, :))), 1, 'omitnan'))); m(abs(diff([m m(end)])) > 180) = NaN;
            plot(tq, m, 'k', 'LineWidth', 2.5); xline(0, 'k:'); yline(0, 'k:')
            ylim([-180 180]); yticks(-180:90:180); xlim(win); xlabel('time from LED onset (s)'); ylabel('bump position - bump(0) (deg)')
            title(sprintf('%s: %d/%d pulses with bump at t=0 (%d in holds)', names{k}, nnz(ok), nP, nnz(ok & inHold)), 'FontSize', 10)
        end
        sgtitle(tag, 'Interpreter', 'none')
        exportgraphics(fig, fullfile(outDir, [tag '_bumpStimAligned.png']), 'Resolution', 130); close(fig);
        fprintf('saved %s\n', tag);
    end
end
