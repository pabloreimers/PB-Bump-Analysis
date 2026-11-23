%% fig_p2x2_empty_example_drift
% Scratch: what is the bump doing in an EMPTY>P2X2 fly with a large mean
% bump displacement per stim (fig_p2x2_displacement_per_stim.m)? Shows every
% trial of one fly across the whole 300 s trial, so the stim-locked movement
% can be judged against the bump's behaviour between stims:
%   row 1: whole-PB atp dF/F (ejections)
%   row 2: EPG im.z heatmap (PB angle x time), bump position overlaid
%          (rho < rho_min masked, then interpolated -- as the displacement
%          script; drawn wrapped, on the heatmap's own PB-angle axis)
%   row 3: unwrapped bump position (masked + interpolated), re-zeroed at the
%          trial start, plus im.rho (right axis)
% Stim onsets = dotted lines; the 0-15 s scoring window after each stim is
% shaded. The command window gives, per stim, the bump's net movement over
% the 15 s BEFORE the stim vs the 15 s after (raw, unflipped sign): if the
% "displacement" is just steady drift, the two match.
%
% Needs the workspace from fig_p2x2_displacement_per_stim.m (all_data,
% genotype, fly_id, trial_is_right, fly_rot); runs it if missing.
% example_fly picks the fly by rank of mean displacement among Empty flies
% (1 = largest).

%% setup
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if ~exist('fly_rot','var')
    run(fullfile(repo_root,'ugly_figures','scripts','fig_p2x2_displacement_per_stim.m')); close all
end
example_rank = 1;
win_s = 15; rho_cut = 0.1;
export_dir = fullfile(repo_root,'ugly_figures','exports','V2','scratch');

e_rows = find(strcmp(genotype,'empty') & ~isnan(trial_is_right));
e_flies = unique(fly_id(e_rows));
[~,ord] = sort(fly_rot{1,1},'descend');
ex_fly = e_flies{ord(example_rank)};
ex_trials = e_rows(strcmp(fly_id(e_rows), ex_fly))';
fprintf('\nexample Empty fly (rank %d): %s, mean displacement %+.2f rad/stim\n', example_rank, ex_fly, 2*pi*fly_rot{1,1}(ord(example_rank)));

%% figure
gcamp_color = [0 .2 1]; atp_color = [.8 0 0];
nt = numel(ex_trials);
fig = figure('Color','w','Position',[50 50 1500 300*nt]);
T = tiledlayout(fig, nt, 1, 'TileSpacing','loose','Padding','compact');
for r = 1:nt
    i = ex_trials(r); tr = all_data(i);
    xb = tr.ft.xb(:); xf = tr.ft.xf(:);
    on = xf(find(diff([false; logical(tr.ft.stims(:))])==1));
    alpha = unwrap(tr.im.alpha);
    mu_m = masked_mu(tr.im.mu(:), tr.im.rho(:), xb, rho_cut);
    side = 'LR'; side = side(trial_is_right(i)+1);

    % per-stim pre vs post (raw sign)
    pre = nan(1,numel(on)); post = pre;
    m = @(a,b) mean(mu_m(xb>=a & xb<b));
    for s = 1:numel(on)
        pre(s)  = m(on(s)-1, on(s))       - m(on(s)-win_s-1, on(s)-win_s);
        post(s) = m(on(s)+win_s-1, on(s)+win_s) - m(on(s)-1, on(s));
    end
    fprintf('  %s (%s-side stim): bump net movement, 15 s BEFORE each stim %s rad | 15 s AFTER %s rad | mean |ball yaw| %.3f\n', ...
        trial_name(tr.meta), side, mat2str(round(pre,1)), mat2str(round(post,1)), mean(abs(tr.ft.r_speed)));

    Tr = tiledlayout(T, 6, 1, 'TileSpacing','none','Padding','none'); Tr.Layout.Tile = r;
    axa = nexttile(Tr, 1); hold(axa,'on')
    plot(axa, xb, mean(tr.atp.d,1), 'Color', atp_color)
    ylabel(axa,'ATP dF/F'); xlim(axa, xb([1 end])); axa.XTick = [];
    title(axa, sprintf('%s, %s-side stims; displacement 0-15 s after stims 1-4 (flipped): %s rad', trial_name(tr.meta), side, ...
        mat2str(round((1-2*trial_is_right(i))*post,1))), 'Interpreter','none')

    axh = nexttile(Tr, 2, [2 1]); hold(axh,'on')
    imagesc(axh, xb, alpha, tr.im.z); colormap(axh, white_to_color(gcamp_color))
    clim(axh, [prctile(min(tr.im.z,[],1),5) prctile(max(tr.im.z,[],1),95)])
    % bump position on both PB halves (alpha spans -pi..pi twice), wrapped; break lines at wraps
    for off = [0 2*pi]
        y = mod(mu_m + pi, 2*pi) - pi + off + (alpha(1) - (-pi));
        y(abs(diff([y; y(end)])) > pi) = nan;
        plot(axh, xb, y, '-', 'Color', [1 .5 0], 'LineWidth', 1)
    end
    ylim(axh, alpha([1 end])); xlim(axh, xb([1 end])); axh.XTick = []; ylabel(axh,'PB angle')

    axu = nexttile(Tr, 4, [3 1]); hold(axu,'on')
    for s = 1:numel(on)
        patch(axu, on(s)+[0 win_s win_s 0], [-1 -1 1 1]*1e3, [.9 .9 .9], 'EdgeColor','none')
    end
    yyaxis(axu,'left'); plot(axu, xb, mu_m - mu_m(1), 'Color', gcamp_color, 'LineWidth', 1.2); ylabel(axu,'bump (rad, unwrapped)')
    yl = [min(mu_m-mu_m(1)) max(mu_m-mu_m(1))] + [-1 1]; ylim(axu, yl)
    yyaxis(axu,'right'); plot(axu, xb, tr.im.rho, 'Color', [.6 .6 .6]); ylabel(axu,'\rho'); ylim(axu,[0 1])
    yline(axu, rho_cut, ':', 'Color', [.4 .4 .4]);
    axu.YAxis(1).Color = gcamp_color; axu.YAxis(2).Color = [.5 .5 .5];
    xlim(axu, xb([1 end])); xlabel(axu, 'time (s)')
    for ax = [axa axh axu]; xline(ax, on, ':k'); end
end
title(T, sprintf('Empty>P2X2 example, %s (rank %d by mean displacement per stim); shading = 0-%d s scoring window', ex_fly, example_rank, win_s), 'Interpreter','none')

if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig, fullfile(export_dir, sprintf('fig_p2x2_empty_example_drift_rank%d.png', example_rank)), 'Resolution', 150)

%% functions
function m = masked_mu(mu, rho, t, rho_min)
    % as fig_p2x2_displacement_per_stim.m: NaN where rho < rho_min, unwrap the
    % rest, interpolate the gaps, nearest value at the ends
    ok = rho >= rho_min & ~isnan(mu);
    m = nan(size(mu)); m(ok) = unwrap(mu(ok));
    m = interp1(t(ok), m(ok), t, 'linear'); m = fillmissing(m, 'nearest');
end
function s = trial_name(meta)
    parts = strsplit(meta,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    if strcmpi(parts{end},'registration'); s = parts{end-1}; else; s = parts{end}; end
end
function cmap = white_to_color(c)
    n = 256; cmap = [linspace(1,c(1),n)', linspace(1,c(2),n)', linspace(1,c(3),n)'];
end
