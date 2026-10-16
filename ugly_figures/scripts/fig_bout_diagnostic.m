%% fig_bout_diagnostic
% Why does Figure 3's bout-slope mobility metric read ~2 for a control fly
% whose bump looks like it tracks? One figure per dark pre-ATP trial of
% the chosen fly:
%   left column (time):  im.z heatmap with im.mu overlaid; the SAME
%                        smoothed, unwrapped bump and heading the metric
%                        uses (walk_trial_bouts: mu upsampled to 60 Hz then
%                        1 s Gaussian; r_speed 1 s Gaussian, integrated),
%                        with every detected walking bout shaded and
%                        numbered; |smoothed r_speed| with the 0.25 rad/s
%                        bout threshold; rho.
%   right column:        the bout scatter (bump path vs heading path per
%                        bout, numbered, sized by sqrt duration) with the
%                        pooled through-origin slope; per-bout ratio vs
%                        heading path (log x) to show which bouts carry
%                        the ratio; and how the slope changes if bouts
%                        below a heading-path floor are dropped.
% Per-bout numbers are printed so individual bouts can be looked up.

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports','bout_diagnostic');
if ~isfolder(export_dir), mkdir(export_dir), end
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% which fly / which trials
fly_tag = '20260728/fly 18';
smooth_window = 60; turn_thresh = .25; max_gap_frames = 30; min_walking_frames = 30; % Figure 3 settings

%% labels
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
idx = find(contains(fly_id, fly_tag)); [~,o] = sort(trial_num(idx)); idx = idx(o);
p = find(has_pulse(idx),1); pre = idx(1:p-1); dark_pre = pre(is_dark(pre));
fprintf('%s: dark pre-ATP trials: %s\n', fly_tag, mat2str(trial_num(dark_pre)'));

%% one figure per dark pre trial
for i = dark_pre'
    tr = all_data(i); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf));
    mu_f = interp1(xb, unwrap(tr.im.mu(:)), xf, 'linear', 'extrap');
    mu_s = smoothdata(mu_f, 'gaussian', smooth_window);
    r_s  = smoothdata(tr.ft.r_speed(:), 'gaussian', smooth_window);
    head_s = cumtrapz(xf, r_s);
    rho_f = interp1(xb, tr.im.rho(:), xf, 'linear', 'extrap');
    is_walking = abs(r_s) > turn_thresh; is_walking = imclose(is_walking, ones(max_gap_frames+1,1)); is_walking = bwareaopen(is_walking, min_walking_frames);
    d = diff([false; is_walking; false]); b0 = find(d==1); b1 = find(d==-1)-1; nb = numel(b0);
    bump_p = nan(nb,1); head_p = bump_p; dur = bump_p; rho_b = bump_p; net_head = bump_p; net_bump = bump_p;
    for b = 1:nb
        seg = b0(b):b1(b);
        bump_p(b) = sum(abs(diff(mu_s(seg)))); head_p(b) = sum(abs(r_s(seg)))*dt; dur(b) = numel(seg)*dt;
        rho_b(b) = mean(rho_f(seg)); net_head(b) = head_s(seg(end))-head_s(seg(1)); net_bump(b) = mu_s(seg(end))-mu_s(seg(1));
    end
    w = sqrt(dur); slope = (w.*head_p) \ (w.*bump_p);
    ratio_b = bump_p ./ head_p;

    % slope with a heading-path floor
    floors = [0 .5 1 2 3 5];
    slope_floor = nan(size(floors)); n_floor = slope_floor;
    for k = 1:numel(floors)
        ok = head_p >= floors(k); n_floor(k) = sum(ok);
        if sum(ok) >= 5, slope_floor(k) = (w(ok).*head_p(ok)) \ (w(ok).*bump_p(ok)); end
    end
    % unweighted and ratio-of-sums alternatives
    slope_unw = head_p \ bump_p; ratio_sums = sum(bump_p)/sum(head_p);
    % bump path outside bouts, for context
    bump_out = sum(abs(diff(mu_s(~is_walking)))); head_out = sum(abs(r_s(~is_walking)))*dt;

    fprintf('\n== trial %d (%.0f s): %d bouts, pooled slope %.2f (unweighted %.2f, ratio of sums %.2f); in-bout fraction %.0f%%; outside bouts: bump %.0f rad vs heading %.0f rad\n', ...
        trial_num(i), xf(end)-xf(1), nb, slope, slope_unw, ratio_sums, 100*mean(is_walking), bump_out, head_out);
    fprintf('   slope with heading-path floor %s rad: %s  (n bouts %s)\n', mat2str(floors), mat2str(round(slope_floor,2)), mat2str(n_floor));
    fprintf('   %-4s %-7s %-6s %-8s %-8s %-6s %-6s %-8s %-8s\n', 'bout','t0(s)','dur','headPath','bumpPath','ratio','rho','netHead','netBump');
    for b = 1:nb, fprintf('   %-4d %-7.1f %-6.1f %-8.2f %-8.2f %-6.2f %-6.2f %-+8.2f %-+8.2f\n', b, xf(b0(b)), dur(b), head_p(b), bump_p(b), ratio_b(b), rho_b(b), net_head(b), net_bump(b)); end

    % ---- figure
    fig = figure(70); clf; set(fig,'Color','w','Units','pixels','Position',[20 20 1900 950])
    tl = tiledlayout(4, 3, 'TileSpacing','compact','Padding','compact');
    % heatmap
    ax1 = nexttile(tl, 1, [1 2]); hold(ax1,'on')
    imagesc(ax1, xb, 1:32, tr.im.z, [-1.5 3]); colormap(ax1, parula)
    alpha = tr.im.alpha(:);
    for h = 1:2, rows = (h-1)*16+(1:16); y = interp1(alpha(rows), rows, tr.im.mu, 'linear','extrap'); y(abs(diff([y(1);y]))>8) = nan; plot(ax1, xb, y, '-', 'Color',[1 1 1 .8], 'LineWidth', .6); end
    set(ax1,'YDir','normal','YTick',[1 16 17 32]); ylabel(ax1,'PB wedge'); xlim(ax1,[xf(1) xf(end)]); ylim(ax1,[.5 32.5])
    title(ax1, sprintf('%s, dark pre-ATP trial %d: im.z with im.mu (white)', fly_tag, trial_num(i)), 'Interpreter','none')
    % traces with bouts
    ax2 = nexttile(tl, 4, [1 2]); hold(ax2,'on')
    yl = [min([mu_s-mu_s(1); head_s]) max([mu_s-mu_s(1); head_s])]; yl = yl + [-1 1]*0.05*diff(yl);
    for b = 1:nb
        col = [.88 .88 .88]; if mod(b,2), col = [.8 .8 .8]; end
        patch(ax2, xf([b0(b) b1(b) b1(b) b0(b)]), yl([1 1 2 2]), col, 'EdgeColor','none', 'HandleVisibility','off')
        text(ax2, mean(xf([b0(b) b1(b)])), yl(2)-0.04*diff(yl), sprintf('%d', b), 'FontSize', 6, 'HorizontalAlignment','center', 'Color', [.3 .3 .3], 'HandleVisibility','off')
    end
    plot(ax2, xf, head_s, 'k-', 'LineWidth', 1, 'DisplayName', 'heading (\int smoothed r\_speed)')
    plot(ax2, xf, mu_s - mu_s(1), '-', 'Color', [.85 .33 .1], 'LineWidth', 1, 'DisplayName', 'bump (im.mu, 1 s Gaussian)')
    plot(ax2, xf, mu_f - mu_f(1), '-', 'Color', [.85 .33 .1 .3], 'LineWidth', .5, 'DisplayName', 'bump, unsmoothed')
    xlim(ax2,[xf(1) xf(end)]); ylim(ax2, yl); ylabel(ax2,'unwrapped angle (rad)'); legend(ax2,'Location','best','FontSize',7,'Box','off')
    title(ax2, sprintf('the traces the metric uses; shaded + numbered = walking bouts (%d)', nb))
    % r_speed
    ax3 = nexttile(tl, 7, [1 2]); hold(ax3,'on')
    plot(ax3, xf, abs(r_s), 'k-', 'LineWidth', .8); yline(ax3, turn_thresh, 'r--');
    area(ax3, xf, double(is_walking)*turn_thresh, 'FaceColor',[.85 .85 .85], 'EdgeColor','none', 'FaceAlpha', .6)
    xlim(ax3,[xf(1) xf(end)]); ylabel(ax3,'|r\_speed| smoothed (rad/s)'); ylim(ax3,[0 max(2, prctile(abs(r_s),99.5))])
    title(ax3, sprintf('bout definition: |smoothed r\\_speed| > %.2f rad/s, gaps <= %.1f s closed, bouts < %.1f s dropped', turn_thresh, max_gap_frames*dt, min_walking_frames*dt))
    % rho
    ax4 = nexttile(tl, 10, [1 2]); hold(ax4,'on')
    plot(ax4, xb, tr.im.rho, '-', 'Color', [.4 .2 .6], 'LineWidth', .8); yline(ax4, .3, ':r')
    xlim(ax4,[xf(1) xf(end)]); ylabel(ax4,'rho'); xlabel(ax4,'time (s)'); ylim(ax4,[0 max(1.2, max(tr.im.rho))])
    linkaxes([ax1 ax2 ax3 ax4],'x')

    % scatter with bout numbers
    ax5 = nexttile(tl, 3); hold(ax5,'on')
    scatter(ax5, head_p, bump_p, 15 + 8*sqrt(dur), rho_b, 'filled', 'MarkerFaceAlpha', .7); clim(ax5, [.2 .8]); colormap(ax5, parula)
    text(ax5, head_p, bump_p, arrayfun(@(b) sprintf(' %d', b), 1:nb, 'UniformOutput', false), 'FontSize', 6)
    xm = max([head_p; bump_p; 1]); plot(ax5, [0 xm], slope*[0 xm], 'r-', 'LineWidth', 1.5); plot(ax5, [0 xm], [0 xm], 'k:')
    axis(ax5,'square'); xlim(ax5,[0 xm*1.05]); ylim(ax5,[0 xm*1.05]); xlabel(ax5,'heading path per bout (rad)'); ylabel(ax5,'bump path per bout (rad)')
    cb = colorbar(ax5); cb.Label.String = 'mean rho in bout';
    title(ax5, sprintf('bout scatter: pooled slope %.2f (size = sqrt dur)', slope))
    % ratio vs heading path
    ax6 = nexttile(tl, 6); hold(ax6,'on')
    scatter(ax6, head_p, ratio_b, 15 + 8*sqrt(dur), rho_b, 'filled', 'MarkerFaceAlpha', .7); clim(ax6, [.2 .8])
    text(ax6, head_p, ratio_b, arrayfun(@(b) sprintf(' %d', b), 1:nb, 'UniformOutput', false), 'FontSize', 6)
    set(ax6,'XScale','log','YScale','log'); yline(ax6, 1, 'k:'); yline(ax6, slope, 'r-')
    xlabel(ax6,'heading path per bout (rad, log)'); ylabel(ax6,'bump / heading per bout (log)'); grid(ax6,'on')
    title(ax6, 'per-bout ratio: small bouts vs large bouts')
    % slope vs floor
    ax7 = nexttile(tl, 9); hold(ax7,'on')
    plot(ax7, floors, slope_floor, 'ro-', 'MarkerFaceColor','r', 'LineWidth', 1.5)
    for k = 1:numel(floors), text(ax7, floors(k), slope_floor(k)+.05, sprintf('n=%d', n_floor(k)), 'FontSize', 7, 'HorizontalAlignment','center'); end
    yline(ax7, 1, 'k:'); xlabel(ax7,'drop bouts with heading path < x rad'); ylabel(ax7,'pooled slope'); ylim(ax7,[0 max(2.5, max(slope_floor)+.3)]); grid(ax7,'on')
    title(ax7, sprintf('slope vs heading-path floor | unweighted %.2f, ratio of sums %.2f', slope_unw, ratio_sums))
    % net displacement view
    ax8 = nexttile(tl, 12); hold(ax8,'on')
    scatter(ax8, abs(net_head), abs(net_bump), 15 + 8*sqrt(dur), rho_b, 'filled', 'MarkerFaceAlpha', .7); clim(ax8, [.2 .8])
    text(ax8, abs(net_head), abs(net_bump), arrayfun(@(b) sprintf(' %d', b), 1:nb, 'UniformOutput', false), 'FontSize', 6)
    xm2 = max([abs(net_head); abs(net_bump); 1]); plot(ax8, [0 xm2], [0 xm2], 'k:'); axis(ax8,'square')
    sl_net = (w.*abs(net_head)) \ (w.*abs(net_bump)); plot(ax8, [0 xm2], sl_net*[0 xm2], 'r-', 'LineWidth', 1.5)
    xlabel(ax8,'|net heading change| per bout (rad)'); ylabel(ax8,'|net bump change| per bout (rad)')
    title(ax8, sprintf('NET displacement instead of path length: slope %.2f', sl_net))
    exportgraphics(fig, fullfile(export_dir, sprintf('bout_diag_%s_trial%d.png', regexprep(fly_tag,'[/ ]','_'), trial_num(i))), 'Resolution', 110)
end

%% helpers
function [fid, tnum] = walk_meta_fly_trial(meta_path)
    parts = strsplit(meta_path,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    if ~(isfield(trial,'atp') && isfield(trial.atp,'d') && ~isempty(trial.atp.d)); tf = false; return; end
    stims = logical(trial.ft.stims(:)); onsets = find(diff([false;stims])==1);
    if isempty(onsets); tf = false; return; end
    xb = trial.ft.xb; xf = trial.ft.xf; t_onset = xf(onsets); atp_sig = sum(trial.atp.d,1);
    t_common = linspace(-peri_win,peri_win,101);
    atp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets), atp_aligned(s,:) = interp1(xb - t_onset(s),atp_sig,t_common,'linear',nan); end
    atp_peri = mean(atp_aligned,1,'omitnan');
    base_idx = t_common>=base_win(1) & t_common<base_win(2); peak_idx = t_common>peak_win(1) & t_common<=peak_win(2);
    tf = (median(atp_peri(peak_idx),'omitnan')-median(atp_peri(base_idx),'omitnan')) > peak_factor*std(atp_peri(base_idx),'omitnan');
end
