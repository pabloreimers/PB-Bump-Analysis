%% fig_flagged_flies_heatmaps
% For the empty>p2x2 control flies whose closed-loop pre-ATP traces looked
% decoupled in fig_cl_pre_traces_sweep.m, show the underlying data: the
% z-scored PB heatmap (im.z, 32 wedges), the extracted bump (im.mu),
% the cue (-ft.cue) and the fly's heading (cumtrapz r_speed), all
% unwrapped and re-zeroed, plus rho. Left column = the whole trial, right
% column = the same 60 s snippet the trace grids used (first closed-loop
% pre-ATP trial, window with the most walking). Nothing is smoothed or
% masked here beyond what the dataset carries.

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports','cl_pre_traces');
if ~isfolder(export_dir), mkdir(export_dir), end
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% which flies
flagged = {'20260728/fly 18', '20260721/fly 2', '20260727/fly 3'};
snip_len = 60;

%% labels
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end

for k = 1:numel(flagged)
    idx = find(contains(fly_id, flagged{k})); [~,o] = sort(trial_num(idx)); idx = idx(o);
    p = find(has_pulse(idx),1); pre = idx(1:p-1); cl = pre(~is_dark(pre)); i0 = cl(1);
    tr = all_data(i0); xf = tr.ft.xf(:); xb = tr.ft.xb(:); dt = median(diff(xf));
    mu_u  = unwrap(tr.im.mu(:));
    cue_u = unwrap(-tr.ft.cue(:));
    r = tr.ft.r_speed(:); r(isnan(r)) = 0; head = cumtrapz(xf, r);
    r_s = smoothdata(tr.ft.r_speed(:),'gaussian',60); w = abs(r_s) > .25; w = imclose(w,ones(31,1)); w = bwareaopen(w,30);
    wlen = round(snip_len/dt); cover = movmean(double(w), wlen, 'Endpoints','discard'); [~,k0] = max(cover);
    t0 = xf(k0); t1 = t0 + snip_len;

    fig = figure(50+k); clf; set(fig,'Color','w','Units','pixels','Position',[30 30 1800 850])
    tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
    for col = 1:2
        if col == 1, tr_lim = [xf(1) xf(end)]; else, tr_lim = [t0 t1]; end
        kb = xb >= tr_lim(1) & xb <= tr_lim(2); kf = xf >= tr_lim(1) & xf <= tr_lim(2);
        % re-zero at first sample of the window
        mu0 = interp1(xb, mu_u, tr_lim(1), 'linear', 'extrap'); c0 = interp1(xf, cue_u, tr_lim(1)); h0 = interp1(xf, head, tr_lim(1));

        % heatmap with bump + cue overlaid as wedge index
        ax = nexttile(tl, col); hold(ax,'on')
        imagesc(ax, xb(kb), 1:32, tr.im.z(:,kb), [-1.5 3]); colormap(ax, parula)
        alpha = tr.im.alpha(:);
        for h = 1:2
            rows = (h-1)*16 + (1:16);
            y = interp1(alpha(rows), rows, tr.im.mu(kb), 'linear', 'extrap'); y(abs(diff([y(1);y]))>8) = nan;
            plot(ax, xb(kb), y, '-', 'Color', [1 1 1 .8], 'LineWidth', .7)
            cw = mod(-tr.ft.cue(kf)+pi, 2*pi)-pi; yc = interp1(alpha(rows), rows, cw, 'linear', 'extrap'); yc(abs(diff([yc(1);yc]))>8) = nan;
            plot(ax, xf(kf), yc, '-', 'Color', [0 1 1 .7], 'LineWidth', .6)
        end
        yline(ax, 16.5, 'w:'); set(ax,'YDir','normal','YTick',[1 16 17 32]); xlim(ax, tr_lim); ylim(ax, [.5 32.5])
        ylabel(ax, 'PB wedge'); if col == 1, title(ax, sprintf('full trial %d (%.0f s): im.z, white = im.mu, cyan = -cue', trial_num(i0), xf(end)-xf(1))); else, title(ax, sprintf('%d s snippet from t = %.0f s', snip_len, t0)); end

        % unwrapped traces
        ax = nexttile(tl, 2+col); hold(ax,'on')
        plot(ax, xf(kf), head(kf)-h0, 'k-', 'LineWidth', 1, 'DisplayName', 'heading (\int r\_speed)')
        plot(ax, xf(kf), cue_u(kf)-c0, '-', 'Color', [0 .7 .85], 'LineWidth', 1, 'DisplayName', 'cue (-ft.cue)')
        plot(ax, xb(kb), mu_u(kb)-mu0, '-', 'Color', [.85 .33 .1], 'LineWidth', 1, 'DisplayName', 'bump (im.mu)')
        xlim(ax, tr_lim); ylabel(ax, 'unwrapped angle (rad, re-zeroed)'); grid(ax,'on')
        legend(ax, 'Location','best', 'FontSize', 8, 'Box','off')

        % rho + walking
        ax = nexttile(tl, 4+col); hold(ax,'on')
        area(ax, xf(kf), double(w(kf))*0.1, 'FaceColor', [.85 .85 .85], 'EdgeColor','none', 'DisplayName', 'walking (bout mask)')
        plot(ax, xb(kb), tr.im.rho(kb), '-', 'Color', [.4 .2 .6], 'LineWidth', .9, 'DisplayName', 'rho')
        yline(ax, .3, ':r', 'HandleVisibility','off')
        xlim(ax, tr_lim); ylim(ax, [0 max(1.2, max(tr.im.rho(kb)))]); ylabel(ax, 'rho (dotted = 0.3)'); xlabel(ax, 'time (s)'); grid(ax,'on')
        legend(ax, 'Location','best', 'FontSize', 8, 'Box','off')
    end
    title(tl, sprintf('%s (empty>p2x2), closed-loop pre-ATP trial %d -- pattern %s', flagged{k}, trial_num(i0), tr.ft.pattern), 'Interpreter','none')
    exportgraphics(fig, fullfile(export_dir, sprintf('heatmap_%s.png', regexprep(flagged{k},'[/ ]','_'))), 'Resolution', 110)
    fprintf('%s: trial %d, %d closed-loop pre trials, snippet from %.0f s; mean rho %.2f; corr(diff heading, diff cue) = %.2f\n', flagged{k}, trial_num(i0), numel(cl), t0, mean(tr.im.rho), corr(diff(head), diff(cue_u)));
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
