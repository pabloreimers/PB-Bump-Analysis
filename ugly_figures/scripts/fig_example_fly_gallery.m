%% fig_example_fly_gallery
% Gallery of every fly in the Figure 3 walking dataset as a candidate
% example for fig_example_flies_pb_heatmaps.m. One panel per fly (one
% figure per genotype): the per-walking-bout scatter of bump path length
% vs heading path length for that fly's dark trials, pre-ATP (gray) and
% post-ATP (red), with Figure 3's sqrt(duration)-weighted through-origin
% slopes. The title carries what you need to pick:
%   r_CL  = closed-loop tracking quality, pre-ATP (circ corr of -cue vs mu
%           on rho > 0.2 AND cue moving; mean over CL-pre trials, with the
%           per-trial values in brackets)
%   mob   = pooled dark bout slope pre -> post (Figure 3's metric)
%   rho   = mean rho over dark trials pre -> post
%   nb    = pooled dark bout counts pre / post
% Panels are sorted by r_CL within genotype. Flies without a detected
% pulse trial are omitted. Bout detection and slope are the same as
% Figure_3_claude.m (walk_trial_bouts, min 15 bouts per pool).

%% paths + load
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root,'circ_stats'))
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);

%% labels
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id); n_flies = numel(fly_list);

%% bouts per trial + per-trial tracking r + mean rho
bout_mu = cell(n,1); bout_cue = bout_mu; bout_dur = bout_mu; trial_rho = nan(n,1); track_r = nan(n,1);
for i = 1:n
    [bout_mu{i}, bout_cue{i}, bout_dur{i}] = walk_trial_bouts(all_data(i), 60, .25, 30, 30);
    trial_rho(i) = mean(all_data(i).im.rho);
    if ~is_dark(i), track_r(i) = trial_tracking_r(all_data(i), .2, .2); end
end

%% per fly summary
F = struct('f',{},'geno',{},'short',{},'r',{},'r_list',{},'pre_idx',{},'post_idx',{},'cl_pre_idx',{},'mob_pre',{},'mob_post',{},'rho_pre',{},'rho_post',{},'nb_pre',{},'nb_post',{});
for f = 1:n_flies
    idx = find(fly_ix==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    p = find(has_pulse(idx),1); if isempty(p), continue, end
    pre = idx(1:p-1); post = idx(p+1:end); post = post(~has_pulse(post));
    s.f = f; s.geno = 'empty>p2x2'; if is_lpsp(idx(1)), s.geno = 'lpsp>p2x2'; end
    s.short = strrep(fly_list{f},'Z:/pablo/lpsp_p2x2_walking/','');
    s.cl_pre_idx = pre(~is_dark(pre)); s.r_list = track_r(s.cl_pre_idx); s.r = mean(s.r_list,'omitnan');
    s.pre_idx = pre(is_dark(pre)); s.post_idx = post(is_dark(post));
    s.mob_pre = pooled_slope(bout_cue(s.pre_idx), bout_mu(s.pre_idx), bout_dur(s.pre_idx), 15);
    s.mob_post = pooled_slope(bout_cue(s.post_idx), bout_mu(s.post_idx), bout_dur(s.post_idx), 15);
    s.rho_pre = mean(trial_rho(s.pre_idx)); s.rho_post = mean(trial_rho(s.post_idx));
    s.nb_pre = numel(vertcat(bout_mu{s.pre_idx})); s.nb_post = numel(vertcat(bout_mu{s.post_idx}));
    F(end+1) = s; %#ok<AGROW>
end

%% gallery, one figure per genotype
genos = {'empty>p2x2','lpsp>p2x2'};
for gi = 1:2
    G = F(strcmp({F.geno}, genos{gi})); [~,o] = sort([G.r], 'descend', 'MissingPlacement','last'); G = G(o);
    nG = numel(G); n_col = 4; n_row = ceil(nG/n_col);
    fig = figure(60+gi); clf; set(fig,'Color','w','Units','pixels','Position',[20 20 1700 330*n_row])
    tl = tiledlayout(n_row, n_col, 'TileSpacing','compact','Padding','compact');
    for k = 1:nG
        s = G(k); ax = nexttile(tl,k); hold(ax,'on')
        cols = {[.35 .35 .35], [.85 .2 .2]}; labs = {'dark PRE','dark POST'}; sets = {s.pre_idx, s.post_idx}; mobs = [s.mob_pre, s.mob_post];
        xmax = 1;
        for c = 1:2
            ii = sets{c}; if isempty(ii), continue, end
            cue = vertcat(bout_cue{ii}); mu = vertcat(bout_mu{ii}); dur = vertcat(bout_dur{ii});
            scatter(ax, cue, mu, 10 + 3*sqrt(dur), cols{c}, 'filled', 'MarkerFaceAlpha', .5)
            xmax = max([xmax; cue]);
            if ~isnan(mobs(c)), plot(ax, [0 xmax], mobs(c)*[0 xmax], '-', 'Color', cols{c}, 'LineWidth', 1.5); end
        end
        plot(ax, [0 xmax], [0 xmax], 'k:')
        xlim(ax, [0 xmax*1.05]); ylim(ax, [0 max(xmax*1.05, 1)]); axis(ax,'square'); set(ax,'FontSize',7)
        title(ax, sprintf('#%d %s\nr_{CL}=%.2f %s\nmob %.2f -> %.2f | rho %.2f -> %.2f | nb %d/%d', k, s.short, s.r, mat2str(round(s.r_list',2)), s.mob_pre, s.mob_post, s.rho_pre, s.rho_post, s.nb_pre, s.nb_post), 'FontSize', 7.5)
        if mod(k-1,n_col)==0, ylabel(ax,'bump path / bout (rad)'); end
        if k > nG-n_col, xlabel(ax,'heading path / bout (rad)'); end
    end
    title(tl, sprintf('%s candidate example flies, sorted by closed-loop tracking r (pre-ATP). gray = dark pre bouts, red = dark post bouts, lines = pooled sqrt(dur)-weighted slopes (Figure 3 metric), dotted = 1:1', genos{gi}), 'Interpreter','none','FontSize',10)
    exportgraphics(fig, fullfile(export_dir, sprintf('fig_example_fly_gallery_%s.png', strrep(genos{gi},'>','_'))), 'Resolution', 120)
    fprintf('\n== %s ==\n%-3s %-18s %-6s %-7s %-7s %-7s %-7s %s\n', genos{gi}, '#','fly','r_CL','mobPre','mobPost','rhoPre','rhoPost','r per CL-pre trial');
    for k = 1:nG, s = G(k); fprintf('%-3d %-18s %-6.2f %-7.2f %-7.2f %-7.2f %-7.2f %s\n', k, s.short, s.r, s.mob_pre, s.mob_post, s.rho_pre, s.rho_post, mat2str(round(s.r_list',2))); end
end

%% helpers
function r = trial_tracking_r(trial, rho_thr, cue_speed_thr)
    xb = trial.ft.xb; xf = trial.ft.xf;
    cue_u = interp1(xf, unwrap(-trial.ft.cue), xb, 'linear', 'extrap');
    cue_sp = abs(gradient(smoothdata(cue_u,'movmean',5), xb));
    good = trial.im.rho > rho_thr & cue_sp > cue_speed_thr;
    r = nan; if sum(good) > 100, r = circ_corrcc(mod(cue_u(good),2*pi), mod(trial.im.mu(good),2*pi)); end
end

function s = pooled_slope(cue_c, mu_c, dur_c, min_bouts)
    c = vertcat(cue_c{:}); m = vertcat(mu_c{:}); d = vertcat(dur_c{:});
    s = nan; if numel(c) >= min_bouts, w = sqrt(d); s = (w.*c) \ (w.*m); end
end

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

function [mov_mu, mov_cue, dur] = walk_trial_bouts(trial, smooth_window, turn_thresh, max_gap_frames, min_walking_frames)
    xf = trial.ft.xf; dt = median(diff(xf));
    mu = interp1(trial.ft.xb,unwrap(trial.im.mu),xf,'linear','extrap');
    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth = smoothdata(mu,'gaussian',smooth_window);
    is_walking = abs(r_speed_smooth) > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);
    d = diff([false;is_walking;false]); bout_starts = find(d==1); bout_ends = find(d==-1)-1;
    mov_mu = nan(numel(bout_starts),1); mov_cue = mov_mu; dur = mov_mu;
    for b = 1:numel(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(r_speed_smooth(bout_starts(b):bout_ends(b))),'omitnan')*dt;
        dur(b)     = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end
