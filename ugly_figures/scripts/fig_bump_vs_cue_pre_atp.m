%% fig_bump_vs_cue_pre_atp
% QC of the bump-extraction method behind Figure_3_claude.m's bump-mobility
% row: for every CLOSED-LOOP trial recorded BEFORE a fly received any ATP
% stimulation, scatter the bump position (im.mu) against the cue position
% (-ft.cue, sign flipped to mu's convention -- see lpsp_p2x2_walking_script.m
% section 1). If extraction is good, each panel should be a tight diagonal
% band (one arbitrary offset per fly, wrapping at +/-pi).
%
% Dataset: .data/lpsp_p2x2_walking_20260728_post0708.mat (the subset Figure
% 3 uses). Trial selection, per fly:
%   * chronological order by trial number (walk_meta_fly_trial)
%   * "pre-ATP" = every trial strictly before that fly's FIRST detected
%     perturbation trial (walk_trial_has_pulse, same detector and settings
%     as Figure_3_claude.m / lpsp_p2x2_walking_script.m section 3)
%   * closed loop = ft.pattern is the brightbar pattern (not 'background')
% Flies with no detected pulse trial are skipped (we can't know what is
% "pre" for them), and that is reported in the command window.
%
% Points are colored by im.rho (bump vector strength) so low-confidence
% samples can be told apart from the real tracking.

%% paths
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
addpath(fullfile(repo_root,'circ_stats'))
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
if ~isfolder(export_dir); mkdir(export_dir); end

%% load (reuse an already-loaded copy in the base workspace if present)
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n = numel(all_data);
fprintf('loaded %d trials from %s\n', n, source_file);

%% label trials: fly, trial number, genotype, lighting, pulse
is_lpsp = arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data);
is_dark = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data);
fly_id = cell(n,1); trial_num = nan(n,1); has_pulse = false(n,1);
for i = 1:n
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), 5, [-5,0], [0,3], 3);
end
[fly_list,~,fly_ix] = unique(fly_id);
n_flies = numel(fly_list);

%% select pre-ATP closed-loop trials
sel = []; % trial indices, grouped by fly in chronological order
sel_fly = [];
n_skipped_flies = 0;
for f = 1:n_flies
    idx = find(fly_ix==f);
    [~,o] = sort(trial_num(idx)); idx = idx(o);
    p = find(has_pulse(idx),1);
    if isempty(p)
        n_skipped_flies = n_skipped_flies+1;
        fprintf('  skipping %s: no detected pulse trial\n', fly_list{f});
        continue
    end
    pre = idx(1:p-1);
    pre = pre(~is_dark(pre));
    sel = [sel; pre(:)]; sel_fly = [sel_fly; repmat(f,numel(pre),1)]; %#ok<AGROW>
end
fprintf('%d pre-ATP closed-loop trials from %d flies (%d flies skipped)\n', numel(sel), numel(unique(sel_fly)), n_skipped_flies);

%% per-trial bump vs cue (on the imaging clock, xb)
rho_thresh = .1; % same as walking script section 1; only used to report a circ corr, all points are plotted
cue_b = cell(numel(sel),1); mu_b = cue_b; rho_b = cue_b; ccorr = nan(numel(sel),1);
for k = 1:numel(sel)
    i = sel(k);
    xb = all_data(i).ft.xb; xf = all_data(i).ft.xf;
    cue = interp1(xf, unwrap(-all_data(i).ft.cue), xb, 'linear', 'extrap'); % cue downsampled onto imaging clock; sign flipped to mu's convention
    cue_b{k} = mod(cue+pi,2*pi)-pi;
    mu_b{k}  = mod(all_data(i).im.mu+pi,2*pi)-pi;
    rho_b{k} = all_data(i).im.rho;
    good = rho_b{k} > rho_thresh;
    ccorr(k) = circ_corrcc(cue_b{k}(good), mu_b{k}(good));
end

%% plot: one subpanel per trial
n_pan = numel(sel);
n_col = ceil(sqrt(n_pan*1.4)); n_row = ceil(n_pan/n_col);
fig = figure(31); clf
set(fig,'Color','w','Units','pixels','Position',[50 50 min(240*n_col,1900) min(220*n_row,1000)])
tl = tiledlayout(n_row,n_col,'TileSpacing','compact','Padding','compact');
for k = 1:n_pan
    i = sel(k); f = sel_fly(k);
    ax = nexttile(tl); hold(ax,'on')
    scatter(ax, cue_b{k}, mu_b{k}, 4, rho_b{k}, 'filled', 'MarkerFaceAlpha', .5)
    xlim(ax,[-pi,pi]); ylim(ax,[-pi,pi]); axis(ax,'square')
    xticks(ax,[-pi,0,pi]); yticks(ax,[-pi,0,pi]);
    xticklabels(ax,{'-\pi','0','\pi'}); yticklabels(ax,{'-\pi','0','\pi'})
    clim(ax,[0,.6])
    geno = 'empty'; if is_lpsp(i); geno = 'lpsp'; end
    date_str = regexp(fly_list{f},'\d{8}','match','once');
    fnum = regexp(fly_list{f},'fly (\d+)','tokens','once');
    title(ax, sprintf('%s fly %s (%s) tr%d  r=%.2f', date_str(5:end), fnum{1}, geno, trial_num(i), ccorr(k)), 'FontSize',8, 'FontWeight','normal')
    if is_lpsp(i); ax.Title.Color = [.8 .2 .2]; end
    set(ax,'FontSize',7)
end
cb = colorbar(ax); cb.Layout.Tile = 'east'; cb.Label.String = 'rho (bump strength)';
xlabel(tl,'cue position (-ft.cue, rad)'); ylabel(tl,'bump position (im.mu, rad)')
title(tl, sprintf('bump vs cue, pre-ATP closed-loop trials (%d trials, %d flies) -- red titles = lpsp>p2x2; r = circ corr (rho>%.1f)', n_pan, numel(unique(sel_fly)), rho_thresh), 'FontSize',10)

%% export
exportgraphics(fig, fullfile(export_dir,'fig_bump_vs_cue_pre_atp.png'), 'Resolution',150)
fprintf('exported to %s\n', fullfile(export_dir,'fig_bump_vs_cue_pre_atp.png'))

%% helpers (ported from Figure_3_claude.m)
function [fid, tnum] = walk_meta_fly_trial(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor)
    if ~(isfield(trial,'atp') && isfield(trial.atp,'d') && ~isempty(trial.atp.d)); tf = false; return; end
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    if isempty(onsets); tf = false; return; end
    xb = trial.ft.xb; xf = trial.ft.xf;
    t_onset = xf(onsets);
    atp_sig = sum(trial.atp.d,1);
    t_common = linspace(-peri_win,peri_win,101);
    atp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets)
        atp_aligned(s,:) = interp1(xb - t_onset(s),atp_sig,t_common,'linear',nan);
    end
    atp_peri = mean(atp_aligned,1,'omitnan');
    base_idx = t_common>=base_win(1) & t_common<base_win(2);
    peak_idx = t_common>peak_win(1) & t_common<=peak_win(2);
    tf = (median(atp_peri(peak_idx),'omitnan')-median(atp_peri(base_idx),'omitnan')) > peak_factor*std(atp_peri(base_idx),'omitnan');
end
