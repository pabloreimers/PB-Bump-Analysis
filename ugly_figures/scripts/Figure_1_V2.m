%% Figure_1_V2
% Figure 1, V2 draft: LPsP background + unilateral LPsP activation rotates
% the EPG bump (dangling fly, open-loop P2X2 iontophoresis).
%
% Layout (3 x 3; exact placement is left to Illustrator):
%   row 1: [compass cartoon (placeholder)] [compass inputs (placeholder)]
%          [LPsP morphology, FlyWire whole brain, posterior view]
%   row 2: [dangling-fly schematic (placeholder)]
%          [example LPsP>P2X2 stim] [paired example empty>P2X2 stim]
%   row 3: [average bump position, LPsP>P2X2 vs empty>P2X2]
%          [net bump displacement per stim, 4-15 s, left-stim convention,
%           one dot per fly] [same, raw sign, split by stimulated side]
%
% Not in this figure, saved for the last (long-term LPsP activation) figure:
% after stims, LPsP>P2X2 bumps follow the open-loop bright-bar sweeps more
% (fig_p2x2_sweep_tracking.m, fig_p2x2_side_effects.m; see
% ugly_figures/V2_NOTES.md).
%
% Conventions for this figure set: time is shown as a SCALE BAR, never as a
% time axis, so it survives resizing in Illustrator. The two example stims
% are oriented VERTICALLY: time runs top -> bottom, PB angle runs left ->
% right. Right of each heatmap: a red bar over the ATP TTL (0-2 s) above the
% 5 s time scale bar; left: a short "Time" arrow pointing down. No x axis
% (the mean image above shows the PB).
%
% Example mean images (EPG blue / ATP red, above each heatmap): made from the
% raw per-plane registered movies by trial_mean_rgb (ported from
% Figure_3_claude.m). The two example trials' files are copied from the Z:
% drive into Dropbox (raw_dir below, one folder per trial, flat:
% imagingData_trial001.mat + mask.mat). !! If a trial's raw files are
% missing, that example FALLS BACK to a crop of V1 Figure 3's exported PNG
% (ugly_figures/all_figs/panels/Fig3/F3_r12.png) with a warning -- low-res,
% glomerulus-1 end trimmed. Don't ship a figure that printed that warning.
%
% Exports go to the repo (ugly_figures/all_figs/V2, exports/V2) AND to
% Dropbox (dropbox_dir below), when that folder exists on this machine.
%
% Morphology panel: PNG rendered by fig_lpsp_morphology_flywire.py (FlyWire
% v783, public data; see that script's header for the env and the
% PB/LAL-only shading rationale). Run it first if the PNG is missing:
%   ~/venvs/flywire/bin/python ugly_figures/scripts/fig_lpsp_morphology_flywire.py
%
% P2X2 data + example picks: ported UNCHANGED from Figure_3_claude.m (V1
% Figure 3) -- dataset lpsp_p2x2_reredo_20260518.mat, per-stim ejection
% confirmation (get_trial_stim_events), stim side called from each trial's
% most confident ejection, LPsP example = rank 3 on the reference side by
% post-ejection bump-vs-heading divergence (8-15 s), empty example = same
% side, atp amplitude within +-25% of the LPsP example's, ranked by atp
% trace smoothness, rank 3. See Figure_3_claude.m for the reasoning behind
% each of these choices. The average-trace panel is V1 Figure 3's row 2
% left panel: every confirmed ejection, flipped onto the example's side,
% averaged per fly first, then mean +/- SEM across flies.

%% paths
clearvars % stale variables from the fig_p2x2_* scratch scripts (same names, different shapes) break this script
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath')))); % ugly_figures/scripts -> repo root
addpath(fullfile(repo_root,'circ_stats'))
data_dir = fullfile(repo_root,'.data'); % Windows machine keeps datasets in .data/, the Mac in data/
if ~isfolder(data_dir); data_dir = fullfile(repo_root,'data'); end
export_dir   = fullfile(repo_root,'ugly_figures','exports','V2');
all_figs_dir = fullfile(repo_root,'ugly_figures','all_figs','V2');
morph_png    = fullfile(repo_root,'ugly_figures','exports','V2','fig1_anatomy','fig1_lpsp_morphology_posterior.png');
lpsp_dropbox = fullfile(getenv('HOME'),'Dropbox-HarvardUniversity','Pablo Reimers','LPsP');
raw_dir      = lpsp_dropbox;                      % example trials' raw movies + masks, one folder per trial
dropbox_dir  = fullfile(lpsp_dropbox,'V2');       % shared copy of the V2 figures

%% load lpsp_p2x2_reredo data
source_file = 'lpsp_p2x2_reredo_20260518.mat';
tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:); clear tmp
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype + fly ID per trial (folder "fly N _ <genotype>" in all_data.meta)
genotype = cell(n_trials,1);
fly_id   = cell(n_trials,1);
for i = 1:n_trials
    fname = trial_fly_folder_name(all_data(i).meta);
    if contains(fname,'lpsp','IgnoreCase',true)
        genotype{i} = 'lpsp>p2x2';
    elseif contains(fname,'empty','IgnoreCase',true)
        genotype{i} = 'empty>p2x2';
    else
        genotype{i} = '';
    end
    fly_id{i} = trial_fly_id(all_data(i).meta);
end

%% confirmed-ejection stim events, every trial (as Figure_3_claude.m)
peri_win    = 5;
base_win    = [-peri_win,0];
peak_win    = [0,3];
peak_factor = 3;
t_grid = linspace(-5,30,351);

trial_events = cell(n_trials,1);
for i = 1:n_trials
    trial_events{i} = get_trial_stim_events(all_data(i), t_grid, peri_win, base_win, peak_win, peak_factor);
end
n_ejections_confirmed_all = cellfun(@numel,trial_events);

%% stim side per trial, from its most confident ejection (as Figure_3_claude.m)
trial_is_right = nan(n_trials,1);
for i = 1:n_trials
    ev = trial_events{i};
    if isempty(ev); continue; end
    rfs = [ev.right_frac];
    [~,most_confident] = max(abs(rfs-0.5));
    trial_is_right(i) = rfs(most_confident) > 0.5;
end

%% LPsP example stim (as Figure_3_claude.m: rank 3 on the reference side)
settle_win = [8,15];
rho_gate   = 0.25;
atp_peak_amp_thresh = 2.0;
settle_idx = t_grid>settle_win(1) & t_grid<=settle_win(2);

lpsp_rows = find(strcmp(genotype,'lpsp>p2x2'));
cand_trial = []; cand_stim = []; cand_score = []; cand_side = [];
for i = lpsp_rows(:)'
    if n_ejections_confirmed_all(i) < 2; continue; end
    ev = trial_events{i};
    for s = 1:numel(ev)
        if ev(s).atp_peak_amp < atp_peak_amp_thresh; continue; end
        if mean(ev(s).rho_grid(settle_idx),'omitnan') <= rho_gate; continue; end
        cand_trial(end+1) = i; %#ok<AGROW>
        cand_stim(end+1)  = s; %#ok<AGROW>
        cand_score(end+1) = abs(median(ev(s).mu_rel(settle_idx) - ev(s).head_rel(settle_idx),'omitnan')); %#ok<AGROW>
        cand_side(end+1)  = trial_is_right(i); %#ok<AGROW>
    end
end
[~,best] = max(cand_score);
reference_is_right = logical(cand_side(best));
if reference_is_right; ref_side_str = 'right'; else; ref_side_str = 'left'; end
ref_pool = find(cand_side==reference_is_right);
[~,o] = sort(cand_score(ref_pool),'descend'); ref_pool = ref_pool(o);
example_rank_ref = 3;
best = ref_pool(example_rank_ref);
example_i = cand_trial(best); example_s = cand_stim(best);

%% paired empty example (as Figure_3_claude.m: same side, amplitude-matched, smoothest, rank 3)
emp_amp_tol = 0.25;
example_rank_empty = 3;
ref_amp = trial_events{example_i}(example_s).atp_peak_amp;
empty_rows = find(strcmp(genotype,'empty>p2x2'));
emp_trial = []; emp_stim = []; emp_amp_err = []; emp_rough = [];
for i = empty_rows(:)'
    if n_ejections_confirmed_all(i) < 2 || trial_is_right(i)~=reference_is_right; continue; end
    ev = trial_events{i};
    for s = 1:numel(ev)
        if ev(s).atp_peak_amp < atp_peak_amp_thresh; continue; end
        if mean(ev(s).rho_grid(settle_idx),'omitnan') <= rho_gate; continue; end
        emp_trial(end+1)   = i; %#ok<AGROW>
        emp_stim(end+1)    = s; %#ok<AGROW>
        emp_amp_err(end+1) = abs(ev(s).atp_peak_amp - ref_amp); %#ok<AGROW>
        emp_rough(end+1)   = atp_trace_roughness(all_data(i), ev(s).onset_t); %#ok<AGROW>
    end
end
emp_in_band  = find(emp_amp_err <= emp_amp_tol*ref_amp);
emp_out_band = find(emp_amp_err >  emp_amp_tol*ref_amp);
[~,o_in]  = sort(emp_rough(emp_in_band),'ascend');
[~,o_out] = sort(emp_amp_err(emp_out_band),'ascend');
emp_pool = [emp_in_band(o_in), emp_out_band(o_out)];
best2 = emp_pool(example_rank_empty);
example_i2 = emp_trial(best2); example_s2 = emp_stim(best2);

example_idx      = [example_i, example_i2];
example_stim_idx = [example_s, example_s2];
example_geno     = {'LPsP>P2X2','Empty>P2X2'};
for c = 1:2
    fprintf('example %s: fly %s, trial "%s", stim %d (%s-side stim), atp peak %.2f dF/F\n', example_geno{c}, ...
        fly_id{example_idx(c)}, meta_display(all_data(example_idx(c)).meta), example_stim_idx(c), ref_side_str, ...
        trial_events{example_idx(c)}(example_stim_idx(c)).atp_peak_amp);
end

%% average bump-position trace per genotype (as Figure_3_claude.m row 2 left)
trace_pre_s = -2; trace_post_s = 20;
trace_idx = t_grid>=trace_pre_s & t_grid<=trace_post_s;
t_trace = t_grid(trace_idx);
geno_order = {'empty>p2x2','lpsp>p2x2'};
geno_label = {'Empty>P2X2','LPsP>P2X2'};

geno_mean = nan(2,numel(t_trace)); geno_sem = nan(2,numel(t_trace)); geno_n = zeros(1,2);
for gi = 1:2
    rows = find(strcmp(genotype,geno_order{gi}));
    these_flies = unique(fly_id(rows));
    M = nan(numel(these_flies),numel(t_trace));
    for f = 1:numel(these_flies)
        stim_traces = [];
        for i = rows(strcmp(fly_id(rows),these_flies{f}))'
            ev = trial_events{i};
            for s = 1:numel(ev)
                mu_rel = ev(s).mu_rel(trace_idx);
                if trial_is_right(i) ~= reference_is_right; mu_rel = -mu_rel; end % flip onto the example's side
                stim_traces(end+1,:) = mu_rel; %#ok<AGROW>
            end
        end
        if ~isempty(stim_traces); M(f,:) = mean(stim_traces,1,'omitnan'); end
    end
    M = M(~all(isnan(M),2),:);
    geno_n(gi) = size(M,1);
    geno_mean(gi,:) = mean(M,1,'omitnan');
    geno_sem(gi,:)  = std(M,[],1,'omitnan') ./ sqrt(sum(~isnan(M),1));
end
fprintf('average trace: n=%d empty flies, n=%d LPsP flies (%s-flip)\n', geno_n(1), geno_n(2), ref_side_str);

%% net bump displacement per stim, 4-15 s after onset (row 3, panels 2-3)
% Settled in fig_p2x2_displacement_per_stim.m (see its header for the
% exploration): per confirmed stim, net displacement = mean unwrapped bump
% position over [disp_t1-1, disp_t1] minus over [disp_t0-1, disp_t0] s re
% onset. Starting at 4 s (chosen 2026-11-23, was 2 s) keeps the window
% outside the ATP-evoked bright spot: during the 2 s pulse a broad
% brightening at the pipette site (both genotypes) pulls the PVA phase into
% it, and as it fades over ~2-4 s the phase can sweep across the PB with no
% bump actually moving. The effect of interest is the bump's movement after
% that. The window ends well before the bar starts sweeping again (18.6 s; stims whose
% window would reach a sweep are dropped -- none are). Bump position is
% NaN'd where rho < disp_rho_min, unwrapped across the rest and linearly
% interpolated through the gaps. Panel 2: signed onto the LEFT-stim
% convention (right-side trials flipped), per fly mean over its stims.
% Panel 3: raw sign, per fly mean over its left-side and right-side trials
% separately (lines join the same fly).
% OPEN ISSUE (ugly_figures/V2_NOTES.md, "Come back to"): in some Empty flies
% the PVA phase still wanders through static bright glomeruli when there's
% no clear bump; rho doesn't flag it. A bump-quality filter may help.
disp_t0 = 4; disp_t1 = 15; disp_rho_min = 0.1;
fly_disp = cell(1,2); fly_disp_side = cell(1,2); n_disp_stims = cell(1,2); n_sweep_drop = 0;
for gi = 1:2
    rows = find(strcmp(genotype,geno_order{gi}) & ~isnan(trial_is_right));
    these_flies = unique(fly_id(rows));
    fd = nan(numel(these_flies),1); fds = nan(numel(these_flies),2); ns = zeros(numel(these_flies),1);
    for f = 1:numel(these_flies)
        d_all = []; d_side = {[],[]};
        for i = rows(strcmp(fly_id(rows),these_flies{f}))'
            tr = all_data(i); xb = tr.ft.xb(:); xf = tr.ft.xf(:); dtf = median(diff(xf));
            m = masked_unwrapped_mu(tr.im.mu(:), tr.im.rho(:), xb, disp_rho_min);
            cue_v = [0; diff(unwrap(tr.ft.cue(:)))] ./ [1; diff(xf)];
            sweeping = bwareaopen(abs(smoothdata(cue_v,'movmean',round(.5/dtf))) > 0.05, round(1/dtf));
            for e = 1:numel(trial_events{i})
                t0 = trial_events{i}(e).onset_t;
                next_sweep = xf(find(xf > t0 & sweeping, 1));
                if ~isempty(next_sweep) && t0 + disp_t1 > next_sweep; n_sweep_drop = n_sweep_drop + 1; continue; end
                a = m(xb >= t0+disp_t0-1 & xb < t0+disp_t0); b = m(xb >= t0+disp_t1-1 & xb < t0+disp_t1);
                if isempty(a) || isempty(b); continue; end
                d = mean(b) - mean(a);
                d_all(end+1) = d * (1 - 2*(trial_is_right(i) ~= reference_is_right)); %#ok<AGROW> % onto the left-stim convention
                d_side{trial_is_right(i)+1}(end+1) = d;
            end
        end
        fd(f) = mean(d_all); ns(f) = numel(d_all);
        for sd = 1:2; if ~isempty(d_side{sd}); fds(f,sd) = mean(d_side{sd}); end; end
    end
    fly_disp{gi} = fd; fly_disp_side{gi} = fds; n_disp_stims{gi} = ns;
end
fprintf('net displacement %g-%g s (rho < %.2f masked; %d stims dropped for a bar sweep):\n', disp_t0, disp_t1, disp_rho_min, n_sweep_drop);
for gi = 1:2
    v = fly_disp{gi}; s2 = fly_disp_side{gi}; b2 = all(~isnan(s2),2);
    fprintf('  %-11s %+.2f +/- %.2f rad (n=%d flies, %d-%d stims/fly; vs 0 p=%.3g) | left-side %+.2f, right-side %+.2f rad (paired p=%.3g, n=%d)\n', ...
        geno_label{gi}, mean(v), std(v)/sqrt(numel(v)), numel(v), min(n_disp_stims{gi}), max(n_disp_stims{gi}), signrank(v), ...
        mean(s2(:,1),'omitnan'), mean(s2(:,2),'omitnan'), signrank(s2(b2,1),s2(b2,2)), sum(b2));
end
fprintf('  LPsP vs Empty (ranksum): p=%.3g\n', ranksum(fly_disp{1}, fly_disp{2}));

%% example mean images: raw data when available, else the V1 PNG crop
% trial_mean_rgb returns each image with glomerulus 1 (lowest PB angle) at
% the BOTTOM, as in V1 Figure 3; here PB angle runs left -> right, so it's
% rotated 90 deg clockwise (bottom -> left).
gcamp_color = [0 .2 1];   % EPG calcium blue (V1 convention)
atp_color   = [.8 0 0];   % atp channel red
mean_img_planes = 3:10;   % as V1: planes 1-2 hold bright cell bodies above the PB

% fallback only: pixel boxes [x0 x1 y0 y1] of the two mean images in
% F3_r12.png (largest non-white region beside each heatmap), bottom ~30 px
% trimmed where V1 drew its EPG/ATP labels over the image
mean_png       = fullfile(repo_root,'ugly_figures','all_figs','panels','Fig3','F3_r12.png');
mean_png_box   = [30 255 441 1031; 2705 2943 363 1109]; % rows: LPsP, Empty
mean_png_trial = {'20260428-15_epg_syt8s_blind_p2x2', '20260429-11_epg_syt8s_blind_p2x2'};
mean_png_stim  = [2 4];

mean_rgb = cell(1,2);
for c = 1:2
    trial_name = meta_display(all_data(example_idx(c)).meta);
    trial_dir  = fullfile(raw_dir, trial_name);
    if isfile(fullfile(trial_dir,'imagingData_trial001.mat')) && isfile(fullfile(trial_dir,'mask.mat'))
        fprintf('mean image %d: raw data, %s\n', c, trial_dir);
        mean_rgb{c} = rot90(trial_mean_rgb(trial_dir, gcamp_color, atp_color, mean_img_planes), -1);
    else
        warning('Figure1:meanImageFallback', ['mean image %d: no raw files in %s -- using the low-res crop ' ...
            'of V1 Figure 3''s PNG. Don''t ship this version.'], c, trial_dir);
        assert(strcmp(trial_name, mean_png_trial{c}) && example_stim_idx(c)==mean_png_stim(c), 'Figure1:meanImage', ...
            'example %d changed and has no raw files -- the PNG crop doesn''t match it either', c);
        I = im2double(imread(mean_png)); b = mean_png_box(c,:);
        mean_rgb{c} = rot90(I(b(3):b(4), b(1):b(2), :), -1); clear I
    end
end

%% plot
geno_color  = [.5 .5 .5; gcamp_color]; % empty = gray, LPsP = EPG blue
zoom_pre_s  = 2; zoom_post_s = 20;     % example window around the stim, s

fig = figure('color','w','Visible','off','Position',[50 50 1500 1900]);
% 7 x 3 grid: row 1 -> grid rows 1-2, row 2 (vertical examples) -> 3-5, row 3 -> 6-7
T = tiledlayout(fig,7,3,'TileSpacing','loose','Padding','compact');
tile = @(r,c) (r-1)*3 + c;

% ---- row 1
placeholder(nexttile(T,tile(1,1),[2 1]), 'compass cartoon')
placeholder(nexttile(T,tile(1,2),[2 1]), 'compass inputs')
axm = nexttile(T,tile(1,3),[2 1]);
if isfile(morph_png)
    image(axm, imread(morph_png)); axis(axm,'image','off')
else
    placeholder(axm, 'LPsP morphology: run fig_lpsp_morphology_flywire.py')
end

% ---- row 2
placeholder(nexttile(T,tile(3,1),[3 1]), 'dangling fly schematic')
% The two example panels are placed by hand inside their grid slot (not as
% nested tiledlayouts): the mean image has to be exactly as wide as the
% heatmap below it, with its own aspect ratio, and MATLAB doesn't report
% reliable positions for axes inside nested layouts. The slot itself comes
% from a throwaway top-level tile (top-level tile positions ARE in figure
% units, since T fills the figure).
fig_px = fig.Position(3:4);
heat_axes = gobjects(1,2);
for c = 1:2
    i = example_idx(c); ex = all_data(i);
    xb = ex.ft.xb;
    stim_t = trial_events{i}(example_stim_idx(c)).onset_t;
    stims_ttl = logical(ex.ft.stims(:));                       % ATP TTL on the behavior clock
    ttl_end = ex.ft.xf(find(ex.ft.xf(:) > stim_t & ~stims_ttl, 1)); % first sample after the TTL goes low
    atp_dur = ttl_end - stim_t;                                % 2.02 s for every stim in this dataset
    zi = xb>=stim_t-zoom_pre_s & xb<=stim_t+zoom_post_s;
    tr = xb(zi) - stim_t;                        % time re: stim onset
    alpha = unwrap(ex.im.alpha);
    z = ex.im.z(:,zi);
    clim_c = [prctile(min(z,[],1),5), prctile(max(z,[],1),95)];

    slot_ax = nexttile(T,tile(3,1+c),[3 1]); drawnow
    b = slot_ax.Position; delete(slot_ax)         % [x y w h], figure-normalized
    a  = size(mean_rgb{c},1) / size(mean_rgb{c},2);
    xh = b(1) + .22*b(3); wh = .70*b(3);           % heatmap (and mean image) x-span
    hi = wh*fig_px(1)*a / fig_px(2);               % mean-image height that keeps its aspect
    yb = b(2) + .10*b(4);                          % heatmap bottom (x labels + colorbar below)
    yt = b(2) + .97*b(4) - hi - .004;              % heatmap top (image + title above)
    pos_heat  = [xh, yb, wh, yt-yb];
    pos_img   = [xh, yt+.004, wh, hi];
    pos_cb    = [xh, yb - .03*b(4), wh, .012*b(4)];      % directly under the heatmap (it has no x axis)

    % EPG calcium heatmap: PB angle across, time down
    axh = axes(fig,'Position',pos_heat); hold(axh,'on')
    imagesc(axh, alpha, tr, z', clim_c)
    colormap(axh, white_to_color(gcamp_color))
    plot(axh, [min(alpha) max(alpha)], [0 0], ':k', 'LineWidth', 1)
    set(axh,'YDir','reverse','Box','off'); xlim(axh,[min(alpha) max(alpha)]); ylim(axh,[-zoom_pre_s zoom_post_s])
    cb = colorbar(axh,'southoutside'); cb.Label.String = 'EPG z-score';
    axh.Position = pos_heat; cb.Position = pos_cb;   % colorbar() shrinks the axes; put both back
    minmax_colorbar(cb)
    add_time_scalebar_vertical(axh, 5)
    % ATP: red bar spanning the TTL (0 to atp_dur s), in the same column as the time scale bar
    xl_h = xlim(axh); x_bar = xl_h(2) + .04*diff(xl_h);
    plot(axh, [x_bar x_bar], [0 atp_dur], '-', 'Color', atp_color, 'LineWidth', 5, 'Clipping','off')
    text(axh, x_bar + .03*diff(xl_h), atp_dur/2, 'ATP', 'Color', atp_color, 'FontWeight','bold', 'FontSize', 10, ...
        'HorizontalAlignment','left', 'VerticalAlignment','middle', 'Clipping','off')
    % "Time" arrow down the left side of the heatmap (time runs top -> bottom)
    % (short: spans the top ~third of the heatmap, starting at its top edge)
    x_arr = xl_h(1) - .06*diff(xl_h); yl_h = ylim(axh); y_arr = yl_h(1) + [0 .3]*diff(yl_h);
    plot(axh, [x_arr x_arr], y_arr, 'k-', 'LineWidth', 1.5, 'Clipping','off')
    plot(axh, x_arr, y_arr(2), 'kv', 'MarkerFaceColor','k', 'MarkerSize', 8, 'Clipping','off')
    text(axh, x_arr - .03*diff(xl_h), mean(y_arr), 'Time', 'Rotation', 90, 'HorizontalAlignment','center', ...
        'VerticalAlignment','bottom', 'FontSize', 11, 'Clipping','off')
    axh.XColor = 'none'; % no PB-angle ticks either: the mean image above shows the PB
    heat_axes(c) = axh;

    % mean image: PB angle left -> right, same width as the heatmap
    axm = axes(fig,'Position',pos_img);
    image(axm, mean_rgb{c}); axis(axm,'off')
    text(axm, 1, 1, 'EPG ', 'Color', gcamp_color, 'FontWeight','bold', 'FontSize', 10, ...
        'HorizontalAlignment','right', 'VerticalAlignment','top')
    text(axm, 1, size(mean_rgb{c},1), 'ATP ', 'Color', atp_color, 'FontWeight','bold', 'FontSize', 10, ...
        'HorizontalAlignment','right', 'VerticalAlignment','bottom')
    title(axm, sprintf('%s (%s-side stim)', example_geno{c}, ref_side_str), 'FontWeight','bold', 'Visible','on')
end

% ---- row 3
ax3 = nexttile(T,tile(6,1),[2 1]); hold(ax3,'on')
h = gobjects(1,2);
for gi = 1:2
    patch(ax3,[t_trace,fliplr(t_trace)],[geno_mean(gi,:)+geno_sem(gi,:),fliplr(geno_mean(gi,:)-geno_sem(gi,:))], ...
        geno_color(gi,:),'FaceAlpha',.15,'EdgeColor','none')
    h(gi) = plot(ax3,t_trace,geno_mean(gi,:),'Color',geno_color(gi,:),'LineWidth',2);
end
xline(ax3,0,':k'); yline(ax3,0,':k');
xlim(ax3,[trace_pre_s trace_post_s])
ylabel(ax3,'bump position (rad, unwrapped)')
legend(ax3,h,arrayfun(@(g) sprintf('%s (n=%d flies)',geno_label{g},geno_n(g)),1:2,'UniformOutput',false), ...
    'Location','northwest','Box','off')
title(ax3,'bump position aligned to ATP onset')
set_pi_yticks(ax3)
add_time_scalebar(ax3, 5)

% row 3, panel 2: net displacement 4-15 s, left-stim convention, one dot per fly
ax4 = nexttile(T,tile(6,2),[2 1]); hold(ax4,'on')
rng(1) % reproducible jitter
for gi = 1:2
    v = fly_disp{gi};
    scatter(ax4, gi + .12*(rand(numel(v),1)-.5), v, 36, geno_color(gi,:), 'filled', 'MarkerFaceAlpha', .6)
    errorbar(ax4, gi+.3, mean(v), std(v)/sqrt(numel(v)), 'o', 'Color','k', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 0)
end
yline(ax4, 0, ':k');
xlim(ax4,[.5 2.7]); xticks(ax4,[1 2]); xticklabels(ax4, geno_label)
ylabel(ax4, sprintf('net bump displacement per stim (rad, %s-stim sign)', ref_side_str))
add_pbar(ax4, 1, 2, ranksum(fly_disp{1}, fly_disp{2}))
title(ax4, sprintf('%g-%g s after ATP onset', disp_t0, disp_t1))

% row 3, panel 3: same, raw sign, split by stimulated side; lines join the same fly
ax5 = nexttile(T,tile(6,3),[2 1]); hold(ax5,'on')
xpos = [1 2; 3.4 4.4]; % [left right] per genotype
for gi = 1:2
    s2 = fly_disp_side{gi}; b2 = all(~isnan(s2),2);
    plot(ax5, xpos(gi,:), s2(b2,:)', '-', 'Color', [geno_color(gi,:) .35])
    for sd = 1:2
        v = s2(~isnan(s2(:,sd)),sd);
        scatter(ax5, xpos(gi,sd)*ones(size(v)), v, 30, geno_color(gi,:), 'filled', 'MarkerFaceAlpha', .6)
    end
    errorbar(ax5, xpos(gi,:) + [-.18 .18], mean(s2,1,'omitnan'), std(s2,0,1,'omitnan')./sqrt(sum(~isnan(s2),1)), ...
        'ok-', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 0)
    add_pbar(ax5, xpos(gi,1), xpos(gi,2), signrank(s2(b2,1), s2(b2,2)))
end
yline(ax5, 0, ':k');
xlim(ax5,[.4 5]); xticks(ax5, reshape(xpos',1,[])); xticklabels(ax5, {'left','right','left','right'})
for gi = 1:2 % genotype names under each pair of side ticks
    text(ax5, (mean(xpos(gi,:))-.4)/4.6, -.1, geno_label{gi}, 'Units','normalized', 'HorizontalAlignment','center', ...
        'VerticalAlignment','top', 'FontSize', 11, 'Color', geno_color(gi,:))
end
ylabel(ax5, 'net bump displacement per stim (rad, raw sign)')
title(ax5, 'by stimulated side (lines join the same fly)')

% text sizes (titles keep their own)
for ax_k = findall(fig,'Type','axes')'
    if strcmp(ax_k.Visible,'off'); continue; end
    fs = ax_k.Title.FontSize; ax_k.FontSize = 11; ax_k.Title.FontSize = fs;
end

%% export
if ~isfolder(export_dir); mkdir(export_dir); end
if ~isfolder(all_figs_dir); mkdir(all_figs_dir); end
exportgraphics(fig, fullfile(export_dir,'Figure_1_V2.png'), 'Resolution', 200)
exportgraphics(fig, fullfile(all_figs_dir,'Fig1_V2.pdf'), 'ContentType', 'vector')
if isfolder(lpsp_dropbox) % shared copy in Dropbox (PDF for Illustrator, PNG for a quick look)
    if ~isfolder(dropbox_dir); mkdir(dropbox_dir); end
    copyfile(fullfile(all_figs_dir,'Fig1_V2.pdf'), dropbox_dir)
    copyfile(fullfile(export_dir,'Figure_1_V2.png'), fullfile(dropbox_dir,'Fig1_V2.png'))
    fprintf('copied Fig1_V2.pdf/.png to %s\n', dropbox_dir);
end
fig.Visible = 'on';

%% ===================== functions =====================
function m = masked_unwrapped_mu(mu, rho, t, rho_min)
    % bump position with weak-bump frames removed: NaN where rho < rho_min,
    % unwrapped across the remaining frames only, linearly interpolated
    % through the gaps, nearest solid value at the ends (as
    % fig_p2x2_displacement_per_stim.m)
    ok = rho >= rho_min & ~isnan(mu);
    m = nan(size(mu));
    if sum(ok) < 2; return; end
    m(ok) = unwrap(mu(ok));
    m = fillmissing(interp1(t(ok), m(ok), t, 'linear'), 'nearest');
end

function minmax_colorbar(cb)
    % label only the colorbar's two ends (figure-set convention)
    lim = cb.Limits;
    cb.Ticks = lim;
    cb.TickLabels = arrayfun(@(v) sprintf('%.1f', v), lim, 'UniformOutput', false);
end

function add_pbar(ax, x1, x2, p)
    % bracket between x1 and x2 above the data with the p-value; raises ylim to fit
    yl = ylim(ax); yb = yl(2) + .04*diff(yl);
    plot(ax, [x1 x1 x2 x2], yb - [.02 0 0 .02]*diff(yl), 'k', 'HandleVisibility','off')
    text(ax, mean([x1 x2]), yb, sprintf('p=%.2g',p), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize', 9)
    ylim(ax, [yl(1) yb + .1*diff(yl)])
end

function placeholder(ax, label)
    % dashed empty box with a label, for panels that will be drawn in Illustrator
    hold(ax,'on')
    rectangle(ax,'Position',[0 0 1 1],'LineStyle','--','EdgeColor',[.6 .6 .6])
    text(ax,.5,.5,label,'HorizontalAlignment','center','Color',[.45 .45 .45],'FontSize',12,'FontAngle','italic')
    xlim(ax,[0 1]); ylim(ax,[0 1]); axis(ax,'off')
end

function add_time_scalebar_vertical(ax, len_s)
    % vertical len_s time scale bar just right of the axes' bottom edge
    % (time on y). Call after the limits are final (freezes them).
    xl = xlim(ax); yl = ylim(ax); xlim(ax,xl); ylim(ax,yl);
    ax.YColor = 'none';
    x = xl(2) + .04*diff(xl);
    if strcmp(ax.YDir,'reverse'); y = yl(2) - [len_s 0]; else; y = yl(1) + [0 len_s]; end
    plot(ax, [x x], y, 'k-', 'LineWidth', 2, 'Clipping','off', 'HandleVisibility','off')
    text(ax, x + .03*diff(xl), mean(y), sprintf('%g s',len_s), 'HorizontalAlignment','left', ...
        'VerticalAlignment','middle', 'FontSize', 9, 'Clipping','off')
end

function add_time_scalebar(ax, len_s)
    % replace a time axis with a len_s scale bar under the axes' bottom-right
    % corner (ported from Figure_3_claude.m). Call after limits are final.
    xl = xlim(ax); yl = ylim(ax); xlim(ax,xl); ylim(ax,yl);
    ax.XColor = 'none';
    if strcmp(ax.YDir,'reverse'); y = yl(2) + .06*diff(yl); dy = .04*diff(yl); else; y = yl(1) - .06*diff(yl); dy = -.04*diff(yl); end
    x = xl(2) - [len_s 0];
    plot(ax, x, [y y], 'k-', 'LineWidth', 2, 'Clipping', 'off', 'HandleVisibility', 'off')
    text(ax, mean(x), y + dy, sprintf('%g s', len_s), 'HorizontalAlignment','center', 'VerticalAlignment', 'middle', 'FontSize', 9, 'Clipping', 'off')
end

function ev = get_trial_stim_events(trial, t_grid, peri_win, base_win, peak_win, peak_factor)
    % one entry per CONFIRMED-ejection stim -- ported unchanged from
    % Figure_3_claude.m (see its header for every field). Fields used here:
    % onset_t, right_frac, mu_rel, head_rel, rho_grid, atp_peak_amp.
    xb = trial.ft.xb; xf = trial.ft.xf;
    stims  = logical(trial.ft.stims(:));
    onsets = find(diff([false;stims])==1);
    atp_sig = sum(trial.atp.d,1);
    t_onset = xf(onsets);
    mu_unwrapped = unwrap(trial.im.mu);
    dt = median(diff(xf));
    heading = cumsum(trial.ft.r_speed)*dt;
    atp_dff_avg = mean(trial.atp.d,1);
    t_common_local = linspace(-peri_win,peri_win,101);
    ev = struct('onset_t',{},'right_frac',{},'mu_rel',{},'head_rel',{},'rho_grid',{},'atp_dff_rel',{},'atp_peak_amp',{});
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s);
        xf_rel = xf - t_onset(s);
        peri = interp1(xb_rel,atp_sig,t_common_local,'linear',nan);
        base_idx = t_common_local >= base_win(1) & t_common_local < base_win(2);
        peak_idx = t_common_local > peak_win(1) & t_common_local <= peak_win(2);
        if (median(peri(peak_idx),'omitnan') - median(peri(base_idx),'omitnan')) <= peak_factor*std(peri(base_idx),'omitnan')
            continue % not a confirmed ejection
        end
        peak_frames = xb_rel>peak_win(1) & xb_rel<=peak_win(2);
        wedge_peak = mean(trial.atp.d(:,peak_frames),2);
        mu_peri   = interp1(xb_rel,mu_unwrapped,t_grid,'linear',nan);
        head_peri = interp1(xf_rel,heading,t_grid,'linear',nan);
        k = numel(ev)+1;
        ev(k).onset_t     = t_onset(s);
        ev(k).right_frac  = mean(wedge_peak(1:16)) / (mean(wedge_peak(1:16)) + mean(wedge_peak(17:32)));
        ev(k).mu_rel      = mu_peri   - interp1(t_grid,mu_peri,0);
        ev(k).head_rel    = head_peri - interp1(t_grid,head_peri,0);
        ev(k).rho_grid    = interp1(xb_rel,trial.im.rho,t_grid,'linear',nan);
        ev(k).atp_dff_rel = interp1(xb_rel,atp_dff_avg,t_grid,'linear',nan);
        ev(k).atp_peak_amp = max(ev(k).atp_dff_rel(t_grid>=0 & t_grid<=5)) - mean(ev(k).atp_dff_rel(t_grid>=-2 & t_grid<0),'omitnan');
    end
end

function name = trial_fly_folder_name(meta_path)
    % "fly N _ <genotype>" folder above the trial folder; splits on both
    % '\' and '/' since .meta paths were saved on Windows
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part); name = ''; else; name = parts{fly_part(1)}; end
end

function fid = trial_fly_id(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part); fid = meta_path; else; fid = strjoin(parts(1:fly_part(1)),'/'); end
end

function s = meta_display(meta_path)
    % trial folder name (skipping a trailing generic "registration" folder)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if strcmpi(parts{end},'registration') && numel(parts)>1; s = parts{end-1}; else; s = parts{end}; end
end

function rough = atp_trace_roughness(ex, onset_t)
    % SD of one stim's whole-PB atp dF/F (-2..20 s) around a 1 s-gaussian
    % smoothed copy, over the smoothed peak above baseline
    xb = ex.ft.xb;
    idx = xb>=onset_t-2 & xb<=onset_t+20;
    tr = mean(ex.atp.d(:,idx),1);
    dt = median(diff(xb));
    sm = smoothdata(tr,'gaussian',round(1/dt));
    rough = std(tr-sm) / (max(sm) - mean(sm(1:round(2/dt))));
end

function set_pi_yticks(ax)
    % y ticks at multiples of pi/2, labeled in pi
    yl = ylim(ax); ylim(ax,yl);
    q = 1; while diff(yl)/(q*pi/2) > 8; q = q*2; end
    k = q*(ceil(yl(1)/(q*pi/2)):floor(yl(2)/(q*pi/2)));
    yticks(ax,k*pi/2);
    yticklabels(ax,arrayfun(@pi_half_label,k,'UniformOutput',false));
end

function s = pi_half_label(k)
    if k == 0; s = '0'; return; end
    if mod(k,2) == 0; n = k/2; d = ''; else; n = k; d = '/2'; end
    if n == 1; s = ['\pi' d]; elseif n == -1; s = ['-\pi' d]; else; s = sprintf('%d\\pi%s',n,d); end
end

function len = nice_scale_len(x)
    % largest 1/2/5 x 10^k value <= x
    p = 10^floor(log10(x));
    m = [1 2 5 10];
    len = p*m(find(m*p <= x,1,'last'));
end

function cmap = white_to_color(max_color)
    n = 256;
    cmap = [linspace(1,max_color(1),n)', linspace(1,max_color(2),n)', linspace(1,max_color(3),n)'];
end

function rgb = trial_mean_rgb(trial_dir, gcamp_color, atp_color, planes)
    % whole-trial mean image of one trial, two-channel RGB on white -- ported
    % from Figure_3_claude.m (see there for the full rationale). Differences:
    % takes the trial FOLDER holding both imagingData_trial001.mat and
    % mask.mat (the Dropbox copies are flat; on Z: the movie sits one level
    % down in registration/). img{1} = gcamp, img{2} = atp, per-plane
    % rows x cols x planes x frames; only `planes` are summed. Cropped to the
    % mask bounding box (+5 px), rotated so glomerulus 1 is at the BOTTOM
    % (and mirrored, as V1). Each channel clipped to its own 1st-99.5th
    % (gcamp: 95th) percentile over PB-mask pixels.
    S = load(fullfile(trial_dir,'imagingData_trial001.mat'),'img');
    g = sum(mean(S.img{1}(:,:,planes,:),4),3);
    a = sum(mean(S.img{2}(:,:,planes,:),4),3);
    clear S
    mk = load(fullfile(trial_dir,'mask.mat'),'mask');
    [ym,xm] = find(mk.mask);
    pad = 5;
    rows = max(1,min(ym)-pad):min(size(g,1),max(ym)+pad);
    cols = max(1,min(xm)-pad):min(size(g,2),max(xm)+pad);
    g = fliplr(rot90(g(rows,cols))); a = fliplr(rot90(a(rows,cols)));
    m = fliplr(rot90(mk.mask(rows,cols)));
    norm01 = @(x,ref,hi) min(max((x-prctile(ref,1))/(prctile(ref,hi)-prctile(ref,1)),0),1);
    rgb = max(1 - norm01(g,g(m),95).*reshape(1-gcamp_color,1,1,3) - norm01(a,a(m),99.5).*reshape(1-atp_color,1,1,3), 0);
end
