%% fig_p2x2_displacement_per_stim
% Scratch (not in Figure 1 yet): how far the bump rotates per ATP stim,
% LPsP>P2X2 vs Empty>P2X2, ONE DOT PER FLY, plotted in RADIANS (signed net
% displacement; the command window also gives rotations = rad / 2*pi). Dataset and stim/side handling are Figure_1_V2.m's
% (lpsp_p2x2_reredo_20260518.mat; every trial has 4 TTL stims at
% 41/101/161/221 s, 60 s apart, in a 300 s trial; ~16 stims per fly).
%
% Per stim: net displacement = unwrapped bump position (im.mu) averaged over
% the avg_s before (onset + end_s) minus the same average over the avg_s
% before onset, over 0-15 s after each stim (windows_s; chosen 2026-11-23).
% Longer windows were tried: over the whole 0-58 s inter-stim interval the
% bump's own drift (~10 rad unsigned per minute in both genotypes -- not the
% fly turning: the ball barely moves in these dangling flies) swamps the
% stim effect; 0-20 s (Figure 1's average-trace window) behaves like 0-15 s.
% Signed, not absolute: flipped onto one side convention (positive = the
% direction a left-side stim pushes the bump, as Figure 1's average trace),
% using each trial's side call from its most confident ejection. Only stims
% with a confirmed ejection (per-stim atp detector, as Figure 1) -- in this
% dataset that is essentially every stim. Per fly = mean over its stims;
% group = mean +/- SEM across flies. Tests: signrank vs 0 per group,
% ranksum between groups.
%
% STIM NUMBER (stim_nums): which of each trial's 4 stims to include, e.g.
% [3 4] for only the later stims (stim_nums_override in the workspace).
%
% BAR STILL ONLY: the visual pattern (bright bar, or the background pattern
% in the other trial type) is swept in open loop on a fixed schedule: it
% holds still from ~0.5 s before each stim until 18.6 s after it, then sweeps
% (~47 rad of cue motion) until just before the next stim. Only the still
% period is scored -- a window is set to NaN for any stim where it would run
% into a sweep (none at 15 s; the count is printed). Sweeps are found from
% ft.cue: |cue velocity| > 0.05 rad/s (0.5 s smoothing), bouts < 1 s dropped
% (the cue channel blips briefly at each stim).
%
% FULL ROTATIONS (command window only; was a panel, dropped 2026-11-23 as it
% didn't separate the groups): % of each fly's stims whose net rotation over
% the window exceeds one full turn (|net| > 2*pi), and the share > +1 turn
% toward the stim side.
%
% LOW-RHO MASKING (rho_min): when the bump is weak (PVA strength im.rho low)
% its phase is mostly noise, and unwrapping turns noisy jumps past +-pi into
% spurious rotation -- inflating apparent wiggle, most of all over long
% windows. So bump position is set to NaN on frames with rho < rho_min (0.1:
% removes the near-zero-rho frames, ~6% of frames, while keeping almost all
% data; a 0-0.3 sweep moved the LPsP effect only within +0.2-0.35 rotations
% but the LPsP-vs-Empty p between 0.01 and 0.4, so the threshold was fixed
% before looking further), unwrapped across the remaining frames only, and linearly
% interpolated through the gaps (each gap is crossed by the shortest path
% between the solid frames on either side; gaps at a trial's ends take the
% nearest solid value). The plot uses the masked position; the command
% window also prints the unmasked numbers and how much of each window was
% masked, per genotype.

%% paths + data
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir = fullfile(repo_root,'.data'); if ~isfolder(data_dir); data_dir = fullfile(repo_root,'data'); end
export_dir = fullfile(repo_root,'ugly_figures','exports','V2','scratch');
if ~exist('all_data','var') || numel(all_data) ~= 101
    tmp = load(fullfile(data_dir,'lpsp_p2x2_reredo_20260518.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
n_trials = numel(all_data);

%% parameters
windows_s = 15;       % score from onset to onset + this (s; next stim at 60 s). A vector plots one panel per window
avg_s = 1;            % average mu over this long at each end
rho_min = 0.1;        % bump position NaN'd below this PVA strength, then interpolated
if exist('rho_min_override','var'); rho_min = rho_min_override; end % set in the workspace to try another threshold
stim_nums = 1:4;      % which stims within each trial (TTL order; every trial has 4) to include
if exist('stim_nums_override','var'); stim_nums = stim_nums_override; end
peri_win = 5; base_win = [-5 0]; peak_win = [0 3]; peak_factor = 3; % ejection detector (as Figure 1)
nw = numel(windows_s);

%% per trial: genotype, fly, side call, per-stim displacement for each window
n_sweep_excluded = zeros(1,nw); % stims whose window would run into a bar sweep
genotype = cell(n_trials,1); fly_id = cell(n_trials,1);
trial_is_right = nan(n_trials,1);
stim_disp = cell(n_trials,1); stim_conf = cell(n_trials,1); % stim_disp{i}: n_stims x nw (rho-masked)
stim_disp_raw = cell(n_trials,1); stim_lowfrac = cell(n_trials,1); % unmasked; fraction of 0-max(window) frames masked
for i = 1:n_trials
    tr = all_data(i);
    fname = fly_folder(tr.meta);
    if contains(fname,'lpsp','IgnoreCase',true); genotype{i} = 'lpsp';
    elseif contains(fname,'empty','IgnoreCase',true); genotype{i} = 'empty';
    else; genotype{i} = ''; end
    fly_id{i} = fly_path(tr.meta);

    xb = tr.ft.xb; xf = tr.ft.xf;
    onsets = xf(find(diff([false; logical(tr.ft.stims(:))])==1));
    cue = unwrap(tr.ft.cue(:)); dtf = median(diff(xf));
    cue_v = [0; diff(cue)] ./ [1; diff(xf(:))];
    sweeping = bwareaopen(abs(smoothdata(cue_v,'movmean',round(.5/dtf))) > 0.05, round(1/dtf));
    mu_raw = unwrap(tr.im.mu(:));
    mu = masked_unwrapped_mu(tr.im.mu(:), tr.im.rho(:), xb(:), rho_min);
    low = tr.im.rho(:) < rho_min;
    atp_sig = sum(tr.atp.d,1);
    tc = linspace(-peri_win,peri_win,101);
    d = nan(numel(onsets),nw); d_raw = d; lowfrac = nan(numel(onsets),1);
    conf = false(numel(onsets),1); rfs = nan(numel(onsets),1);
    for s = 1:numel(onsets)
        t0 = onsets(s);
        pre = xb>=t0-avg_s & xb<t0;
        for w = 1:nw
            post = xb>=t0+windows_s(w)-avg_s & xb<t0+windows_s(w);
            next_sweep = xf(find(xf(:) > t0 & sweeping, 1)); % bar starts moving again
            bar_still = isempty(next_sweep) || t0 + windows_s(w) <= next_sweep;
            n_sweep_excluded(w) = n_sweep_excluded(w) + ~bar_still;
            if any(pre) && any(post) && bar_still
                d(s,w)     = mean(mu(post))     - mean(mu(pre));
                d_raw(s,w) = mean(mu_raw(post)) - mean(mu_raw(pre));
            end
        end
        lowfrac(s) = mean(low(xb>=t0 & xb<t0+max(windows_s)));
        peri = interp1(xb-t0, atp_sig, tc, 'linear', nan);
        bi = tc>=base_win(1) & tc<base_win(2); pk = tc>peak_win(1) & tc<=peak_win(2);
        conf(s) = (median(peri(pk),'omitnan') - median(peri(bi),'omitnan')) > peak_factor*std(peri(bi),'omitnan');
        if conf(s)
            wp = mean(tr.atp.d(:, xb-t0>peak_win(1) & xb-t0<=peak_win(2)),2);
            rfs(s) = mean(wp(1:16)) / (mean(wp(1:16)) + mean(wp(17:32)));
        end
    end
    if any(conf)
        [~,k] = max(abs(rfs-.5)); trial_is_right(i) = rfs(k) > .5; % most confident ejection decides the side
    end
    stim_disp{i} = d; stim_conf{i} = conf; stim_disp_raw{i} = d_raw; stim_lowfrac{i} = lowfrac;
end

%% pool per fly, in rotations (sign-flipped to the LEFT-stim convention, as Figure 1)
ref_is_right = false;
geno_order = {'empty','lpsp'}; geno_label = {'Empty>P2X2','LPsP>P2X2'};
fly_rot = cell(2,nw); fly_rot_raw = cell(2,nw); fly_n = cell(2,1); fly_low = cell(2,1);
fly_full = cell(2,nw); fly_full_dir = cell(2,nw); % % of stims with |net| > 1 turn; > +1 turn toward the stim side
for g = 1:2
    rows = find(strcmp(genotype,geno_order{g}) & ~isnan(trial_is_right));
    flies = unique(fly_id(rows));
    fly_n{g} = zeros(numel(flies),1);
    fly_low{g} = nan(numel(flies),1);
    for w = 1:nw
        [fly_rot{g,w}, fly_rot_raw{g,w}, fly_full{g,w}, fly_full_dir{g,w}] = deal(nan(numel(flies),1));
    end
    for f = 1:numel(flies)
        dd = []; dr = []; lf = [];
        for i = rows(strcmp(fly_id(rows),flies{f}))'
            sgn = 1 - 2*(trial_is_right(i) ~= ref_is_right);
            use = stim_conf{i} & ismember((1:numel(stim_conf{i}))', stim_nums); % confirmed, and a requested stim number
            dd = [dd; sgn*stim_disp{i}(use,:)];     %#ok<AGROW>
            dr = [dr; sgn*stim_disp_raw{i}(use,:)]; %#ok<AGROW>
            lf = [lf; stim_lowfrac{i}(use)];        %#ok<AGROW>
        end
        fly_n{g}(f) = size(dd,1); fly_low{g}(f) = mean(lf,'omitnan');
        for w = 1:nw
            fly_rot{g,w}(f)     = mean(dd(:,w),'omitnan') / (2*pi);
            fly_rot_raw{g,w}(f) = mean(dr(:,w),'omitnan') / (2*pi);
            rot = dd(~isnan(dd(:,w)),w) / (2*pi);
            fly_full{g,w}(f)     = 100*mean(abs(rot) > 1);
            fly_full_dir{g,w}(f) = 100*mean(rot > 1);
        end
    end
end
fprintf('\nstims excluded because the window ran into a bar sweep: %s (windows %s s)\n', mat2str(n_sweep_excluded), mat2str(windows_s));
fprintf('frames with rho < %.2f in the 0-%d s windows (masked): Empty %.1f%% (fly range %.1f-%.1f%%), LPsP %.1f%% (%.1f-%.1f%%)\n', ...
    rho_min, max(windows_s), 100*mean(fly_low{1}), 100*min(fly_low{1}), 100*max(fly_low{1}), ...
    100*mean(fly_low{2}), 100*min(fly_low{2}), 100*max(fly_low{2}));
for w = 1:nw
    fprintf('\n[0-%d s after each stim, rotations per stim; rho < %.2f masked + interpolated]\n', windows_s(w), rho_min);
    for g = 1:2
        v = fly_rot{g,w};
        fprintf('  %-11s n=%2d flies (stims/fly %d-%d): %+.3f +/- %.3f rot (%+.2f rad), signrank vs 0 p=%.3g\n', ...
            geno_label{g}, numel(v), min(fly_n{g}), max(fly_n{g}), mean(v), std(v)/sqrt(numel(v)), 2*pi*mean(v), signrank(v));
    end
    fprintf('  ranksum LPsP vs Empty: p=%.3g\n', ranksum(fly_rot{1,w}, fly_rot{2,w}));
    for g = 1:2
        fprintf('  %-11s stims with >1 full turn: %.1f +/- %.1f%% of stims per fly (%d/%d flies have any); toward stim side: %.1f%%\n', ...
            geno_label{g}, mean(fly_full{g,w}), std(fly_full{g,w})/sqrt(numel(fly_full{g,w})), sum(fly_full{g,w}>0), ...
            numel(fly_full{g,w}), mean(fly_full_dir{g,w}));
    end
    fprintf('  ranksum >1 turn %%, LPsP vs Empty: p=%.3g\n', ranksum(fly_full{1,w}, fly_full{2,w}));
    fprintf('  (unmasked: Empty %+.3f +/- %.3f, LPsP %+.3f +/- %.3f rot, ranksum p=%.3g)\n', ...
        mean(fly_rot_raw{1,w}), std(fly_rot_raw{1,w})/sqrt(numel(fly_rot_raw{1,w})), ...
        mean(fly_rot_raw{2,w}), std(fly_rot_raw{2,w})/sqrt(numel(fly_rot_raw{2,w})), ranksum(fly_rot_raw{1,w}, fly_rot_raw{2,w}));
end

%% plot: one panel per window, per-fly mean net displacement in radians
gcol = [.5 .5 .5; 0 .2 1];
fig = figure('Color','w','Position',[100 100 120+380*nw 480]);
tl = tiledlayout(fig,1,nw,'TileSpacing','loose');
rng(1) % reproducible jitter
for w = 1:nw
    ax = nexttile(tl); hold(ax,'on')
    for g = 1:2
        v = 2*pi*fly_rot{g,w}; % rotations -> rad
        scatter(ax, g + .12*(rand(numel(v),1)-.5), v, 36, gcol(g,:), 'filled', 'MarkerFaceAlpha', .6)
        errorbar(ax, g+.3, mean(v), std(v)/sqrt(numel(v)), 'o', 'Color','k', 'MarkerFaceColor','k', 'LineWidth', 1.5, 'CapSize', 0)
    end
    yline(ax, 0, ':k');
    xlim(ax,[.5 2.7]); xticks(ax,[1 2]); xticklabels(ax, geno_label)
    ylabel(ax, 'bump displacement per stim (rad, left-stim sign)')
    p = ranksum(fly_rot{1,w}, fly_rot{2,w});
    yl = ylim(ax); yb = yl(2) + .05*diff(yl);
    plot(ax,[1 1 2 2], yb - [.02 0 0 .02]*diff(yl), 'k'); text(ax,1.5,yb,sprintf('p=%.2g',p),'HorizontalAlignment','center','VerticalAlignment','bottom')
    ylim(ax,[yl(1) yb+.1*diff(yl)])
    title(ax, sprintf('0-%d s after each stim', windows_s(w)))
end
title(tl, {'net bump displacement per ATP stim, one dot per fly', ...
    sprintf('n=%d Empty, %d LPsP flies; stims %s; \\rho < %.2f masked', numel(fly_rot{1,1}), numel(fly_rot{2,1}), ...
    strjoin(arrayfun(@num2str,stim_nums,'UniformOutput',false),','), rho_min)}, 'FontSize', 11)

if ~isfolder(export_dir); mkdir(export_dir); end
exportgraphics(fig, fullfile(export_dir,sprintf('fig_p2x2_displacement_per_stim_rho%03d_stims%s.png',round(100*rho_min),sprintf('%d',stim_nums))), 'Resolution', 200)

%% functions
function name = fly_folder(meta)
    parts = strsplit(meta,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    k = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')),1);
    if isempty(k); name = ''; else; name = parts{k}; end
end
function fid = fly_path(meta)
    parts = strsplit(meta,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    k = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')),1);
    if isempty(k); fid = meta; else; fid = strjoin(parts(1:k),'/'); end
end
function m = masked_unwrapped_mu(mu, rho, t, rho_min)
    % bump position with low-rho frames removed: NaN where rho < rho_min,
    % unwrapped across the remaining (solid) frames only, then linearly
    % interpolated through the gaps; gaps at the ends take the nearest solid
    % value. NaN throughout if fewer than 2 solid frames.
    ok = rho >= rho_min & ~isnan(mu);
    m = nan(size(mu));
    if sum(ok) < 2; return; end
    m(ok) = unwrap(mu(ok));
    m = interp1(t(ok), m(ok), t, 'linear');
    m = fillmissing(m, 'nearest');
end
