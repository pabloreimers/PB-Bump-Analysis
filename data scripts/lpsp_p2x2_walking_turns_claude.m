%% lpsp_p2x2_walking_turns_claude
% Turn-size-resolved bump motion for the LPsP>P2X2 ATP-iontophoresis
% walking dataset (.data/lpsp_p2x2_walking_20260728_post0708.mat, the same
% subset Figure_3_claude.m's bump-mobility row uses). A port of QUESTION 1
% of data scripts/lpsp_rnai_turns_claude.m to this dataset:
%
%   Every discrete turn the fly makes is detected from its own rotational
%   velocity (ft.r_speed), sized by how far the fly rotated over the turn
%   (integrated r_speed), and binned by that size (15, 30, 45, ... deg).
%   Within each bin the fly's heading and the bump's position are aligned
%   to turn onset, sign-flipped so every turn is positive, and overlaid
%   (mean +/- SEM). A summary then plots how far the bump moved vs. how far
%   the fly moved per turn-size bin, per fly.
%
% What is different from the RNAi version:
%   * CONDITIONS. Each fly's trials are split into PRE (before its first
%     detected ATP-pulse trial) and POST (after it, pulse trials excluded)
%     -- same detector/logic as Figure_3_claude.m / lpsp_p2x2_walking_
%     script.m section 3 -- and into closed loop vs dark by ft.pattern. So
%     the grid is genotype (empty>p2x2, lpsp>p2x2) x phase (pre, post) x
%     light (closed loop, dark). Flies with no detected pulse trial are
%     dropped (no pre/post for them).
%   * BUMP ESTIMATE. The RNAi script recomputed the PVA from im.f (that
%     dataset was built unsmoothed). Here we use the dataset's own im.mu /
%     im.rho, i.e. exactly the bump the Figure 3 mobility metric is built
%     on (imgData_reg -> 10-frame Gaussian -> wedge dF/F -> z -> PVA, see
%     lpsp_p2x2_script_minimal.m), plus the same 0.2 s Gaussian on the
%     unwrapped trace as the RNAi turn analysis. Set recompute_pva = true
%     to use the RNAi script's im.f movmean-3 PVA instead.
%   * HEADING. This dataset has no ft.heading; heading = cumtrapz(r_speed),
%     which is what Figure 3 uses too. Sign checked directly: integrated
%     r_speed and -ft.cue move together (corr of increments 0.7-0.9 on the
%     example closed-loop trials), and im.mu tracks -ft.cue, so bump and
%     fly displacements are directly comparable with no sign flip.
%   * CLOSED-LOOP CUE GAIN is not assumed to be 0.8 -- it is estimated per
%     closed-loop trial (slope of cue increments vs heading increments) and
%     the median is drawn as the dashed reference line.
%
%   * LATERALIZATION. Each ATP ejection targets one PB hemisphere. Per fly
%     the stimulated side is called from the atp channel the same way as
%     lpsp_p2x2_walking_script.m section 16a: mean atp.d over centroids
%     1:16 vs 17:32 in the 0-3 s post-stim window of its pulse trials;
%     1:16 brighter = "right" (the codebase's right_idx convention -- not
%     independently verified against anatomy, but applied consistently).
%     RIGHT-stimulated flies are then FLIPPED into the left-stimulated
%     frame by negating ft.r_speed (hence heading and turn direction),
%     ft.cue and im.mu. Nothing else in this analysis depends on side:
%     rho, f_speed, bout detection and the rho gate are sign-blind. (If
%     im.z heatmaps are ever drawn in this frame, glomerulus rows would
%     need mirroring too -- not done here.) Turns are then split by
%     direction: "left" = positive r_speed in the aligned frame (CCW under
%     FicTrac's velYaw convention, i.e. toward the stimulated side IF that
%     convention holds on this rig -- same unverified assumption as the
%     walking script's 16b), "right" = negative. Empty flies get the same
%     side call and flip so both genotypes sit in the same frame.
%
% !!! RHO / MOBILITY CAVEAT (see Figure_3_claude.m, "RHO AND MOBILITY ARE
% COUPLED"): lpsp post-ATP dark trials have lower bump vector strength, and
% lower rho goes with higher measured bump movement. Every turn here is
% gated on mean rho over the turn window > turn_rho_thresh (same as the
% RNAi analysis), which partly protects per-turn displacements, but any
% pre->post change in bump displacement should be reported alongside the
% rho change. The summary figure's third panel shows per-condition rho.
%
% Run from the repo root one %% section at a time or top to bottom.
% Figures export as PNG to ugly_figures/exports/walking_turns/.
repo_root = fileparts(fileparts(mfilename('fullpath'))); % data scripts -> repo root; paths below are anchored here so run() from anywhere works
set(groot,'defaultAxesToolbarVisible','off');

%% 0) parameters
data_dir    = fullfile(repo_root,'.data');
source_file = 'lpsp_p2x2_walking_20260728_post0708.mat';
export_dir  = fullfile(repo_root,'ugly_figures','exports','walking_turns');
if ~exist(export_dir,'dir'), mkdir(export_dir); end
export_prefix = 'p2x2turns';

geno_keep   = {'empty>p2x2','lpsp>p2x2'};
geno_colors = [0.00,0.45,0.74;   % empty  blue
               0.85,0.33,0.10];  % lpsp   orange
phase_name  = {'pre','post'};
dir_name    = {'all turns','left turns (+r_speed, aligned frame)','right turns (-r_speed, aligned frame)'}; dir_short = {'all','left','right'}; dir_sign = [0, 1, -1];
flip_to_left = true; % flip right-stimulated flies into the left-stimulated frame (false = analyze unflipped, as before)
cond_name   = {'closed loop','dark'};

recompute_pva   = false; % false = dataset im.mu/im.rho (Figure 3's bump); true = RNAi-style movmean-3 PVA from im.f
im_f_frames     = 3;
bump_smooth_s   = 0.20; % s, Gaussian on unwrapped mu (imaging timebase)
bump_lag_s      = 0.13; % s, bump lags behavior (gain_scratch_claude.m)

% turn detection -- identical to lpsp_rnai_turns_claude.m
turn_vel_smooth_s = 0.50; turn_vel_thresh = 0.30; turn_max_gap_s = 0.20; turn_min_dur_s = 0.15;
turn_rho_thresh   = 0.20; turn_win_s = [-1, 3]; turn_grid_dt = 0.02; turn_quiet_pre_s = 0.50;
turn_norm_tau     = (-0.5:0.02:1.5)';
turn_bin_edges_deg = [7.5:15:97.5, 135, 180, 360];
min_turns_per_bin  = 3;

% pulse-trial detection -- identical to Figure_3_claude.m / walking script section 3
pulse_peri_win = 5; pulse_base_win = [-5,0]; pulse_peak_win = [0,3]; pulse_peak_factor = 3;

%% 1) load + labels
if ~exist('all_data','var') || ~exist('walk_src','var') || ~strcmp(walk_src,'post0708')
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:); clear tmp
    walk_src = 'post0708';
end
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

genotype = repmat({'empty>p2x2'},n_trials,1);
genotype(arrayfun(@(x)contains(x.meta,'_lpsp_'), all_data)) = {'lpsp>p2x2'};
is_dark  = arrayfun(@(x)contains(x.ft.pattern,'background'), all_data); is_dark = is_dark(:);
fly_id = cell(n_trials,1); trial_num = nan(n_trials,1); has_pulse = false(n_trials,1);
for i = 1:n_trials
    [fly_id{i}, trial_num(i)] = walk_meta_fly_trial(all_data(i).meta);
    has_pulse(i) = walk_trial_has_pulse(all_data(i), pulse_peri_win, pulse_base_win, pulse_peak_win, pulse_peak_factor);
end
[fly_list,~,fly_num] = unique(fly_id);
n_flies = numel(fly_list);

% phase per trial: 1 = pre (before first pulse trial), 2 = post (after it, non-pulse), NaN = pulse trial or fly without pulse
phase = nan(n_trials,1);
fly_geno = cell(n_flies,1);
for f = 1:n_flies
    idx = find(fly_num==f); [~,o] = sort(trial_num(idx)); idx = idx(o);
    fly_geno{f} = genotype{idx(1)};
    p = find(has_pulse(idx),1);
    if isempty(p), continue, end
    phase(idx(1:p-1)) = 1;
    post = idx(p+1:end); phase(post(~has_pulse(post))) = 2;
end
use_trials = find(~isnan(phase));
fprintf('%d flies (%d empty, %d lpsp); %d/%d trials usable (pre %d, post %d); %d flies have no pulse trial and are dropped\n', ...
    n_flies, sum(strcmp(fly_geno,'empty>p2x2')), sum(strcmp(fly_geno,'lpsp>p2x2')), numel(use_trials), n_trials, ...
    sum(phase==1), sum(phase==2), sum(arrayfun(@(f) ~any(has_pulse(fly_num==f)), 1:n_flies)));

% stimulated side per fly (walking script 16a logic) -> flip sign per trial
fly_side = nan(n_flies,1); side_margin = nan(n_flies,1); % +1 = "right" (centroids 1:16 brighter), -1 = "left"
for f = 1:n_flies
    pt = find(fly_num==f & has_pulse);
    if isempty(pt), continue, end
    h1 = []; h2 = [];
    for i = pt'
        stims = logical(all_data(i).ft.stims(:)); onsets = find(diff([false;stims])==1);
        xb = all_data(i).ft.xb(:); t_on = all_data(i).ft.xf(onsets);
        for k = 1:numel(onsets)
            in_peak = (xb-t_on(k)) > pulse_peak_win(1) & (xb-t_on(k)) <= pulse_peak_win(2);
            h1(end+1) = mean(all_data(i).atp.d(1:16,in_peak),'all'); h2(end+1) = mean(all_data(i).atp.d(17:32,in_peak),'all'); %#ok<AGROW>
        end
    end
    fly_side(f) = 2*(mean(h1,'omitnan') > mean(h2,'omitnan')) - 1;
    side_margin(f) = (mean(h1,'omitnan') - mean(h2,'omitnan')) / (abs(mean(h1,'omitnan')) + abs(mean(h2,'omitnan')) + eps); % -1..1, 0 = ambiguous
end
trial_flip = ones(n_trials,1);
if flip_to_left, trial_flip(fly_side(fly_num)==1) = -1; end
fprintf('stim side: %d right (flipped), %d left (reference), %d undetermined; |side margin| median %.2f, min %.2f\n', ...
    sum(fly_side==1), sum(fly_side==-1), sum(isnan(fly_side)), median(abs(side_margin),'omitnan'), min(abs(side_margin)));
for gi = 1:2, fprintf('  %s: %d right, %d left\n', geno_keep{gi}, sum(fly_side==1 & strcmp(fly_geno,geno_keep{gi})), sum(fly_side==-1 & strcmp(fly_geno,geno_keep{gi}))); end

% closed-loop cue gain: slope of cue increments vs heading increments, per closed-loop trial
cue_gain_trial = nan(n_trials,1);
for i = find(~is_dark)'
    xf = all_data(i).ft.xf(:); dt = median(diff(xf));
    h = cumtrapz(xf, fillmissing(all_data(i).ft.r_speed(:),'constant',0));
    c = unwrap(-all_data(i).ft.cue(:));
    dh = diff(smoothdata(h,'gaussian',round(.5/dt))); dc = diff(smoothdata(c,'gaussian',round(.5/dt)));
    ok = ~isnan(dh) & ~isnan(dc) & abs(dh) > 1e-4;
    if sum(ok) > 1000, cue_gain_trial(i) = dh(ok) \ dc(ok); end
end
cue_gain = median(cue_gain_trial,'omitnan');
fprintf('closed-loop cue gain (cue increments vs heading increments): median %.2f, IQR %.2f-%.2f over %d trials\n', ...
    cue_gain, prctile(cue_gain_trial,25), prctile(cue_gain_trial,75), sum(~isnan(cue_gain_trial)));

%% 2) bump position per trial
pva_mu = cell(n_trials,1); pva_rho = cell(n_trials,1);
for ti = use_trials'
    if recompute_pva
        [pva_mu{ti},pva_rho{ti}] = trial_pva_movmean(all_data(ti),im_f_frames);
    else
        pva_mu{ti} = all_data(ti).im.mu(:); pva_rho{ti} = all_data(ti).im.rho(:);
    end
end

%% 3) detect turns, extract aligned traces + per-turn displacements
turn_t = (turn_win_s(1):turn_grid_dt:turn_win_s(2))';
turns = [];
for ti = use_trials'
    trial = all_data(ti); % flipped copy for right-stimulated flies: r_speed (-> heading, turn direction), cue, mu
    trial.ft.r_speed = trial_flip(ti)*trial.ft.r_speed; trial.ft.cue = trial_flip(ti)*trial.ft.cue;
    tr = trial_turns(trial,trial_flip(ti)*pva_mu{ti},pva_rho{ti},bump_smooth_s,bump_lag_s, ...
        turn_vel_smooth_s,turn_vel_thresh,turn_max_gap_s,turn_min_dur_s,turn_rho_thresh,turn_win_s,turn_t,turn_quiet_pre_s,turn_norm_tau);
    if isempty(tr), continue, end
    [tr.trial] = deal(ti); [tr.fly] = deal(fly_num(ti)); [tr.geno] = deal(genotype{ti}); [tr.dark] = deal(is_dark(ti)); [tr.phase] = deal(phase(ti));
    turns = [turns; tr(:)]; %#ok<AGROW>
end
assert(~isempty(turns), 'no turns detected');
turn_size  = [turns.size_deg]';
turn_bin   = discretize(turn_size,turn_bin_edges_deg);
n_bins     = numel(turn_bin_edges_deg)-1;
bin_labels = cell(1,n_bins);
for b = 1:n_bins
    if diff(turn_bin_edges_deg(b:b+1)) <= 15.01, bin_labels{b} = sprintf('%g%s',mean(turn_bin_edges_deg(b:b+1)),char(176));
    else, bin_labels{b} = sprintf('%g-%g%s',turn_bin_edges_deg(b),turn_bin_edges_deg(b+1),char(176)); end
end
turn_fly = [turns.fly]'; turn_dark = [turns.dark]'; turn_phase = [turns.phase]'; turn_quiet = [turns.quiet_pre]'; turn_rho = [turns.rho_mean]'; turn_dir = [turns.dir]';
fprintf('%d turns detected (%d within the size bins); per bin: %s\n', numel(turns), sum(~isnan(turn_bin)), ...
    strjoin(arrayfun(@(b) sprintf('%s=%d',bin_labels{b},sum(turn_bin==b)),1:n_bins,'UniformOutput',false),', '));
for ci = 1:2, for ph = 1:2, for gi = 1:2
    sel = turn_dark==(ci==2) & turn_phase==ph & strcmp({turns.geno}',geno_keep{gi});
    fprintf('  %-11s %-4s %-10s: %5d turns (%d left / %d right) from %2d flies, mean rho over turn %.2f\n', cond_name{ci}, phase_name{ph}, geno_keep{gi}, sum(sel), sum(sel&turn_dir==1), sum(sel&turn_dir==-1), numel(unique(turn_fly(sel))), mean(turn_rho(sel)));
end, end, end

%% 4) per fly x light x phase x direction x bin: mean displacements + through-origin slope
n_dir = numel(dir_sign);
n_turns_fb = zeros(n_flies,2,2,n_dir,n_bins); dfly_fb = nan(n_flies,2,2,n_dir,n_bins); dbump_fb = dfly_fb; ratio_fb = dfly_fb; rho_fb = dfly_fb;
for f = 1:n_flies
    for ci = 1:2
        for ph = 1:2
            for di = 1:n_dir
                sel_fc = turn_fly==f & turn_dark==(ci==2) & turn_phase==ph;
                if dir_sign(di) ~= 0, sel_fc = sel_fc & turn_dir==dir_sign(di); end
                for b = 1:n_bins
                    sel = find(sel_fc & turn_bin==b);
                    n_turns_fb(f,ci,ph,di,b) = numel(sel);
                    if numel(sel) < min_turns_per_bin, continue, end
                    dfly = [turns(sel).dfly_deg]'; dbump = [turns(sel).dbump_deg]';
                    dfly_fb(f,ci,ph,di,b) = mean(dfly); dbump_fb(f,ci,ph,di,b) = mean(dbump); ratio_fb(f,ci,ph,di,b) = dfly \ dbump;
                    rho_fb(f,ci,ph,di,b) = mean(turn_rho(sel));
                end
            end
        end
    end
end

%% figure 1: group aligned overlays -- rows = genotype x phase, cols = turn-size bin; one figure per light condition x turn direction
row_labels = {}; row_sel = {}; row_col = [];
for gi = 1:2, for ph = 1:2
    row_labels{end+1} = sprintf('%s %s', geno_keep{gi}, upper(phase_name{ph})); %#ok<AGROW>
    row_sel{end+1} = struct('gi',gi,'ph',ph); %#ok<AGROW>
    row_col(end+1,:) = geno_colors(gi,:) * (1 - .45*(ph==1)) + .45*(ph==1)*[1 1 1]; %#ok<AGROW> pre = lighter tint, post = full color
end, end
flip_str = ''; if flip_to_left, flip_str = ' [right-stim flies flipped to left-stim frame]'; end
for ci = 1:2
    for di = 1:n_dir
        figure(10*ci+di); clf
        set(gcf,'Name',sprintf('turn overlays, %s, %s',cond_name{ci},dir_name{di}),'Position',[40,40,210*n_bins,200*numel(row_labels)],'Color','w')
        tl = tiledlayout(numel(row_labels),n_bins,'TileSpacing','compact','Padding','compact');
        for r = 1:numel(row_labels)
            gi = row_sel{r}.gi; ph = row_sel{r}.ph;
            sel = turn_dark==(ci==2) & turn_phase==ph & strcmp({turns.geno}',geno_keep{gi});
            if dir_sign(di) ~= 0, sel = sel & turn_dir==dir_sign(di); end
            FM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,'fly_trace',min_turns_per_bin);
            BM = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,'bump_trace',min_turns_per_bin);
            for b = 1:n_bins
                ax = nexttile(tl,(r-1)*n_bins+b); hold(ax,'on')
                FT = reshape(FM(:,b,:),n_flies,[]); BT = reshape(BM(:,b,:),n_flies,[]);
                n_f = sum(~all(isnan(FT),2));
                if n_f >= 1
                    plot_sem(turn_t,FT,[0,0,0]); plot_sem(turn_t,BT,row_col(r,:));
                    if ci == 1, plot(turn_t,cue_gain*mean(FT,1,'omitnan'),'--','Color',[.5,.5,.5]); end
                end
                xline(ax,0,':k'); xlim(ax,turn_win_s)
                yl = turn_bin_edges_deg(b+1)*1.3; if b == n_bins, yl = 250; end; ylim(ax,[-yl*.3, yl])
                title(ax,sprintf('%s (n=%d flies)',bin_labels{b},n_f),'FontSize',8)
                if b == 1, ylabel(ax,sprintf('%s\nrotation (deg)',row_labels{r}),'Interpreter','none'), end
                if r == numel(row_labels), xlabel(ax,'time from onset (s)'), end
                set(ax,'FontSize',7)
            end
        end
        ref_str = ''; if ci == 1, ref_str = sprintf('; dashed = %.2f x fly (measured cue gain)', cue_gain); end
        title(tl,sprintf('%s, %s%s: fly heading (black) vs bump (color; pre = light, post = dark), mean +/- SEM across flies (each fly = mean of its >= %d turns per bin)%s', ...
            cond_name{ci},dir_name{di},flip_str,min_turns_per_bin,ref_str),'Interpreter','none','FontSize',10)
        exportgraphics(gcf,fullfile(export_dir,sprintf('%s_fig1_overlay_%s_%s.png',export_prefix,strrep(cond_name{ci},' ',''),dir_short{di})),'Resolution',120)
    end
end

%% figure 2: gain per turn-size bin, pre vs post per genotype -- rows = light, cols = all / left / right turns
figure(30); clf; set(gcf,'Name','gain per bin by turn direction','Position',[40,40,1700,800],'Color','w')
tl = tiledlayout(2,n_dir,'TileSpacing','compact','Padding','compact');
ph_style = {'--','-'};
for ci = 1:2
    for di = 1:n_dir
        ax = nexttile(tl,(ci-1)*n_dir+di); hold(ax,'on')
        for gi = 1:2
            flies_g = find(strcmp(fly_geno,geno_keep{gi}));
            for ph = 1:2
                R = reshape(ratio_fb(flies_g,ci,ph,di,:),numel(flies_g),n_bins);
                xj = (1:n_bins) + ((gi-1)*2+ph - 2.5)*0.18;
                plot(ax,repmat(xj,size(R,1),1)',R','.','Color',[geno_colors(gi,:),.3],'MarkerSize',8,'HandleVisibility','off')
                hE = errorbar(ax,xj,mean(R,1,'omitnan'),sem(R,1),[ph_style{ph} 'o'],'Color',geno_colors(gi,:),'LineWidth',1.8,'DisplayName',sprintf('%s %s (n=%d)',geno_keep{gi},phase_name{ph},sum(any(~isnan(R),2))));
                if ph == 2, hE.MarkerFaceColor = geno_colors(gi,:); end
            end
        end
        yline(ax,1,':k','HandleVisibility','off'); if ci == 1, yline(ax,cue_gain,'--','Color',[.5,.5,.5],'HandleVisibility','off'); end; yline(ax,0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
        xticks(ax,1:n_bins); xticklabels(ax,bin_labels); xlim(ax,[.5,n_bins+.5]); ylim(ax,[-1, 3])
        xlabel(ax,'turn size (fly rotation)'); if di == 1, ylabel(ax,'bump / fly displacement (per-fly slope)'), end
        legend(ax,'Location','northwest','Interpreter','none','FontSize',7,'Box','off')
        title(ax,sprintf('%s: %s',cond_name{ci},dir_name{di}),'Interpreter','none','FontSize',9)
    end
end
title(tl,['gain per turn-size bin, pre (dashed/open) vs post (solid/filled) ATP, by genotype and turn direction' flip_str],'Interpreter','none')
exportgraphics(gcf,fullfile(export_dir,sprintf('%s_fig2_gain_by_direction.png',export_prefix)),'Resolution',120)

%% figure 2b: mean rho over the turn per bin, same layout (the caveat panel)
figure(32); clf; set(gcf,'Name','rho per bin by turn direction','Position',[40,40,1700,700],'Color','w')
tl = tiledlayout(2,n_dir,'TileSpacing','compact','Padding','compact');
for ci = 1:2
    for di = 1:n_dir
        ax = nexttile(tl,(ci-1)*n_dir+di); hold(ax,'on')
        for gi = 1:2
            flies_g = find(strcmp(fly_geno,geno_keep{gi}));
            for ph = 1:2
                R = reshape(rho_fb(flies_g,ci,ph,di,:),numel(flies_g),n_bins);
                xj = (1:n_bins) + ((gi-1)*2+ph - 2.5)*0.18;
                hE = errorbar(ax,xj,mean(R,1,'omitnan'),sem(R,1),[ph_style{ph} 'o'],'Color',geno_colors(gi,:),'LineWidth',1.8,'DisplayName',sprintf('%s %s',geno_keep{gi},phase_name{ph}));
                if ph == 2, hE.MarkerFaceColor = geno_colors(gi,:); end
            end
        end
        yline(ax,turn_rho_thresh,':r','HandleVisibility','off')
        xticks(ax,1:n_bins); xticklabels(ax,bin_labels); xlim(ax,[.5,n_bins+.5]); ylim(ax,[.2 .9])
        xlabel(ax,'turn size'); if di == 1, ylabel(ax,'mean rho over turn'), end
        legend(ax,'Location','northwest','Interpreter','none','FontSize',7,'Box','off')
        title(ax,sprintf('%s: %s',cond_name{ci},dir_name{di}),'Interpreter','none','FontSize',9)
    end
end
title(tl,'bump vector strength (rho) over the turn, per bin -- CAVEAT: rho and bump movement are coupled (see Figure_3_claude.m)','Interpreter','none')
exportgraphics(gcf,fullfile(export_dir,sprintf('%s_fig2b_rho_by_direction.png',export_prefix)),'Resolution',120)

%% figure 3: paired per-fly change in gain (post - pre) per bin, by direction
figure(31); clf; set(gcf,'Name','paired gain change by direction','Position',[40,40,1700,700],'Color','w')
tl = tiledlayout(2,n_dir,'TileSpacing','compact','Padding','compact');
for ci = 1:2
    for di = 1:n_dir
        ax = nexttile(tl,(ci-1)*n_dir+di); hold(ax,'on')
        for gi = 1:2
            flies_g = find(strcmp(fly_geno,geno_keep{gi}));
            D = reshape(ratio_fb(flies_g,ci,2,di,:) - ratio_fb(flies_g,ci,1,di,:),numel(flies_g),n_bins);
            xj = (1:n_bins) + (gi-1.5)*0.25;
            plot(ax,repmat(xj,size(D,1),1)',D','.','Color',[geno_colors(gi,:),.35],'MarkerSize',9,'HandleVisibility','off')
            errorbar(ax,xj,mean(D,1,'omitnan'),sem(D,1),'o-','Color',geno_colors(gi,:),'MarkerFaceColor',geno_colors(gi,:),'LineWidth',1.8, ...
                'DisplayName',sprintf('%s (n=%d paired flies)',geno_keep{gi},max(sum(~isnan(D),1))))
            for b = 1:n_bins
                d = D(~isnan(D(:,b)),b);
                if numel(d) >= 5, p = signrank(d); if p < .05, text(ax,xj(b),max(d)+.1,sprintf('p=%.3f',p),'Color',geno_colors(gi,:),'FontSize',7,'HorizontalAlignment','center'); end, end
            end
        end
        yline(ax,0,'-k','HandleVisibility','off'); xticks(ax,1:n_bins); xticklabels(ax,bin_labels); xlim(ax,[.5,n_bins+.5]); ylim(ax,[-2,2.5])
        xlabel(ax,'turn size (fly rotation)'); if di == 1, ylabel(ax,'post - pre gain'), end
        legend(ax,'Location','northwest','Interpreter','none','FontSize',7,'Box','off')
        title(ax,sprintf('%s: %s',cond_name{ci},dir_name{di}),'Interpreter','none','FontSize',9)
    end
end
title(tl,'per-fly change in bump/fly gain after the ATP pulse trial, per turn-size bin and direction (p: signrank vs 0, shown when < .05; 9 bins x 12 panels, so expect a few by chance)','Interpreter','none')
exportgraphics(gcf,fullfile(export_dir,sprintf('%s_fig3_paired_gain_change_by_direction.png',export_prefix)),'Resolution',120)

%% figure 4: the lateralization test -- per fly, gain(left turns) - gain(right turns) per bin, pre vs post, by genotype
% if the ATP ejection biases the bump toward/away from the stimulated
% side, the post-ATP left-right asymmetry should differ from pre in lpsp
% but not empty flies. Pre = dashed/open, post = solid/filled.
figure(33); clf; set(gcf,'Name','left-right gain asymmetry','Position',[40,40,1200,480],'Color','w')
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for ci = 1:2
    ax = nexttile(tl,ci); hold(ax,'on')
    for gi = 1:2
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        for ph = 1:2
            A = reshape(ratio_fb(flies_g,ci,ph,2,:) - ratio_fb(flies_g,ci,ph,3,:),numel(flies_g),n_bins);
            xj = (1:n_bins) + ((gi-1)*2+ph - 2.5)*0.18;
            plot(ax,repmat(xj,size(A,1),1)',A','.','Color',[geno_colors(gi,:),.3],'MarkerSize',8,'HandleVisibility','off')
            hE = errorbar(ax,xj,mean(A,1,'omitnan'),sem(A,1),[ph_style{ph} 'o'],'Color',geno_colors(gi,:),'LineWidth',1.8,'DisplayName',sprintf('%s %s (n=%d)',geno_keep{gi},phase_name{ph},max(sum(~isnan(A),1))));
            if ph == 2, hE.MarkerFaceColor = geno_colors(gi,:); end
        end
        dA = reshape((ratio_fb(flies_g,ci,2,2,:) - ratio_fb(flies_g,ci,2,3,:)) - (ratio_fb(flies_g,ci,1,2,:) - ratio_fb(flies_g,ci,1,3,:)),numel(flies_g),n_bins);
        for b = 1:n_bins
            d = dA(~isnan(dA(:,b)),b);
            if numel(d) >= 5, p = signrank(d); if p < .05, text(ax,b+(gi-1.5)*.3,2.2-(gi-1)*.25,sprintf('p=%.3f',p),'Color',geno_colors(gi,:),'FontSize',7,'HorizontalAlignment','center'); end, end
        end
    end
    yline(ax,0,'-k','HandleVisibility','off'); xticks(ax,1:n_bins); xticklabels(ax,bin_labels); xlim(ax,[.5,n_bins+.5]); ylim(ax,[-2.5,2.5])
    xlabel(ax,'turn size (fly rotation)'); ylabel(ax,'gain(left turns) - gain(right turns)')
    legend(ax,'Location','southwest','Interpreter','none','FontSize',7,'Box','off')
    title(ax,sprintf('%s: left-right gain asymmetry per fly',cond_name{ci}))
end
title(tl,'left vs right turn gain asymmetry in the left-stim-aligned frame, pre (dashed/open) vs post (solid/filled); p = signrank on the per-fly post-pre change in asymmetry (shown when < .05)','Interpreter','none')
exportgraphics(gcf,fullfile(export_dir,sprintf('%s_fig4_left_right_asymmetry.png',export_prefix)),'Resolution',120)

%% report
for di = 1:n_dir
    fprintf('\n== per-fly gain (bump/fly slope), mean over flies -- %s ==\n', dir_name{di})
    fprintf('%-11s %-4s %-10s', 'light','phase','geno'); fprintf(' %7s', bin_labels{:}); fprintf('\n')
    for ci = 1:2, for gi = 1:2, for ph = 1:2
        flies_g = find(strcmp(fly_geno,geno_keep{gi}));
        R = reshape(ratio_fb(flies_g,ci,ph,di,:),numel(flies_g),n_bins);
        fprintf('%-11s %-4s %-10s', cond_name{ci}, phase_name{ph}, geno_keep{gi}); fprintf(' %7.2f', mean(R,1,'omitnan')); fprintf('\n')
    end, end, end
end

%% functions -- ported from lpsp_rnai_turns_claude.m unless noted
function [fid, tnum] = walk_meta_fly_trial(meta_path) % from Figure_3_claude.m
    parts = strsplit(meta_path,{'\','/'}); parts(cellfun(@isempty,parts)) = [];
    fly_pos = find(startsWith(parts,'fly '),1);
    fid = strjoin(parts(1:fly_pos),'/');
    tnum = str2double(regexp(parts{fly_pos+1},'-(\d+)[_-]','tokens','once'));
end

function tf = walk_trial_has_pulse(trial, peri_win, base_win, peak_win, peak_factor) % from Figure_3_claude.m
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

function [mu_new,rho_new] = trial_pva_movmean(trial,smooth_frames)
    f_smooth = movmean(trial.im.f,smooth_frames,2);
    f_z = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);
    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp] = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function h = trial_heading(trial)
    % no ft.heading in this dataset: integrate r_speed (what Figure 3 does)
    r = trial.ft.r_speed(:); r(isnan(r)) = 0;
    h = cumtrapz(trial.ft.xf(:),r);
end

function [mu_f,rho_f,xb] = bump_on_fictrac(trial,mu_raw,rho_raw,mu_smooth_s)
    xf = trial.ft.xf(:); n_im = numel(mu_raw);
    if isfield(trial.ft,'xb') && numel(trial.ft.xb) == n_im, xb = trial.ft.xb(:); else, xb = linspace(xf(1),xf(end),n_im)'; end
    dt_im = median(diff(xb)); win_im = max(1,round(mu_smooth_s/dt_im));
    mu_s = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    mu_f  = interp1(xb,mu_s,xf,'linear');
    rho_f = interp1(xb,rho_raw(:),xf,'linear');
end

function turns = trial_turns(trial,mu_raw,rho_raw,mu_smooth_s,lag_s,vel_smooth_s,vel_thresh,max_gap_s,min_dur_s,rho_thresh,win_s,t_grid,quiet_pre_s,tau_grid)
    xf = trial.ft.xf(:); dt = median(diff(xf));
    heading = trial_heading(trial);
    [mu_f,rho_f] = bump_on_fictrac(trial,mu_raw,rho_raw,mu_smooth_s);
    win = max(1,round(vel_smooth_s/dt));
    r_sm = smoothdata(trial.ft.r_speed(:),'gaussian',win); r_sm(isnan(r_sm)) = 0;
    mu_vel = smoothdata(gradient(mu_f)/dt,'gaussian',win);
    max_gap = max(1,round(max_gap_s/dt)); min_dur = max(1,round(min_dur_s/dt)); n_quiet = max(1,round(quiet_pre_s/dt));
    turns = struct('onset_s',{},'dur_s',{},'size_deg',{},'dir',{},'quiet_pre',{},'fly_trace',{},'bump_trace',{}, ...
                   'fly_vel_trace',{},'bump_vel_trace',{},'fly_norm',{},'bump_norm',{},'dfly_deg',{},'dbump_deg',{},'mu_onset',{},'rho_mean',{});
    for s = [1,-1]
        active = s*r_sm > vel_thresh;
        active = imclose(active,ones(max_gap+1,1));
        active = bwareaopen(active,min_dur);
        d = diff([false;active;false]); on = find(d==1); off = find(d==-1)-1;
        for k = 1:numel(on)
            t_on = xf(on(k)); t_off = xf(off(k));
            if t_on + win_s(1) < xf(1) || max(t_on + win_s(2), t_off + lag_s) > xf(end), continue, end
            rng = on(k):off(k);
            if mean(rho_f(rng),'omitnan') < rho_thresh, continue, end
            dfly = heading(off(k)) - heading(on(k));
            if abs(dfly) < 1e-6, continue, end
            mu_on  = interp1(xf,mu_f,t_on + lag_s); mu_off = interp1(xf,mu_f,t_off + lag_s);
            tq = t_on + t_grid; bump_tr = interp1(xf,mu_f,tq);
            if isnan(mu_on) || isnan(mu_off) || any(isnan(bump_tr)), continue, end
            tau_q = t_on + tau_grid*(t_off - t_on);
            mu_at_on = interp1(xf,mu_f,t_on);
            T.onset_s = t_on; T.dur_s = t_off - t_on; T.size_deg = rad2deg(abs(dfly)); T.dir = s;
            T.quiet_pre = on(k) > n_quiet && all(abs(r_sm(on(k)-n_quiet:on(k)-1)) < vel_thresh);
            T.fly_trace  = single(s*rad2deg(interp1(xf,heading,tq) - heading(on(k))));
            T.bump_trace = single(s*rad2deg(bump_tr - mu_at_on));
            T.fly_vel_trace  = single(s*rad2deg(interp1(xf,r_sm,tq)));
            T.bump_vel_trace = single(s*rad2deg(interp1(xf,mu_vel,tq)));
            T.fly_norm  = single(s*rad2deg(interp1(xf,heading,tau_q) - heading(on(k))));
            T.bump_norm = single(s*rad2deg(interp1(xf,mu_f,tau_q) - mu_at_on));
            T.dfly_deg  = s*rad2deg(dfly); T.dbump_deg = s*rad2deg(mu_off - mu_on);
            T.mu_onset  = mod(mu_at_on+pi,2*pi)-pi; T.rho_mean = mean(rho_f(rng),'omitnan');
            turns(end+1) = T; %#ok<AGROW>
        end
    end
    [~,order] = sort([turns.onset_s]); turns = turns(order);
end

function M = fly_bin_means(turns,sel,turn_bin,turn_fly,n_flies,n_bins,field,min_turns)
    n_t = numel(turns(1).(field)); M = nan(n_flies,n_bins,n_t);
    for f = unique(turn_fly(sel))'
        for b = 1:n_bins
            idx = find(sel & turn_fly==f & turn_bin==b);
            if numel(idx) < min_turns, continue, end
            M(f,b,:) = mean(cat(2,turns(idx).(field))',1,'omitnan');
        end
    end
end

function s = sem(X,dim)
    s = std(X,0,dim,'omitnan') ./ sqrt(max(1,sum(~isnan(X),dim)));
end

function plot_sem(t,X,c)
    t = t(:)'; m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m); if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end
