%% fig_rnai_th_lapping_examples
% "Lapping": the bump over-rotates the fly by a full circle or more -- it
% keeps sweeping around the PB, in the direction the fly turned, while the
% fly turns only a fraction of that (e.g. 20250923 fly 4 trial 11, ~290-310
% s: fly ~+100 deg, bump ~+600 deg). This finds lapping episodes in closed-
% loop (light_cond = 'closed') or dark (light_cond = 'dark') trials and shows
% one per main-genotype fly (the named example first): z-scored PB heatmap with
% fly heading and bump PVA drawn on each PB half, and unwrapped rotation below.
% The same detector is used in both conditions; it measures the bump against
% the fly's heading (FicTrac), which in closed loop also drives the bar.
%
% LAP = the bump over-rotates the fly by >= lap_min_excess_deg (150), plus the
% quality checks below. (Earlier versions also required the bump to go >= 1.5x
% the fly's turn; dropped -- the absolute 150 deg minimum is what separates
% laps from the small post-turn over-rotations seen in controls.)
%
% Detector, per trial, sliding windows of each length in lap_win_s (step
% lap_step_s; 5 s catches fast laps, 20 s slow ones):
%   bump  = PVA of im.f (movmean 3 frames, z-scored per glomerulus),
%           unwrapped, 0.2 s Gaussian (same as fig_rnai_turns_examples.m)
%   fly   = unwrapped ft.heading on the imaging timebase
%   a window counts if
%     - the fly turns >= fly_min_deg and the bump moves the same way
%     - |bump| - |fly| >= lap_min_excess_deg (150)
%     - median PVA strength rho >= rho_min, no raw frame-to-frame PVA
%       jump > jump_max_deg (unwrap errors)
%     - the traces DIVERGE STEADILY: with d(t) = bump - fly (signed so a
%       lap is d rising), no step_win_s (0.3 s) stretch carries more than
%       step_frac_max of the net divergence (a relocation jump lands within a
%       few frames; a fast lap sweeps over >= ~1 s), and d's net change is >= straight_min of its
%       path length (0.5 s grid, so frame noise doesn't count as path).
%       This rejects a single bump jump with otherwise parallel traces,
%       which the same checks on the bump trace alone let through.
%     - the bump outpaces the fly THROUGHOUT the turn: over 1 s bins where
%       the fly turns >= bin_fly_min_deg (same direction), the median of
%       bump change / fly change is >= bin_ratio_min, with at least
%       min_turn_frac of the window's 1 s bins turning (>= 1), so a brief
%       turn can count in a short window (one jump bin can't carry the
%       median, and drift while the fly is still doesn't count)
%     - the bump still outpaces the fly (fly >= fly_min_deg, same direction)
%       counting only 0.5 s bins where
%       the heading moves >= move_min_deg: divergence that piles up while
%       the bar is stationary (fly still, or FicTrac stalled) can't make a lap
%   overlapping passing windows, of any length -> one episode (shown by its
%   max-excess window; win_s in the table says which length that was).
% The counts per genotype are printed; the controls are scanned the same way.
% A second figure shows control flies' episodes from the same detector, a
% third summarizes: an example lap, laps per fly and fraction of flies with a
% lap per genotype (and LPsP>TH high / low gain), with p-values.
% Trials: closed loop -- the integrative-gain inclusion (walking + working cue
% + bump); dark -- walking + bump (mean raw rho > 0.2), and windows where the
% fly's mounting bias (leak of forward walking into measured rotation, fit per
% trial) is removed from the heading (bias_correct; the heading is the
% integral of ft.r_speed, so the integrated leak slope x forward speed is
% subtracted) -- a fly mounted at an angle reads straight walking as turning,
% which the bump rightly ignores. With bias_correct = false, windows where the
% leak adds more than bias_max_deg are skipped instead.
% Exports are tagged _dark in the dark.
% Needs all_data (loaded if absent) and the integrative-gain cache (for
% trial labels and per-fly gain). Run from the repo root.

%% parameters + labels
lap_win_s = [5,20]; lap_step_s = 1; % window lengths (s) scanned; fast laps show in short windows, slow ones in long
fly_min_deg = 30; rho_min = 0.4; jump_max_deg = 90;
rho_min_dark = 0.4; % median PVA strength threshold used in the dark (0.3 tried: added mostly noisy weak-bump hits, more in controls)
light_cond = 'closed'; % 'closed' (closed-loop trials) or 'dark' (dark trials)
lap_min_excess_deg = 150; % LAP: the bump goes at least this far past the fly (|bump| - |fly|)
bias_correct = true; % dark only: subtract the fitted mounting bias (forward-speed leak into turning) from the heading
bias_max_deg = 20;   % dark, bias_correct = false: instead skip windows where the bias adds more than this much heading
step_win_s = 0.3; step_frac_max = 0.35; % no step_win_s stretch may carry more than this share of the divergence
%   (0.3 s: a relocation jump lands within a few frames -- share 0.42-0.91 -- while a
%   fast lap sweeping over ~1.5 s spreads out -- 0.22-0.33; confirmed laps <= 0.20)
straight_min = 0.4; % divergence net / path length
bin_fly_min_deg = 5; min_turn_frac = 0.2;   % 1 s bins counted as turning; at least this fraction of the window's bins (>= 1)
bin_ratio_min = 1;                          % median bump/fly ratio over those turning bins: the bump outpaces the fly during the turn
move_min_deg = 10; % per 0.5 s bin (20 deg/s): heading change below this = bar stationary (FicTrac jitter at rest can exceed 5 deg/s)
named_example = {'20250923','4','11'}; % date, fly, trial of the example to show first
pad_s = 3;       % s shown before / after each episode window
n_show = 8;      % flies shown per figure (one episode each, largest excess)
geno_set = 'th'; % 'th' or 'vglut': which set the example galleries and the summary show (all genotypes are always scanned)
switch geno_set
    case 'th'
        main_geno = 'lpsp>th'; main_name = 'LPsP > TH-RNAi'; ctl_geno = {'empty>th','lpsp>mcherry'}; set_tag = 'th';
        grp_def = {'LPsP > TH-RNAi','lpsp>th',[]; 'Empty > TH-RNAi','empty>th',[]; 'LPsP > mCherry','lpsp>mcherry',[]; ...
                   'LPsP>TH high gain','lpsp>th',true; 'LPsP>TH low gain','lpsp>th',false};
        grp_col = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.93,0.69,0.13; 0.60,0.15,0.00; 0.95,0.62,0.45];
        grp_x = [1,2,3,4.6,5.6]; cmp = [1,2; 1,3; 4,5];
    case 'vglut'
        main_geno = 'lpsp>vglut'; main_name = 'LPsP > VGlut-RNAi'; ctl_geno = {'empty>vglut','lpsp>mcherry'}; set_tag = 'vglut';
        grp_def = {'LPsP > VGlut-RNAi','lpsp>vglut',[]; 'Empty > VGlut-RNAi','empty>vglut',[]; 'LPsP > mCherry','lpsp>mcherry',[]};
        grp_col = [0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13];
        grp_x = 1:3; cmp = [1,2; 1,3];
end
gain_split = 0.85;
summary_example = {'20250905','5','12'}; % date, fly, trial of the lap shown in the summary figure (its largest episode)
scan_geno = {'lpsp>th','empty>th','lpsp>mcherry','lpsp>vglut','empty>vglut'};
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
if ~exist('all_data','var')
    tmp = load(fullfile(data_dir,'lpsp_rnai_joint_nosmooth_reg_20260831.mat'),'all_data'); all_data = tmp.all_data(:); clear tmp
end
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
g_fly = nan(max(IG.fly_num),1); % closed-loop integrative gain per fly
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
out_dir = fullfile('ugly_figures','rnai_turns','laps_stalls');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% detect lapping episodes
if strcmp(light_cond,'dark')
    cond_lbl = 'dark'; set_tag = [set_tag '_dark']; rho_min = rho_min_dark;
    cand = find(ismember(IG.genotype,scan_geno) & IG.is_dark); ok_t = false(size(cand));
    for j = 1:numel(cand), tr = all_data(cand(j)); ok_t(j) = sum((tr.ft.f_speed > .1) | (abs(tr.ft.r_speed)*2 > .1)) > numel(tr.ft.f_speed)*.2 && mean(tr.im.rho) > .2; end
    scan_trials = cand(ok_t);                                  % dark: walking + bump
else
    cond_lbl = 'closed loop';
    scan_trials = find(ismember(IG.genotype,scan_geno) & ~IG.is_dark & IG.inc); % closed loop: walking + working cue + bump (integrative-gain inclusion)
end
ev = zeros(0,8); % trial, t0, t1, dbump, dfly, excess, median rho, window length
for i = scan_trials(:)'
    [xb,mu_raw,mu_u,rho,hd] = trial_bump_heading(all_data(i));
    dt = median(diff(xb)); st = max(1,round(lap_step_s/dt)); n1 = max(1,round(1/dt)); n_half_s = max(1,round(0.5/dt)); n_step = max(1,round(step_win_s/dt));
    jump = [0; abs(diff(mu_raw))] > deg2rad(jump_max_deg);
    slope = 0; fb = zeros(size(xb));
    if strcmp(light_cond,'dark') % mounting bias: per-trial slope of turning on forward speed during walking
        trx = all_data(i); xf = trx.ft.xf(:); wf = max(1,round(0.5/median(diff(xf))));
        r = smoothdata(fillmissing(trx.ft.r_speed(:),'constant',0),'gaussian',wf); f = smoothdata(fillmissing(trx.ft.f_speed(:),'constant',0),'gaussian',wf);
        m = f > 1; if sum(m) > 200, bb = [ones(sum(m),1) f(m)] \ r(m); slope = bb(2); end
        fb = interp1(xf,f,xb);
        if bias_correct % heading = integral of r_speed, so remove the integrated leak slope x forward speed
            hd = hd - interp1(xf,cumtrapz(xf,slope*f),xb); slope = 0; % corrected: nothing left for the window filter below
        end
    end
    hits = zeros(0,7);
    for w_s = lap_win_s % every window length; hits are pooled and merged below
    nW = round(w_s/dt);
    for a = 1:st:numel(xb)-nW
        b = a + nW;
        db = rad2deg(mu_u(b)-mu_u(a)); df = rad2deg(hd(b)-hd(a));
        if isnan(df) || isnan(db) || abs(df) < fly_min_deg || sign(db) ~= sign(df), continue, end
        exc = abs(db) - abs(df);
        if exc < lap_min_excess_deg || median(rho(a:b)) < rho_min || any(jump(a:b)), continue, end
        if rad2deg(abs(slope*mean(fb(a:b))))*w_s > bias_max_deg, continue, end % dark: mounting bias could fake the mismatch
        % divergence d(t) = bump - fly, signed so a lap is d rising: it must
        % ramp, not step (no step_win_s stretch > step_frac_max of the net
        % divergence; net >= straight_min of its path on a 0.5 s grid) -- a bump
        % jump with otherwise parallel traces fails, a fast continuous sweep passes
        d = sign(df)*rad2deg((mu_u(a:b)-mu_u(a)) - (hd(a:b)-hd(a)));
        if any(isnan(d)), continue, end
        step_max = max(abs(d(1+n_step:end) - d(1:end-n_step)));
        d_ds = d(1:n_half_s:end);
        if step_max > step_frac_max*exc || exc/sum(abs(diff(d_ds))) < straight_min, continue, end
        % over-rotation throughout the turn: in 1 s bins where the fly turns
        % >= bin_fly_min_deg (at least min_turn_frac of the window's bins, so a
        % brief turn can count in a short window), median(bump / fly change) >= bin_ratio_min
        k_b = a:n1:b; dB1 = diff(rad2deg(mu_u(k_b))); dF1 = diff(rad2deg(hd(k_b)));
        turning = abs(dF1) >= bin_fly_min_deg & sign(dF1) == sign(df);
        if sum(turning) < max(1,ceil(min_turn_frac*numel(dF1))) || median(dB1(turning)./dF1(turning)) < bin_ratio_min, continue, end
        % the over-rotation must also hold counting only 0.5 s bins where the
        % heading moves (>= move_min_deg): divergence that piles up while the
        % bar is stationary (fly still, or FicTrac stalled) can't make a lap
        k_h = a:n_half_s:b; dBh = diff(rad2deg(mu_u(k_h))); dFh = diff(rad2deg(hd(k_h)));
        mov = abs(dFh) >= move_min_deg; df_m = sum(dFh(mov)); db_m = sum(dBh(mov));
        if abs(df_m) < fly_min_deg || sign(db_m) ~= sign(df) || abs(db_m) <= abs(df_m), continue, end
        hits(end+1,:) = [xb(a), xb(b), db, df, exc, median(rho(a:b)), w_s]; %#ok<AGROW>
    end
    end
    % overlapping hits (any window length) -> one episode, represented by its max-excess window
    while ~isempty(hits)
        [~,k] = max(hits(:,5)); e = hits(k,:);
        ev(end+1,:) = [i, e]; %#ok<AGROW>
        hits(hits(:,1) < e(2) & hits(:,2) > e(1),:) = [];
    end
end
E = array2table(ev,'VariableNames',{'trial','t0','t1','dbump','dfly','excess','rho','win_s'});
E.geno = IG.genotype(E.trial); E.fly = IG.fly_num(E.trial); E.gain = g_fly(E.fly);
[E.date,E.fly_n,E.trial_n] = cellfun(@meta_labels,{all_data(E.trial).meta}','UniformOutput',false);
fprintf('\n%s episodes (bump >= %g deg past the fly, fly >= %g deg, median rho >= %g, windows %s s):\n',cond_lbl,lap_min_excess_deg,fly_min_deg,rho_min,mat2str(lap_win_s));
fprintf('%-13s %8s %10s %14s %16s\n','genotype','episodes','per 10 min','flies with >=1','episodes/fly');
for gi = 1:numel(scan_geno)
    tr_g = scan_trials(strcmp(IG.genotype(scan_trials),scan_geno{gi}));
    mins = sum(arrayfun(@(i) all_data(i).ft.xf(end)-all_data(i).ft.xf(1),tr_g))/60;
    s = strcmp(E.geno,scan_geno{gi}); fl = unique(IG.fly_num(tr_g));
    n_f = accumarray(E.fly(s),1,[max(IG.fly_num),1]); n_f = n_f(fl);
    fprintf('%-13s %8d %10.2f %8d of %3d %9.2f (median %g)\n',scan_geno{gi},sum(s),10*sum(s)/mins, ...
        sum(n_f > 0),numel(fl),mean(n_f),median(n_f));
end
for gi = 1:numel(scan_geno) % which window length represents each episode
    s = strcmp(E.geno,scan_geno{gi});
    fprintf('  %-13s by window: %s\n',scan_geno{gi},strjoin(arrayfun(@(w) sprintf('%g s: %d',w,sum(s & E.win_s==w)),lap_win_s,'UniformOutput',false),', '));
end
% LPsP>TH, high- vs low-gain flies
tr_th = scan_trials(strcmp(IG.genotype(scan_trials),'lpsp>th')); fl = unique(IG.fly_num(tr_th)); fl = fl(~isnan(g_fly(fl)));
n_f = accumarray(E.fly(strcmp(E.geno,'lpsp>th')),1,[max(IG.fly_num),1]); n_f = n_f(fl); hi = g_fly(fl) > gain_split;
fprintf('LPsP>TH high gain: %d of %d flies with >=1 episode, %.2f episodes/fly; low gain: %d of %d, %.2f/fly; ranksum p = %.3g\n', ...
    sum(n_f(hi) > 0),sum(hi),mean(n_f(hi)),sum(n_f(~hi) > 0),sum(~hi),mean(n_f(~hi)),ranksum(n_f(hi),n_f(~hi)));
E = sortrows(E,'excess','descend');
save(fullfile(data_dir,sprintf('lpsp_rnai_lap_episodes_%s.mat',strrep(cond_lbl,' ','_'))),'E','lap_min_excess_deg','lap_win_s','light_cond'); % for fig_rnai_th_lap_onset.m
% disp(E(:,{'geno','date','fly_n','trial_n','t0','t1','dbump','dfly','rho','gain'})) % full episode list

%% pick one episode per main-genotype fly: the named example first (TH set), then by excess
Eth = E(strcmp(E.geno,main_geno),:);
is_named = strcmp(Eth.date,named_example{1}) & strcmp(Eth.fly_n,named_example{2}) & strcmp(Eth.trial_n,named_example{3});
if any(is_named), Eth = [Eth(find(is_named,1),:); Eth(~is_named,:)]; end
[~,first] = unique(Eth.fly,'stable'); pick = Eth(first(1:min(n_show,numel(first))),:);

%% figure: main-genotype examples (heatmap with heading / bump overlay above, unwrapped rotation below)
if isempty(pick), disp(['no ' main_name ' episodes']), else
plot_episodes(pick,all_data,gain_split,pad_s,any(is_named), ...
    {sprintf(['laps in ' main_name ' (' cond_lbl '): bump >= %g%s past a >= %g%s fly turn within %s s windows, median \\rho >= %g, one per fly'],lap_min_excess_deg,char(176),fly_min_deg,char(176),mat2str(lap_win_s),rho_min), ...
    'top: PB activity (z-scored) with bump (red) and fly heading (blue) on each half; bottom: unwrapped rotation, dotted every 360\circ; shaded / dotted = detected window'})
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_examples.png',set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_examples.pdf',set_tag)),'ContentType','image','Resolution',300)
end

%% figure: control examples, same detector, one episode per control fly (largest excess first)
Ec = E(ismember(E.geno,ctl_geno),:);
[~,first] = unique(Ec.fly,'stable'); pick_ctl = Ec(first(1:min(n_show,numel(first))),:);
if isempty(pick_ctl)
    disp('no control episodes')
else
    plot_episodes(pick_ctl,all_data,gain_split,pad_s,false, ...
        {sprintf('laps in control flies (%s, %s), same detector, one per fly, largest over-rotation first',strjoin(ctl_geno,', '),cond_lbl), ...
        'top: PB activity (z-scored) with bump (red) and fly heading (blue) on each half; bottom: unwrapped rotation, dotted every 360\circ; shaded / dotted = detected window'})
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_examples_controls.png',set_tag)),'Resolution',200)
    exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_examples_controls.pdf',set_tag)),'ContentType','image','Resolution',300)
end

%% figure: summary -- an example lap, laps per fly, fraction of flies with a lap, stats
% (a) one lap: PB heatmap with bump / cue (bar) position overlay, and
%     unwrapped rotation of bump and cue (closed loop: cue = 0.8 x heading;
%     detection itself uses heading).
% (b) laps per 10 min of closed loop, one dot per fly (rate, since flies have
%     1-4 closed-loop trials); (c) fraction of flies with >= 1 lap. Groups:
%     the TH set, then LPsP>TH split into high- / low-gain flies (gain_split).
%     p: ranksum on per-fly rates (b), Fisher exact on flies with >= 1 lap (c),
%     LPsP>TH vs each control and high vs low gain.
% groups (grp_def, grp_col, grp_x, cmp) come from geno_set in the parameters
n_fly_laps = accumarray(E.fly,1,[max(IG.fly_num),1]);
fly_min = accumarray(IG.fly_num(scan_trials),arrayfun(@(i) all_data(i).ft.xf(end)-all_data(i).ft.xf(1),scan_trials)/60,[max(IG.fly_num),1]);
G = struct('rate',{},'has',{},'n',{});
for gi = 1:size(grp_def,1)
    fl = unique(IG.fly_num(scan_trials(strcmp(IG.genotype(scan_trials),grp_def{gi,2}))));
    if ~isempty(grp_def{gi,3}), fl = fl(~isnan(g_fly(fl)) & (g_fly(fl) > gain_split) == grp_def{gi,3}); end
    G(gi).rate = 10*n_fly_laps(fl)./fly_min(fl); G(gi).has = n_fly_laps(fl) > 0; G(gi).n = numel(fl);
end
p_rate = nan(size(cmp,1),1); p_frac = p_rate;
for c = 1:size(cmp,1)
    a = G(cmp(c,1)); b = G(cmp(c,2));
    p_rate(c) = ranksum(a.rate,b.rate);
    [~,p_frac(c)] = fishertest([sum(a.has),sum(~a.has); sum(b.has),sum(~b.has)]);
    fprintf('%s vs %s: laps / 10 min ranksum p = %.2g; flies with >= 1 lap Fisher p = %.2g\n',grp_def{cmp(c,1),1},grp_def{cmp(c,2),1},p_rate(c),p_frac(c));
end
k = find(strcmp(E.date,summary_example{1}) & strcmp(E.fly_n,summary_example{2}) & strcmp(E.trial_n,summary_example{3}) & strcmp(E.geno,main_geno));
if isempty(k), k = find(strcmp(E.geno,main_geno)); end % example not in this set: the set's largest episode
figure('Color','w','Position',[40,40,1500,760]); clf
tl = tiledlayout(2,4,'TileSpacing','compact','Padding','compact');
if isempty(k)
    nexttile(tl,1,[2,2]); axis off; text(0.5,0.5,['no ' main_name ' laps'],'HorizontalAlignment','center','FontSize',12)
else
[~,j] = max(E.excess(k)); ex = E(k(j),:);
t0 = ex.t0; win = [-pad_s, ex.t1-t0+pad_s]; S = ex_snippet(all_data(ex.trial),t0,win);
nexttile(tl,1,[1,2]); hold on
n_half = size(S.f_z,1)/2;
imagesc(S.xb,1:2*n_half,S.f_z); clim([-1.5,3]); colormap(gca,flipud(gray))
if strcmp(light_cond,'dark'), ref = S.fly; ref_abs = S.fly_abs; ref_lbl = 'fly heading'; else, ref = S.cue; ref_abs = S.cue_abs; ref_lbl = 'cue (bar) position'; end
for h = 1:2 % bump and reference (bar in closed loop, heading in the dark) on each PB half, started at that half's bump at t0
    [tb,rb] = angle_to_row(S.xb,S.bump_abs(:,h),S.alpha); [tc,rc] = angle_to_row(S.xf,ref_abs(:,h),S.alpha);
    plot(tb,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',3); plot(tc,rc+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',3)
end
yline(n_half+0.5,'k'); xline([0,ex.t1-t0],'k:'); xlim(win); ylim([0.5,2*n_half+0.5]); yticks([]); set(gca,'YDir','reverse')
ylabel('PB glomeruli'); title(sprintf('(a) example lap: %s, %s fly %s trial %s (gain %.2f)',main_name,ex.date{1},ex.fly_n{1},ex.trial_n{1},ex.gain),'FontWeight','normal')
nexttile(tl,5,[1,2]); hold on
patch([0,ex.t1-t0,ex.t1-t0,0],[-1e4,-1e4,1e4,1e4],[.93,.93,.93],'EdgeColor','none','HandleVisibility','off')
plot(S.xf,ref,'Color',[0,0.45,0.74],'LineWidth',1.5,'DisplayName',ref_lbl)
plot(S.xb,S.bump,'Color',[0.85,0.2,0.1],'LineWidth',1.5,'DisplayName','bump (PVA)')
yl = [min([ref;S.bump]),max([ref;S.bump])]; ylim(yl + [-1,1]*0.08*diff(yl)); xlim(win)
if ex.dbump > 0, yt = 0:360:yl(2); else, yt = -fliplr(0:360:-yl(1)); end
for y = yt, yline(y,':','Color',[.6,.6,.6],'HandleVisibility','off'), end
xlabel(sprintf('time from %.0f s (s)',t0)); ylabel('rotation (deg)'); legend('Location','best')
dref = interp1(S.xf,ref,ex.t1-t0) - interp1(S.xf,ref,0);
title(sprintf('bump %+.0f%s vs %s %+.0f%s in the detected window (shaded); dotted every 360%s',ex.dbump,char(176),strtok(ref_lbl,' '),dref,char(176),char(176)),'FontWeight','normal')
end
nexttile(tl,3,[2,1]); hold on
for gi = 1:numel(G)
    y = G(gi).rate;
    scatter(grp_x(gi) + (rand(size(y))-.5)*.4,y,22,grp_col(gi,:),'filled','MarkerFaceAlpha',.55)
    errorbar(grp_x(gi),mean(y),std(y)/sqrt(numel(y)),'o','Color','k','MarkerFaceColor',grp_col(gi,:),'MarkerSize',8,'LineWidth',1.2)
end
add_p(gca,grp_x,cmp,p_rate)
xticks(grp_x); xticklabels(arrayfun(@(gi) sprintf('%s (n=%d)',grp_def{gi,1},G(gi).n),1:numel(G),'UniformOutput',false)); xtickangle(35)
xlim([0.4,max(grp_x)+0.6]); ylabel(['laps per 10 min of ' cond_lbl ' (per fly)']); title('(b) lap rate, dots = flies, mean \pm SEM','FontWeight','normal')
nexttile(tl,4,[2,1]); hold on
for gi = 1:numel(G)
    bar(grp_x(gi),mean(G(gi).has),0.7,'FaceColor',grp_col(gi,:),'EdgeColor','none')
    text(grp_x(gi),mean(G(gi).has)+0.02,sprintf('%d/%d',sum(G(gi).has),G(gi).n),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9)
end
ylim([0,1]); add_p(gca,grp_x,cmp,p_frac); yticks(0:0.25:1)
xticks(grp_x); xticklabels(arrayfun(@(gi) grp_def{gi,1},1:numel(G),'UniformOutput',false)); xtickangle(35)
xlim([0.4,max(grp_x)+0.6]); ylabel('fraction of flies with >= 1 lap'); title('(c) flies with a lap','FontWeight','normal')
title(tl,{['lapping: the bump over-rotates the fly (' cond_lbl ')'], ...
    sprintf('lap: in a %s s window the fly turns >= %g%s and the bump goes >= %g%s further, same direction, steadily diverging, median \\rho >= %g', ...
    strjoin(arrayfun(@num2str,lap_win_s,'UniformOutput',false),' or '),fly_min_deg,char(176),lap_min_excess_deg,char(176),rho_min), ...
    'p: ranksum on per-fly lap rates (b), Fisher exact on flies with >= 1 lap (c)'},'FontSize',11)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_summary.png',set_tag)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lapping_summary.pdf',set_tag)),'ContentType','image','Resolution',300)

%% functions
function add_p(ax,x,cmp,p)
    % significance brackets above the data, stacked
    yl = ylim(ax); top = yl(2); h = 0.06*diff(yl);
    for c = 1:size(cmp,1)
        y = top + (c-0.4)*h*1.6;
        plot(ax,x(cmp(c,[1,1,2,2])),[y-h*0.4,y,y,y-h*0.4],'k','LineWidth',0.8,'HandleVisibility','off')
        if p(c) < 0.001, ps = sprintf('p = %.1e',p(c)); else, ps = sprintf('p = %.3f',p(c)); end
        text(ax,mean(x(cmp(c,:))),y,ps,'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8)
    end
    ylim(ax,[yl(1), top + (size(cmp,1)+0.6)*h*1.6])
end

function plot_episodes(pick,all_data,gain_split,pad_s,named_first,ttl_lines)
    % one column per episode, 4 per row: heatmap with bump / heading overlay, then unwrapped rotation
    nc = 4; nr = ceil(height(pick)/nc);
    figure('Color','w','Position',[40,40,430*nc,400*nr]); clf
    tl = tiledlayout(2*nr,nc,'TileSpacing','compact','Padding','compact');
    for k = 1:height(pick)
        i = pick.trial(k); t0 = pick.t0(k); win = [-pad_s, pick.t1(k)-t0+pad_s];
        S = ex_snippet(all_data(i),t0,win);
        r = floor((k-1)/nc); c = mod(k-1,nc)+1;
        nexttile(tl,2*r*nc+c); hold on
        n_half = size(S.f_z,1)/2;
        imagesc(S.xb,1:2*n_half,S.f_z); clim([-1.5,3]); colormap(gca,flipud(gray))
        for h = 1:2
            [tb,rb] = angle_to_row(S.xb,S.bump_abs(:,h),S.alpha); [tf,rf] = angle_to_row(S.xf,S.fly_abs(:,h),S.alpha);
            plot(tb,rb+(h-1)*n_half,'Color',[0.85,0.2,0.1],'LineWidth',1.5); plot(tf,rf+(h-1)*n_half,'Color',[0,0.45,0.74],'LineWidth',1.5)
        end
        yline(n_half+0.5,'k'); xline([0,pick.t1(k)-t0],'k:')
        xlim(win); ylim([0.5,2*n_half+0.5]); yticks([]); set(gca,'YDir','reverse')
        if isnan(pick.gain(k)), g_lbl = 'gain n/a'; elseif pick.gain(k) > gain_split, g_lbl = sprintf('gain %.2f, high',pick.gain(k)); else, g_lbl = sprintf('gain %.2f, low',pick.gain(k)); end
        ttl = sprintf('%s | %s fly %s trial %s (%s)',pick.geno{k},pick.date{k},pick.fly_n{k},pick.trial_n{k},g_lbl);
        if named_first && k == 1, ttl = ['named example: ' ttl]; end
        title(ttl,'FontSize',8,'FontWeight','normal','Interpreter','none')
        if c == 1, ylabel('PB glomeruli'), end
        nexttile(tl,(2*r+1)*nc+c); hold on
        patch([0,pick.t1(k)-t0,pick.t1(k)-t0,0],[-1e4,-1e4,1e4,1e4],[.93,.93,.93],'EdgeColor','none','HandleVisibility','off')
        plot(S.xf,S.fly,'Color',[0,0.45,0.74],'LineWidth',1.5,'DisplayName','fly heading')
        plot(S.xb,S.bump,'Color',[0.85,0.2,0.1],'LineWidth',1.5,'DisplayName','bump (PVA)')
        yl = [min([S.fly;S.bump]),max([S.fly;S.bump])]; ylim(yl + [-1,1]*0.08*diff(yl)); xlim(win)
        if pick.dbump(k) > 0, yt = 0:360:yl(2); else, yt = -fliplr(0:360:-yl(1)); end
        for y = yt, yline(y,':','Color',[.6,.6,.6],'HandleVisibility','off'), end
        title(sprintf('window: bump %+.0f%s, fly %+.0f%s',pick.dbump(k),char(176),pick.dfly(k),char(176)),'FontSize',9,'FontWeight','normal')
        xlabel(sprintf('time from %.0f s (s)',t0)); if c == 1, ylabel('rotation (deg)'), end
        if k == 1, legend('Location','best'), end
    end
    title(tl,ttl_lines)
end

function [xb,mu_raw,mu_u,rho,hd] = trial_bump_heading(tr)
    % PVA (raw unwrapped and 0.2 s smoothed) and unwrapped heading on the imaging timebase
    xf = tr.ft.xf(:); nb = size(tr.im.f,2);
    if numel(tr.ft.xb) == nb, xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),nb)'; end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2))./std(f_s,0,2);
    [x_,y_] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_,2),mean(y_,2));
    mu_raw = unwrap(mu);
    mu_u = smoothdata(mu_raw,'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = interp1(xf,unwrap(fillmissing(tr.ft.heading(:),'linear','EndValues','nearest')),xb);
end

function [d,f,t] = meta_labels(meta)
    % date, fly number and trial number from the ScanImage path
    tok = regexp(meta,'\\(\d{8})\\fly\s*(\d+)[^\\]*\\\d{8}-(\d+)_','tokens','once'); % fly folder may carry a suffix ('fly 6 _ empty _ vglut')
    if isempty(tok), d = ''; f = ''; t = ''; else, [d,f,t] = deal(tok{:}); end
end

function S = ex_snippet(tr,t0,win_s)
    % as fig_rnai_turns_examples.m: z-scored PB activity, bump and heading in
    % t0+win_s, times relative to t0, angles in deg relative to t0, plus
    % per-PB-half absolute PVA / heading for the heatmap overlay
    S = [];
    xf = tr.ft.xf(:);
    if isfield(tr.ft,'xb') && numel(tr.ft.xb) == size(tr.im.f,2), xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),size(tr.im.f,2))'; end
    win_s = [max(win_s(1),max(xf(1),xb(1))-t0), min(win_s(2),min(xf(end),xb(end))-t0)];
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2)) ./ std(f_s,0,2);
    [x_tmp,y_tmp] = pol2cart(tr.im.alpha(:)',f_z');
    [mu,rho] = cart2pol(mean(x_tmp,2),mean(y_tmp,2)); mu = unwrap(mu);
    mu = smoothdata(mu,'gaussian',max(1,round(0.2/median(diff(xb)))));
    hd = fillmissing(tr.ft.heading(:),'linear','EndValues','nearest'); hd = unwrap(hd);
    in_im = xb >= t0+win_s(1) & xb <= t0+win_s(2);
    in_ft = xf >= t0+win_s(1) & xf <= t0+win_s(2);
    S.xb = xb(in_im)-t0; S.f_z = f_z(:,in_im); S.bump = rad2deg(mu(in_im)-interp1(xb,mu,t0)); S.rho = rho(in_im);
    S.xf = xf(in_ft)-t0; S.fly = rad2deg(hd(in_ft)-interp1(xf,hd,t0));
    cu = unwrap_cue(tr.ft.cue); % bar position (closed loop: 0.8 x heading)
    S.cue = rad2deg(cu(in_ft)-interp1(xf,cu,t0));
    n_half = size(tr.im.f,1)/2;
    S.alpha = tr.im.alpha(1:n_half); S.alpha = S.alpha(:)';
    S.bump_abs = nan(sum(in_im),2); S.fly_abs = nan(sum(in_ft),2); S.cue_abs = S.fly_abs;
    for h = 1:2
        rows = (h-1)*n_half + (1:n_half);
        [x_h,y_h] = pol2cart(tr.im.alpha(rows),f_z(rows,:)');
        mu_h = smoothdata(unwrap(cart2pol(mean(x_h,2),mean(y_h,2))),'gaussian',max(1,round(0.2/median(diff(xb)))));
        S.bump_abs(:,h) = mu_h(in_im); S.fly_abs(:,h) = interp1(xb,mu_h,t0) + deg2rad(S.fly);
        S.cue_abs(:,h) = interp1(xb,mu_h,t0) + deg2rad(S.cue);
    end
end

function [t,row] = angle_to_row(t,ang,alpha)
    % PVA angle -> fractional row within one PB half (as fig_rnai_turns_examples.m)
    n = numel(alpha); d = 2*pi/n;
    row = mod(ang(:) - alpha(1) + d/2, 2*pi)/d + 0.5;
    jump = [false; abs(diff(row)) > n/2];
    t = t(:); k = find(jump);
    for j = numel(k):-1:1
        t = [t(1:k(j)-1); NaN; t(k(j):end)]; row = [row(1:k(j)-1); NaN; row(k(j):end)];
    end
end

function cu = unwrap_cue(cue)
    % -ft.cue unwrapped over the valid samples, THEN gaps filled (filling first
    % interpolates across the +-pi wrap the long way round: a fake 360 deg jump)
    cu = -cue(:); ok = ~isnan(cu);
    cu(ok) = unwrap(cu(ok)); cu = fillmissing(cu,'linear','EndValues','nearest');
end
