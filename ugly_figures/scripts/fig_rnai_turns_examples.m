%% fig_rnai_turns_examples
% Example turns per genotype -- z-scored PB heatmap with heading and bump
% overlaid, and heading vs. bump traces below -- for one genotype set and
% one 15 deg turn-size bin (centers 15, 30, ..., same bins as the extended
% size overlays in fig_rnai_turns_overlays.m). Generalizes the 30 deg
% TH-set examples in fig_rnai_turns_quiet_cl.m. Closed loop; turn_set
% picks quiet-baseline turns ('quiet') or all turns ('all'), matching the
% overlay figure the examples should represent. The pick per genotype is
% the turn whose aligned fly and bump traces (-1..3 s) are closest to that
% genotype's group-mean traces for the bin (same turn set), among turns
% with a solid bump; with 'quiet', candidates must also be clean of other
% turns > 15 deg over -2..4 s. With 'all' that cleanliness rule is off:
% the all-turns mean itself includes the preceding counter-turn and the
% continued turning after the turn, so excluding those would exclude the
% typical turn. Set geno_set / turn_set / ex_bin_deg below; ex_override
% pins specific turns (indices from the printed ranking).
% Needs the turn cache and all_data (loaded if absent). Run from the repo root.

%% parameters + load turn cache
geno_set = 'vglut';
turn_set = 'all';  % 'quiet' or 'all'
switch geno_set
    case 'vglut'
        plot_geno  = {'lpsp>vglut','empty>vglut','lpsp>mcherry'};
        plot_label = {'LPsP > VGlut-RNAi','Empty > VGlut-RNAi','LPsP > mCherry'};
    case 'th'
        plot_geno  = {'lpsp>th','empty>th','lpsp>mcherry'};
        plot_label = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'};
end
cache_file = fullfile('.data','lpsp_rnai_turns_cache.mat');
load(cache_file,'turns','turn_quiet','turn_fly','turn_dark', ...
    'turn_t','min_turns_per_bin','fly_geno','n_flies','geno_keep','geno_colors');
n_g = numel(plot_geno);
sel = turn_dark==0;
if strcmp(turn_set,'quiet'), sel = sel & turn_quiet; set_lbl = 'quiet-baseline'; set_tag = '';
else, set_lbl = 'all'; set_tag = '_allturns'; end
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

%% example turns: candidates -- every turn of the set in the bin, per genotype
% The pick is made on the raw data below, and (with no fixed window) the
% plotted window is then widened as far as all three picks allow.
ex_bin_deg     = 135;         % turn-size bin center (deg); bin = center +/- 7.5 deg
ex_core_s      = [-2, 4];     % s around onset: scoring window (and, for 'quiet', must be clean)
ex_max_win_s   = [-7.5, 7.5]; % s, the widest shared window considered
iso_max_deg    = 15;          % "other turns" = detected turns bigger than this
ex_rho_p10_min = 0.3;         % bump gate: 10th percentile of PVA strength over the core window
ex_override    = nan(1,3);    % turn index per genotype (from the printed ranking) to override the automatic pick; NaN = automatic
ex_fixed_win_s = [-2, 5];     % s; plot this window instead of the derived shared one (empty = use the shared window). Other turns may then fall inside it.

turn_trial = [turns.trial]'; turn_on = [turns.onset_s]'; turn_size = [turns.size_deg]'; turn_off = turn_on + [turns.dur_s]';
in_bin = abs(turn_size - ex_bin_deg) < 7.5;
ex_b = 1; bin_lbl = sprintf('%g%s',ex_bin_deg,char(176));
% per-fly mean traces for this bin (group means for the pick below), as in the overlay figures
FM = nan(n_flies,1,numel(turn_t)); BM = FM;
for f = unique(turn_fly(sel & in_bin))'
    idx = find(sel & turn_fly==f & in_bin);
    if numel(idx) < min_turns_per_bin, continue, end
    FM(f,1,:) = mean(cat(2,turns(idx).fly_trace),2,'omitnan');
    BM(f,1,:) = mean(cat(2,turns(idx).bump_trace),2,'omitnan');
end
ex_cand = cell(1,n_g);
for gi = 1:n_g
    flies_g = strcmp(fly_geno,plot_geno{gi});
    idx = find(sel & in_bin & flies_g(turn_fly));
    clean = true(size(idx));
    if strcmp(turn_set,'quiet')
        for k = 1:numel(idx)
            i = idx(k);
            clean(k) = ~any(turn_trial==turn_trial(i) & (1:numel(turns))'~=i & turn_size > iso_max_deg & ...
                turn_off > turn_on(i)+ex_core_s(1) & turn_on < turn_on(i)+ex_core_s(2));
        end
    end
    ex_cand{gi} = idx(clean);
    fprintf('%s: %d %s-turn candidates in the %s bin (%d flies in the group mean)\n',plot_geno{gi},sum(clean),set_lbl,bin_lbl, ...
        sum(~all(isnan(reshape(FM(flies_g,1,:),sum(flies_g),[])),2)));
end

%% example turns: load raw trials (~3.7 GB, slow -- skipped if all_data is already in the workspace)
if ~exist('all_data','var')
    load(cache_file,'source_file');
    tmp = load(fullfile('.data',source_file),'all_data'); all_data = tmp.all_data(:); clear tmp
end

%% example turns: pick one per genotype, then the widest window all three allow
% Score on the core window. Gate: a solid bump (10th percentile of PVA
% strength >= ex_rho_p10_min). Then smallest
%   RMS(fly trace - group-mean fly trace) + RMS(bump trace - group-mean bump trace)
% over the aligned -1..3 s traces, against this genotype's group-mean
% traces for the bin (as in the size-binned overlay figures):
% the turn that looks most like
% the genotype's mean. (Also tried: step sizes (1..3 s minus -2..0 s) vs.
% targets plus pre/post SDs -> thrown off by pre-onset wiggles; bump step
% vs. g x fly step; free-gain fit residual -> favored bumps that barely
% moved; bump-heading correlation -> favored bumps that overshot
% several-fold; whole-window RMS vs. g x heading -> favored windows where
% nothing moved.) The top few are printed; set ex_override(gi) to a turn
% index from that list to show a different one.
% Shared window: for each pick, how far before/after onset the trace
% stays free of other turns > iso_max_deg (and inside the recording);
% the plotted window is the tightest of the three, capped at ex_max_win_s.
step = @(t,y) mean(y(t>=1 & t<=3),'omitnan') - mean(y(t>=-2 & t<=0),'omitnan');
ex_pick = nan(1,n_g); free_pre = nan(1,n_g); free_post = nan(1,n_g);
for gi = 1:n_g
    flies_g = strcmp(fly_geno,plot_geno{gi});
    f_mean = mean(reshape(FM(flies_g,ex_b,:),sum(flies_g),[]),1,'omitnan')'; % group-mean traces, as in the overlay figure
    b_mean = mean(reshape(BM(flies_g,ex_b,:),sum(flies_g),[]),1,'omitnan')';
    cand = ex_cand{gi}; err = nan(size(cand)); info = cell(size(cand));
    for k = 1:numel(cand)
        S_c = ex_snippet(all_data(turns(cand(k)).trial),turns(cand(k)).onset_s,ex_core_s);
        if isempty(S_c) || prctile(S_c.rho,10) < ex_rho_p10_min, continue, end
        dB = turns(cand(k)).dir*step(S_c.xb,S_c.bump); dF = turns(cand(k)).dir*step(S_c.xf,S_c.fly);
        err(k) = sqrt(mean((turns(cand(k)).fly_trace - f_mean).^2,'omitnan')) + sqrt(mean((turns(cand(k)).bump_trace - b_mean).^2,'omitnan'));
        info{k} = sprintf('bump step %.0f, fly step %.0f deg; rho p10 %.2f',dB,dF,prctile(S_c.rho,10));
    end
    [e_sorted,o] = sort(err); o = o(~isnan(e_sorted));
    fprintf('%s (group mean at 3 s: fly %.0f, bump %.0f deg), %d candidates pass the bump gate; best:\n',plot_geno{gi},f_mean(end),b_mean(end),numel(o));
    for k = o(1:min(5,numel(o)))'
        fprintf('  turn %6d: fly %3d, trial %3d, %2.0f deg at %4.0f s; score %5.1f; %s\n',cand(k),turns(cand(k)).fly,turns(cand(k)).trial, ...
            turns(cand(k)).size_deg,turns(cand(k)).onset_s,err(k),info{k});
    end
    if ~isnan(ex_override(gi)), i = ex_override(gi); else, i = cand(o(1)); end
    ex_pick(gi) = i;
    % clean extent around the pick
    tr = all_data(turns(i).trial); xf = tr.ft.xf(:);
    if isfield(tr.ft,'xb') && numel(tr.ft.xb) == size(tr.im.f,2), xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),size(tr.im.f,2))'; end
    others = turn_trial==turn_trial(i) & (1:numel(turns))'~=i & turn_size > iso_max_deg;
    prev_end = max([max(xf(1),xb(1)); turn_off(others & turn_off <= turn_on(i))]);
    next_on  = min([min(xf(end),xb(end)); turn_on(others & turn_on >= turn_on(i))]);
    free_pre(gi) = turn_on(i) - prev_end; free_post(gi) = next_on - turn_on(i);
    fprintf('  -> picked turn %d: clean from %.1f s before to %.1f s after onset\n',i,free_pre(gi),free_post(gi));
end
ex_win_s = [-floor(min([free_pre,-ex_max_win_s(1)])*2)/2, floor(min([free_post,ex_max_win_s(2)])*2)/2]; % rounded inward to 0.5 s
fprintf('shared window: %.1f..%.1f s around onset (%.1f s)\n',ex_win_s,diff(ex_win_s));
if ~isempty(ex_fixed_win_s), ex_win_s = ex_fixed_win_s; fprintf('using fixed window %g..%g s instead\n',ex_win_s); end

%% example turns: PB heatmap (z-scored) + heading/bump, one column per genotype
figure('Color','w','Position',[60,60,420*n_g,560]); clf
tl3 = tiledlayout(2,n_g,'TileSpacing','compact','Padding','compact');
for gi = 1:n_g
    c = geno_colors(strcmp(geno_keep,plot_geno{gi}),:);
    i = ex_pick(gi);
    S = ex_snippet(all_data(turns(i).trial),turns(i).onset_s,ex_win_s);
    t0 = turns(i).onset_s; s = turns(i).dir; % sign-flip like the overlays so the turn is positive

    nexttile(tl3,gi)
    n_half = size(S.f_z,1)/2;
    imagesc(S.xb,1:2*n_half,S.f_z); hold on
    % white -> genotype color -> darker shade, so each heatmap matches its genotype's color elsewhere
    cmap = interp1([0;0.65;1],[1,1,1; c; c*0.45],linspace(0,1,256)');
    clim([-1.5,3]); colormap(gca,cmap)
    % heading (black) and bump PVA (white, black-edged) mapped onto glomerulus rows, each PB half with its own PVA
    h_ov = gobjects(1,2);
    for h = 1:2
        [tb,rb] = angle_to_row(S.xb,S.bump_abs(:,h),S.alpha);
        [tf,rf] = angle_to_row(S.xf,S.fly_abs(:,h),S.alpha);
        plot(tb,rb+(h-1)*n_half,'k','LineWidth',3,'HandleVisibility','off')
        h_ov(2) = plot(tb,rb+(h-1)*n_half,'w','LineWidth',1.5);
        h_ov(1) = plot(tf,rf+(h-1)*n_half,'k','LineWidth',1.5);
    end
    yline(n_half+0.5,'k','LineWidth',1); xline(0,':k','LineWidth',1)
    xlim(ex_win_s); ylim([0.5,2*n_half+0.5]); yticks([]);
    title(sprintf('%s\nfly %d, trial %d, %.0f%s turn at %.0f s',plot_label{gi},turns(i).fly,turns(i).trial,turns(i).size_deg,char(176),t0),'FontWeight','normal')
    if gi == 1, ylabel('PB glomeruli (both halves)'), legend(h_ov,{'fly heading','bump (PVA)'},'Location','southwest','Color',[.9,.9,.9]), end
    cb = colorbar; if gi == n_g, cb.Label.String = 'z-scored \DeltaF'; end

    nexttile(tl3,n_g+gi); hold on
    yl = [min([s*S.fly;s*S.bump]), max([s*S.fly;s*S.bump])]; yl = yl + [-1,1]*0.1*diff(yl);
    patch([0,turns(i).dur_s,turns(i).dur_s,0],yl([1,1,2,2]),[.9,.9,.9],'EdgeColor','none','HandleVisibility','off')
    plot(S.xf,s*S.fly,'k','LineWidth',1.5,'DisplayName','fly heading')
    plot(S.xb,s*S.bump,'Color',c,'LineWidth',1.5,'DisplayName','bump (PVA)')
    xline(0,':k','HandleVisibility','off'); xlim(ex_win_s); ylim(yl)
    xlabel('time from turn onset (s)')
    if gi == 1, ylabel('rotation from onset (deg)'), legend('Location','northwest'), end
end
title(tl3,sprintf('closed loop, %s %g%s turns: example per genotype closest to its group-mean traces (shaded = detected turn)',set_lbl,ex_bin_deg,char(176)))


%% export examples
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_examples_%gdeg%s.pdf',geno_set,ex_bin_deg,set_tag)),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_turns_%s_examples_%gdeg%s.png',geno_set,ex_bin_deg,set_tag)),'Resolution',200)

%% functions
function S = ex_snippet(tr,t0,win_s)
    % z-scored PB activity, bump and fly heading for the window t0+win_s,
    % times relative to t0 and angles in deg relative to their value at t0.
    % Same activity/bump estimate as lpsp_rnai_turns_claude.m (movmean 3
    % frames, z-score over the whole trial, PVA, 0.2 s Gaussian on the
    % unwrapped mu). Empty if the window runs outside fictrac or imaging.
    S = [];
    xf = tr.ft.xf(:);
    if isfield(tr.ft,'xb') && numel(tr.ft.xb) == size(tr.im.f,2), xb = tr.ft.xb(:); else, xb = linspace(xf(1),xf(end),size(tr.im.f,2))'; end
    if t0+win_s(1) < max(xf(1),xb(1)) || t0+win_s(2) > min(xf(end),xb(end)), return, end
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2)) ./ std(f_s,0,2);
    [x_tmp,y_tmp] = pol2cart(tr.im.alpha(:)',f_z');
    [mu,rho] = cart2pol(mean(x_tmp,2),mean(y_tmp,2)); mu = unwrap(mu);
    mu = smoothdata(mu,'gaussian',max(1,round(0.2/median(diff(xb)))));
    if isfield(tr.ft,'heading') && numel(tr.ft.heading) == numel(xf), hd = fillmissing(tr.ft.heading(:),'linear','EndValues','nearest');
    else, r = tr.ft.r_speed(:); r(isnan(r)) = 0; hd = cumtrapz(xf,r); end
    in_im = xb >= t0+win_s(1) & xb <= t0+win_s(2);
    in_ft = xf >= t0+win_s(1) & xf <= t0+win_s(2);
    S.xb = xb(in_im)-t0; S.f_z = f_z(:,in_im); S.bump = rad2deg(mu(in_im)-interp1(xb,mu,t0)); S.rho = rho(in_im);
    S.xf = xf(in_ft)-t0; S.fly = rad2deg(hd(in_ft)-interp1(xf,hd,t0));
    % absolute angles (rad) for the heatmap overlay, PER PB HALF: each
    % half's own PVA (same smoothing), and heading started at that half's
    % bump position at onset (closed loop: bump moves with the same sign as
    % heading, offset arbitrary). Per half because the combined PVA above
    % averages both halves, and in some flies the two halves' bumps sit at
    % different alpha -- the combined PVA then lands between the bands.
    n_half = size(tr.im.f,1)/2;
    S.alpha = tr.im.alpha(1:n_half); S.alpha = S.alpha(:)';
    S.bump_abs = nan(sum(in_im),2); S.fly_abs = nan(sum(in_ft),2);
    for h = 1:2
        rows = (h-1)*n_half + (1:n_half);
        [x_h,y_h] = pol2cart(tr.im.alpha(rows),f_z(rows,:)');
        mu_h = smoothdata(unwrap(cart2pol(mean(x_h,2),mean(y_h,2))),'gaussian',max(1,round(0.2/median(diff(xb)))));
        S.bump_abs(:,h) = mu_h(in_im); S.fly_abs(:,h) = interp1(xb,mu_h,t0) + deg2rad(S.fly);
    end
end

function [t,row] = angle_to_row(t,ang,alpha)
    % PVA angle -> fractional row within one PB half, whose rows sit at
    % evenly spaced angles alpha (row k centered on alpha(k)); NaN-breaks
    % the line where it wraps from one edge of the half to the other
    n = numel(alpha); d = 2*pi/n;
    row = mod(ang(:) - alpha(1) + d/2, 2*pi)/d + 0.5;
    jump = [false; abs(diff(row)) > n/2];
    t = t(:); k = find(jump);
    for j = numel(k):-1:1
        t = [t(1:k(j)-1); NaN; t(k(j):end)]; row = [row(1:k(j)-1); NaN; row(k(j):end)];
    end
end

function plot_sem(t,X,c)
    % mean +/- SEM shading across rows of X, line in color c
    t = t(:)';
    m = mean(X,1,'omitnan'); s = std(X,0,1,'omitnan')./sqrt(max(1,sum(~isnan(X),1)));
    ok = ~isnan(m);
    if ~any(ok), return, end
    patch([t(ok),fliplr(t(ok))],[m(ok)+s(ok),fliplr(m(ok)-s(ok))],c,'EdgeColor','none','FaceAlpha',.2)
    plot(t,m,'Color',c,'LineWidth',1.5)
end
