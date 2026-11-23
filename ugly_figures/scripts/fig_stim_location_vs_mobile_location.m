%% fig_stim_location_vs_mobile_location
% Scratch: is the post-ATP increase in dark bump mobility SPATIALLY
% specific -- larger at the PB positions where the bump sat while LPsP was
% being stimulated? Walking dataset (lpsp_p2x2_walking_20260728_post0708).
%
% Two REFERENCE angles per fly, each run through the identical pipeline:
%  (A) theta_stim -- where the BUMP was during stimulation. In the fly's
%      first pulse trial (closed loop, 29 TTL pulses, 10 s apart), for every
%      CONFIRMED pulse (whole-PB atp.d median 0-3 s minus median -2..0 s >
%      3 baseline SDs), collect im.mu over 0..stim_win s after onset;
%      theta_stim = their circular mean, R = resultant length (in closed
%      loop the bump follows the cue, so R is low: the bump was spread out).
%  (B) theta_hot -- where the ATP LANDED (per explicit request). Average
%      the confirmed pulses' baseline-subtracted 0-3 s atp.d profile (32
%      wedges); keep the PB half (wedges 1:16 or 17:32) holding its peak;
%      theta_hot = the population-vector angle of that half's positive
%      profile over im.alpha -- the same wedge-angle frame im.mu is a PVA
%      in, so a bump at mu = theta_hot has its peak ON the ATP-injected
%      glomeruli. Fixed per fly (unlike A), so a sharper spatial reference.
% For each reference:
%  POSITION-RESOLVED MOBILITY. In the fly's dark trials before (pre) and
%    after (post) the pulse trial -- the same pools as Figure 3 -- take
%    every sample inside a walking bout (Figure 3 bout definition), bump
%    smoothed as in walk_trial_bouts (60-sample Gaussian on ft.xf). Bin each
%    sample by bump position RELATIVE to the reference (n_bin bins); per
%    bin, mobility = sum |d mu| / sum |r_speed|*dt (Figure 3's path-length
%    ratio, made local). A bin needs >= min_head rad of heading path in both
%    pre and post, else NaN. d_mob(bin) = post - pre.
%  TESTS. (a) mean d_mob vs relative position per genotype: a spatially
%    specific effect is a peak at 0 in lpsp only. (b) per fly, Spearman
%    across bins between d_mob and closeness to the reference,
%    cos(relative angle) (+ = more mobile near it); signrank vs 0 per
%    genotype, ranksum between. (c) near (2 bins around 0) minus opposite
%    (2 bins at +/-pi). empty>P2X2 runs the identical pipeline as the
%    control for any bias the alignment itself introduces.
%
% REQUIRES the Figure_3_claude.m workspace (walk_* variables, circ_stats on path).

%% check workspace
need = {'walk_data','walk_fly_ix','walk_fly_list','walk_trial_num','walk_has_pulse','walk_is_dark','walk_fly_group','walk_delta','gcamp_color', ...
        'walk_smooth_window','walk_turn_thresh','walk_max_gap_frames','walk_min_walking_frames'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');
wrap = @(x) angle(exp(1i*x));

%% parameters
n_bin = 8; edges = linspace(-pi,pi,n_bin+1); ctr = edges(1:end-1)+pi/n_bin;
stim_win = 1; base_win = [-2 0]; peak_win = [0 3]; peak_factor = 3;
min_head = pi; % rad of heading path per bin, in pre AND post

%% per fly: both reference angles
nF = numel(walk_fly_list);
ref_names = {'theta_stim','theta_hot'};
ref_label = {'bump position during stimulation (\theta_{stim})','ATP hotspot (\theta_{hot})'};
F = struct();
for a = 1:2
    F(a).occ = nan(nF,n_bin); F(a).d_mob = nan(nF,n_bin); F(a).r_cos = nan(nF,1); F(a).near_far = nan(nF,1);
end
R_stim = nan(nF,1); hot_minus_stim = nan(nF,1);
for f = 1:nF
    if isnan(walk_delta(f)); continue; end
    ti = find(walk_fly_ix==f); [~,o] = sort(walk_trial_num(ti)); ti = ti(o);
    p = find(walk_has_pulse(ti),1);
    pre  = ti(1:p-1); pre = pre(walk_is_dark(pre));
    post = ti(p+1:end); post = post(~walk_has_pulse(post) & walk_is_dark(post));

    % pulse trial: confirmed pulses -> bump samples during stim + atp profile
    tr = walk_data(ti(p)); xb = tr.ft.xb(:); xf = tr.ft.xf(:); mu_k = tr.im.mu(:);
    st = logical(tr.ft.stims(:)); t_on = xf(diff([false;st])==1);
    atp_m = mean(tr.atp.d,1)'; mu_s = []; prof = zeros(size(tr.atp.d,1),1);
    for k = 1:numel(t_on)
        bwin = xb>=t_on(k)+base_win(1) & xb<t_on(k)+base_win(2);
        q = xb>t_on(k)+peak_win(1) & xb<=t_on(k)+peak_win(2);
        if median(atp_m(q)) - median(atp_m(bwin)) <= peak_factor*std(atp_m(bwin)); continue; end
        mu_s = [mu_s; mu_k(xb>=t_on(k) & xb<=t_on(k)+stim_win)]; %#ok<AGROW>
        prof = prof + mean(tr.atp.d(:,q),2) - mean(tr.atp.d(:,bwin),2);
    end
    if numel(mu_s) < 50; continue; end
    theta_stim = circ_mean(mu_s); R_stim(f) = circ_r(mu_s);
    [~,pk] = max(prof); half = (1:16) + 16*(pk>16);
    w = max(prof(half),0); alpha = tr.im.alpha(:);
    theta_hot = angle(sum(w .* exp(1i*alpha(half))));
    hot_minus_stim(f) = wrap(theta_hot - theta_stim);
    thetas = [theta_stim, theta_hot];

    for a = 1:2
        th = thetas(a);
        F(a).occ(f,:) = histcounts(wrap(mu_s-th), edges, 'Normalization','probability');
        [dm_pre, dh_pre]   = binned_paths(walk_data(pre),  th, edges, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames);
        [dm_post, dh_post] = binned_paths(walk_data(post), th, edges, walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames);
        ok = dh_pre >= min_head & dh_post >= min_head;
        d = nan(1,n_bin); d(ok) = dm_post(ok)./dh_post(ok) - dm_pre(ok)./dh_pre(ok);
        F(a).d_mob(f,:) = d;
        if sum(ok) >= 5; F(a).r_cos(f) = corr(cos(ctr(ok))', d(ok)', 'type','Spearman'); end
        near = ceil(n_bin/2)+[0 1]; far = [1 n_bin];
        F(a).near_far(f) = mean(d(near),'omitnan') - mean(d(far),'omitnan');
    end
end

%% stats
geno_label = {'empty>P2X2','lpsp>P2X2'}; gcol = [.5 .5 .5; gcamp_color];
for gi = 1:2
    g = walk_fly_group==gi & ~isnan(R_stim);
    fprintf('%-11s bump-during-stim concentration R median %.2f [%.2f-%.2f]; theta_hot - theta_stim: |median| %.0f deg, circ R %.2f\n', ...
        geno_label{gi}, median(R_stim(g)), min(R_stim(g)), max(R_stim(g)), abs(circ_median(hot_minus_stim(g)))*180/pi, circ_r(hot_minus_stim(g)));
end
for a = 1:2
    fprintf('--- aligned to %s ---\n', ref_names{a});
    for gi = 1:2
        g = walk_fly_group==gi & ~isnan(F(a).r_cos);
        fprintf('%-11s n=%2d  Spearman(d_mob, cos rel. angle) median %+.2f, signrank p=%.3f | near - opposite median %+.2f, signrank p=%.3f\n', ...
            geno_label{gi}, sum(g), median(F(a).r_cos(g)), signrank(F(a).r_cos(g)), median(F(a).near_far(g),'omitnan'), signrank(F(a).near_far(g)));
    end
    g1 = walk_fly_group==1 & ~isnan(F(a).r_cos); g2 = walk_fly_group==2 & ~isnan(F(a).r_cos);
    fprintf('between genotypes: Spearman ranksum p=%.3f, near-opposite ranksum p=%.3f\n', ranksum(F(a).r_cos(g1),F(a).r_cos(g2)), ranksum(F(a).near_far(g1),F(a).near_far(g2)));
end

%% plot: one row per reference angle
fig = figure('color','w','Visible','off','Position',[50 50 1700 1000]);
tl = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
deg = ctr*180/pi; all_dmob_axes = gobjects(0);
for a = 1:2
    ax = nexttile(tl); hold(ax,'on')
    for gi = 1:2
        g = walk_fly_group==gi & ~isnan(F(a).r_cos); X = F(a).occ(g,:);
        errorbar(ax, deg, mean(X,1), std(X,[],1)/sqrt(sum(g)), '-o', 'Color', gcol(gi,:), 'MarkerFaceColor', gcol(gi,:), 'LineWidth', 1.5, 'DisplayName', geno_label{gi})
    end
    yline(ax,1/n_bin,':k','HandleVisibility','off')
    xlim(ax,[-180 180]); xticks(ax,-180:90:180); xlabel(ax,'bump position re: reference (deg)'); ylabel(ax,'fraction of stim-window samples')
    legend(ax,'Location','northwest','FontSize',8,'Box','off'); title(ax,{sprintf('reference: %s',ref_label{a}),'where the bump was during pulses'},'FontSize',9)
    for gi = 1:2
        ax = nexttile(tl); hold(ax,'on')
        g = walk_fly_group==gi & ~isnan(F(a).r_cos);
        plot(ax, deg, F(a).d_mob(g,:)', '-', 'Color', [gcol(gi,:) .25])
        m = mean(F(a).d_mob(g,:),1,'omitnan'); se = std(F(a).d_mob(g,:),[],1,'omitnan')./sqrt(sum(~isnan(F(a).d_mob(g,:)),1));
        errorbar(ax, deg, m, se, '-o', 'Color', gcol(gi,:), 'MarkerFaceColor', gcol(gi,:), 'LineWidth', 2)
        yline(ax,0,':k'); xline(ax,0,':k')
        xlim(ax,[-180 180]); xticks(ax,-180:90:180); xlabel(ax,'bump position re: reference (deg)'); ylabel(ax,'\Delta local mobility (dark post - pre)')
        title(ax,sprintf('%s: mobility change by position (n=%d)',geno_label{gi},sum(g)),'FontSize',9)
        all_dmob_axes(end+1) = ax; %#ok<SAGROW>
    end
    ax = nexttile(tl); hold(ax,'on')
    for gi = 1:2
        g = walk_fly_group==gi & ~isnan(F(a).r_cos);
        scatter(ax, gi + .08*randn(sum(g),1), F(a).r_cos(g), 40, gcol(gi,:), 'filled', 'MarkerFaceAlpha', .7)
        errorbar(ax, gi+.3, mean(F(a).r_cos(g)), std(F(a).r_cos(g))/sqrt(sum(g)), 'ok', 'MarkerFaceColor','k', 'LineWidth', 2)
        text(ax, gi, 1.07, sprintf('p=%.2g', signrank(F(a).r_cos(g))), 'HorizontalAlignment','center','FontSize',9)
    end
    yline(ax,0,':k'); xlim(ax,[.5 2.5]); ylim(ax,[-1 1.15]); xticks(ax,[1 2]); xticklabels(ax,geno_label)
    ylabel(ax,'Spearman: \Delta local mobility vs closeness to reference'); title(ax,'per fly, across position bins (+ = more mobile near it)','FontSize',9)
end
linkaxes(all_dmob_axes,'y')
title(tl,sprintf('Is the post-ATP mobility increase larger near where LPsP was stimulated? (dark trials; %d bins; stim window 0-%d s after each confirmed pulse)',n_bin,stim_win),'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_stim_location_vs_mobile_location.png'), 'Resolution', 150)
fig.Visible = 'on';

%% local functions
function [dm, dh] = binned_paths(trials, theta, edges, sm, thr, max_gap, min_walk)
    % per position bin (bump position re: theta), summed |d mu| and
    % |r_speed|*dt over walking-bout samples -- Figure 3's per-bout path
    % lengths, accumulated by where the bump was instead of by bout
    dm = zeros(1,numel(edges)-1); dh = dm;
    for k = 1:numel(trials)
        tr = trials(k); xf = tr.ft.xf(:); dt = median(diff(xf));
        r_s = smoothdata(tr.ft.r_speed(:),'gaussian',sm);
        iw = abs(r_s) > thr; iw = imclose(iw,ones(max_gap+1,1)); iw = bwareaopen(iw,min_walk);
        mu = smoothdata(interp1(tr.ft.xb(:), unwrap(tr.im.mu(:)), xf, 'linear','extrap'),'gaussian',sm);
        in = iw(1:end-1) & iw(2:end); % sample pairs inside one bout
        d_mu = abs(diff(mu)); d_h = abs(r_s(1:end-1))*dt;
        b = discretize(angle(exp(1i*(mu(1:end-1)-theta))), edges);
        for j = 1:numel(dm)
            s = in & b==j;
            dm(j) = dm(j) + sum(d_mu(s),'omitnan'); dh(j) = dh(j) + sum(d_h(s),'omitnan');
        end
    end
end
