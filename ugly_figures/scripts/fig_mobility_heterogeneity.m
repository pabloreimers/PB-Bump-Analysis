%% fig_mobility_heterogeneity
% Scratch: is post-ATP dark bump mobility uneven ACROSS THE BRIDGE, and if
% so, do bump position during stimulation, the ATP hotspot, or their
% coincidence predict WHERE? Walking dataset
% (lpsp_p2x2_walking_20260728_post0708); everything in the absolute bump-
% angle frame (im.mu wrapped, n_bin bins), per fly, same pre/post dark
% pools and walking-bout definition as Figure 3. Follows
% fig_stim_location_vs_mobile_location.m (same theta_hot definition and
% confirmed-pulse test; bump-during-stim window 0..stim_win s).
%
% Q1 HETEROGENEITY. Position-resolved mobility per bin = sum |d mu| /
%   sum |r_speed|*dt over walking-bout samples whose bump sits in that bin
%   (bins need >= min_head rad of heading path). Unevenness = CV of the
%   bin mobilities (std/mean over valid bins). NULL: rotate each walking
%   bout's bump positions by its own random angle, re-bin, recompute CV
%   (n_null times) -- keeps every bout's own path lengths, destroys any
%   consistent position -> mobility mapping. het_z = (CV - mean null)/sd
%   null, het_p = fraction of null CVs >= observed. Run on pre AND post: a
%   fixed imaging/PVA artifact (e.g. a dim glomerulus making mu jumpy
%   there) would show up in both; what's new is post-only structure.
%   het_z grows with the amount of data (a tighter null), and flies have
%   ~2 dark pre trials vs ~1 post, so pre vs post is compared on an effect
%   size instead: excess CV = observed CV - mean null CV, with pre
%   SUBSAMPLED to the post bout count (n_sub random bout subsets).
% Q2 PREDICTORS (within fly, pooled per genotype). y = d_mob(bin) = post -
%   pre, z-scored across that fly's valid bins (removes the fly's overall
%   increase -- only WHERE is tested). Predictors per bin, z-scored within
%   fly: occ = fraction of the bump's stim-window samples in that bin;
%   hot = cos(bin centre - theta_hot); occ x hot = "bump sat on the
%   hotspot". OLS on all flies' bins; NULL: circularly shift each fly's
%   d_mob profile by a random nonzero number of bins (keeps the profile's
%   shape, breaks its alignment to the predictors), n_perm times.
% Q3 COINCIDENCE, per fly: overlap = fraction of stim-window bump samples
%   within +/-45 deg of theta_hot. Spearman vs the fly's overall dark
%   mobility change (Figure 3's walk_delta) and vs post het_z.
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
min_head = pi; n_null = 500; n_perm = 2000; overlap_half = pi/4; n_sub = 30; n_null_sub = 100;
rng(1); % reproducible nulls

%% per fly
nF = numel(walk_fly_list);
D = struct('mob_pre',nan(nF,n_bin),'mob_post',nan(nF,n_bin),'d_mob',nan(nF,n_bin),'occ',nan(nF,n_bin), ...
    'theta_hot',nan(nF,1),'overlap',nan(nF,1),'het_z',nan(nF,2),'het_p',nan(nF,2),'null_band',nan(nF,n_bin,2), ...
    'exc',nan(nF,2),'exc_pre_matched',nan(nF,1),'n_bout',nan(nF,2));
for f = 1:nF
    if isnan(walk_delta(f)); continue; end
    ti = find(walk_fly_ix==f); [~,o] = sort(walk_trial_num(ti)); ti = ti(o);
    p = find(walk_has_pulse(ti),1);
    pre  = ti(1:p-1); pre = pre(walk_is_dark(pre));
    post = ti(p+1:end); post = post(~walk_has_pulse(post) & walk_is_dark(post));

    % pulse trial: bump during stim (0..stim_win s) + ATP hotspot
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
    [~,pk] = max(prof); half = (1:16) + 16*(pk>16); w = max(prof(half),0); alpha = tr.im.alpha(:);
    D.theta_hot(f) = angle(sum(w .* exp(1i*alpha(half))));
    D.occ(f,:) = histcounts(wrap(mu_s), edges, 'Normalization','probability');
    D.overlap(f) = mean(abs(wrap(mu_s - D.theta_hot(f))) <= overlap_half);

    % position-resolved mobility (absolute frame) + rotation null, pre and post
    S = {bout_samples(walk_data(pre), walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames), ...
         bout_samples(walk_data(post), walk_smooth_window, walk_turn_thresh, walk_max_gap_frames, walk_min_walking_frames)};
    mob = nan(2,n_bin);
    for c = 1:2
        [m, ok] = bin_mobility(S{c}, zeros(max(S{c}.bout),1), edges, min_head);
        mob(c,ok) = m(ok);
        cv_obs = std(m(ok))/mean(m(ok));
        cv_null = nan(n_null,1); null_prof = nan(n_null,n_bin);
        for r = 1:n_null
            [mn, okn] = bin_mobility(S{c}, 2*pi*rand(max(S{c}.bout),1), edges, min_head);
            cv_null(r) = std(mn(okn))/mean(mn(okn)); null_prof(r,okn) = mn(okn);
        end
        D.het_z(f,c) = (cv_obs - mean(cv_null,'omitnan'))/std(cv_null,'omitnan');
        D.het_p(f,c) = mean(cv_null >= cv_obs);
        D.null_band(f,:,c) = prctile(null_prof,97.5) - prctile(null_prof,2.5);
        D.exc(f,c) = cv_obs - mean(cv_null,'omitnan');
        D.n_bout(f,c) = numel(unique(S{c}.bout));
    end
    % pre subsampled to the post bout count -> data-matched excess CV
    ub = unique(S{1}.bout); e = nan(n_sub,1);
    for d = 1:n_sub
        keep = ismember(S{1}.bout, ub(randperm(numel(ub), min(D.n_bout(f,2), numel(ub)))));
        Ss = structfun(@(v) v(keep), S{1}, 'UniformOutput', false);
        [m, ok] = bin_mobility(Ss, zeros(max(S{1}.bout),1), edges, min_head);
        cvn = nan(n_null_sub,1);
        for r = 1:n_null_sub
            [mn, okn] = bin_mobility(Ss, 2*pi*rand(max(S{1}.bout),1), edges, min_head); cvn(r) = std(mn(okn))/mean(mn(okn));
        end
        e(d) = std(m(ok))/mean(m(ok)) - mean(cvn,'omitnan');
    end
    D.exc_pre_matched(f) = mean(e,'omitnan');
    D.mob_pre(f,:) = mob(1,:); D.mob_post(f,:) = mob(2,:);
    D.d_mob(f,:) = mob(2,:) - mob(1,:);
end

%% Q1: heterogeneity
geno_label = {'empty>P2X2','lpsp>P2X2'}; gcol = [.5 .5 .5; gcamp_color]; cond = {'pre','post'};
fprintf('Q1 HETEROGENEITY (CV of bin mobility vs per-bout rotation null, %d rotations)\n', n_null);
for gi = 1:2
    g = walk_fly_group==gi & ~isnan(D.het_z(:,2));
    for c = 1:2
        fprintf('  %-11s %-4s n=%2d  het_z median %+.2f  flies p<.05: %d/%d  | signrank het_z vs 0 p=%.3f\n', ...
            geno_label{gi}, cond{c}, sum(g), median(D.het_z(g,c)), sum(D.het_p(g,c)<.05), sum(g), signrank(D.het_z(g,c)));
    end
    fprintf('  %-11s bouts pre/post median %d/%d | excess CV: pre %.2f, pre bout-matched %.2f, post %.2f | post vs matched pre signrank p=%.3f\n', ...
        geno_label{gi}, median(D.n_bout(g,1)), median(D.n_bout(g,2)), median(D.exc(g,1)), median(D.exc_pre_matched(g)), median(D.exc(g,2)), signrank(D.exc(g,2), D.exc_pre_matched(g)));
    % is the spatial pattern stable pre -> post? (a fixed artifact would be)
    r_pp = arrayfun(@(f) corr(D.mob_pre(f,:)', D.mob_post(f,:)', 'type','Spearman','rows','complete'), find(g));
    fprintf('  %-11s pre vs post bin-mobility profile, per-fly Spearman median %+.2f, signrank p=%.3f\n', geno_label{gi}, median(r_pp,'omitnan'), signrank(r_pp));
end

%% Q2: within-fly predictors of WHERE mobility changed
pred_names = {'bump during stim (occ)','ATP hotspot (cos)','bump on hotspot (occ x cos)'};
fprintf('Q2 PREDICTORS of z-scored d_mob across bins (pooled OLS, %d circular-shift permutations)\n', n_perm);
B = nan(2,3); Bp = nan(2,3); Bnull = cell(1,2);
for gi = 1:2
    flies = find(walk_fly_group==gi & ~isnan(D.theta_hot) & sum(~isnan(D.d_mob),2) >= 5);
    [Y, X] = q2_design(D, flies, ctr, 0);
    B(gi,:) = (X \ Y)';
    bn = nan(n_perm,3);
    for r = 1:n_perm
        [Yr, Xr] = q2_design(D, flies, ctr, randi(n_bin-1, numel(flies), 1));
        bn(r,:) = (Xr \ Yr)';
    end
    Bnull{gi} = bn;
    Bp(gi,:) = mean(abs(bn) >= abs(B(gi,:)), 1);
    for k = 1:3
        fprintf('  %-11s n=%2d flies  beta %-28s %+.2f  perm p=%.3f\n', geno_label{gi}, numel(flies), pred_names{k}, B(gi,k), Bp(gi,k));
    end
end

%% Q3: fly-level coincidence
fprintf('Q3 COINCIDENCE (fraction of stim-window bump samples within +/-%d deg of the hotspot)\n', round(overlap_half*180/pi));
for gi = 1:2
    g = walk_fly_group==gi & ~isnan(D.overlap) & ~isnan(D.het_z(:,2));
    [r1,p1] = corr(D.overlap(g), walk_delta(g), 'type','Spearman');
    [r2,p2] = corr(D.overlap(g), D.het_z(g,2), 'type','Spearman');
    fprintf('  %-11s n=%2d overlap median %.2f (chance %.2f) | vs overall d_mobility r_s %+.2f p=%.2f | vs post het_z r_s %+.2f p=%.2f\n', ...
        geno_label{gi}, sum(g), median(D.overlap(g)), 2*overlap_half/(2*pi), r1, p1, r2, p2);
end

%% plot
fig = figure('color','w','Visible','off','Position',[40 40 1800 950]);
tl = tiledlayout(fig,2,4,'TileSpacing','compact','Padding','compact');
deg = ctr*180/pi;
% 1: example lpsp fly, pre/post bin mobility (largest post het_z)
gl = find(walk_fly_group==2 & ~isnan(D.het_z(:,2))); [~,k] = max(D.het_z(gl,2)); fx = gl(k);
ax = nexttile(tl); hold(ax,'on')
plot(ax, deg, D.mob_pre(fx,:), '-o', 'Color', [.5 .5 .5], 'MarkerFaceColor', [.5 .5 .5], 'LineWidth', 1.5)
plot(ax, deg, D.mob_post(fx,:), '-o', 'Color', gcamp_color, 'MarkerFaceColor', gcamp_color, 'LineWidth', 1.5)
xline(ax, D.theta_hot(fx)*180/pi, '--', 'Color', [.8 0 0], 'LineWidth', 1.2)
xlim(ax,[-180 180]); xticks(ax,-180:90:180); xlabel(ax,'bump position (deg, PB frame)'); ylabel(ax,'local mobility (bump/heading path)')
pstr0 = @(q) ifelse_p(q, n_null);
legend(ax,{sprintf('dark pre (het z %.1f, %s)',D.het_z(fx,1),pstr0(D.het_p(fx,1))), sprintf('dark post (het z %.1f, %s)',D.het_z(fx,2),pstr0(D.het_p(fx,2))),'ATP hotspot'},'Location','best','FontSize',8,'Box','off')
title(ax,{'1. most uneven lpsp fly (post)', strrep(walk_fly_list{fx},'Z:/pablo/lpsp_p2x2_walking/','')},'FontSize',9,'Interpreter','none')
% 2: excess CV (effect size), bout-matched pre vs post, per fly, by genotype
ax = nexttile(tl); hold(ax,'on')
for gi = 1:2
    g = find(walk_fly_group==gi & ~isnan(D.exc(:,2)));
    xo = (gi-1)*3; E = [D.exc_pre_matched(g), D.exc(g,2)];
    plot(ax, xo+[1 2], E', '-', 'Color', [gcol(gi,:) .35])
    scatter(ax, xo+ones(numel(g),1), E(:,1), 30, gcol(gi,:), 'filled'); scatter(ax, xo+2*ones(numel(g),1), E(:,2), 30, gcol(gi,:), 'filled')
    text(ax, xo+1.5, -.12, sprintf('post vs pre* p=%.2g', signrank(E(:,2),E(:,1))), 'HorizontalAlignment','center','FontSize',8)
end
yline(ax,0,':k'); xticks(ax,[1 2 4 5]); xticklabels(ax,{'pre*','post','pre*','post'}); xlim(ax,[.5 5.5])
xlabel(ax,'empty>P2X2           lpsp>P2X2    (*pre subsampled to post bout count)')
ylabel(ax,'unevenness across PB (CV - null CV)'); title(ax,'2. is mobility uneven across the bridge? (0 = no more than chance)','FontSize',9)
% 3: Q2 coefficients with permutation null
ax = nexttile(tl,3,[1 2]); hold(ax,'on')
for gi = 1:2
    xs = (1:3) + (gi-1.5)*.25;
    lo = prctile(Bnull{gi},2.5); hi = prctile(Bnull{gi},97.5);
    for k = 1:3; plot(ax, [xs(k) xs(k)], [lo(k) hi(k)], '-', 'Color', [gcol(gi,:) .4], 'LineWidth', 6, 'HandleVisibility','off'); end
    scatter(ax, xs, B(gi,:), 70, gcol(gi,:), 'filled', 'DisplayName', geno_label{gi})
    for k = 1:3; text(ax, xs(k), B(gi,k), sprintf('  p=%.2g',Bp(gi,k)), 'FontSize', 8, 'Color', gcol(gi,:)); end
end
yline(ax,0,':k','HandleVisibility','off'); xticks(ax,1:3); xticklabels(ax,pred_names); xlim(ax,[.5 3.5])
ylabel(ax,'\beta (z-scored \Delta local mobility per z of predictor)'); legend(ax,'Location','best','FontSize',8,'Box','off')
title(ax,'3. what predicts WHERE mobility changed (within fly; bars = 95% of circular-shift null)','FontSize',9)
% 4: overlap vs overall d_mobility; 5: overlap vs post het_z
yl_list = {walk_delta, D.het_z(:,2)}; ylab = {'\Delta bump mobility (Figure 3, dark post - pre)','post unevenness (het z)'};
for q = 1:2
    ax = nexttile(tl); hold(ax,'on'); h = gobjects(1,2); leg = cell(1,2);
    for gi = 1:2
        g = walk_fly_group==gi & ~isnan(D.overlap) & ~isnan(D.het_z(:,2));
        [rs,ps] = corr(D.overlap(g), yl_list{q}(g), 'type','Spearman');
        h(gi) = scatter(ax, D.overlap(g), yl_list{q}(g), 50, gcol(gi,:), 'filled', 'MarkerFaceAlpha', .7);
        leg{gi} = sprintf('%s: r_s %+.2f, p=%.2g', geno_label{gi}, rs, ps);
    end
    xline(ax, 2*overlap_half/(2*pi), ':k', 'chance', 'HandleVisibility','off')
    xlabel(ax,sprintf('bump-on-hotspot overlap during stim (\\pm%d deg)',round(overlap_half*180/pi))); ylabel(ax,ylab{q})
    legend(ax,h,leg,'Location','southoutside','FontSize',8); title(ax,sprintf('%d. coincidence vs %s',3+q,ylab{q}),'FontSize',9)
end
% 6: mean pre/post bin-mobility profile aligned to hotspot, lpsp (context)
ax = nexttile(tl,[1 2]); hold(ax,'on')
for gi = 1:2
    g = find(walk_fly_group==gi & ~isnan(D.theta_hot));
    Z = nan(numel(g),n_bin);
    for j = 1:numel(g)
        y = D.d_mob(g(j),:); okb = ~isnan(y); if sum(okb) < 5; continue; end
        y(okb) = (y(okb)-mean(y(okb)))/std(y(okb));
        rel = wrap(ctr - D.theta_hot(g(j))); [~,ord] = sort(rel); Z(j,:) = y(ord); % reorder bins by angle re: hotspot
    end
    rel_c = sort(wrap(ctr - mean(ctr))); % nominal relative bin centres
    errorbar(ax, rel_c*180/pi + (gi-1.5)*6, mean(Z,1,'omitnan'), std(Z,[],1,'omitnan')./sqrt(sum(~isnan(Z),1)), '-o', 'Color', gcol(gi,:), 'MarkerFaceColor', gcol(gi,:), 'LineWidth', 1.5, 'DisplayName', geno_label{gi})
end
yline(ax,0,':k','HandleVisibility','off'); xline(ax,0,'--','Color',[.8 0 0],'HandleVisibility','off'); xlim(ax,[-180 180]); xticks(ax,-180:45:180)
xlabel(ax,'bump position re: ATP hotspot (deg, nearest bin)'); ylabel(ax,'within-fly z of \Delta local mobility'); legend(ax,'Location','best','FontSize',8,'Box','off')
title(ax,'6. shape of the change re: hotspot (overall increase removed per fly)','FontSize',9)
title(tl,sprintf('Is post-ATP dark mobility uneven across the PB, and does bump position during stim (0-%d s) and/or the ATP hotspot predict where?',stim_win),'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_mobility_heterogeneity.png'), 'Resolution', 150)
fig.Visible = 'on';

%% local functions
function s = ifelse_p(q, n)
    % permutation p: report the resolution limit instead of p=0
    if q == 0; s = sprintf('p<%.3g',1/n); else; s = sprintf('p=%.2g',q); end
end

function S = bout_samples(trials, sm, thr, max_gap, min_walk)
    % per sample pair inside a walking bout: bump position, |d mu|, |r_speed|*dt, bout id
    S = struct('mu',[],'dm',[],'dh',[],'bout',[]); nb = 0;
    for k = 1:numel(trials)
        tr = trials(k); xf = tr.ft.xf(:); dt = median(diff(xf));
        r_s = smoothdata(tr.ft.r_speed(:),'gaussian',sm);
        iw = abs(r_s) > thr; iw = imclose(iw,ones(max_gap+1,1)); iw = bwareaopen(iw,min_walk);
        mu = smoothdata(interp1(tr.ft.xb(:), unwrap(tr.im.mu(:)), xf, 'linear','extrap'),'gaussian',sm);
        lab = bwlabel(iw); in = iw(1:end-1) & iw(2:end); % sample pairs inside one bout
        id = lab(1:end-1) + nb; nb = nb + max(lab);      % bout ids unique across trials
        p0 = mu(1:end-1); dmu = abs(diff(mu)); dh = abs(r_s(1:end-1))*dt;
        S.mu = [S.mu; p0(in)]; S.dm = [S.dm; dmu(in)]; S.dh = [S.dh; dh(in)]; S.bout = [S.bout; id(in)]; %#ok<AGROW>
    end
end

function [m, ok] = bin_mobility(S, rot, edges, min_head)
    % mobility per position bin, with each bout's positions rotated by rot(bout)
    if isempty(S.bout); m = nan(1,numel(edges)-1); ok = false(size(m)); return; end
    b = discretize(angle(exp(1i*(S.mu + rot(S.bout)))), edges);
    dm = accumarray(b, S.dm, [numel(edges)-1 1])'; dh = accumarray(b, S.dh, [numel(edges)-1 1])';
    ok = dh >= min_head; m = dm./dh;
end

function [Y, X] = q2_design(D, flies, ctr, shift)
    % stacked within-fly z-scored d_mob and predictors; shift(j) circularly shifts fly j's d_mob (null)
    Y = []; X = [];
    for j = 1:numel(flies)
        f = flies(j);
        y = circshift(D.d_mob(f,:), shift(min(j,numel(shift))));
        occ = D.occ(f,:); hot = cos(ctr - D.theta_hot(f));
        okb = ~isnan(y); if sum(okb) < 5; continue; end
        z = @(v) (v(okb)-mean(v(okb)))/std(v(okb));
        x1 = z(occ); x2 = z(hot); x3 = x1.*x2; x3 = (x3-mean(x3))/std(x3);
        Y = [Y; z(y)']; X = [X; [x1' x2' x3']]; %#ok<AGROW>
    end
end
