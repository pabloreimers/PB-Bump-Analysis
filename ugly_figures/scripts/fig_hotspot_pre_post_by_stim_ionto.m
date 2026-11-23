%% fig_hotspot_pre_post_by_stim_ionto
% Scratch: iontophoresis dataset (lpsp_p2x2_reredo_20260518, Figure 3 rows
% 1-2). Bump occupancy and EPG fluorescence at the ATP puff location, ~20 s
% after each stim vs BEFORE THE FLY WAS EVER PUFFED, split by stim number.
% Companion to fig_hotspot_occupancy_fluorescence_ionto.m (same measures).
%
% TRIALS: per fly in acquisition order (folder name "YYYYMMDD-N_"). Every
% fly's first trial is open-loop bar (O), its second dark (D), then the
% pipette is MOVED (trial 3 on: the other PB half, 85-179 deg away in
% almost every fly -- checked 2026-11-16). Per explicit request, trials
% after the pipette moves are excluded: the SITE is set by the first
% trial with >= min_stims confirmed ejections (whole-PB atp.d median 0-3 s
% minus median -2..0 s > 3 baseline SDs); later trials are kept while they
% stay in the same puffed PB half and within same_site deg of it, stopping
% at the first that doesn't. In practice: trials 1-2, stims 1-4 open loop
% then 5-8 dark. (20260514 fly 1 and fly 5 had no confirmed ejection in
% trial 1, so their site/numbering starts at trial 2.)
% PRE: base_win of the fly's FIRST trial (before its first TTL) -- the only
%   never-puffed period. It is OPEN LOOP, so stims 5-8 (dark) also carry
%   the open-loop -> dark change; shaded in the summary panels.
% POST: from post_start s after each confirmed stim until the next TTL
%   (any TTL, confirmed or not) -- 20-60 s for stims 1-3 of a trial -- or,
%   for a trial's last stim, until the trial ends (20-79 s). Stim number counts TTLs from
%   the site's first trial (4 per trial); unconfirmed stims are skipped.
% HOTSPOT: theta_hot = PVA of the confirmed stims' pooled 0-3 s atp.d
%   profile over all kept trials, in the puffed half (as the companion
%   script); every kept trial is aligned to it.
% MEASURES (n_bin bins re: theta_hot): occupancy (bump = im.mu, the
%   both-halves PVA, rho > rho_gate); raw F (window-mean im.f per wedge /
%   that window's all-wedge mean -- comparable across trials); dF/F
%   (im.d: the dataset's PER-TRIAL dF/F, so for stims 5-8 it compares a
%   different trial's normalization with the pre -- read raw F first).
% ABSOLUTE RAW F (second figure): the relative measures above can't tell
%   "hotspot brighter" from "rest of the PB dimmer". So, un-normalized
%   im.f averaged over 4 wedge groups -- hotspot (puffed half, +/-f_half),
%   mirror (other half, same angles), rest of the puffed half (> 2*f_half
%   from the hotspot), whole PB (all 32) -- as % change from the same
%   group's never-puffed pre. Bleaching pulls every group down; a real
%   local brightening shows as hotspot > mirror and hotspot > whole PB.
% STATS at the hotspot, per fly and stim: occupancy change within
%   +/-occ_half; fluorescence change in the puffed half minus the other
%   half within +/-f_half. Signrank vs 0 per genotype, ranksum between.
%
% REQUIRES the Figure_3_claude.m workspace (all_data, genotype, fly_id).

%% check workspace
need = {'all_data','genotype','fly_id','gcamp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');
wrap = @(x) angle(exp(1i*x));

%% parameters
n_bin = 16; edges = linspace(-pi,pi,n_bin+1); ctr = edges(1:end-1)+pi/n_bin;
conf_base = [-2 0]; conf_peak = [0 3]; peak_factor = 3; min_stims = 2; rho_gate = .2;
base_win = [1 39];  % s into the fly's first trial (first TTL at 41 s)
post_start = 20;   % s after each confirmed stim; the window runs from here to the NEXT TTL (or the trial end, for a trial's last stim)
same_site = 60*pi/180; n_stim_max = 8;
site = 2; % 1 = first pipette site (trials 1-2 in practice), 2 = the site after the pipette moved (trials 3-4)
% SITE 2 (per explicit request): its trials start at the first trial with
% confirmed ejections after the site-1 run, extended by the same same-site
% rule. Its PRE is the first trial AT THAT SITE, before its first TTL --
% NOT never-puffed: the fly already had ~8 puffs at site 1.
if site == 1
    site_tag = ''; base_desc = sprintf('the fly''s never-puffed baseline (first %d-%d s of trial 1)', base_win);
else
    site_tag = sprintf('_site%d', site); base_desc = sprintf('the pre-first-puff baseline of the first trial at site %d (%d-%d s; fly already puffed at site 1)', site, base_win);
end
occ_half = pi/4; f_half = pi/8;
meas_names = {'bump occupancy','raw F (wedge / all-wedge mean)','dF/F'};
geno_order = {'empty>p2x2','lpsp>p2x2'}; geno_label = {'empty>P2X2','lpsp>P2X2'};

%% per fly
n = numel(all_data); acq = nan(n,1);
for i = 1:n; tok = regexp(all_data(i).meta,'\d{8}-(\d+)_','tokens','once'); acq(i) = str2double(tok{1}); end
is_dark = arrayfun(@(x) contains(x.ft.pattern,'background'), all_data);
[fly_list,~,fly_ix] = unique(fly_id); nF = numel(fly_list);
fly_geno = nan(nF,1); Ldark = nan(nF,n_stim_max); % per fly: was stim sn in a dark trial
Dm = repmat({nan(nF,n_bin,n_stim_max,2)},1,3); % (fly, bin, stim number, half)
kept = cell(nF,1);
grp_names = {'hotspot (puffed half)','mirror (other half)','rest of puffed half','whole PB'};
Fabs = nan(nF,n_stim_max,4); % % change of absolute raw F from the fly's never-puffed pre, per wedge group
for f = 1:nF
    r = find(fly_ix==f); [~,o] = sort(acq(r)); r = r(o);
    fly_geno(f) = find(strcmp(geno_order, genotype{r(1)}));
    % per trial: confirmed stims, profile, half, angle
    P = struct('t_on',{},'prof',{},'half',{},'th',{},'ok',{});
    for k = 1:numel(r); P(k) = trial_atp(all_data(r(k)), conf_base, conf_peak, peak_factor, min_stims); end
    k0 = find([P.ok],1); if isempty(k0); continue; end
    ks = k0;
    for k = k0+1:numel(r)
        if ~P(k).ok || P(k).half ~= P(k0).half || abs(wrap(P(k).th - P(k0).th)) > same_site; break; end
        ks(end+1) = k; %#ok<AGROW>
    end
    for st_ = 2:site % advance to the next site: first confirmed-ejection trial after the previous site's run
        k0 = find([P.ok] & (1:numel(r)) > ks(end), 1); if isempty(k0); ks = []; break; end
        ks = k0;
        for k = k0+1:numel(r)
            if ~P(k).ok || P(k).half ~= P(k0).half || abs(wrap(P(k).th - P(k0).th)) > same_site; break; end
            ks(end+1) = k; %#ok<AGROW>
        end
    end
    if isempty(ks); continue; end
    kept{f} = acq(r(ks))';
    % site hotspot from all kept trials
    prof = sum([P(ks).prof],2); puffed = (1:16) + 16*(P(k0).half==2); other = setdiff(1:32,puffed);
    alpha = all_data(r(k0)).im.alpha(:); w = max(prof(puffed),0);
    h = struct('th', angle(sum(w .* exp(1i*alpha(puffed)))), 'puffed', puffed, 'other', other);
    % pre: first trial, before its first TTL
    if site == 1; tr1 = all_data(r(1)); else; tr1 = all_data(r(k0)); end % pre: never-puffed (site 1) / first trial at this site (site 2)
    xb1 = tr1.ft.xb(:);
    [oB, FB, dB] = window_profiles(tr1, xb1>=base_win(1) & xb1<=base_win(2), h, edges, rho_gate);
    G_pre = abs_groups(tr1, xb1>=base_win(1) & xb1<=base_win(2), h, f_half);
    % post: each confirmed stim of the kept trials, numbered by TTL from the site's first trial
    for j = 1:numel(ks)
        tr = all_data(r(ks(j))); xb = tr.ft.xb(:); xf = tr.ft.xf(:);
        st = logical(tr.ft.stims(:)); t_all = xf(diff([false;st])==1);
        for s = 1:numel(t_all)
            sn = (j-1)*4 + s; if sn > n_stim_max || ~any(abs(P(ks(j)).t_on - t_all(s)) < 1e-6); continue; end
            t_end = xb(end); if s < numel(t_all); t_end = t_all(s+1); end
            kw = xb > t_all(s)+post_start & xb < t_end;
            [oW, FW, dW] = window_profiles(tr, kw, h, edges, rho_gate);
            Dm{1}(f,:,sn,1) = oW - oB;
            Fabs(f,sn,:) = 100*(abs_groups(tr, kw, h, f_half) ./ G_pre - 1);
            Dm{2}(f,:,sn,:) = reshape((FW - FB).', 1, n_bin, 1, 2);
            Dm{3}(f,:,sn,:) = reshape((dW - dB).', 1, n_bin, 1, 2);
            Ldark(f,sn) = is_dark(r(ks(j)));
        end
    end
end
stim_dark = round(mean(Ldark,1,'omitnan')); % majority lighting per stim number (titles / shading)
odd = find(any(Ldark ~= stim_dark & ~isnan(Ldark),2));
if ~isempty(odd); fprintf('flies whose lighting differs from the majority at some stim number: %s\n', strjoin(fly_list(odd),', ')); end
fprintf('site %d: kept trials per fly (acquisition numbers):\n', site); for f = 1:nF; fprintf('  %-38s %s\n', fly_list{f}, mat2str(kept{f})); end

%% stats at the hotspot, per stim number
hot_occ = abs(ctr) <= occ_half; hot_f = abs(ctr) <= f_half;
S = nan(nF,n_stim_max,3);
S(:,:,1) = squeeze(mean(Dm{1}(:,hot_occ,:,1),2,'omitnan'));
for m = 2:3; S(:,:,m) = squeeze(mean(Dm{m}(:,hot_f,:,1),2,'omitnan') - mean(Dm{m}(:,hot_f,:,2),2,'omitnan')); end
P_g = nan(2,n_stim_max,3); P_b = nan(n_stim_max,3);
for m = 1:3
    fprintf('site %d, %s (post = +%d s after stim until the next stim)\n', site, meas_names{m}, post_start);
    for sn = 1:n_stim_max
        line = sprintf('  stim %d (%s)', sn, ifelse(stim_dark(sn),'dark','OL bar'));
        for gi = 1:2
            x = S(fly_geno==gi,sn,m); x = x(~isnan(x));
            if numel(x) >= 4; P_g(gi,sn,m) = signrank(x); end
            line = [line sprintf(' | %s n=%2d med %+.3f p=%.3f', geno_label{gi}(1:4), numel(x), median(x), P_g(gi,sn,m))]; %#ok<AGROW>
        end
        x1 = S(fly_geno==1,sn,m); x2 = S(fly_geno==2,sn,m); P_b(sn,m) = ranksum(x1(~isnan(x1)), x2(~isnan(x2)));
        fprintf('%s | between p=%.3f\n', line, P_b(sn,m));
    end
end

%% plot: rows = measures; columns = stim 1..n_stim_max profiles + a summary column
gcol = [.5 .5 .5; gcamp_color]; deg = ctr*180/pi;
fig = figure('color','w','Visible','off','Position',[10 10 2700 1050]);
tl = tiledlayout(fig,3,n_stim_max+2,'TileSpacing','compact','Padding','compact');
ylabs = {'\Delta occupancy (fraction)','\Delta raw F (wedge / all-wedge mean)','\Delta dF/F'};
for m = 1:3
    row_ax = gobjects(1,n_stim_max);
    for sn = 1:n_stim_max
        ax = nexttile(tl); hold(ax,'on')
        if stim_dark(sn); set(ax,'Color',[.95 .95 .95]); end
        halves = 1; if m > 1; halves = 1:2; end
        for gi = 1:2
            for hf = halves
                Y = Dm{m}(fly_geno==gi,:,sn,hf); Y = Y(~all(isnan(Y),2),:); if isempty(Y); continue; end
                ls = '-'; if hf == 2; ls = '--'; end
                errorbar(ax, deg + (gi-1.5)*4, mean(Y,1,'omitnan'), std(Y,[],1,'omitnan')./sqrt(sum(~isnan(Y),1)), ls, 'Color', gcol(gi,:), 'LineWidth', 1.2, 'CapSize', 2)
            end
        end
        yline(ax,0,':k'); xline(ax,0,'--','Color',[.8 0 0])
        xlim(ax,[-180 180]); xticks(ax,[-180 0 180])
        text(ax,.03,.97,sprintf('e p=%.2g\nl p=%.2g',P_g(1,sn,m),P_g(2,sn,m)),'Units','normalized','VerticalAlignment','top','FontSize',7)
        if m == 1; title(ax, sprintf('stim %d (%s)', sn, ifelse(stim_dark(sn),'dark','OL bar')), 'FontSize', 9); end
        if m == 3; xlabel(ax,'re: hotspot (deg)'); end
        if sn == 1; ylabel(ax, ylabs{m}); end
        row_ax(sn) = ax;
    end
    linkaxes(row_ax,'y')
    % summary: hotspot statistic vs stim number
    ax = nexttile(tl,[1 2]); hold(ax,'on')
    patch(ax,[4.5 n_stim_max+.5 n_stim_max+.5 4.5],[-1 -1 1 1]*10,[.93 .93 .93],'EdgeColor','none','HandleVisibility','off')
    for gi = 1:2
        X = S(fly_geno==gi,:,m);
        errorbar(ax, (1:n_stim_max)+(gi-1.5)*.15, mean(X,1,'omitnan'), std(X,[],1,'omitnan')./sqrt(sum(~isnan(X),1)), '-o', ...
            'Color', gcol(gi,:), 'MarkerFaceColor', gcol(gi,:), 'LineWidth', 1.5, 'DisplayName', geno_label{gi})
        sig = find(P_g(gi,:,m) < .05);
        if ~isempty(sig); text(ax, sig+(gi-1.5)*.15, mean(X(:,sig),1,'omitnan') + std(X(:,sig),[],1,'omitnan')./sqrt(sum(~isnan(X(:,sig)),1)), '*', 'Color', gcol(gi,:), 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment','center', 'VerticalAlignment','bottom'); end
    end
    yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[.5 n_stim_max+.5]); xticks(ax,1:n_stim_max)
    mu_all = []; for gi = 1:2; X = S(fly_geno==gi,:,m); mu_all = [mu_all, mean(X,1,'omitnan')+[-1;1].*std(X,[],1,'omitnan')./sqrt(sum(~isnan(X),1))]; end %#ok<AGROW>
    yl = [min([mu_all(:);0]) max([mu_all(:);0])]; ylim(ax, yl + [-.15 .25]*diff(yl)) % scaled to the means +/- SEM, not single flies
    xlabel(ax,'stim number (grey = dark trial)'); ylabel(ax,'hotspot statistic (mean \pm SEM)')
    if m == 1; legend(ax,'Location','best','FontSize',8); end
    title(ax, sprintf('%s at hotspot (* p<.05 vs pre)', meas_names{m}), 'FontSize', 9)
end
title(tl,{sprintf('Iontophoresis dataset, pipette site %d: +%d s after each stim until the next stim (or trial end) vs %s, by stim number; same pipette site only', site, post_start, base_desc), ...
    'solid = puffed PB half, dashed = other half (same angles); fluorescence hotspot statistic = puffed minus other half within \pm22.5 deg; occupancy within \pm45 deg'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,['fig_hotspot_pre_post_by_stim_ionto' site_tag '.png']), 'Resolution', 120)
fig.Visible = 'on';

%% absolute raw F: is the hotspot actually brighter, or the rest dimmer?
fprintf('site %d ABSOLUTE raw F, %% change from %s (median over flies); hotspot-mirror and hotspot-whole PB signrank\n', site, base_desc);
for gi = 1:2
    g = fly_geno==gi;
    for sn = 1:n_stim_max
        X = squeeze(Fabs(g,sn,:)); X = X(~any(isnan(X),2),:);
        fprintf('  %-11s stim %d (%-6s) n=%2d | hot %+5.1f  mirror %+5.1f  rest %+5.1f  PB %+5.1f | hot-mirror %+5.1f p=%.3f | hot-PB %+5.1f p=%.3f\n', geno_label{gi}, sn, ifelse(stim_dark(sn),'dark','OL'), size(X,1), ...
            median(X(:,1)), median(X(:,2)), median(X(:,3)), median(X(:,4)), median(X(:,1)-X(:,2)), signrank(X(:,1),X(:,2)), median(X(:,1)-X(:,4)), signrank(X(:,1),X(:,4)));
    end
end
fig2 = figure('color','w','Visible','off','Position',[30 30 1700 900]);
tl2 = tiledlayout(fig2,2,3,'TileSpacing','compact','Padding','compact');
gcols = [gcamp_color; .55 .55 .55; .3 .75 .3; 0 0 0]; gls = {'-','--','-','-'};
ax2 = gobjects(2,3);
for gi = 1:2
    g = fly_geno==gi;
    % (a) each group's % change
    ax = nexttile(tl2); hold(ax,'on'); ax2(gi,1) = ax;
    patch(ax,[4.5 n_stim_max+.5 n_stim_max+.5 4.5],[-1 -1 1 1]*100,[.93 .93 .93],'EdgeColor','none','HandleVisibility','off')
    for q = 1:4
        X = Fabs(g,:,q);
        errorbar(ax, (1:n_stim_max)+(q-2.5)*.08, mean(X,1,'omitnan'), std(X,[],1,'omitnan')./sqrt(sum(~isnan(X),1)), gls{q}, 'Color', gcols(q,:), 'Marker','o', 'MarkerFaceColor', gcols(q,:), 'MarkerSize', 4, 'LineWidth', 1.4, 'DisplayName', grp_names{q})
    end
    yline(ax,0,':k','HandleVisibility','off'); xlim(ax,[.5 n_stim_max+.5]); xticks(ax,1:n_stim_max)
    ylabel(ax,'absolute raw F, % change from pre'); xlabel(ax,'stim number (grey = dark)')
    legend(ax,'Location','southwest','FontSize',7); title(ax,sprintf('%s: absolute F by wedge group',geno_label{gi}),'FontSize',9)
    % (b,c) hotspot minus mirror, hotspot minus whole PB, per fly
    comps = {[1 2],[1 4]}; clab = {'hotspot - mirror','hotspot - whole PB'};
    for c = 1:2
        ax = nexttile(tl2); hold(ax,'on'); ax2(gi,c+1) = ax;
        patch(ax,[4.5 n_stim_max+.5 n_stim_max+.5 4.5],[-1 -1 1 1]*100,[.93 .93 .93],'EdgeColor','none','HandleVisibility','off')
        D = Fabs(g,:,comps{c}(1)) - Fabs(g,:,comps{c}(2));
        plot(ax, 1:n_stim_max, D', '-', 'Color', [gcamp_color*(gi==2)+[.5 .5 .5]*(gi==1) .25])
        errorbar(ax, 1:n_stim_max, mean(D,1,'omitnan'), std(D,[],1,'omitnan')./sqrt(sum(~isnan(D),1)), '-ok', 'MarkerFaceColor','k', 'LineWidth', 1.6)
        for sn = 1:n_stim_max
            x = D(:,sn); x = x(~isnan(x)); if numel(x) < 4; continue; end
            pv = signrank(x); if pv < .05; text(ax, sn, mean(x)+std(x)/sqrt(numel(x)), '*', 'FontSize', 16, 'FontWeight','bold', 'HorizontalAlignment','center','VerticalAlignment','bottom'); end
        end
        yline(ax,0,':k'); xlim(ax,[.5 n_stim_max+.5]); xticks(ax,1:n_stim_max)
        ylabel(ax,['\Delta% change: ' clab{c}]); xlabel(ax,'stim number (grey = dark)')
        title(ax,sprintf('%s: %s (* p<.05)',geno_label{gi},clab{c}),'FontSize',9)
    end
end
% y-limits from the data (the +/-100 dark-trial shading must not drive them): column 1 from the
% group means +/- SEM, columns 2-3 from the per-fly differences; shared across genotypes per column
lim1 = [0 0]; lim23 = [0 0];
for gi = 1:2
    X = Fabs(fly_geno==gi,:,:); M = squeeze(mean(X,1,'omitnan')); E = squeeze(std(X,[],1,'omitnan'))./sqrt(squeeze(sum(~isnan(X),1)));
    lim1 = [min([lim1(1); M(:)-E(:)]) max([lim1(2); M(:)+E(:)])];
    D = [X(:,:,1)-X(:,:,2), X(:,:,1)-X(:,:,4)]; lim23 = [min([lim23(1); D(:)]) max([lim23(2); D(:)])];
end
set(ax2(:,1),'YLim',lim1 + [-.1 .1]*diff(lim1)); set(ax2(:,2:3),'YLim',lim23 + [-.05 .1]*diff(lim23));
title(tl2,{sprintf('Pipette site %d. Absolute (un-normalized) raw F, +%d s after each stim until the next, as %% change from %s',site,post_start,base_desc), ...
    'bleaching pulls everything down; a local brightening = hotspot falls less than its mirror / the whole PB'},'FontSize',10)
exportgraphics(fig2, fullfile(export_dir,['fig_hotspot_pre_post_by_stim_ionto_absF' site_tag '.png']), 'Resolution', 130)
fig2.Visible = 'on';

%% local functions
function s = ifelse(c, a, b)
    if c; s = a; else; s = b; end
end

function G = abs_groups(tr, keep, h, f_half)
    % un-normalized window-mean im.f over 4 wedge groups: hotspot (puffed
    % half, |rel| <= f_half), mirror (other half, same), rest of the puffed
    % half (|rel| > 2*f_half), whole PB
    wrap = @(x) angle(exp(1i*x));
    fm = mean(tr.im.f(:,keep),2); rel = abs(wrap(tr.im.alpha(:) - h.th));
    inP = false(32,1); inP(h.puffed) = true;
    G = [mean(fm(inP & rel <= f_half)), mean(fm(~inP & rel <= f_half)), mean(fm(inP & rel > 2*f_half)), mean(fm)];
end

function P = trial_atp(tr, base_win, peak_win, peak_factor, min_stims)
    % confirmed stims of one trial + their pooled baseline-subtracted 0-3 s
    % atp.d profile, puffed half and PVA angle
    xb = tr.ft.xb(:); xf = tr.ft.xf(:);
    st = logical(tr.ft.stims(:)); t_all = xf(diff([false;st])==1);
    atp_m = mean(tr.atp.d,1)'; prof = zeros(size(tr.atp.d,1),1); keep = false(numel(t_all),1);
    for k = 1:numel(t_all)
        b = xb>=t_all(k)+base_win(1) & xb<t_all(k)+base_win(2);
        q = xb>t_all(k)+peak_win(1) & xb<=t_all(k)+peak_win(2);
        keep(k) = median(atp_m(q)) - median(atp_m(b)) > peak_factor*std(atp_m(b));
        if keep(k); prof = prof + mean(tr.atp.d(:,q),2) - mean(tr.atp.d(:,b),2); end
    end
    [~,pk] = max(prof); half = 1 + (pk>16); hh = (1:16) + 16*(half==2);
    w = max(prof(hh),0); alpha = tr.im.alpha(:);
    P = struct('t_on', t_all(keep), 'prof', prof, 'half', half, 'th', angle(sum(w .* exp(1i*alpha(hh)))), 'ok', sum(keep) >= min_stims);
end

function [occ, F, d] = window_profiles(tr, keep, h, edges, rho_gate)
    % samples 'keep' (imaging clock): occupancy re: h.th; per-half [2 x bin]
    % raw F (normalized by that window's all-wedge mean) and dF/F
    nb = numel(edges)-1; wrap = @(x) angle(exp(1i*x));
    mu = tr.im.mu(:); ok = keep(:) & tr.im.rho(:) > rho_gate;
    occ = histcounts(wrap(mu(ok)-h.th), edges, 'Normalization','probability');
    fm = mean(tr.im.f(:,keep),2); fm = fm / mean(fm); dm = mean(tr.im.d(:,keep),2);
    F = nan(2,nb); d = nan(2,nb); halves = {h.puffed, h.other};
    for hf = 1:2
        wv = halves{hf}; b = discretize(wrap(tr.im.alpha(wv)'-h.th), edges);
        F(hf,:) = accumarray(b(:), fm(wv), [nb 1], @mean, NaN)';
        d(hf,:) = accumarray(b(:), dm(wv), [nb 1], @mean, NaN)';
    end
end
