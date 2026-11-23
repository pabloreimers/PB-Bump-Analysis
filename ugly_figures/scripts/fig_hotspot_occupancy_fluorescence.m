%% fig_hotspot_occupancy_fluorescence
% Scratch: does ATP puffing change (1) where the bump sits and (2) EPG
% fluorescence AT THE PUFF LOCATION -- during the puffs and in the trials
% after? Walking dataset (lpsp_p2x2_walking_20260728_post0708), closed
% loop and dark, everything aligned to each fly's ATP hotspot.
%
% HOTSPOT, per PULSE TRIAL (as fig_stim_location_vs_mobile_location.m):
% average that trial's confirmed pulses' (whole-PB atp.d median 0-3 s minus
% median -2..0 s > 3 baseline SDs) baseline-subtracted 0-3 s atp.d
% profile; the PB half (wedges 1:16 / 17:32) holding its peak is the
% PUFFED half; theta_hot = population-vector angle of that half's positive
% profile over im.alpha (the frame im.mu lives in).
% MOST FLIES HAVE 2-3 PULSE TRIALS whose hotspots differ (first vs second:
% median ~40 deg, up to 120 deg; one fly switched halves). Per explicit
% request, every trial is aligned to the hotspot of the MOST RECENT pulse
% trial at or before it: each pulse trial to its own, each post trial to
% the last pulse trial preceding it. The pre baseline it is compared with
% (same lighting) is re-referenced to that SAME hotspot, so every
% post - pre difference is measured against one location. Within a fly,
% differences are averaged over pulse trials / post trials.
% CONDITIONS, per fly (trial order from the folder name; Figure 3 pools):
%   during : peri-pulse, 0..3 s after each confirmed pulse vs -2..0 s before
%   pulse  : the whole pulse trial (closed loop) vs pre-ATP closed loop
%   post CL: closed-loop trials after the pulse trial vs pre-ATP closed loop
%   post DK: dark trials after vs dark trials before (Figure 3's pools)
% MEASURES, in n_bin angle bins relative to theta_hot:
%   occupancy : fraction of samples with the bump (rho > rho_gate) there.
%               Bump position = im.mu, which is ALREADY the population vector
%               of im.z over all 32 wedges, i.e. BOTH PB halves (checked
%               2026-11-16: 0 deg median difference vs a recomputed 32-wedge
%               PVA; either half alone differs by ~12 deg) -- so occupancy is
%               not tied to the puffed half.
%   raw F     : time-mean im.f per wedge / that trial's mean over all 32
%               wedges (removes bleaching / global drift: only LOCAL change)
%   dF/F      : time-mean im.d per wedge (dataset's own per-trial dF/F)
%   For fluorescence, the PUFFED half and the OTHER half are kept apart:
%   the bump lights the same angle in both halves, so a change in the
%   puffed half that the other half's same angle doesn't share is a local
%   ATP effect, not just more bump time there. During: (F after - F
%   before)/F before per wedge for raw F; dF/F after - before.
% STATS at the hotspot: occupancy within +/-occ_half; fluorescence in the
%   wedges within +/-f_half (puffed half) minus the same angles in the other
%   half. Signrank vs 0 per genotype, ranksum between.
%
% REQUIRES the Figure_3_claude.m workspace (walk_* variables, circ_stats on path).

%% check workspace
need = {'walk_data','walk_fly_ix','walk_fly_list','walk_trial_num','walk_has_pulse','walk_is_dark','walk_fly_group','walk_delta','gcamp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');
wrap = @(x) angle(exp(1i*x));

%% parameters
n_bin = 16; edges = linspace(-pi,pi,n_bin+1); ctr = edges(1:end-1)+pi/n_bin;
base_win = [-2 0]; peak_win = [0 3]; peak_factor = 3; rho_gate = .2;
occ_half = pi/4; f_half = pi/8; % hotspot windows: occupancy +/-45 deg, fluorescence +/-22.5 deg (~nearest 2 wedges)
comp_names = {'during puffs (0-3 s vs -2-0 s)','pulse trials vs pre CL','post CL vs pre CL','post dark vs pre dark'};
meas_names = {'bump occupancy','raw F (wedge / trial mean)','dF/F'};

%% per fly
nF = numel(walk_fly_list);
% Dm{m}(f,bin,comp,half): change in measure m; half 1 = puffed, 2 = other (occupancy uses half 1 only)
Dm = repmat({nan(nF,n_bin,4,2)},1,3); n_pulse_trials = nan(nF,1);
for f = 1:nF
    if isnan(walk_delta(f)); continue; end
    ti = find(walk_fly_ix==f); [~,o] = sort(walk_trial_num(ti)); ti = ti(o);
    pp = find(walk_has_pulse(ti));
    pre = ti(1:pp(1)-1); preCL = pre(~walk_is_dark(pre)); preDK = pre(walk_is_dark(pre));

    % every pulse trial: its own hotspot + confirmed-pulse windows
    H = struct('th',{},'puffed',{},'other',{},'in_a',{},'in_b',{},'pos',{});
    for j = 1:numel(pp)
        [h, ok] = pulse_hotspot(walk_data(ti(pp(j))), base_win, peak_win, peak_factor);
        if ok; h.pos = pp(j); H(end+1) = h; end %#ok<SAGROW>
    end
    if isempty(H); continue; end
    n_pulse_trials(f) = numel(H);

    acc = repmat({nan(n_bin,2,0)},3,4); % acc{m,c}: stacked [bin x half x instance]
    for j = 1:numel(H)
        tr = walk_data(ti(H(j).pos)); h = H(j);
        % 1) during puffs, this pulse trial vs its own hotspot
        [oa, Fa, da] = profiles(tr, h.in_a, h.th, edges, rho_gate, h.puffed, h.other);
        [ob, Fb, db] = profiles(tr, h.in_b, h.th, edges, rho_gate, h.puffed, h.other);
        acc{1,1}(:,:,end+1) = [(oa-ob)' nan(n_bin,1)]; acc{2,1}(:,:,end+1) = ((Fa-Fb)./Fb)'; acc{3,1}(:,:,end+1) = (da-db)';
        % 2) whole pulse trial vs pre CL, at this hotspot
        if ~isempty(preCL)
            [oA, FA, dA] = cond_profiles(tr, h.th, edges, rho_gate, h.puffed, h.other);
            [oB, FB, dB] = cond_profiles(walk_data(preCL), h.th, edges, rho_gate, h.puffed, h.other);
            acc{1,2}(:,:,end+1) = [(oA-oB)' nan(n_bin,1)]; acc{2,2}(:,:,end+1) = (FA-FB)'; acc{3,2}(:,:,end+1) = (dA-dB)';
        end
    end
    % 3-4) each post (non-pulse) trial vs pre of the same lighting, at the most recent pulse trial's hotspot
    for k = pp(1)+1:numel(ti)
        if walk_has_pulse(ti(k)); continue; end
        j = find([H.pos] < k, 1, 'last'); if isempty(j); continue; end
        h = H(j); is_dk = walk_is_dark(ti(k));
        base = preCL; c = 3; if is_dk; base = preDK; c = 4; end
        if isempty(base); continue; end
        [oA, FA, dA] = cond_profiles(walk_data(ti(k)), h.th, edges, rho_gate, h.puffed, h.other);
        [oB, FB, dB] = cond_profiles(walk_data(base), h.th, edges, rho_gate, h.puffed, h.other);
        acc{1,c}(:,:,end+1) = [(oA-oB)' nan(n_bin,1)]; acc{2,c}(:,:,end+1) = (FA-FB)'; acc{3,c}(:,:,end+1) = (dA-dB)';
    end
    for m = 1:3
        for c = 1:4
            if ~isempty(acc{m,c}); Dm{m}(f,:,c,:) = reshape(mean(acc{m,c},3,'omitnan'), 1, n_bin, 1, 2); end
        end
    end
end
fprintf('pulse trials with >= 3 confirmed pulses per fly: %s\n', mat2str(n_pulse_trials(~isnan(n_pulse_trials))'));

%% stats at the hotspot
geno_label = {'empty>P2X2','lpsp>P2X2'}; gcol = [.5 .5 .5; gcamp_color];
hot_occ = abs(ctr) <= occ_half; hot_f = abs(ctr) <= f_half;
S = nan(nF,4,3); % fly x comparison x measure: hotspot statistic
for c = 1:4
    S(:,c,1) = mean(Dm{1}(:,hot_occ,c,1),2,'omitnan');
    for m = 2:3
        S(:,c,m) = mean(Dm{m}(:,hot_f,c,1),2,'omitnan') - mean(Dm{m}(:,hot_f,c,2),2,'omitnan'); % puffed minus other half, same angles
    end
end
fprintf('hotspot statistic: occupancy = change in fraction within +/-%d deg; fluorescence = change in puffed half minus other half, within +/-%.1f deg\n', round(occ_half*180/pi), f_half*180/pi);
for m = 1:3
    for c = 1:4
        line = sprintf('%-27s %-32s', meas_names{m}, comp_names{c}); pg = nan(1,2);
        for gi = 1:2
            x = S(walk_fly_group==gi,c,m); x = x(~isnan(x));
            if numel(x) >= 4; pg(gi) = signrank(x); end
            line = [line sprintf(' | %s n=%2d med %+.3f p=%.3f', geno_label{gi}(1:4), numel(x), median(x), pg(gi))]; %#ok<AGROW>
        end
        x1 = S(walk_fly_group==1,c,m); x2 = S(walk_fly_group==2,c,m);
        fprintf('%s | between p=%.3f\n', line, ranksum(x1(~isnan(x1)), x2(~isnan(x2))));
    end
end

%% plot: rows = measures, columns = comparisons
fig = figure('color','w','Visible','off','Position',[30 30 1900 1150]);
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
deg = ctr*180/pi;
ylabs = {'\Delta occupancy (fraction)','\Delta raw F (wedge / trial mean; during: %)','\Delta dF/F'};
for m = 1:3
    row_ax = gobjects(1,4);
    for c = 1:4
        ax = nexttile(tl); hold(ax,'on'); h = gobjects(0); lab = {};
        halves = 1; if m > 1; halves = 1:2; end
        for gi = 1:2
            g = walk_fly_group==gi;
            for hf = halves
                Y = Dm{m}(g,:,c,hf); Y = Y(~all(isnan(Y),2),:); if isempty(Y); continue; end
                ls = '-'; if hf == 2; ls = '--'; end
                h(end+1) = errorbar(ax, deg + (gi-1.5)*3, mean(Y,1,'omitnan'), std(Y,[],1,'omitnan')./sqrt(sum(~isnan(Y),1)), ls, ...
                    'Color', gcol(gi,:), 'LineWidth', 1.4, 'CapSize', 2); %#ok<SAGROW>
                if m == 1; lab{end+1} = sprintf('%s (n=%d)', geno_label{gi}, size(Y,1)); %#ok<SAGROW>
                elseif hf == 1; lab{end+1} = sprintf('%s puffed half (n=%d)', geno_label{gi}, size(Y,1)); %#ok<SAGROW>
                else; lab{end+1} = sprintf('%s other half', geno_label{gi}); end %#ok<SAGROW>
            end
        end
        yline(ax,0,':k','HandleVisibility','off'); xline(ax,0,'--','Color',[.8 0 0],'HandleVisibility','off')
        xlim(ax,[-180 180]); xticks(ax,-180:90:180)
        % hotspot stat p-values in the corner
        ptxt = '';
        for gi = 1:2
            x = S(walk_fly_group==gi,c,m); x = x(~isnan(x));
            if numel(x) >= 4; ptxt = [ptxt sprintf('%s hotspot p=%.2g  ', geno_label{gi}(1:4), signrank(x))]; end %#ok<AGROW>
        end
        text(ax, .02, .98, ptxt, 'Units','normalized','VerticalAlignment','top','FontSize',7)
        if m == 3; xlabel(ax,'position re: ATP hotspot (deg)'); end
        if c == 1; ylabel(ax, ylabs{m}); end
        if m == 1; title(ax, comp_names{c}, 'FontSize', 10); end
        if c == 4 && ~isempty(h); legend(ax, h, lab, 'Location','southoutside','FontSize',7); end
        row_ax(c) = ax;
    end
    linkaxes(row_ax(2:4),'y') % during-puff column is on its own scale (a transient, not a trial mean)
end
title(tl,{'Bump occupancy and EPG fluorescence across the PB, aligned to the ATP hotspot of the most recent pulse trial (red dashed)', ...
    'rows: occupancy | raw F | dF/F;  solid = puffed PB half, dashed = other half (same angles); corner p = hotspot statistic vs 0'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_hotspot_occupancy_fluorescence.png'), 'Resolution', 140)
fig.Visible = 'on';

%% local functions
function [h, ok] = pulse_hotspot(tr, base_win, peak_win, peak_factor)
    % this pulse trial's ATP hotspot (puffed half + PVA angle) and its
    % confirmed pulses' before (base_win) / after (peak_win) sample masks
    xb = tr.ft.xb(:); xf = tr.ft.xf(:);
    st = logical(tr.ft.stims(:)); t_on = xf(diff([false;st])==1);
    atp_m = mean(tr.atp.d,1)'; prof = zeros(size(tr.atp.d,1),1);
    in_b = false(size(xb)); in_a = false(size(xb)); nk = 0;
    for k = 1:numel(t_on)
        bwin = xb>=t_on(k)+base_win(1) & xb<t_on(k)+base_win(2);
        q = xb>t_on(k)+peak_win(1) & xb<=t_on(k)+peak_win(2);
        if median(atp_m(q)) - median(atp_m(bwin)) <= peak_factor*std(atp_m(bwin)); continue; end
        prof = prof + mean(tr.atp.d(:,q),2) - mean(tr.atp.d(:,bwin),2);
        in_b = in_b | bwin; in_a = in_a | q; nk = nk + 1;
    end
    ok = nk >= 3;
    [~,pk] = max(prof); puffed = (1:16) + 16*(pk>16); other = setdiff(1:32, puffed);
    w = max(prof(puffed),0); alpha = tr.im.alpha(:);
    h = struct('th', angle(sum(w .* exp(1i*alpha(puffed)))), 'puffed', puffed, 'other', other, 'in_a', in_a, 'in_b', in_b, 'pos', NaN);
end

function [occ, F, d] = profiles(tr, keep, th, edges, rho_gate, puffed, other)
    % one trial, samples 'keep' (imaging clock): occupancy re: th; per-half
    % fluorescence profiles binned by wedge angle re: th. F is raw im.f
    % (NOT trial-normalized) -- used for the peri-pulse % change.
    nb = numel(edges)-1; wrap = @(x) angle(exp(1i*x));
    mu = tr.im.mu(:); ok = keep(:) & tr.im.rho(:) > rho_gate;
    occ = histcounts(wrap(mu(ok)-th), edges, 'Normalization','probability');
    F = nan(2,nb); d = nan(2,nb);
    fm = mean(tr.im.f(:,keep),2); dm = mean(tr.im.d(:,keep),2);
    halves = {puffed, other};
    for hf = 1:2
        wv = halves{hf}; b = discretize(wrap(tr.im.alpha(wv)'-th), edges);
        F(hf,:) = accumarray(b, fm(wv), [nb 1], @mean, NaN)';
        d(hf,:) = accumarray(b, dm(wv), [nb 1], @mean, NaN)';
    end
end

function [occ, F, d] = cond_profiles(trials, th, edges, rho_gate, puffed, other)
    % mean over trials: occupancy, raw F normalized by the trial's all-wedge
    % mean (removes bleaching/global drift), dF/F
    nb = numel(edges)-1; occ = zeros(1,nb); F = zeros(2,nb); d = zeros(2,nb); n = 0;
    for k = 1:numel(trials)
        tr = trials(k); keep = true(size(tr.im.mu(:)));
        [o, Fk, dk] = profiles(tr, keep, th, edges, rho_gate, puffed, other);
        Fk = Fk / mean(tr.im.f(:)); % trial-normalized raw F
        occ = occ + o; F = F + Fk; d = d + dk; n = n + 1;
    end
    occ = occ/n; F = F/n; d = d/n;
end
