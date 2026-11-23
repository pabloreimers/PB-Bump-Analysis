%% fig_hotspot_occupancy_fluorescence_ionto
% Scratch: the fig_hotspot_occupancy_fluorescence.m question on Figure 3's
% IONTOPHORESIS dataset (lpsp_p2x2_reredo_20260518, rows 1-2 of Figure 3):
% does ATP change (1) where the bump sits and (2) EPG fluorescence AT THE
% PUFF LOCATION, during and after each ejection? Here the visual cue is a
% bright bar moving in OPEN LOOP (checked: fly yaw vs -cue velocity r =
% 0.00), or the trial is dark (background pattern). Every trial is 300 s
% with 4 ATP stims at 41/101/161/221 s, so there are no separate pre/post
% trials: the baseline is each trial's own pre-stim period.
%
% HOTSPOT, per trial: average that trial's confirmed stims' (whole-PB
% atp.d median 0-3 s minus median -2..0 s > 3 baseline SDs; >= min_stims
% of them) baseline-subtracted 0-3 s atp.d profile; the PB half holding its
% peak is the PUFFED half; theta_hot = population-vector angle of that
% half's positive profile over im.alpha (the frame im.mu lives in). Each
% trial is aligned to its own hotspot.
% WINDOWS (s after each confirmed stim), each vs the trial's BASELINE
% (base_win, before the first stim):
%   during 0-3 | early after 3-10 | late after 10-58 (the settled part of
%   the 60 s inter-stim interval)
% MEASURES in n_bin angle bins re: theta_hot:
%   occupancy : fraction of samples with the bump (rho > rho_gate) there;
%               bump position = im.mu, the 32-wedge (both-halves) PVA
%   raw F     : window-mean im.f per wedge / that window's mean over all 32
%               wedges (only LOCAL change; bleaching/global drift removed)
%   dF/F      : window-mean im.d per wedge
%   Fluorescence kept per half: puffed vs other (same angles) -- the bump
%   lights both halves at the same angle, so puffed-minus-other isolates a
%   local ATP effect from a change in bump occupancy.
% STATS at the hotspot (per fly: mean over its trials of that lighting):
%   occupancy within +/-occ_half; fluorescence within +/-f_half, puffed
%   minus other half. Signrank vs 0 per genotype, ranksum between.
%
% REQUIRES the Figure_3_claude.m workspace (all_data, genotype, fly_id).

%% check workspace
need = {'all_data','genotype','fly_id','gcamp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');

%% parameters
n_bin = 16; edges = linspace(-pi,pi,n_bin+1); ctr = edges(1:end-1)+pi/n_bin;
conf_base = [-2 0]; conf_peak = [0 3]; peak_factor = 3; min_stims = 2; rho_gate = .2;
base_win = [1 39]; % s, before the first stim (41 s)
wins = {[0 3],[3 10],[10 58]}; win_names = {'during (0-3 s)','early after (3-10 s)','late after (10-58 s)'};
light_names = {'open-loop bar','dark'};
occ_half = pi/4; f_half = pi/8;
meas_names = {'bump occupancy','raw F (wedge / all-wedge mean)','dF/F'};

%% per trial -> per fly
geno_order = {'empty>p2x2','lpsp>p2x2'}; geno_label = {'empty>P2X2','lpsp>P2X2'};
[fly_list,~,fly_ix] = unique(fly_id); nF = numel(fly_list);
fly_geno = nan(nF,1);
n_trial = numel(all_data); is_dark = arrayfun(@(x) contains(x.ft.pattern,'background'), all_data);
% T{m}(trial,bin,window,half)
T = repmat({nan(n_trial,n_bin,3,2)},1,3); used = false(n_trial,1);
for i = 1:n_trial
    tr = all_data(i);
    [h, ok, t_on] = trial_hotspot(tr, conf_base, conf_peak, peak_factor, min_stims);
    if ~ok; continue; end
    xb = tr.ft.xb(:);
    kb = xb >= base_win(1) & xb <= base_win(2);
    [oB, FB, dB] = window_profiles(tr, kb, h, edges, rho_gate);
    for w = 1:3
        kw = false(size(xb));
        for k = 1:numel(t_on); kw = kw | (xb > t_on(k)+wins{w}(1) & xb <= t_on(k)+wins{w}(2)); end
        [oW, FW, dW] = window_profiles(tr, kw, h, edges, rho_gate);
        T{1}(i,:,w,1) = oW - oB;
        T{2}(i,:,w,:) = reshape((FW - FB).', 1, n_bin, 1, 2);
        T{3}(i,:,w,:) = reshape((dW - dB).', 1, n_bin, 1, 2);
    end
    used(i) = true;
end
% per fly x lighting: mean over trials
Dm = repmat({nan(nF,n_bin,3,2,2)},1,3); % (fly,bin,window,half,lighting)
for f = 1:nF
    rows = find(fly_ix==f);
    fly_geno(f) = find(strcmp(geno_order, genotype{rows(1)}));
    for L = 1:2
        r = rows(used(rows) & is_dark(rows) == (L==2));
        if isempty(r); continue; end
        for m = 1:3; Dm{m}(f,:,:,:,L) = mean(T{m}(r,:,:,:),1,'omitnan'); end
    end
end
fprintf('%d/%d trials used (>= %d confirmed stims); flies: %d empty, %d lpsp\n', sum(used), n_trial, min_stims, sum(fly_geno==1), sum(fly_geno==2));

%% stats at the hotspot
hot_occ = abs(ctr) <= occ_half; hot_f = abs(ctr) <= f_half;
S = nan(nF,3,2,3); % fly x window x lighting x measure
for w = 1:3
    for L = 1:2
        S(:,w,L,1) = mean(Dm{1}(:,hot_occ,w,1,L),2,'omitnan');
        for m = 2:3
            S(:,w,L,m) = mean(Dm{m}(:,hot_f,w,1,L),2,'omitnan') - mean(Dm{m}(:,hot_f,w,2,L),2,'omitnan');
        end
    end
end
fprintf('hotspot statistic: occupancy = change within +/-%d deg; fluorescence = change in puffed half minus other half within +/-%.1f deg; vs trial baseline %d-%d s\n', round(occ_half*180/pi), f_half*180/pi, base_win);
for m = 1:3
    for L = 1:2
        for w = 1:3
            line = sprintf('%-27s %-13s %-21s', meas_names{m}, light_names{L}, win_names{w});
            for gi = 1:2
                x = S(fly_geno==gi,w,L,m); x = x(~isnan(x)); p = NaN; if numel(x) >= 4; p = signrank(x); end
                line = [line sprintf(' | %s n=%2d med %+.3f p=%.3f', geno_label{gi}(1:4), numel(x), median(x), p)]; %#ok<AGROW>
            end
            x1 = S(fly_geno==1,w,L,m); x2 = S(fly_geno==2,w,L,m);
            fprintf('%s | between p=%.3f\n', line, ranksum(x1(~isnan(x1)), x2(~isnan(x2))));
        end
    end
end

%% plot: rows = measures, columns = lighting x window
gcol = [.5 .5 .5; gcamp_color]; deg = ctr*180/pi;
fig = figure('color','w','Visible','off','Position',[20 20 2200 1100]);
tl = tiledlayout(fig,3,6,'TileSpacing','compact','Padding','compact');
ylabs = {'\Delta occupancy (fraction)','\Delta raw F (wedge / all-wedge mean)','\Delta dF/F'};
for m = 1:3
    row_ax = gobjects(1,6);
    for L = 1:2
        for w = 1:3
            col = (L-1)*3 + w; ax = nexttile(tl); hold(ax,'on'); h = gobjects(0); lab = {};
            halves = 1; if m > 1; halves = 1:2; end
            for gi = 1:2
                g = fly_geno==gi;
                for hf = halves
                    Y = Dm{m}(g,:,w,hf,L); Y = Y(~all(isnan(Y),2),:); if isempty(Y); continue; end
                    ls = '-'; if hf == 2; ls = '--'; end
                    h(end+1) = errorbar(ax, deg + (gi-1.5)*3, mean(Y,1,'omitnan'), std(Y,[],1,'omitnan')./sqrt(sum(~isnan(Y),1)), ls, ...
                        'Color', gcol(gi,:), 'LineWidth', 1.3, 'CapSize', 2); %#ok<SAGROW>
                    if m == 1; lab{end+1} = sprintf('%s (n=%d)', geno_label{gi}, size(Y,1)); %#ok<SAGROW>
                    elseif hf == 1; lab{end+1} = sprintf('%s puffed half (n=%d)', geno_label{gi}, size(Y,1)); %#ok<SAGROW>
                    else; lab{end+1} = sprintf('%s other half', geno_label{gi}); end %#ok<SAGROW>
                end
            end
            yline(ax,0,':k','HandleVisibility','off'); xline(ax,0,'--','Color',[.8 0 0],'HandleVisibility','off')
            xlim(ax,[-180 180]); xticks(ax,-180:90:180)
            ptxt = '';
            for gi = 1:2
                x = S(fly_geno==gi,w,L,m); x = x(~isnan(x));
                if numel(x) >= 4; ptxt = [ptxt sprintf('%s p=%.2g  ', geno_label{gi}(1:4), signrank(x))]; end %#ok<AGROW>
            end
            x1 = S(fly_geno==1,w,L,m); x2 = S(fly_geno==2,w,L,m);
            ptxt = [ptxt sprintf('| e vs l p=%.2g', ranksum(x1(~isnan(x1)), x2(~isnan(x2))))]; %#ok<AGROW>
            text(ax, .02, .98, ptxt, 'Units','normalized','VerticalAlignment','top','FontSize',7)
            if m == 3; xlabel(ax,'position re: ATP hotspot (deg)'); end
            if col == 1; ylabel(ax, ylabs{m}); end
            if m == 1; title(ax, sprintf('%s: %s', light_names{L}, win_names{w}), 'FontSize', 9); end
            if col == 6 && ~isempty(h); legend(ax, h, lab, 'Location','southoutside','FontSize',7); end
            row_ax(col) = ax;
        end
    end
    linkaxes(row_ax,'y')
end
title(tl,{'Iontophoresis dataset (Figure 3 rows 1-2): bump occupancy and EPG fluorescence across the PB re: each trial''s ATP hotspot (red dashed), vs the trial''s pre-stim baseline (1-39 s)', ...
    'rows: occupancy | raw F | dF/F;  solid = puffed PB half, dashed = other half;  corner: hotspot statistic vs 0 per genotype (signrank) and empty vs lpsp (ranksum)'},'FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_hotspot_occupancy_fluorescence_ionto.png'), 'Resolution', 130)
fig.Visible = 'on';

%% local functions
function [h, ok, t_on] = trial_hotspot(tr, base_win, peak_win, peak_factor, min_stims)
    % this trial's ATP hotspot (puffed half + PVA angle) from its confirmed stims; t_on = confirmed onsets
    xb = tr.ft.xb(:); xf = tr.ft.xf(:);
    st = logical(tr.ft.stims(:)); t_all = xf(diff([false;st])==1);
    atp_m = mean(tr.atp.d,1)'; prof = zeros(size(tr.atp.d,1),1); keep = false(numel(t_all),1);
    for k = 1:numel(t_all)
        b = xb>=t_all(k)+base_win(1) & xb<t_all(k)+base_win(2);
        q = xb>t_all(k)+peak_win(1) & xb<=t_all(k)+peak_win(2);
        keep(k) = median(atp_m(q)) - median(atp_m(b)) > peak_factor*std(atp_m(b));
        if keep(k); prof = prof + mean(tr.atp.d(:,q),2) - mean(tr.atp.d(:,b),2); end
    end
    t_on = t_all(keep); ok = sum(keep) >= min_stims;
    [~,pk] = max(prof); puffed = (1:16) + 16*(pk>16); other = setdiff(1:32, puffed);
    w = max(prof(puffed),0); alpha = tr.im.alpha(:);
    h = struct('th', angle(sum(w .* exp(1i*alpha(puffed)))), 'puffed', puffed, 'other', other);
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
