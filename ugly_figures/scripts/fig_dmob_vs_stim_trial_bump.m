%% fig_dmob_vs_stim_trial_bump
% Scratch: does the per-fly increase in dark bump mobility (Figure 3's
% walk_delta, post - pre) track what the bump did DURING that fly's ATP
% pulse trial? Walking dataset (lpsp_p2x2_walking_20260728_post0708).
%
% The pulse trial is the fly's FIRST detected perturbation trial
% (walk_has_pulse; the one that splits pre/post in Figure 3). In this
% dataset every pulse trial is CLOSED LOOP, 300 s, 29 TTL pulses 10 s
% apart -- so "bump position" is read as the bump-cue OFFSET
% (mu - (-cue), wrapped), and per-pulse windows are short: baseline
% -2..0 s, response 2..8 s (ends before the next pulse).
%
% Per pulse, kept only if the whole-PB atp.d mean confirms an ejection
% (median 0-3 s minus median -2..0 s > 3 baseline SDs -- the
% get_trial_stim_events / walk_trial_has_pulse test, per pulse):
%   d_off      : circular change in mean offset, response - baseline
%   |d_off|    : its magnitude (how far the pulse knocked the bump off the cue)
%   toward     : d_off signed + if it moved the bump TOWARD the ATP hotspot
%   hot_dist   : |circular distance| bump -> hotspot at onset (-0.5..0 s)
% Hotspot = the wedge with the largest trial-averaged 0-3 s atp.d profile
% over the confirmed pulses, as an angle in the same frame as im.mu
% (im.alpha mod 2pi -- mu is the PVA over those same wedge angles).
% Per trial:
%   r_stim     : circ corr(mu, -cue) in the pulse trial (rho > .2 and cue
%                moving > .2 rad/s, as Figure_3_scratch.m trial_tracking_r)
%   d_r        : r_stim minus the mean r over the fly's pre-ATP CL trials
%   off_sd     : circular SD of the offset across the whole pulse trial
% Each per-fly value is the mean over confirmed pulses. Spearman per
% genotype vs walk_delta.
%
% REQUIRES the Figure_3_claude.m workspace (walk_* variables, circ_stats on path).

%% check workspace
need = {'walk_data','walk_fly_ix','walk_fly_list','walk_trial_num','walk_has_pulse','walk_is_dark','walk_fly_group','walk_delta','gcamp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW>
if ~isempty(missing); error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', ')); end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports');
wrap = @(x) angle(exp(1i*x));

%% per-fly metrics from the pulse trial
base_win = [-2 0]; resp_win = [2 8]; peak_win = [0 3]; onset_win = [-.5 0]; peak_factor = 3;
nF = numel(walk_fly_list);
M = struct('abs_doff',nan(nF,1),'toward',nan(nF,1),'hot_dist',nan(nF,1),'r_stim',nan(nF,1),'d_r',nan(nF,1),'off_sd',nan(nF,1),'n_pulse',zeros(nF,1));
for f = 1:nF
    if isnan(walk_delta(f)); continue; end
    ti = find(walk_fly_ix==f); [~,o] = sort(walk_trial_num(ti)); ti = ti(o);
    p = find(walk_has_pulse(ti),1); i = ti(p);
    tr = walk_data(i); xb = tr.ft.xb(:); xf = tr.ft.xf(:);
    mu  = tr.im.mu(:);
    cue = wrap(interp1(xf, unwrap(-tr.ft.cue(:)), xb, 'linear', 'extrap')); % heading (-cue) on the imaging clock
    off = wrap(mu - cue);
    st = logical(tr.ft.stims(:)); t_on = xf(diff([false;st])==1);
    atp_m = mean(tr.atp.d,1)';

    % confirmed pulses + trial hotspot
    keep = false(numel(t_on),1); prof = zeros(size(tr.atp.d,1),1);
    for k = 1:numel(t_on)
        b = xb>=t_on(k)+base_win(1) & xb<t_on(k)+base_win(2);
        q = xb>t_on(k)+peak_win(1) & xb<=t_on(k)+peak_win(2);
        keep(k) = median(atp_m(q)) - median(atp_m(b)) > peak_factor*std(atp_m(b));
        if keep(k); prof = prof + mean(tr.atp.d(:,q),2) - mean(tr.atp.d(:,b),2); end
    end
    if sum(keep) < 3; continue; end
    [~,hw] = max(prof); hot = wrap(tr.im.alpha(hw));

    d_off = nan(numel(t_on),1); toward = d_off; hot_dist = d_off;
    for k = find(keep)'
        b = xb>=t_on(k)+base_win(1) & xb<t_on(k)+base_win(2);
        r = xb>=t_on(k)+resp_win(1) & xb<=t_on(k)+resp_win(2);
        z = xb>=t_on(k)+onset_win(1) & xb<=t_on(k)+onset_win(2);
        d_off(k) = wrap(circ_mean(off(r)) - circ_mean(off(b)));
        mu0 = circ_mean(mu(z));
        hot_dist(k) = abs(wrap(hot - mu0));
        toward(k) = d_off(k) * sign(wrap(hot - mu0)); % + = offset moved toward the hotspot
    end
    M.abs_doff(f) = mean(abs(d_off),'omitnan'); M.toward(f) = mean(toward,'omitnan');
    M.hot_dist(f) = mean(hot_dist,'omitnan'); M.n_pulse(f) = sum(keep);
    M.off_sd(f) = circ_std(off(tr.im.rho(:) > .2));

    % tracking in the pulse trial vs pre-ATP closed loop
    M.r_stim(f) = cl_track_r(tr);
    pre_cl = ti(1:p-1); pre_cl = pre_cl(~walk_is_dark(pre_cl));
    r_pre = arrayfun(@(j) cl_track_r(walk_data(j)), pre_cl);
    M.d_r(f) = M.r_stim(f) - mean(r_pre,'omitnan');
end

%% correlate with delta mobility, per genotype
metrics = {'abs_doff','toward','hot_dist','r_stim','d_r','off_sd'};
mlab = {'|\Delta offset| per pulse (rad)','\Delta offset toward ATP hotspot (rad)','bump-hotspot distance at onset (rad)', ...
        'bump-cue tracking r, pulse trial','\Delta tracking r (pulse trial - pre CL)','offset circular SD, pulse trial'};
geno_label = {'empty>P2X2','lpsp>P2X2'}; gcol = [.5 .5 .5; gcamp_color];
fig = figure('color','w','Visible','off','Position',[50 50 1500 900]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');
fprintf('%-10s %-11s %3s %8s %7s\n','metric','genotype','n','Spearman','p');
for m = 1:numel(metrics)
    ax = nexttile(tl); hold(ax,'on'); h = gobjects(1,2); leg = cell(1,2);
    x = M.(metrics{m});
    for gi = 1:2
        ok = walk_fly_group==gi & ~isnan(x) & ~isnan(walk_delta);
        [rs,ps] = corr(x(ok), walk_delta(ok), 'type','Spearman');
        h(gi) = scatter(ax, x(ok), walk_delta(ok), 50, gcol(gi,:), 'filled', 'MarkerFaceAlpha', .7);
        leg{gi} = sprintf('%s: r_s %+.2f, p=%.2g (n=%d)', geno_label{gi}, rs, ps, sum(ok));
        fprintf('%-10s %-11s %3d %+8.2f %7.3f\n', metrics{m}, geno_label{gi}, sum(ok), rs, ps);
    end
    yline(ax,0,':k','HandleVisibility','off')
    xlabel(ax, mlab{m}); ylabel(ax,'\Delta bump mobility (dark, post - pre)')
    legend(ax,h,leg,'Location','southoutside','FontSize',8)
end
title(tl,'Per-fly increase in dark bump mobility vs. what the bump did during the ATP pulse trial (closed loop, 29 pulses / 10 s; confirmed pulses only)','FontSize',10)
exportgraphics(fig, fullfile(export_dir,'fig_dmob_vs_stim_trial_bump.png'), 'Resolution', 150)
fig.Visible = 'on';

%% local functions
function r = cl_track_r(tr)
    % circ corr of bump vs -cue, rho > .2 and cue moving > .2 rad/s (Figure_3_scratch.m trial_tracking_r)
    xb = tr.ft.xb(:); cue_u = interp1(tr.ft.xf(:), unwrap(-tr.ft.cue(:)), xb, 'linear', 'extrap');
    sp = abs(gradient(smoothdata(cue_u,'movmean',5), xb));
    good = tr.im.rho(:) > .2 & sp > .2;
    r = nan; if sum(good) > 100; r = circ_corrcc(mod(cue_u(good),2*pi), mod(tr.im.mu(good),2*pi)); end
end
