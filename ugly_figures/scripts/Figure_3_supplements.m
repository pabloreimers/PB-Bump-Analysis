%% Figure_3_supplements
% Figure 3 supplements S1-S4. Exported as vector PDFs to
% ugly_figures/all_figs/supplements/Fig3_S<N>_<name>_V1.pdf (PNG previews in
% ugly_figures/exports/supplements/):
%   S1  bump position by stim number        (its lpsp stim-number panel is also Figure 3 row 3 col 2)
%   S2  whole-PB atp by stim number
%   S3  EPG calcium, stim-side wedges, by stim number
%   S4  EPG calcium, non-stim-side wedges, by stim number
%   S5  the two Figure 3 example stims: EPG calcium heatmap with the atp-channel
%       heatmap of the SAME trial/stim directly below (documentation; was
%       Figure 3 row 2)
%
% Figure 3's average bump-position row (Figure_3_claude.m,
% "row 4: average bump-position trace") split out by stim number, to ask
% whether the bump becomes MORE mobile at successive ATP ejections.
%
% Every lpsp_p2x2_reredo trial has exactly 4 TTL stims (~60 s apart), and
% every trial with any confirmed ejection has all 4 confirmed (checked
% 2026-10-20: 99/101 trials 4/4, 2 trials 0/4), so "stim N" below is the
% trial's own Nth TTL -- recovered by matching each confirmed event's
% onset_t back to the trial's TTL onsets rather than assuming ev(N) = stim N.
%
% Differences from the Figure 3 row:
%  - BASELINE: mu is unwrapped and baseline-subtracted using the mean over
%    -2..0 s before onset (Figure 3 zeroes at the single onset sample). Less
%    sensitive to one noisy frame at t=0; otherwise the same trace.
%  - Same side flip (onto the Figure 3 example's side, reference_is_right),
%    same per-fly-then-per-genotype pooling, both light conditions pooled.
%
% Row 1: stim number WITHIN a trial (1-4), one line per stim number.
% Row 2: stim order across the fly's WHOLE SESSION (trials ordered by their
%        acquisition number in the trial folder name, e.g. "20260428-15"),
%        binned by which trial in the session (trial 1 = stims 1-4, ...).
% Right columns: per-stim summary metrics, averaged per fly, mean +/- SEM
%        over flies, vs. stim number / session trial.
%
% Two figures, same layout:
%  BUMP: side-flipped bump position. Metrics: late-window (late_win) mean
%        displacement, signed (positive = the direction the Figure 3
%        example's side pushes the bump) and unsigned (|.| per stim, THEN
%        averaged -- so opposite-direction stims don't cancel).
%  ATP:  whole-PB average atp.d dF/F (ev.atp_dff_rel, no side flip --
%        it's a whole-PB mean). Asks whether later stims simply eject more
%        ATP. Metrics: peak (max over 0-5 s of the baseline-subtracted
%        trace -- same definition as atp_peak_amp) and the raw pre-stim
%        baseline (-2..0 s), which rises if the previous stim's ATP hasn't
%        cleared ~60 s later.
%  EPG:  mean EPG calcium (im.z) over the STIM-SIDE and NON-STIM-SIDE 16
%        wedges (ev.im_z_half1/2, picked per trial_is_right exactly as
%        Figure_3_claude.m's side_pick), baseline-subtracted. Metrics: peak
%        (max 0-3 s), sustained mean (3-10 s, where lpsp and empty separate
%        in Figure 3) and the raw pre-stim baseline, to see whether calcium
%        has recovered by the next stim.
%
% REQUIRES the Figure_3_claude.m workspace (trial_events, genotype, fly_id,
% trial_is_right, reference_is_right, all_data, t_grid, ...). Run that
% script first; this one only reads its variables.

%% check workspace
need = {'example_disp','example_stim_idx','example_col_label','example_geno_label','trial_events','genotype','fly_id','trial_is_right','reference_is_right','all_data','t_grid','ref_side_str','gcamp_color','atp_color'};
missing = {}; for v = need; if exist(v{1},'var')~=1; missing{end+1} = v{1}; end; end %#ok<AGROW> -- not cellfun: exist() inside an anonymous fn checks that fn's workspace
if ~isempty(missing)
    error('run Figure_3_claude.m first -- missing: %s', strjoin(missing,', '))
end
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
export_dir = fullfile(repo_root,'ugly_figures','exports','supplements');
supp_dir   = fullfile(repo_root,'ugly_figures','all_figs','supplements');
if ~isfolder(export_dir); mkdir(export_dir); end
if ~isfolder(supp_dir);   mkdir(supp_dir);   end

%% per-stim traces with stim number and session order
pre_s = -2; post_s = 20; base_win = [-2 0]; late_win = [10 20]; peak_win = [0 5]; epg_peak_win = [0 3]; epg_sus_win = [3 10];
tr_idx  = t_grid>=pre_s & t_grid<=post_s; t_tr = t_grid(tr_idx);
base_ix = t_grid>=base_win(1) & t_grid<base_win(2);
late_ix = t_tr>=late_win(1) & t_tr<=late_win(2);
peak_ix = t_tr>=peak_win(1) & t_tr<=peak_win(2);
epg_peak_ix = t_tr>=epg_peak_win(1) & t_tr<=epg_peak_win(2);
epg_sus_ix  = t_tr>epg_sus_win(1) & t_tr<=epg_sus_win(2);

n_trials = numel(all_data);
acq_num = nan(n_trials,1); % acquisition number from the trial folder name, "YYYYMMDD-N_..."
for i = 1:n_trials
    tok = regexp(all_data(i).meta,'\d{8}-(\d+)_','tokens','once');
    if ~isempty(tok); acq_num(i) = str2double(tok{1}); end
end

S = struct('fly',{},'geno',{},'stim_num',{},'sess_trial',{},'bump',{},'atp',{},'atp_base',{},'epg_stim',{},'epg_stim_base',{},'epg_non',{},'epg_non_base',{});
[fly_list,~,fly_ix] = unique(fly_id);
for f = 1:numel(fly_list)
    tr = find(fly_ix==f);
    [~,o] = sort(acq_num(tr)); tr = tr(o); % chronological within fly
    for k = 1:numel(tr)
        i = tr(k);
        st = logical(all_data(i).ft.stims(:));
        ttl_t = all_data(i).ft.xf(diff([false;st])==1);
        ev = trial_events{i};
        for e = 1:numel(ev)
            [dmin,sn] = min(abs(ttl_t - ev(e).onset_t));
            if dmin > 1e-6; continue; end
            mu = ev(e).mu_rel;                         % unwrapped im.mu, zeroed at onset
            mu = mu - mean(mu(base_ix),'omitnan');     % -> baseline-subtracted (-2..0 s mean)
            if trial_is_right(i) ~= reference_is_right; mu = -mu; end
            a_base = mean(ev(e).atp_dff_rel(base_ix),'omitnan');
            a = ev(e).atp_dff_rel - a_base;
            if trial_is_right(i) % wedges 1:16 = right half = stim side for a RIGHT-stim trial (Figure_3_claude.m side_pick)
                zs = ev(e).im_z_half1; zn = ev(e).im_z_half2;
            else
                zs = ev(e).im_z_half2; zn = ev(e).im_z_half1;
            end
            zs_base = mean(zs(base_ix),'omitnan'); zn_base = mean(zn(base_ix),'omitnan');
            zs = zs - zs_base; zn = zn - zn_base;
            S(end+1) = struct('fly',f,'geno',genotype{i},'stim_num',sn,'sess_trial',k, ...
                'bump',mu(tr_idx),'atp',a(tr_idx),'atp_base',a_base, ...
                'epg_stim',zs(tr_idx),'epg_stim_base',zs_base,'epg_non',zn(tr_idx),'epg_non_base',zn_base); %#ok<SAGROW>
        end
    end
end
fprintf('%d confirmed stims from %d flies\n', numel(S), numel(unique([S.fly])));

%% signal definitions: trace field, two per-stim scalar metrics, labels
sig = struct();
sig(1).name = 'bump'; sig(1).fld = 'bump'; sig(1).color = gcamp_color;
sig(1).ylab = sprintf('bump pos. (rad, %s-flip, baseline-sub.)',ref_side_str);
sig(1).metric_fn = {@(x,s) mean(x(late_ix),'omitnan'), @(x,s) mean(abs(x(late_ix)),'omitnan')};
sig(1).metric_name = {'signed late displacement', 'unsigned late displacement'};
sig(1).metric_ylab = {sprintf('mean bump pos. %d-%d s (rad, signed)',late_win), sprintf('mean |bump pos.| %d-%d s (rad)',late_win)};
sig(1).shade = late_win; sig(1).zero_line = true;
sig(1).title = 'Bump-position response to ATP ejection, split by stim number';
sig(1).file = 'Fig3_S1_bump_position_by_stim_number';

sig(2).name = 'atp'; sig(2).fld = 'atp'; sig(2).color = atp_color;
sig(2).ylab = 'whole-PB atp dF/F (baseline-sub.)';
sig(2).metric_fn = {@(x,s) max(x(peak_ix)), @(x,s) s.atp_base};
sig(2).metric_name = {'peak ejection', 'pre-stim baseline'};
sig(2).metric_ylab = {sprintf('peak atp dF/F %d-%d s (baseline-sub.)',peak_win), sprintf('raw atp dF/F, %d to %d s',base_win)};
sig(2).shade = peak_win; sig(2).zero_line = false;
sig(2).title = 'Whole-PB atp response to ATP ejection, split by stim number';
sig(2).file = 'Fig3_S2_atp_by_stim_number';

side_fld = {'epg_stim','epg_non'}; side_lab = {'STIM-SIDE','NON-STIM-SIDE'}; side_file = {'stimside','nonstimside'};
for k = 1:2
    q = 2+k; f0 = side_fld{k};
    sig(q).name = f0; sig(q).fld = f0; sig(q).color = gcamp_color;
    sig(q).ylab = sprintf('EPG z, %s wedges (baseline-sub.)',lower(side_lab{k}));
    sig(q).metric_fn = {@(x,s) max(x(epg_peak_ix)), @(x,s) mean(x(epg_sus_ix),'omitnan'), @(x,s) s.([f0 '_base'])};
    sig(q).metric_name = {'peak', 'sustained', 'pre-stim baseline'};
    sig(q).metric_ylab = {sprintf('peak EPG z %d-%d s (baseline-sub.)',epg_peak_win), sprintf('mean EPG z %d-%d s (baseline-sub.)',epg_sus_win), sprintf('raw EPG z, %d to %d s',base_win)};
    sig(q).shade = epg_sus_win; sig(q).zero_line = true;
    sig(q).title = sprintf('EPG calcium (im.z), %s wedges, response to ATP ejection split by stim number',side_lab{k});
    sig(q).file = sprintf('Fig3_S%d_epg_%s_by_stim_number',q,side_file{k});
end

%% pool: per fly x bin, then across flies
geno_order = {'empty>p2x2','lpsp>p2x2'};
geno_label = {'empty>P2X2','lpsp>P2X2'};
n_stim = 4;
n_sess = 4; % session trials 1..4 (a few flies have a 5th; folded into 4 as "4+")
groupings = {'stim_num', n_stim, 'stim %d in trial', 'stim number within trial'; 'sess_trial', n_sess, 'trial %d in session', 'trial number within session'};

for si = 1:numel(sig)
    for gr = 1:2
        bfld = groupings{gr,1}; nb = groupings{gr,2};
        for gi = 1:2
            Sg = S(strcmp({S.geno},geno_order{gi}));
            flies = unique([Sg.fly]);
            nm = numel(sig(si).metric_fn);
            M = nan(numel(flies),numel(t_tr),nb); X = nan(numel(flies),nb,nm);
            for a = 1:numel(flies)
                Sf = Sg([Sg.fly]==flies(a));
                b_of = min([Sf.(bfld)], nb);
                for b = 1:nb
                    Sb = Sf(b_of==b); if isempty(Sb); continue; end
                    M(a,:,b) = mean(vertcat(Sb.(sig(si).fld)),1,'omitnan');
                    for m = 1:nm
                        X(a,b,m) = mean(arrayfun(@(q) sig(si).metric_fn{m}(q.(sig(si).fld),q), Sb),'omitnan');
                    end
                end
            end
            sig(si).P(gr,gi).M = M; sig(si).P(gr,gi).X = X;

            % stats: does the metric change with bin? Friedman (flies complete in all bins) + per-fly Spearman trend
            for m = 1:nm
                Xm = X(:,:,m); Xm = Xm(all(~isnan(Xm),2),:);
                p = nan; if size(Xm,1) >= 3; p = friedman(Xm,1,'off'); end
                r = diag(corr(repmat(1:nb,size(Xm,1),1)',Xm','type','Spearman'));
                fprintf('%-5s %-11s %-28s %-27s n=%2d  Friedman p=%.3f  median per-fly Spearman=%+.2f  bin means: %s\n', ...
                    sig(si).name, geno_label{gi}, groupings{gr,4}, sig(si).metric_name{m}, size(Xm,1), p, median(r), mat2str(mean(Xm,1),3));
            end
        end
    end
end

%% plot: one figure per signal
for si = 1:numel(sig)
    nm = numel(sig(si).metric_fn); nc = 2+nm;
    fig = figure('color','w','Visible','off','Position',[50 50 450*nc 950]);
    tl = tiledlayout(fig,2,nc,'TileSpacing','compact','Padding','compact');
    for gr = 1:2
        nb = groupings{gr,2};
        trace_axes = gobjects(1,2);
        for gi = 1:2
            base_col = [.35 .35 .35; sig(si).color]; base_col = base_col(gi,:); % empty = grays, lpsp = channel color, as in Figure 3
            cmap = white_to_color_local(base_col, nb); % light -> dark with bin
            ax = nexttile(tl,(gr-1)*nc+gi); hold(ax,'on')
            M = sig(si).P(gr,gi).M; h = gobjects(1,nb); lab = cell(1,nb);
            for b = 1:nb
                Mb = M(:,:,b); Mb = Mb(~all(isnan(Mb),2),:);
                mu_b = mean(Mb,1,'omitnan'); se_b = std(Mb,[],1,'omitnan')./sqrt(sum(~isnan(Mb),1));
                patch(ax,[t_tr,fliplr(t_tr)],[mu_b+se_b,fliplr(mu_b-se_b)],cmap(b,:),'FaceAlpha',.12,'EdgeColor','none')
                h(b) = plot(ax,t_tr,mu_b,'Color',cmap(b,:),'LineWidth',2);
                lab{b} = sprintf(groupings{gr,3},b); if gr==2 && b==nb; lab{b} = [lab{b} '+']; end
                lab{b} = sprintf('%s (n=%d)',lab{b},size(Mb,1));
            end
            xline(ax,0,':k'); yline(ax,0,':k'); xregion(ax,sig(si).shade(1),sig(si).shade(2),'FaceColor',[.9 .9 .9],'FaceAlpha',.3);
            xlim(ax,[pre_s post_s]); xlabel(ax,'time from stim onset (s)')
            if gi==1; ylabel(ax,sig(si).ylab); end
            legend(ax,h,lab,'Location','best','FontSize',8)
            title(ax,sprintf('%s: by %s',geno_label{gi},groupings{gr,4}))
            trace_axes(gi) = ax;
        end
        linkaxes(trace_axes,'y')
        for m = 1:nm
            ax = nexttile(tl,(gr-1)*nc+2+m); hold(ax,'on')
            gcol = [.5 .5 .5; sig(si).color];
            for gi = 1:2
                X = sig(si).P(gr,gi).X(:,:,m); xo = (gi-1.5)*.15;
                plot(ax,(1:nb)+xo,X','Color',[gcol(gi,:) .25],'HandleVisibility','off')
                errorbar(ax,(1:nb)+xo,mean(X,1,'omitnan'),std(X,[],1,'omitnan')./sqrt(sum(~isnan(X),1)),'-o', ...
                    'Color',gcol(gi,:),'MarkerFaceColor',gcol(gi,:),'LineWidth',2,'DisplayName',geno_label{gi})
            end
            if sig(si).zero_line; yline(ax,0,':k','HandleVisibility','off'); end
            xlim(ax,[.5 nb+.5]); xticks(ax,1:nb)
            xlabel(ax,groupings{gr,4}); ylabel(ax,sig(si).metric_ylab{m})
            legend(ax,'Location','best','FontSize',8)
            title(ax,sprintf('%s (per fly)',sig(si).metric_name{m}))
        end
    end
    title(tl,sprintf('Figure 3 - supplement S%d: %s (grey band = %g-%g s metric window)',si,sig(si).title,sig(si).shade))
    exportgraphics(fig, fullfile(export_dir,[sig(si).file '.png']), 'Resolution', 200)
    exportgraphics(fig, fullfile(supp_dir,[sig(si).file '_V1.pdf']), 'ContentType', 'vector')
    fig.Visible = 'on';
end

%% S5: example-stim heatmaps, EPG calcium over atp channel (same trial, same stim)
% Exactly the two example stims Figure 3 row 1 shows (example_disp from
% Figure_3_claude.m: same zoom window and per-window color limits), each
% column an EPG calcium (im.z) heatmap with the atp-channel (atp.d, dF/F)
% heatmap of the same trial underneath on a shared time axis -- so the
% ejection's spatial footprint can be read against the bump.
fig = figure('color','w','Visible','off','Position',[50 50 1700 900]);
tl = tiledlayout(fig,2,2,'TileSpacing','compact','Padding','compact');
for c = 1:2
    d = example_disp(c);
    axg = nexttile(tl,c); hold(axg,'on')
    imagesc(axg,d.xb,d.im_alpha,d.ex.im.z,d.im_clim)
    colormap(axg,white_to_color_local(gcamp_color,256)); cb = colorbar(axg); cb.Label.String = 'EPG z-score';
    xline(axg,d.stim_t,':k','LineWidth',1)
    xlim(axg,[d.zoom_t0,d.zoom_t1]); ylim(axg,[min(d.im_alpha),max(d.im_alpha)]); xticks(axg,[])
    if c==1; ylabel(axg,'PB angle (rad)'); end
    title(axg,sprintf('%s-stim example (%s), stim %d: EPG calcium (im.z)\n%s', example_col_label{c}, example_geno_label{c}, example_stim_idx(c), d.ex.meta),'Interpreter','none','FontSize',9)
    axa = nexttile(tl,2+c); hold(axa,'on')
    imagesc(axa,d.xb,d.atp_alpha,d.ex.atp.d,d.atp_clim)
    colormap(axa,white_to_color_local(atp_color,256)); cb = colorbar(axa); cb.Label.String = 'atp dF/F';
    xline(axa,d.stim_t,':k','LineWidth',1)
    xlim(axa,[d.zoom_t0,d.zoom_t1]); ylim(axa,[min(d.atp_alpha),max(d.atp_alpha)])
    if c==1; ylabel(axa,'PB angle (rad)'); end
    xlabel(axa,'time (s)'); title(axa,'atp channel (atp.d, dF/F), same trial','FontSize',9)
    linkaxes([axg,axa],'x')
end
title(tl,'Figure 3 - supplement S5: example stims, EPG calcium (top) with the atp-channel heatmap of the same trial (bottom); dotted line = stim onset')
exportgraphics(fig, fullfile(export_dir,'Fig3_S5_example_atp_heatmaps.png'), 'Resolution', 200)
exportgraphics(fig, fullfile(supp_dir,'Fig3_S5_example_atp_heatmaps_V1.pdf'), 'ContentType', 'vector')
fig.Visible = 'on';

%% local functions
function cmap = white_to_color_local(c, n)
    % n colors from near-white to c (skips pure white so the lightest line stays visible)
    % (n=256 for a heatmap: starts near-white at .75 mix -- use the full ramp then)
    if n >= 64; w = linspace(1,0,n)'; else; w = linspace(.75,0,n)'; end
    cmap = (1-w).*c + w.*[1 1 1];
end
