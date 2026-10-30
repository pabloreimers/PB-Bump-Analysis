%% fig_figure3_p2x2_traces
% LPsP stimulation (ATP onto lpsp>P2X2) vs empty>P2X2, on the two
% lpsp_p2x2 datasets that Figure_3_claude.m uses:
%   LEFT  - lpsp_p2x2_reredo_20260518.mat (EPG syt8s, 2 s ATP pulses, 4 per
%           trial; genotype from the "fly N _ lpsp/empty" folder)
%   MIDDLE - lpsp_p2x2_walking_20260728_post0708.mat (EPG 8m, walking, 0.5 s
%           pulses every 10 s, 29 or 59 per trial; genotype "_lpsp_" in the
%           trial name, as in Figure_3_claude.m). Stims are binned by their
%           number within the trial (1-5, 6-10, ...), one compressed panel
%           per bin (shared y axis, scale bar, arrowhead = stim onset).
% Signals come from the EJECTION-SIDE half of the PB only (16 glomeruli):
% each trial's side = the half with more evoked red, summed over all its
% stims. Top row: EPG green dF/F (mean of those 16 glomeruli). Bottom row:
% red (ATP dye) dF/F summed over the same 16. Far right: per fly and bin,
% the mean of the fly's bin-averaged trace over the whole 0-10 s
% inter-pulse interval, lpsp>P2X2 top and empty>P2X2 bottom, GCaMP left
% axis and dye right axis.
% lpsp>P2X2 faded orange, empty>P2X2 faded blue, lpsp - empty black
% (SEMs combined in quadrature). Each stim is baseline-subtracted
% (mean of base_win before its own onset), stims are averaged per fly, then
% mean +/- SEM across flies. Bin panels share their row's y axis.
%
% The walking set has two trial types: 199 x 300 s trials with 29 pulses
% and 57 x 600 s trials with 59 pulses. The 600 s trials show NO dye
% response at any stim (mean evoked red -0.01 in every 10-stim bin vs
% +0.65 to +0.95 in the 300 s trials; checked 2026-10-30), i.e. nothing
% visibly left the pipette, so they are excluded (max_stims below).
%
% Every delivered stim (ft.stims rising edge) is used, not only detected
% ejections: the 0.5 s walking pulses are mostly below the per-stim
% ejection detector (~15%), and the averaged dye trace shows what was
% actually ejected. In the reredo set 396/400 stims are detected ejections,
% so the left side is effectively the same either way. All flies in both
% sets are called ATP (fig_ionto_celltypes.m, fig06).
%
% Run from the repo root: addpath('ugly_figures/ionto_celltypes'); fig_figure3_p2x2_traces
% -> fig08_figure3_p2x2_lpsp_minus_empty.pdf (5-stim bins, 0-10 s mean).
% To explore, set bin_size and/or pk_metric ('mean' | 'max' = peak in 0-8 s)
% before running; non-default settings are saved with a _binN / _max suffix.
% ej_only = true keeps only CONFIRMED ejections (same detector as
% fig_ionto_celltypes.m: whole-PB summed dye, median over 0-3 s more than 3
% baseline SDs above the -5..0 s median). Walking ejections are then split
% into rounds: within each fly's same-side run of trials (restarting when
% the ejection side switches, as in fig07), round 1 = the run's first
% trial and round 2 = its second (~5 min later; later trials left out).
% Ejections are numbered 1..n within each trial and binned by bin_size;
% round-1 bins, a dotted line, then round-2 bins. Saved with an _ej suffix.
% ej_trials = true keeps every commanded stim but only from trials with >= 1
% confirmed ejection, split into round 1 | round 2 (..._ejtrials); add
% split_rounds = false to pool all such trials by stim # (..._ejtrials_pooled).

data_dir = '.data';
out_dir  = fullfile('ugly_figures','ionto_celltypes');
cl = [0.90 0.45 0.05]; ce = [0.10 0.40 0.80]; % lpsp>P2X2, empty>P2X2

%% left: reredo, 2 s pulses
L.t = -10:0.1:30; L.base = [-5 0]; L.side_win = 3;
tmp = load(fullfile(data_dir,'lpsp_p2x2_reredo_20260518.mat'),'all_data'); ad = tmp.all_data(:); clear tmp
if ~exist('ej_only','var'); ej_only = false; end
% ej_trials = true: keep ALL commanded stims, but only from trials with at
% least one confirmed ejection (saved with an _ejtrials suffix)
if ~exist('ej_trials','var'); ej_trials = false; end
% split_rounds (ej_trials mode only): true = round 1 | round 2 panels;
% false = all kept trials pooled by stim # (saved as ..._ejtrials_pooled)
if ~exist('split_rounds','var'); split_rounds = true; end
L = collect(L,ad,@reredo_labels,inf);
n_all_L = numel(L.fly);
if ej_only; L = keep_rows(L,L.ej); end
if ej_trials && ~ej_only; L = keep_rows(L,trial_has_ej(L)); end
fprintf('reredo: %d stims, lpsp %d flies, empty %d flies\n',numel(L.fly),numel(unique(L.fly(L.is_lpsp))),numel(unique(L.fly(~L.is_lpsp))));

%% right: walking post0708, 0.5 s pulses, binned by stim number within trial
R.t = -2:0.1:10; R.base = [-2 0]; R.side_win = 2; % 10 s = next pulse
if ~exist('bin_size','var'); bin_size = 5; end
% summary metric per fly x bin: 'mean' = mean of the bin-averaged trace over
% the whole 0-10 s inter-pulse interval; 'max' = its peak in 0-8 s
if ~exist('pk_metric','var'); pk_metric = 'mean'; end
out_name = 'fig08_figure3_p2x2_lpsp_minus_empty';
if bin_size ~= 5; out_name = sprintf('%s_bin%d',out_name,bin_size); end
if ~strcmp(pk_metric,'mean'); out_name = [out_name '_' pk_metric]; end
if ej_only; out_name = [out_name '_ej']; elseif ej_trials; out_name = [out_name '_ejtrials']; end
if ej_trials && ~ej_only && ~split_rounds; out_name = [out_name '_pooled']; end
tmp = load(fullfile(data_dir,'lpsp_p2x2_walking_20260728_post0708.mat'),'all_data'); ad = tmp.all_data(:); clear tmp
max_stims = 30; % drop the 600 s / 59-pulse trials (no dye response, see header)
R = collect(R,ad,@walking_labels,max_stims);
n_all = numel(R.fly);
new_trial_bin = []; % first panel of round 2 (ejection mode)
if ej_only || (ej_trials && split_rounds)
    % ROUNDS: within each fly's same-side run of trials (count restarts at a
    % side switch), round 1 = the run's first trial, round 2 = its second
    % (~5 min later); later trials are left out. Round-1 bins are followed
    % by round-2 bins, separated by a dotted line.
    %  ej_only  : confirmed ejections only, numbered 1..n WITHIN each trial
    %  ej_trials: every commanded stim (stim # within the trial), but only
    %             trials with >= 1 confirmed ejection -- rounds are counted
    %             over those trials
    if ej_only
        R = keep_rows(R,R.ej);
        fprintf('walking: %d/%d stims are confirmed ejections\n',numel(R.fly),n_all);
        R.cnt = nan(numel(R.fly),1);
        [gid,~] = findgroups(R.fly,R.tkey);
        for g = 1:max(gid)
            j = find(gid==g); [~,o] = sort(R.snum(j)); R.cnt(j(o)) = 1:numel(j);
        end
        unit = 'ej'; cnt_x = 'confirmed ejection # within trial (round 1 | round 2)';
    else
        n_tr = numel(unique(findgroups(R.fly,R.tkey)));
        R = keep_rows(R,trial_has_ej(R));
        fprintf('walking: kept %d/%d trials with >=1 confirmed ejection (%d of %d commanded stims)\n', ...
            numel(unique(findgroups(R.fly,R.tkey))),n_tr,numel(R.fly),n_all);
        R.cnt = R.snum;
        unit = 'stims'; cnt_x = 'stim # within trial (round 1 | round 2)';
    end
    [~,R.tord] = ejection_rank(R);
    nbr = ceil(max_stims/bin_size); % bins per round
    rb  = ceil(R.cnt/bin_size);
    col = (R.tord-1)*nbr + rb; col(R.tord > 2) = nan;
    fprintf('round 1: %d %s, round 2: %d, later rounds (left out): %d\n',sum(R.tord==1),unit,sum(R.tord==2),sum(R.tord>2));
    % keep a bin only if both genotypes have >= 2 flies in it
    ok = arrayfun(@(c) numel(unique(R.fly(col==c & R.is_lpsp)))>=2 && numel(unique(R.fly(col==c & ~R.is_lpsp)))>=2, 1:2*nbr);
    cols = find(ok); nb = numel(cols);
    R.bin = nan(numel(R.fly),1);
    for b = 1:nb; R.bin(col==cols(b)) = b; end
    col_round = ceil(cols/nbr); col_rb = cols - (col_round-1)*nbr;
    col_lbl = arrayfun(@(b) sprintf('%d-%d',(col_rb(b)-1)*bin_size+1,min(col_rb(b)*bin_size,max(R.cnt(R.tord==col_round(b))))),1:nb,'UniformOutput',false);
    col_ttl = arrayfun(@(b) sprintf('R%d %s',col_round(b),col_lbl{b}),1:nb,'UniformOutput',false);
    new_trial_bin = find(col_round==2,1);
else
    if ej_trials % pooled: all kept trials together, binned by stim # within the trial
        n_tr = numel(unique(findgroups(R.fly,R.tkey)));
        R = keep_rows(R,trial_has_ej(R));
        fprintf('walking: kept %d/%d trials with >=1 confirmed ejection (%d of %d commanded stims)\n', ...
            numel(unique(findgroups(R.fly,R.tkey))),n_tr,numel(R.fly),n_all);
    end
    R.cnt = R.snum; cnt_x = 'stim # in trial';
    R.bin = ceil(R.cnt/bin_size);
    nb = max(R.bin);
    keep_bin = arrayfun(@(b) numel(unique(R.fly(R.bin==b & R.is_lpsp)))>=2 && numel(unique(R.fly(R.bin==b & ~R.is_lpsp)))>=2, 1:nb);
    nb = find(keep_bin,1,'last'); % bins past this have <2 flies of a genotype
    col_lbl = arrayfun(@(b) sprintf('%d-%d',(b-1)*bin_size+1,min(b*bin_size,max(R.cnt))),1:nb,'UniformOutput',false);
    col_ttl = strcat({'stims '},col_lbl);
    col_round = ones(1,nb); col_rb = 1:nb;
end
fprintf('walking: %d stims, lpsp %d flies, empty %d flies, %d bins of %d (%d long trials skipped)\n',numel(R.fly), ...
    numel(unique(R.fly(R.is_lpsp))),numel(unique(R.fly(~R.is_lpsp))),nb,bin_size,R.n_skipped);

%% figure
figure(8); clf; set(gcf,'Position',[20 40 min(1900,1500+110*nb) 720],'Color','w') % capped at screen width: wider canvases get squeezed
wL = 3; wS = 3; wB = 2; % grid columns for the reredo, summary and slope panels; each bin panel gets 1 (compressed)
nc = wL+nb+wS+wB;
tl = tiledlayout(2,nc,'TileSpacing','compact','Padding','compact');
sig = {'g','r'}; ylab = {'EPG green \DeltaF/F (ejection-side 16 glom.)','red (ATP dye) \Sigma\DeltaF/F (ejection side)'};
axR = gobjects(2,nb);
for row = 1:2
    % left: reredo (2 tiles wide)
    ax = nexttile(tl,(row-1)*nc+1,[1 wL]); hold(ax,'on')
    nf = three_traces(ax,L.t,L.(sig{row}),true(size(L.fly)),L.is_lpsp,L.fly,cl,ce);
    xline(ax,0,':k'); yline(ax,0,':k'); xlim(ax,L.t([1 end])); ylabel(ax,ylab{row})
    if row == 1
        title(ax,sprintf('reredo: EPG syt8s, 2 s ATP pulses\nlpsp %d flies / empty %d flies',nf(1),nf(2)),'FontSize',9)
    else
        xlabel(ax,'time from stim onset (s)')
    end
    % middle: walking, one compressed panel per bin of stim numbers. Shared y
    % axis shown on the first panel only; no x axes -- a time scale bar sits on
    % the last bottom panel and an arrowhead under t = 0 marks each stim onset.
    for b = 1:nb
        ax = nexttile(tl,(row-1)*nc+wL+b); hold(ax,'on'); axR(row,b) = ax;
        nf = three_traces(ax,R.t,R.(sig{row}),R.bin==b,R.is_lpsp,R.fly,cl,ce);
        yline(ax,0,':k'); xlim(ax,R.t([1 end]))
        ax.XAxis.Visible = 'off';
        if b > 1; ax.YAxis.Visible = 'off'; end
        if row == 1
            title(ax,sprintf('%s\nn = %d/%d',col_ttl{b},nf(1),nf(2)),'FontSize',8) % flies lpsp/empty in this bin
        end
    end
    linkaxes(axR(row,:),'y')
    yl = ylim(axR(row,1)); h = diff(yl);
    for b = 1:nb   % stim-onset arrowhead just above the bottom edge
        plot(axR(row,b),0,yl(1)+0.04*h,'^k','MarkerFaceColor','k','MarkerSize',5,'Clipping','off')
    end
    if row == 2    % time scale bar on the last bin panel
        sb = 2; x0 = R.t(end)-sb; y0 = yl(1)+0.12*h;
        plot(axR(row,nb),[x0 x0+sb],[y0 y0],'k-','LineWidth',2)
        text(axR(row,nb),x0+sb/2,y0,sprintf('%g s',sb),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',8)
        text(axR(row,1),0,yl(1),'stim','HorizontalAlignment','center','VerticalAlignment','top','FontSize',7)
    end
end

% far right: per-fly peaks, lpsp>P2X2 (top) and empty>P2X2 (bottom). For
% each fly and walking stim bin, the fly's stims in that bin are averaged
% into one trace and its maximum within pk_win is taken -> one dot per fly
% per bin, GCaMP on the left axis (black), dye on the right axis (red);
% thin lines join each fly across bins, thick = mean +/- SEM across flies.
% Both panels share axis limits so the genotypes compare directly.
if strcmp(pk_metric,'mean')
    pk_win = [0 10]; w = R.t >= 0 & R.t < 10; stat = @(x) mean(x,'omitnan'); met_lbl = 'mean'; % up to, not including, the next pulse
else
    pk_win = [0 8];  w = R.t >= 0 & R.t <= 8; stat = @max; met_lbl = 'peak';
end
geno_lbl = {'lpsp>P2X2','empty>P2X2'}; geno_sel = {R.is_lpsp,~R.is_lpsp};
cdye = [0.80 0.15 0.15]; dx = 0.1;
P = cell(2,2); % {genotype, signal}: flies x bins of peak values
for gi = 1:2
    flies = unique(R.fly(geno_sel{gi}));
    for row = 1:2
        M = nan(numel(flies),nb);
        for k = 1:numel(flies)
            for b = 1:nb
                j = geno_sel{gi} & R.bin==b & strcmp(R.fly,flies{k});
                if any(j); M(k,b) = stat(mean(R.(sig{row})(j,w),1,'omitnan')); end
            end
        end
        P{gi,row} = M;
    end
end
axP = gobjects(2,1);
for gi = 1:2
    ax = nexttile(tl,(gi-1)*nc+wL+nb+1,[1 wS]); hold(ax,'on'); axP(gi) = ax;
    for row = 1:2
        if row == 1; yyaxis(ax,'left'); c = [0 0 0]; x = (1:nb)-dx; else; yyaxis(ax,'right'); c = cdye; x = (1:nb)+dx; end
        M = P{gi,row};
        plot(ax,x,M','-','Color',[c .15],'LineWidth',.5)
        scatter(ax,repmat(x,size(M,1),1),M,14,c,'filled','MarkerFaceAlpha',.4)
        n = sum(~isnan(M),1);
        errorbar(ax,x,mean(M,1,'omitnan'),std(M,0,1,'omitnan')./sqrt(n),'-o','Color',c,'MarkerFaceColor',c,'LineWidth',2,'CapSize',6)
    end
    yyaxis(ax,'left');  ylabel(ax,[met_lbl ' GCaMP \DeltaF/F (ejection side)']); ax.YColor = 'k';
    yyaxis(ax,'right'); ylabel(ax,[met_lbl ' dye \Sigma\DeltaF/F (ejection side)']); ax.YColor = cdye;
    xlim(ax,[0.5 nb+0.5]); set(ax,'XTick',1:nb,'XTickLabel',col_lbl,'XTickLabelRotation',90*(nb>8))
    xlabel(ax,cnt_x)
    if ~isempty(new_trial_bin) && new_trial_bin > 1 && new_trial_bin <= nb
        xline(ax,new_trial_bin-0.5,':k','LineWidth',1.2,'HandleVisibility','off')
    end
    title(ax,sprintf('%s: per-fly %s (%s %g-%g s), %d flies',geno_lbl{gi},met_lbl,pk_metric,pk_win,size(P{gi,1},1)),'FontSize',9)
    for row = 1:2
        fprintf('%-11s %-5s %s by bin (mean +/- SEM over flies): %s\n',geno_lbl{gi},sig{row},met_lbl, ...
            strjoin(arrayfun(@(b) sprintf('%.3f+/-%.3f',mean(P{gi,row}(:,b),'omitnan'),std(P{gi,row}(:,b),'omitnan')/sqrt(sum(~isnan(P{gi,row}(:,b))))),1:nb,'UniformOutput',false),'  '));
    end
end
% shared limits across the two genotype panels, per axis side
for side = {'left','right'}
    lims = zeros(2,2);
    for gi = 1:2; yyaxis(axP(gi),side{1}); lims(gi,:) = ylim(axP(gi)); end
    for gi = 1:2; yyaxis(axP(gi),side{1}); ylim(axP(gi),[min(0,min(lims(:,1))) max(lims(:,2))]); end
end
% far right: build-up slope from a linear mixed model on the per-fly bin
% values of the summary panels (one row per fly x stim bin; stim # = the
% bin's middle stim within the trial, so both rounds share one axis):
%   resp ~ dye + stim*lpsp + (1|fly)
% dye has one shared slope; 'stim' = empty's GCaMP change per stim at
% matched dye, 'stim:lpsp' = how much steeper lpsp is. Plotted: each
% genotype's stim slope (does it still rise at matched ejection strength?),
% with 95% CIs from the fixed-effect covariance.
x_stim = (col_rb-1)*bin_size + (bin_size+1)/2;
T = table;
for gi = 1:2
    [fi,bi] = ndgrid(1:size(P{gi,1},1),1:nb); okm = ~isnan(P{gi,1}(:)) & ~isnan(P{gi,2}(:));
    T = [T; table(strcat(geno_lbl{gi},'_',string(fi(okm))),x_stim(bi(okm))',P{gi,2}(okm),P{gi,1}(okm),repmat(double(gi==1),sum(okm),1), ...
        'VariableNames',{'fly','stim','dye','resp','lpsp'})]; %#ok<AGROW>
end
T.fly = categorical(T.fly);
lm = fitlme(T,'resp ~ dye + stim*lpsp + (1|fly)');
cn = lm.CoefficientNames;
cv = @(nm) double(strcmp(cn,nm));
c_int = double(ismember(cn,{'stim:lpsp','lpsp:stim'}));
C = {cv('stim')+c_int, cv('stim')}; % lpsp slope, empty slope
c_lbl = {'lpsp>P2X2','empty>P2X2'}; c_col = {cl,ce};
tq = tinv(0.975,lm.DFE);
ax = nexttile(tl,wL+nb+wS+1,[2 wB]); hold(ax,'on')
for k = 1:2
    est = C{k}*lm.Coefficients.Estimate; se = sqrt(C{k}*lm.CoefficientCovariance*C{k}');
    errorbar(ax,k,1e3*est,1e3*tq*se,'o','Color',c_col{k},'MarkerFaceColor',c_col{k},'MarkerSize',8,'LineWidth',2,'CapSize',10) % 10^-3 dF/F per stim
    p0 = coefTest(lm,C{k}); % H0: slope = 0
    text(ax,k+0.12,1e3*(est+tq*se),sprintf('p = %.2g',p0),'VerticalAlignment','bottom','FontSize',9)
    fprintf('LME stim slope %-12s %+.3f x10^-3 dF/F per stim [%+.3f, %+.3f], p = %.3g\n',c_lbl{k},1e3*est,1e3*(est-tq*se),1e3*(est+tq*se),p0);
end
fprintf('LME dye slope %+.4f dF/F per dye unit, p = %.3g; %d rows, %d flies\n',lm.Coefficients.Estimate(strcmp(cn,'dye')),lm.Coefficients.pValue(strcmp(cn,'dye')),height(T),numel(categories(T.fly)));
yline(ax,0,':k'); xlim(ax,[0.4 2.6]); set(ax,'XTick',1:2,'XTickLabel',c_lbl,'XTickLabelRotation',30)
ylabel(ax,'GCaMP \DeltaF/F per stim (\times10^{-3}, dye-adjusted; 95% CI)')
title(ax,{'build-up slope (mixed model)','resp ~ dye + stim\timeslpsp + (1|fly)',sprintf('(%d-stim bins per fly, %d flies)',bin_size,numel(categories(T.fly)))},'FontSize',9)
if ej_only
    ej_txt = sprintf(['confirmed ejections only (walking: %d of %d stims shown); ejections numbered within each trial, ' ...
        'round 1 / round 2 = first / second trial of a same-side run'],sum(~isnan(R.bin)),n_all);
elseif ej_trials && ~split_rounds
    ej_txt = sprintf(['all commanded stims, binned by stim # within the trial, from trials with >=1 confirmed ejection only ' ...
        '(walking: %d of %d stims; reredo: %d of %d)'],sum(~isnan(R.bin)),n_all,numel(L.fly),n_all_L);
elseif ej_trials
    ej_txt = sprintf(['all commanded stims (stim # within trial) from trials with >=1 confirmed ejection only; ' ...
        'round 1 / round 2 = first / second such trial of a same-side run (walking: %d of %d stims shown; reredo: %d of %d)'], ...
        sum(~isnan(R.bin)),n_all,numel(L.fly),n_all_L);
else
    ej_txt = 'all delivered stims, binned by stim # within the trial';
end
title(tl,{'LPsP stimulation: lpsp>P2X2 (orange), empty>P2X2 (blue), lpsp - empty (black); mean \pm SEM across flies', ej_txt, ...
    'left: lpsp\_p2x2\_reredo (Figure 3 main set)', ...
    sprintf('middle: lpsp\\_p2x2\\_walking post-7/8, 300 s trials (29 x 0.5 s pulses every 10 s; lpsp %d / empty %d flies; shared y axis; \\^ = stim onset)', ...
        numel(unique(R.fly(R.is_lpsp))),numel(unique(R.fly(~R.is_lpsp)))), ...
    sprintf('%d walking trials of 600 s / 59 pulses excluded: no dye response at any stim',R.n_skipped)})
% dotted line through both trace rows in the gap before new_trial_bin: from
% there on, ejections mostly come from a later trial (~5 min after the
% previous trial's last pulse)
if ~isempty(new_trial_bin) && new_trial_bin > 1 && new_trial_bin <= nb
    drawnow
    fp = getpixelposition(gcf);
    pa = getpixelposition(axR(1,new_trial_bin-1),true); pb = getpixelposition(axR(1,new_trial_bin),true);
    pbot = getpixelposition(axR(2,new_trial_bin),true);
    xg = ((pa(1)+pa(3)) + pb(1))/2 / fp(3);
    y0 = pbot(2)/fp(4); y1 = (pb(2)+pb(4))/fp(4);
    annotation(gcf,'line',[xg xg],[y0-0.01 y1],'LineStyle',':','LineWidth',1.2,'Color','k')
    annotation(gcf,'textbox',[xg-0.06 y0-0.045 0.12 0.03],'String','round 2: next trial (~5 min later)','EdgeColor','none', ...
        'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7) % under the bottom row
end
exportgraphics(gcf,fullfile(out_dir,[out_name '.pdf']))
fprintf('figure written to %s\n',out_dir);

%% functions
function D = collect(D,ad,label_fun,max_stims)
    % one row per delivered stim: baseline-subtracted green / red traces on
    % D.t, the stim's number within its trial, fly id and genotype; trials
    % with more than max_stims stims are skipped
    D.g = []; D.r = []; D.snum = []; D.fly = {}; D.is_lpsp = []; D.n_skipped = 0;
    D.ej = []; D.tkey = []; D.side = [];
    bidx = D.t >= D.base(1) & D.t < D.base(2);
    for i = 1:numel(ad)
        tr = ad(i);
        if ~isfield(tr,'atp') || ~isstruct(tr.atp) || ~isfield(tr.atp,'d') || isempty(tr.atp.d); continue; end
        if ~isfield(tr.ft,'stims') || isempty(tr.ft.stims) || max(tr.ft.stims) <= 0; continue; end
        [is_lpsp,fid] = label_fun(tr.meta);
        if isnan(is_lpsp); continue; end
        st = tr.ft.stims(:) > 0.5*max(tr.ft.stims);
        on = tr.ft.xf(diff([false;st])==1);
        if numel(on) > max_stims; D.n_skipped = D.n_skipped + 1; continue; end
        n  = min(numel(tr.ft.xb),size(tr.im.d,2));
        xb = tr.ft.xb(1:n); xb = xb(:);
        % ejection side for the whole trial: evoked red per glomerulus (0..side_win s
        % after onset minus baseline, positive part), summed over all the trial's
        % stims -- per-stim side is unreliable for the weak 0.5 s pulses
        A  = tr.atp.d(:,1:n); ev = zeros(size(A,1),1);
        for s = 1:numel(on)
            p = xb > on(s) & xb <= on(s)+D.side_win; b = xb >= on(s)+D.base(1) & xb < on(s);
            if any(p) && any(b); ev = ev + max(mean(A(:,p),2) - mean(A(:,b),2),0); end
        end
        half = size(A,1)/2;
        if sum(ev(1:half)) >= sum(ev(half+1:end)); rows = 1:half; side = 1; else; rows = half+1:2*half; side = 2; end
        tk = regexp(tr.meta,'(\d{8})-(\d+)','tokens','once'); tkey = str2double(tk{1})*1e4 + str2double(tk{2});
        asum = sum(A,1)'; % whole-PB dye, for the confirmed-ejection call (as in fig_ionto_celltypes.m)
        gs = mean(tr.im.d(rows,1:n),1)'; rs = sum(A(rows,:),1)';
        for s = 1:numel(on)
            g = interp1(xb-on(s),gs,D.t,'linear',nan); r = interp1(xb-on(s),rs,D.t,'linear',nan);
            if all(isnan(g(bidx))); continue; end
            D.g(end+1,:) = g - mean(g(bidx),'omitnan'); D.r(end+1,:) = r - mean(r(bidx),'omitnan'); %#ok<AGROW>
            D.snum(end+1,1) = s; D.fly{end+1,1} = fid; D.is_lpsp(end+1,1) = is_lpsp; %#ok<AGROW>
            bw = xb >= on(s)-5 & xb < on(s); pw = xb > on(s) & xb <= on(s)+3;
            bm = median(asum(bw)); bs = std(asum(bw));
            D.ej(end+1,1) = (median(asum(pw)) - bm) > 3*bs; D.tkey(end+1,1) = tkey; D.side(end+1,1) = side; %#ok<AGROW>
        end
    end
    D.is_lpsp = logical(D.is_lpsp);
end

function keep = trial_has_ej(D)
    % true for every stim whose trial (fly + trial key) has >= 1 confirmed ejection
    g = findgroups(D.fly,D.tkey);
    keep = ismember(g,unique(g(logical(D.ej))));
end

function D = keep_rows(D,m)
    % subset every per-stim field of D to the rows in logical mask m
    m = logical(m);
    for f = {'g','r','snum','fly','is_lpsp','ej','tkey','side'}
        D.(f{1}) = D.(f{1})(m,:);
    end
end

function [rk,tord] = ejection_rank(D)
    % confirmed-ejection number within each fly: trials in chronological
    % order, stims in order within a trial; restarts at 1 whenever the
    % trial's ejection side differs from the previous trial's (fig07 rule).
    % tord = which trial of that same-side run (1st, 2nd, ...) each stim is from
    rk = nan(numel(D.fly),1); tord = nan(numel(D.fly),1);
    for f = unique(D.fly)'
        j = find(strcmp(D.fly,f{1}));
        [~,~,it] = unique(D.tkey(j));
        tside = accumarray(it,D.side(j),[],@(x) x(1));
        seg = cumsum([1; diff(tside) ~= 0]);
        for s = 1:seg(end)
            js = j(seg(it) == s);
            [~,o] = sortrows([D.tkey(js) D.snum(js)]);
            rk(js(o)) = 1:numel(js);
            [~,~,tt] = unique(D.tkey(js)); tord(js) = tt;
        end
    end
end

function [v,lo,hi,tp,dots] = diff_peak(t,X,sel,is_lpsp,fly,win,n_boot)
    % maximum of (lpsp fly-mean trace - empty fly-mean trace) within win, a
    % bootstrap CI over flies, and per-lpsp-fly dots at the peak time
    w = t >= win(1) & t <= win(2); tw = t(w);
    ML = fly_matrix(X(:,w),sel & is_lpsp,fly); ME = fly_matrix(X(:,w),sel & ~is_lpsp,fly);
    d = mean(ML,1,'omitnan') - mean(ME,1,'omitnan');
    [v,k] = max(d); tp = tw(k);
    dots = ML(:,k) - mean(ME(:,k),'omitnan');
    bv = nan(n_boot,1);
    for i = 1:n_boot
        db = mean(ML(randi(size(ML,1),size(ML,1),1),:),1,'omitnan') - mean(ME(randi(size(ME,1),size(ME,1),1),:),1,'omitnan');
        bv(i) = max(db);
    end
    lo = prctile(bv,2.5); hi = prctile(bv,97.5);
end

function M = fly_matrix(X,idx,fly)
    % one row per fly: that fly's mean trace over the selected stims
    f = unique(fly(idx)); M = nan(numel(f),size(X,2));
    for k = 1:numel(f); M(k,:) = mean(X(idx & strcmp(fly,f{k}),:),1,'omitnan'); end
end

function nf = three_traces(ax,t,X,sel,is_lpsp,fly,cl,ce)
    % per-fly means within sel, then lpsp / empty mean +/- SEM and their difference
    [mL,sL,nL] = fly_avg(X,sel & is_lpsp,fly); [mE,sE,nE] = fly_avg(X,sel & ~is_lpsp,fly);
    % genotype traces faded toward white (thin, light SEM) so the black
    % lpsp - empty difference stands out
    fade = @(c) 0.45*c + 0.55;
    shade(ax,t,mE,sE,fade(ce),0.8,0.12); shade(ax,t,mL,sL,fade(cl),0.8,0.12)
    shade(ax,t,mL-mE,sqrt(sL.^2+sE.^2),[0 0 0],2.2,0.25)
    nf = [nL nE];
end

function [m,se,nf] = fly_avg(X,idx,fly)
    f = unique(fly(idx)); nf = numel(f); M = nan(nf,size(X,2));
    for k = 1:nf; M(k,:) = mean(X(idx & strcmp(fly,f{k}),:),1,'omitnan'); end
    m = mean(M,1,'omitnan'); se = std(M,0,1,'omitnan')/sqrt(max(nf,1));
end

function shade(ax,t,m,se,col,lw,fa)
    if nargin < 7; fa = .2; end
    ok = ~isnan(m); se(isnan(se)) = 0;
    fill(ax,[t(ok) fliplr(t(ok))],[m(ok)+se(ok) fliplr(m(ok)-se(ok))],col,'FaceAlpha',fa,'EdgeColor','none')
    plot(ax,t,m,'Color',col,'LineWidth',lw)
end

function [is_lpsp,fid] = reredo_labels(meta)
    % genotype from the "fly N _ <empty|lpsp>" folder (as Figure_3_claude.m); fly = date + fly number
    p = strsplit(meta,'\'); p(cellfun(@isempty,p)) = [];
    fi = find(startsWith(p,'fly'),1,'last');
    fly_folder = lower(p{fi});
    if contains(fly_folder,'lpsp'); is_lpsp = 1; elseif contains(fly_folder,'empty'); is_lpsp = 0; else; is_lpsp = nan; end
    n = regexp(fly_folder,'^fly\s*(\d+)','tokens','once');
    fid = sprintf('%s|fly %s',p{fi-1},n{1});
end

function [is_lpsp,fid] = walking_labels(meta)
    % genotype "_lpsp_" in the trial name (as Figure_3_claude.m); fly = date + fly number
    is_lpsp = double(contains(meta,'_lpsp_'));
    p = strsplit(meta,'\'); p(cellfun(@isempty,p)) = [];
    fi = find(startsWith(p,'fly'),1,'last');
    n = regexp(lower(p{fi}),'^fly\s*(\d+)','tokens','once');
    fid = sprintf('%s|fly %s',p{fi-1},n{1});
end
