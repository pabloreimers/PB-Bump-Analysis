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

data_dir = '.data';
out_dir  = fullfile('ugly_figures','ionto_celltypes');
cl = [0.90 0.45 0.05]; ce = [0.10 0.40 0.80]; % lpsp>P2X2, empty>P2X2

%% left: reredo, 2 s pulses
L.t = -10:0.1:30; L.base = [-5 0]; L.side_win = 3;
tmp = load(fullfile(data_dir,'lpsp_p2x2_reredo_20260518.mat'),'all_data'); ad = tmp.all_data(:); clear tmp
L = collect(L,ad,@reredo_labels,inf);
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
tmp = load(fullfile(data_dir,'lpsp_p2x2_walking_20260728_post0708.mat'),'all_data'); ad = tmp.all_data(:); clear tmp
max_stims = 30; % drop the 600 s / 59-pulse trials (no dye response, see header)
R = collect(R,ad,@walking_labels,max_stims);
R.bin = ceil(R.snum/bin_size);
nb = max(R.bin);
keep_bin = arrayfun(@(b) numel(unique(R.fly(R.bin==b & R.is_lpsp)))>=2 && numel(unique(R.fly(R.bin==b & ~R.is_lpsp)))>=2, 1:nb);
nb = find(keep_bin,1,'last'); % bins past this have <2 flies of a genotype (only the 59-stim trials reach them)
fprintf('walking: %d stims, lpsp %d flies, empty %d flies, %d bins of %d (%d long trials skipped)\n',numel(R.fly), ...
    numel(unique(R.fly(R.is_lpsp))),numel(unique(R.fly(~R.is_lpsp))),nb,bin_size,R.n_skipped);

%% figure
figure(8); clf; set(gcf,'Position',[20 40 1300+110*nb 680],'Color','w')
wL = 3; wS = 3;    % grid columns for the reredo and summary panels; each bin panel gets 1 (compressed)
nc = wL+nb+wS;
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
            stim_lbl = sprintf('%d-%d',(b-1)*bin_size+1,min(b*bin_size,max(R.snum)));
            title(ax,sprintf('stims %s',stim_lbl),'FontSize',8)
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
    xlim(ax,[0.5 nb+0.5]); set(ax,'XTick',1:nb,'XTickLabel',arrayfun(@(b) sprintf('%d-%d',(b-1)*bin_size+1,min(b*bin_size,max(R.snum))),1:nb,'UniformOutput',false),'XTickLabelRotation',0)
    xlabel(ax,'stim # in trial')
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
title(tl,{'LPsP stimulation: lpsp>P2X2 (orange), empty>P2X2 (blue), lpsp - empty (black); mean \pm SEM across flies, all delivered stims', ...
    'left: lpsp\_p2x2\_reredo (Figure 3 main set)', ...
    sprintf('middle: lpsp\\_p2x2\\_walking post-7/8, 300 s trials (29 x 0.5 s pulses every 10 s; lpsp %d / empty %d flies), binned by stim number within the trial (shared y axis; \\^ = stim onset)', ...
        numel(unique(R.fly(R.is_lpsp))),numel(unique(R.fly(~R.is_lpsp)))), ...
    sprintf('%d walking trials of 600 s / 59 pulses excluded: no dye response at any stim',R.n_skipped)})
exportgraphics(gcf,fullfile(out_dir,[out_name '.pdf']))
fprintf('figure written to %s\n',out_dir);

%% functions
function D = collect(D,ad,label_fun,max_stims)
    % one row per delivered stim: baseline-subtracted green / red traces on
    % D.t, the stim's number within its trial, fly id and genotype; trials
    % with more than max_stims stims are skipped
    D.g = []; D.r = []; D.snum = []; D.fly = {}; D.is_lpsp = []; D.n_skipped = 0;
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
        if sum(ev(1:half)) >= sum(ev(half+1:end)); rows = 1:half; else; rows = half+1:2*half; end
        gs = mean(tr.im.d(rows,1:n),1)'; rs = sum(A(rows,:),1)';
        for s = 1:numel(on)
            g = interp1(xb-on(s),gs,D.t,'linear',nan); r = interp1(xb-on(s),rs,D.t,'linear',nan);
            if all(isnan(g(bidx))); continue; end
            D.g(end+1,:) = g - mean(g(bidx),'omitnan'); D.r(end+1,:) = r - mean(r(bidx),'omitnan'); %#ok<AGROW>
            D.snum(end+1,1) = s; D.fly{end+1,1} = fid; D.is_lpsp(end+1,1) = is_lpsp; %#ok<AGROW>
        end
    end
    D.is_lpsp = logical(D.is_lpsp);
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
