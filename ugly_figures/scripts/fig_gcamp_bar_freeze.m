%% fig_gcamp_bar_freeze
% End-of-trial cue freeze in the two EPG>GCaMP datasets used by
% Figure_1_claude.m, and how well the frozen cue holds the bump still while
% the fly keeps turning.
%
%   .data/epg_7f_20261007.mat    gain_change rig, EPG>GCaMP7f, 300 s trials
%                                (one 600 s), bar (closed loop, gain 0.7) or dark
%   .data/hackathon_20250729.mat hackathon rig, EPG>syt8m. Only the G4-pattern
%                                trials (4 px bright bar, 8 px vertical
%                                grating, dark background) have the freeze --
%                                the berg4 gain-staircase trials keep
%                                integrating the cue to the very end (checked:
%                                cue changes up to the last sample, and no
%                                >5 s frozen-cue run with turning anywhere
%                                except two 10 s holds in one trial), so they
%                                are excluded here.
%
% ==FREEZE DETECTION==
% Freeze onset = time of the last cue step larger than freeze_step_rad.
% The epg7f cue is despiked first (same despike_cue as Figure_1_claude.m:
% the fictrac-socket cuePos log has isolated 1-2 sample jumps, which would
% otherwise look like late cue movement). Dark trials have the same freeze
% in their logged cue, but there is nothing on the screen to freeze, so they
% serve as the control: there the bump should keep integrating rotation.
%
% =="HOLD" METRIC==
% Per trial and window (pre = the win_pre_s seconds before the freeze, freeze
% = freeze onset + settle_s to trial end), bin the imaging frames into
% bin_s-second bins (circular mean of im.mu over non-flash frames with
% im.rho > rho_thresh) and regress the bin-to-bin bump displacement on the
% bin-to-bin fly heading change (through the origin). Heading is taken
% lag_s earlier than the bump bin (GCaMP/bump lag behind behavior; lag_s is
% picked per dataset to maximize the pooled closed-loop pre-freeze slope).
% The slope's sign is normalized per dataset so closed-loop pre-freeze is
% positive. Slope ~ closed-loop gain before the freeze; ~0 during the freeze
% if the frozen cue holds the bump; ~1 in the dark (path integration).
% The same regression with the cue in place of the bump confirms the cue
% itself stopped moving.

%% paths + constants
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');

freeze_step_rad = 0.02; % rad, a cue step larger than this counts as cue movement
win_pre_s       = 30;   % s, pre-freeze comparison window
settle_s        = 1;    % s, skip the first second of the freeze
bin_s           = 0.5;  % s, bin width for the displacement regression
rho_thresh      = 0.2;  % im.rho, discernable bump (as Figure_1_claude.m)
flash_mad_thresh = 8;
min_turn_rad    = pi;   % need at least this much total turning in a window for a slope
min_pairs       = 10;   % and at least this many valid bin pairs
lag_grid_s      = 0:0.1:1;
max_cue_dev_rad = 0.5;  % rad, drop trials whose cue strays further than this during the freeze

%% load + standardize both datasets into one trial list
T = struct('dataset',{},'fly',{},'cond',{},'xf',{},'cue',{},'heading',{},'xb',{},'mu',{},'rho',{},'keep',{},'label',{});

S = load(fullfile(data_dir,'epg_7f_20261007.mat'),'all_data'); src = S.all_data(:);
for i = 1:numel(src)
    s = src(i);
    parts = split_path(s.meta.flyDir);
    tr.dataset = 'epg7f';
    if startsWith(parts{end},'fly')
        tr.fly = sprintf('%s_%s',parts{end-1},strrep(parts{end},' ',''));
    else
        tr.fly = parts{end};
    end
    tr.cond    = cond_from_pattern(s.ft.pattern);
    tr.xf      = s.ft.xf(:);
    tr.cue     = despike_cue(s.ft.cue_raw(:));
    tr.heading = unwrap(s.ft.heading_raw(:));
    tr.xb      = s.ft.xb(:);
    tr.mu      = s.im.mu(:);
    tr.rho     = s.im.rho(:);
    bad = false(numel(tr.mu),1);
    if isfield(s.im,'badVol') && ~isempty(s.im.badVol), bad = logical(s.im.badVol(:)); end
    tr.keep    = ~(reshape(detect_flash_frames(s.im.f,flash_mad_thresh),[],1) | bad);
    tr.label   = sprintf('epg7f #%d (%s)',i,tr.fly);
    T(end+1) = tr; %#ok<SAGROW>
end

S = load(fullfile(data_dir,'hackathon_20250729.mat'),'all_data'); src = S.all_data(:);
n_berg4 = 0;
for i = 1:numel(src)
    s = src(i);
    if strcmpi(s.ft.pattern,'berg4'), n_berg4 = n_berg4+1; continue, end
    % bump recomputed from smoothed im.f, same recipe as Figure_1_claude.m / hackathon_claude.m
    f   = smoothdata(s.im.f,2,'movmean',5);
    f0  = prctile(f,7,2);
    z   = zscore((f-f0)./f0,[],2);
    [x,y]    = pol2cart(s.im.alpha(:)',z');
    [mu,rho] = cart2pol(mean(x,2),mean(y,2));
    tok = regexp(s.meta,'\d+_[a-z]{2}_\d+','match','once');
    if isempty(tok), parts = split_path(s.meta); tok = parts{end-2}; end
    tr.dataset = 'hackathon';
    tr.fly     = tok;
    tr.cond    = cond_from_pattern(s.ft.pattern);
    tr.xf      = s.ft.xf(:);
    tr.cue     = s.ft.cue(:);
    tr.heading = unwrap(s.ft.heading(:));
    tr.xb      = get_xb(s.ft,numel(mu));
    tr.mu      = mu;
    tr.rho     = rho;
    tr.keep    = ~reshape(detect_flash_frames(s.im.f,flash_mad_thresh),[],1);
    tr.label   = sprintf('hackathon #%d (%s)',i,tok);
    T(end+1) = tr; %#ok<SAGROW>
end
fprintf('%d trials (%d berg4 hackathon trials excluded: no end-of-trial freeze)\n',numel(T),n_berg4);

%% detect the freeze in every trial
for k = 1:numel(T)
    c  = T(k).cue; t = T(k).xf;
    dc = abs(angle(exp(1i*diff(c))));
    j  = find(dc > freeze_step_rad,1,'last');
    T(k).t_freeze   = t(j+1);
    T(k).freeze_dur = t(end) - t(j+1);
    in = t > t(j+1);
    T(k).freeze_cue_dev = max(abs(angle(exp(1i*(c(in)-c(j+1))))));      % how far the cue strays during the freeze
    T(k).freeze_turn    = sum(abs(diff(T(k).heading(in))));              % total turning while frozen
end

fprintf('\n%-30s %-8s %8s %8s %10s %10s\n','trial','cond','onset s','dur s','cue dev','turned rad');
for k = 1:numel(T)
    fprintf('%-30s %-8s %8.1f %8.1f %10.3f %10.1f\n',T(k).label,T(k).cond,T(k).t_freeze,T(k).freeze_dur,T(k).freeze_cue_dev,T(k).freeze_turn);
end
% a few hackathon trials have a cue that creeps by sub-threshold steps
% during the "freeze" (hackathon #32 is the clear case: onset found at
% 347.5 s of a 600 s trial, cue strays 1.06 rad) -- not a clean freeze, so drop them
not_clean = [T.freeze_cue_dev] > max_cue_dev_rad;
fprintf('excluding %d trials whose cue strays > %.1f rad during the freeze: %s\n', sum(not_clean), max_cue_dev_rad, strjoin({T(not_clean).label},', '));
T = T(~not_clean);

%% per-dataset lag (pooled closed-loop pre-freeze bump slope)
datasets = {'epg7f','hackathon'};
lag_s = nan(1,2);
for di = 1:2
    sel = find(strcmp({T.dataset},datasets{di}) & ~strcmp({T.cond},'dark'));
    sw = nan(numel(lag_grid_s),1);
    for L = 1:numel(lag_grid_s)
        sl = nan(numel(sel),1);
        for q = 1:numel(sel)
            tr = T(sel(q));
            sl(q) = window_slope(tr,tr.t_freeze-win_pre_s,tr.t_freeze,lag_grid_s(L),'bump',bin_s,rho_thresh,min_turn_rad,min_pairs);
        end
        sw(L) = abs(median(sl,'omitnan'));
    end
    [~,b] = max(sw); lag_s(di) = lag_grid_s(b);
    fprintf('%s: bump lag %.1f s\n',datasets{di},lag_s(di));
end

%% slopes: bump and cue vs heading, pre vs freeze
for k = 1:numel(T)
    tr = T(k); di = find(strcmp(datasets,tr.dataset));
    w_pre = [tr.t_freeze-win_pre_s, tr.t_freeze];
    w_frz = [tr.t_freeze+settle_s,  tr.xf(end)];
    T(k).bump_pre = window_slope(tr,w_pre(1),w_pre(2),lag_s(di),'bump',bin_s,rho_thresh,min_turn_rad,min_pairs);
    T(k).bump_frz = window_slope(tr,w_frz(1),w_frz(2),lag_s(di),'bump',bin_s,rho_thresh,min_turn_rad,min_pairs);
    T(k).cue_pre  = window_slope(tr,w_pre(1),w_pre(2),0,'cue',bin_s,rho_thresh,min_turn_rad,min_pairs);
    T(k).cue_frz  = window_slope(tr,w_frz(1),w_frz(2),0,'cue',bin_s,rho_thresh,min_turn_rad,min_pairs);
end
% sign-normalize per dataset so closed-loop pre-freeze slopes are positive
for di = 1:2
    in = strcmp({T.dataset},datasets{di});
    cl = in & ~strcmp({T.cond},'dark');
    sb = sign(median([T(cl).bump_pre],'omitnan'));
    sc = sign(median([T(cl).cue_pre],'omitnan'));
    for k = find(in)
        T(k).bump_pre = sb*T(k).bump_pre; T(k).bump_frz = sb*T(k).bump_frz;
        T(k).cue_pre  = sc*T(k).cue_pre;  T(k).cue_frz  = sc*T(k).cue_frz;
    end
end

fprintf('\nmedian slopes (n trials with both windows valid):\n');
for di = 1:2
    for cc = {'bar','grating','dark'}
        in = strcmp({T.dataset},datasets{di}) & strcmp({T.cond},cc{1});
        ok = in & ~isnan([T.bump_pre]) & ~isnan([T.bump_frz]);
        if ~any(in), continue, end
        fprintf('  %-9s %-8s bump pre %5.2f -> freeze %5.2f (n=%d/%d trials, %d flies, paired signrank p=%.3g) | cue pre %5.2f -> freeze %5.2f\n', datasets{di}, cc{1}, ...
            median([T(ok).bump_pre]), median([T(ok).bump_frz]), sum(ok), sum(in), numel(unique({T(ok).fly})), ...
            signrank([T(ok).bump_pre],[T(ok).bump_frz]), ...
            median([T(in).cue_pre],'omitnan'), median([T(in).cue_frz],'omitnan'));
    end
end

%% plot
cond_list  = {'bar','grating','dark'};
cond_color = [0 .2 1; .2 .65 .9; .5 .5 .5];
ds_title   = {'EPG > GCaMP7f (epg7f)','EPG > syt8m (hackathon, G4 trials)'};

figure('color','w','Position',[50 50 1100 1150]); clf
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');

for di = 1:2
    % row 1: example closed-loop bar trial = most turning during the freeze
    % among trials with valid pre and freeze slopes
    cand = find(strcmp({T.dataset},datasets{di}) & strcmp({T.cond},'bar') & ~isnan([T.bump_pre]) & ~isnan([T.bump_frz]));
    [~,b] = max([T(cand).freeze_turn]); tr = T(cand(b));
    t0 = tr.t_freeze - 30; t1 = tr.xf(end);
    ax = nexttile(tl,di); hold(ax,'on')
    patch(ax,[tr.t_freeze t1 t1 tr.t_freeze],[-1 -1 1 1]*1e3,[.92 .92 .92],'EdgeColor','none')
    kf = tr.xf>=t0 & tr.xf<=t1;
    h_ref = interp1(tr.xf,tr.heading,tr.t_freeze);
    c_ref = interp1(tr.xf,unwrap(tr.cue),tr.t_freeze);
    % bump: unwrap over good frames only, referenced to the freeze onset
    kb = find(tr.xb>=t0 & tr.xb<=t1 & tr.keep & tr.rho>rho_thresh);
    mu_u = unwrap(tr.mu(kb));
    [~,r0] = min(abs(tr.xb(kb)-tr.t_freeze));
    raw_sign = raw_bump_sign(tr,lag_s(di),bin_s,rho_thresh,win_pre_s);
    cue_sign = raw_cue_sign(tr,bin_s,win_pre_s);
    plot(ax,tr.xf(kf),tr.heading(kf)-h_ref,'-','Color',[.6 .6 .6],'LineWidth',1)
    plot(ax,tr.xf(kf),cue_sign*(unwrap(tr.cue(kf))-c_ref),'-k','LineWidth',1.2)
    scatter(ax,tr.xb(kb),raw_sign*(mu_u-mu_u(r0)),6,cond_color(1,:),'filled')
    xline(ax,tr.t_freeze,'--k')
    yl = [min(tr.heading(kf)-h_ref) max(tr.heading(kf)-h_ref)]; yl = yl + [-1 1]*0.1*diff(yl);
    ylim(ax,yl); xlim(ax,[t0 t1])
    xlabel(ax,'time (s)'); ylabel(ax,'displacement from freeze onset (rad)')
    title(ax,sprintf('%s\nexample: %s, %.0f rad turned during freeze',ds_title{di},tr.label,tr.freeze_turn),'Interpreter','none','FontSize',9)
    if di == 1
        legend(ax,{'freeze','fly heading','cue (sign-aligned)','bump (sign-aligned)'},'Location','northwest','FontSize',8)
    end

    % rows 2-3: cue and bump slopes, pre vs freeze, per trial
    for row = 2:3
        ax = nexttile(tl,(row-1)*2+di); hold(ax,'on')
        if row == 2, f_pre = 'cue_pre';  f_frz = 'cue_frz';  yl_str = 'cue displacement / heading change';
        else,        f_pre = 'bump_pre'; f_frz = 'bump_frz'; yl_str = 'bump displacement / heading change'; end
        xt = []; xtl = {};
        for ci = 1:3
            in = find(strcmp({T.dataset},datasets{di}) & strcmp({T.cond},cond_list{ci}));
            if isempty(in), continue, end
            x0 = 3*(ci-1)+1; xt = [xt x0 x0+1]; xtl = [xtl {'pre','freeze'}]; %#ok<AGROW>
            a = [T(in).(f_pre)]; z = [T(in).(f_frz)];
            ok = ~isnan(a) & ~isnan(z); a = a(ok); z = z(ok); % paired trials only
            plot(ax,[x0 x0+1],[a; z],'-','Color',[cond_color(ci,:) .35])
            scatter(ax,x0+zeros(size(a)),a,14,cond_color(ci,:),'filled')
            scatter(ax,x0+1+zeros(size(z)),z,14,cond_color(ci,:),'filled')
            plot(ax,[x0 x0+1],[median(a) median(z)],'-k','LineWidth',2.5)
            text(ax,x0+0.5,1.5,cond_list{ci},'HorizontalAlignment','center','FontWeight','bold','Color',cond_color(ci,:))
            text(ax,x0+0.5,-0.5,sprintf('n=%d trials\n%d flies',numel(a),numel(unique({T(in(ok)).fly}))),'HorizontalAlignment','center','FontSize',7)
            if row == 3
                text(ax,x0+0.5,1.35,sprintf('%.2f \\rightarrow %.2f\np=%.2g',median(a),median(z),signrank(a,z)),'HorizontalAlignment','center','FontSize',7)
            end
        end
        yline(ax,0,':k')
        xticks(ax,xt); xticklabels(ax,xtl); xlim(ax,[0.5 max(xt)+0.5]); ylim(ax,[-0.7 1.6])
        ylabel(ax,yl_str)
        if row == 2, title(ax,'cue vs. heading (confirms the freeze)','FontSize',9)
        else,        title(ax,sprintf('bump vs. heading (bump lag %.1f s)',lag_s(di)),'FontSize',9); end
    end
end
title(tl,{'End-of-trial cue freeze: does the frozen cue hold the EPG bump?', ...
    sprintf('slope of %.1f s bin displacements vs. heading change; pre = %d s before freeze, freeze = onset+%d s to trial end; one line per trial (bold = median)',bin_s,win_pre_s,settle_s)},'FontSize',10)

exportgraphics(gcf,fullfile(export_dir,'fig_gcamp_bar_freeze.png'),'Resolution',200)

%% ===================== local functions =====================
function cond = cond_from_pattern(p)
    p = char(p);
    if contains(p,'background'), cond = 'dark';
    elseif contains(p,'grating'), cond = 'grating';
    else, cond = 'bar'; end
end

function s = window_slope(tr, ta, tb, lag, which, bin_s, rho_thresh, min_turn, min_pairs)
    % regression slope (through the origin) of bin-to-bin displacement of
    % the bump ('bump') or cue ('cue') on bin-to-bin heading change, within
    % [ta,tb]. Heading is read lag s before each bin centre.
    edges = ta:bin_s:tb;
    s = nan;
    if numel(edges) < 3, return, end
    ctr = edges(1:end-1) + bin_s/2;
    y = nan(numel(ctr),1);
    for b = 1:numel(ctr)
        if strcmp(which,'bump')
            k = tr.xb>=edges(b) & tr.xb<edges(b+1) & tr.keep & tr.rho>rho_thresh;
            if any(k), y(b) = angle(sum(exp(1i*tr.mu(k)))); end
        else
            k = tr.xf>=edges(b) & tr.xf<edges(b+1);
            if any(k), y(b) = angle(sum(exp(1i*tr.cue(k)))); end
        end
    end
    h  = interp1(tr.xf,tr.heading,ctr(:)-lag);
    dy = angle(exp(1i*diff(y)));
    dh = diff(h);
    ok = ~isnan(dy) & ~isnan(dh);
    if sum(ok) < min_pairs || sum(abs(dh(ok))) < min_turn, return, end
    s = dh(ok) \ dy(ok);
end

function sg = raw_bump_sign(tr, lag, bin_s, rho_thresh, win_pre_s)
    sg = sign(window_slope(tr,tr.t_freeze-win_pre_s,tr.t_freeze,lag,'bump',bin_s,rho_thresh,0,2));
    if isnan(sg) || sg == 0, sg = 1; end
end

function sg = raw_cue_sign(tr, bin_s, win_pre_s)
    sg = sign(window_slope(tr,tr.t_freeze-win_pre_s,tr.t_freeze,0,'cue',bin_s,0,0,2));
    if isnan(sg) || sg == 0, sg = 1; end
end

function is_flash = detect_flash_frames(f, mad_k)
    % as in Figure_1_claude.m
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;
    is_bright  = frame_mean > median(frame_mean) + mad_k*mad(frame_mean,1);
    is_uniform = frame_cv <= median(frame_cv);
    is_flash   = is_bright & is_uniform;
end

function xb = get_xb(ft, n_im)
    if isfield(ft,'xb') && numel(ft.xb) == n_im
        xb = ft.xb(:);
    else
        xb = linspace(ft.xf(1),ft.xf(end),n_im)';
    end
end

function parts = split_path(p)
    parts = regexp(p,'[\\/]','split');
    parts(cellfun(@isempty,parts)) = [];
end

function cue_clean = despike_cue(cue, max_dev_rad, max_run, win)
    % as in Figure_1_claude.m: replace isolated 1-4 sample spikes in the
    % wrapped cue trace by circular interpolation
    if nargin < 2, max_dev_rad = 0.3; end
    if nargin < 3, max_run = 4; end
    if nargin < 4, win = 7; end
    cue = cue(:);
    cu  = unwrap(cue);
    dev = angle(exp(1i*(cu - movmedian(cu,win,'Endpoints','shrink'))));
    bad = abs(dev) > max_dev_rad;
    d = diff([false; bad; false]); s = find(d==1); e = find(d==-1)-1;
    for r = 1:numel(s)
        if e(r)-s(r)+1 > max_run, bad(s(r):e(r)) = false; end
    end
    cue_clean = cue;
    good = find(~bad);
    if ~any(bad) || numel(good) < 2, return, end
    cue_clean(bad) = angle(interp1(good,exp(1i*cue(good)),find(bad),'linear','extrap'));
end
