%% fig_grab_vs_syt7f_wedge_acf
% Per-wedge version of fig_grab_vs_syt7f_peak_acf.m: does EPG>GRAB(DA2m)
% still look sharper than LPsP>syt7f when each PB wedge is analyzed on its
% own?
%
% Why: the max-across-wedges "peak" trace can jump between wedges as the
% bump moves, so it can change faster than any one wedge does. A single
% wedge's trace has no such jumps, so its ACF is the fairer comparison with
% a sensor's off time constant. (It is still not a pure sensor measurement:
% for EPG>GRAB, a wedge's ACF also includes how long the bump dwells there.)
%
% Same trials and filters as fig_grab_vs_syt7f_peak_acf.m. Per trial, each
% of the 32 wedges' im.d traces gets its own mean-subtracted, biased ACF;
% wedges are then averaged within trial and trials within fly.
%
% Noise: uncorrelated noise shortens an ACF, and single wedges are noisy,
% so a width difference between indicators could be a noise difference.
% Extrapolating the ACF back to lag 0 can't remove this, because im.d in
% lpsp_cl is already temporally smoothed: the lag-1 autocorrelation of its
% frame-to-frame differences is +0.46 to +0.58 (white noise would be -0.5),
% so noise leaks into the first few lags instead of sitting only at lag 0.
% All widths here therefore include that smoothing too.
%
% Instead, the cross-hemisphere correlation pairs each left wedge with its
% homologous right wedge (wedge w with w+16). These see the same signal
% with independent shot noise, so noise drops out of the numerator at every
% lag. It is symmetrized over +/- lag, pooled across the 16 pairs, and
% normalized to its own lag-0 value. It assumes the two hemispheres' signals
% are in phase, which the script checks by printing the lag of the peak.

%% paths + constants
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(repo_root), repo_root = pwd; end % evaluated as a cell from repo root
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','grab_gcamp_tau');

flash_mad_thresh      = 8;
min_forward_speed     = 0.5;
min_forward_time_frac = 0.20;
max_lag_s             = 10;
n_wedge               = 32;

indicator_order = {'LPsP>syt7f','EPG>GRAB(DA2m)'};
ind_color       = {[.7 0 .7],[0 .7 0]};

%% load lpsp_cl
S = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data');
D = S.all_data(:);
clear S

%% per-trial flash-filled wedge traces and peak trace
n = numel(D);
T = struct('indic',repmat({''},n,1),'fly','','ok',false,'d',[],'peak',[],'fs',nan);
for i = 1:n
    meta = D(i).meta;
    if isempty(meta) || ~ischar(meta), continue, end
    ctype = regexp(meta,'_(LPsP|EPG)_','tokens','once');
    if isempty(ctype), continue, end
    T(i).indic = [ctype{1} '>' trial_indicator(meta)];
    T(i).fly   = [ctype{1} ':' trial_fly_id(meta)];
    ft = D(i).ft;
    T(i).ok = ismember(T(i).indic,indicator_order) && ...
              mean(ft.f_speed > min_forward_speed) >= min_forward_time_frac;
    if ~T(i).ok, continue, end

    flash = detect_flash_frames(D(i).im.f,flash_mad_thresh);
    T(i).d    = fill_flash(D(i).im.d,flash);
    T(i).peak = max(T(i).d,[],1);
    T(i).fs   = (size(D(i).im.d,2)-1)/(ft.xf(end)-ft.xf(1)); % lpsp_cl has no ft.xb
end
T = T([T.ok]);
fs = median([T.fs]);
max_lag = round(max_lag_s*fs);
lags_s  = (0:max_lag)/fs;
fprintf('%d trials kept; imaging fs = %.2f Hz\n',numel(T),fs);

%% ACF per wedge -> per fly -> per indicator
R = struct();
for k = 1:numel(indicator_order)
    rows  = find(strcmp({T.indic},indicator_order{k}));
    flies = unique({T(rows).fly});
    nf    = numel(flies);
    R(k).peak     = nan(nf,max_lag+1);
    R(k).wedge    = nan(nf,max_lag+1);
    R(k).cross    = nan(nf,max_lag+1);
    R(k).by_wedge = nan(nf,n_wedge); % width of each wedge position, per fly
    R(k).r0       = nan(nf,1);       % lag-0 left/right correlation
    R(k).pk_lag   = nan(nf,1);       % lag of the unsymmetrized left/right peak, s
    for ff = 1:nf
        tr = rows(strcmp({T(rows).fly},flies{ff}));
        a_w = zeros(n_wedge,max_lag+1);
        a_p = zeros(1,max_lag+1);
        a_c = zeros(1,max_lag+1);
        r0  = zeros(1,numel(tr));
        pk  = zeros(1,numel(tr));
        for j = 1:numel(tr)
            t = tr(j);
            a_w = a_w + row_acf(T(t).d,max_lag);
            a_p = a_p + row_acf(T(t).peak,max_lag);
            [c,r0(j),pk(j)] = cross_hemi(T(t).d,max_lag);
            a_c = a_c + c;
        end
        a_w = a_w/numel(tr);

        R(k).peak(ff,:)  = a_p/numel(tr);
        R(k).wedge(ff,:) = mean(a_w,1,'omitnan');
        R(k).cross(ff,:) = a_c/numel(tr);
        R(k).r0(ff)      = mean(r0);
        R(k).pk_lag(ff)  = mean(pk)/fs;
        for w = 1:n_wedge
            R(k).by_wedge(ff,w) = width_1e(a_w(w,:),lags_s);
        end
    end
    R(k).peak_w  = arrayfun(@(ff) width_1e(R(k).peak(ff,:),lags_s),(1:nf)');
    R(k).wedge_w = arrayfun(@(ff) width_1e(R(k).wedge(ff,:),lags_s),(1:nf)');
    R(k).cross_w = arrayfun(@(ff) width_1e(R(k).cross(ff,:),lags_s),(1:nf)');
    R(k).n_fly = nf; R(k).n_trial = numel(rows);
    fprintf('%-15s %2d flies / %2d trials | 1/e width (s, median [IQR]): peak %.2f [%.2f %.2f], wedge %.2f [%.2f %.2f], cross-hemi %.2f [%.2f %.2f]\n', ...
        indicator_order{k},nf,numel(rows), ...
        median(R(k).peak_w),prctile(R(k).peak_w,[25 75]), ...
        median(R(k).wedge_w),prctile(R(k).wedge_w,[25 75]), ...
        median(R(k).cross_w,'omitnan'),prctile(R(k).cross_w,[25 75]));
    fprintf('   left/right lag-0 r: median %.2f; L/R peak lag: median %.3f s; signrank peak vs wedge: p = %.3g; wedge vs cross-hemi: p = %.3g\n', ...
        median(R(k).r0),median(R(k).pk_lag),signrank(R(k).peak_w,R(k).wedge_w),signrank(R(k).wedge_w,R(k).cross_w));
end
metric_names = {'peak_w','wedge_w','cross_w'};
for s = 1:numel(metric_names)
    fprintf('ranksum %s vs %s, %s: p = %.3g\n',indicator_order{:},metric_names{s}, ...
        ranksum(R(1).(metric_names{s}),R(2).(metric_names{s})));
end

%% figure
fig = figure('Color','w','Position',[50 50 1500 850]);
tl = tiledlayout(fig,2,3,'TileSpacing','compact','Padding','compact');

% mean ACFs: single wedge vs. peak trace
ax = nexttile(tl,1); hold(ax,'on')
h = gobjects(1,2*numel(indicator_order));
for k = 1:numel(indicator_order)
    h(2*k-1) = plot_mean_sem(ax,lags_s,R(k).wedge,ind_color{k},'-');
    h(2*k)   = plot_mean_sem(ax,lags_s,R(k).peak,ind_color{k},'--');
    h(2*k-1).DisplayName = sprintf('%s, single wedge (%d flies)',indicator_order{k},R(k).n_fly);
    h(2*k).DisplayName   = sprintf('%s, peak trace',indicator_order{k});
end
yline(ax,exp(-1),'k:'); yline(ax,0,'k-');
xlim(ax,[0 max_lag_s]); ylim(ax,[-.2 1]);
xlabel(ax,'lag (s)'); ylabel(ax,'autocorrelation'); title(ax,'single wedge vs. peak trace')
legend(ax,h,'Location','northeast','FontSize',8)

% zoom: within-wedge vs. noise-robust cross-hemisphere
ax = nexttile(tl,2); hold(ax,'on')
h3 = gobjects(1,2*numel(indicator_order));
for k = 1:numel(indicator_order)
    h3(2*k-1) = plot_mean_sem(ax,lags_s,R(k).wedge,ind_color{k},'-');
    h3(2*k)   = plot_mean_sem(ax,lags_s,R(k).cross,ind_color{k},':');
    h3(2*k-1).DisplayName = sprintf('%s, within wedge',indicator_order{k});
    h3(2*k).DisplayName   = sprintf('%s, cross-hemisphere',indicator_order{k});
end
yline(ax,exp(-1),'k:');
xlim(ax,[0 4]); ylim(ax,[0 1]);
xlabel(ax,'lag (s)'); ylabel(ax,'correlation'); title(ax,'within wedge vs. cross-hemisphere (noise robust)')
legend(ax,h3,'Location','northeast','FontSize',8)

% widths per fly
ax = nexttile(tl,3); hold(ax,'on')
for s = 1:numel(metric_names)
    for k = 1:numel(indicator_order)
        w = R(k).(metric_names{s});
        x = (s-1)*3 + k;
        swarmchart(ax,x*ones(size(w)),w,18,ind_color{k},'filled','MarkerFaceAlpha',.6,'XJitterWidth',.4);
        plot(ax,x+[-.3 .3],median(w,'omitnan')*[1 1],'k-','LineWidth',2);
    end
end
ax.XTick = [1.5 4.5 7.5]; ax.XTickLabel = {'peak','single wedge','cross-hemisphere'};
xlim(ax,[0 9]); ylim(ax,[0 inf]);
ylabel(ax,'1/e width (s), per fly'); title(ax,'width (narrower = sharper)')

% paired per-fly: peak trace vs. single wedge
ax = nexttile(tl,4); hold(ax,'on')
for k = 1:numel(indicator_order)
    x = (k-1)*3 + [1 2];
    plot(ax,x,[R(k).peak_w R(k).wedge_w]','-o','Color',[ind_color{k} .5],'MarkerSize',4, ...
        'MarkerFaceColor',ind_color{k},'MarkerEdgeColor',ind_color{k});
    plot(ax,x,[median(R(k).peak_w) median(R(k).wedge_w)],'k-','LineWidth',2.5);
end
ax.XTick = [1 2 4 5]; ax.XTickLabel = {'peak','wedge','peak','wedge'};
xlim(ax,[.5 5.5]); ylim(ax,[0 inf]);
ylabel(ax,'ACF 1/e width (s), per fly'); title(ax,'same fly: peak trace vs. single wedge (magenta = syt7f, green = GRAB)')

% width by wedge position
ax = nexttile(tl,5,[1 2]); hold(ax,'on')
h2 = gobjects(1,numel(indicator_order));
for k = 1:numel(indicator_order)
    h2(k) = plot_mean_sem(ax,1:n_wedge,R(k).by_wedge,ind_color{k},'-');
    h2(k).DisplayName = indicator_order{k};
end
xline(ax,16.5,'k:');
xlim(ax,[1 n_wedge]); ylim(ax,[0 inf]);
xlabel(ax,'wedge (dotted line = between hemispheres)'); ylabel(ax,'ACF 1/e width (s)')
title(ax,'width by wedge position (mean \pm SEM across flies)')
legend(ax,h2,'Location','best')

title(tl,'Per-wedge ACF: LPsP>syt7f vs. EPG>GRAB(DA2m), Figure 2 trial set (closed loop + dark)')

if ~exist(export_dir,'dir'), mkdir(export_dir); end
exportgraphics(fig,fullfile(export_dir,'fig_grab_vs_syt7f_wedge_acf.png'),'Resolution',150);

%% functions
function d = fill_flash(d,flash)
    % linearly interpolate each wedge over flash frames, keeping frame spacing uniform
    if any(flash) && sum(~flash) >= 2
        t = 1:size(d,2);
        d(:,flash) = interp1(t(~flash),d(:,~flash)',t(flash),'linear','extrap')';
    end
end

function a = row_acf(x,max_lag)
    % biased ACF of each row, normalized so lag 0 = 1
    x = x - mean(x,2,'omitnan');
    x(isnan(x)) = 0;
    n = size(x,2);
    v = sum(x.^2,2);
    a = nan(size(x,1),max_lag+1);
    for L = 0:max_lag
        a(:,L+1) = sum(x(:,1:n-L).*x(:,1+L:n),2)./v;
    end
end

function [c,r0,pk_lag] = cross_hemi(d,max_lag)
    % correlation between homologous wedges (w, w+16), pooled across pairs,
    % symmetrized over +/- lag and normalized to its lag-0 value. Also returns
    % the lag-0 left/right correlation and the lag (frames, left leading =
    % positive) of the unsymmetrized peak.
    d = d - mean(d,2,'omitnan');
    d(isnan(d)) = 0;
    x = d(1:16,:); y = d(17:32,:);
    n = size(d,2);
    fwd = nan(1,max_lag+1); bwd = nan(1,max_lag+1);
    for L = 0:max_lag
        fwd(L+1) = sum(x(:,1:n-L).*y(:,1+L:n),'all');
        bwd(L+1) = sum(y(:,1:n-L).*x(:,1+L:n),'all');
    end
    c  = (fwd+bwd)/(fwd(1)+bwd(1));
    r0 = fwd(1)/sqrt(sum(x.^2,'all')*sum(y.^2,'all'));
    two_sided = [fliplr(bwd(2:end)) fwd];
    [~,j] = max(two_sided);
    pk_lag = j - (max_lag+1);
end

function w = width_1e(a,lags)
    j = find(a < exp(-1),1);
    if isempty(j) || j == 1, w = nan; return, end
    w = interp1(a([j-1 j]),lags([j-1 j]),exp(-1));
end

function h = plot_mean_sem(ax,x,A,col,style)
    mu = mean(A,1,'omitnan');
    se = std(A,0,1,'omitnan')./sqrt(sum(~isnan(A),1));
    patch(ax,[x fliplr(x)],[mu-se fliplr(mu+se)],col,'FaceAlpha',.15,'EdgeColor','none');
    h = plot(ax,x,mu,style,'Color',col,'LineWidth',1.8);
end

function ind = trial_indicator(meta_path)
    % same rules as Figure_2_claude.m
    if contains(meta_path,'_EPG_syt7f','IgnoreCase',true)
        ind = 'GRAB(DA2m)'; % mislabeled at acquisition
    elseif contains(meta_path,'syt7f','IgnoreCase',true)
        ind = 'syt7f';
    elseif contains(meta_path,'grab','IgnoreCase',true) || contains(meta_path,'DA2m','IgnoreCase',true)
        ind = 'GRAB(DA2m)';
    else
        ind = '';
    end
end

function fid = trial_fly_id(meta_path)
    % Figure_2_claude.m's lpsp_cl rule: <date>_<trailing fly suffix>
    parts = strsplit(meta_path,{'\\','/'});
    parts(cellfun(@isempty,parts)) = [];
    tname    = parts{end};
    date_str = regexp(tname,'^\d{8}','match','once');
    suffix   = regexp(tname,'_(\d+)$','tokens','once');
    if isempty(suffix), fid = date_str; else, fid = [date_str,'_',suffix{1}]; end
end

function is_flash = detect_flash_frames(f, mad_k)
    % ported from Figure_2_claude.m
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;
    is_flash   = frame_mean > median(frame_mean) + mad_k*mad(frame_mean,1) & frame_cv <= median(frame_cv);
end
