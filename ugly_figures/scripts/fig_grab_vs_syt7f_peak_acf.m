%% fig_grab_vs_syt7f_peak_acf
% Is the EPG>GRAB(DA2m) PB signal temporally "sharper" than LPsP>syt7f?
%
% Uses the Figure_2_claude.m trial set (lpsp_cl dataset -- the only one of
% Figure 2's three sources that holds syt7f/GRAB trials -- with the same
% indicator, fly-ID, flash-frame, and walked-forward-enough filters).
% The 5 lpsp_cl trials whose folders say "EPG_syt7f" are really
% EPG>GRAB(DA2m) (mislabeled at acquisition) and are counted as GRAB.
%
% Per trial, on the imaging timebase (~9.1 Hz):
%   peak  = max across the 32 PB wedges of im.d (dF/F) per frame
%   p90z  = 90th percentile across wedges of im.z (Figure 2's response variable)
% flash frames are linearly interpolated over, not dropped, so the frame
% spacing (and hence the ACF lag axis) stays uniform.
%
% ACF: each trace is mean-subtracted and normalized (biased xcorr,
% 'coeff'); trials are averaged within fly, flies within indicator. Width
% summary per fly = lag where the ACF first drops below 1/e (linearly
% interpolated between frames). A narrower ACF = faster fluctuations =
% "sharper" signal.
%
% Control: the same ACF on |r_speed| resampled to the imaging frames. If
% the two indicators' behavior ACFs match but their imaging ACFs differ,
% the difference is in the signal (indicator kinetics / cell type), not in
% how the flies happened to walk.

%% paths + constants
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(repo_root), repo_root = pwd; end % evaluated as a cell from repo root
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');

flash_mad_thresh      = 8;
min_forward_speed     = 0.5;
min_forward_time_frac = 0.20;
max_lag_s             = 10;
example_s             = 60;

indicator_order = {'LPsP>syt7f','EPG>GRAB(DA2m)'};
ind_color       = {[.7 0 .7],[0 .7 0]};

%% load lpsp_cl
S = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data');
D = S.all_data(:);
clear S

%% per-trial peak / p90z / |r_speed| on the imaging timebase
n = numel(D);
T = struct('indic',repmat({''},n,1),'fly','','dark',false,'ok',false, ...
           'peak',[],'p90z',[],'rsp',[],'fs',nan);
for i = 1:n
    meta = D(i).meta;
    if isempty(meta) || ~ischar(meta), continue, end
    ctype = regexp(meta,'_(LPsP|EPG)_','tokens','once');
    if isempty(ctype), continue, end
    T(i).indic = [ctype{1} '>' trial_indicator(meta)];
    T(i).fly   = [ctype{1} ':' trial_fly_id(meta)];
    T(i).dark  = contains(meta,'_dark','IgnoreCase',true);
    ft = D(i).ft;
    T(i).ok = ismember(T(i).indic,indicator_order) && ...
              mean(ft.f_speed > min_forward_speed) >= min_forward_time_frac;
    if ~T(i).ok, continue, end

    n_im = size(D(i).im.d,2);
    xb   = linspace(ft.xf(1),ft.xf(end),n_im)'; % lpsp_cl has no ft.xb (same fallback as Figure 2's get_xb)
    flash = detect_flash_frames(D(i).im.f,flash_mad_thresh)';
    T(i).peak = fill_flash(max(D(i).im.d,[],1)',flash);
    T(i).p90z = fill_flash(prctile(D(i).im.z,90,1)',flash);
    T(i).rsp  = interp1(ft.xf(:),abs(ft.r_speed(:)),xb,'linear','extrap');
    T(i).fs   = 1/median(diff(xb));
end
T = T([T.ok]);
fs = median([T.fs]);
max_lag = round(max_lag_s*fs);
lags_s  = (0:max_lag)/fs;
fprintf('%d trials kept; imaging fs = %.2f Hz\n',numel(T),fs);

%% ACF per trial -> per fly -> per indicator
sig_names = {'peak','p90z','rsp'};
sig_label = {'peak dF/F (max across wedges)','p90 z-score across wedges','|r\_speed| (behavior control)'};
R = struct();
for k = 1:numel(indicator_order)
    rows  = find(strcmp({T.indic},indicator_order{k}));
    flies = unique({T(rows).fly});
    for s = 1:numel(sig_names)
        acf_fly = nan(numel(flies),max_lag+1);
        for ff = 1:numel(flies)
            tr = rows(strcmp({T(rows).fly},flies{ff}));
            a  = cell2mat(arrayfun(@(t) trace_acf(T(t).(sig_names{s}),max_lag),tr,'UniformOutput',false)');
            acf_fly(ff,:) = mean(a,1);
        end
        R(k).(sig_names{s})       = acf_fly;
        R(k).([sig_names{s} '_w']) = arrayfun(@(ff) width_1e(acf_fly(ff,:),lags_s),1:numel(flies))';
    end
    R(k).n_fly = numel(flies); R(k).n_trial = numel(rows);
    fprintf('%-11s %2d flies / %2d trials | 1/e width (s, median [IQR]): peak %.2f [%.2f %.2f], p90z %.2f [%.2f %.2f], |r_speed| %.2f [%.2f %.2f]\n', ...
        indicator_order{k},R(k).n_fly,R(k).n_trial, ...
        median(R(k).peak_w),prctile(R(k).peak_w,[25 75]), ...
        median(R(k).p90z_w),prctile(R(k).p90z_w,[25 75]), ...
        median(R(k).rsp_w),prctile(R(k).rsp_w,[25 75]));
end
for s = 1:numel(sig_names)
    p = ranksum(R(1).([sig_names{s} '_w']),R(2).([sig_names{s} '_w']));
    fprintf('ranksum %s vs %s, %s 1/e width: p = %.3g\n',indicator_order{:},sig_names{s},p);
end

%% figure
fig = figure('Color','w','Position',[50 50 1500 900]);
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');

% row 1-2, cols 1-2: example traces (one trial per indicator, highest-variance 60s window)
for k = 1:numel(indicator_order)
    rows = find(strcmp({T.indic},indicator_order{k}) & ~[T.dark]);
    [~,j] = max(arrayfun(@(t) numel(T(t).peak),rows)); t = rows(j);
    nwin = round(example_s*fs);
    idx  = 1:min(nwin,numel(T(t).peak));
    ax = nexttile(tl,(k-1)*3+1,[1 2]); hold(ax,'on')
    tt = (idx-1)/fs;
    yyaxis(ax,'left');  plot(ax,tt,T(t).peak(idx),'-','Color',ind_color{k},'LineWidth',1); ylabel(ax,'peak dF/F')
    ax.YColor = ind_color{k};
    yyaxis(ax,'right'); plot(ax,tt,T(t).rsp(idx),'-','Color',[.5 .5 .5 .6]); ylabel(ax,'|r\_speed| (rad/s)')
    ax.YColor = [.4 .4 .4];
    title(ax,sprintf('%s example: fly %s, first %d s',indicator_order{k},T(t).fly,example_s),'Interpreter','none')
    xlabel(ax,'time (s)')
end

% col 3, rows 1-2: width summary per fly
ax = nexttile(tl,3,[2 1]); hold(ax,'on')
for s = 1:numel(sig_names)
    for k = 1:numel(indicator_order)
        w = R(k).([sig_names{s} '_w']);
        x = (s-1)*3 + k;
        swarmchart(ax,x*ones(size(w)),w,18,ind_color{k},'filled','MarkerFaceAlpha',.6,'XJitterWidth',.4);
        plot(ax,x+[-.3 .3],median(w)*[1 1],'k-','LineWidth',2);
    end
end
ax.XTick = [1.5 4.5 7.5]; ax.XTickLabel = {'peak dF/F','p90 z','|r\_speed|'};
ylabel(ax,'ACF 1/e width (s), per fly'); title(ax,'ACF width (narrower = sharper)')
xlim(ax,[0 9])

% row 3: mean ACFs
for s = 1:numel(sig_names)
    ax = nexttile(tl,6+s); hold(ax,'on')
    h = gobjects(1,numel(indicator_order));
    for k = 1:numel(indicator_order)
        A  = R(k).(sig_names{s});
        mu = mean(A,1); se = std(A,0,1)/sqrt(size(A,1));
        patch(ax,[lags_s fliplr(lags_s)],[mu-se fliplr(mu+se)],ind_color{k},'FaceAlpha',.2,'EdgeColor','none');
        h(k) = plot(ax,lags_s,mu,'-','Color',ind_color{k},'LineWidth',2, ...
            'DisplayName',sprintf('%s (%d flies)',indicator_order{k},R(k).n_fly));
    end
    yline(ax,exp(-1),'k:'); yline(ax,0,'k-');
    xlabel(ax,'lag (s)'); ylabel(ax,'autocorrelation'); title(ax,sig_label{s})
    xlim(ax,[0 max_lag_s]); ylim(ax,[-.2 1])
    if s == 1
        legend(ax,h,'Location','northeast')
    end
end
title(tl,'PB peak activity: LPsP>syt7f vs. EPG>GRAB(DA2m), Figure 2 trial set (closed loop + dark)')

if ~exist(export_dir,'dir'), mkdir(export_dir); end
exportgraphics(fig,fullfile(export_dir,'fig_grab_vs_syt7f_peak_acf.png'),'Resolution',150);

%% functions
function x = fill_flash(x,flash)
    if any(flash) && sum(~flash) >= 2
        t = (1:numel(x))';
        x(flash) = interp1(t(~flash),x(~flash),t(flash),'linear','extrap');
    end
end

function a = trace_acf(x,max_lag)
    x = x(:) - mean(x,'omitnan');
    x(isnan(x)) = 0;
    n = numel(x);
    a = nan(1,max_lag+1);
    v = sum(x.^2);
    for L = 0:max_lag
        a(L+1) = sum(x(1:n-L).*x(1+L:n))/v; % biased estimator, a(1) = 1
    end
end

function w = width_1e(a,lags)
    j = find(a < exp(-1),1);
    if isempty(j), w = nan; return, end
    w = interp1(a([j-1 j]),lags([j-1 j]),exp(-1));
end

function ind = trial_indicator(meta_path)
    % same rules as Figure_2_claude.m
    if contains(meta_path,'_EPG_syt7f','IgnoreCase',true)
        ind = 'GRAB(DA2m)'; % mislabeled at acquisition, see header
    elseif contains(meta_path,'syt7f','IgnoreCase',true)
        ind = 'syt7f';
    elseif contains(meta_path,'grab','IgnoreCase',true) || contains(meta_path,'DA2m','IgnoreCase',true)
        ind = 'GRAB(DA2m)';
    else
        ind = '';
    end
end

function fid = trial_fly_id(meta_path)
    % Figure_2_claude.m's lpsp_cl rule: <date>_<trailing trial suffix>
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
