%% fig_dlight_freeze_vs_dark_methods
% EPG > dLight: is fluorescence during a cue freeze below the speed-matched
% DARK level? epg_dlight_claude_script.m (data scripts/) said yes;
% fig_freeze_fluor_vs_rotation.m says freeze ~= dark (both below moving
% cue). This script reproduces the old analysis and then swaps its choices
% one at a time for the new figure's, to show which choice decides it.
%
% Row "original" = epg_dlight_claude_script.m sections 'detect freeze
% blocks' + 'fluorescence vs. rotational speed' as written:
%   - .data/epg_dlight_20260415.mat
%   - fluorescence = im.z (z-scored per trial and glomerulus) interpolated
%     onto ft.xf; summed over the 32 glomeruli, and peak = max of the
%     hemisphere-averaged 16-glomerulus profile
%   - freeze = samples where the smoothed fly speed cr > 0.1 rad/s AND the
%     smoothed |cue speed| < 0.01 rad/s for >= 0.5 s ('orig' detector).
%     Frozen-cue samples where the fly is NOT turning stay in "closed loop".
%   - no margin around freeze on/offsets; per-fly speed bins need > 5 s
% Other rows change, in turn:
%   fluor 'dff'     - im.d (dF/F, one F0 per trial and glomerulus) instead of im.z
%   detector 'new'  - fig_freeze_fluor_vs_rotation.m's freeze blocks: the
%                     whole cue-still run (heading spans >= 0.5 rad),
%                     ambiguous long still runs dropped, 1.5 s dropped after
%                     every freeze on/offset
%   dataset 0417    - the newer 153-trial build used by the new figure
%
% Columns: summed and peak fluorescence vs. speed (mean +/- sem across
% flies, as in the old script), and the paired per-fly frozen - dark
% difference in each speed bin. Units differ between z and dF/F rows, so
% the difference column is divided by each variant's closed-loop level.
%
% ==WHAT IT SHOWS (2026-10-30)==
% The old "frozen below dark" result is real in the old analysis (summed z
% 10-27% of the closed-loop level below dark over 0.25-1.5 rad/s, paired
% per fly) and survives the new freeze detector, but in dF/F it shrinks to
% ~0 over 0.25-1 rad/s, leaving only the two fastest bins (1-1.5 rad/s,
% 7-11 flies) below dark. im.z z-scores each trial separately, so freeze
% frames (inside closed-loop trials, which carry the larger moving-cue
% signal) are normalized by a larger mean/std than dark frames (normalized
% within dark trials); that pushes freeze below dark in z units.
% Within-trial comparisons (freeze vs. moving cue) are unaffected.
% In every variant the slowest bin (0-0.25 rad/s) has frozen ABOVE dark,
% which is why averaging the difference over all speed bins hides the
% faster-bin deficit.

%% paths + constants
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');

speed_edges = 0:.25:3;
speed_x     = speed_edges(1:end-1) + diff(speed_edges)/2;
min_bin_s   = 5;    % as in epg_dlight_claude_script.m
pre_smooth  = 90; post_smooth = 90;
min_flies_diff = 5; % flies needed to draw a speed bin of the paired difference curve

variants = struct( ...
    'name',     {'original (0415, z, orig detector)','0415, dF/F, orig detector','0415, z, new detector','0415, dF/F, new detector','0417, dF/F, new detector (new figure)'}, ...
    'file',     {'epg_dlight_20260415.mat','epg_dlight_20260415.mat','epg_dlight_20260415.mat','epg_dlight_20260415.mat','epg_dlight_20260417.mat'}, ...
    'fluor',    {'z','dff','z','dff','dff'}, ...
    'detector', {'orig','orig','new','new','new'});
conds    = {'closed loop','frozen','dark'};
cond_rgb = [0 0 0; .85 .1 .1; .2 .4 1];

%% run every variant
loaded = containers.Map();
V = struct();
for v = 1:numel(variants)
    if ~isKey(loaded,variants(v).file)
        S = load(fullfile(data_dir,variants(v).file),'all_data');
        loaded(variants(v).file) = S.all_data(:);
    end
    A = loaded(variants(v).file);
    [V(v).sum, V(v).peak, V(v).flies] = run_variant(A, variants(v).fluor, variants(v).detector, speed_edges, min_bin_s, pre_smooth, post_smooth);
    for m = {'sum','peak'}
        B = V(v).(m{1});
        V(v).([m{1} '_fd']) = matched_diff(B(:,2,:), B(:,3,:)); % freeze - dark
        V(v).([m{1} '_fm']) = matched_diff(B(:,2,:), B(:,1,:)); % freeze - moving
    end
    fprintf('%-42s', variants(v).name);
    for m = {'sum','peak'}
        fd = V(v).([m{1} '_fd']); fm = V(v).([m{1} '_fm']);
        fprintf(' | %s: frz-dark %+.3f (n=%d p=%.2g), frz-moving %+.3f (n=%d p=%.2g)', m{1}, ...
            median(fd,'omitnan'), sum(~isnan(fd)), signrank(fd(~isnan(fd))), median(fm,'omitnan'), sum(~isnan(fm)), signrank(fm(~isnan(fm))));
    end
    fprintf('\n');
end

%% plot: one row per variant
nv = numel(variants);
figure('color','w','Position',[50 50 1250 330*nv]); clf
tl = tiledlayout(nv,3,'TileSpacing','compact','Padding','compact');
for v = 1:nv
    unit = 'z'; if strcmp(variants(v).fluor,'dff'), unit = 'dF/F'; end
    for mi = 1:2
        m = {'sum','peak'}; m = m{mi};
        ax = nexttile(tl,(v-1)*3+mi); hold(ax,'on')
        B = V(v).(m); h = gobjects(1,3);
        for c = 1:3
            y = squeeze(B(:,c,:));
            mu = mean(y,2,'omitnan'); se = std(y,0,2,'omitnan')./sqrt(sum(~isnan(y),2));
            ok = ~isnan(mu) & ~isnan(se);
            patch(ax,[speed_x(ok) fliplr(speed_x(ok))],[mu(ok)'+se(ok)' fliplr(mu(ok)'-se(ok)')],cond_rgb(c,:),'FaceAlpha',.2,'EdgeColor','none')
            h(c) = plot(ax,speed_x,mu,'-o','Color',cond_rgb(c,:),'MarkerFaceColor',cond_rgb(c,:),'MarkerSize',3,'LineWidth',1.5);
        end
        xlim(ax,[0 3]); xlabel(ax,'rotational speed (rad/s)')
        if strcmp(m,'sum'), ylabel(ax,sprintf('summed %s (32 glomeruli)',unit)), else, ylabel(ax,sprintf('peak %s (hemisphere-avg max)',unit)), end
        if mi == 1
            title(ax,sprintf('%s  (n=%d flies)',variants(v).name,numel(V(v).flies)),'Interpreter','none','FontSize',9,'HorizontalAlignment','left','Units','normalized','Position',[0 1.02 0])
        end
        if v == 1 && mi == 1, legend(ax,h,conds,'Location','northwest'), end
    end
    % paired per-fly frozen - dark difference in each speed bin (as in the
    % old script's 'frozen - dark' consistency figures), normalized per
    % metric by that variant's median closed-loop level so sum and peak can
    % share one axis; bins with < min_flies_diff flies are not drawn.
    % Stars: signrank across flies in that bin, p < 0.05.
    ax = nexttile(tl,(v-1)*3+3); hold(ax,'on')
    mcol = [.85 .1 .1; .45 .1 .6];
    hh = gobjects(1,2);
    for mi = 1:2
        m = {'sum','peak'}; m = m{mi};
        B = V(v).(m);
        scale = abs(median(squeeze(mean(B(:,1,:),1,'omitnan')),'omitnan'));
        d = squeeze(B(:,2,:) - B(:,3,:)) / scale;
        n_ok = sum(~isnan(d),2);
        mu = mean(d,2,'omitnan'); se = std(d,0,2,'omitnan')./sqrt(n_ok);
        mu(n_ok < min_flies_diff) = nan; se(n_ok < min_flies_diff) = nan;
        ok = ~isnan(mu);
        patch(ax,[speed_x(ok) fliplr(speed_x(ok))],[mu(ok)'+se(ok)' fliplr(mu(ok)'-se(ok)')],mcol(mi,:),'FaceAlpha',.2,'EdgeColor','none')
        hh(mi) = plot(ax,speed_x,mu,'-o','Color',mcol(mi,:),'MarkerFaceColor',mcol(mi,:),'MarkerSize',3,'LineWidth',1.5);
        for b = find(ok)'
            if signrank(d(b,~isnan(d(b,:)))) < 0.05
                text(ax,speed_x(b),mu(b)+se(b)+0.02,'*','Color',mcol(mi,:),'HorizontalAlignment','center','FontSize',12)
            end
        end
    end
    yline(ax,0,':k'); xlim(ax,[0 1.75])
    xlabel(ax,'rotational speed (rad/s)')
    ylabel(ax,sprintf('frozen - dark, paired per fly\n(/ median closed-loop level)'))
    if v == 1, legend(ax,hh,{'summed','peak'},'Location','southwest'), end
end
title(tl,{'EPG > dLight: freeze vs. dark, old analysis (top row) and one change at a time', ...
    sprintf('right column: per-fly frozen - dark in each speed bin, mean \\pm sem across flies (>= %d flies), * signrank p < 0.05',min_flies_diff)},'FontSize',10)
exportgraphics(gcf,fullfile(export_dir,'fig_dlight_freeze_vs_dark_methods.png'),'Resolution',180)

%% ===================== local functions =====================
function [Bsum, Bpeak, flies] = run_variant(A, fluor, detector, edges, min_bin_s, pre_smooth, post_smooth)
    n = numel(A);
    is_dark = false(n,1); fly = cell(n,1);
    for i = 1:n
        is_dark(i) = contains(A(i).ft.pattern,'background');
        tmp = strsplit(A(i).meta,'\');
        ff = tmp(cellfun(@(x) startsWith(x,'fly '),tmp));
        if isempty(ff), fly{i} = A(i).meta; else, fly{i} = strjoin(tmp(1:find(strcmp(tmp,ff{1}))),'\'); end
    end
    [flies,~,fid] = unique(fly);
    lab = cell(n,1); cr_all = cell(n,1); fs_ = cell(n,1); fp_ = cell(n,1);
    keep_trial = true(n,1);
    for i = 1:n
        ft = A(i).ft; xf = ft.xf(:); dt = mean(diff(xf));
        r  = ft.r_speed(:);
        cr = smoothdata(abs(smoothdata(r,'gaussian',pre_smooth)),'movmean',[post_smooth,0]);
        cr_all{i} = cr;
        L = zeros(numel(xf),1);                  % 0 closed loop, 1 frozen, 2 dark, nan dropped
        if is_dark(i)
            L(:) = 2;
        elseif strcmp(detector,'orig')
            c  = ft.cue(:);
            dc = smoothdata([diff(unwrap(c))/dt;0],'gaussian',pre_smooth);
            cc = smoothdata(abs(dc),'movmean',[post_smooth,0]);
            frz = cr > .1 & bwareaopen(cc < 1e-2, max(1,round(.5/dt)));
            L(frz) = 1;
        else
            [frz, amb] = find_freezes(ft.cue(:), r, 1/dt, 1e-3, 6, 3, 0.5, 8);
            if mean(frz) > 0.9, keep_trial(i) = false; end
            L(frz) = 1; L(amb) = nan;
            nset = round(1.5/dt);
            for e = find(diff([0; frz; 0]) ~= 0)'
                L(e:min(numel(L),e+nset-1)) = nan;
            end
        end
        lab{i} = L;
        if strcmp(fluor,'z'), X = A(i).im.z; else, X = A(i).im.d; end
        Xf = interp1(A(i).ft.xb, X', xf);         % xf x 32
        half = size(Xf,2)/2;
        fs_{i} = sum(Xf,2);
        fp_{i} = max((Xf(:,1:half) + Xf(:,half+1:end))/2,[],2);
    end
    nominal_dt = mean(diff(A(1).ft.xf));
    nb = numel(edges)-1; nf = numel(flies);
    Bsum = nan(nb,3,nf); Bpeak = nan(nb,3,nf);
    for f = 1:nf
        tr = find(fid==f & keep_trial);
        if isempty(tr), continue, end
        cr = cat(1,cr_all{tr}); L = cat(1,lab{tr}); ys = cat(1,fs_{tr}); yp = cat(1,fp_{tr});
        for c = 1:3
            k = L == c-1;
            Bsum(:,c,f)  = bin_by_speed(cr(k), ys(k), edges, nominal_dt, min_bin_s);
            Bpeak(:,c,f) = bin_by_speed(cr(k), yp(k), edges, nominal_dt, min_bin_s);
        end
    end
end

function d = matched_diff(a, b)
    % per fly: mean over speed bins where both curves exist (>= 2 bins)
    a = squeeze(a); b = squeeze(b);
    if isvector(a), a = a(:); b = b(:); end
    ok = ~isnan(a) & ~isnan(b);
    d = nan(size(a,2),1);
    for f = 1:size(a,2)
        if sum(ok(:,f)) >= 2, d(f) = mean(a(ok(:,f),f) - b(ok(:,f),f)); end
    end
end

function m = bin_by_speed(cr, y, edges, dt, min_sec)
    % as in epg_dlight_claude_script.m (strictly more than min_sec)
    m = nan(numel(edges)-1,1);
    for j = 1:numel(edges)-1
        idx = cr >= edges(j) & cr < edges(j+1);
        if sum(idx)*dt > min_sec, m(j) = mean(y(idx),'omitnan'); end
    end
end

function [frz, amb] = find_freezes(cue, r, fs, still_step, bridge_n, min_s, min_range, amb_s)
    % as in fig_freeze_fluor_vs_rotation.m
    n = numel(cue);
    still = [abs(angle(exp(1i*diff(cue)))) < still_step; true];
    d = diff([1; still; 1]); s = find(d==-1); e = find(d==1)-1;
    for k = 1:numel(s)
        if e(k)-s(k)+1 <= bridge_n && s(k) > 1 && e(k) < n && abs(angle(exp(1i*(cue(e(k)+1)-cue(s(k)))))) < still_step
            still(s(k):e(k)) = true;
        end
    end
    heading = cumsum(r)/fs;
    frz = false(n,1); amb = false(n,1);
    d = diff([0; still; 0]); s = find(d==1); e = find(d==-1)-1;
    for k = 1:numel(s)
        dur = (e(k)-s(k)+1)/fs;
        if dur < min_s, continue, end
        h = heading(s(k):e(k));
        if max(h)-min(h) >= min_range
            frz(s(k):e(k)) = true;
        elseif dur >= amb_s
            amb(s(k):e(k)) = true;
        end
    end
end
