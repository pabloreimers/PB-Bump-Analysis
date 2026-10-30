%% fig_freeze_fluor_vs_rotation
% Mean and peak PB fluorescence vs. the fly's rotational speed, for three
% rotation-scaling indicators, split by visual condition:
%   moving cue - closed loop, cue following the fly
%   freeze     - closed loop, cue held still while the fly keeps turning
%   dark       - background pattern, no visual cue
%
%   LPsP > syt7f      .data/lpsp_cl_data_20240206.mat ('LPsP_syt7f' trials)
%   EPG > GRAB(DA2m)  .data/lpsp_cl_data_20240206.mat ('GRAB' trials + the
%                     5 mislabeled 'EPG_syt7f' trials, see trial_indicator)
%   EPG > dLight      .data/epg_dlight_20260417.mat (newest build, 153 trials)
%
% ==CONDITIONS==
% lpsp_cl has no ft.pattern; dark vs. closed loop comes from the trial's
% trialSettings.csv on Z:\pablo\lpsp_cl (patternPath contains 'background'),
% falling back to the folder name ('_dark') when no CSV is found -- same
% resolution as grab_freeze_scratch.m. epg_dlight stores ft.pattern.
%
% Freezes are found from the cue itself: runs where the cue does not move
% (|dcue| < still_step_rad per 60 Hz sample, bridging cue-log glitches of up
% to bridge_n samples) lasting >= min_freeze_s, during which the fly's
% integrated heading spans >= min_turn_range_rad. In closed loop a heading
% excursion that large would have moved the bar many pixels, so such a run
% can only be a freeze. lpsp_cl trials have one ~20 s freeze at the end;
% many dLight trials also have 10 s freezes every 60 s. Cue-still runs that
% do NOT qualify (the fly barely turned) but last >= ambiguous_s could be a
% freeze during which the fly stood still, so they are dropped from all
% conditions; shorter still runs are ordinary closed loop with a standing
% fly. Trials whose cue is "frozen" for > 90% of the trial (cue not logged /
% not closed loop; dLight 20260417-2) are dropped.
%
% ==FLUORESCENCE METRICS== (per imaging frame, flash frames excluded)
%   mean - mean dF/F (im.d) over all 32 glomeruli
%   peak - hemisphere-averaged 16-glomerulus ring, best 3 circularly
%          adjacent glomeruli averaged, max over the ring (best3_neighbor_peak
%          from grab_freeze_scratch.m; less sensitive to one noisy glomerulus
%          than a raw max)
% dF/F rather than im.z: z-scoring per trial would rescale dark trials
% separately from closed-loop trials and hide a condition difference.
%
% ==ROTATIONAL SPEED== |r_speed| smoothed as in grab_freeze_scratch.m
% (gaussian, then a causal ~1.5 s running mean), so each frame is matched to
% the fly's RECENT rotation -- these indicators rise and decay over ~1 s
% after a turn. Frames within settle_s of a freeze onset or offset are
% dropped from every condition so the slow signal from one condition doesn't
% leak into the other.
%
% ==TIME-IN-TRIAL CAVEAT==
% lpsp_cl freezes all sit in the last ~20 s of a 600 s trial, where dF/F has
% drifted (bleaching / baseline), so freeze vs. moving cue is partly
% confounded with time in trial there. Checked 2026-10-30 against
% moving-cue frames from only the 60 s before each freeze: for syt7f about
% half of the freeze-vs-moving level drop is this drift (0.87 instead of
% ~0.8, n.s.); for dLight, whose freezes recur through the trial, it changes
% nothing. That reference was shown in an earlier version and then dropped
% from the figure at the user's request.
%
% Per fly: bin frames by speed, mean per bin (needs >= min_bin_s of data);
% curves are mean +/- sem across flies. Slope rows: each fly's least-squares
% slope of fluorescence on speed per condition (frames pooled, >= min_fit_s
% of data spanning >= min_speed_span rad/s); freeze is tested against moving
% cue and against dark (paired by fly). Bottom row: each fly's speed-matched
% level ratio, freeze / moving cue and freeze / dark (mean over the speed
% bins both have) -- catches an offset that the slope alone would miss.

%% paths + constants
repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root,'.data');
export_dir = fullfile(repo_root,'ugly_figures','exports');
addpath(fullfile(repo_root,'circ_stats'))

still_step_rad     = 1e-3; % rad per sample, cue counted as not moving
bridge_n           = 6;    % samples, bridge moving-cue glitches this short inside a still run
min_freeze_s       = 3;    % s
min_turn_range_rad = 0.5;  % rad, heading excursion during a still run that proves a freeze
ambiguous_s        = 8;    % s, non-qualifying still runs this long are dropped
settle_s           = 1.5;  % s, dropped after each freeze onset/offset
pre_smooth         = 90;   % samples (60 Hz), gaussian on r_speed
post_smooth        = 90;   % samples, causal running mean of |r_speed|
flash_mad_thresh   = 8;
speed_edges        = 0:0.25:3;
min_bin_s          = 2;    % s of a fly's data needed in a speed bin
min_flies_per_bin  = 3;    % draw the across-fly mean only where >= this many flies
min_fit_s          = 10;   % s of data needed for a per-fly slope
min_speed_span     = 1;    % rad/s, speed range needed for a per-fly slope

sensors     = {'syt7f','GRAB(DA2m)','dLight'};
sensor_lbl  = {'LPsP > syt7f','EPG > GRAB(DA2m)','EPG > dLight'};
% condition codes in T.cond: 0 moving cue, 1 freeze, 2 dark
conds       = {'moving cue','freeze','dark'};
cond_rgb    = [0 0 0; .85 .1 .1; .2 .4 1];

%% load lpsp_cl + epg_dlight
S = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); A1 = S.all_data(:);
S = load(fullfile(data_dir,'epg_dlight_20260417.mat'),'all_data');   A2 = S.all_data(:);
csv_idx = build_csv_index('Z:\pablo\lpsp_cl\');

T = struct('sensor',{},'fly',{},'dark',{},'label',{},'xb',{},'spd',{},'cond',{},'fmean',{},'fpeak',{},'n_freeze',{},'freeze_s',{});
n_drop_long = 0; dropped = {};
for src = 1:2
    if src == 1, A = A1; ds = 'lpsp_cl'; else, A = A2; ds = 'epg_dlight'; end
    for i = 1:numel(A)
        meta = A(i).meta;
        if isempty(meta) || ~ischar(meta) || isempty(A(i).ft) || isempty(A(i).ft.cue), continue, end
        ind = trial_indicator(meta);
        si  = find(strcmp(sensors,ind));
        if isempty(si), continue, end
        ft = A(i).ft; im = A(i).im;

        % dark vs. closed loop
        if isfield(ft,'pattern') && ~isempty(ft.pattern)
            is_dark = contains(char(ft.pattern),'background','IgnoreCase',true);
        else
            res = resolve_csv_for_trial(meta, trial_folder_name(meta), csv_idx);
            if any(strcmp(res.status,{'exact','fallback'})) && ~isnan(res.info.is_dark)
                is_dark = res.info.is_dark;
            else
                is_dark = contains(meta,'_dark','IgnoreCase',true);
            end
        end

        % behavior timebase labels: 0 moving cue, 1 freeze, 2 dark, nan dropped
        xf = ft.xf(:); fs = 1/median(diff(xf)); r = ft.r_speed(:);
        lab = zeros(numel(xf),1); n_frz = 0; frz_s = 0;
        if is_dark
            lab(:) = 2;
        else
            [frz, amb] = find_freezes(ft.cue(:), r, fs, still_step_rad, bridge_n, min_freeze_s, min_turn_range_rad, ambiguous_s);
            if mean(frz) > 0.9
                dropped{end+1} = sprintf('%s (cue still for %.0f%% of the trial)', trial_folder_name(meta), 100*mean(frz)); %#ok<SAGROW>
                n_drop_long = n_drop_long + 1;
                continue
            end
            lab(frz) = 1;
            lab(amb) = nan;
            % settle margins after every freeze onset and offset
            edges_i = find(diff([0; frz; 0]) ~= 0);
            n_set = round(settle_s*fs);
            for e = edges_i(:)'
                lab(e:min(numel(lab),e+n_set-1)) = nan;
            end
            d = diff([0; frz; 0]); n_frz = sum(d==1); frz_s = sum(frz)/fs;
        end
        cr = smoothdata(abs(smoothdata(r,'gaussian',pre_smooth)),'movmean',[post_smooth 0]);

        % fluorescence at imaging frames
        n_im = size(im.d,2);
        if isfield(ft,'xb') && numel(ft.xb) == n_im, xb = ft.xb(:); else, xb = linspace(xf(1),xf(end),n_im)'; end
        keep = ~reshape(detect_flash_frames(im.f,flash_mad_thresh),[],1);
        fmean = mean(im.d,1)';
        fpeak = best3_neighbor_peak(im.d)';
        spd   = interp1(xf,cr,xb);
        cnd   = interp1(xf,lab,xb,'nearest');
        cnd(~keep) = nan;

        tr.sensor = sensors{si};
        tr.fly    = [ds '|' trial_fly_id(ds,meta)];
        tr.dark   = is_dark;
        tr.label  = trial_folder_name(meta);
        tr.xb = xb; tr.spd = spd; tr.cond = cnd; tr.fmean = fmean; tr.fpeak = fpeak;
        tr.n_freeze = n_frz; tr.freeze_s = frz_s;
        T(end+1) = tr; %#ok<SAGROW>
    end
end
fprintf('dropped %d trials: %s\n', n_drop_long, strjoin(dropped,'; '));
for si = 1:3
    in = strcmp({T.sensor},sensors{si});
    fprintf('%-11s %3d trials (%d CL, %d dark), %d flies; %d freezes, %.0f s of freeze in %d trials\n', sensors{si}, sum(in), ...
        sum(in & ~[T.dark]), sum(in & [T.dark]), numel(unique({T(in).fly})), sum([T(in).n_freeze]), sum([T(in).freeze_s]), sum(in & [T.n_freeze]>0));
end

%% per fly: speed-binned curves + slopes, per sensor x condition
speed_x = speed_edges(1:end-1) + diff(speed_edges)/2;
% the comparisons that are tested: [condition, its reference]
tests = [2 1; 2 3]; % freeze vs. moving cue; freeze vs. dark
R = struct();
for si = 1:3
    in = find(strcmp({T.sensor},sensors{si}));
    flies = unique({T(in).fly});
    nf = numel(flies);
    R(si).flies = flies;
    for m = {'mean','peak'}
        R(si).(['curve_' m{1}]) = nan(numel(speed_x),3,nf);
        R(si).(['slope_' m{1}]) = nan(nf,3);
        R(si).(['ratio_' m{1}]) = nan(nf,2); % freeze / moving cue, freeze / dark
    end
    for f = 1:nf
        tf = in(strcmp({T(in).fly},flies{f}));
        spd = cat(1,T(tf).spd); cnd = cat(1,T(tf).cond);
        F.mean = cat(1,T(tf).fmean); F.peak = cat(1,T(tf).fpeak);
        dt  = median(diff(T(tf(1)).xb));
        for c = 1:3
            k = cnd == c-1 & ~isnan(spd);
            for m = {'mean','peak'}
                y = F.(m{1});
                R(si).(['curve_' m{1}])(:,c,f) = bin_by_speed(spd(k),y(k),speed_edges,dt,min_bin_s);
                if sum(k)*dt >= min_fit_s && range(spd(k)) >= min_speed_span
                    b = [ones(sum(k),1) spd(k)] \ y(k);
                    R(si).(['slope_' m{1}])(f,c) = b(2);
                end
            end
        end
        for m = {'mean','peak'}
            Cv = R(si).(['curve_' m{1}])(:,:,f);
            for rc = 1:2
                ref = tests(rc,2);
                ok = ~isnan(Cv(:,2)) & ~isnan(Cv(:,ref));
                if sum(ok) >= 2, R(si).(['ratio_' m{1}])(f,rc) = mean(Cv(ok,2))/mean(Cv(ok,ref)); end
            end
        end
    end
end

fprintf('\nper-fly slope of fluorescence on rotational speed (dF/F per rad/s), median [n flies]\n');
for si = 1:3
    for m = {'mean','peak'}
        Sl = R(si).(['slope_' m{1}]);
        s = sprintf('%-11s %-4s', sensors{si}, m{1});
        for c = 1:3
            s = [s sprintf('  %s %.3f [%d]', conds{c}, median(Sl(:,c),'omitnan'), sum(~isnan(Sl(:,c))))]; %#ok<AGROW>
        end
        for t = 1:2
            ok = ~isnan(Sl(:,tests(t,1))) & ~isnan(Sl(:,tests(t,2)));
            s = [s sprintf(' | %s vs %s p=%.2g', conds{tests(t,1)}, conds{tests(t,2)}, signrank(Sl(ok,tests(t,1)),Sl(ok,tests(t,2))))]; %#ok<AGROW>
        end
        for rc = 1:2
            rt = R(si).(['ratio_' m{1}])(:,rc);
            s = [s sprintf(' | speed-matched freeze/%s level %.2f (n=%d, p=%.2g)', conds{tests(rc,2)}, median(rt,'omitnan'), sum(~isnan(rt)), signrank(rt(~isnan(rt))-1))]; %#ok<AGROW>
        end
        fprintf('%s\n', s);
    end
end

%% plot: rows = mean curve, peak curve, mean slope, peak slope, slope delta, level ratio; cols = sensors
figure('color','w','Position',[50 50 1300 1950]); clf
tl = tiledlayout(6,3,'TileSpacing','compact','Padding','compact');
metric = {'mean','peak'};
ylab   = {'mean dF/F (all glomeruli)','peak dF/F (best 3 adjacent)'};
for row = 1:2
    for si = 1:3
        ax = nexttile(tl,(row-1)*3+si); hold(ax,'on')
        Y = R(si).(['curve_' metric{row}]);
        h = gobjects(1,3);
        for c = 1:3
            y = squeeze(Y(:,c,:)); if isvector(y), y = y(:); end
            n_ok = sum(~isnan(y),2);
            mu = mean(y,2,'omitnan'); se = std(y,0,2,'omitnan')./sqrt(n_ok);
            mu(n_ok < min_flies_per_bin) = nan; se(n_ok < min_flies_per_bin) = nan;
            plot(ax,speed_x,y,'-','Color',[cond_rgb(c,:) .12])
            v = ~isnan(mu);
            patch(ax,[speed_x(v) fliplr(speed_x(v))],[mu(v)'+se(v)' fliplr(mu(v)'-se(v)')],cond_rgb(c,:),'FaceAlpha',.2,'EdgeColor','none')
            h(c) = plot(ax,speed_x,mu,'-o','Color',cond_rgb(c,:),'MarkerFaceColor',cond_rgb(c,:),'MarkerSize',4,'LineWidth',2);
        end
        n_fly_c = arrayfun(@(c) sum(any(~isnan(squeeze(Y(:,c,:))),1)), 1:3);
        xlabel(ax,'rotational speed (rad/s)'); ylabel(ax,ylab{row})
        xlim(ax,[0 speed_edges(end)])
        if row == 1
            title(ax,sprintf('%s\nflies: moving %d, freeze %d, dark %d',sensor_lbl{si},n_fly_c),'FontSize',10)
        end
        if row == 1 && si == 1, legend(ax,h,conds,'Location','northwest'), end
    end
end
order = [1 2 3]; % x positions: moving, freeze, dark
for row = 3:4
    m = metric{row-2};
    for si = 1:3
        ax = nexttile(tl,(row-1)*3+si); hold(ax,'on')
        Sl = R(si).(['slope_' m])(:,order);
        plot(ax,1:3,Sl','-','Color',[0 0 0 .15])
        for x = 1:3
            c = order(x);
            scatter(ax,x*ones(size(Sl,1),1),Sl(:,x),16,cond_rgb(c,:),'filled')
            med = median(Sl(:,x),'omitnan');
            plot(ax,x+[-.3 .3],[med med],'-','Color',cond_rgb(c,:),'LineWidth',3)
        end
        yl = ylim(ax);
        for t = 1:2
            xa = find(order==tests(t,1)); xr = find(order==tests(t,2));
            ok = ~isnan(Sl(:,xa)) & ~isnan(Sl(:,xr));
            if sum(ok) >= 3
                % label sits over the reference column (both tests share the freeze column)
                text(ax,xr,yl(2),sprintf('freeze vs %s\np=%.2g (n=%d)',conds{tests(t,2)},signrank(Sl(ok,xa),Sl(ok,xr)),sum(ok)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7)
            end
        end
        yline(ax,0,':k')
        xticks(ax,1:3); xticklabels(ax,conds(order)); xlim(ax,[.5 3.5])
        ylabel(ax,sprintf('slope of %s dF/F\non speed (per rad/s)',m))
        if si == 1, title(ax,sprintf('per-fly %s-fluorescence slope (one line per fly)',m),'FontSize',9), end
    end
end
% row 5: per-fly slope difference, freeze minus reference (same 4 groups as
% the level-ratio row below); the signrank against 0 is the same paired
% test as in rows 3-4
for si = 1:3
    ax = nexttile(tl,12+si); hold(ax,'on')
    xl = {}; all_d = [];
    for mi = 1:2
        Sl = R(si).(['slope_' metric{mi}]);
        for rc = 1:2
            x = 2*(mi-1) + rc + 0.4*(mi-1);
            ref = tests(rc,2);
            dlt = Sl(:,2) - Sl(:,ref);
            scatter(ax,x+0.07*randn(size(dlt)),dlt,16,cond_rgb(ref,:),'filled','MarkerFaceAlpha',.6)
            plot(ax,x+[-.3 .3],median(dlt,'omitnan')*[1 1],'-k','LineWidth',3)
            all_d = [all_d; dlt]; %#ok<AGROW>
            xl{end+1} = sprintf('vs %s',conds{ref}); %#ok<SAGROW>
        end
    end
    % y-limits from the data, with headroom for the stats text
    lim = max(abs(all_d),[],'omitnan'); if isempty(lim) || isnan(lim) || lim == 0, lim = 1; end
    ylim(ax,[-1.1 1.6]*lim)
    for mi = 1:2
        Sl = R(si).(['slope_' metric{mi}]);
        for rc = 1:2
            x = 2*(mi-1) + rc + 0.4*(mi-1);
            dlt = Sl(:,2) - Sl(:,tests(rc,2)); ok = ~isnan(dlt);
            if sum(ok) >= 3
                text(ax,x,1.55*lim,sprintf('%.3f\np=%.2g\n(n=%d)',median(dlt(ok)),signrank(dlt(ok)),sum(ok)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7)
            end
        end
    end
    yline(ax,0,':k')
    xticks(ax,[1 2 3.4 4.4]); xticklabels(ax,xl); xlim(ax,[.4 5])
    ylabel(ax,sprintf('\\Delta slope, freeze - reference\n(dF/F per rad/s)'))
    title(ax,'mean                                   peak','FontSize',9)
end

% row 6: speed-matched level ratio
for si = 1:3
    ax = nexttile(tl,15+si); hold(ax,'on')
    % 4 groups: mean vs moving cue, mean vs dark, peak vs moving cue, peak vs dark
    xl = {};
    for mi = 1:2
        for rc = 1:2
            x = 2*(mi-1) + rc + 0.4*(mi-1); % small gap between the mean and peak pairs
            ref = tests(rc,2);
            rt = R(si).(['ratio_' metric{mi}])(:,rc);
            scatter(ax,x+0.07*randn(size(rt)),rt,16,cond_rgb(ref,:),'filled','MarkerFaceAlpha',.6)
            plot(ax,x+[-.3 .3],median(rt,'omitnan')*[1 1],'-k','LineWidth',3)
            ok = ~isnan(rt);
            if sum(ok) >= 3
                text(ax,x,1.95,sprintf('%.2f\np=%.2g\n(n=%d)',median(rt(ok)),signrank(rt(ok)-1),sum(ok)),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',7)
            end
            xl{end+1} = sprintf('vs %s',conds{ref}); %#ok<SAGROW>
        end
        text(ax,1.5 + 2.4*(mi-1),2.15,metric{mi},'HorizontalAlignment','center','FontWeight','bold','FontSize',9)
    end
    yline(ax,1,':k')
    xticks(ax,[1 2 3.4 4.4]); xticklabels(ax,xl); xlim(ax,[.4 5]); ylim(ax,[0.3 2.0])
    ylabel(ax,sprintf('freeze / reference\nspeed-matched level'))
end
title(tl,{'PB fluorescence vs. rotational speed: moving cue vs. cue freeze vs. dark', ...
    sprintf('speed = causal ~1.5 s smoothed |r\\_speed|; per-fly bins (>= %d s), mean \\pm sem across flies (>= %d flies); %.1f s after freeze on/offsets excluded',min_bin_s,min_flies_per_bin,settle_s)},'FontSize',10)
exportgraphics(gcf,fullfile(export_dir,'fig_freeze_fluor_vs_rotation.png'),'Resolution',200)

%% companion figure: speed dependence of the freeze effect
% The level-ratio row above averages over all shared speed bins, which can
% hide a speed-dependent effect: for dLight, freeze sits ABOVE dark in the
% slowest bin and BELOW it at 1-1.5 rad/s (see
% fig_dlight_freeze_vs_dark_methods.m). Here the paired per-fly difference
% freeze - reference is shown in every speed bin instead (dF/F units, same
% per-fly curves as rows 1-2), mean +/- sem across flies, drawn where >=
% min_flies_per_bin flies have both conditions; * = signrank p < 0.05 in
% that bin (uncorrected). Numbers under each curve = flies per bin.
figure('color','w','Position',[80 80 1300 700]); clf
tl2 = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
for row = 1:2
    for si = 1:3
        ax = nexttile(tl2,(row-1)*3+si); hold(ax,'on')
        Y = R(si).(['curve_' metric{row}]);
        hd = gobjects(1,2);
        y_lo = inf; y_hi = -inf; curves = cell(1,2);
        for rc = 1:2
            ref = tests(rc,2);
            d = squeeze(Y(:,2,:) - Y(:,ref,:)); if isvector(d), d = d(:); end
            n_ok = sum(~isnan(d),2);
            mu = mean(d,2,'omitnan'); se = std(d,0,2,'omitnan')./sqrt(n_ok);
            mu(n_ok < min_flies_per_bin) = nan; se(n_ok < min_flies_per_bin) = nan;
            v = ~isnan(mu);
            patch(ax,[speed_x(v) fliplr(speed_x(v))],[mu(v)'+se(v)' fliplr(mu(v)'-se(v)')],cond_rgb(ref,:),'FaceAlpha',.2,'EdgeColor','none')
            hd(rc) = plot(ax,speed_x,mu,'-o','Color',cond_rgb(ref,:),'MarkerFaceColor',cond_rgb(ref,:),'MarkerSize',4,'LineWidth',2);
            curves{rc} = struct('d',d,'mu',mu,'se',se,'n',n_ok,'v',v);
            if any(v), y_lo = min(y_lo,min(mu(v)-se(v))); y_hi = max(y_hi,max(mu(v)+se(v))); end
        end
        pad = 0.25*(y_hi-y_lo);
        ylim(ax,[y_lo-2*pad, y_hi+pad])
        for rc = 1:2
            C = curves{rc}; ref = tests(rc,2);
            for b = find(C.v)'
                if signrank(C.d(b,~isnan(C.d(b,:)))) < 0.05
                    % freeze - moving cue stars below its curve, freeze - dark above, so they don't land on the other curve
                    if rc == 1
                        text(ax,speed_x(b),C.mu(b)-C.se(b),'*','Color',cond_rgb(ref,:),'HorizontalAlignment','center','VerticalAlignment','top','FontSize',13)
                    else
                        text(ax,speed_x(b),C.mu(b)+C.se(b),'*','Color',cond_rgb(ref,:),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',13)
                    end
                end
                % flies per bin, two text rows at the bottom of the panel
                text(ax,speed_x(b),y_lo-(0.3+0.6*rc)*pad,sprintf('%d',C.n(b)),'Color',cond_rgb(ref,:),'HorizontalAlignment','center','FontSize',7)
            end
        end
        yline(ax,0,':k'); xlim(ax,[0 2.5])
        xlabel(ax,'rotational speed (rad/s)')
        ylabel(ax,sprintf('freeze - reference\n(%s dF/F)',metric{row}))
        if row == 1, title(ax,sensor_lbl{si},'FontSize',10), end
        if row == 1 && si == 1, legend(ax,hd,{'freeze - moving cue','freeze - dark'},'Location','northeast'), end
    end
end
title(tl2,{'Speed dependence of the freeze effect: paired per-fly difference in each speed bin', ...
    sprintf('mean \\pm sem across flies (>= %d flies with both conditions); * signrank p < 0.05 in that bin, uncorrected; numbers = flies per bin',min_flies_per_bin)},'FontSize',10)
exportgraphics(gcf,fullfile(export_dir,'fig_freeze_fluor_vs_rotation_by_speed.png'),'Resolution',200)

%% ===================== local functions =====================
function [frz, amb] = find_freezes(cue, r, fs, still_step, bridge_n, min_s, min_range, amb_s)
    % frz: samples inside a cue-still run that is a confirmed freeze (fly's
    % heading spans >= min_range during it); amb: still runs that don't
    % qualify but last >= amb_s
    n = numel(cue);
    still = [abs(angle(exp(1i*diff(cue)))) < still_step; true];
    % bridge cue-log spikes inside still runs: a short moving run after which
    % the cue is back where it started. (A plain length rule would bridge
    % everything -- in closed loop the cue moves in isolated single-sample
    % one-pixel steps.)
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

function m = bin_by_speed(cr, y, edges, dt, min_sec)
    m = nan(numel(edges)-1,1);
    for j = 1:numel(edges)-1
        idx = cr >= edges(j) & cr < edges(j+1);
        if sum(idx)*dt >= min_sec, m(j) = mean(y(idx),'omitnan'); end
    end
end

function peak3 = best3_neighbor_peak(x)
    % as in grab_freeze_scratch.m: hemisphere-average the 32 glomeruli into a
    % 16-ring, average each 3 circularly adjacent glomeruli, take the max
    nh = size(x,1)/2;
    xh = (x(1:nh,:) + x(nh+1:end,:))/2;
    a3 = (xh + circshift(xh,-1,1) + circshift(xh,-2,1))/3;
    peak3 = max(a3,[],1);
end

function is_flash = detect_flash_frames(f, mad_k)
    frame_mean = mean(f,1);
    frame_cv   = std(f,0,1) ./ frame_mean;
    is_bright  = frame_mean > median(frame_mean) + mad_k*mad(frame_mean,1);
    is_uniform = frame_cv <= median(frame_cv);
    is_flash   = is_bright & is_uniform;
end

function ind = trial_indicator(meta_path)
    if contains(meta_path,'_EPG_syt7f','IgnoreCase',true)
        ind = 'GRAB(DA2m)'; % mislabeled at acquisition (20221117-7/8/9, 20221121-2/3)
    elseif contains(meta_path,'syt7f','IgnoreCase',true)
        ind = 'syt7f';
    elseif contains(meta_path,'syt8m','IgnoreCase',true)
        ind = 'syt8m';
    elseif contains(meta_path,'dlight','IgnoreCase',true)
        ind = 'dLight';
    elseif contains(meta_path,'grab','IgnoreCase',true) || contains(meta_path,'DA2m','IgnoreCase',true)
        ind = 'GRAB(DA2m)';
    else
        ind = '';
    end
end

function name = trial_folder_name(meta_path)
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        name = parts{end-1};
    else
        name = parts{end};
    end
end

function fid = trial_fly_id(dataset_name, meta_path)
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    if strcmp(dataset_name,'epg_dlight')
        fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+$','once')));
        if isempty(fly_part), fid = meta_path; else, fid = strjoin(parts(1:fly_part(1)),filesep); end
    else % lpsp_cl: date + trailing suffix digit
        tname    = trial_folder_name(meta_path);
        date_str = regexp(tname,'^\d{8}','match','once');
        suffix   = regexp(tname,'_(\d+)$','tokens','once');
        if isempty(suffix), fid = date_str; else, fid = [date_str,'_',suffix{1}]; end
    end
end

function idx = build_csv_index(base_dir)
    % from grab_freeze_scratch.m
    d = dir(fullfile(base_dir,'**','csv','trialSettings.csv'));
    idx.names = cell(numel(d),1); idx.paths = cell(numel(d),1);
    for i = 1:numel(d)
        [~,trial_name] = fileparts(fileparts(d(i).folder));
        idx.names{i} = trial_name;
        idx.paths{i} = fullfile(d(i).folder,d(i).name);
    end
    is_test = cellfun(@(n) contains(n,'flash_test','IgnoreCase',true), idx.names);
    idx.names = idx.names(~is_test); idx.paths = idx.paths(~is_test);
    idx.by_name = containers.Map('KeyType','char','ValueType','any');
    for i = 1:numel(idx.names)
        key = lower(idx.names{i});
        if isKey(idx.by_name,key), idx.by_name(key) = [idx.by_name(key), idx.paths(i)];
        else, idx.by_name(key) = idx.paths(i); end
    end
end

function res = resolve_csv_for_trial(trial_dir, tname, idx)
    % from grab_freeze_scratch.m: exact folder-name match, else the
    % ficTracData_DAQ base id; 'ambiguous' if candidates disagree
    res.info = read_trial_condition(''); res.status = 'not_found';
    key = lower(tname); used_fallback = false; candidates = {};
    if isKey(idx.by_name,key)
        candidates = idx.by_name(key);
    else
        used_fallback = true;
        dd = dir(fullfile(trial_dir,'*ficTracData_DAQ*'));
        if ~isempty(dd)
            base_id = regexp(dd(1).name,'^(\d{8}-\d+)','match','once');
            if ~isempty(base_id)
                hit = strcmp(idx.names,base_id) | startsWith(idx.names,[base_id,'_']);
                candidates = idx.paths(hit);
            end
        end
    end
    if isempty(candidates), return, end
    infos = cellfun(@read_trial_condition, candidates, 'UniformOutput', false);
    agree = all(cellfun(@(x) isequal(x.is_dark,infos{1}.is_dark) && strcmp(x.pattern,infos{1}.pattern), infos));
    res.info = infos{1};
    if ~agree, res.status = 'ambiguous'; elseif used_fallback, res.status = 'fallback'; else, res.status = 'exact'; end
end

function info = read_trial_condition(csv_path)
    info.pattern = ''; info.is_dark = nan;
    if isempty(csv_path) || ~isfile(csv_path), return, end
    T = readtable(csv_path,'TextType','string');
    if ismember('patternPath',T.Properties.VariableNames) && height(T) > 0 && strlength(T.patternPath(1)) > 0
        info.pattern = char(T.patternPath(1));
        info.is_dark = contains(info.pattern,'background','IgnoreCase',true);
    end
end
