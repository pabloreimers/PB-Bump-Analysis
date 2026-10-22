%% fig_rnai_integrative_gain
% "Integrative gain" -- the gain that best maps the fly's heading onto the
% bump over a sliding 20 s window, estimated every second through the
% trial -- for the TH set (lpsp>th, empty>th, lpsp>mcherry) and the VGlut
% set (lpsp>vglut, empty>vglut, lpsp>mcherry), one panel each.
%
% Method: data scripts/lpsp_rnai.m ("calculate integrative gain") +
% lpsp_rnai_joint.m ("extract the integrative gain of each trial, group
% identical trials"). Those stored the result as all_data(i).gain in the
% older per-project datasets; the current joint dataset
% (lpsp_rnai_joint_nosmooth_reg_20260831.mat, same one the turn analyses
% use) doesn't carry it, so it is recomputed here, same recipe:
%   bump: im.mu with rho < 0.1 set to NaN, unwrapped, interpolated onto
%         the fictrac timebase (ft.xb -> ft.xf), 60-sample Gaussian
%   lag:  bump shifted 20 fictrac samples (~0.33 s) behind heading
%   per second t: window t-10..t+10 s; gain = argmin over k = 0:0.05:5 of
%         circ_var(circ_dist(bump, k x heading)); also the loss at the
%         minimum (v) and circ_var(heading) in the window (hv)
%   keep: hv >= 0.1 (fly rotated enough to constrain the gain) and v <= 1;
%         samples with |r_speed| < 0.5 AND rho < 0.5 dropped; t > 530 s
%         dropped; a trial needs > 500 surviving fictrac samples
%   trials: walking (> 20% of samples with f_speed > 0.1 or 2|r_speed| >
%         0.1), working cue (cue moved > 5 rad and tracks -0.8 x integrated
%         r_speed, circ corr > 0.5), and a bump (mean rho > 0.2)
%   per fly x light condition: mean over all surviving samples of its trials
% Vectorized (all candidate gains at once per window) but otherwise the
% same computation; the only deliberate difference is that each window is
% centered on xt(j) -- the original centered window j on j while storing
% it at xt(j) = j-1, a 1 s offset. Results are cached to
% .data/lpsp_rnai_intgain_cache.mat (gitignored) since the full pass takes
% a few minutes.
%
% Figure: top row, per-fly mean gain (closed loop | dark) by genotype, one
% panel per set; bottom row, closed-loop gain vs. time in trial (60 s
% bins), mean +/- SEM across flies. Run from the repo root.

%% parameters
addpath(fullfile(pwd,'circ_stats')); % circ_corrcc
data_file  = fullfile('.data','lpsp_rnai_joint_nosmooth_reg_20260831.mat');
cache_file = fullfile('.data','lpsp_rnai_intgain_cache.mat');
sets = struct('label',{'TH set','VGlut set'}, ...
              'geno',{{'lpsp>th','empty>th','lpsp>mcherry'},{'lpsp>vglut','empty>vglut','lpsp>mcherry'}}, ...
              'names',{{'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'},{'LPsP > VGlut-RNAi','Empty > VGlut-RNAi','LPsP > mCherry'}});
geno_all    = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};
geno_colors = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13]; % as lpsp_rnai_turns_claude.m
win = [-10,10]; lag = 20; k_grid = 0:.05:5; mu_rho_thresh = .1; mu_smooth = 60;
hv_thresh = .1; v_thresh = 1; r_thresh = .5; rho_thresh_s = .5; t_max = 530; min_samples = 500;
for_thresh = .1; for_length = .2; cue_thresh = 5; rho_thresh_trial = .2;
time_bin_s = 60;

%% compute (or load) per-trial integrative gain + trial labels
if exist(cache_file,'file')
    load(cache_file)
else
    if ~exist('all_data','var')
        tmp = load(data_file,'all_data'); all_data = tmp.all_data(:); clear tmp
    end
    n_trials = numel(all_data);
    genotype = cell(n_trials,1); is_dark = false(n_trials,1); fly_id = cell(n_trials,1);
    inc = false(n_trials,1); g_trial = cell(n_trials,1); xf_trial = cell(n_trials,1);
    tic
    for i = 1:n_trials
        tr = all_data(i);
        [fly_seg,trial_seg] = meta_fly_and_trial_seg(tr.meta);
        combo = lower(regexprep([fly_seg,trial_seg],'[_\s]',''));
        if contains(combo,'lpsp'), drv = 'lpsp'; elseif contains(combo,'empty'), drv = 'empty'; else, drv = ''; end
        if contains(combo,'vglut'), tgt = 'vglut'; elseif contains(combo,'mcherry'), tgt = 'mcherry'; else, tgt = 'th'; end
        genotype{i} = [drv '>' tgt];
        is_dark(i) = contains(tr.ft.pattern,'background');
        fly_id{i} = meta_fly_id(tr.meta);

        % trial inclusion, as lpsp_rnai_joint.m
        walk = sum((tr.ft.f_speed > for_thresh) | (abs(tr.ft.r_speed)*2 > for_thresh)) > numel(tr.ft.f_speed)*for_length;
        tmp2 = -cumsum(tr.ft.r_speed)*.8/60; ok = ~isnan(tr.ft.cue) & ~isnan(tmp2);
        cue_ok = sum(abs(diff(unwrap(tr.ft.cue))),'omitnan') > cue_thresh && circ_corrcc(tr.ft.cue(ok),tmp2(ok)) > .5;
        bump_ok = mean(tr.im.rho) > rho_thresh_trial;
        inc(i) = walk && cue_ok && bump_ok;

        % integrative gain, as lpsp_rnai.m
        m = tr.im.mu(:); m(tr.im.rho(:) < mu_rho_thresh) = nan;
        nf = 1:numel(m);
        m = interp1(tr.ft.xb(nf),unwrap(m),tr.ft.xf(:));
        m = smoothdata(m,1,'gaussian',mu_smooth);
        h = tr.ft.heading(:);
        m = m(lag+1:end); h = h(1:end-lag); xf = tr.ft.xf(1:end-lag); xf = xf(:);
        xt = (0:ceil(max(tr.ft.xf)))';
        g = nan(numel(xt),1); v = inf(numel(xt),1); hv = nan(numel(xt),1);
        % circ_var(circ_dist(m, k h)) = 1 - |mean exp(i(m - k h))| over the
        % window; the per-sample terms for every k are computed once, and
        % window means come from cumulative sums (NaN samples omitted, as
        % circ_var(...,'omitnan'))
        Z = exp(1i*(m - h*k_grid)); okZ = ~isnan(Z); Z(~okZ) = 0;
        cZ = [zeros(1,numel(k_grid)); cumsum(Z,1)]; cN = [zeros(1,numel(k_grid)); cumsum(okZ,1)];
        H = exp(1i*h); cH = [0; cumsum(H)];
        lo = arrayfun(@(t) find(xf > t+win(1),1,'first'),xt,'UniformOutput',false);
        hi = arrayfun(@(t) find(xf < t+win(2),1,'last'),xt,'UniformOutput',false);
        for j = 1:numel(xt)
            if isempty(lo{j}) || isempty(hi{j}) || hi{j} < lo{j}, continue, end
            a = lo{j}; b = hi{j};
            hv(j) = 1 - abs((cH(b+1)-cH(a))/(b-a+1));          % circ_var(heading in window)
            if isnan(m(a)), continue, end                       % original: m_tmp - m_tmp(1) all NaN -> no gain
            R = abs((cZ(b+1,:)-cZ(a,:)) ./ (cN(b+1,:)-cN(a,:)));
            [vmin,kk] = min(1 - R);
            if ~isnan(vmin), g(j) = k_grid(kk); v(j) = vmin; end
        end
        % sample-level filters, as lpsp_rnai_joint.m
        g(hv < hv_thresh | v > v_thresh) = nan;
        gs = interp1(xt,g,tr.ft.xf(:));
        rho = interp1(tr.ft.xb(nf),tr.im.rho(:),tr.ft.xf(:));
        gs(abs(tr.ft.r_speed(:)) < r_thresh & rho < rho_thresh_s) = nan;
        gs(tr.ft.xf(:) > t_max) = nan;
        if sum(~isnan(gs)) > min_samples, g_trial{i} = gs; xf_trial{i} = tr.ft.xf(:); end
        if mod(i,50) == 0, fprintf('%d/%d trials, %.0f s\n',i,n_trials,toc); end
    end
    [fly_list,~,fly_num] = unique(fly_id);
    save(cache_file,'genotype','is_dark','fly_num','fly_list','inc','g_trial','xf_trial','-v7.3');
end
use = inc & ~cellfun(@isempty,g_trial);
fprintf('%d of %d trials pass inclusion (walking, working cue, bump) and have > %d gain samples\n',sum(use),numel(use),min_samples);

%% per fly x light condition: mean gain, and closed-loop gain vs. time in trial
n_flies = max(fly_num);
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = genotype{find(fly_num==f,1)}; end
g_fly = nan(n_flies,2); % columns: closed loop, dark
t_edges = 0:time_bin_s:t_max+time_bin_s; n_tb = numel(t_edges)-1;
g_time = nan(n_flies,n_tb);
for f = 1:n_flies
    for d = 0:1
        tr_i = find(use & fly_num==f & is_dark==d);
        if isempty(tr_i), continue, end
        g_fly(f,d+1) = mean(vertcat(g_trial{tr_i}),'omitnan');
        if d == 0
            G = vertcat(g_trial{tr_i}); T = vertcat(xf_trial{tr_i});
            tb = discretize(T,t_edges);
            for b = 1:n_tb, g_time(f,b) = mean(G(tb==b),'omitnan'); end
        end
    end
end

%% figure
figure('Color','w','Position',[60,60,1250,780]); clf
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
ax_top = gobjects(1,2); ax_bot = gobjects(1,2);
for s = 1:2
    gs = sets(s).geno; n_g = numel(gs);
    gcol = zeros(n_g,3); for gi = 1:n_g, gcol(gi,:) = geno_colors(strcmp(geno_all,gs{gi}),:); end

    % top: per-fly mean gain, closed loop (left group) and dark (right group)
    ax_top(s) = nexttile(tl,s); hold on
    xt_lbl = {};
    for d = 1:2
        for gi = 1:n_g
            x0 = (d-1)*(n_g+1) + gi;
            y = g_fly(strcmp(fly_geno,gs{gi}),d); y = y(~isnan(y));
            if ~isempty(y)
                scatter(x0 + (rand(size(y))-.5)*.3,y,18,gcol(gi,:),'filled','MarkerFaceAlpha',.4)
                errorbar(x0,mean(y),std(y)/sqrt(numel(y)),'o','Color',gcol(gi,:)*.6,'MarkerFaceColor',gcol(gi,:)*.6,'LineWidth',2,'MarkerSize',7)
            end
            xt_lbl{end+1} = sprintf('%s (n=%d)',sets(s).names{gi},numel(y)); %#ok<AGROW>
        end
    end
    yline(0.8,':k')
    xticks([1:n_g, n_g+1+(1:n_g)]); xticklabels(xt_lbl); xtickangle(30)
    xlim([0, 2*n_g+2]) % normalized label positions above assume this range
    text(mean(1:n_g)/(2*n_g+2),0.99,'closed loop','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
    text((n_g+1+mean(1:n_g))/(2*n_g+2),0.99,'dark','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
    ylabel('integrative gain (per-fly mean)'); title(sets(s).label)

    % bottom: closed-loop gain vs. time in trial
    ax_bot(s) = nexttile(tl,2+s); hold on
    tc = t_edges(1:end-1) + time_bin_s/2;
    for gi = 1:n_g
        Y = g_time(strcmp(fly_geno,gs{gi}),:); Y = Y(any(~isnan(Y),2),:);
        m = mean(Y,1,'omitnan'); se = std(Y,0,1,'omitnan')./sqrt(max(1,sum(~isnan(Y),1)));
        ok = ~isnan(m);
        patch([tc(ok),fliplr(tc(ok))],[m(ok)+se(ok),fliplr(m(ok)-se(ok))],gcol(gi,:),'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off')
        plot(tc,m,'o-','Color',gcol(gi,:),'MarkerFaceColor',gcol(gi,:),'LineWidth',1.5,'DisplayName',sprintf('%s (n=%d flies)',sets(s).names{gi},size(Y,1)))
    end
    yline(0.8,':k','HandleVisibility','off')
    xlim([0,t_max]); xlabel('time in trial (s)'); ylabel('integrative gain (closed loop)')
    legend('Location','southeast'); title(sprintf('%s: closed-loop gain over the trial (%d s bins)',sets(s).label,time_bin_s))
end
linkaxes(ax_top,'y'); linkaxes(ax_bot,'y')
ylim(ax_top(1),[0,max(1.6,max(g_fly(:))*1.05)]);
title(tl,{'integrative gain: best-fit bump/heading gain in sliding 20 s windows (lpsp\_rnai.m method), trials with walking + working cue + bump', ...
    'dots = flies, circles = mean \pm SEM; dotted = 0.8 cue gain'})

%% export
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_integrative_gain.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_integrative_gain.png'),'Resolution',200)

%% functions
function [fly_seg,trial_seg] = meta_fly_and_trial_seg(meta_path)
    % same as lpsp_rnai_turns_claude.m
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        trial_seg = parts{end-1}; fly_seg = parts{end-2};
    else
        trial_seg = parts{end};   fly_seg = parts{end-1};
    end
end

function fid = meta_fly_id(meta_path)
    % same as lpsp_rnai_turns_claude.m
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        fly_part_end = numel(parts)-2;
    else
        fly_part_end = numel(parts)-1;
    end
    fid = strjoin(parts(1:fly_part_end),'\');
end
