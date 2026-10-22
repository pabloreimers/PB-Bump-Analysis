%% fig_rnai_velocity_gain
% "Velocity gain" from data scripts/lpsp_rnai_claude_v2.m step 8 (figure
% 30) -- an instantaneous-velocity regression, as opposed to the sliding-
% window "integrative gain" of fig_rnai_integrative_gain.m -- for the TH
% set (lpsp>th, empty>th, lpsp>mcherry) and the VGlut set (lpsp>vglut,
% empty>vglut, lpsp>mcherry), one column each.
%
% Method, exactly as v2 step 8 (settled in gain_scratch_claude.m):
%   bump:  im.f movmean 3 frames -> z-score per glomerulus -> PVA (mu,rho)
%          -> unwrap -> 0.2 s Gaussian -> onto fictrac timebase (imaging
%          frames spread evenly over ft.xf, as v2) -> d/dt  = bump_vel
%   fly:   ft.r_speed, 1.0 s Gaussian                     = fly_vel
%   cue:   -ft.cue (NaN dropouts set to pi), unwrapped, 0.75 s Gaussian,
%          d/dt                                           = cue_vel
%   lag:   fly_vel/cue_vel shifted 8 fictrac frames (~0.13 s) earlier
%   samples kept: 0.2 < |fly_vel| < 5 rad/s, |bump_vel| < 10 rad/s, rho > 0.2
%   per fly x light condition, pooling that fly's trials: OLS slope of
%          bump_vel ~ 1 + fly_vel (gain_fly, closed-loop target 0.8) and
%          bump_vel ~ 1 + cue_vel (gain_cue, closed-loop target 1, ~0
%          expected in the dark); NaN if < 50 kept samples or the fly
%          turned above threshold in < 1% of all samples.
% All trials of each fly count (v2 applies no walking/cue/bump trial
% filter, unlike the integrative-gain recipe).
% Run from the repo root (needs all_data; loads it if absent).

%% parameters
data_file = fullfile('.data','lpsp_rnai_joint_nosmooth_reg_20260831.mat');
sets = struct('label',{'TH set','VGlut set'}, ...
              'geno',{{'lpsp>th','empty>th','lpsp>mcherry'},{'lpsp>vglut','empty>vglut','lpsp>mcherry'}}, ...
              'names',{{'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > mCherry'},{'LPsP > VGlut-RNAi','Empty > VGlut-RNAi','LPsP > mCherry'}});
geno_all    = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};
geno_colors = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13]; % as lpsp_rnai_turns_claude.m
im_f_frames = 3; bump_smooth_s = 0.20; cue_smooth_s = 0.75; vel_smooth_s = 1.00; lag_frames = 8;
vel_thresh = 0.2; bump_thresh = 10; rho_thresh = 0.2; vel_max = 5;

%% labels + per fly x light condition gain
if ~exist('all_data','var')
    tmp = load(data_file,'all_data'); all_data = tmp.all_data(:); clear tmp
end
n_trials = numel(all_data);
genotype = cell(n_trials,1); is_dark = false(n_trials,1); fly_id = cell(n_trials,1);
for i = 1:n_trials
    [fly_seg,trial_seg] = meta_fly_and_trial_seg(all_data(i).meta);
    combo = lower(regexprep([fly_seg,trial_seg],'[_\s]',''));
    if contains(combo,'lpsp'), drv = 'lpsp'; elseif contains(combo,'empty'), drv = 'empty'; else, drv = ''; end
    if contains(combo,'vglut'), tgt = 'vglut'; elseif contains(combo,'mcherry'), tgt = 'mcherry'; else, tgt = 'th'; end
    genotype{i} = [drv '>' tgt];
    is_dark(i) = contains(all_data(i).ft.pattern,'background');
    fly_id{i} = meta_fly_id(all_data(i).meta);
end
[~,~,fly_num] = unique(fly_id);
n_flies = max(fly_num);
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = genotype{find(fly_num==f,1)}; end

pva_mu = cell(n_trials,1); pva_rho = cell(n_trials,1);
for i = find(ismember(genotype,geno_all))'
    [pva_mu{i},pva_rho{i}] = trial_pva_movmean(all_data(i),im_f_frames);
end
gain_fly = nan(n_flies,2); gain_cue = nan(n_flies,2); % columns: closed loop, dark
for f = 1:n_flies
    for d = 0:1
        tr_i = find(fly_num==f & is_dark==d);
        if isempty(tr_i), continue, end
        [gain_cue(f,d+1),gain_fly(f,d+1)] = fly_gain_v3(all_data,tr_i,pva_mu,pva_rho,bump_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max,true);
    end
end

%% stats printout
for d = 1:2
    for k = 1:numel(geno_all)
        y = gain_fly(strcmp(fly_geno,geno_all{k}),d); y = y(~isnan(y));
        fprintf('%-13s %-11s gain_fly %.2f +/- %.2f (n=%d)\n',geno_all{k},char(string(d==1).replace("true","closed loop").replace("false","dark")),mean(y),std(y)/sqrt(numel(y)),numel(y));
    end
end
pairs = {'lpsp>th','empty>th'; 'lpsp>th','lpsp>mcherry'; 'lpsp>vglut','empty>vglut'; 'lpsp>vglut','lpsp>mcherry'};
for p = 1:size(pairs,1)
    [~,pv] = ttest2(gain_fly(strcmp(fly_geno,pairs{p,1}),1),gain_fly(strcmp(fly_geno,pairs{p,2}),1));
    fprintf('closed-loop gain_fly t-test %s vs %s: p = %.3f\n',pairs{p,:},pv);
end

%% figure: rows = gain_fly, gain_cue; columns = TH set, VGlut set
figure('Color','w','Position',[60,60,1250,780]); clf
tl = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
metrics = {gain_fly, 'gain_{fly} = slope(bump vel ~ fly vel)', 0.8;
           gain_cue, 'gain_{cue} = slope(bump vel ~ cue vel)', 1};
ax = gobjects(2,2);
for r = 1:2
    for s = 1:2
        gs = sets(s).geno; n_g = numel(gs);
        ax(r,s) = nexttile(tl,(r-1)*2+s); hold on
        xt_lbl = {};
        for d = 1:2
            for gi = 1:n_g
                c = geno_colors(strcmp(geno_all,gs{gi}),:);
                x0 = (d-1)*(n_g+1) + gi;
                y = metrics{r,1}(strcmp(fly_geno,gs{gi}),d); y = y(~isnan(y));
                if ~isempty(y)
                    scatter(x0 + (rand(size(y))-.5)*.3,y,18,c,'filled','MarkerFaceAlpha',.4)
                    errorbar(x0,mean(y),std(y)/sqrt(numel(y)),'o','Color',c*.6,'MarkerFaceColor',c*.6,'LineWidth',2,'MarkerSize',7)
                end
                xt_lbl{end+1} = sprintf('%s (n=%d)',sets(s).names{gi},numel(y)); %#ok<AGROW>
            end
        end
        yline(metrics{r,3},':k'); yline(0,'-','Color',[.85,.85,.85])
        xticks([1:n_g, n_g+1+(1:n_g)]); xticklabels(xt_lbl); xtickangle(30); xlim([0, 2*n_g+2])
        text(mean(1:n_g)/(2*n_g+2),0.99,'closed loop','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
        text((n_g+1+mean(1:n_g))/(2*n_g+2),0.99,'dark','Units','normalized','HorizontalAlignment','center','VerticalAlignment','top','FontWeight','bold')
        ylabel(metrics{r,2}); if r == 1, title(sets(s).label), end
    end
    linkaxes(ax(r,:),'y')
end
title(tl,{'velocity gain (lpsp\_rnai\_claude\_v2.m step 8): per-fly OLS slope of bump velocity on fly / cue velocity, |fly vel| > 0.2 rad/s, rho > 0.2, 0.13 s lag', ...
    'dots = flies, circles = mean \pm SEM; dotted = target (0.8 fly, 1 cue)'})

%% export
out_dir = fullfile('ugly_figures','exports');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_velocity_gain.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_velocity_gain.png'),'Resolution',200)

%% functions (copied unchanged from lpsp_rnai_claude_v2.m / lpsp_rnai_turns_claude.m)
function [mu_new,rho_new,f_z] = trial_pva_movmean(trial,smooth_frames)
    f_smooth = movmean(trial.im.f,smooth_frames,2);
    f_z = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);
    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp] = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function [gain_cue,gain_fly] = fly_gain_v3(all_data,trial_list,mu_cache,rho_cache,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max,use_intercept)
    fly_vel = []; cue_vel = []; bump_vel = []; valid = logical([]);
    for k = 1:numel(trial_list)
        ti = trial_list(k);
        [fv,cv,bv,vd] = trial_gain_vectors_v3(all_data(ti),mu_cache{ti},rho_cache{ti},mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
            vel_thresh,bump_thresh,rho_thresh,vel_max);
        fly_vel = [fly_vel; fv]; cue_vel = [cue_vel; cv]; bump_vel = [bump_vel; bv]; valid = [valid; vd]; %#ok<AGROW>
    end
    min_frac_moving = 0.01;
    n_valid = sum(valid);
    if n_valid < 50 || mean(abs(fly_vel) > vel_thresh) < min_frac_moving
        gain_cue = nan; gain_fly = nan; return
    end
    if use_intercept
        b_cue = [ones(n_valid,1),cue_vel(valid)] \ bump_vel(valid); gain_cue = b_cue(2);
        b_fly = [ones(n_valid,1),fly_vel(valid)] \ bump_vel(valid); gain_fly = b_fly(2);
    else
        if std(cue_vel(valid)) < 0.01 || std(fly_vel(valid)) < 0.01
            gain_cue = nan; gain_fly = nan; return
        end
        gain_cue = cue_vel(valid) \ bump_vel(valid);
        gain_fly = fly_vel(valid) \ bump_vel(valid);
    end
end

function [fly_vel,cue_vel,bump_vel,valid] = trial_gain_vectors_v3(trial,mu_raw,rho_raw,mu_smooth_s,cue_smooth_s,vel_smooth_s,lag_frames, ...
        vel_thresh,bump_thresh,rho_thresh,vel_max)
    xf = trial.ft.xf;
    dt = median(diff(xf));
    n_im = numel(mu_raw);
    xb = linspace(xf(1),xf(end),n_im)';
    dt_im = (xf(end)-xf(1)) / (n_im-1);
    win_im = max(1,round(mu_smooth_s/dt_im));
    mu_smoothed = smoothdata(unwrap(mu_raw(:)),'gaussian',win_im);
    bump_vel_full = gradient(interp1(xb,mu_smoothed,xf,'linear','extrap'))/dt;
    rho_full = interp1(xb,rho_raw(:),xf,'linear','extrap');
    win_vel = max(1,round(vel_smooth_s/dt));
    fly_vel_full = smoothdata(trial.ft.r_speed(:),'gaussian',win_vel);
    win_cue = max(1,round(cue_smooth_s/dt));
    cue_filled = trial.ft.cue(:); cue_filled(isnan(cue_filled)) = pi; % fill_nan_gaps_pi
    cue_smoothed = smoothdata(unwrap(-cue_filled),'gaussian',win_cue);
    cue_vel_full = gradient(cue_smoothed)/dt;
    fly_vel  = fly_vel_full(1:end-lag_frames);
    cue_vel  = cue_vel_full(1:end-lag_frames);
    bump_vel = bump_vel_full(lag_frames+1:end);
    rho_i    = rho_full(lag_frames+1:end);
    valid = abs(fly_vel) > vel_thresh & abs(fly_vel) < vel_max & abs(bump_vel) < bump_thresh & rho_i > rho_thresh;
end

function [fly_seg,trial_seg] = meta_fly_and_trial_seg(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        trial_seg = parts{end-1}; fly_seg = parts{end-2};
    else
        trial_seg = parts{end};   fly_seg = parts{end-1};
    end
end

function fid = meta_fly_id(meta_path)
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        fly_part_end = numel(parts)-2;
    else
        fly_part_end = numel(parts)-1;
    end
    fid = strjoin(parts(1:fly_part_end),'\');
end
