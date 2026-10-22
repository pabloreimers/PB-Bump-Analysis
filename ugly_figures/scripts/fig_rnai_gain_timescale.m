%% fig_rnai_gain_timescale
% Why do the two gain metrics disagree? The velocity gain
% (fig_rnai_velocity_gain.m; lpsp_rnai_claude_v2.m step 8) is a
% moment-to-moment slope of bump velocity on fly velocity; the integrative
% gain (fig_rnai_integrative_gain.m; lpsp_rnai.m) is a position fit over
% 20 s windows. This bridges them with one family of estimates:
%
%   DISPLACEMENT GAIN AT TIMESCALE dt: OLS slope (with intercept) of bump
%   displacement mu(t+dt) - mu(t) on heading displacement h(t+dt) - h(t),
%   pooled per fly over its closed-loop trials, for dt from 0.1 s
%   (~ the velocity gain) to 20 s (~ the integrative gain's window).
%   Bump/heading exactly as the velocity gain (im.f movmean 3 -> z-score
%   -> PVA, 0.2 s Gaussian, 8-frame lag; heading = ft.heading), samples
%   every dt/2, kept if rho > 0.2 at both ends and the fly's mean turning
%   speed over the interval is 0.2..5 rad/s (the velocity gain's gate,
%   applied to the interval). Trials restricted to the integrative gain's
%   inclusion set (walking + working cue + bump), so trial selection is
%   held fixed across timescales.
%
% Plus two checks on the velocity gain itself: the same trial inclusion,
% and longer lags (0.13 / 0.25 / 0.4 s).
% Needs all_data and the integrative-gain cache (run
% fig_rnai_integrative_gain.m once). Run from the repo root.

%% parameters
data_file  = fullfile('.data','lpsp_rnai_joint_nosmooth_reg_20260831.mat');
cache_file = fullfile('.data','lpsp_rnai_intgain_cache.mat');
geno_all    = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};
geno_names  = {'LPsP > TH-RNAi','Empty > TH-RNAi','LPsP > VGlut-RNAi','Empty > VGlut-RNAi','LPsP > mCherry'};
geno_colors = [0.85,0.33,0.10; 0.00,0.45,0.74; 0.49,0.18,0.56; 0.47,0.67,0.19; 0.93,0.69,0.13];
set_geno = {[1,2,5],[3,4,5]}; set_label = {'TH set','VGlut set'};
dts = [0.1, 0.25, 0.5, 1, 2, 5, 10, 20]; % s
lag_frames = 8; vel_thresh = 0.2; vel_max = 5; rho_thresh = 0.2;

if ~exist('all_data','var')
    tmp = load(data_file,'all_data'); all_data = tmp.all_data(:); clear tmp
end
C = load(cache_file,'genotype','is_dark','fly_num','inc');
n_flies = max(C.fly_num);
fly_geno = cell(n_flies,1); for f = 1:n_flies, fly_geno{f} = C.genotype{find(C.fly_num==f,1)}; end
use = C.inc & ~C.is_dark & ismember(C.genotype,geno_all);

%% per-trial bump / heading on the fictrac timebase (velocity-gain pipeline)
bump = cell(numel(all_data),1); head = bump; rhoF = bump; dtf = nan(numel(all_data),1);
for i = find(use)'
    tr = all_data(i); xf = tr.ft.xf(:); dtf(i) = median(diff(xf));
    f_s = movmean(tr.im.f,3,2); f_z = (f_s - mean(f_s,2)) ./ std(f_s,0,2);
    [x_tmp,y_tmp] = pol2cart(tr.im.alpha(:)',f_z'); [mu,rho] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    n_im = numel(mu); xb = linspace(xf(1),xf(end),n_im)';
    mu = smoothdata(unwrap(mu),'gaussian',max(1,round(0.2/((xf(end)-xf(1))/(n_im-1)))));
    m = interp1(xb,mu,xf,'linear','extrap'); r = interp1(xb,rho,xf,'linear','extrap');
    h = tr.ft.heading(:); h = fillmissing(h,'linear','EndValues','nearest');
    bump{i} = m(lag_frames+1:end); rhoF{i} = r(lag_frames+1:end); head{i} = h(1:end-lag_frames); % bump lags behavior
end

%% displacement gain per fly x timescale
G = nan(n_flies,numel(dts));
for f = unique(C.fly_num(use))'
    tr_i = find(use & C.fly_num==f);
    for k = 1:numel(dts)
        dB = []; dH = [];
        for i = tr_i'
            n = round(dts(k)/dtf(i)); step = max(1,round(n/2));
            a = (1:step:numel(head{i})-n)'; b = a + n;
            dh = head{i}(b) - head{i}(a); db = bump{i}(b) - bump{i}(a);
            spd = abs(dh)/dts(k);
            ok = spd > vel_thresh & spd < vel_max & rhoF{i}(a) > rho_thresh & rhoF{i}(b) > rho_thresh & ~isnan(db);
            dB = [dB; db(ok)]; dH = [dH; dh(ok)]; %#ok<AGROW>
        end
        if numel(dH) < 30, continue, end
        beta = [ones(numel(dH),1),dH] \ dB; G(f,k) = beta(2);
    end
end

%% figure: displacement gain vs. timescale, one panel per set
figure('Color','w','Position',[60,60,1150,430]); clf
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for s = 1:2
    nexttile; hold on
    for gi = set_geno{s}
        Y = G(strcmp(fly_geno,geno_all{gi}),:); Y = Y(any(~isnan(Y),2),:);
        errorbar(dts,mean(Y,1,'omitnan'),std(Y,0,1,'omitnan')./sqrt(sum(~isnan(Y),1)),'o-','Color',geno_colors(gi,:), ...
            'MarkerFaceColor',geno_colors(gi,:),'LineWidth',1.5,'DisplayName',sprintf('%s (n=%d)',geno_names{gi},size(Y,1)))
    end
    yline(0.8,':k','HandleVisibility','off'); set(gca,'XScale','log'); xticks(dts); xlim([0.08,25])
    xlabel('timescale \Deltat (s)'); ylabel('slope of \Deltabump on \Deltaheading'); title(set_label{s}); legend('Location','southeast')
end
title(tl,{'closed loop: bump vs. heading displacement gain as a function of timescale (same trials as the integrative gain)', ...
    'short \Deltat \approx velocity gain, long \Deltat \approx integrative gain; mean \pm SEM across flies'})
out_dir = fullfile('ugly_figures','exports');
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_gain_timescale.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_gain_timescale.png'),'Resolution',200)

for gi = 1:numel(geno_all)
    Y = G(strcmp(fly_geno,geno_all{gi}),:);
    fprintf('%-13s %s\n',geno_all{gi},strjoin(arrayfun(@(k) sprintf('%gs %.2f',dts(k),mean(Y(:,k),'omitnan')),1:numel(dts),'UniformOutput',false),'  '));
end

%% gain by size of heading excursion (dt = 1 s)
% Same 1 s displacements, split by how far the fly turned in that second;
% through-origin slope of (sign-flipped) bump displacement on |heading
% displacement| per bin, per fly. Least-squares slopes like the velocity
% gain are dominated by the big-excursion bins (each sample weighted by
% its squared excursion); the integrative gain's bounded circular loss
% over 20 s windows weights every moment ~equally, so it mostly reflects
% the small-excursion bins.
exc_edges = [0.2, 0.5, 1, 2, 5]; % rad of heading change within 1 s
min_flies_mean = 7;              % group mean +/- SEM drawn only for bins with more flies than this (dots always shown)
n_eb = numel(exc_edges)-1;
Gb = nan(n_flies,n_eb);
for f = unique(C.fly_num(use))'
    dB = []; dH = [];
    for i = find(use & C.fly_num==f)'
        n = round(1/dtf(i)); a = (1:max(1,round(n/2)):numel(head{i})-n)'; b = a + n;
        dh = head{i}(b) - head{i}(a); db = bump{i}(b) - bump{i}(a);
        ok = rhoF{i}(a) > rho_thresh & rhoF{i}(b) > rho_thresh & ~isnan(db);
        dB = [dB; db(ok)]; dH = [dH; dh(ok)]; %#ok<AGROW>
    end
    x = abs(dH); y = dB.*sign(dH); eb = discretize(x,exc_edges);
    for k = 1:n_eb, w = eb==k; if sum(w) >= 20, Gb(f,k) = x(w)\y(w); end, end
end
exc_labels = arrayfun(@(k) sprintf('%.0f-%.0f%s',rad2deg(exc_edges(k)),rad2deg(exc_edges(k+1)),char(176)),1:n_eb,'UniformOutput',false);

figure('Color','w','Position',[60,60,1150,430]); clf
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
for s = 1:2
    nexttile; hold on
    n_g = numel(set_geno{s});
    for j = 1:n_g
        gi = set_geno{s}(j);
        Y = Gb(strcmp(fly_geno,geno_all{gi}),:); Y = Y(any(~isnan(Y),2),:);
        xj = (1:n_eb) + (j-(n_g+1)/2)*(0.7/n_g);
        plot(repmat(xj,size(Y,1),1)',Y','.','Color',[geno_colors(gi,:),.3],'MarkerSize',8,'HandleVisibility','off')
        n_b = sum(~isnan(Y),1); m = mean(Y,1,'omitnan'); se = std(Y,0,1,'omitnan')./sqrt(n_b);
        m(n_b <= min_flies_mean) = NaN; se(n_b <= min_flies_mean) = NaN; % too few flies for a group mean: dots only
        errorbar(xj,m,se,'o-','Color',geno_colors(gi,:), ...
            'MarkerFaceColor',geno_colors(gi,:),'LineWidth',2,'DisplayName',sprintf('%s (n=%d)',geno_names{gi},size(Y,1)))
    end
    yline(0.8,':k','HandleVisibility','off'); yline(0,'-','Color',[.85,.85,.85],'HandleVisibility','off')
    xticks(1:n_eb); xticklabels(exc_labels); xlim([0.5,n_eb+0.5]); ylim([-0.5,2])
    xlabel('heading change within 1 s'); ylabel('bump / heading displacement gain'); title(set_label{s}); legend('Location','northeast')
end
title(tl,{'closed loop: gain over 1 s displacements, split by how far the fly turned (same trials as the integrative gain)', ...
    sprintf('faint = flies, thick = mean \\pm SEM (bins with > %d flies only); dotted = 0.8 cue gain',min_flies_mean)})
for gi = 1:numel(geno_all), fprintf('%-13s flies per excursion bin: %s\n',geno_all{gi},mat2str(sum(~isnan(Gb(strcmp(fly_geno,geno_all{gi}),:)),1))); end
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_gain_by_excursion.pdf'),'ContentType','vector')
exportgraphics(gcf,fullfile(out_dir,'fig_rnai_gain_by_excursion.png'),'Resolution',200)
