%% fig_pb_folded_profile_by_speed
% Average PB profile (z-scored dF/F, left/right halves folded together) as a
% function of signed rotational speed, LPsP>syt7f vs EPG>GRAB(DA2m)
% (lpsp_cl_data_20240206.mat, closed loop + dark pooled).
%
% Per trial:
%   - z = each bin's dF/F z-scored over time (as im.z in process_im), but
%     recomputed with brief global artifact flashes excluded (they inflate
%     the SD). Because z is per bin, each glomerulus's resting offset is
%     removed: profiles show which glomeruli RISE more with speed, in units
%     of their own variability -- not absolute brightness.
%   - fold: bins 1-16 (one half) and 17-32 (other half) both run -pi..pi in
%     the same direction (alpha), so the halves are averaged non-mirrored.
%     The second half is shifted by the integer k in -3..3 bins that best
%     matches the first half's spatial structure (per-frame spatial mean
%     removed), WITHOUT wrapping (PB ends aren't adjacent; overhanging bins
%     keep a single half). Mirrored folding was also tested and loses in
%     every trial; the best k is consistently -1 to -2 (a ~1 glomerulus
%     offset between halves, possibly anatomical rather than mask error).
%   - signed |vr| (FicTrac r_speed, as recorded) shifted by the sensor's
%     lag (xcorr peak of PB-mean vs |vr|), binned; mean folded profile per
%     bin (>= 20 frames).
% Then trials -> flies (date_flyN) -> mean across flies.
%
% Result when written: both sensors rise across the PB with speed. GRAB is
% ~3x more compartmentalized at fast turns (spatial CV of the folded
% profile ~0.30 vs ~0.10; two humps at folded glomeruli ~3-4 and ~12-13,
% falling at the ends), LPsP calcium is nearly flat. Both are lowest at the
% mask ends (folded 15-16, and 1 for GRAB) -- partial end bins or real.
% A weak turn-direction gradient (-vr favours folded 2-8, +vr 9-15) is
% marginal (LPsP p~0.05, GRAB p~0.13); no left-vs-right hemisphere
% asymmetry with turn direction in either sensor.

%% load
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_data_20240206.mat')), data_dir = '.data'; else, data_dir = 'data'; end
d = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); % 'EPG_syt7f' folders are a naming typo for GRAB(DA2m)
isLpsp = contains(meta,'LPsP_syt7f');
grp    = {find(isLpsp), find(isGrab)};
names  = {'LPsP > syt7f','EPG > GRAB(DA2m)'};
c      = [0.85 0.33 0.10; 0.00 0.45 0.74];
fr     = 9.1;
shifts = -3:3;
spE = [-inf -3 -2 -1.25 -0.5 -0.15 0.15 0.5 1.25 2 3 inf]; nb = numel(spE)-1;
spc = [-3.5 -2.5 -1.6 -0.9 -0.3 0 0.3 0.9 1.6 2.5 3.5];   % nominal bin centres for plotting
min_frames = 20;
flyid = @(m) regexprep(regexp(m,'\d{8}-\d+_.*$','match','once'),'^(\d{8})-\d+_.*_(\d+)$','$1_$2');

%% per trial: fold + bin by speed
PR = cell(1,2); FLY = cell(1,2); SH = cell(1,2); HF = cell(1,2);
for g = 1:2
    N = numel(grp{g}); P = nan(16,nb,N); H = nan(N,nb); fl = cell(N,1); sh = nan(N,1);
    for n = 1:N
        i = grp{g}(n); s = d(i); D = s.im.d; T = size(D,2);
        m = mean(D,1); r = m - movmedian(m,21); bad = r > 6*1.4826*mad(r,1); bad = conv(double(bad),ones(1,7),'same') > 0;
        Z = (D - mean(D(:,~bad),2)) ./ std(D(:,~bad),0,2); Z(:,bad) = nan;

        L = Z(1:16,:); Rr = Z(17:32,:); sc = nan(size(shifts));
        for j = 1:numel(shifts)
            Rk = shiftpad(Rr,shifts(j)); ov = any(~isnan(Rk),2);
            A = L(ov,~bad); B = Rk(ov,~bad); A = A - mean(A,1); B = B - mean(B,1);
            sc(j) = corr(A(:),B(:));
        end
        [~,j] = max(sc); sh(n) = shifts(j);
        Zf = mean(cat(3,L,shiftpad(Rr,sh(n))),3,'omitnan');

        xf = s.ft.xf; xb = linspace(xf(1),xf(end),T);
        sp = interp1(xf,abs(smoothdata(s.ft.r_speed,'gaussian',30)),xb);
        lg = -round(1*fr):round(3*fr); x = mean(D,1); x(bad) = mean(x(~bad));
        xc = xcorr(x-mean(x), sp-mean(sp), lg(end), 'coeff'); [~,jl] = max(xc(lg+lg(end)+1));
        vr = circshift(interp1(xf,smoothdata(s.ft.r_speed,'gaussian',30),xb), lg(jl));
        b = discretize(vr,spE);
        LR = mean(Z(1:16,:),1) - mean(Z(17:32,:),1);   % unfolded hemisphere difference
        for k = 1:nb
            sel = b==k & ~any(isnan(Zf),1);
            if sum(sel) >= min_frames, P(:,k,n) = mean(Zf(:,sel),2); H(n,k) = mean(LR(sel),'omitnan'); end
        end
        fl{n} = flyid(meta{i});
    end
    [uf,~,fi] = unique(fl);
    PF = nan(16,nb,numel(uf)); HH = nan(numel(uf),nb);
    for f = 1:numel(uf)
        PF(:,:,f) = mean(P(:,:,fi==f),3,'omitnan'); HH(f,:) = mean(H(fi==f,:),1,'omitnan');
    end
    PR{g} = PF; FLY{g} = uf; SH{g} = sh; HF{g} = HH;
    fprintf('%-18s %d trials, %d flies; best fold shift counts (k = %s): %s\n', names{g}, N, numel(uf), mat2str(shifts), mat2str(histcounts(sh,[shifts-0.5 shifts(end)+0.5])));
end

%% stats
CVF = cell(1,2);
for g = 1:2
    PF = PR{g}; S = PF - mean(PF,1);
    fast = [1 2 nb-1 nb];
    CVF{g} = squeeze(mean(std(PF(:,fast,:),0,1)./mean(PF(:,fast,:),1),2,'omitnan'));
    Dp  = squeeze(mean(S(:,nb-2:nb,:),2,'omitnan') - mean(S(:,1:3,:),2,'omitnan'));
    idx = mean(Dp(9:15,:),1) - mean(Dp(2:8,:),1);
    asym = mean(HF{g}(:,nb-2:nb),2,'omitnan') - mean(HF{g}(:,1:3),2,'omitnan');
    fprintf('%-18s spatial CV at fast turns %.3f +/- %.3f | direction gradient (+fast - -fast, folded 9-15 minus 2-8) %.3f +/- %.3f, %d/%d flies >0, signrank p=%.2g | hemisphere asymmetry %.3f, p=%.2g\n', ...
        names{g}, mean(CVF{g},'omitnan'), std(CVF{g},'omitnan')/sqrt(sum(~isnan(CVF{g}))), mean(idx,'omitnan'), std(idx,'omitnan')/sqrt(sum(~isnan(idx))), ...
        sum(idx>0), sum(~isnan(idx)), signrank(idx), mean(asym,'omitnan'), signrank(asym));
end
fprintf('spatial CV at fast turns, ranksum p = %.2g\n', ranksum(CVF{1},CVF{2}));

%% figure
cmap = [linspace(0,.85,5)' linspace(.3,.85,5)' ones(5,1); .3 .3 .3; ones(5,1) linspace(.85,.3,5)' linspace(.85,0,5)'];
figure(1); clf; set(gcf,'Color','w','Position',[50 50 1600 820]); tl = tiledlayout(2,5,'TileSpacing','compact');
for g = 1:2
    M = mean(PR{g},3,'omitnan'); S = M - mean(M,1); nF = numel(FLY{g});
    nexttile; imagesc(1:nb,1:16,M); colorbar; xticks(1:nb); xticklabels(compose('%.1f',spc)); xtickangle(60)
    ylabel({names{g},'folded glomerulus'}); xlabel('signed rotational speed (rad/s)'); title(sprintf('mean z (%d flies)',nF))
    nexttile; hold on
    for k = 1:nb, plot(1:16,M(:,k),'-','Color',cmap(k,:),'LineWidth',1.5); end
    xlim([1 16]); xlabel('folded glomerulus'); ylabel('mean z'); title('profile per speed bin (blue -vr, grey still, red +vr)')
    nexttile; imagesc(1:nb,1:16,S); colorbar; xticks(1:nb); xticklabels(compose('%.1f',spc)); xtickangle(60)
    xlabel('signed rotational speed (rad/s)'); title('shape: each column minus its mean')
    nexttile; sdf = squeeze(std(PR{g} - mean(PR{g},1),0,1));
    errorbar(spc, mean(sdf,2,'omitnan'), std(sdf,0,2,'omitnan')./sqrt(sum(~isnan(sdf),2)), 'o-','Color',c(g,:),'LineWidth',1.5)
    xlabel('signed rotational speed (rad/s)'); ylabel('spatial SD across glomeruli (z)'); title('compartment contrast vs speed')
    nexttile; hold on
    for f = 1:nF, plot(1:16, PR{g}(:,1,f) - mean(PR{g}(:,1,f)) , '-', 'Color',[cmap(1,:) .4]); plot(1:16, PR{g}(:,nb,f) - mean(PR{g}(:,nb,f)), '-', 'Color',[cmap(nb,:) .4]); end
    plot(1:16, S(:,1),'-','Color',cmap(1,:),'LineWidth',2.5); plot(1:16, S(:,nb),'-','Color',cmap(nb,:),'LineWidth',2.5)
    xlim([1 16]); xlabel('folded glomerulus'); ylabel('shape (z)'); title('fastest bins per fly (thin) and mean (thick)')
end
title(tl, sprintf('folded PB profile vs signed rotational speed -- spatial CV at fast turns: LPsP %.2f, GRAB %.2f (p = %.1g)', ...
    mean(CVF{1},'omitnan'), mean(CVF{2},'omitnan'), ranksum(CVF{1},CVF{2})))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_folded_profile_by_speed.png'),'Resolution',150);

%% functions
function Y = shiftpad(X,k)
% shift rows of X by k without wrapping; vacated rows are NaN
Y = nan(size(X)); n = size(X,1);
if k >= 0, Y(1+k:n,:) = X(1:n-k,:); else, Y(1:n+k,:) = X(1-k:n,:); end
end
