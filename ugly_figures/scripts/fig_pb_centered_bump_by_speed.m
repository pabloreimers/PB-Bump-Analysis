%% fig_pb_centered_bump_by_speed
% Bump SHAPE (regardless of bump position) as a function of signed
% rotational speed: LPsP>syt7f, EPG>GRAB(DA2m) (lpsp_cl_data_20240206.mat)
% and, as a reference for a real EPG bump, EPG>GCaMP7f (epg_7f_20261007.mat).
% Closed loop + dark pooled.
%
% Per trial:
%   - z = each bin's dF/F z-scored over time (artifact frames excluded).
%   - the two PB halves (bins 1-16, 17-32; both span -pi..pi in alpha) are
%     each centered on the bump, then averaged (= fold + center).
%   - CROSS-CENTERING: each half is centered on the bump phase measured from
%     the OTHER half (PVA), corrected by that trial's left-right phase offset
%     (circular mean of phase_L - phase_R). Centering a profile on its own
%     PVA aligns noise peaks too and manufactures a bump even from pure
%     noise; the halves' noise is independent, so cross-centering doesn't.
%     (A self-centered version is computed too, for comparison.)
%   - centering is a fractional circular shift (linear interpolation on the
%     16-bin ring), so the bump peak sits at 0 rad without +/-0.5 bin jitter.
%   - all frames are used (no bump-strength threshold): frames with a weak
%     or absent bump contribute a ~flat profile, so the average reflects
%     how bump-like the signal is at that speed, not a pre-selected subset.
%   - signed rotational speed (r_speed as recorded), shifted by the sensor's
%     lag (xcorr peak of PB-mean vs |vr|), binned.
% Trials -> flies -> mean across flies. Amplitude = peak - trough of the
% centered profile; width = full width at half (peak - trough) height.

%% load
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_data_20240206.mat')), data_dir = '.data'; else, data_dir = 'data'; end
d = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); % 'EPG_syt7f' folders are a naming typo for GRAB(DA2m)
isLpsp = contains(meta,'LPsP_syt7f');
E = load(fullfile(data_dir,'epg_7f_20261007.mat'),'all_data'); E = E.all_data(:);

names = {'LPsP > syt7f','EPG > GRAB(DA2m)','EPG > GCaMP7f (reference)'};
c     = [0.85 0.33 0.10; 0.00 0.45 0.74; 0.3 0.3 0.3];
spE = [-inf -3 -2 -1.25 -0.5 -0.15 0.15 0.5 1.25 2 3 inf]; nb = numel(spE)-1;
spc = [-3.5 -2.5 -1.6 -0.9 -0.3 0 0.3 0.9 1.6 2.5 3.5];
min_frames = 20;
alpha16 = linspace(-pi, pi - 2*pi/16, 16)';
flyid_cl = @(m) regexprep(regexp(m,'\d{8}-\d+_.*$','match','once'),'^(\d{8})-\d+_.*_(\d+)$','$1_$2');

%% unify trials: {D (32 x T dF/F), vr (signed, imaging clock), fly}
U = {[],[],[]};
for g = 1:2
    ids = find((g==1 & isLpsp) | (g==2 & isGrab));
    for i = ids
        s = d(i); T = size(s.im.d,2); xf = s.ft.xf; xb = linspace(xf(1),xf(end),T);
        U{g} = [U{g}; struct('D',s.im.d,'vr',interp1(xf,smoothdata(s.ft.r_speed,'gaussian',30),xb),'fly',flyid_cl(meta{i}))]; %#ok<AGROW>
    end
end
for k = 1:numel(E)
    s = E(k); f = smoothdata(s.im.f,2,'gaussian',5); f0 = prctile(f,7,2);
    parts = strsplit(strrep(s.meta.flyDir,'\','/'),'/'); parts = parts(~cellfun(@isempty,parts));
    if startsWith(parts{end},'fly'), fl = [parts{end-1} '_' strrep(parts{end},' ','')]; else, fl = parts{end}; end
    vr = interp1(s.ft.xf, smoothdata(s.ft.r_speed_raw,'gaussian',30), s.ft.xb);
    D = (f - f0)./f0; D(:,logical(s.im.badVol)) = nan;
    U{3} = [U{3}; struct('D',D,'vr',vr,'fly',fl)]; %#ok<AGROW>
end

%% per trial: center + bin
PRx = cell(1,3); PRs = cell(1,3); FLY = cell(1,3);
for g = 1:3
    N = numel(U{g}); Px = nan(16,nb,N); Ps = Px; fl = {U{g}.fly}';
    for n = 1:N
        D = U{g}(n).D; T = size(D,2);
        m = mean(D,1,'omitnan'); r = m - movmedian(m,21,'omitnan');
        bad = r > 6*1.4826*mad(r,1) | isnan(m); bad = conv(double(bad),ones(1,7),'same') > 0;
        Z = (D - mean(D(:,~bad),2)) ./ std(D(:,~bad),0,2);
        L = Z(1:16,:); R = Z(17:32,:);
        phL = angle(exp(1i*alpha16).' * L); phR = angle(exp(1i*alpha16).' * R);
        dLR = angle(mean(exp(1i*(phL(~bad) - phR(~bad)))));          % per-trial left-right offset
        Lx = center_ring(L, phR + dLR); Rx = center_ring(R, phL - dLR);  % cross-centered
        Ls = center_ring(L, phL);        Rs = center_ring(R, phR);         % self-centered (biased)
        Cx = (Lx + Rx)/2; Cs = (Ls + Rs)/2;

        sp = abs(U{g}(n).vr); x = m; x(bad) = mean(x(~bad)); fr_ = 9.1;
        lg = -round(1*fr_):round(3*fr_); sp0 = sp; sp0(isnan(sp0)) = 0;
        xc = xcorr(x-mean(x), sp0-mean(sp0), lg(end), 'coeff'); [~,jl] = max(xc(lg+lg(end)+1));
        b = discretize(circshift(U{g}(n).vr, lg(jl)), spE);
        for k = 1:nb
            sel = b==k & ~bad;
            if sum(sel) >= min_frames
                Px(:,k,n) = mean(Cx(:,sel),2); Ps(:,k,n) = mean(Cs(:,sel),2);
            end
        end
    end
    [uf,~,fi] = unique(fl);
    PRx{g} = cell2mat(reshape(arrayfun(@(f) mean(Px(:,:,fi==f),3,'omitnan'),1:numel(uf),'UniformOutput',false),1,1,[]));
    PRs{g} = cell2mat(reshape(arrayfun(@(f) mean(Ps(:,:,fi==f),3,'omitnan'),1:numel(uf),'UniformOutput',false),1,1,[]));
    FLY{g} = uf;
    fprintf('%-26s %d trials, %d flies\n', names{g}, N, numel(uf));
end

%% shape metrics per fly: amplitude (peak - trough) and FWHM
AMP = cell(1,3); WID = cell(1,3);
for g = 1:3
    PF = PRx{g}; nF = size(PF,3); AMP{g} = nan(nb,nF); WID{g} = nan(nb,nF);
    for f = 1:nF
        for k = 1:nb
            p = PF(:,k,f); if any(isnan(p)), continue, end
            AMP{g}(k,f) = max(p) - min(p);
            WID{g}(k,f) = fwhm_ring(p);
        end
    end
    fast = [1 2 nb-1 nb]; slow = [5 6 7];
    fprintf('%-26s centered bump amplitude (peak-trough, z): still %.2f, fast %.2f | FWHM (deg): still %.0f, fast %.0f\n', names{g}, ...
        mean(mean(AMP{g}(slow,:),1),'omitnan'), mean(mean(AMP{g}(fast,:),1),'omitnan'), ...
        mean(mean(WID{g}(slow,:),1),'omitnan'), mean(mean(WID{g}(fast,:),1),'omitnan'));
end

%% figure
cmap = [linspace(0,.85,5)' linspace(.3,.85,5)' ones(5,1); .3 .3 .3; ones(5,1) linspace(.85,.3,5)' linspace(.85,0,5)'];
th = rad2deg(alpha16);
figure(1); clf; set(gcf,'Color','w','Position',[50 50 1700 1050]); tl = tiledlayout(3,5,'TileSpacing','compact');
for g = 1:3
    M = mean(PRx{g},3,'omitnan'); Ms = mean(PRs{g},3,'omitnan'); nF = numel(FLY{g});
    nexttile; imagesc(1:nb, th, M); colorbar; xticks(1:nb); xticklabels(compose('%.1f',spc)); xtickangle(60)
    ylabel({names{g},'angle from bump centre (deg)'}); xlabel('signed rotational speed (rad/s)'); title(sprintf('cross-centered mean z (%d flies)',nF))
    nexttile; hold on
    for k = 1:nb, plot(th, M(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 1.5); end
    xlim([-180 157.5]); xlabel('angle from bump centre (deg)'); ylabel('mean z'); title('bump shape per speed bin (blue -vr, grey still, red +vr)')
    nexttile; hold on
    for k = 1:nb
        p = M(:,k); pn = (p - min(p))/(max(p) - min(p));
        plot(th, pn, '-', 'Color', cmap(k,:), 'LineWidth', 1.5)
    end
    xlim([-180 157.5]); ylim([-0.05 1.05]); xlabel('angle from bump centre (deg)'); ylabel('normalized (min 0, max 1)'); title('shape only (each profile min-max normalized)')
    nexttile; hold on
    errorbar(spc, mean(AMP{g},2,'omitnan'), std(AMP{g},0,2,'omitnan')./sqrt(sum(~isnan(AMP{g}),2)), 'o-', 'Color', c(g,:), 'LineWidth', 1.5)
    as = squeeze(max(PRs{g},[],1) - min(PRs{g},[],1));
    plot(spc, mean(as,2,'omitnan'), ':', 'Color', c(g,:), 'LineWidth', 1.2)
    xlabel('signed rotational speed (rad/s)'); ylabel('bump amplitude, peak - trough (z)'); legend({'cross-centered','self-centered (biased)'},'Box','off','Location','north')
    title('amplitude vs speed')
    nexttile
    errorbar(spc, mean(WID{g},2,'omitnan'), std(WID{g},0,2,'omitnan')./sqrt(sum(~isnan(WID{g}),2)), 'o-', 'Color', c(g,:), 'LineWidth', 1.5)
    xlabel('signed rotational speed (rad/s)'); ylabel('FWHM (deg)'); title('width vs speed'); ylim([0 360])
end
title(tl,{'bump shape vs rotational speed: each PB half centered on the OTHER half''s bump phase (fractional shift), halves averaged', ...
          'per-bin z of dF/F; trials -> flies -> mean; closed loop + dark pooled'})
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_by_speed.png'),'Resolution',150);

%% can a static nonlinearity map the LPsP curves onto the GRAB curves?
% Treat the fly-averaged cross-centered profiles as 16 angles x 11 speeds
% points. A static monotonic f needs GRAB(theta,v) = f(Ca(theta,v)) for all
% of them, i.e. one non-decreasing curve through the scatter. Best possible
% f = isotonic regression. Compared with a 2-input model where GRAB depends
% on calcium's global level (PB mean at that speed, Cg) and local deviation
% (Cl = Ca - Cg) with a gain interaction: a + b*Cg + (c + d*Cg)*Cl.
% Caveats: these are averaged curves (f applied per frame then averaged is
% not f of the average -- Jensen), sensors are in different flies, and each
% is z-scored in its own units (any monotonic rescaling is absorbed by f).
% The frame-level rank-correlation test in fig_grab_vs_lpsp_ca_nonlinearity.m
% doesn't have the averaging caveat.
Mc = mean(PRx{1},3,'omitnan'); Mg = mean(PRx{2},3,'omitnan');
x = Mc(:); y = Mg(:); [~,Vg] = ndgrid(alpha16, spc);
r2 = @(yy,yh) 1 - sum((yy-yh).^2)/sum((yy-mean(yy)).^2);
yiso = isofit(x,y);
Cg = repmat(mean(Mc,1),16,1); Cl = Mc - Cg;
Xg = [ones(size(y)) Cg(:) Cl(:) Cg(:).*Cl(:)]; yg = Xg*(Xg\y);
Xa = [ones(size(y)) Cg(:) Cl(:)];              ya = Xa*(Xa\y);
fprintf('\nGRAB curves from LPsP curves: best monotonic f R2 %.3f | additive global+local R2 %.3f | gain model R2 %.3f\n', r2(y,yiso), r2(y,ya), r2(y,yg));
beta = arrayfun(@(k) polyfit(Cl(:,k), Mg(:,k)-mean(Mg(:,k)), 1)*[1;0], 1:nb);
fprintf('within-speed slope dGRAB/dCa_local: %s ; corr with Ca global level %.2f\n', mat2str(beta,2), corr(beta', mean(Mc,1)'));
rng(1); nB = 500; R2b = nan(nB,2);
for bb = 1:nB
    mc = mean(PRx{1}(:,:,randi(size(PRx{1},3),size(PRx{1},3),1)),3,'omitnan');
    mg = mean(PRx{2}(:,:,randi(size(PRx{2},3),size(PRx{2},3),1)),3,'omitnan');
    cg = repmat(mean(mc,1),16,1); cl = mc - cg; Xb = [ones(numel(mg),1) cg(:) cl(:) cg(:).*cl(:)];
    R2b(bb,:) = [r2(mg(:),isofit(mc(:),mg(:))) r2(mg(:),Xb*(Xb\mg(:)))];
end
fprintf('fly bootstrap: monotonic f R2 %.2f [%.2f %.2f], gain model %.2f [%.2f %.2f], gain model better in %.0f%%\n', ...
    median(R2b(:,1)), prctile(R2b(:,1),[2.5 97.5]), median(R2b(:,2)), prctile(R2b(:,2),[2.5 97.5]), 100*mean(R2b(:,2) > R2b(:,1)));

figure(2); clf; set(gcf,'Color','w','Position',[50 50 1700 450]); tiledlayout(1,4,'TileSpacing','compact');
nexttile; hold on
scatter(x, y, 30, Vg(:), 'filled'); [xs,o] = sort(x); plot(xs, yiso(o), 'k-', 'LineWidth', 2)
colormap(gca,cmap); cb = colorbar; cb.Label.String = 'signed rotational speed bin (rad/s)';
xlabel('LPsP syt7f mean z'); ylabel('GRAB mean z'); title(sprintf('all angle x speed points + best monotonic f (R^2 = %.2f)', r2(y,yiso)))
nexttile; hold on
for k = 1:nb, plot(Mc(:,k), Mg(:,k), '-o', 'Color', cmap(k,:), 'MarkerSize', 3, 'LineWidth', 1.2); end
xlabel('LPsP z'); ylabel('GRAB z'); title('one strand per speed bin (loop over angle)')
nexttile; hold on
Yi = reshape(yiso,16,nb); Yg = reshape(yg,16,nb);
for k = [1 3 6 9 11]
    plot(th, Mg(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 2)
    plot(th, Yi(:,k), '--', 'Color', cmap(k,:), 'LineWidth', 1.2)
    plot(th, Yg(:,k), ':', 'Color', cmap(k,:), 'LineWidth', 1.8)
end
xlabel('angle from bump centre (deg)'); ylabel('z'); title({'GRAB (solid), best f(Ca) (dashed),', sprintf('gain model (dotted, R^2 = %.2f)', r2(y,yg))})
nexttile
plot(mean(Mc,1), beta, 'o', 'MarkerFaceColor', c(2,:), 'Color', c(2,:)); hold on
for k = 1:nb, text(mean(Mc(:,k)), beta(k), sprintf('  %.1f', spc(k)), 'FontSize', 8); end
xlabel('LPsP global level at that speed (mean z)'); ylabel('dGRAB / dCa_{local} (within-speed slope)')
title(sprintf('local gain grows with global calcium (r = %.2f)', corr(beta', mean(Mc,1)')))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_static_nl_test.png'),'Resolution',150);

%% normalization check: is the static-NL result a z-score artifact?
% Same bump centering (phases from z, cross-half), but the averaged VALUES
% are (1) per-bin z, (2) dF/F, (3) dF/F / one scalar per trial (PB-mean
% 5-95% range; same for every bin, so within-trial ordering = dF/F, but
% expression differences across flies are evened out).
% Result when written: the no-single-static-f conclusion holds for all three
% (best monotonic f reproduces only ~35-60% of GRAB's fast-turn bump; many
% order-violating point pairs). Normalization-dependent: in dF/F the LPsP
% bump DOES grow with speed (~1.8x), just less than its global level
% (~3.6x), so LPsP relative contrast (bump / mean) falls ~0.40 -> 0.19
% while GRAB's rises ~0.34 -> 0.85; and the simple 4-parameter gain model's
% advantage over the best f is large in z units but not in dF/F.
normN = {'z (per bin)','dF/F','dF/F / trial scalar'}; PRv = cell(3,2);
for g = 1:2
    N = numel(U{g}); Pv = nan(16,nb,N,3); fl = {U{g}.fly}';
    for n = 1:N
        D = U{g}(n).D; m = mean(D,1,'omitnan'); r = m - movmedian(m,21,'omitnan');
        bad = r > 6*1.4826*mad(r,1) | isnan(m); bad = conv(double(bad),ones(1,7),'same') > 0;
        Z = (D - mean(D(:,~bad),2)) ./ std(D(:,~bad),0,2);
        phL = angle(exp(1i*alpha16).'*Z(1:16,:)); phR = angle(exp(1i*alpha16).'*Z(17:32,:));
        dLR = angle(mean(exp(1i*(phL(~bad)-phR(~bad)))));
        V = {Z, D, D/(prctile(m(~bad),95) - prctile(m(~bad),5))};
        x_ = m; x_(bad) = mean(x_(~bad)); sp_ = abs(U{g}(n).vr); sp_(isnan(sp_)) = 0; lg = -9:27;
        xc = xcorr(x_-mean(x_), sp_-mean(sp_), lg(end), 'coeff'); [~,jl] = max(xc(lg+lg(end)+1));
        b = discretize(circshift(U{g}(n).vr, lg(jl)), spE);
        for v = 1:3
            C = (center_ring(V{v}(1:16,:), phR+dLR) + center_ring(V{v}(17:32,:), phL-dLR))/2;
            for k = 1:nb, sel = b==k & ~bad; if sum(sel) >= min_frames, Pv(:,k,n,v) = mean(C(:,sel),2); end, end
        end
    end
    [uf,~,fi] = unique(fl);
    for v = 1:3, PRv{v,g} = cell2mat(reshape(arrayfun(@(f) mean(Pv(:,:,fi==f,v),3,'omitnan'),1:numel(uf),'UniformOutput',false),1,1,[])); end
end
figure(3); clf; set(gcf,'Color','w','Position',[50 50 1600 900]); tl3 = tiledlayout(2,3,'TileSpacing','compact');
for v = 1:3
    mc = mean(PRv{v,1},3,'omitnan'); mg = mean(PRv{v,2},3,'omitnan'); yi = reshape(isofit(mc(:),mg(:)),16,nb);
    cg = repmat(mean(mc,1),16,1); cl = mc - cg; Xg = [ones(numel(mg),1) cg(:) cl(:) cg(:).*cl(:)];
    [I,J] = meshgrid(1:numel(mc)); viol = mc(I) < mc(J) & mg(I) - mg(J) > 0.25*range(mg(:));
    fprintf('%-20s best f R2 %.3f, gain model R2 %.3f, strongly order-violating pairs %d | fast-turn bump: GRAB %.3f vs best f(Ca) %.3f | bump/mean contrast still->fast: Ca %.2f->%.2f, GRAB %.2f->%.2f\n', ...
        normN{v}, r2(mg(:),yi(:)), r2(mg(:),Xg*(Xg\mg(:))), sum(viol(:)), mean(range(mg(:,[1 nb]))), mean(range(yi(:,[1 nb]))), ...
        range(mc(:,6))/mean(mc(:,6)), mean(range(mc(:,[1 nb]))./mean(mc(:,[1 nb]))), range(mg(:,6))/mean(mg(:,6)), mean(range(mg(:,[1 nb]))./mean(mg(:,[1 nb]))));
    nexttile(tl3,v); hold on
    for k = 1:nb, plot(mc(:,k), mg(:,k), '-o', 'Color', cmap(k,:), 'MarkerSize', 3, 'LineWidth', 1.4); end
    [xs,o] = sort(mc(:)); yy = yi(:); plot(xs, yy(o), 'k-', 'LineWidth', 1)
    xlabel('LPsP calcium'); ylabel('GRAB'); title(sprintf('%s: strands + best monotonic f', normN{v}))
    nexttile(tl3,v+3); hold on
    for k = [1 3 6 9 11], plot(th, mg(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 2); plot(th, yi(:,k), '--', 'Color', cmap(k,:), 'LineWidth', 1.3); end
    xlabel('angle from bump centre (deg)'); title('GRAB (solid) vs best f(Ca) (dashed)')
end
title(tl3,'same bump centering, three normalizations of the averaged values')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_normalization_check.png'),'Resolution',150);

%% power law per speed bin: GRAB = b * Ca^n (dF/F, both positive)
% Per speed bin: fit (b, n) on log GRAB vs log Ca across the 16 angles --
% n = how strongly the calcium bump is expanded WITHIN that speed. Shared:
% one (b, n) for all angle x speed points -- what a static power law would
% be. If a static power law held, per-bin n would be constant and equal the
% shared n. Fly-bootstrap CIs (resample flies of each sensor separately).
% Result when written: per-bin n rises from ~0.85 at rest to ~3-3.5 in fast
% turns, while the shared n is ~0.9 (outside the per-bin CIs for |vr| >~
% 0.9 rad/s): an exponent near 1 explains how GRAB scales ACROSS speeds,
% but the bump shape WITHIN a fast turn needs ~3.
vP = 2;   % dF/F (use 3 for dF/F / trial scalar -- same picture)
mc = mean(PRv{vP,1},3,'omitnan'); mg = mean(PRv{vP,2},3,'omitnan');
pfit = @(xc,yg) polyfit(log(xc(:)), log(yg(:)), 1);   % [n, log b]
nk = nan(1,nb); bk = nk;
for k = 1:nb, p_ = pfit(mc(:,k),mg(:,k)); nk(k) = p_(1); bk(k) = exp(p_(2)); end
ps = pfit(mc(:),mg(:)); ns = ps(1); bs = exp(ps(2));
rng(2); B = 500; nkb = nan(B,nb); nsb = nan(B,1);
for bb = 1:B
    a1 = mean(PRv{vP,1}(:,:,randi(size(PRv{vP,1},3),size(PRv{vP,1},3),1)),3,'omitnan');
    a2 = mean(PRv{vP,2}(:,:,randi(size(PRv{vP,2},3),size(PRv{vP,2},3),1)),3,'omitnan');
    if any(a1(:) <= 0) || any(a2(:) <= 0), continue, end
    for k = 1:nb, p_ = pfit(a1(:,k),a2(:,k)); nkb(bb,k) = p_(1); end
    p_ = pfit(a1(:),a2(:)); nsb(bb) = p_(1);
end
fprintf('\npower law GRAB = b*Ca^n (%s):\n  per-bin n %s\n  95%% CI low %s\n  95%% CI high %s\n  shared n %.2f [%.2f %.2f]\n', normN{vP}, ...
    mat2str(nk,2), mat2str(prctile(nkb,2.5),2), mat2str(prctile(nkb,97.5),2), ns, prctile(nsb,[2.5 97.5]));

figure(4); clf; set(gcf,'Color','w','Position',[50 50 1800 700]); tl4 = tiledlayout(2,7,'TileSpacing','compact');
for k = 1:nb
    nexttile(tl4,k); hold on
    plot(th, mg(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 2.5)
    plot(th, bk(k)*mc(:,k).^nk(k), 'k--', 'LineWidth', 1.3)
    plot(th, bs*mc(:,k).^ns, 'k:', 'LineWidth', 1.5)
    xlim([-180 157.5]); title(sprintf('%+.1f rad/s: n = %.1f', spc(k), nk(k)))
    if k == 1, ylabel('GRAB dF/F'); legend({'GRAB','b_k Ca^{n_k} (this bin)',sprintf('b Ca^{%.2f} (shared)',ns)},'Box','off','Location','north','FontSize',8); end
    if k > 7, xlabel('angle from bump (deg)'); end
end
nexttile(tl4,[1 3]); hold on
lo = prctile(nkb,2.5); hi = prctile(nkb,97.5);
fill([spc fliplr(spc)], [lo fliplr(hi)], c(2,:), 'FaceAlpha', .15, 'EdgeColor', 'none')
plot(spc, nk, 'o-', 'Color', c(2,:), 'MarkerFaceColor', c(2,:), 'LineWidth', 1.5)
yline(ns, 'k:', sprintf('shared n = %.2f', ns), 'LineWidth', 1.5); yline(1, ':', 'Color', [.6 .6 .6])
ylim([-1 8]); xlabel('signed rotational speed (rad/s)'); ylabel('exponent n'); title('within-bin exponent (fly bootstrap 95% CI) vs shared exponent')
title(tl4, sprintf('GRAB vs (LPsP calcium)^n, bump-centered %s profiles: n fit separately per speed bin (dashed) or shared across bins (dotted)', normN{vP}))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_powerlaw.png'),'Resolution',150);

%% straight rescale per speed bin -> moving threshold + global gain
% Within each speed bin GRAB is ~ a straight-line rescale of calcium,
% GRAB = a(v) + b(v)*Ca = b(v)*(Ca - th(v)). (Linearity alone is weak
% evidence -- Ca spans a narrow range within a bin -- but HOW a, b move with
% speed is informative.) Result when written (dF/F): th(v) ~ Ca_mean(v) -
% 0.15 (slope ~1.0, r ~1.00 across speeds) and b(v) rises with Ca_mean
% (r ~0.98). That compresses to 3 parameters:
%   GRAB(theta,v) ~ (b0 + b1*Ca_mean(v)) * (Ca(theta,v) - Ca_mean(v) + c)
% which matches the 22-parameter per-bin lines (bump-shape R2 ~0.94) and
% beats the best static monotonic f (~0.82) and a static power law (~0.44).
% The gain is nearly proportional to Ca_mean, so the DA bump ~ global
% calcium x local calcium deviation. Fits exclude the +3.5 rad/s bin (its
% averaged calcium bump is off-centre).
vL = 2; ok = 1:nb-1;
mc = mean(PRv{vL,1},3,'omitnan'); mg = mean(PRv{vL,2},3,'omitnan');
Cm = repmat(mean(mc,1),16,1); dev = mc - Cm; Y = mg(:,ok);
lin = nan(nb,2); for k = 1:nb, lin(k,:) = polyfit(mc(:,k), mg(:,k), 1); end
thv = -lin(:,2)'./lin(:,1)';
fC = @(q) sum((Y - (q(1) + q(2)*Cm(:,ok)).*(dev(:,ok) + q(3))).^2,'all');
qC = fminsearch(fC, [0.1 1 0.15], optimset('MaxFunEvals',2e4,'MaxIter',2e4));
yC = (qC(1) + qC(2)*Cm).*(dev + qC(3));
yA = nan(size(mg)); for k = 1:nb, yA(:,k) = polyval(lin(k,:), mc(:,k)); end
ps = polyfit(log(mc(:,ok)), log(Y), 1); yP = exp(ps(2))*mc.^ps(1);
yI = reshape(isofit(reshape(mc(:,ok),[],1), Y(:)), 16, numel(ok));
Wb = @(yh) 1 - sum(((Y-mean(Y,1)) - (yh-mean(yh,1))).^2,'all')/sum((Y-mean(Y,1)).^2,'all');
pth = polyfit(mean(mc(:,ok),1), thv(ok), 1);
fprintf('\nmoving threshold: th ~ %.2f + %.2f*Ca_mean (r = %.2f); gain b vs Ca_mean r = %.2f\n', pth(2), pth(1), corr(thv(ok)', mean(mc(:,ok),1)'), corr(lin(ok,1), mean(mc(:,ok),1)'));
fprintf('bump-shape R2: per-bin lines (22 p) %.3f | 3-param global x local %.3f [gain = %.3f + %.3f*Ca_mean, c = %.3f] | best static f %.3f | static power law %.3f\n', ...
    Wb(yA(:,ok)), Wb(yC(:,ok)), qC, Wb(yI), Wb(yP(:,ok)));

figure(5); clf; set(gcf,'Color','w','Position',[50 50 1800 700]); tl5 = tiledlayout(2,7,'TileSpacing','compact');
for k = 1:nb
    nexttile(tl5,k); hold on
    plot(th, mg(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 2.5)
    plot(th, yC(:,k), 'k--', 'LineWidth', 1.3)
    xlim([-180 157.5]); title(sprintf('%+.1f rad/s', spc(k)))
    if k == 1, ylabel('GRAB dF/F'); legend({'GRAB','3-param global x local'},'Box','off','Location','north','FontSize',8); end
    if k > 7, xlabel('angle from bump (deg)'); end
end
nexttile(tl5,[1 3]); hold on
camean = mean(mc,1);
yyaxis left
plot(camean(ok), thv(ok), 'o', 'MarkerFaceColor', c(2,:), 'Color', c(2,:)); plot(camean(ok), camean(ok), 'k:'); plot(camean(ok), polyval(pth, camean(ok)), '-', 'Color', c(2,:))
ylabel('threshold th(v): Ca where the line hits GRAB = 0')
yyaxis right
plot(camean(ok), lin(ok,1), 's', 'MarkerFaceColor', c(1,:), 'Color', c(1,:)); ylabel('gain b(v) (slope of GRAB vs Ca)')
xlabel('LPsP global calcium at that speed (PB-mean dF/F)')
title(sprintf('per-bin line GRAB = b(v)(Ca - th(v)): th = %.2f + %.2f Ca_{mean} (dotted = identity), gain rises with Ca_{mean}', pth(2), pth(1)))
title(tl5, sprintf('GRAB ~ (%.2f + %.2f Ca_{mean}) x (Ca - Ca_{mean} + %.2f): 3 parameters, bump-shape R^2 = %.2f (best static f %.2f, static power law %.2f)', ...
    qC, Wb(yC(:,ok)), Wb(yI), Wb(yP(:,ok))))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_global_x_local.png'),'Resolution',150);

%% predict LPsP calcium from EPG bump + rotational speed
% LPsP (dF/F) ~ G(|v|) + k * EPG bump (z, each speed bin's mean removed).
% EPG is in z because epg_7f's f is offset-subtracted (negative raw values),
% so its dF/F isn't meaningful; k absorbs the units. Caveat: averaged
% cross-centered bumps are ~97-100% first harmonic (cosine) for all three
% sensors at 16-bin resolution, so beyond amplitude the "shape" carries
% little information -- the bump-amplitude ratios below are the real test.
% Result when written: LPsP bump / EPG bump ~0.10-0.13 at all speeds (fixed,
% additive coupling; rest slightly lower), while GRAB bump / EPG bump rises
% ~4x; GRAB bump / (LPsP bump x LPsP global) is the most constant ratio,
% linking to the global x local GRAB rule. 4 parameters (G = g0 +
% g1*tanh(|v|/s), k) give LPsP total R2 ~0.97, bump R2 ~0.92.
Me = mean(PRx{3},3,'omitnan'); Edev = Me - mean(Me,1);
Ml = mean(PRv{2,1},3,'omitnan'); Mg2 = mean(PRv{2,2},3,'omitnan');
ok = 1:nb-1; Y = Ml(:,ok); S_ = repmat(abs(spc(ok)),16,1);
Gp = @(q) q(1) + q(2)*tanh(S_/max(q(3),0.05)); o = optimset('MaxFunEvals',4e4,'MaxIter',4e4);
q0 = fminsearch(@(q) sum((Y - Gp(q)).^2,'all'), [min(Y(:)) range(Y(:)) 1], o);
q1 = fminsearch(@(q) sum((Y - Gp(q(1:3)) - q(4)*Edev(:,ok)).^2,'all'), [q0 0.1], o);
q3 = fminsearch(@(q) sum((Y - Gp(q(1:3)).*(1 + q(4)*Edev(:,ok))).^2,'all'), [q0 0.1], o);
tot = @(yh) 1 - sum((Y-yh).^2,'all')/sum((Y-mean(Y,'all')).^2,'all');
bmp = @(yh) 1 - sum(((Y-mean(Y,1)) - (yh-mean(yh,1))).^2,'all')/sum((Y-mean(Y,1)).^2,'all');
yL = Gp(q1(1:3)) + q1(4)*Edev(:,ok); yM = Gp(q3(1:3)).*(1 + q3(4)*Edev(:,ok));
fprintf('\nLPsP from speed + EPG bump (bins excl +3.5): speed only total R2 %.3f | + k*EPG (additive, 4 p) total %.3f, bump %.3f | multiplicative (4 p) total %.3f, bump %.3f\n', ...
    tot(Gp(q0)), tot(yL), bmp(yL), tot(yM), bmp(yM));
fprintf('  G(|v|) = %.3f + %.3f*tanh(|v|/%.2f), k = %.3f\n', q1);
A1 = @(X) arrayfun(@(k) 2*abs(sum((X(:,k)-mean(X(:,k))).*exp(-1i*alpha16)))/16, 1:size(X,2));
Ae = A1(Me); Al = A1(Ml); Ag = A1(Mg2); Gl = mean(Ml,1);
fprintf('  bump-amplitude ratios by speed %s:\n    LPsP/EPG %s\n    GRAB/EPG %s\n    GRAB/(LPsP x LPsP global) %s\n', mat2str(spc,2), mat2str(Al./Ae,2), mat2str(Ag./Ae,2), mat2str(Ag./(Al.*Gl),2));

figure(6); clf; set(gcf,'Color','w','Position',[50 50 1800 700]); tl6 = tiledlayout(2,7,'TileSpacing','compact');
yLall = q1(1) + q1(2)*tanh(repmat(abs(spc),16,1)/q1(3)) + q1(4)*Edev;
for k = 1:nb
    nexttile(tl6,k); hold on
    plot(th, Ml(:,k), '-', 'Color', cmap(k,:), 'LineWidth', 2.5); plot(th, yLall(:,k), 'k--', 'LineWidth', 1.3)
    xlim([-180 157.5]); title(sprintf('%+.1f rad/s', spc(k)))
    if k == 1, ylabel('LPsP dF/F'); legend({'LPsP','G(|v|) + k EPG bump'},'Box','off','Location','south','FontSize',8); end
    if k > 7, xlabel('angle from bump (deg)'); end
end
nexttile(tl6,[1 3]); hold on
plot(spc, (Al./Ae)/mean(Al(ok)./Ae(ok)), 'o-', 'Color', c(1,:), 'MarkerFaceColor', c(1,:), 'LineWidth', 1.5)
plot(spc, (Ag./Ae)/mean(Ag(ok)./Ae(ok)), 's-', 'Color', c(2,:), 'MarkerFaceColor', c(2,:), 'LineWidth', 1.5)
plot(spc, (Ag./(Al.*Gl))/mean(Ag(ok)./(Al(ok).*Gl(ok))), 'd--', 'Color', [.4 .4 .4], 'LineWidth', 1.2)
yline(1,':'); xlabel('signed rotational speed (rad/s)'); ylabel('ratio / its mean (1 = constant across speeds)')
legend({'LPsP bump / EPG bump','GRAB bump / EPG bump','GRAB bump / (LPsP bump x LPsP global)'},'Box','off','Location','north')
title('bump-amplitude ratios across speed (+3.5 rad/s bin unreliable)')
title(tl6, sprintf('LPsP calcium ~ G(|v|) + k x EPG bump: G = %.2f + %.2f tanh(|v|/%.1f), k = %.2f; total R^2 %.2f, bump R^2 %.2f (multiplicative alternative: bump R^2 %.2f)', ...
    q1, tot(yL), bmp(yL), bmp(yM)))
fontsize(gcf,11,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_pb_centered_bump_lpsp_from_epg.png'),'Resolution',150);

%% functions
function yhat = isofit(x, y)
% isotonic (non-decreasing) least-squares fit of y on x, pool adjacent violators
[~, o] = sort(x(:)); ys = y(o); ys = ys(:);
v = ys; w = ones(size(v)); idx = num2cell((1:numel(v))'); k = 1;
while k < numel(v)
    if v(k) > v(k+1)
        v(k) = (v(k)*w(k) + v(k+1)*w(k+1))/(w(k)+w(k+1)); w(k) = w(k)+w(k+1); idx{k} = [idx{k}; idx{k+1}];
        v(k+1) = []; w(k+1) = []; idx(k+1) = []; k = max(k-1,1);
    else
        k = k+1;
    end
end
yh = zeros(size(ys)); for j = 1:numel(v), yh(idx{j}) = v(j); end
yhat = zeros(numel(y),1); yhat(o) = yh; yhat = reshape(yhat, size(y));
end

function C = center_ring(X, ph)
% circularly resample each column of X (16 bins on the ring, alpha -pi..pi)
% so that angle ph(t) moves to 0 rad (bin 9); linear interpolation
nb = size(X,1); step = 2*pi/nb; T = size(X,2);
nanph = isnan(ph(:)'); ph(nanph) = 0;          % dropped frames -> NaN output columns
a = linspace(-pi, pi - step, nb)';
pos = (a + ph(:)' + pi)/step;                  % fractional 0-based source index, nb x T
pos = mod(pos, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
cols = repmat(0:T-1, nb, 1)*nb;
C = (1-w).*X(lo+1+cols) + w.*X(hi+1+cols);
C(:,nanph) = nan;
end

function wdeg = fwhm_ring(p)
% full width at half (peak - trough) height of a centered 16-bin ring
% profile, by linear interpolation outward from the peak bin
nb = numel(p); p = p(:); half = min(p) + (max(p)-min(p))/2; [~,k0] = max(p);
step = 360/nb; w = 0;
for dirn = [-1 1]
    k = k0; dist = 0;
    for s = 1:nb/2
        kn = mod(k-1+dirn, nb)+1;
        if p(kn) < half
            dist = dist + (p(k)-half)/(p(k)-p(kn)); break
        end
        dist = dist + 1; k = kn;
    end
    w = w + dist;
end
wdeg = w*step;
end
