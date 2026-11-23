%% fig_summary_lpsp_grab_bump_model
% Summary figure: how bump shapes were measured, and the two-step model
%   LPsP calcium  ~ G(|vr|) + k * EPG bump                      (additive)
%   GRAB (DA)     ~ (g0 + g1 * LPsP global) * (1 + k2 * EPG bump)  (multiplicative)
% Standalone (recomputes everything). Data:
%   LPsP>syt7f, EPG>GRAB(DA2m): data/lpsp_cl_data_20240206.mat (closed loop +
%     dark; '*_EPG_syt7f_*' folders are a naming typo for GRAB trials)
%   EPG>GCaMP7f: data/epg_7f_20261007.mat -- its f is offset-subtracted
%     (raw values go negative), so per-bin dF/F (dividing by F0) is
%     meaningless there.
%
% NORMALIZATION (norm_mode below) -- all three sensors start from raw f and
% get the SAME treatment (a stopgap until all data are reprocessed
% identically):
%   'dF_scalar' (default): dF = f - F0 (per-bin 7th pct), divided by ONE
%       positive number per trial (5-95% range of the PB-mean dF). The
%       subtraction cancels EPG's background offset; the single divisor
%       avoids per-bin F0 division. Trade-off: brighter glomeruli weigh
%       more (no per-bin F0 normalization).
%   'z': per-bin z-score of f, all three sensors.
%   'dff_legacy': per-bin dF/F for LPsP/GRAB, z for EPG (previous version;
%       mixed units).
%
% Bump-shape pipeline (row 1), per trial:
%   1) per-bin z-score over time (artifact frames excluded) -> 32 bins, two
%      PB halves (1-16, 17-32), both spanning -pi..pi.
%   2) bump phase of each half = PVA angle of that half's z profile.
%   3) CROSS-CENTERING: each half is circularly shifted (fractional, linear
%      interpolation) so the OTHER half's phase (+/- the trial's mean
%      left-right offset) sits at 0 deg; halves averaged. Centering on a
%      half's own phase would align noise peaks and manufacture a bump.
%   4) signed rotational speed, shifted by the sensor's lag (xcorr of
%      PB-mean vs |vr|), binned; per-bin mean of the centered profile, in
%      norm_mode units (centering phases always come from z). Trials ->
%      flies -> mean of flies.
%   5) +/- speed bins folded to |vr| = 0, 0.3, 0.9, 1.6, 2.5 rad/s. (The
%      +/-3.5 bins are dropped: the +3.5 LPsP average is off-centre.)
% Result when written (data_src = 'reproc_nobg', norm_mode = 'dff', PMT offset
% removed, bg_static = false, snr_frac = 0.5; 9 EPG trials excluded for a
% meanImg/im.f scale mismatch; no bin needed flooring): LPsP ~ 0.36 + 1.01
% tanh(|vr|/2.0) + 0.27 x EPG bump (total R2 0.97, bump R2 0.95;
% multiplicative 0.86). GRAB ~ (0.08 + 0.06 |vr|) x (1 + 1.2 EPG bump)
% (bump R2 0.97; additive 0.85; LPsP-global version 0.97; best static
% f(LPsP Ca) 0.74; static power law 0.32). Bump / EPG bump, rest -> 2.5
% rad/s: LPsP 1.1x, GRAB 5.2x. bg_static = true gives the same fits but 8
% EPG bins with F0 <= 0.
% Result (data_src = 'reproc', norm_mode = 'dF_scalar',
% snr_frac = 0.5): LPsP ~ 0.28 + 0.75 tanh(|vr|/1.9) + 0.10 x EPG bump
% (total R2 0.97, bump R2 0.93; multiplicative alternative 0.89). GRAB ~
% (0.06 + 1.06 LPsP_global) x (1 + 0.52 EPG bump) (bump R2 0.97; additive
% 0.84; best static f(LPsP Ca) 0.81; static power law 0.36). Bump amplitude
% / EPG bump amplitude rest -> 2.5 rad/s: LPsP 1.4x, GRAB 5.2x. Dropping
% low-SNR bins doesn't bias the EPG bump (shape r = 0.999 vs no dropping).
% Earlier result (norm_mode = 'dff_legacy', legacy data): LPsP ~ 0.17 + 0.39 tanh(|vr|/1.7) + 0.05 x EPG bump
% (total R2 0.97, bump R2 0.91; multiplicative alternative 0.87). LPsP bump /
% EPG bump changes ~1.4x from rest to 2.5 rad/s, GRAB's ~4x. GRAB ~ (0.08 +
% 0.85 LPsP_global) x (0.18 + 0.06 EPG bump) (bump R2 0.98; additive 0.89;
% best static f(LPsP Ca) 0.83; best static power law 0.40).
% Note: earlier exploratory versions z-scored EPG dF/F built on negative
% baselines (sign-flipping some bins); here EPG is z-scored f directly.
% Caveats: three different sets of flies (EPG on another rig, CL gain 0.7);
% these are fly-averaged curves, not frame-by-frame predictions; each
% signal is centered on its own bump (assumes LPsP and EPG bumps in phase);
% at 16 bins the averaged bumps are ~cosines, so "shape" ~ amplitude.

%% load
if ~exist('data_src','var'), data_src = 'reproc'; end       % 'reproc' (lpsp_cl_reproc_bgsub, all sensors, one pipeline) | 'legacy'
if ~exist('snr_frac','var'), snr_frac = 0.5; end           % drop bins with SNR < snr_frac x trial median (0 = keep all)
if ~exist('bg_static','var'), bg_static = false; end      % reproc_nobg only: also subtract the tissue-background LEVEL (bg_pb - offset), not just its fluctuations
if ~exist('grab_driver','var'), grab_driver = 'lpsp_global'; end  % GRAB gain term: 'lpsp_global' (measured LPsP PB-mean) | 'speed' (g0 + g1 tanh(|vr|/s))
if exist('.data','dir') && isfile(fullfile('.data','lpsp_cl_data_20240206.mat')), data_dir = '.data'; else, data_dir = 'data'; end

alpha16 = linspace(-pi, pi - 2*pi/16, 16)'; th = rad2deg(alpha16);
spE = [-inf -3 -2 -1.25 -0.5 -0.15 0.15 0.5 1.25 2 3 inf]; nb = numel(spE)-1;
sv  = [0 0.3 0.9 1.6 2.5];                                  % folded |vr| levels (nominal)
fold = @(M) [M(:,6,:), (M(:,5,:)+M(:,7,:))/2, (M(:,4,:)+M(:,8,:))/2, (M(:,3,:)+M(:,9,:))/2, (M(:,2,:)+M(:,10,:))/2];
min_frames = 20;
if ~exist('norm_mode','var'), norm_mode = 'dF_scalar'; end   % 'dF_scalar' | 'z' | 'dff_legacy' (set before running to override)
names = {'EPG > GCaMP7f','LPsP > syt7f','EPG > GRAB(DA2m)'};      % order used throughout: 1 EPG, 2 LPsP, 3 GRAB
switch norm_mode
    case 'dF_scalar',  units = repmat({'\DeltaF / trial scale'},1,3);
    case 'z',          units = repmat({'z'},1,3);
    case 'dff_legacy', units = {'z','dF/F','dF/F'};
    case 'dff',        units = repmat({'dF/F'},1,3);
end
scol  = [0.30 0.30 0.30; 0.85 0.33 0.10; 0.00 0.45 0.74];
vcol  = [0.15 0.15 0.15; 0.35 0.35 0.70; 0.30 0.55 0.90; 0.95 0.55 0.25; 0.80 0.15 0.10];   % |vr| levels
flyid = @(m) regexprep(regexp(m,'\d{8}-\d+_.*$','match','once'),'^(\d{8})-\d+_.*_(\d+)$','$1_$2');

%% unify trials: V = f (32 x T), bad = excluded frames, vr on imaging clock, fr, fly
% U{g}(n): V, bad, vr, xb, fr, fly. Frame rate and timestamps are per trial.
U = {[],[],[]};
switch data_src
    case 'reproc'
        % One pipeline for all sensors (background-subtracted; leak flashes
        % corrected -> only badVol frames are excluded; no ad-hoc flash
        % detector). lpsp_cl ft.xb are evenly spaced placeholders spanning the
        % FicTrac record; trials with a volume clock (epg_7f) start imaging
        % ~1.6 ms before FicTrac, so that offset is applied to the
        % placeholders for consistency (negligible in practice). The ~0.5-0.8 s
        % lag of LPsP/GRAB bumps behind the cue is NOT a clock problem (cue vs
        % FicTrac agree on both rigs; same-frame co-imaging shows GRAB trailing
        % EPG calcium by ~0.5 s = sensor kinetics); per-trial lags absorb it.
        A = load(fullfile(data_dir,'lpsp_cl_reproc_bgsub_20261117.mat'),'all_data'); A = A.all_data(:);
        info = [A.info]; ind = {info.indicator}; dsn = {info.dataset};
        isV = strcmp(dsn,'epg_7f'); x0 = median(arrayfun(@(s) s.ft.xb(1) - s.ft.xf(1), A(isV)));
        gidx = {find(strcmp(ind,'EPG>GCaMP7f')), find(strcmp(ind,'LPsP>syt7f')), find(strcmp(ind,'EPG>GRAB(DA2m)'))};
        for g = 1:3
            for i = gidx{g}
                s = A(i); xf = s.ft.xf(:); xb = s.ft.xb(:)';
                if ~strcmp(s.info.dataset,'epg_7f'), xb = xb + x0; end
                if strcmp(s.info.dataset,'epg_7f'), fl = sprintf('epg7f_fly%02d', s.info.fly_num); else, fl = flyid(s.meta); end
                U{g} = [U{g}; struct('V',s.im.f,'bad',logical(s.info.badVol(:)'),'vr',interp1(xf,smoothdata(s.ft.r_speed(:),'gaussian',30),xb,'linear','extrap'), ...
                    'xb',xb,'fr',s.info.fps,'fly',fl)]; %#ok<AGROW>
            end
        end
    case 'reproc_nobg'
        % Same trials/labels/timing as 'reproc', but im.f WITHOUT the
        % pipeline's background subtraction (which over-subtracted: ~200
        % counts on the lpsp_cl rig vs a measured background region of ~55,
        % ~1400 vs ~890 on the EPG rig -> negative/floored baselines). Here:
        %   f = f_nobg - PMT offset - (background fluctuation)
        % PMT offset per trial = darkest-pixel level (0.1th pct of meanImg):
        % ~3-5 counts on the lpsp_cl rig, ~320-640 on most EPG trials (hard
        % floor in the image histogram). Background fluctuation = smoothed
        % bg_pb minus its median: removes leak flashes and drift without
        % lowering the baseline. bg_static = true also removes the tissue
        % background level (bg_pb - offset). Trials whose meanImg/bg_pb are
        % on a different scale from im.f (dark level >= the PB's median
        % baseline; EPG flies 12-15, May-Jun 2023) are excluded.
        A = load(fullfile(data_dir,'lpsp_cl_reproc_nobg_20261117.mat'),'all_data'); A = A.all_data(:);
        info = [A.info]; ind = {info.indicator}; dsn = {info.dataset};
        isV = strcmp(dsn,'epg_7f'); x0 = median(arrayfun(@(s) s.ft.xb(1) - s.ft.xf(1), A(isV)));
        gidx = {find(strcmp(ind,'EPG>GCaMP7f')), find(strcmp(ind,'LPsP>syt7f')), find(strcmp(ind,'EPG>GRAB(DA2m)'))};
        n_excl = 0;
        for g = 1:3
            for i = gidx{g}
                s = A(i); xf = s.ft.xf(:); xb = s.ft.xb(:)'; bad = logical(s.info.badVol(:)');
                off = prctile(double(s.flash.meanImg(:)), 0.1);
                if off >= median(prctile(s.im.f(:,~bad),7,2))
                    n_excl = n_excl + 1; fprintf('  excluded (meanImg/bg on a different scale from im.f): %s\n', s.info.movie_src); continue
                end
                bgs = smoothdata(double(s.flash.bg_pb(:)'), 'gaussian', s.flash.params.smoothWindow);
                if bg_static, V = s.im.f - bgs; else, V = s.im.f - off - (bgs - median(bgs(~bad))); end
                if ~strcmp(s.info.dataset,'epg_7f'), xb = xb + x0; end
                if strcmp(s.info.dataset,'epg_7f'), fl = sprintf('epg7f_fly%02d', s.info.fly_num); else, fl = flyid(s.meta); end
                U{g} = [U{g}; struct('V',V,'bad',bad,'vr',interp1(xf,smoothdata(s.ft.r_speed(:),'gaussian',30),xb,'linear','extrap'), ...
                    'xb',xb,'fr',s.info.fps,'fly',fl)]; %#ok<AGROW>
            end
        end
        fprintf('reproc_nobg: %d trials excluded for scale mismatch; bg_static = %d\n', n_excl, bg_static);
    case 'legacy'
        d = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
        d = d(~cellfun(@isempty,{d.meta})); meta = {d.meta};
        isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); isLpsp = contains(meta,'LPsP_syt7f');
        E = load(fullfile(data_dir,'epg_7f_20261007.mat'),'all_data'); E = E.all_data(:);
        for k = 1:numel(E)
            s = E(k); f = smoothdata(s.im.f,2,'gaussian',5);
            parts = strsplit(strrep(s.meta.flyDir,'\','/'),'/'); parts = parts(~cellfun(@isempty,parts));
            if startsWith(parts{end},'fly'), fl = [parts{end-1} '_' strrep(parts{end},' ','')]; else, fl = parts{end}; end
            U{1} = [U{1}; struct('V',f,'bad',legacy_flash(f) | logical(s.im.badVol(:)'),'vr',interp1(s.ft.xf, smoothdata(s.ft.r_speed_raw,'gaussian',30), s.ft.xb), ...
                'xb',s.ft.xb(:)','fr',9.1,'fly',fl)]; %#ok<AGROW>
        end
        for g = 2:3
            for i = find((g==2 & isLpsp) | (g==3 & isGrab))
                s = d(i); T = size(s.im.d,2); xf = s.ft.xf; xb = linspace(xf(1),xf(end),T);
                U{g} = [U{g}; struct('V',s.im.f,'bad',legacy_flash(s.im.f),'vr',interp1(xf,smoothdata(s.ft.r_speed,'gaussian',30),xb), ...
                    'xb',xb,'fr',9.1,'fly',flyid(meta{i}))]; %#ok<AGROW>
            end
        end
end

%% per trial: drop low-SNR bins, cross-center, bin by signed speed; fly averages
% Low-SNR bins: per-bin SNR = slow (signal) / fast (noise) variance of f,
% from autocovariance at lags 4-6 frames; bins below snr_frac x the trial's
% median SNR are set to NaN (mostly the dim PB midline and ends). Bump
% phases are then computed from the remaining bins, and the two centered
% halves are averaged where available (each half covers the other's gaps).
PR = cell(1,3); FLYS = cell(1,3); EX = struct(); NDROP = cell(1,3);
for g = 1:3
    N = numel(U{g}); P = nan(16,nb,N); ndrop = zeros(N,1);
    for n = 1:N
        u = U{g}(n); V = u.V; fr_n = u.fr;
        bad = u.bad | all(isnan(V),1);
        if snr_frac > 0
            x = V(:,~bad); x = x - mean(x,2);
            ac = mean(cell2mat(arrayfun(@(L) mean(x(:,1:end-L).*x(:,1+L:end),2), 4:6, 'UniformOutput', false)),2);
            snr = ac ./ max(var(x,0,2) - ac, eps);
            drop = snr < snr_frac*median(snr); V(drop,:) = nan; ndrop(n) = sum(drop);
        end
        Z = (V - mean(V(:,~bad),2,'omitnan')) ./ std(V(:,~bad),0,2,'omitnan');
        Zp = Z; Zp(isnan(Zp)) = 0;
        phL = angle(exp(1i*alpha16).'*Zp(1:16,:)); phR = angle(exp(1i*alpha16).'*Zp(17:32,:));
        dLR = angle(mean(exp(1i*(phL(~bad)-phR(~bad)))));
        F0 = prctile(V(:,~bad),7,2); dF = V - F0;
        switch norm_mode
            case 'dF_scalar'
                mdF = mean(dF(:,~bad),1,'omitnan'); X = dF / (prctile(mdF,95) - prctile(mdF,5));
            case 'z'
                X = Z;
            case 'dff_legacy'
                if g == 1, X = Z; else, X = dF ./ F0; end
            case 'dff'      % plain per-bin dF/F, no flooring; any non-positive baseline -> bin dropped
                nonpos = F0 <= 0; if any(nonpos), fprintf('  %s trial %d: %d bins with F0 <= 0 dropped\n', names{g}, n, sum(nonpos)); end
                F0(nonpos) = nan; X = dF ./ F0;
        end
        C = mean(cat(3, center_ring(X(1:16,:), phR + dLR), center_ring(X(17:32,:), phL - dLR)), 3, 'omitnan');
        m = mean(V,1,'omitnan'); x_ = m; x_(bad) = mean(x_(~bad)); sp_ = abs(u.vr); sp_(isnan(sp_)) = 0;
        lg = -round(fr_n):round(3*fr_n);
        xc = xcorr(x_-mean(x_), sp_-mean(sp_), lg(end), 'coeff'); [~,jl] = max(xc(lg+lg(end)+1));
        b = discretize(circshift(u.vr, lg(jl)), spE);
        for k = 1:nb
            sel = b==k & ~bad;
            if sum(sel) >= min_frames
                Pk = mean(C(:,sel),2,'omitnan'); Pk(sum(~isnan(C(:,sel)),2) < min_frames) = nan; P(:,k,n) = Pk;
            end
        end
        if g == 1   % keep per-trial pieces to pick a method example later
            EX(n).Z = Z; EX(n).phL = phL; EX(n).phR = phR; EX(n).dLR = dLR; EX(n).vr = u.vr; EX(n).bad = bad; EX(n).xb = u.xb; EX(n).fr = fr_n;
            EX(n).lock = mean(abs(exp(1i*alpha16).'*Zp(1:16,~bad)))/mean(sum(abs(Zp(1:16,~bad)),1));
        end
    end
    [uf,~,fi] = unique({U{g}.fly}');
    PR{g} = cell2mat(reshape(arrayfun(@(f) mean(P(:,:,fi==f),3,'omitnan'), 1:numel(uf), 'UniformOutput', false),1,1,[]));
    FLYS{g} = uf; NDROP{g} = ndrop;
    fprintf('%-18s %d trials, %d flies; low-SNR bins dropped per trial: median %d (max %d)\n', names{g}, N, numel(uf), median(ndrop), max(ndrop));
end
PF = cellfun(@(x) fold(x), PR, 'UniformOutput', false);            % 16 x 5 x flies, folded
M  = cellfun(@(x) mean(x,3,'omitnan'), PF, 'UniformOutput', false); % 16 x 5 fly means

%% models on the fly-averaged folded curves
Ef = M{1}; Lf = M{2}; Gf = M{3};
Ed = Ef - mean(Ef,1);                         % EPG bump (z, each speed's mean removed)
Lg = repmat(mean(Lf,1),16,1);                 % LPsP global level per speed
S_ = repmat(sv,16,1);
o  = optimset('MaxFunEvals',4e4,'MaxIter',4e4);
R2 = @(Y,yh) 1 - sum((Y-yh).^2,'all')/sum((Y-mean(Y,'all')).^2,'all');
B2 = @(Y,yh) 1 - sum(((Y-mean(Y,1))-(yh-mean(yh,1))).^2,'all')/sum((Y-mean(Y,1)).^2,'all');
Gp = @(q,S) q(1) + q(2)*tanh(S/max(q(3),0.05));
% LPsP = G(|v|) + k*EPG bump
qL = fminsearch(@(q) sum((Lf - Gp(q(1:3),S_) - q(4)*Ed).^2,'all'), [0.1 0.4 1 0.1], o);
yL = Gp(qL(1:3),S_) + qL(4)*Ed;
qLm = fminsearch(@(q) sum((Lf - Gp(q(1:3),S_).*(1 + q(4)*Ed)).^2,'all'), [qL(1:3) 0.1], o);
yLm = Gp(qLm(1:3),S_).*(1 + qLm(4)*Ed);
% GRAB = (g0 + g1*LPsP global) * (1 + k2*EPG bump)   [3 params; a free constant in the
% second factor would be redundant with the scale of the first]
g0_ = [min(Gf(:)) range(Gf(:))/max(range(Lg(:)),eps) 0.3];
qGL = fminsearch(@(q) sum((Gf - (q(1)+q(2)*Lg).*(1+q(3)*Ed)).^2,'all'), g0_, o);
yGL = (qGL(1)+qGL(2)*Lg).*(1+qGL(3)*Ed);
qGaL = fminsearch(@(q) sum((Gf - (q(1)+q(2)*Lg+q(3)*Ed)).^2,'all'), [0 0.2 0.02], o);
yGaL = qGaL(1) + qGaL(2)*Lg + qGaL(3)*Ed;
% same model with the gain driven by rotational speed directly (no LPsP data):
% GRAB = (g0 + g1 |vr|) * (1 + k2*EPG bump). A saturating tanh(|vr|/s) gain
% runs to s -> infinity (i.e. linear) over 0-2.5 rad/s, unlike LPsP's
% global calcium, which saturates (s ~ 1.9 rad/s) -- so linear is used.
qGS = fminsearch(@(q) sum((Gf - (q(1)+q(2)*S_).*(1+q(3)*Ed)).^2,'all'), [min(Gf(:)) range(Gf(:))/2.5 0.3], o);
yGS = (qGS(1)+qGS(2)*S_).*(1+qGS(3)*Ed);
qGaS = fminsearch(@(q) sum((Gf - (q(1)+q(2)*S_+q(3)*Ed)).^2,'all'), [min(Gf(:)) range(Gf(:))/2.5 0.1], o);
yGaS = qGaS(1) + qGaS(2)*S_ + qGaS(3)*Ed;
fprintf('GRAB gain from speed: (%.3f + %.3f |v|)(1 + %.3f EPG): total R2 %.3f, bump R2 %.3f | additive alt bump %.3f\n', qGS, R2(Gf,yGS), B2(Gf,yGS), B2(Gf,yGaS));
fprintf('  speed-driven vs LPsP-global-driven GRAB predictions: corr %.3f; corr(LPsP global, |v|) across speeds %.3f\n', corr(yGS(:),yGL(:)), corr(mean(Lf,1)', sv'));
switch grab_driver
    case 'lpsp_global', qG = qGL; yG = yGL; qGa = qGaL; yGa = yGaL;
    case 'speed',       qG = qGS; yG = yGS; qGa = qGaS; yGa = yGaS;
end
if all(Lf(:) > 0) && all(Gf(:) > 0)
    ps = polyfit(log(Lf(:)), log(Gf(:)), 1); yP = exp(ps(2))*Lf.^ps(1);   % best static power law of LPsP Ca
else
    ps = [nan nan]; yP = nan(size(Gf));                                   % undefined for non-positive units (z)
end
yI = reshape(isofit(Lf(:), Gf(:)), size(Gf));                           % best static monotonic f (any shape)
fprintf('LPsP = G(|v|) + k*EPG: total R2 %.3f, bump R2 %.3f [G = %.3f + %.3f tanh(|v|/%.2f), k = %.3f] | multiplicative alt bump R2 %.3f\n', R2(Lf,yL), B2(Lf,yL), qL, B2(Lf,yLm));
fprintf('GRAB = (g0+g1*LPsPglob)(1+k2*EPG): total R2 %.3f, bump R2 %.3f [%s] | additive alt bump %.3f | static power law bump %.3f (n=%.2f) | best static f bump %.3f\n', ...
    R2(Gf,yGL), B2(Gf,yGL), mat2str(qGL,3), B2(Gf,yGaL), B2(Gf,yP), ps(1), B2(Gf,yI));

%% bump-amplitude ratios with fly bootstrap
A1 = @(X) squeeze(2*abs(sum((X - mean(X,1)).*exp(-1i*alpha16),1))/16);   % first-harmonic amplitude per column
rng(1); B = 500; RL = nan(B,5); RG = RL;
for bb = 1:B
    e_ = mean(PF{1}(:,:,randi(size(PF{1},3),size(PF{1},3),1)),3,'omitnan');
    l_ = mean(PF{2}(:,:,randi(size(PF{2},3),size(PF{2},3),1)),3,'omitnan');
    g_ = mean(PF{3}(:,:,randi(size(PF{3},3),size(PF{3},3),1)),3,'omitnan');
    RL(bb,:) = A1(l_)./A1(e_); RG(bb,:) = A1(g_)./A1(e_);
end
rl = A1(Lf)./A1(Ef); rg = A1(Gf)./A1(Ef);
fprintf('bump amplitude / EPG bump amplitude by |vr| %s: LPsP %s, GRAB %s\n', mat2str(sv), mat2str(rl',2), mat2str(rg',3));

%% method example: EPG trial with a clear bump, a turning stretch
lk = [EX.lock]; lk(NDROP{1}(:)' > 0) = -inf;            % example trial: no dropped bins (clean illustration)
[~,ex] = max(lk); X0 = EX(ex);
spd = abs(X0.vr); spd(isnan(spd)) = 0; w = round(20*X0.fr);
score = movsum(spd, w) .* (movsum(~X0.bad, w) == w);
[~,cmid] = max(score); a = max(1, cmid - floor(w/2)); bnd = min(numel(spd), a + w - 1); win = a:bnd;
% example frame: both halves carry a clear bump that cross-centers to within one bin of 0 deg
amp_ = abs(exp(1i*alpha16).'*X0.Z(1:16,win)) + abs(exp(1i*alpha16).'*X0.Z(17:32,win));
c1w = center_ring(X0.Z(1:16,win), X0.phR(win) + X0.dLR); c2w = center_ring(X0.Z(17:32,win), X0.phL(win) - X0.dLR);
[~,k1] = max(c1w,[],1); [~,k2] = max(c2w,[],1);
t1 = th(k1); t2 = th(k2);
good_ = abs(t1(:))' <= 22.5 & abs(t2(:))' <= 22.5 & abs(X0.phL(win)) < pi/2 & abs(X0.phR(win)) < pi/2;
amp_(~good_) = -inf; [~,jf] = max(amp_); tf = win(jf);

%% figure
fig = figure(1); clf; set(fig,'Color','w','Position',[40 40 1750 1150]);
tl = tiledlayout(6,12,'TileSpacing','compact','Padding','compact');
toBin = @(p) mod(p+pi,2*pi)/(2*pi)*16 + 0.5;
lbl = @(ax,s) text(ax,-0.12,1.08,s,'Units','normalized','FontWeight','bold','FontSize',14);

% A: example EPG trial, two PB halves, cross-estimated bump phases; speed below
axA = nexttile(tl,1,[1 4]); t = X0.xb(win) - X0.xb(win(1));
hA = imagesc(axA, t, 1:32, X0.Z(:,win)); set(hA,'AlphaData',~isnan(X0.Z(:,win))); colormap(axA, parula); clim(axA, prctile(X0.Z(:,win),[5 99.5],'all')); hold(axA,'on')
yline(axA, 16.5, 'w-', 'LineWidth', 1.5)
pL = toBin(X0.phR(win) + X0.dLR); pR = toBin(X0.phL(win) - X0.dLR) + 16;
pL(abs(diff([pL pL(end)])) > 8) = nan; pR(abs(diff([pR pR(end)])) > 8) = nan;
plot(axA, t, pL, 'w-', 'LineWidth', 1.3); plot(axA, t, pR, 'w-', 'LineWidth', 1.3)
xline(axA, X0.xb(tf) - X0.xb(win(1)), 'r-', 'LineWidth', 1.2)
set(axA,'XTick',[],'YTick',[8 24],'YTickLabel',{'half 1','half 2'}); ylabel(axA,'PB bin')
title(axA,'per-bin z (EPG example); white = bump phase taken from the OTHER half','FontWeight','normal')
lbl(axA,'A')
axA2 = nexttile(tl,13,[1 4]); plot(axA2, t, X0.vr(win), 'k-', 'LineWidth', 1); hold(axA2,'on'); yline(axA2,0,':')
xline(axA2, X0.xb(tf) - X0.xb(win(1)), 'r-', 'LineWidth', 1.2); xlim(axA2,[t(1) t(end)])
xlabel(axA2,'time (s)'); ylabel(axA2,'v_r (rad/s)'); title(axA2,'rotational speed (lag-shifted, then binned)','FontWeight','normal')

% B: one frame -> two halves -> cross-centered -> averaged
% B top: the two halves of one frame are the same bump, offset by a fixed
% per-trial delta (circular mean over frames of half-1 phase - half-2 phase)
gC = [0.2 0.6 0.2]; pC = [0.6 0.2 0.6];
axB = nexttile(tl,5,[1 2]); hold(axB,'on')
plot(axB, th, X0.Z(1:16,tf), '-', 'Color', gC, 'LineWidth', 1.6)
plot(axB, th, X0.Z(17:32,tf), '-', 'Color', pC, 'LineWidth', 1.6)
yTop = max(X0.Z(:,tf)) + 0.6; pkL = rad2deg(X0.phL(tf)); pkR = rad2deg(X0.phR(tf));
plot(axB, pkL, yTop, 'v', 'Color', gC, 'MarkerFaceColor', gC, 'MarkerSize', 6)
plot(axB, pkR, yTop, 'v', 'Color', pC, 'MarkerFaceColor', pC, 'MarkerSize', 6)
plot(axB, [pkR pkL], (yTop + 0.55)*[1 1], 'k-', 'LineWidth', 1.1)
if pkL < pkR, hd = '<'; else, hd = '>'; end
plot(axB, pkL, yTop + 0.55, hd, 'Color', 'k', 'MarkerFaceColor', 'k', 'MarkerSize', 5)
text(axB, (pkL + pkR)/2, yTop + 1.4, sprintf('\\delta = %.0f%s', rad2deg(X0.dLR), char(176)), 'HorizontalAlignment','center','FontSize',9)
xlim(axB,[-180 157.5]); ylim(axB,[min(X0.Z(:,tf)) - 0.3, yTop + 2.0]); set(axB,'XTick',[-90 0 90]); xlabel(axB,'position in its half (deg)'); ylabel(axB,'z')
title(axB,{'one frame (red line in A): {\color[rgb]{0.2,0.6,0.2}half 1}, {\color[rgb]{0.6,0.2,0.6}half 2}', ...
    'same bump, offset by a fixed \delta per trial'},'FontWeight','normal'); lbl(axB,'B')
% B bottom: each half centred on the OTHER half's phase, corrected by delta
axB2 = nexttile(tl,17,[1 2]); hold(axB2,'on')
c1 = center_ring(X0.Z(1:16,tf), X0.phR(tf) + X0.dLR); c2 = center_ring(X0.Z(17:32,tf), X0.phL(tf) - X0.dLR);
plot(axB2, th, c1, '-', 'Color', gC); plot(axB2, th, c2, '-', 'Color', pC); plot(axB2, th, (c1+c2)/2, 'k-', 'LineWidth', 2)
xline(axB2, 0, ':'); xlim(axB2,[-180 157.5]); xlabel(axB2,'angle from bump (deg)'); ylabel(axB2,'z')
title(axB2,{'{\color[rgb]{0.2,0.6,0.2}half 1} shifted by (half-2 phase + \delta),', ...
    '{\color[rgb]{0.6,0.2,0.6}half 2} by (half-1 phase - \delta); black = average'},'FontWeight','normal')

% C: fly-averaged centered bumps per |vr|, three sensors
for g = 1:3
    axC = nexttile(tl, 7 + 2*(g-1), [2 2]); hold(axC,'on')
    for j = 1:5, plot(axC, th, M{g}(:,j), '-', 'Color', vcol(j,:), 'LineWidth', 1.8); end
    xlim(axC,[-180 157.5]); xlabel(axC,'angle from bump (deg)'); ylabel(axC, units{g})
    title(axC, sprintf('%s (%d flies)', names{g}, numel(FLYS{g})), 'Color', scol(g,:))
    if g == 1, lbl(axC,'C'); end
end

% D: LPsP global vs |vr|
axD = nexttile(tl,25,[2 4]); hold(axD,'on')
gl_f = squeeze(mean(PF{2},1));                         % 5 x flies
for f = 1:size(gl_f,2), plot(axD, sv, gl_f(:,f), '-', 'Color', [scol(2,:) 0.25]); end
plot(axD, sv, mean(Lf,1), 'o', 'Color', scol(2,:), 'MarkerFaceColor', scol(2,:), 'MarkerSize', 7)
ss = linspace(0,2.6,100); plot(axD, ss, Gp(qL(1:3),ss), 'k-', 'LineWidth', 1.5)
xlabel(axD,'|rotational speed| (rad/s)'); ylabel(axD,['LPsP global (PB-mean, ' units{2} ')'])
title(axD, sprintf('global: G(|v_r|) = %.2f + %.2f tanh(|v_r| / %.1f)', qL(1:3)), 'FontWeight','normal')
legend(axD, {'flies'}, 'Box','off','Location','southeast'); lbl(axD,'D')

% E: LPsP data vs G + k*EPG bump
axE = nexttile(tl,29,[2 4]); hold(axE,'on')
for j = 1:5
    plot(axE, th, Lf(:,j), '-', 'Color', vcol(j,:), 'LineWidth', 2.2)
    plot(axE, th, yL(:,j), '--', 'Color', vcol(j,:)*0.7, 'LineWidth', 1.2)
end
xlim(axE,[-180 157.5]); xlabel(axE,'angle from bump (deg)'); ylabel(axE,['LPsP (' units{2} ')'])
title(axE, {sprintf('LPsP \\approx G(|v_r|) + %.2f \\times EPG bump', qL(4)), ...
    sprintf('solid data, dashed model: R^2 %.2f, bump R^2 %.2f (\\times instead of +: %.2f)', R2(Lf,yL), B2(Lf,yL), B2(Lf,yLm))}, 'FontWeight','normal')
lbl(axE,'E')

% F: bump amplitude relative to EPG bump amplitude
axF = nexttile(tl,33,[2 4]); hold(axF,'on')
nl = mean(rl); ng = mean(rg);
fill(axF, [sv fliplr(sv)], [prctile(RL,2.5) fliplr(prctile(RL,97.5))]/nl, scol(2,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none')
fill(axF, [sv fliplr(sv)], [prctile(RG,2.5) fliplr(prctile(RG,97.5))]/ng, scol(3,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none')
h1 = plot(axF, sv, rl/nl, 'o-', 'Color', scol(2,:), 'MarkerFaceColor', scol(2,:), 'LineWidth', 1.8);
h2 = plot(axF, sv, rg/ng, 's-', 'Color', scol(3,:), 'MarkerFaceColor', scol(3,:), 'LineWidth', 1.8);
yline(axF, 1, ':'); xlabel(axF,'|rotational speed| (rad/s)'); ylabel(axF,'bump amp. / EPG bump amp. (/ mean)')
legend(axF, [h1 h2], {sprintf('LPsP Ca: changes %.1fx (rest -> %.1f rad/s)', rl(end)/rl(1), sv(end)), ...
    sprintf('GRAB: changes %.1fx', rg(end)/rg(1))}, 'Box','off','Location','northwest')
title(axF,'bump amplitude relative to EPG bump (fly bootstrap 95% CI)','FontWeight','normal'); lbl(axF,'F')

% G: GRAB vs LPsP -- one strand per speed, no single static function
axG = nexttile(tl,49,[2 4]); hold(axG,'on')
hv = gobjects(1,5);
for j = 1:5, hv(j) = plot(axG, Lf(:,j), Gf(:,j), '-o', 'Color', vcol(j,:), 'MarkerSize', 3, 'MarkerFaceColor', vcol(j,:), 'LineWidth', 1.3); end
xx = linspace(min(Lf(:)), max(Lf(:)), 100); if ~isnan(ps(1)), plot(axG, xx, exp(ps(2))*xx.^ps(1), 'k:', 'LineWidth', 1.5); end
legend(axG, hv, compose('|v_r| %.1f rad/s', sv), 'Box','off','Location','northwest','FontSize',9)
xlabel(axG,['LPsP calcium (' units{2} ')']); ylabel(axG,['GRAB (' units{3} ')'])
title(axG, {'GRAB vs LPsP at each angle: one steep strand per speed', ...
    sprintf('no single f(Ca): best power law (dotted) bump R^2 %.2f, best monotonic f %.2f', B2(Gf,yP), B2(Gf,yI))}, 'FontWeight','normal')
lbl(axG,'G')

% H: GRAB data vs (g0 + g1*LPsP global)(c + k*EPG bump)
axH = nexttile(tl,53,[2 4]); hold(axH,'on')
for j = 1:5
    plot(axH, th, Gf(:,j), '-', 'Color', vcol(j,:), 'LineWidth', 2.2)
    plot(axH, th, yG(:,j), '--', 'Color', vcol(j,:)*0.7, 'LineWidth', 1.2)
end
xlim(axH,[-180 157.5]); xlabel(axH,'angle from bump (deg)'); ylabel(axH,['GRAB (' units{3} ')'])
if strcmp(grab_driver,'speed')
    hEq = sprintf('GRAB \\approx (%.3g + %.3g |v_r|) \\times (1 + %.3g EPG bump)', qG);
else
    hEq = sprintf('GRAB \\approx (%.3g + %.3g LPsP_{global}) \\times (1 + %.3g EPG bump)', qG);
end
title(axH, {hEq, ...
    sprintf('solid data, dashed model: R^2 %.2f, bump R^2 %.2f (+ instead of \\times: %.2f)', R2(Gf,yG), B2(Gf,yG), B2(Gf,yGa))}, 'FontWeight','normal')
lbl(axH,'H')

% I: schematic -- inputs on top (speed, EPG bump), outputs below (LPsP Ca, GRAB)
axI = nexttile(tl,57,[2 4]); hold(axI,'on'); axis(axI,'off'); W = 15.2;   % ~ tile aspect, so equal data aspect fills it
daspect(axI,[1 1 1]); xlim(axI,[0 W]); ylim(axI,[0 10])
box_ = @(x,y,wd,ht,str,col) [rectangle(axI,'Position',[x y wd ht],'Curvature',0.25,'EdgeColor',col,'LineWidth',1.6), ...
    text(axI, x+wd/2, y+ht/2, str, 'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',10,'Color',col)];
arrow_ = @(p1,p2,col) arrow_draw(axI,p1,p2,col);
bw = 4.2; bh = 1.6; xL = 1.0; xR = W - 1.0 - bw; yT = 8.0; yB = 1.9;
box_(xL,yT,bw,bh,{'rotational','speed'},[0 0 0]);
box_(xR,yT,bw,bh,{'EPG bump','(heading)'},scol(1,:));
box_(xL,yB,bw,bh,{'LPsP syt7f','(Pre Ca)'},scol(2,:));
box_(xR,yB,bw,bh,{'EPG GRAB','(DA Release)'},scol(3,:));
cL = xL + bw/2; cR = xR + bw/2;
% operator circles where each output's inputs combine: + for calcium (panel E), x for GRAB (panel H)
r0 = 0.5; yO = yB + bh + 1.0; oL = [cL yO]; oR = [cR yO];
op_ = @(c,str,col) [rectangle(axI,'Position',[c(1)-r0 c(2)-r0 2*r0 2*r0],'Curvature',[1 1],'EdgeColor',col,'LineWidth',1.6,'FaceColor','w'), ...
    text(axI, c(1), c(2), str, 'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',15,'FontWeight','bold','Color',col)];
toCirc = @(p1,c) c - r0*(c-p1)/norm(c-p1);                      % arrow end on the circle's edge
lab_ = @(p1,p2,t,str,col,side,dy) text(axI, p1(1)+t*(p2(1)-p1(1)) + side*0.25, p1(2)+t*(p2(2)-p1(2)) + dy, str, ...
    'HorizontalAlignment', ifelse_h(side), 'FontSize', 9, 'Color', col);
% calcium: G(|vr|) + 0.27 * EPG bump
pS = [cL yT]; arrow_(pS, toCirc(pS,oL), [0 0 0]); lab_(pS, oL, 0.45, 'G(|v_r|)', [0 0 0], -1, 0);
pE = [xR+0.5 yT]; arrow_(pE, toCirc(pE,oL), scol(1,:)); text(axI, oL(1)+r0+0.25, oL(2)-0.35, sprintf('%.2f\\cdotEPG bump', qL(4)), 'HorizontalAlignment','left','FontSize',9,'Color',scol(1,:));
arrow_([cL yO-r0], [cL yB+bh], scol(2,:)); op_(oL, '+', scol(2,:));
% GRAB: (g0 + g1 |vr|) * (1 + k2 * EPG bump)
pE2 = [cR yT]; arrow_(pE2, toCirc(pE2,oR), scol(1,:)); lab_(pE2, oR, 0.45, {'1+k_2\cdot','EPG bump'}, scol(1,:), 1, 0);
if strcmp(grab_driver,'speed')
    pS2 = [xL+bw-0.5 yT]; arrow_(pS2, toCirc(pS2,oR), [0 0 0]); text(axI, oR(1)-r0-0.25, oR(2)-0.35, 'g_0+g_1|v_r|', 'HorizontalAlignment','right','FontSize',9,'Color',[0 0 0]);
    eqG = sprintf('GRAB = (%.2f + %.2f|v_r|)(1 + %.2f EPG)', qG(1), qG(2), qG(3));
else
    pC = [xL+bw yB+bh*0.75]; arrow_(pC, toCirc(pC,oR), scol(2,:)); lab_(pC, oR, 0.5, 'g_0 + g_1 Ca_{glob}', scol(2,:), -1, 0.3);
    eqG = sprintf('GRAB = (%.2f + %.2f Ca_{glob})(1 + %.2f EPG)', qG(1), qG(2), qG(3));
end
arrow_([cR yO-r0], [cR yB+bh], scol(3,:)); op_(oR, '\times', scol(3,:));
text(axI, cL, yB-0.55, sprintf('Ca = G(|v_r|) + %.2f EPG', qL(4)), 'HorizontalAlignment','center','FontSize',8.5,'Color',scol(2,:))
text(axI, cR, yB-0.55, eqG, 'HorizontalAlignment','center','FontSize',8.5,'Color',scol(3,:))
text(axI, cL, yB-1.3, sprintf('bump R^2 %.2f (\\times instead: %.2f)', B2(Lf,yL), B2(Lf,yLm)), 'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4])
text(axI, cR, yB-1.3, sprintf('bump R^2 %.2f (+ instead: %.2f)', B2(Gf,yG), B2(Gf,yGa)), 'HorizontalAlignment','center','FontSize',8,'Color',[0.4 0.4 0.4])
title(axI,{'summary: calcium ADDS its inputs (E),','dopamine MULTIPLIES them (H)'},'FontWeight','normal'); lbl(axI,'I')

title(tl, {'LPsP presynaptic calcium vs. dopamine release across the PB', ...
    sprintf('data: %s; units: %s (norm\\_mode = %s); bins with SNR < %.1f x trial median dropped; fly-averaged, bump-centered; different flies per sensor; |v_r| folded', strrep(data_src,'_','\_'), units{2}, strrep(norm_mode,'_','\_'), snr_frac)}, 'FontSize', 13)
fontsize(fig, 11, 'pixels')
exportgraphics(fig, fullfile('ugly_figures','exports',['fig_summary_lpsp_grab_bump_model_' data_src repmat('_bgstatic',1,bg_static & strcmp(data_src,'reproc_nobg')) '_' norm_mode '_grabFrom_' grab_driver '.png']), 'Resolution', 200);
exportgraphics(fig, fullfile('ugly_figures','exports',['fig_summary_lpsp_grab_bump_model_' data_src repmat('_bgstatic',1,bg_static & strcmp(data_src,'reproc_nobg')) '_' norm_mode '_grabFrom_' grab_driver '.pdf']), 'ContentType', 'vector');

%% functions
function h = ifelse_h(side)
% text alignment for a label placed left (-1) or right (+1) of an arrow
if side < 0, h = 'right'; else, h = 'left'; end
end

function arrow_draw(ax, p1, p2, col)
% straight arrow from p1 to p2 (data units, equal aspect axes): line + filled head
d = p2 - p1; L = norm(d); u = d/L; nrm = [-u(2) u(1)]; hl = 0.45; hw = 0.22;
plot(ax, [p1(1) p2(1)-u(1)*hl], [p1(2) p2(2)-u(2)*hl], '-', 'Color', col, 'LineWidth', 1.4)
base = p2 - u*hl;
patch(ax, [p2(1) base(1)+nrm(1)*hw base(1)-nrm(1)*hw], [p2(2) base(2)+nrm(2)*hw base(2)-nrm(2)*hw], col, 'EdgeColor', col)
end

function bad = legacy_flash(f)
% legacy-data artifact flashes: fast residual of the PB-mean > 6 robust SDs, dilated +/-3 frames
m = mean(f,1,'omitnan'); r = m - movmedian(m,21,'omitnan');
bad = r > 6*1.4826*mad(r,1) | isnan(m); bad = conv(double(bad),ones(1,7),'same') > 0;
end

function C = center_ring(X, ph)
% circularly resample each column of X (16-bin ring, alpha -pi..pi) so that
% angle ph(t) moves to 0 rad; linear interpolation; NaN phase -> NaN column
nb = size(X,1); step = 2*pi/nb; T = size(X,2);
nanph = isnan(ph(:)'); ph(nanph) = 0;
a = linspace(-pi, pi - step, nb)';
pos = mod((a + ph(:)' + pi)/step, nb); lo = floor(pos); w = pos - lo; hi = mod(lo+1, nb);
cols = repmat(0:T-1, nb, 1)*nb;
C = (1-w).*X(lo+1+cols) + w.*X(hi+1+cols);
C(:,nanph) = nan;
end

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
