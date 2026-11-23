%% fig_grab_stalls_and_sustained_turns
% Q: are there moments where the fly rotates but the EPG bump is stuck
% ("stalls"), and during those is GRAB(DA2m) more transient than LPsP
% calcium?
%
% Which data can answer it:
%   - LPsP>syt7f / EPG>GRAB (lpsp_cl_data_20240206.mat): no EPG bump channel.
%     LPsP calcium has no usable bump phase (its own PVA is pinned by the
%     static compartment profile), so stalls can't be defined there at all.
%     GRAB's own phase is too unreliable when DA is low.
%   - EPG jRGECO1a + EPG GRAB(DA2m) co-imaging (epg_coimage_20260120.mat):
%     the only set with an EPG bump and GRAB in the same fly, same frames.
%     Section 1 finds stalls from the jRGECO1a bump with the core criteria
%     of fig_rnai_vglut_stall_examples.m and shows them with both channels.
%     (Many of these trials have fewer frames than ft.xb -- timestamps use
%     Figure_1_claude.m's get_xb fallback, otherwise cue and bump misalign.)
%   - No dataset has LPsP calcium and an EPG bump in the same fly
%     (lpsp_joint_20260126.mat is EPG calcium in the RNAi sets).
%
% Result when written: candidate stalls exist (32 in 15 well-tracking
% co-imaged closed-loop trials) but none is convincing by eye -- the
% jRGECO1a bump is faint/diffuse in most trials and several "stuck" phases
% coincide with GRAB blobs that DO move with the cue (phase-estimate
% failures, not stuck bumps). Section 2 therefore asks the global half of
% the question without a bump estimate: during sustained turns, does the
% PB-mean signal decay faster for GRAB than for LPsP calcium? It doesn't
% (both stay up while the fly keeps turning).

%% load
if exist('.data','dir') && isfile(fullfile('.data','epg_coimage_20260120.mat')), data_dir = '.data'; else, data_dir = 'data'; end
C = load(fullfile(data_dir,'epg_coimage_20260120.mat'),'all_data'); C = C.all_data(:);
d = load(fullfile(data_dir,'lpsp_cl_data_20240206.mat'),'all_data'); d = d.all_data(:);
d = d(~cellfun(@isempty,{d.meta}));
meta   = {d.meta};
isGrab = contains(meta,'GRAB','IgnoreCase',true) | contains(meta,'EPG_syt7f'); % 'EPG_syt7f' folders are a naming typo for GRAB(DA2m)
isLpsp = contains(meta,'LPsP_syt7f');
fr_lpsp = 9.1;

%% 1) stalls from the co-imaged EPG jRGECO1a bump
P = struct('win_s',[5 20],'step_s',1,'cue_min_deg',180,'bump_frac',0.25,'ratio_max',0.25,'gain_range',[0.4 1.3], ...
           'follow_max',0.4,'rho_pct',25,'jump_max_deg',90,'bin_cue_min_deg',5,'min_turn_frac',0.2);
lock_min = 0.4;   % trial must track the cue overall (own-phase lock R, best lag in [-2,+1] s)
isCL = arrayfun(@(s) contains(string(s.ft.pattern),"bar"), C);
CO = struct('mu',{},'rho',{},'cue',{},'head',{},'zr',{},'zg',{},'bad',{},'sp',{},'fr',{},'id',{},'R',{});
for i = find(isCL(:))'
    s = C(i); n = size(s.grab.f,2); xf = s.ft.xf(:);
    if isfield(s.ft,'xb') && numel(s.ft.xb) == n, xb = s.ft.xb(:)'; else, xb = linspace(xf(1),xf(end),n); end
    frr = (n-1)/(xb(end)-xb(1));
    bad = flash_frames(s.grab.f,8) | flash_frames(s.geco.f,8); bad = conv(double(bad),ones(1,5),'same') > 0;
    u = struct('mu',s.geco.mu(:)','rho',s.geco.rho(:)','cue',angle(interp1(xf,exp(1i*s.ft.cue(:)),xb)), ...
        'head',interp1(xf,cumtrapz(xf,s.ft.r_speed(:)),xb),'zr',s.geco.z,'zg',s.grab.z,'bad',bad, ...
        'sp',interp1(xf,abs(s.ft.r_speed(:)),xb),'fr',frr,'id',i,'R',nan);
    u.R = cue_lock(u.mu,u.cue,frr);
    CO(end+1) = u; %#ok<SAGROW>
end
good = [CO.R] >= lock_min;
EPc = [];
for k = find(good)
    ep = find_stalls(CO(k).mu,CO(k).rho,CO(k).cue,CO(k).head,CO(k).fr,P);
    if ~isempty(ep), EPc = [EPc; repmat(k,size(ep,1),1) ep]; end %#ok<AGROW>
end
dur_min = sum(arrayfun(@(u) numel(u.mu)/u.fr, CO(good)))/60;
fprintf('co-imaged closed loop: %d trials, %d track the cue (R >= %.1f, %.0f min); candidate stalls %d (%.2f / 10 min)\n', ...
    numel(CO), sum(good), lock_min, dur_min, size(EPc,1), size(EPc,1)/dur_min*10);

% gallery: largest-gap candidates, jRGECO1a (top) and GRAB (bottom), z-scored per bin
EPc = sortrows(EPc,-4); nshow = min(16,size(EPc,1)); chn = {'EPG jRGECO1a','EPG GRAB(DA2m)'};
toBin = @(p) mod(p+pi,2*pi)/(2*pi)*16 + 0.5;
figure(1); clf; set(gcf,'Color','w','Position',[50 50 1800 2000]);
tl = tiledlayout(2*ceil(nshow/4),4,'TileSpacing','tight','Padding','compact');
for r = 1:nshow
    u = CO(EPc(r,1)); pd = round(6*u.fr); a = max(1,EPc(r,2)-pd); b = min(numel(u.mu),EPc(r,3)+pd); t = (a:b)/u.fr;
    row = floor((r-1)/4); col = mod(r-1,4)+1;
    for ch = 1:2
        nexttile(tl,(2*row+ch-1)*4+col);
        if ch == 1, Z = u.zr(:,a:b); else, Z = u.zg(:,a:b); end
        Z(:,u.bad(a:b)) = nan; imagesc(t,1:32,Z); caxis(prctile(Z(:),[5 99])); hold on
        mb = toBin(u.mu(a:b)); cb = toBin(-u.cue(a:b));
        mb(abs(diff([mb mb(end)])) > 8) = nan; cb(abs(diff([cb cb(end)])) > 8) = nan;
        if ch == 1, plot(t,mb,'w','LineWidth',1); plot(t,mb+16,'w','LineWidth',1); end
        plot(t,cb,'r','LineWidth',.8); plot(t,cb+16,'r','LineWidth',.8);
        xline(EPc(r,2)/u.fr,'y','LineWidth',1.5); xline(EPc(r,3)/u.fr,'y','LineWidth',1.5); set(gca,'XTick',[],'YTick',[])
        if ch == 1, title(sprintf('#%d %s, gap %.0f%s, |vr| %.1f', r, trial_tag(C(u.id).meta), EPc(r,4), char(176), ...
                mean(u.sp(EPc(r,2):EPc(r,3)))),'FontSize',8,'Interpreter','none'); end
        if col == 1, ylabel(chn{ch}); end
    end
end
title(tl,{'candidate stalls from the EPG jRGECO1a bump (white = jRGECO1a phase, red = -cue, yellow = stall window)', ...
          'z-scored per bin; 6 s padding; none convincing by eye -- see header'},'FontSize',11)
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_coimage_stall_candidates.png'),'Resolution',100);

%% 2) sustained turns: PB-mean time course, LPsP calcium vs GRAB (lpsp_cl, closed loop + dark)
names = {'LPsP > syt7f','EPG > GRAB(DA2m)'}; c = [0.85 0.33 0.10; 0 0.45 0.74];
turn_thr = 0.5; turn_min_s = 3; quiet_s = 1;
prepad = round(3*fr_lpsp); postpad = round(8*fr_lpsp);
TT = {[],[]}; SPt = TT;
grpL = {find(isLpsp), find(isGrab)};
for g = 1:2
    for i = grpL{g}
        s = d(i); D = s.im.d; T = size(D,2); xf = s.ft.xf; xb = linspace(xf(1),xf(end),T);
        m = mean(D,1); r = m - movmedian(m,21); bad = r > 6*1.4826*mad(r,1); bad = conv(double(bad),ones(1,7),'same') > 0; m(bad) = nan;
        sp = interp1(xf,abs(smoothdata(s.ft.r_speed,'gaussian',30)),xb);
        on = sp > turn_thr; dd = diff([0 on 0]); st = find(dd==1); en = find(dd==-1)-1;
        lo = prctile(m,5); hi = prctile(m,95);
        for e = 1:numel(st)
            if en(e)-st(e)+1 < round(turn_min_s*fr_lpsp) || st(e)-prepad < 1 || st(e)+postpad > T, continue, end
            if e > 1 && st(e)-en(e-1) < round(quiet_s*fr_lpsp), continue, end
            w = st(e)-prepad:st(e)+postpad;
            TT{g}(end+1,:) = (m(w)-lo)/(hi-lo); SPt{g}(end+1,:) = sp(w);
        end
    end
end
tt = (-prepad:postpad)/fr_lpsp;
wpk = prepad + (1:round(1.5*fr_lpsp)); wl = prepad + (round((turn_min_s-1)*fr_lpsp):round(turn_min_s*fr_lpsp)); wb = 1:round(fr_lpsp);
rng(0);
for g = 1:2
    X = TT{g}; n = size(X,1);
    r0 = @(I) (mean(mean(X(I,wl),1,'omitnan')) - mean(mean(X(I,wb),1,'omitnan'))) / ...
              (max(mean(X(I,wpk),1,'omitnan')) - mean(mean(X(I,wb),1,'omitnan')));
    bs = arrayfun(@(~) r0(randi(n,n,1)),1:1000);
    fprintf('%-18s turns %3d: PB-mean late (%g-%g s) / early peak (0-1.5 s) = %.2f [%.2f %.2f]; speed late/early %.2f\n', names{g}, n, ...
        turn_min_s-1, turn_min_s, r0(1:n), prctile(bs,[2.5 97.5]), mean(mean(SPt{g}(:,wl)))/mean(mean(SPt{g}(:,wpk))));
end
figure(2); clf; set(gcf,'Color','w','Position',[100 100 1000 380]);
subplot(1,2,1); hold on
for g = 1:2
    mu = mean(TT{g},'omitnan'); se = std(TT{g},'omitnan')/sqrt(size(TT{g},1));
    fill([tt fliplr(tt)],[mu+se fliplr(mu-se)],c(g,:),'FaceAlpha',.2,'EdgeColor','none'); plot(tt,mu,'Color',c(g,:),'LineWidth',2)
end
xline(0); xline(turn_min_s,'--'); xlabel('time from turn onset (s)'); ylabel('PB-mean dF/F (trial 5-95% range)')
title(sprintf('turns: |vr| > %.1f rad/s for \\geq %g s', turn_thr, turn_min_s))
subplot(1,2,2); hold on
for g = 1:2, plot(tt,mean(SPt{g},'omitnan'),'Color',c(g,:),'LineWidth',2); end
xline(0); xline(turn_min_s,'--'); xlabel('time from turn onset (s)'); ylabel('|rotational speed| (rad/s)')
legend(arrayfun(@(g) sprintf('%s (n=%d)',names{g},size(TT{g},1)),1:2,'UniformOutput',false),'Box','off')
fontsize(gcf,12,'pixels')
exportgraphics(gcf, fullfile('ugly_figures','exports','fig_grab_sustained_turns.png'),'Resolution',150);

%% functions
function bad = flash_frames(f, k)
% Figure_1_claude.m's detect_flash_frames: bright AND spatially uniform frames
fm = mean(f,1); cv = std(f,0,1)./fm;
bad = fm > median(fm) + k*mad(fm,1) & cv <= median(cv);
end

function R = cue_lock(mu, cue, fr)
T = numel(mu); R = 0;
for t = -round(2*fr):round(fr)
    idx = 1:T; idx = idx(idx+t >= 1 & idx+t <= T);
    for sg = [-1 1], R = max(R, abs(mean(exp(1i*(mu(idx) + sg*cue(idx+t))),'omitnan'))); end
end
end

function s = trial_tag(meta)
tok = regexp(meta,'\d{8}-\d+','match','once'); if isempty(tok), tok = meta; end
s = tok;
end

function ep = find_stalls(mu, rho, cue, fly_head, fr, P)
% core criteria of fig_rnai_vglut_stall_examples.m on the imaging clock:
% cue turns >= cue_min_deg in the window, bump moves <= bump_frac of it and
% stays near its start, cue/fly rotation agree, solid bump (median rho >=
% the trial's rho_pct percentile, no big phase jumps), stuck throughout the
% turn (1 s bins). Overlapping windows merged; returns [start end gap_deg].
mu_u = unwrap(mu(:)'); cu = unwrap(cue(:)'); fh = unwrap(fly_head(:)');
T = numel(mu_u); n1 = max(1,round(fr)); st = max(1,round(P.step_s*fr));
rho_thr = prctile(rho, P.rho_pct); cand = [];
for w_s = P.win_s
    w = round(w_s*fr);
    for a = 1:st:T-w
        b = a+w;
        dc = rad2deg(cu(b)-cu(a)); db = rad2deg(mu_u(b)-mu_u(a)); dfh = rad2deg(fh(b)-fh(a));
        if abs(dc) < P.cue_min_deg || abs(db) > P.bump_frac*abs(dc), continue, end
        g = abs(dc)/max(abs(dfh),eps); if g < P.gain_range(1) || g > P.gain_range(2), continue, end
        if max(abs(rad2deg(mu_u(a:b)-mu_u(a)))) > P.follow_max*abs(dc), continue, end
        r = rho(a:b); if median(r) < rho_thr || any(isnan(r)), continue, end
        if any(rad2deg(abs(angle(exp(1i*diff(mu(a:b)))))) > P.jump_max_deg), continue, end
        e = a:n1:b; dC = rad2deg(diff(cu(e))); dB = rad2deg(diff(mu_u(e)));
        turning = abs(dC) >= P.bin_cue_min_deg & sign(dC) == sign(dc);
        if sum(turning) < max(1,ceil(P.min_turn_frac*numel(dC))) || median(dB(turning)./dC(turning)) > P.ratio_max, continue, end
        cand(end+1,:) = [a b abs(dc-db)]; %#ok<AGROW>
    end
end
ep = [];
if isempty(cand), return, end
cand = sortrows(cand,1); cur = cand(1,:); grpEnd = cand(1,2);
for k = 2:size(cand,1)
    if cand(k,1) <= grpEnd
        grpEnd = max(grpEnd,cand(k,2)); if cand(k,3) > cur(3), cur = cand(k,:); end
    else
        ep(end+1,:) = cur; cur = cand(k,:); grpEnd = cand(k,2); %#ok<AGROW>
    end
end
ep(end+1,:) = cur;
end
