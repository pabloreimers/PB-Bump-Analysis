%% fig_rnai_th_lap_onset
% In LAPS, does the bump start turning before the fly? Laps from the window
% detector (fig_rnai_th_lapping_examples.m, which saves the episode tables
% data/lpsp_rnai_lap_episodes_{closed_loop,dark}.mat -- run it first, in both
% conditions). Each lap is matched to the turn-cache turn it contains (same
% direction, overlapping the lap window; the largest if several). All other
% turns >= min_turn_deg of the same flies are the non-lap reference, so lap and
% non-lap turns are measured identically.
%
%   fly onset  = the turn's onset (turn cache: smoothed turning speed crossing
%                0.3 rad/s), t = 0
%   bump onset = from the bump's peak velocity in the turn (onset - 0.5 s to
%                offset + 1 s), the last time going back that bump velocity (same
%                0.5 s smoothing as the fly's, signed in the turn direction) was
%                below the same 0.3 rad/s threshold. NEGATIVE = the bump started
%                turning first. Onsets earlier than the cached window (-1 s) are
%                censored at -1 s and counted.
% Stats: within each fly that has lap turns, median onset of lap vs non-lap
% turns (signrank across flies), plus turn-level fractions with the bump first.
% Needs the turn cache, the integrative-gain cache and the episode tables.
% Run from the repo root.

%% parameters
min_turn_deg = 30; thr_dps = rad2deg(0.3); gain_split = 0.85;
geno_set = 'th'; % 'th' or 'vglut'
if strcmp(geno_set,'vglut'), genos = {'lpsp>vglut','empty>vglut'}; main_name = 'LPsP > VGlut-RNAi'; else, genos = {'lpsp>th','empty>th'}; main_name = 'LPsP > TH-RNAi'; end
main_fld = strrep(genos{1},'>','_');
if exist('.data','dir'), data_dir = '.data'; else, data_dir = 'data'; end
IG = load(fullfile(data_dir,'lpsp_rnai_intgain_cache.mat'),'genotype','is_dark','fly_num','inc','g_trial');
TC = load(fullfile(data_dir,'lpsp_rnai_turns_cache.mat'),'turns','turn_t');
use_g = IG.inc & ~cellfun(@isempty,IG.g_trial) & ~IG.is_dark;
g_fly = nan(max(IG.fly_num),1);
for f = unique(IG.fly_num(use_g))', g_fly(f) = mean(vertcat(IG.g_trial{use_g & IG.fly_num==f}),'omitnan'); end
out_dir = fullfile('ugly_figures','rnai_turns','laps_stalls');
t = TC.turn_t(:); tt = [TC.turns.trial]'; on = [TC.turns.onset_s]'; du = [TC.turns.dur_s]'; dir_t = [TC.turns.dir]'; sz = [TC.turns.size_deg]';

%% per-turn bump onset (relative to fly onset)
n = numel(TC.turns); b_on = nan(n,1); cens = false(n,1);
for k = 1:n
    bv = double(TC.turns(k).bump_vel_trace); m = t >= -0.5 & t <= du(k) + 1;
    [~,ip] = max(bv.*m - 1e9*~m);
    j = find(bv(1:ip) < thr_dps,1,'last');
    if isempty(j), b_on(k) = t(1); cens(k) = true; else, b_on(k) = t(min(j+1,ip)); end
end

%% label lap turns, per condition
res = struct(); conds = {'closed_loop','dark'};
for ci = 1:2
    L = load(fullfile(data_dir,sprintf('lpsp_rnai_lap_episodes_%s.mat',conds{ci})),'E'); E = L.E;
    dk = ci == 2; ok_t = IG.is_dark(tt) == dk & sz >= min_turn_deg; if ~dk, ok_t = ok_t & IG.inc(tt); end
    is_lap = false(n,1); n_unmatched = 0;
    for e = 1:height(E)
        if ~ismember(E.geno{e},genos), continue, end
        c = find(ok_t & tt == E.trial(e) & on < E.t1(e) & on + du > E.t0(e) & dir_t == sign(E.dfly(e)));
        if isempty(c), n_unmatched = n_unmatched + 1; continue, end
        [~,j] = max(sz(c)); is_lap(c(j)) = true;
    end
    res(ci).is_lap = is_lap; res(ci).ok_t = ok_t; res(ci).n_unmatched = n_unmatched;
    fprintf('\n%s: %d laps matched to a turn, %d unmatched\n',strrep(conds{ci},'_',' '),sum(is_lap),n_unmatched);
    for gi = 1:numel(genos)
        g = genos{gi}; ing = strcmp(IG.genotype(tt),g); lapfl = unique(IG.fly_num(tt(is_lap & ing)));
        mL = is_lap & ing; mN = ok_t & ~is_lap & ing & ismember(IG.fly_num(tt),lapfl);
        fprintf('  %s: lap turns %d (bump first %.0f%%, median bump onset %+.2f s, %d censored) | non-lap turns of the same %d flies %d (bump first %.0f%%, median %+.2f s)\n', ...
            g,sum(mL),100*mean(b_on(mL) < 0),median(b_on(mL)),sum(cens(mL)),numel(lapfl),sum(mN),100*mean(b_on(mN) < 0),median(b_on(mN)));
        pf = nan(numel(lapfl),2);
        for q = 1:numel(lapfl)
            mf = IG.fly_num(tt) == lapfl(q); pf(q,:) = [median(b_on(mL & mf)), median(b_on(mN & mf))];
        end
        pf = pf(all(~isnan(pf),2),:);
        if size(pf,1) >= 3, p = signrank(pf(:,1),pf(:,2)); else, p = NaN; end
        fprintf('    per fly (n=%d): lap median %+.2f s vs non-lap %+.2f s, signrank p = %.3g\n',size(pf,1),median(pf(:,1)),median(pf(:,2)),p);
        res(ci).(strrep(g,'>','_')) = struct('mL',mL,'mN',mN,'pf',pf,'p',p);
    end
end

%% figure
figure('Color','w','Position',[40,40,1650,820]); clf
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
cl = [0.60,0.15,0.00]; cn = [0.6,0.6,0.6];
for ci = 1:2
    R = res(ci).(main_fld);
    % (a) onset-aligned mean velocity, lap vs non-lap turns (normalized per turn to its peak)
    nexttile(tl,(ci-1)*3+1); hold on
    for typ = 1:2
        if typ == 1, m = R.mL; c = cl; lbl = 'lap turns'; else, m = R.mN; c = cn; lbl = 'non-lap turns, same flies'; end
        B = cat(2,TC.turns(m).bump_vel_trace); F = cat(2,TC.turns(m).fly_vel_trace); B = double(B)./max(double(B),[],1); F = double(F)./max(double(F),[],1);
        mb = mean(B,2,'omitnan'); sb = std(B,0,2,'omitnan')/sqrt(size(B,2)); mf = mean(F,2,'omitnan');
        patch([t; flipud(t)],[mb+sb; flipud(mb-sb)],c,'EdgeColor','none','FaceAlpha',.25,'HandleVisibility','off')
        plot(t,mb,'Color',c,'LineWidth',2,'DisplayName',sprintf('bump, %s (n=%d)',lbl,sum(m)))
        plot(t,mf,'--','Color',c,'LineWidth',1.5,'DisplayName',sprintf('fly, %s',lbl))
    end
    xline(0,':k','HandleVisibility','off'); xlim([-1,2.5]); xlabel('time from fly turn onset (s)'); ylabel('angular velocity (per turn / its peak)')
    legend('Location','northeast','FontSize',7); title({sprintf(['(%s) %s, ' main_name ': onset-aligned velocity'],char('a'+(ci-1)*3),strrep(conds{ci},'_',' ')),'mean \pm SEM across turns'},'FontWeight','normal','FontSize',9)
    % (b) distribution of bump onset relative to fly onset
    nexttile(tl,(ci-1)*3+2); hold on
    edges = -1:0.1:1.5;
    histogram(b_on(R.mN),edges,'Normalization','probability','FaceColor',cn,'EdgeColor','none','FaceAlpha',.6)
    histogram(b_on(R.mL),edges,'Normalization','probability','FaceColor',cl,'EdgeColor','none','FaceAlpha',.6)
    xline(0,'k'); xlabel('bump onset - fly onset (s; < 0 = bump first, -1 = earlier/censored)'); ylabel('fraction of turns')
    legend(sprintf('non-lap (%d%% bump first)',round(100*mean(b_on(R.mN) < 0))),sprintf('lap (%d%% bump first)',round(100*mean(b_on(R.mL) < 0))),'Location','northeast')
    title({sprintf('(%s) when does the bump start?',char('b'+(ci-1)*3)),['turns of ' main_name ' flies that lap']},'FontWeight','normal','FontSize',9)
    % (c) per fly: lap vs non-lap median onset
    nexttile(tl,(ci-1)*3+3); hold on
    pf = R.pf;
    if ~isempty(pf)
        plot([1,2],pf','-','Color',[cl,.35]); scatter(ones(size(pf,1),1),pf(:,1),30,cl,'filled'); scatter(2*ones(size(pf,1),1),pf(:,2),30,cn,'filled')
    end
    yline(0,':k'); xticks([1,2]); xticklabels({'lap turns','non-lap turns'}); xlim([0.6,2.4]); ylabel('median bump onset - fly onset (s)')
    title({sprintf('(%s) per fly (n=%d)',char('c'+(ci-1)*3),size(pf,1)),sprintf('signrank p = %.2g',R.p)},'FontWeight','normal','FontSize',9)
end
title(tl,{['Do laps start with the bump? Lap turns vs other turns of the same ' main_name ' flies'], ...
    sprintf('bump onset = last crossing of %.0f%s/s (same threshold and smoothing as the fly''s turn onset) before the bump''s peak velocity',thr_dps,char(176))},'FontSize',11)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lap_onset.png',geno_set)),'Resolution',200)
exportgraphics(gcf,fullfile(out_dir,sprintf('fig_rnai_%s_lap_onset.pdf',geno_set)),'ContentType','vector')
