%% fig_ionto_celltypes
% Mean effect of dopamine iontophoresis on each imaged cell type, pooled
% across every dopamine_ionto dataset in .data, plus the effect of LPsP
% stimulation (lpsp>P2X2 + ATP iontophoresis, releasing DA + glutamate)
% on EPG, computed as lpsp>P2X2 minus a stim-matched empty>P2X2 control.
%
% Per-stim logic (same as lpsp_p2x2_claude.m): ft.stims rising edge = stim
% sent; the red (dye) channel's 0-3 s median must clear 3 baseline SDs
% (-5..0 s) for the stim to count as an ejection; ejection side = PB half
% with the larger evoked red response. Response = green dF/F (mean of 32
% glomeruli), 0-3 s minus -5..0 s. Pooled one value per fly, then mean +/-
% SEM across flies. Stims < 9.5 s after the previous onset are dropped.
%
% DOPAMINE datasets: cell type / drug from the trial folder name.
% dopamine_ionto_pilot (GRAB-DA2m, no ft/stim trace) is skipped. walking
% and redo are treated as DA (folder names), though their claude scripts'
% headers say ATP.
%
% P2X2 datasets (every .data/lpsp_p2x2*.mat, incl. the redo/syt8s trials
% built by lpsp_p2x2_unsaved_build.m in this folder): checked against
% a scan of every Z:\pablo\lpsp_p2x2* trial folder (2026-10-29). Each trial
% exists once on Z:, but the .mat pulls captured many trials before AND
% after folders were reorganized (todo\, blind _ low\, reredo\blind\,
% "fly N" -> "fly N _ control"), so the same trial can sit under two
% paths. Trials are therefore deduplicated by root + trial folder name
% (newest .mat copy wins), and a fly is root + date + fly NUMBER (genotype
% suffix on the fly folder ignored).
% Genotype: "empty"/"control" in fly or trial folder -> empty>P2X2;
% "lpsp" or 2024 "th_p2x2" -> lpsp>P2X2 (th>P2X2 = LPsP stimulation,
% confirmed by Pablo 2026-10-29); blinded trials ("blind" in the trial
% name) whose fly folder carries no tag -> lpsp>P2X2 (unblinding convention:
% only one genotype gets tagged).
% Excluded: anything under "open loop", "exclusions" or "pilot (dont use)"
% (the folders lpsp_p2x2_revisit.m also skipped), and 2024 "_low" (low-ATP)
% trials. The 2024 trials have no ft.stims, so their onsets are detected
% from the red channel (findpeaks, prominence 9 on a 20-frame forward mean,
% as in lpsp_p2x2_script_minimal.m) -- every detected event is an ejection
% by construction, so they have no failed-stim comparison.
% Stim matching: one stratum per experiment (dataset root) x stim length,
% since each experiment interleaved both genotypes; a stratum is shown only
% if both genotypes have >= 2 flies with ejections.
% Effect = lpsp per-fly mean minus empty per-fly mean; SEM combined in
% quadrature; p = two-sample t-test across flies. The trials used are
% written to tables/p2x2_trials_used.csv.
%
% Each fly is also called ATP vs dopamine (dye SNR, EPG-green response,
% name label; see the "call each fly's ejected substance" section), and a
% stim-order analysis asks whether later ejections (same fly, same side)
% have a bigger dye-adjusted effect.
%
% Run from the repo root (batch-friendly, ~4 min with the cache):
%   addpath('ugly_figures/ionto_celltypes'); fig_ionto_celltypes
% The P2X2 extraction loads ~13 GB, so it is cached to
% .data/fig_dopamine_ionto_p2x2_cache.mat (rebuild_p2x2 = true to redo, ~20 min).
%
% Figures (-> ugly_figures/ionto_celltypes/*.pdf; tables -> tables/):
%  1 peri-stim green + red time courses per condition
%  2 per-fly summary of the 0-3 s / 3-15 s / ipsi / contra effects
%  3 ejection QC: per-stim z + ejection fraction
%  4 artifact checks: fast (0.5 s) response; kinetics; red amp vs green
%  5 glomerulus heatmaps, ejection side aligned to rows 1-16
%  6 ATP vs dopamine evidence and calls per fly
%  7 dye and dye-adjusted response vs ejection number (same fly, same side)

if ~exist('rebuild_p2x2','var'); rebuild_p2x2 = false; end
data_dir = '.data';
out_dir  = fullfile('ugly_figures','ionto_celltypes');    % this folder: scripts, PDFs
tab_dir  = fullfile(out_dir,'tables');                      % ledgers, per-fly calls, stats
if ~isfolder(tab_dir); mkdir(tab_dir); end

t_common = -10:0.1:30;
win.base = t_common>=-5 & t_common<0;
win.pk   = t_common>0  & t_common<=3;
win.late = t_common>3  & t_common<=15;

%% dopamine iontophoresis: load + per-stim extraction
ds_names = {'dopamine_ionto_epg_pilot','dopamine_ionto_lpsp_pilot','dopamine_ionto_20250722', ...
            'dopamine_ionto_redo_20260518','dopamine_ionto_walking_reg_20260831'};
S = empty_stims();
S_trial = {}; S_snum = []; % trial path + order within trial, parallel to S (for the stim-order analysis)
for k = 1:numel(ds_names)
    tmp = load(fullfile(data_dir,[ds_names{k} '.mat']),'all_data');
    ad  = tmp.all_data;
    fprintf('loaded %s (%d trials)\n',ds_names{k},numel(ad));
    for i = 1:numel(ad)
        [ct,dr] = da_trial_labels(ad(i).meta,ds_names{k});
        if ~has_stims(ad(i)) || ~has_atp(ad(i)); continue; end
        [on_t,slen] = stim_onsets(ad(i).ft);
        s = extract_stims(ad(i),on_t,slen,t_common,win,true);
        [s.ds] = deal(ds_names{k}); [s.fly] = deal([ds_names{k} ':' fly_folder(ad(i).meta)]);
        [s.cond] = deal([ct ' | ' dr]); [s.geno] = deal('');
        S = [S s]; %#ok<AGROW>
        S_trial = [S_trial repmat({ad(i).meta},1,numel(s))]; S_snum = [S_snum 1:numel(s)]; %#ok<AGROW>
    end
end
clear tmp ad

%% LPsP>P2X2 vs empty>P2X2: load (deduplicated) + per-stim extraction, cached
cache_file = fullfile(data_dir,'fig_dopamine_ionto_p2x2_cache.mat');
if isfile(cache_file) && ~rebuild_p2x2
    load(cache_file,'P','p2x2_trials')
    fprintf('loaded P2X2 stim cache (%d stims, %d trials)\n',numel(P),numel(p2x2_trials));
else
    fs = dir(fullfile(data_dir,'lpsp_p2x2*.mat'));
    [~,ord] = sort([fs.datenum],'descend'); fs = fs(ord); % newest pull first: its copy of a trial wins
    seen = containers.Map();
    p2x2_trials = struct('key',{},'path',{},'file',{},'fly',{},'geno',{},'cond',{},'used',{},'why',{},'n_stims',{},'n_ej',{});
    P = empty_stims();
    for k = 1:numel(fs)
        tmp = load(fullfile(fs(k).folder,fs(k).name),'all_data'); ad = tmp.all_data; clear tmp
        n_new = 0;
        for i = 1:numel(ad)
            tr  = regexprep(ad(i).meta,'\\registration[^\\]*\\?$','');
            key = p2x2_trial_key(tr);
            if isKey(seen,key); continue; end
            seen(key) = fs(k).name; n_new = n_new+1;
            [geno,ind,why] = p2x2_labels(tr);
            fid = p2x2_fly_id(tr);
            rec = struct('key',key,'path',tr,'file',fs(k).name,'fly',fid,'geno',geno,'cond','', ...
                'used',false,'why',why,'n_stims',0,'n_ej',0);
            if isempty(why) && ~has_atp(ad(i)); rec.why = 'no red channel'; end
            if isempty(rec.why)
                if has_stims(ad(i))
                    [on_t,slen] = stim_onsets(ad(i).ft); from_stims = true;
                else
                    on_t = atp_detect_onsets(ad(i)); slen = nan(size(on_t)); from_stims = false;
                end
                s = extract_stims(ad(i),on_t,slen,t_common,win,from_stims);
                if isempty(s)
                    rec.why = 'no stims found';
                else
                    [s.ds] = deal(fs(k).name); [s.fly] = deal(fid); [s.geno] = deal(geno);
                    for j = 1:numel(s); s(j).cond = p2x2_stratum(tr,ind,s(j).stimlen,from_stims); end
                    P = [P s]; %#ok<AGROW>
                    rec.used = true; rec.cond = strjoin(unique({s.cond}),';');
                    rec.n_stims = numel(s); rec.n_ej = sum([s.ej]);
                end
            end
            p2x2_trials(end+1) = rec; %#ok<AGROW>
        end
        fprintf('%s: %d trials, %d new\n',fs(k).name,numel(ad),n_new);
        clear ad
    end
    save(cache_file,'P','p2x2_trials','-v7.3')
end
% exclusion status and genotype from where each trial lives on Z: TODAY
% (tables/p2x2_z_trial_status.csv, from the 2026-10-29 scan of every
% Z:\pablo\lpsp_p2x2* trial folder): many .mat copies were saved before
% trials were moved into "pilot (dont use)"/"exclusions" or their fly
% folders were unblinded, so the stored path can be stale. P was appended
% in p2x2_trials order (n_stims per used trial), which links stims to trials.
zs = readtable(fullfile(tab_dir,'p2x2_z_trial_status.csv'),'TextType','char','Delimiter',',');
zmap = containers.Map(zs.key,num2cell(1:height(zs)));
stim_trial = repelem(find([p2x2_trials.used]),[p2x2_trials([p2x2_trials.used]).n_stims]);
P_trial = {p2x2_trials(stim_trial).path};
P_snum  = cell2mat(arrayfun(@(n) 1:n,[p2x2_trials([p2x2_trials.used]).n_stims],'UniformOutput',false));
keep = true(1,numel(P));
for k = 1:numel(p2x2_trials)
    if ~isKey(zmap,p2x2_trials(k).key); p2x2_trials(k).z_status = 'not on Z'; continue; end
    z = zmap(p2x2_trials(k).key);
    p2x2_trials(k).z_status = zs.z_status{z}; p2x2_trials(k).z_path = zs.z_path{z};
    if ~isempty(zs.z_geno{z}) && ~strcmp(zs.z_geno{z},p2x2_trials(k).geno) && p2x2_trials(k).used
        fprintf('genotype from current Z folder: %s -> %s (%s)\n',p2x2_trials(k).geno,zs.z_geno{z},zs.z_path{z});
        p2x2_trials(k).geno = zs.z_geno{z}; [P(stim_trial==k).geno] = deal(zs.z_geno{z});
    end
    if ~strcmp(zs.z_status{z},'ok') && p2x2_trials(k).used
        p2x2_trials(k).used = false; p2x2_trials(k).why = zs.z_status{z};
        keep(stim_trial==k) = false;
    end
end
fprintf('dropped %d stims from %d trials now under excluded Z folders\n',sum(~keep),numel(unique(stim_trial(~keep))));
P = P(keep); P_trial = P_trial(keep); P_snum = P_snum(keep);
writetable(struct2table(p2x2_trials),fullfile(tab_dir,'p2x2_trials_used.csv'));

%% pool everything, define conditions
A    = [S P];
trial_all = [S_trial P_trial]'; snum_all = [S_snum P_snum]';
G    = vertcat(A.g); Rd = vertcat(A.r);
fly  = {A.fly}'; cond = {A.cond}'; geno = {A.geno}';
ej   = [A.ej]' & [A.isi]'>=9.5;
fl   = ~[A.ej]' & [A.z]'<1 & [A.isi]'>=9.5; % clearly failed ejections
fprintf('%d stims total (%d DA, %d P2X2); %d ejections\n',numel(A),numel(S),numel(P),sum(ej));

%% per-fly ejection features (for calling ATP vs dopamine per fly)
% ATP ejection reliably raises EPG GCaMP fluorescence, and the ATP pipette
% dye is much brighter than the dopamine one, so per fly: median red
% evoked amplitude / z over its stims, and mean green 0-3 s response over
% its ejection stims. Name label: trial folder says da/dopamine.
[~,i_g] = min(abs(t_common-1)); % green at ~1 s, before slow effects
uf = unique(fly);
FF = table('Size',[numel(uf) 15],'VariableTypes',{'string','string','string','string','double','double','double','double','double','double','double','logical','double','double','double'}, ...
    'VariableNames',{'fly','ds','cond','geno','stimlen','n_stims','n_ej','red_amp','red_z','green_0_3','green_1s','name_da','red_amp_ej','red_z_ej','red_z_p90'});
p2x2_paths = containers.Map();
if exist('p2x2_trials','var') && isstruct(p2x2_trials)
    for k = 1:numel(p2x2_trials)
        f = p2x2_trials(k).fly;
        if isKey(p2x2_paths,f); p2x2_paths(f) = [p2x2_paths(f) {p2x2_trials(k).path}]; else; p2x2_paths(f) = {p2x2_trials(k).path}; end
    end
end
for k = 1:numel(uf)
    j  = strcmp(fly,uf{k}); je = j & ej;
    sl = [A(j).stimlen];
    FF.fly(k) = uf{k}; FF.ds(k) = A(find(j,1)).ds; FF.cond(k) = strjoin(unique(cond(j)),';'); FF.geno(k) = A(find(j,1)).geno;
    FF.stimlen(k) = median(sl,'omitnan'); FF.n_stims(k) = sum(j); FF.n_ej(k) = sum(je);
    FF.red_amp(k) = median([A(j).atp_amp]); FF.red_z(k) = median([A(j).z]);
    FF.green_0_3(k) = mean([A(je).early]); FF.green_1s(k) = mean(G(je,i_g),'omitnan');
    FF.red_amp_ej(k) = median([A(je).atp_amp]); FF.red_z_ej(k) = median([A(je).z]);
    FF.red_z_p90(k)  = prctile([A(j).z],90); % best ejections, robust to many failed pulses
    if isKey(p2x2_paths,uf{k})
        FF.name_da(k) = any(~cellfun(@isempty,regexpi(p2x2_paths(uf{k}),'(_da_|_da$|dopamine)','once')));
    else
        FF.name_da(k) = contains(FF.cond(k),'DA'); % dopamine_ionto sets: drug parsed from the trial name
    end
end
%% call each fly's ejected substance: ATP vs dopamine, with a confidence
% (all trials of one fly used the same substance). Three votes, +1 = ATP,
% -1 = dopamine, 0 = uninformative, NaN = not available:
%  dye SNR  : the ATP pipette dye is far brighter. ATP if ejection-stim red
%             amplitude >= 8 and z >= 30; dopamine if < 5 and < 30.
%             Calibrated on the named dopamine sets (amp 0.5-5, z 5-29) vs
%             the p2x2 sets (amp 10-60, z 40-230).
%  EPG green: ATP ejection reliably raises EPG GCaMP. Strong ATP if green at
%             ~1 s >= 0.05, weak ATP (+0.5) if 0.015-0.05, dopamine if <= 0.
%             EPG-GCaMP flies only (LPsP and GRAB-DA also rise with dopamine).
%  label    : -1 if the trial/dataset name says da/dopamine.
% high = >= 2 votes, all one sign, at least one strong; medium = one
% strong vote alone, or agreeing votes where the strong one is missing;
% low = votes of both signs (or only a weak vote). Flies with no ejection
% stim get no call. ACh-only flies are not called.
% Known ground truth (Pablo, 2026-10-29): the walking p2x2 0.5 s flies were
% ATP at a low ejection strength (so their dye looks dopamine-dim and the
% EPG increase is weak); they are called "ATP / confirmed" regardless of
% votes, which are still computed and reported.
atp_confirmed = startsWith(FF.fly,"lpsp_p2x2_walking|") & contains(FF.cond,"| 0.5 s");
FF.vote_snr = nan(height(FF),1); FF.vote_green = nan(height(FF),1); FF.vote_label = nan(height(FF),1);
FF.call = repmat("none",height(FF),1); FF.conf = repmat("",height(FF),1);
for k = 1:height(FF)
    if FF.n_ej(k) == 0 || (contains(FF.cond(k),'ACh') && ~contains(FF.cond(k),'DA')); continue; end
    a = FF.red_amp_ej(k); z = FF.red_z_ej(k);
    if a >= 8 && z >= 30; FF.vote_snr(k) = 1; elseif a < 5 && z < 30; FF.vote_snr(k) = -1; else; FF.vote_snr(k) = 0; end
    if contains(FF.cond(k),'EPG') && ~contains(FF.cond(k),'GRAB')
        g = FF.green_1s(k);
        if g >= 0.05; FF.vote_green(k) = 1; elseif g >= 0.015; FF.vote_green(k) = 0.5; elseif g <= 0; FF.vote_green(k) = -1; else; FF.vote_green(k) = 0; end
    end
    if FF.name_da(k) || (FF.geno(k)=="" && contains(FF.cond(k),'DA')); FF.vote_label(k) = -1; end
    v = [FF.vote_snr(k) FF.vote_green(k) FF.vote_label(k)]; v = v(~isnan(v) & v~=0);
    s = sum(v);
    if s > 0; FF.call(k) = "ATP"; elseif s < 0; FF.call(k) = "dopamine"; else; FF.call(k) = "uncertain"; end
    if isempty(v) || (any(v>0) && any(v<0)) || all(abs(v) < 1)
        FF.conf(k) = "low";
    elseif numel(v) >= 2 && any(abs(v)==1)
        FF.conf(k) = "high";
    else
        FF.conf(k) = "medium";
    end
    if atp_confirmed(k); FF.call(k) = "ATP"; FF.conf(k) = "confirmed"; end
end
writetable(FF,fullfile(tab_dir,'per_fly_ejection_features.csv'));
j = atp_confirmed & FF.n_ej>0;
fprintf('\nwalking 0.5 s p2x2 (confirmed ATP): EPG green at ~1 s > 0 in %d/%d flies (empty %d/%d, lpsp %d/%d); median %.3f\n', ...
    sum(FF.green_1s(j)>0),sum(j),sum(FF.green_1s(j&FF.geno=="empty")>0),sum(j&FF.geno=="empty"), ...
    sum(FF.green_1s(j&FF.geno=="lpsp")>0),sum(j&FF.geno=="lpsp"),median(FF.green_1s(j)));
fprintf('\n=== substance calls per group (flies with >=1 ejection) ===\n');
grp_cond = FF.cond; grp_cond(FF.geno~="") = FF.cond(FF.geno~="") + " | " + FF.geno(FF.geno~="");
[ug,~,ig] = unique(grp_cond);
for k = 1:numel(ug)
    j = ig==k & FF.call~="none";
    fprintf('%-58s ATP %2d (h%d m%d l%d)  DA %2d (h%d m%d l%d)  unc %d\n',ug(k), ...
        sum(FF.call(j)=="ATP"),sum(FF.call(j)=="ATP"&FF.conf(j)=="high"),sum(FF.call(j)=="ATP"&FF.conf(j)=="medium"),sum(FF.call(j)=="ATP"&FF.conf(j)=="low"), ...
        sum(FF.call(j)=="dopamine"),sum(FF.call(j)=="dopamine"&FF.conf(j)=="high"),sum(FF.call(j)=="dopamine"&FF.conf(j)=="medium"),sum(FF.call(j)=="dopamine"&FF.conf(j)=="low"), ...
        sum(FF.call(j)=="uncertain"));
end
% p2x2 flies called dopamine with HIGH confidence are moved out of the
% lpsp - empty ATP comparisons into their own dopamine column (genotype
% cleared so they pool like a dopamine set); medium/low dopamine calls stay
% in their ATP stratum and only show up in the title counts. Dopamine-set
% flies called ATP with high confidence (none so far) would be dropped.
call_of = containers.Map(cellstr(FF.fly),cellstr(FF.call)); conf_of = containers.Map(cellstr(FF.fly),cellstr(FF.conf));
da_p2x2 = FF.fly(FF.geno~="" & FF.call=="dopamine" & FF.conf=="high");
% one column per indicator x genotype, e.g. "EPG syt8s, lpsp>P2X2 flies called DA"
da_p2x2_conds = {};
for k = 1:numel(da_p2x2)
    j  = strcmp(fly,da_p2x2(k));
    c0 = strjoin(unique(cond(j)),';');
    if contains(c0,'syt8s'); ind_k = 'syt8s'; elseif contains(c0,'8m'); ind_k = '8m'; else; ind_k = '7f'; end
    ck = sprintf('EPG %s, %s>P2X2 flies called DA',ind_k,geno{find(j,1)});
    fprintf('p2x2 fly called dopamine (high) -> own column: %s (%s)\n',da_p2x2(k),c0);
    cond(j) = {ck}; geno(j) = {''}; [A(j).cond] = deal(ck); [A(j).geno] = deal('');
    da_p2x2_conds{end+1} = ck; %#ok<SAGROW>
end
da_p2x2_conds = flip(unique(da_p2x2_conds)); % syt8s first
bad_fly = FF.fly(FF.geno=="" & FF.call=="ATP" & FF.conf=="high");
if ~isempty(bad_fly); fprintf('dropping %d dopamine-set flies called ATP (high):\n',numel(bad_fly)); fprintf('  %s\n',bad_fly); end
drop = ismember(fly,cellstr(bad_fly)); ej(drop) = false; fl(drop) = false;

% each condition: {label, cond string, type}; type 'da' = single group,
% 'diff' = lpsp>P2X2 minus empty>P2X2 within that stim-matched stratum
C = {'LPsP DA',                 'LPsP (syt8m) | DA',            'da'
     'GRAB-DA DA',              'EPG>GRAB-DA2m | DA',           'da'
     'EPG 7f DA',               'EPG (7f) | DA',                'da'
     'EPG 8m DA',               'EPG (syt8m) | DA',             'da'
     'EPG 8s DA (empty>P2X2)',  'EPG (syt8s), empty>P2X2 | DA', 'da'};
% dopamine ejections found inside the p2x2 data, next to the matching
% dopamine_ionto column (e.g. syt8s lpsp>P2X2 beside syt8s empty>P2X2)
for k = 1:numel(da_p2x2_conds)
    lbl = regexprep(da_p2x2_conds{k},'^EPG (\S+), (\S+)>P2X2 flies called DA$','EPG $1 DA ($2>P2X2, called DA)');
    lbl = strrep(lbl,'syt8s','8s');
    C(end+1,:) = {lbl,da_p2x2_conds{k},'da'}; %#ok<SAGROW>
end
C = [C; {'EPG 8m ACh','EPG (syt8m) | ACh','da'; 'EPG 8m DA+ACh','EPG (syt8m) | DA+ACh','da'}];
% LPsP-stim strata: every experiment x stim length with >= 2 flies (with
% ejections) of BOTH genotypes, in a fixed experiment order
strata = unique(cond(~cellfun(@isempty,geno)));
exp_order = {'2024','redo','reredo','syt8s','walking'};
[~,o] = sort(cellfun(@(s) find(startsWith(s,strcat(exp_order,' ')),1), strata)); strata = strata(o);
for k = 1:numel(strata)
    nL = numel(unique(fly(ej & strcmp(cond,strata{k}) & strcmp(geno,'lpsp'))));
    nE = numel(unique(fly(ej & strcmp(cond,strata{k}) & strcmp(geno,'empty'))));
    fprintf('stratum %-34s lpsp %2d flies, empty %2d flies%s\n',strata{k},nL,nE,repmat(' -> dropped',1,nL<2||nE<2));
    if nL >= 2 && nE >= 2; C(end+1,:) = {['LPsP stim: ' strata{k}],strata{k},'diff'}; end %#ok<SAGROW>
end
nC = size(C,1); is_diff = strcmp(C(:,3),'diff');
sel_of = @(c,g,base) base & strcmp(cond,C{c,2}) & strcmp(geno,g);

cg = [0.13 0.55 0.13]; cr = [0.80 0.15 0.15]; cf = [0.55 0.55 0.55];
cl = [0.90 0.45 0.05]; ce = [0.10 0.40 0.80]; % lpsp>P2X2, empty>P2X2

fprintf('\n%-26s %12s %12s %9s %7s\n','condition','n lpsp/DA','n empty','dFF 0-3s','p');
for c = 1:nC
    [v,se,p,nf] = cond_effect(c,'early',A,C,is_diff,ej,cond,geno,fly);
    fprintf('%-26s %12s %12s %9.3f %7.3g  (sem %.3f)\n',C{c,1},sprintf('%d flies',nf(1)),sprintf('%d flies',nf(2)),v,p,se);
end

%% figure 1: time courses
figure(1); clf; set(gcf,'Position',[10 20 2700 760],'Color','w')
tl = tiledlayout(2,nC,'TileSpacing','compact','Padding','compact');
for c = 1:nC
    ax = nexttile(tl,c); hold(ax,'on')
    patch(ax,[0 3 3 0],[-1 -1 1 1],'k','FaceAlpha',.05,'EdgeColor','none')
    if ~is_diff(c)
        a = sel_of(c,'',ej); b = sel_of(c,'',fl);
        [m2,se2,nf2] = fly_mean(G,b,fly); shade(ax,t_common,m2,se2,cf,1)
        [m,se,nf]    = fly_mean(G,a,fly); shade(ax,t_common,m,se,cg,2)
        ttl = sprintf('%s\nejection: %d flies/%d stims\nfailed: %d flies/%d stims\n%s',C{c,1},nf,sum(a),nf2,sum(b), ...
            subst_line(unique(fly(a)),call_of,conf_of));
    else
        aL = sel_of(c,'lpsp',ej); aE = sel_of(c,'empty',ej);
        [mL,seL,nL] = fly_mean(G,aL,fly); [mE,seE,nE] = fly_mean(G,aE,fly);
        shade(ax,t_common,mE,seE,ce,1); shade(ax,t_common,mL,seL,cl,1)
        shade(ax,t_common,mL-mE,sqrt(seL.^2+seE.^2),[0 0 0],2)
        ttl = sprintf('%s\nlpsp: %d flies/%d stims\nempty: %d flies/%d stims\n%s',C{c,1},nL,sum(aL),nE,sum(aE), ...
            subst_line(unique(fly(aL|aE)),call_of,conf_of));
    end
    xline(ax,0,':k'); yline(ax,0,':k'); xlim(ax,[-10 30]); ylim(ax,[-0.25 0.3])
    title(ax,ttl,'FontSize',7)
    if c==1; ylabel(ax,'green \DeltaF/F (mean of 32 glom.)'); end
    ax = nexttile(tl,c+nC); hold(ax,'on')
    if ~is_diff(c)
        [m2,se2] = fly_mean(Rd,b,fly); shade(ax,t_common,m2,se2,cf,1)
        [m,se]   = fly_mean(Rd,a,fly); shade(ax,t_common,m,se,cr,2)
    else
        [m,se] = fly_mean(Rd,aE,fly); shade(ax,t_common,m,se,ce,1.5)
        [m,se] = fly_mean(Rd,aL,fly); shade(ax,t_common,m,se,cl,1.5)
    end
    xline(ax,0,':k'); yline(ax,0,':k'); xlim(ax,[-10 30])
    xlabel(ax,'time from stim onset (s)')
    if c==1; ylabel(ax,'red (dye) \Sigma\DeltaF/F'); end
end
title(tl,{'Peri-stim responses, mean \pm SEM across flies', ...
    'DA columns: green = ejection stims, gray = failed ejections.   LPsP-stim columns: orange = lpsp>P2X2, blue = empty>P2X2, black = lpsp - empty.   Shaded box = 0-3 s window'})
exportgraphics(gcf,fullfile(out_dir,'fig01_timecourses.pdf'))

%% figure 2: per-fly summary
metrics = {'early','late','ipsi','contra'};
mlabels = {'whole PB, 0-3 s','whole PB, 3-15 s','ejection side, 0-3 s','opposite side, 0-3 s'};
figure(2); clf; set(gcf,'Position',[10 20 2200 560],'Color','w')
tl = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for q = 1:4
    ax = nexttile(tl); hold(ax,'on')
    xline(ax,sum(~is_diff)+.5,'-','Color',[.6 .6 .6])
    for c = 1:nC
        [v,se,p,~,dots] = cond_effect(c,metrics{q},A,C,is_diff,ej,cond,geno,fly);
        if is_diff(c); bc = [1 .85 .7]; dc = cl; else; bc = [.85 .85 .85]; dc = [0 0 0]; end
        bar(ax,c,v,.6,'FaceColor',bc,'EdgeColor','none')
        scatter(ax,c+0.15*(rand(numel(dots),1)-.5),dots,18,dc,'filled','MarkerFaceAlpha',.6)
        if ~isnan(se); errorbar(ax,c,v,se,'k','LineWidth',1.5,'CapSize',5); end
        if ~isnan(p); text(ax,c,0.38,sprintf('p=%.2g',p),'HorizontalAlignment','center','FontSize',6.5); end
    end
    yline(ax,0,':k'); ylim(ax,[-0.3 0.42]); xlim(ax,[.4 nC+.6])
    set(ax,'XTick',1:nC,'XTickLabel',C(:,1),'XTickLabelRotation',45,'FontSize',8)
    ylabel(ax,'\DeltaF/F (vs -5..0 s)'); title(ax,mlabels{q})
end
title(tl,{'Effect per fly, bar = mean \pm SEM; ejection stims only', ...
    'left of line: DA iontophoresis (one-sample t-test across flies).  right: LPsP stimulation = lpsp>P2X2 minus stim-matched empty>P2X2 (orange dots = lpsp flies minus empty mean; two-sample t-test)'})
exportgraphics(gcf,fullfile(out_dir,'fig02_summary_per_fly.pdf'))

%% figure 3: ejection QC
% groups: each DA dataset, then each P2X2 stratum x genotype
z = [A.z]'; ds_all = {A.ds}'; isi_ok = [A.isi]'>=9.5; is_p2x2 = ~cellfun(@isempty,geno);
grp_lbl = {'DA: EPG 7f pilot','DA: LPsP pilot','DA: 2025 repilot','DA: redo','DA: walking'};
grp_idx = cellfun(@(x) strcmp(ds_all,x) & isi_ok, ds_names,'UniformOutput',false);
for c = find(is_diff)'
    for g = {'lpsp','empty'}
        grp_lbl{end+1} = sprintf('%s (%s)',strrep(C{c,1},'LPsP stim: ',''),g{1}); %#ok<SAGROW>
        grp_idx{end+1} = is_p2x2 & strcmp(cond,C{c,2}) & strcmp(geno,g{1}) & isi_ok; %#ok<SAGROW>
    end
end
figure(3); clf; set(gcf,'Position',[10 20 1900 520],'Color','w')
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl); hold(ax,'on')
for k = 1:numel(grp_idx)
    idx = find(grp_idx{k});
    if numel(idx) > 600; idx = idx(randperm(numel(idx),600)); end % thin dense groups for plotting
    e = [A(idx).ej]'; col = repmat(cf,numel(idx),1); col(e,:) = repmat(cr,sum(e),1);
    scatter(ax,k+0.5*(rand(numel(idx),1)-.5),max(min(z(idx),200),-10),6,col,'filled','MarkerFaceAlpha',.5)
end
yline(ax,3,'--k','z = 3'); ylim(ax,[-10 60]); xline(ax,numel(ds_names)+.5,'-','Color',[.6 .6 .6])
set(ax,'XTick',1:numel(grp_lbl),'XTickLabel',grp_lbl,'XTickLabelRotation',40,'FontSize',8)
ylabel(ax,'red-channel z (clipped)'); title(ax,'per-stim ejection z (red = called ejection; dense groups thinned to 600 points)')
ax = nexttile(tl); hold(ax,'on')
frac = cellfun(@(x) mean([A(x).ej]), grp_idx); nst = cellfun(@sum, grp_idx);
b = bar(ax,frac,'FaceColor','flat'); b.CData = repmat(cr,numel(frac),1);
b.CData(numel(ds_names)+1:2:end,:) = repmat(cl,numel(frac(numel(ds_names)+1:2:end)),1);
b.CData(numel(ds_names)+2:2:end,:) = repmat(ce,numel(frac(numel(ds_names)+2:2:end)),1);
for k = 1:numel(frac); text(ax,k,frac(k)+.03,sprintf('%d/%d',round(frac(k)*nst(k)),nst(k)),'HorizontalAlignment','center','FontSize',7); end
xline(ax,numel(ds_names)+.5,'-','Color',[.6 .6 .6])
set(ax,'XTick',1:numel(grp_lbl),'XTickLabel',grp_lbl,'XTickLabelRotation',40,'FontSize',8); ylim(ax,[0 1.1])
ylabel(ax,'fraction of stims with detected ejection')
title(ax,'ejection success (stims \geq9.5 s apart); 2024 7f events are detected FROM the red channel, so 100% by construction')
exportgraphics(gcf,fullfile(out_dir,'fig03_ejection_qc.pdf'))

%% figure 4: artifact checks
[~,i05] = min(abs(t_common-0.5));
figure(4); clf; set(gcf,'Position',[10 20 2500 800],'Color','w')
tl = tiledlayout(2,nC,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,[1 6]); hold(ax,'on')
for c = 1:nC
    if is_diff(c)
        grp = {sel_of(c,'lpsp',ej),cl; sel_of(c,'empty',ej),ce};
    else
        grp = {sel_of(c,'',ej),cg; sel_of(c,'',fl),cf};
    end
    for e = 1:2
        f = unique(fly(grp{e,1})); if isempty(f); continue; end
        pf = cellfun(@(x) mean(G(grp{e,1} & strcmp(fly,x),i05),'omitnan'),f);
        x = c + (e-1.5)*0.3;
        bar(ax,x,mean(pf),.28,'FaceColor',grp{e,2},'EdgeColor','none','FaceAlpha',.6)
        scatter(ax,x+0.08*(rand(numel(pf),1)-.5),pf,10,'k','filled')
    end
end
yline(ax,0,':k'); xline(ax,sum(~is_diff)+.5,'-','Color',[.6 .6 .6]); ylim(ax,[-0.3 0.3])
set(ax,'XTick',1:nC,'XTickLabel',C(:,1),'XTickLabelRotation',35,'FontSize',8)
ylabel(ax,'green \DeltaF/F at t = 0.5 s')
title(ax,'fast response at 0.5 s.  DA: green = ejection, gray = failed.  LPsP stim: orange = lpsp>P2X2, blue = empty>P2X2')
ax = nexttile(tl,7,[1 nC-6]); hold(ax,'on')
for c = [1 2 4 6 find(is_diff)']
    if is_diff(c)
        m = fly_mean(G,sel_of(c,'lpsp',ej),fly) - fly_mean(G,sel_of(c,'empty',ej),fly); ls = '-';
    else
        m = fly_mean(G,sel_of(c,'',ej),fly); ls = ':';
    end
    m = smoothdata(m,'movmean',3);
    plot(ax,t_common,m./max(abs(m(t_common>0 & t_common<15))),'LineStyle',ls,'LineWidth',1.5,'DisplayName',C{c,1})
end
m = fly_mean(Rd,sel_of(1,'',ej),fly);
plot(ax,t_common,m./max(m),'k--','LineWidth',2,'DisplayName','red dye (LPsP DA trials)')
xline(ax,0,':k','HandleVisibility','off'); yline(ax,0,':k','HandleVisibility','off')
xlim(ax,[-5 20]); ylim(ax,[-1.2 1.2]); legend(ax,'Location','eastoutside','FontSize',7)
xlabel(ax,'time from stim onset (s)'); ylabel(ax,'normalized to own |peak| (0-15 s)')
title(ax,'kinetics (dotted = DA groups, solid = LPsP stim, lpsp - empty)')
for c = 1:nC
    ax = nexttile(tl,nC+c); hold(ax,'on')
    if is_diff(c)
        grp = {sel_of(c,'empty',ej),ce; sel_of(c,'lpsp',ej),cl};
    else
        grp = {sel_of(c,'',ej),cg};
    end
    ttl = C{c,1};
    for e = 1:size(grp,1)
        x = [A(grp{e,1}).atp_amp]'; y = [A(grp{e,1}).early]';
        if numel(x) > 1500; k = randperm(numel(x),1500); x = x(k); y = y(k); end
        scatter(ax,x,y,6,grp{e,2},'filled','MarkerFaceAlpha',.4)
        if numel(x)>2; ttl = sprintf('%s\nr=%.2f',ttl,corr(x,y,'rows','complete')); end
    end
    title(ax,ttl,'FontSize',7); yline(ax,0,':k'); xlabel(ax,'red amp (0-3 s)','FontSize',7)
    ylim(ax,[-0.4 0.6])
    if c==1; ylabel(ax,'green \DeltaF/F (0-3 s)'); end
end
title(tl,'Artifact checks')
exportgraphics(gcf,fullfile(out_dir,'fig04_artifact_checks.pdf'))

%% figure 5: glomerulus heatmaps, ejection side on top (rows 1-16)
figure(5); clf; set(gcf,'Position',[10 20 2600 440],'Color','w')
tl = tiledlayout(1,nC,'TileSpacing','compact','Padding','compact');
cmap = [linspace(0,1,128)',linspace(0,1,128)',ones(128,1); ones(128,1),linspace(1,0,128)',linspace(1,0,128)'];
for c = 1:nC
    if is_diff(c)
        [ML,nL] = fly_heatmap(A,sel_of(c,'lpsp',ej),fly); [ME,nE] = fly_heatmap(A,sel_of(c,'empty',ej),fly);
        M = ML - ME; ttl = sprintf('%s\nlpsp-empty (%d/%d flies)',C{c,1},nL,nE);
    else
        [M,nf] = fly_heatmap(A,sel_of(c,'',ej),fly); ttl = sprintf('%s\n(%d flies)',C{c,1},nf);
    end
    ax = nexttile(tl); imagesc(ax,t_common,1:32,M); hold(ax,'on')
    yline(ax,16.5,'k-','LineWidth',1); xline(ax,0,':k')
    colormap(ax,cmap); lim = max(0.03,prctile(abs(M(:)),99)); clim(ax,[-lim lim]); colorbar(ax)
    xlim(ax,[-5 20]); set(ax,'YTick',[8.5 24.5],'YTickLabel',{'ejection side','opposite'},'FontSize',7)
    title(ax,ttl,'FontSize',8); xlabel(ax,'time (s)')
end
title(tl,'Glomerulus \DeltaF/F around ejections, ejection hemisphere aligned to top (mean across flies, own color scale; LPsP-stim panels = lpsp>P2X2 minus empty>P2X2)')
exportgraphics(gcf,fullfile(out_dir,'fig05_glomerulus_heatmaps.pdf'))
%% stim order: does the response grow over successive ejections, at matched dye?
% Each ejection stim is ranked within its fly: trials in chronological
% order (YYYYMMDD-N from the trial folder), then stims within the trial --
% i.e. cumulative ejection count across the whole session. Red (dye)
% amplitude 0-3 s is the proxy for how much ATP/DA was ejected. Per
% condition, a linear mixed model with a random intercept per fly:
%   DA columns  : resp ~ dye + log2(rank) + (1|fly)
%   LPsP stim   : resp ~ dye + log2(rank)*is_lpsp + (1|fly)
%                 (the interaction = does the lpsp - empty effect grow)
% log2 because ranks run 1..8 in some sets and 1..800 in walking. The same
% model with within-trial stim number (1..n) instead of rank is reported as
% a second check, and dye ~ log2(rank) tests whether the ejections
% themselves drifted (pipette clogging/depletion). resp = 0-3 s and
% 3-15 s whole-PB green dF/F.
tnum = cellfun(@trial_sort_key,trial_all); % sortable YYYYMMDD*1e4 + trial number
% The count restarts at 1 whenever the ejection side switches: each trial's
% side = majority side of its ejection stims (evoked red per PB half, as
% in extract_stims); within a fly, trials are walked chronologically and a
% new segment starts every time a trial's side differs from the previous
% trial's (including switching back to an earlier side).
rank_all = nan(numel(A),1); seg_all = nan(numel(A),1); right_all = [A.right]';
n_switch = 0; n_fly_switch = 0;
for f = unique(fly(ej))'
    j = find(ej & strcmp(fly,f{1}));
    [ut,~,it] = unique(tnum(j));                        % trials of this fly, chronological
    tside = accumarray(it,double(right_all(j)),[],@(x) double(mean(x) > 0.5));
    seg = cumsum([1; diff(tside) ~= 0]);                % new segment at every side switch
    n_switch = n_switch + seg(end) - 1; n_fly_switch = n_fly_switch + (seg(end) > 1);
    seg_all(j) = seg(it);
    for s = 1:seg(end)
        js = j(seg(it) == s);
        [~,o] = sortrows([tnum(js) snum_all(js)]);
        rank_all(js(o)) = 1:numel(js);
    end
end
fprintf('\nejection side switches: %d, in %d flies (ejection count restarts at each)\n',n_switch,n_fly_switch);
dye_all = [A.atp_amp]'; early_all = [A.early]'; late_all = [A.late]';

ord_res = table('Size',[nC 12],'VariableTypes',[{'string'},repmat({'double'},1,11)], ...
    'VariableNames',{'cond','n_fly','n_stim','max_rank','b_rank_early','p_rank_early','b_rank_late','p_rank_late', ...
                     'b_snum_early','p_snum_early','b_dye_vs_rank','p_dye_vs_rank'});
ord_res{:,2:end} = nan(nC,11);
bin_edges = [1 2 3 4 5 6 11 21 51 inf]; bin_lbl = {'1','2','3','4','5','6-10','11-20','21-50','51+'};
ord_bins = cell(nC,1);
for c = 1:nC
    if is_diff(c)
        idx = ej & strcmp(cond,C{c,2}) & ismember(geno,{'lpsp','empty'});
    else
        idx = ej & strcmp(cond,C{c,2}) & strcmp(geno,'');
    end
    T = table(early_all(idx),late_all(idx),dye_all(idx),log2(rank_all(idx)),snum_all(idx),categorical(fly(idx)), ...
        double(strcmp(geno(idx),'lpsp')),'VariableNames',{'early','late','dye','lrank','snum','fly','is_lpsp'});
    ord_res.cond(c) = C{c,1}; ord_res.n_fly(c) = numel(unique(T.fly)); ord_res.n_stim(c) = height(T);
    ord_res.max_rank(c) = max(2.^T.lrank);
    if height(T) < 6 || ord_res.n_fly(c) < 2 || range(T.lrank) == 0; continue; end
    if is_diff(c); term = 'lrank:is_lpsp'; rhs = 'dye + lrank*is_lpsp + (1|fly)';
    else;          term = 'lrank';         rhs = 'dye + lrank + (1|fly)'; end
    [ord_res.b_rank_early(c),ord_res.p_rank_early(c)] = lme_term(T,['early ~ ' rhs],term);
    [ord_res.b_rank_late(c), ord_res.p_rank_late(c)]  = lme_term(T,['late ~ '  rhs],term);
    [ord_res.b_snum_early(c),ord_res.p_snum_early(c)] = lme_term(T,['early ~ ' strrep(rhs,'lrank','snum')],strrep(term,'lrank','snum'));
    [ord_res.b_dye_vs_rank(c),ord_res.p_dye_vs_rank(c)] = lme_term(T,'dye ~ lrank + (1|fly)','lrank');

    % binned, dye-adjusted response (pooled dye slope removed) per genotype, mean +/- SEM across flies
    bd = T.dye - mean(T.dye); b = [ones(height(T),1) bd] \ T.early; adj = T.early - b(2)*bd;
    bins = discretize(2.^T.lrank,bin_edges);
    grp = {''}; if is_diff(c); grp = {'lpsp','empty'}; end
    out = struct();
    for g = 1:numel(grp)
        if is_diff(c); gi = T.is_lpsp == strcmp(grp{g},'lpsp'); else; gi = true(height(T),1); end
        [out(g).adj_m,out(g).adj_se,out(g).dye_m,out(g).dye_se] = deal(nan(1,numel(bin_lbl)));
        for bb = 1:numel(bin_lbl)
            jj = gi & bins==bb; fl_b = unique(T.fly(jj));
            if isempty(fl_b); continue; end
            pa = arrayfun(@(f) mean(adj(jj & T.fly==f)),fl_b); pd = arrayfun(@(f) mean(T.dye(jj & T.fly==f)),fl_b);
            out(g).adj_m(bb) = mean(pa); out(g).dye_m(bb) = mean(pd);
            if numel(fl_b) > 1; out(g).adj_se(bb) = std(pa)/sqrt(numel(fl_b)); out(g).dye_se(bb) = std(pd)/sqrt(numel(fl_b)); end
        end
    end
    ord_bins{c} = out;
end
writetable(ord_res,fullfile(tab_dir,'stim_order_effects.csv'));
fprintf('\n=== stim order (log2 cumulative ejection rank within fly), dye-adjusted ===\n');
fprintf('%-34s %4s %5s %6s | %-16s %-16s | %-16s | %-16s\n','condition','fly','stim','maxrk','early/doubling','late/doubling','early/stim# in trial','dye/doubling');
for c = 1:nC
    fprintf('%-34s %4d %5d %6d | %+.3f (p=%.2g) | %+.3f (p=%.2g) | %+.3f (p=%.2g) | %+.2f (p=%.2g)\n',ord_res.cond(c),ord_res.n_fly(c),ord_res.n_stim(c),ord_res.max_rank(c), ...
        ord_res.b_rank_early(c),ord_res.p_rank_early(c),ord_res.b_rank_late(c),ord_res.p_rank_late(c), ...
        ord_res.b_snum_early(c),ord_res.p_snum_early(c),ord_res.b_dye_vs_rank(c),ord_res.p_dye_vs_rank(c));
end

%% figure 7: dye and dye-adjusted response vs stim order
figure(7); clf; set(gcf,'Position',[10 20 2800 700],'Color','w')
tl = tiledlayout(2,nC,'TileSpacing','compact','Padding','compact');
for c = 1:nC
    out = ord_bins{c};
    for row = 1:2
        ax = nexttile(tl,c+(row-1)*nC); hold(ax,'on')
        if isempty(out); axis(ax,'off'); if row==1; title(ax,{C{c,1},'too few stims'},'FontSize',7); end; continue; end
        if is_diff(c); cols = {cl,ce}; else; cols = {cg}; if contains(C{c,1},'LPsP DA'); cols = {cg}; end; end
        for g = 1:numel(out)
            if row == 1; m = out(g).dye_m; se = out(g).dye_se; colr = cols{g}*0.6 + [cr]*0.4*(~is_diff(c));
            else;        m = out(g).adj_m; se = out(g).adj_se; colr = cols{g}; end
            ok = ~isnan(m); x = find(ok);
            errorbar(ax,x,m(ok),se(ok),'-o','Color',colr,'MarkerFaceColor',colr,'MarkerSize',3,'LineWidth',1.2,'CapSize',2)
        end
        xlim(ax,[0.5 numel(bin_lbl)+0.5]); set(ax,'XTick',1:numel(bin_lbl),'XTickLabel',bin_lbl,'FontSize',6.5)
        if row == 1
            title(ax,sprintf('%s\n%d flies, %d stims\ndye/doubling %+.2f (p=%.2g)',C{c,1},ord_res.n_fly(c),ord_res.n_stim(c), ...
                ord_res.b_dye_vs_rank(c),ord_res.p_dye_vs_rank(c)),'FontSize',7)
            if c==1; ylabel(ax,'red (dye) amp, 0-3 s'); end
        else
            yline(ax,0,':k')
            if is_diff(c); lab = '\Delta(lpsp-empty)'; else; lab = 'resp'; end
            title(ax,sprintf('%s/doubling\n0-3s %+.3f (p=%.2g)\n3-15s %+.3f (p=%.2g)',lab,ord_res.b_rank_early(c),ord_res.p_rank_early(c), ...
                ord_res.b_rank_late(c),ord_res.p_rank_late(c)),'FontSize',7)
            xlabel(ax,'ejection # (same fly, same side)')
            if c==1; ylabel(ax,'dye-adjusted green \DeltaF/F, 0-3 s'); end
        end
    end
end
title(tl,{'Do later ejections have a bigger effect? Top: dye per ejection (is the pipette drifting?). Bottom: green response with the pooled dye effect removed, by cumulative ejection number within the fly, restarting at 1 whenever the ejection side switches', ...
    'mean \pm SEM across flies per bin; slopes from mixed models (1|fly) per doubling of ejection number, dye as covariate; LPsP-stim columns: orange lpsp, blue empty, slope = lpsp x order interaction'})
exportgraphics(gcf,fullfile(out_dir,'fig07_stim_order.pdf'))

%% figure 6: ATP vs dopamine evidence per fly
FFc = FF(FF.call~="none",:);
is_p2 = FFc.geno~="";
mk_col = containers.Map({'ATP','dopamine','uncertain'},{[0.90 0.45 0.05],[0.45 0.20 0.70],[0.5 0.5 0.5]});
conf_sz = containers.Map({'confirmed','high','medium','low'},{45,60,30,12});
figure(6); clf; set(gcf,'Position',[10 20 2000 760],'Color','w')
tl = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
ax = nexttile(tl,1); hold(ax,'on')
patch(ax,[8 1e3 1e3 8],[30 30 3e3 3e3],[1 .9 .8],'EdgeColor','none','FaceAlpha',.5)
patch(ax,[1e-2 5 5 1e-2],[1e-1 1e-1 30 30],[.9 .85 1],'EdgeColor','none','FaceAlpha',.5)
for k = 1:height(FFc)
    mrk = 'o'; if ~is_p2(k); mrk = 's'; end
    scatter(ax,max(FFc.red_amp_ej(k),1e-2),max(FFc.red_z_ej(k),1e-1),conf_sz(char(FFc.conf(k))),mk_col(char(FFc.call(k))),mrk,'filled','MarkerFaceAlpha',.6)
end
set(ax,'XScale','log','YScale','log'); xlabel(ax,'red amplitude on ejection stims (median, \Sigma\DeltaF/F)'); ylabel(ax,'red z on ejection stims (median)')
title(ax,{'dye SNR vote (orange box = ATP, purple box = dopamine)','circles = p2x2 datasets, squares = dopamine\_ionto datasets'})
ax = nexttile(tl,2); hold(ax,'on')
xg = FFc.green_1s; isepg = contains(FFc.cond,'EPG') & ~contains(FFc.cond,'GRAB');
patch(ax,[1e-1 3e3 3e3 1e-1],[0.05 0.05 0.6 0.6],[1 .9 .8],'EdgeColor','none','FaceAlpha',.5)
patch(ax,[1e-1 3e3 3e3 1e-1],[0.015 0.015 0.05 0.05],[1 .96 .9],'EdgeColor','none','FaceAlpha',.5)
patch(ax,[1e-1 3e3 3e3 1e-1],[-0.3 -0.3 0 0],[.9 .85 1],'EdgeColor','none','FaceAlpha',.5)
for k = find(isepg)'
    mrk = 'o'; if ~is_p2(k); mrk = 's'; end
    scatter(ax,max(FFc.red_z_ej(k),1e-1),xg(k),conf_sz(char(FFc.conf(k))),mk_col(char(FFc.call(k))),mrk,'filled','MarkerFaceAlpha',.6)
end
set(ax,'XScale','log'); ylim(ax,[-0.3 0.6]); xlabel(ax,'red z on ejection stims'); ylabel(ax,'EPG green \DeltaF/F at ~1 s')
title(ax,{'EPG-green vote (EPG-GCaMP flies only)','orange = strong ATP (\geq0.05), pale = weak ATP, purple = dopamine (\leq0)'})
ax = nexttile(tl,3); axis(ax,'off')
text(ax,0,1,{'\bfColor\rm = call:  orange ATP, purple dopamine, gray uncertain', ...
    '\bfSize\rm = confidence:  large high, medium, small low','', ...
    'confirmed: walking p2x2 0.5 s flies, ATP per Pablo (low ejection strength)', ...
    'high: \geq2 votes agree, none conflict','medium: one strong vote, or agreeing weak votes', ...
    'low: votes conflict','', ...
    'p2x2 flies called dopamine (high) get their own column in figures 1-5;', ...
    'medium/low dopamine calls stay in their ATP group (see title counts).'},'VerticalAlignment','top','FontSize',10)
ax = nexttile(tl,[1 3]); hold(ax,'on')
[ug2,~,ig2] = unique(grp_cond(FF.call~="none")); calls = FF.call(FF.call~="none"); confs = FF.conf(FF.call~="none");
cats = {'ATP','confirmed';'ATP','high';'ATP','medium';'ATP','low';'uncertain','low';'dopamine','low';'dopamine','medium';'dopamine','high'};
cc  = [0.60 0.25 0.00; 0.90 0.45 0.05; 0.97 0.65 0.35; 1 .85 .7; .75 .75 .75; .85 .78 .95; .65 .5 .85; 0.45 0.20 0.70];
M = zeros(numel(ug2),size(cats,1));
for k = 1:size(cats,1); for g = 1:numel(ug2); M(g,k) = sum(ig2==g & calls==cats{k,1} & confs==cats{k,2}); end; end
b = barh(ax,M,'stacked','EdgeColor','none'); for k = 1:numel(b); b(k).FaceColor = cc(k,:); end
set(ax,'YTick',1:numel(ug2),'YTickLabel',ug2,'TickLabelInterpreter','none','FontSize',7,'YDir','reverse')
xlabel(ax,'flies'); legend(ax,strcat(cats(:,1),{' / '},cats(:,2)),'Location','eastoutside','FontSize',8)
title(ax,'substance call per fly, by dataset group')
title(tl,'Which substance was ejected? Per-fly evidence (dye SNR, EPG green response, file-name label)')
exportgraphics(gcf,fullfile(out_dir,'fig06_atp_vs_dopamine_calls.pdf'))

fprintf('figures written to %s\n',out_dir);

%% functions
function k = trial_sort_key(p)
    % YYYYMMDD*1e4 + trial number from a "...\YYYYMMDD-N_..." trial path, NaN if absent
    t = regexp(p,'(\d{8})-(\d+)','tokens','once');
    if isempty(t); k = nan; else; k = str2double(t{1})*1e4 + str2double(t{2}); end
end

function [b,p] = lme_term(T,formula,term)
    % fixed-effect estimate and p of one term of a linear mixed model; NaN if it can't be fit
    b = nan; p = nan;
    try
        m = fitlme(T,formula);
        k = strcmp(m.Coefficients.Name,term);
        if any(k); b = m.Coefficients.Estimate(k); p = m.Coefficients.pValue(k); end
    catch
    end
end

function s = subst_line(flies,call_of,conf_of)
    % "ATP 12 (h10 m2 l0) | DA 0 | unc 1" for the flies in one panel
    c = cellfun(@(f) call_of(f),flies,'UniformOutput',false);
    q = cellfun(@(f) conf_of(f),flies,'UniformOutput',false);
    n = @(x,y) sum(strcmp(c,x) & strcmp(q,y));
    s = sprintf('ATP %d (c%d h%d m%d l%d)\nDA %d (h%d m%d l%d), unc %d', ...
        sum(strcmp(c,'ATP')),n('ATP','confirmed'),n('ATP','high'),n('ATP','medium'),n('ATP','low'), ...
        sum(strcmp(c,'dopamine')),n('dopamine','high'),n('dopamine','medium'),n('dopamine','low'),sum(strcmp(c,'uncertain')));
end

function S = empty_stims()
    S = struct('ds',{},'fly',{},'cond',{},'geno',{},'stimlen',{},'z',{},'ej',{},'isi',{},'early',{},'late',{}, ...
        'ipsi',{},'contra',{},'atp_amp',{},'g',{},'r',{},'hm',{},'right',{});
end

function tf = has_stims(tr)
    tf = isfield(tr,'ft') && isstruct(tr.ft) && isfield(tr.ft,'stims') && ~isempty(tr.ft.stims) && max(tr.ft.stims) > 0;
end

function tf = has_atp(tr)
    tf = isfield(tr,'atp') && isstruct(tr.atp) && isfield(tr.atp,'d') && ~isempty(tr.atp.d) && ~all(isnan(tr.atp.d(:)));
end

function [on_t,slen] = stim_onsets(ft)
    st = ft.stims(:); s = st > 0.5*max(st);
    on = find(diff([false;s])==1); off = find(diff([s;false])==-1);
    xf = ft.xf(:); on_t = xf(on); slen = xf(off) - xf(on);
end

function on_t = atp_detect_onsets(tr)
    % stim onsets for trials with no ft.stims: peaks of the red channel
    % (same smoothing/prominence as lpsp_p2x2_script_minimal.m), then walk
    % back from each peak to where the signal first left baseline
    xb  = tr.ft.xb(:)'; n = min(numel(xb),size(tr.atp.d,2)); xb = xb(1:n);
    a   = sum(tr.atp.d(:,1:n),1);
    fr  = 1/median(diff(xb));
    fwd = smoothdata(a,'movmean',[0,20]);
    [~,loc] = findpeaks(fwd,'MinPeakProminence',9,'MinPeakDistance',round(20*fr));
    a3 = smoothdata(a,'movmean',3);
    on_t = nan(numel(loc),1);
    for k = 1:numel(loc)
        p0 = loc(k); if p0 <= round(5*fr); continue; end
        base = median(a3(max(1,p0-round(5*fr)):p0-1));
        [mx,pk] = max(a3(p0:min(n,p0+20))); pk = pk + p0 - 1;
        thr = base + 0.2*(mx - base);
        j = pk; while j > 1 && a3(j) > thr; j = j-1; end
        on_t(k) = xb(j);
    end
    on_t = on_t(~isnan(on_t));
end

function s = extract_stims(tr,on_t,slen,t_common,win,score_ejection)
    s = empty_stims();
    n  = size(tr.im.d,2);
    xb = tr.ft.xb(:)'; xb = xb(1:min(end,n));
    A  = tr.atp.d(:,1:numel(xb)); D = tr.im.d(:,1:numel(xb));
    for j = 1:numel(on_t)
        t0 = on_t(j);
        a  = interp1(xb-t0,A',t_common,'linear',nan)';
        d  = interp1(xb-t0,D',t_common,'linear',nan)';
        asum = sum(a,1);
        bm = median(asum(win.base),'omitnan'); bs = std(asum(win.base),'omitnan');
        pa = median(asum(win.pk),'omitnan');
        ev = mean(a(:,win.pk),2,'omitnan') - mean(a(:,win.base),2,'omitnan'); ev(ev<0) = 0;
        is_right = sum(ev(1:16))/sum(ev) > 0.5;
        d  = d - mean(d(:,win.base),2,'omitnan');
        if ~is_right; d = d([17:32 1:16],:); end % ejection side -> rows 1-16
        g  = mean(d,1,'omitnan');
        hi = mean(d(1:16,:),1,'omitnan'); hc = mean(d(17:32,:),1,'omitnan');
        if j > 1; isi = t0 - on_t(j-1); else; isi = inf; end
        if score_ejection; is_ej = (pa-bm) > 3*bs; else; is_ej = true; end
        if is_ej; hm = single(d); else; hm = []; end % heatmaps only kept for ejections (memory)
        s(end+1) = struct('ds','','fly','','cond','','geno','','stimlen',slen(j),'z',(pa-bm)/bs, ...
            'ej',is_ej,'isi',isi,'early',mean(g(win.pk)),'late',mean(g(win.late)), ...
            'ipsi',mean(hi(win.pk)),'contra',mean(hc(win.pk)),'atp_amp',mean(asum(win.pk))-bm, ...
            'g',g,'r',asum-bm,'hm',hm,'right',is_right); %#ok<AGROW>
    end
end

function [geno,ind,why] = p2x2_labels(tr)
    % genotype from the fly + trial folders, indicator from the trial
    % folder, and why (if at all) the trial is excluded -- see header
    p = strsplit(tr,'\'); p(cellfun(@isempty,p)) = [];
    two = lower(strjoin(p(end-1:end),'\')); tname = lower(p{end});
    if contains(two,'empty') || contains(two,'control')
        geno = 'empty';
    elseif contains(two,'lpsp') || contains(tname,'th_p2x2') || contains(tname,'blind')
        geno = 'lpsp'; % blinded + untagged fly folder = lpsp (only controls were tagged)
    else
        geno = '';
    end
    if contains(tname,'syt8s'); ind = 'syt8s';
    elseif contains(tname,'8m');  ind = '8m';
    else; ind = '7f'; end % 2024 trials (named syt7f/7f, or blinded)
    low = lower(tr);
    if contains(low,'open loop');         why = 'open loop';
    elseif contains(low,'exclusion');     why = 'exclusions folder';
    elseif contains(low,'dont use');      why = 'pilot (dont use)';
    elseif startsWith(p{end},'2024') && contains(tname,'_low'); why = '2024 low ATP';
    elseif isempty(geno);                 why = 'no genotype';
    else;                                 why = '';
    end
end

function c = p2x2_stratum(tr,ind,slen,from_stims)
    % experiment (dataset root) x indicator x stim length
    p = strsplit(tr,'\'); root = p{3};
    names = containers.Map({'lpsp_p2x2','lpsp_p2x2_redo','lpsp_p2x2_reredo','lpsp_p2x2_syt8s','lpsp_p2x2_walking'}, ...
                           {'2024','redo','reredo','syt8s','walking'});
    ex = names(root);
    if ~from_stims; c = sprintf('%s EPG %s | red-detected',ex,ind); return; end
    L = [0.5 2 15 30]; [~,k] = min(abs(L - slen));
    c = sprintf('%s EPG %s | %g s',ex,ind,L(k));
end

function key = p2x2_trial_key(tr)
    % root + trial folder name, normalized ("20260714-4-epg" == "20260714-4_epg");
    % stable across the folder reorganizations that moved trials between .mat pulls
    p = strsplit(tr,'\'); p(cellfun(@isempty,p)) = [];
    key = [p{3} '|' regexprep(lower(p{end}),'^(\d{8}-\d+)[-_]','$1_')];
end

function fid = p2x2_fly_id(tr)
    % root + date + fly NUMBER, so "fly 1" and "fly 1 _ control" are one fly;
    % trials with no fly folder fall back to their own trial name
    p = strsplit(tr,'\'); p(cellfun(@isempty,p)) = [];
    d  = regexp(p{end},'^(\d{8})','tokens','once');
    fi = find(~cellfun(@isempty,regexpi(p,'^fly\s*\d+','once')),1,'last');
    if isempty(fi)
        fid = [p{3} '|' p{end}];
    else
        n = regexpi(p{fi},'^fly\s*(\d+)','tokens','once');
        fid = sprintf('%s|%s|fly %s',p{3},d{1},n{1});
    end
end

function [m,se,nf,M] = fly_mean(X,idx,fly)
    f = unique(fly(idx)); nf = numel(f); M = nan(nf,size(X,2));
    for k = 1:nf; M(k,:) = mean(X(idx & strcmp(fly,f{k}),:),1,'omitnan'); end
    if nf==0; m = nan(1,size(X,2)); se = m; return; end
    m = mean(M,1,'omitnan'); se = std(M,0,1,'omitnan')/sqrt(nf);
end

function [v,se,p,nf,dots] = cond_effect(c,metric,A,C,is_diff,ej,cond,geno,fly)
    % per-fly values of one scalar metric; for a diff condition, lpsp minus empty
    x = [A.(metric)]';
    if ~is_diff(c)
        idx = ej & strcmp(cond,C{c,2}) & strcmp(geno,'');
        [v,se,nf1,dots] = fly_mean(x,idx,fly); nf = [nf1 0];
        if nf1 > 1; [~,p] = ttest(dots); else; p = nan; end
        if nf1 < 2; se = nan; end
    else
        iL = ej & strcmp(cond,C{c,2}) & strcmp(geno,'lpsp');
        iE = ej & strcmp(cond,C{c,2}) & strcmp(geno,'empty');
        [mL,seL,nL,dL] = fly_mean(x,iL,fly); [mE,seE,nE,dE] = fly_mean(x,iE,fly);
        v = mL - mE; se = sqrt(seL^2 + seE^2); nf = [nL nE]; dots = dL - mE;
        if nL > 1 && nE > 1; [~,p] = ttest2(dL,dE); else; p = nan; end
    end
end

function [M,nf] = fly_heatmap(A,idx,fly)
    idx = find(idx); f = unique(fly(idx)); nf = numel(f);
    M = nan(32,numel(A(1).g),nf);
    for k = 1:nf
        j = idx(strcmp(fly(idx),f{k}));
        M(:,:,k) = mean(cat(3,A(j).hm),3,'omitnan');
    end
    M = mean(M,3,'omitnan'); % omitnan: stims near trial start leave NaN pre-stim samples
end

function shade(ax,t,m,se,col,lw)
    if all(isnan(m)); return; end
    se(isnan(se)) = 0; ok = ~isnan(m);
    fill(ax,[t(ok) fliplr(t(ok))],[m(ok)+se(ok) fliplr(m(ok)-se(ok))],col,'FaceAlpha',.2,'EdgeColor','none')
    plot(ax,t,m,'Color',col,'LineWidth',lw)
end

function [ct,dr] = da_trial_labels(meta,ds)
    tn = regexp(meta,'\\(\d{8}-\d+_[^\\]*)\\registration','tokens','once'); tn = lower(tn{1});
    if contains(tn,'lpsp'); ct = 'LPsP (syt8m)';
    elseif contains(tn,'grab'); ct = 'EPG>GRAB-DA2m';
    elseif contains(tn,'syt8s'); ct = 'EPG (syt8s), empty>P2X2';
    elseif contains(tn,'7f'); ct = 'EPG (7f)';
    else; ct = 'EPG (syt8m)'; end
    has_da  = ~isempty(regexp(tn,'(_da_|_da$|dopamine)','once')) || contains(ds,'walking');
    has_ach = contains(tn,'ach');
    if has_da && has_ach; dr = 'DA+ACh'; elseif has_ach; dr = 'ACh'; elseif has_da; dr = 'DA'; else; dr = 'untagged'; end
end

function fid = fly_folder(meta)
    % date folder + "fly N" (or just the date folder when there's no fly level)
    p = strsplit(meta,'\'); p(cellfun(@isempty,p)) = [];
    fi = find(startsWith(p,'fly'),1);
    if isempty(fi); fi = find(~cellfun(@isempty,regexp(p,'^\d{8}$')),1); end
    if isempty(fi); fid = meta; return; end
    fid = strjoin(p(fi-1:fi),'\');
end
