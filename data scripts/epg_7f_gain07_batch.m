%% epg_7f_gain07_batch
% All EPG_7f flies in Z:\pablo\gain_change whose FIRST trial ran at closed-loop
% gain 0.7 (bar or dark), keeping each fly's trials in acquisition order up
% to (not including) the first trial at a different gain. Each kept trial is
% run through pb_gain_trial (raw tif -> registered -> auto mask -> glomeruli
% -> PVA + fictrac) and the results are stored as an all_data struct array in
% .data\epg_7f_20261007.mat, plus a per-trial QC png in
% ugly_figures\exports\epg_7f_gain07\.
%
% Fly selection follows data scripts/gain_change_fly_inventory.csv (2026-10-07):
% plain EPG_7f only (no _split, no _empty_kir, no GRAB, nothing in "to do"),
% and only the single 4-px bright bar / dark patterns -- flies whose first
% trial used a starfield (June-Sept 2023) are excluded.
% Gain per trial is the empirical bar-vs-heading slope from the fictrac/DAQ
% data (the gain is not logged anywhere); 0.7 means |gain| in [0.6 0.8].
%
% Run in batch (several hours; see the time estimate printed at the start):
%   matlab -batch "cd('<repo>'); addpath('data scripts'); epg_7f_gain07_batch"

%% 0. setup
repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(fullfile(repoRootDir, 'claude')); addpath(fullfile(repoRootDir, 'circ_stats')); addpath(repoRootDir);
base = 'Z:\pablo\gain_change';
outFile = fullfile(repoRootDir, '.data', 'epg_7f_20261007.mat');
qcDir = fullfile(repoRootDir, 'ugly_figures', 'exports', 'epg_7f_gain07');
if ~isfolder(qcDir), mkdir(qcDir); end
if ~isfolder(fileparts(outFile)), mkdir(fileparts(outFile)); end

gain_target = 0.7; gain_tol = 0.1;

% fly folders. Dates without "fly n" subfolders hold one fly each (trials
% 3:03-3:49 PM etc., see inventory); their raw *_gainchange duplicates are
% skipped in favour of the renamed copies.
flyDirs = {};
% 20230214 is excluded: different genotype (fly.csv: w+;UAS-7A;60d05GAL4, not the
% SS52578 split), no bridge-shaped signal in the mean image, never preprocessed.
dates = {'20221221','20230106','20230109','20230110','20230112','20230113','20230504','20230607','20230609', ...
         '20230614','20230615','20230616','20230719','20230721','20230731','20230801','20230802'};
for i = 1:numel(dates)
    sub = dir(fullfile(base, dates{i}, 'fly *'));
    sub = sub([sub.isdir]);
    if isempty(sub)
        flyDirs{end+1} = fullfile(base, dates{i}); %#ok<SAGROW>
    else
        for j = 1:numel(sub), flyDirs{end+1} = fullfile(sub(j).folder, sub(j).name); end %#ok<SAGROW>
    end
end
% 20230614 fly 1 is EPG_7f_empty_kir -> drop
flyDirs = flyDirs(~contains(flyDirs, fullfile('20230614', 'fly 1')));
% 20230616 fly 1: no bump in either trial (noise heatmap), fly did not turn,
% bar parked most of the time, no hand mask and the auto mask found no bridge
% -> unusable, drop (see ugly_figures/exports/epg_7f_gain07/fly17_* from the 2026-10-07 run)
flyDirs = flyDirs(~contains(flyDirs, fullfile('20230616', 'fly 1')));

%% 1. enumerate trials per fly, gate on gain, estimate run time
plan = struct('flyDir', {}, 'trialDir', {}, 'tifIdx', {}, 'expNum', {}, 'gain', {}, 'dark', {}, 'tifBytes', {});
fprintf('scanning %d flies...\n', numel(flyDirs));
for i = 1:numel(flyDirs)
    td = dir(fullfile(flyDirs{i}, '*_EPG_7f*')); td = td([td.isdir]);
    td = td(cellfun(@isempty, regexp({td.name}, '_gainchange$|_split|_empty_kir', 'once')));
    expNum = cellfun(@(n) str2double(regexp(n, '^\d{8}-(\d+)_', 'tokens', 'once')), {td.name});
    % one folder per experiment number: 20230109 has duplicate copies of the
    % same trial (e.g. -2_EPG_7f_.7 and -2_EPG_7f_dark_.7, identical tifs).
    % Prefer the folder whose name carries 'dark' if the trial's pattern is
    % background, else the one without.
    keep = true(size(td));
    for e = unique(expNum)
        ii = find(expNum == e);
        if numel(ii) > 1
            isDarkName = contains({td(ii).name}, 'dark');
            ts = fullfile(td(ii(1)).folder, td(ii(1)).name, 'csv', 'trialSettings.csv'); darkPat = false;
            if isfile(ts), try, T = readtable(ts, 'Delimiter', ','); darkPat = contains(char(string(T.patternPath(end))), 'background'); catch, end, end
            pick = find(isDarkName == darkPat, 1); if isempty(pick), pick = 1; end
            keep(ii(setdiff(1:numel(ii), pick))) = false;
        end
    end
    td = td(keep); expNum = expNum(keep);
    [~, ord] = sort(expNum); td = td(ord); expNum = expNum(ord);
    [~, flyName] = fileparts(flyDirs{i}); if ~startsWith(flyName, 'fly'), flyName = ''; end
    seq = {}; stop = false;
    for j = 1:numel(td)
        tDir = fullfile(td(j).folder, td(j).name);
        tifs = dir(fullfile(tDir, '*_trial_*_*.tif'));
        tifTrial = cellfun(@(n) str2double(regexp(n, '_trial_(\d+)_', 'tokens', 'once')), {tifs.name});
        [tifTrial, o] = sort(tifTrial); tifs = tifs(o);
        for k = 1:numel(tifs)
            try
                ft = pb_gain_trial_ft(tDir, tifTrial(k)); % quick gain/dark read, no imaging
            catch ME
                fprintf('  %s trial %d: cannot read fictrac/DAQ (%s) -> stop here for this fly\n', td(j).name, tifTrial(k), ME.message);
                stop = true; break
            end
            is07 = abs(abs(ft.gain) - gain_target) <= gain_tol;
            % only the plain single bright bar (0003_4px_brightbar) or dark
            % (0001_background); any starfield pattern (0006 bar+starfield,
            % 0008 starfield_bar) ends the fly's sequence -- flies that START
            % with a starfield contribute nothing.
            plainPattern = contains(ft.pattern, '4px_brightbar.mat') || contains(ft.pattern, 'background');
            seq{end+1} = sprintf('%s%.2f%s', char('B' + 2*ft.dark), abs(ft.gain), repmat('*', 1, ~plainPattern)); %#ok<SAGROW>
            if ~is07 || ~plainPattern, stop = true; break; end
            plan(end+1) = struct('flyDir', flyDirs{i}, 'trialDir', tDir, 'tifIdx', tifTrial(k), 'expNum', expNum(j), ...
                'gain', ft.gain, 'dark', ft.dark, 'tifBytes', tifs(k).bytes); %#ok<SAGROW>
        end
        if stop, break; end
    end
    nKept = nnz(strcmp({plan.flyDir}, flyDirs{i}));
    fprintf('  %-45s %s -> keep %d trial(s)\n', strrep(flyDirs{i}, [base filesep], ''), strjoin(seq, ' '), nKept); % * = starfield pattern
end
if isempty(plan), error('epg_7f_gain07_batch:empty', 'no trials selected'); end
% ~ 1 s per 50 MB of tif for the read + ~3.5 min registration per 2700 volumes on this machine
estMin = sum([plan.tifBytes]) / 2.2e9 * 5.5;
fprintf('\n%d trials from %d flies selected (%.1f GB of tif). Estimated run time ~%.0f min (~%.1f h).\n\n', ...
    numel(plan), numel(unique({plan.flyDir})), sum([plan.tifBytes])/1e9, estMin, estMin/60);
if exist('gain07_dry_run', 'var') && gain07_dry_run
    fprintf('dry run: stopping before processing.\n'); return
end

%% 2. process, one fly at a time
% Mask priority (pb_gain_trial): 1) the hand-drawn mask.mat in the trial
% folder when it exists (15 of the 18 flies have one per trial); 2) else the
% automatic CV-image mask, computed on the fly's reference trial (first kept
% bar trial, since the CV image needs a visible moving bump -- in the dark the
% bump can be too weak and the mask comes out as a blob) and shifted onto the
% fly's other trials by aligning mean images. out.meta.mask_source says which.
all_data = struct('ft', {}, 'im', {}, 'bump', {}, 'meta', {});
flyList = unique({plan.flyDir}, 'stable');
t0 = tic; nDone = 0;
for f = 1:numel(flyList)
    ip = find(strcmp({plan.flyDir}, flyList{f})); % in acquisition order
    iRef = ip(find(~[plan(ip).dark], 1));
    if isempty(iRef)
        % dark-only fly: try every kept trial's own mask, keep the one that
        % passed sanity on its own CV image (else the widest)
        widths = nan(size(ip)); own = false(size(ip));
        for q = 1:numel(ip)
            try
                o = pb_gain_trial(plan(ip(q)).trialDir, plan(ip(q)).tifIdx, 'verbose', false);
                widths(q) = o.meta.mask_width_frac; own(q) = strcmp(o.meta.mask_source, 'own CV mask');
            catch, end
        end
        score = widths + own; [~, qBest] = max(score); iRef = ip(qBest);
        fprintf('fly %d (%s) has no bar trial kept: reference = trial %d (mask width %.0f%% of frame, own CV mask: %d)\n', f, strrep(flyList{f}, [base filesep], ''), qBest, 100*widths(qBest), own(qBest));
    end
    order = [iRef, setdiff(ip, iRef, 'stable')];
    refMask = [];
    for i = order
        nDone = nDone + 1;
        [~, tName] = fileparts(plan(i).trialDir);
        fprintf('[%d/%d] fly %d: %s trial %d (elapsed %.0f min)%s\n', nDone, numel(plan), f, tName, plan(i).tifIdx, toc(t0)/60, repmat(' [reference]', 1, i == iRef));
        try
            if i == iRef
                out = pb_gain_trial(plan(i).trialDir, plan(i).tifIdx);
                refMask = struct('mask', out.im.mask, 'maskCore', out.im.maskCore, 'projMean', out.im.projMean);
            else
                out = pb_gain_trial(plan(i).trialDir, plan(i).tifIdx, 'refMask', refMask);
            end
        catch ME
            fprintf('  FAILED: %s\n', ME.message); continue
        end
        out.meta.fly_num = f; out.meta.trial_num = find(ip == i); out.meta.flyDir = plan(i).flyDir; out.meta.expNum = plan(i).expNum;
        out.meta.gain_nominal = gain_target; out.meta.is_reference_trial = (i == iRef);
        all_data(end+1) = out; %#ok<SAGROW>
        local_qc_png(out, fullfile(qcDir, sprintf('fly%02d_t%d_%s_trial%d.png', f, out.meta.trial_num, tName, plan(i).tifIdx)));
        save(outFile, 'all_data', 'plan', '-v7.3'); % checkpoint after every trial
    end
end
% acquisition order within each fly (the reference was processed first);
% drop any duplicate (same fly, expNum, tif) that slipped through and renumber
keyFly = arrayfun(@(s) s.meta.fly_num, all_data); keyExp = arrayfun(@(s) s.meta.expNum, all_data); keyTif = arrayfun(@(s) s.meta.tifIdx, all_data);
[~, iu] = unique([keyFly(:) keyExp(:) keyTif(:)], 'rows', 'stable'); all_data = all_data(sort(iu));
keyFly = arrayfun(@(s) s.meta.fly_num, all_data); keyExp = arrayfun(@(s) s.meta.expNum, all_data); keyTif = arrayfun(@(s) s.meta.tifIdx, all_data);
[~, ord] = sortrows([keyFly(:) keyExp(:) keyTif(:)]); all_data = all_data(ord);
for i = 1:numel(all_data), all_data(i).meta.trial_num = nnz(keyFly(ord(1:i)) == all_data(i).meta.fly_num); end
save(outFile, 'all_data', 'plan', '-v7.3');
fprintf('\ndone: %d trials saved to %s (%.0f min)\n', numel(all_data), outFile, toc(t0)/60);

%% local functions
function ft = pb_gain_trial_ft(trialDir, tifIdx)
% gain + dark flag only (reads fictrac/DAQ, no imaging) -- same logic as pb_gain_trial's local_load_ft
ft.pattern = '';
ts = fullfile(trialDir, 'csv', 'trialSettings.csv');
if isfile(ts), try, T = readtable(ts, 'Delimiter', ','); ft.pattern = char(string(T.patternPath(end))); catch, end, end
ft.dark = contains(ft.pattern, 'background');
f = dir(fullfile(trialDir, '*ficTracData_DAQ.mat'));
if contains(ft.pattern, 'with_starfield'), f = []; end % stored cuePos is wrong for that pattern -> raw DAQ (see pb_gain_trial)
if numel(f) == 1
    S = load(fullfile(f(1).folder, f(1).name)); T = S.ftData_DAQ;
    row = find(T.trialNum == tifIdx, 1); if isempty(row), row = min(tifIdx, height(T)); end
    cue = unwrap(T.cuePos{row}(:)' / 192 * 2*pi); head = unwrap(T.intHD{row}(:)'); rate = 60;
else
    d = dir(fullfile(trialDir, sprintf('*daqData*trial_%03d.mat', tifIdx)));
    if numel(d) ~= 1, error('no ficTracData_DAQ.mat and no daqData trial %d', tifIdx); end
    S = load(fullfile(d(1).folder, d(1).name), 'trialData'); D = S.trialData;
    vn = D.Properties.VariableNames;
    cue = unwrap(D.(vn{find(contains(lower(vn), 'panel'), 1)})(:)' / 10 * 2*pi);
    head = unwrap(D.(vn{find(contains(lower(vn), 'yaw'), 1)})(:)' / 10 * 2*pi); rate = 10000;
end
ib = 1:rate:numel(cue); dc = diff(cue(ib)); dh = diff(head(ib)); mv = abs(dh) > 0.15;
if nnz(mv) >= 5, ft.gain = median(dc(mv) ./ dh(mv)); else, ft.gain = NaN; end
end

function local_qc_png(out, pngFile)
fig = figure('Visible', 'off', 'Position', [50 50 1600 1000], 'Color', 'w');
t = out.ft.xb; z = out.im.z; nC = size(z, 1);
ax1 = subplot(4,1,1); imagesc(t, 1:nC, z, 'AlphaData', double(~isnan(z))); clim([-3 3]); set(ax1, 'Color', [.7 .7 .7])
half = 128; colormap(ax1, [[linspace(0,1,half)' linspace(0,1,half)' ones(half,1)]; [ones(half,1) linspace(1,0,half)' linspace(1,0,half)']]);
cond = 'closed loop, bar'; if out.ft.dark, cond = 'closed loop, DARK'; end
ylabel('glomerulus'); title(sprintf('%s trial %d | %s | empirical gain %.2f | fictrac from %s', strrep(out.meta.trialDir, 'Z:\pablo\gain_change\', ''), ...
    out.meta.tifIdx, cond, out.ft.gain_empirical, out.ft.cue_src), 'Interpreter', 'none')
ax2 = subplot(4,1,2); hold on
for k = 1:size(out.bump.exclSegs, 1), patch(out.bump.exclSegs(k,[1 2 2 1]), [-180 -180 180 180], [.75 .75 .75], 'EdgeColor', 'none', 'FaceAlpha', .5); end
b = rad2deg(angle(exp(1i*out.bump.bar_sign*out.bump.cue_lag))); b(abs(diff([b b(end)])) > 180) = NaN; plot(t, b, '-', 'Color', [.3 .3 .3], 'LineWidth', 1.2)
m = rad2deg(out.im.mu); m(~out.bump.ok) = NaN; m(abs(diff([m m(end)])) > 180) = NaN; plot(t, m, 'b.', 'MarkerSize', 4)
ylim([-180 180]); ylabel('deg'); title(sprintf('bump (blue) vs bar (%s, lag %.2f s): offset circ std %.0f deg, %.0f%% volumes used', ...
    out.bump.signNote, out.bump.lag_sec, out.bump.offset_circstd_deg, 100*mean(out.bump.ok)))
ax3 = subplot(4,1,3); hold on
mu_u = out.bump.bar_sign * out.im.mu; mu_u(~out.bump.ok) = NaN; i = ~isnan(mu_u); u = nan(size(mu_u)); u(i) = unwrap(mu_u(i));
hs = out.ft.gain_empirical * (out.ft.heading - out.ft.heading(1)); hs = hs - median(hs(i) - u(i), 'omitnan');
plot(t, rad2deg(hs), '-', 'Color', [.85 .33 .1], 'LineWidth', 1.2); plot(t, rad2deg(u), 'b.', 'MarkerSize', 4)
ylabel('unwrapped deg'); title(sprintf('bump (blue) vs %.2f x heading (orange)', out.ft.gain_empirical))
ax4 = subplot(4,1,4); plot(t, out.im.rho, 'k'); yline(out.bump.rho_thresh, 'r:'); ylabel('rho'); xlabel('trial time (s)')
linkaxes([ax1 ax2 ax3 ax4], 'x'); xlim([t(1) t(end)])
% mask inset
axI = axes('Position', [0.72 0.90 0.27 0.09]); imagesc(axI, out.im.projCV); axis(axI, 'image', 'off'); colormap(axI, bone); hold(axI, 'on')
contour(axI, out.im.mask, [.5 .5], 'r', 'LineWidth', 1); title(axI, sprintf('%s (J=%.2f vs hand)', out.meta.mask_source, out.meta.mask_jaccard_vs_hand), 'FontSize', 8)
pb_export_png(fig, pngFile, 110); close(fig)
end
