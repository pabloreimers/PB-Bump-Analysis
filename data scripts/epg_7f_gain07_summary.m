%% epg_7f_gain07_summary
% One-line-per-trial table of .data\epg_7f_20261007.mat (built by
% epg_7f_gain07_batch.m): fly, trial, condition, gain, mask source, bump
% tracking metrics, walking fraction. Written to
% ugly_figures\exports\epg_7f_gain07\summary.csv and printed.

repoRootDir = fullfile(fileparts(mfilename('fullpath')), '..');
S = load(fullfile(repoRootDir, '.data', 'epg_7f_20261007.mat'), 'all_data');
all_data = S.all_data; clear S
n = numel(all_data);
T = table('Size', [n 14], 'VariableTypes', {'double','double','string','string','double','double','string','double','double','double','double','double','double','string'}, ...
    'VariableNames', {'fly','trial','flyDir','trialFolder','tif','gain','condition','fps','nVol','pct_walking','offset_circstd_deg','lag_s','pct_vol_used','mask_source'});
for i = 1:n
    d = all_data(i); [~, tf] = fileparts(d.meta.trialDir);
    T.fly(i) = d.meta.fly_num; T.trial(i) = d.meta.trial_num;
    T.flyDir(i) = strrep(d.meta.flyDir, 'Z:\pablo\gain_change\', ''); T.trialFolder(i) = tf; T.tif(i) = d.meta.tifIdx;
    T.gain(i) = round(abs(d.ft.gain_empirical), 2);
    if d.ft.dark, T.condition(i) = "dark"; else, T.condition(i) = "bar"; end
    T.fps(i) = round(d.meta.fps, 2); T.nVol(i) = d.meta.nVol;
    T.pct_walking(i) = round(100 * mean(abs(d.ft.f_speed) > 0.5 | abs(d.ft.r_speed) > 0.3)); % same spirit as the lab's walk_idx
    T.offset_circstd_deg(i) = round(d.bump.offset_circstd_deg); T.lag_s(i) = round(d.bump.lag_sec, 2);
    T.pct_vol_used(i) = round(100 * mean(d.bump.ok)); T.mask_source(i) = string(d.meta.mask_source);
end
disp(T)
outCsv = fullfile(repoRootDir, 'ugly_figures', 'exports', 'epg_7f_gain07', 'summary.csv');
writetable(T, outCsv);
fprintf('%d trials, %d flies; %d bar / %d dark; mask sources: %s\nwrote %s\n', n, numel(unique(T.fly)), nnz(T.condition == "bar"), nnz(T.condition == "dark"), ...
    strjoin(cellstr(unique(T.mask_source)), ' | '), outCsv);
