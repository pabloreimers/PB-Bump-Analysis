%% dopamine_ionto_walking_claude
% Builds the all_data struct for the dopamine_ionto_walking dataset (ATP
% iontophoresis near the PB during closed-loop walking, EPG-syt8m/8m
% imaging of the bump, some trials with a simultaneous ATP/dye channel)
% from the raw session tree at base_dir, following the same build
% convention as dopamine_ionto_script.m (process_ft/process_im, one
% all_data(i) per trial with .im/.atp/.ft/.meta), then reports how many
% flies are in the dataset and shows one example trial per fly.
%
% Investigated directly before writing this script (all against the raw
% tree at base_dir, 82 trials total across 20260320, 20260325, and a
% piloting/ subtree with 3 more session dates):
%
% 1) REGISTRATION ALREADY DONE, MASKS ALREADY DRAWN: every one of the 82
%    trials already has registration\imagingData.mat and a trial-level
%    mask.mat (confirmed directly: 0/82 missing either). So unlike
%    dopamine_ionto_script.m's own mask section (which draws a mask via
%    roipoly for anything missing), the mask step below is just a
%    presence check -- it only prompts if a future trial is added without
%    one.
%
% 2) img{1}/img{2} ARE ALREADY Y x X x T, NOT Y x X x Z x T: confirmed
%    directly (size = [100,256,n_frames] for every trial checked, class
%    double/single, no singleton or Z dimension to sum over). This is
%    different from dopamine_ionto_script.m's own assumption
%    (imgData = squeeze(sum(regProduct,3)), squeeze(sum(img{1},3)) before
%    calling process_im) -- that sum(...,3) was collapsing a Z dimension
%    that doesn't exist in this dataset's saved registration output.
%    Applying it here would silently collapse the TIME dimension instead
%    (summing all frames into one), which process_im cannot use as a
%    movie. So below, img{1}/img{2} are passed to process_im directly,
%    with no extra sum() step.
%
% 3) TWO CHANNELS ON MOST, ONE CHANNEL ON A FEW EARLY TRIALS: confirmed
%    directly -- numel(img) is 2 (EPG + ATP/dye) for 80/82 trials, but only
%    1 (EPG only) for the first two trials of 20260320 fly 1 (the
%    iontophoresis electrode probably wasn't in place yet that session).
%    Same guard as dopamine_ionto_script.m: all_data(i).atp is only
%    populated when numel(img)>1.
%
% 4) FICTRAC/DAQ FIELD NAMES MATCH process_ft'S EXPECTATIONS EXACTLY:
%    confirmed directly against one trial's ftData_DAQ/ftData_dat --
%    ftData_DAQ has trialTime (duration), volClock (duration, length
%    matches n_frames exactly), stim (double, ATP-ejection command), and
%    cuePos (range [0,192], matching process_ft's /192*2*pi-pi
%    conversion); ftData_dat has its own trialTime (double, seconds) and
%    velFor. process_ft (copied unchanged from dopamine_ionto_script.m)
%    already expects exactly this. ft.stims = ftData_DAQ.stim{1} is added
%    below (dopamine_ionto_script.m left this commented out) since every
%    later analysis in this codebase (e.g. dopamine_ionto_redo_claude.m)
%    expects it.
%
% 5) ft.pattern (0003_4px_brightbar.mat on the trial checked) comes from
%    each trial's own csv\trialSettings.csv, same as dopamine_ionto_script.m.
%
% 6) do_build/do_save_data (below): the raw build (find trials -> load
%    every registered movie -> process_ft/process_im) is the slow part and
%    doesn't need to change just because a downstream parameter does, so
%    it's split from everything after it. do_build=false loads the
%    already-built .mat instead of re-touching the raw tree at all.
%    do_save_data is separate from do_build (only matters when do_build is
%    true) so re-running the slow build to sanity-check something doesn't
%    also force an equally slow overwrite of the saved .mat.
%
% 7) BUMP RE-ESTIMATION AND FULL-PB PLOTTING (below, after all_data is
%    available either way): mu/rho are re-derived from a 5-frame moving
%    average of each glomerulus's own raw f (im.f/atp.f -- untouched,
%    still the original per-trial raw trace), with heading (ft.cue)
%    smoothed the same way. This only touches already-extracted
%    per-glomerulus traces, so it's cheap enough to always re-run, even
%    when do_build is false. The example-trial figure's y-axis is left to
%    auto-scale to unwrap(im.alpha)'s actual range instead of being
%    clamped to [-pi,pi] -- that clamp was cutting the plot off after the
%    first 16 (of 32) glomeruli, i.e. showing only one PB hemisphere.

%% run configuration
do_build     = false; % true: rebuild all_data from the raw tree at base_dir (slow). false: load the .mat below instead
do_save_data = false; % only consulted when do_build is true -- whether to overwrite the saved .mat with the rebuilt data

data_dir    = fullfile('.data');
source_file = 'dopamine_ionto_walking_reg_20260831.mat';
f0_pct      = 7; % baseline percentile -- used both when building and when re-estimating mu/rho below

%% build all_data from the raw tree, or load the already-built .mat
if do_build
    base_dir = 'Z:\pablo\dopamine_ionto_walking\';
    assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

    all_files = dir(fullfile(base_dir,'**','registration','imagingData.mat'));
    all_files = natsortfiles(all_files);
    n_trials  = numel(all_files);
    fprintf('found %d trials under %s\n', n_trials, base_dir);

    % make sure every trial has a mask (no-op today -- see header point 1
    % -- but kept so a future trial added without one gets one drawn
    % interactively)
    for i = 1:n_trials
        trial_dir = fileparts(all_files(i).folder);
        mask_path = fullfile(trial_dir,'mask.mat');
        if isfile(mask_path); continue; end

        fprintf('drawing mask (missing): %s\n', trial_dir)
        clear img
        load(fullfile(all_files(i).folder,'imagingData.mat'),'img')
        imgData = img{1};

        top_pct = prctile(imgData,98,'all');
        bot_pct = prctile(imgData,5,'all');
        imgData(imgData>top_pct) = top_pct;
        imgData(imgData<bot_pct) = bot_pct;

        figure(1); clf; imagesc(mean(imgData,3)); colormap(bone); axis equal tight; drawnow;
        mask = roipoly(); %#ok<NASGU>
        save(mask_path,'mask')
    end

    % process and store all values
    ft_type = 'movmean'; % smoothing type for fictrac data
    ft_win  = 10;         % smoothing window for fictrac data (samples); gaussian windows have std = win/5
    im_type = {'movmean','movmean'}; % two smoothing steps for im data: summed z-stack, then estimated mu/rho
    im_win  = {1,1};
    n_centroid = 16;      % per hemisphere -- process_im doubles this to 32 total glomeruli

    all_data = struct();

    tic
    for i = 1:n_trials
        clear img

        trial_dir = fileparts(all_files(i).folder);
        fprintf('processing %d/%d: %s ', i, n_trials, trial_dir)

        load(fullfile(all_files(i).folder,'imagingData.mat'),'img')
        load(fullfile(all_files(i).folder,'imgData_reg.mat'),'imgData_reg')
        load(fullfile(trial_dir,'mask.mat'),'mask')

        tmp2 = dir(fullfile(trial_dir,'*ficTracData_DAQ.mat'));
        load(fullfile(tmp2(1).folder,tmp2(1).name))
        tmp2 = dir(fullfile(trial_dir,'*ficTracData_dat.mat'));
        load(fullfile(tmp2(1).folder,tmp2(1).name))
        tmp2 = readtable(fullfile(trial_dir,'csv','trialSettings.csv'));

        all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type); %#ok<NODEF>
        all_data(i).ft.pattern = tmp2.patternPath{1};
        all_data(i).ft.stims   = ftData_DAQ.stim{1};

        all_data(i).im = process_im(imgData_reg, im_win, im_type, mask, n_centroid, f0_pct);
        if numel(img) > 1
            all_data(i).atp = process_im(img{2}, im_win, im_type, mask, n_centroid, f0_pct);
        end
        all_data(i).meta = all_files(i).folder;

        fprintf('ETR: %.2f hours\n', toc/i * (n_trials-i) / 60 / 60)
    end

    if do_save_data
        if ~isfolder(data_dir); mkdir(data_dir); end
        out_file = fullfile(data_dir,source_file);
        save(out_file,'all_data','-v7.3')
        fprintf('saved %d trials to %s\n', n_trials, out_file);
    end
else
    tmp = load(fullfile(data_dir,source_file),'all_data');
    all_data = tmp.all_data(:)';
    n_trials = numel(all_data);
    fprintf('loaded %d trials from %s\n', n_trials, source_file);
end

%% re-estimate bump position (mu/rho) from a smoothed raw f, and smooth heading the same way
% raw f (im.f/atp.f) is left untouched -- only the derived d/z/mu/rho are
% recomputed, from a 5-frame moving average of that raw per-glomerulus
% trace (see mu_rho_from_f). ft.cue gets the same kind of moving-average
% smoothing (5 samples, on its own native fictrac clock xf), applied
% circularly (unwrap -> smooth -> re-wrap, same as process_ft's own cue
% smoothing) since it's a wrapped angle, not a linear signal.
bump_smooth_win = 5; % frames, imaging clock (xb)
cue_smooth_win  = 5; % samples, fictrac clock (xf)

for i = 1:n_trials
    r = mu_rho_from_f(all_data(i).im.f, bump_smooth_win, f0_pct);
    all_data(i).im.d   = r.d;
    all_data(i).im.z   = r.z;
    all_data(i).im.mu  = r.mu;
    all_data(i).im.rho = r.rho;

    if isfield(all_data(i),'atp') && ~isempty(all_data(i).atp)
        r = mu_rho_from_f(all_data(i).atp.f, bump_smooth_win, f0_pct);
        all_data(i).atp.d   = r.d;
        all_data(i).atp.z   = r.z;
        all_data(i).atp.mu  = r.mu;
        all_data(i).atp.rho = r.rho;
    end

    all_data(i).ft.cue = smooth_circular(all_data(i).ft.cue, cue_smooth_win);
end
fprintf('re-estimated mu/rho from a %d-frame moving average of raw f, and re-smoothed heading with a %d-sample moving average\n', bump_smooth_win, cue_smooth_win);

%% fly ID per trial (date folder + "fly N" folder), and fly count
fly_id = cell(1,n_trials);
for i = 1:n_trials
    fly_id{i} = trial_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
fly_num = fly_num(:)'; % row vector, matching other per-trial arrays
n_flies = numel(fly_list);

trials_per_fly = arrayfun(@(f) sum(fly_num==f), 1:n_flies);
fprintf('\n%d trials -> %d flies\n', n_trials, n_flies);
fprintf('\n=== trials per fly ===\n');
for f = 1:n_flies
    fprintf('  %-30s n=%d trials\n', fly_short_label(fly_list{f}), trials_per_fly(f));
end

%% plot: number of trials per fly (answers "how many flies are in the dataset")
figure(1); clf
set(gcf,'Name','flies in dataset','Position',[100,100,900,400])
bar(trials_per_fly,'FaceColor',[0.00,0.45,0.74])
set(gca,'XTick',1:n_flies,'XTickLabel',arrayfun(@(f) fly_short_label(fly_list{f}),1:n_flies,'UniformOutput',false),'XTickLabelRotation',45)
ylabel('number of trials')
title(sprintf('dopamine\\_ionto\\_walking: n=%d flies, %d trials total',n_flies,n_trials))
for f = 1:n_flies
    text(f,trials_per_fly(f)+0.1,num2str(trials_per_fly(f)),'HorizontalAlignment','center')
end

%% pick one example trial per fly: the trial with the highest mean bump
% vector strength (im.rho) among that fly's own trials, i.e. its cleanest/
% most confident bump estimate -- same "pick the clean one" convention
% lpsp_p2x2_walking_script.m uses for its own representative-trial figures.
example_trial = nan(1,n_flies);
for f = 1:n_flies
    idx = find(fly_num==f);
    mean_rho = arrayfun(@(i) mean(all_data(i).im.rho,'omitnan'), idx);
    [~,best] = max(mean_rho);
    example_trial(f) = idx(best);
end

%% figure: example trial per fly -- im.z heatmap across the PB, heading + mu overlaid
n_col_ex = ceil(sqrt(n_flies));
n_row_ex = ceil(n_flies/n_col_ex);

figure(2); clf
set(gcf,'Name','example trial per fly','Position',[50,50,320*n_col_ex,220*n_row_ex])
tl = tiledlayout(n_row_ex,n_col_ex,'TileSpacing','compact','Padding','compact');

for f = 1:n_flies
    i = example_trial(f);
    ax = nexttile(tl); hold(ax,'on')

    imagesc(ax,all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
    colormap(ax,'parula')

    a = plot(ax,all_data(i).ft.xf,-all_data(i).ft.cue,'c','LineWidth',1); a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(ax,all_data(i).ft.xb,all_data(i).im.mu,'w','LineWidth',1);   a.YData(abs(diff(a.YData))>pi) = nan;

    axis(ax,'tight') % NOT clamped to [-pi,pi] -- unwrap(im.alpha) spans both PB hemispheres (32 glomeruli); a fixed [-pi,pi] ylim would cut off the second one
    xlabel(ax,'time (s)')
    title(ax,sprintf('%s (trial %d, rho=%.2f)',fly_short_label(fly_list{f}),i,mean(all_data(i).im.rho,'omitnan')),'FontSize',8,'Interpreter','none')
end
title(tl,'im.z across the PB (color) with heading (cyan) and extracted bump position mu (white) overlaid -- one example trial per fly (highest mean rho)','Interpreter','none')

%% offset variability between mu and heading (circ_dist), by light condition
% same metric/convention as lpsp_p2x2_walking_script.m's own "how
% accurately does mu track the fly's actual heading" check: mu is
% upsampled from the imaging clock (xb) onto the faster fictrac clock
% (xf), restricted to timepoints where the bump estimate is confident
% (im.rho > rho_thresh), then this trial's own mu-heading OFFSET
% VARIABILITY is circ_var(circ_dist(cue,mu)) -- circular variance of the
% offset, not its mean, since the offset's absolute size is arbitrary
% (depends on where mu's zero-point happens to fall); what matters is how
% much it drifts. Light condition is read the same way as elsewhere in
% this codebase (e.g. lpsp_kir_claude.m): ft.pattern containing
% "background" = dark, otherwise closed loop. One point per trial (not
% pooled per fly), matching how lpsp_p2x2_walking_script.m itself plots
% this metric.
rho_thresh = .1; % minimum bump vector strength (im.rho) to trust a mu estimate

is_dark = arrayfun(@(x) contains(x.ft.pattern,'background','IgnoreCase',true), all_data);
light_labels = {'closed loop','dark'};
light_idx    = is_dark(:) + 1; % 1 = closed loop, 2 = dark

mu_offset_var = nan(n_trials,1);
for i = 1:n_trials
    xb  = all_data(i).ft.xb;
    xf  = all_data(i).ft.xf;
    cue = -all_data(i).ft.cue; % cue is stored with the opposite sign convention of mu

    mu  = interp1(xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');
    rho = interp1(xb,all_data(i).im.rho,xf,'linear','extrap');

    idx = rho > rho_thresh;
    if ~any(idx); continue; end
    mu_offset_var(i) = circ_var(circ_dist(cue(idx),mu(idx)));
end

fprintf('\n=== mu-heading offset variability (circ_var of circ_dist(cue,mu)), by light condition ===\n');
for li = 1:2
    valid = light_idx==li & ~isnan(mu_offset_var);
    fprintf('  %-12s n=%d trials, mean=%.3f\n', light_labels{li}, sum(valid), mean(mu_offset_var(valid),'omitnan'));
end

figure(3); clf
set(gcf,'Name','mu-cue offset variability by light condition','Position',[100,100,500,450])
hold on
swarmchart(light_idx,mu_offset_var,'filled','MarkerFaceAlpha',.5,'XJitterWidth',.3);
for li = 1:2
    valid = light_idx==li & ~isnan(mu_offset_var);
    errorbar(li,mean(mu_offset_var(valid),'omitnan'),std(mu_offset_var(valid),'omitnan')/sqrt(sum(valid)),'ok','LineWidth',1.5)
end
xticks(1:2); xticklabels(light_labels); xlim([.5,2.5])
ylabel('circular variance of mu-heading offset')
title(sprintf('mu-heading offset variability by light condition (n=%d trials, rho>%.2g)',sum(~isnan(mu_offset_var)),rho_thresh))

%% detect ATP ejection per trial -- two definitions, permissive and restrictive
% definition 1 (permissive, has_pulse_atp): same convention as
% lpsp_p2x2_walking_script.m's own section 3 -- average this trial's own
% stims (ft.stims rising edges, ALL of them) into one peri-stim atp trace
% on a common relative-time axis, call it "an ejection" if the post-stim
% peak clears peak_factor baseline-SDs above this trial's own pre-stim
% baseline. This only tells you the electrode fired and atp dye moved
% near the PB, not that it got IN.
%
% definition 2 (restrictive, has_pulse_notch): definition 1 AND the EPG
% calcium signal (im.d, summed across glomeruli, same convention as
% atp_sig/gcamp_sig in the figure below) shows a temporally defined
% post-stim dip -- a trough within notch_win that clears notch_factor
% baseline-SDs BELOW this trial's own pre-stim EPG baseline. This is
% evidence the ejected ATP actually reached the PB and perturbed EPG
% activity, not just that the electrode fired. notch_win defaults to the
% same post-stim window as the atp peak check (the notch, if it happens
% at all, isn't expected to lag the ejection itself by much) -- kept as
% its own parameter so it can be tuned independently.
%
% Trials without an atp channel at all (see header point 3) or with zero
% stims are left out of both definitions entirely, not counted as "no
% ejection".
peri_win     = 5;      % seconds before/after each stim to examine
t_common     = linspace(-peri_win,peri_win,101); % common relative-time axis stims are interpolated onto
peak_win     = [0,3];  % seconds after onset checked for an atp response peak
notch_win    = [0,3];  % seconds after onset checked for an EPG trough
base_win     = [-peri_win,0]; % seconds before onset used as this trial's own baseline
peak_factor  = 3;      % post-stim atp peak must clear this many baseline-SDs above baseline to count as "an ejection"
notch_factor = 3;      % post-stim EPG trough must clear this many baseline-SDs below baseline to count as "a notch"

has_pulse_atp   = false(n_trials,1); % definition 1: stim-locked atp peak
has_pulse_notch = false(n_trials,1); % definition 2: definition 1 AND a stim-locked EPG trough
atp_peri   = cell(n_trials,1); % this trial's stim-averaged atp trace (on t_common), kept for plotting/QC
gcamp_peri = cell(n_trials,1); % same, EPG calcium (im.d)

for i = 1:n_trials
    if ~(isfield(all_data(i),'atp') && isfield(all_data(i).atp,'d') && ~isempty(all_data(i).atp.d)); continue; end

    stims  = logical(all_data(i).ft.stims(:));
    onsets = find(diff([false;stims])==1); % rising edge of each stim
    if isempty(onsets); continue; end

    xb      = all_data(i).ft.xb;
    xf      = all_data(i).ft.xf;
    t_onset = xf(onsets);
    atp_sig   = sum(all_data(i).atp.d,1);
    gcamp_sig = sum(all_data(i).im.d,1);

    atp_aligned   = nan(numel(onsets),numel(t_common));
    gcamp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s); % imaging clock, re-centered on this trial's stims
        atp_aligned(s,:)   = interp1(xb_rel,atp_sig,  t_common,'linear',nan);
        gcamp_aligned(s,:) = interp1(xb_rel,gcamp_sig,t_common,'linear',nan);
    end
    atp_peri{i}   = mean(atp_aligned,  1,'omitnan');
    gcamp_peri{i} = mean(gcamp_aligned,1,'omitnan');

    base_idx  = t_common >= base_win(1)  & t_common <  base_win(2);
    peak_idx  = t_common >  peak_win(1)  & t_common <= peak_win(2);
    notch_idx = t_common >  notch_win(1) & t_common <= notch_win(2);

    atp_base_mean = median(atp_peri{i}(base_idx),'omitnan');
    atp_base_std  = std(atp_peri{i}(base_idx),'omitnan');
    atp_peak_amp  = median(atp_peri{i}(peak_idx),'omitnan');
    has_pulse_atp(i) = (atp_peak_amp - atp_base_mean) > peak_factor*atp_base_std;

    im_base_mean = median(gcamp_peri{i}(base_idx),'omitnan');
    im_base_std  = std(gcamp_peri{i}(base_idx),'omitnan');
    im_trough    = min(gcamp_peri{i}(notch_idx));
    has_pulse_notch(i) = has_pulse_atp(i) && (im_base_mean - im_trough) > notch_factor*im_base_std;
end

n_scored = sum(~cellfun(@isempty,atp_peri));
fprintf('\n=== ejection detection (per trial) ===\n');
fprintf('%d/%d trials have an atp channel and >=1 stim\n', n_scored, n_trials);
fprintf('definition 1 (atp peak only):        %d/%d trials (%.0f%%)\n', sum(has_pulse_atp),   n_scored, 100*sum(has_pulse_atp)/max(n_scored,1));
fprintf('definition 2 (atp peak + EPG notch): %d/%d trials (%.0f%%)\n', sum(has_pulse_notch), n_scored, 100*sum(has_pulse_notch)/max(n_scored,1));

%% figure: per-fly overview -- trials, light condition, and which trials are labelled an ejection (either definition)
% one row per fly, one column per that fly's OWN trial (chronological,
% left-aligned -- flies have different trial counts, so shorter flies
% just leave later columns blank). Pixel color = light condition (closed
% loop vs dark). Overlay: a small black star marks a definition-1 ejection
% (atp peak only); a larger gold star marks definition 2 (atp peak + EPG
% notch, a subset of definition 1). Meant to make it visually obvious how
% few trials survive each definition and how unevenly they're spread
% across flies/light conditions -- directly the "are we throwing away a
% lot of data" question the pre/post figure above raised.
max_trials = max(trials_per_fly);
light_grid = nan(n_flies,max_trials);   % 1 = closed loop, 2 = dark, nan = no trial at this column for this fly
def1_grid  = false(n_flies,max_trials);
def2_grid  = false(n_flies,max_trials);

for f = 1:n_flies
    idx = find(fly_num==f); % chronological within this fly -- see header point 1 (natsortfiles)
    light_grid(f,1:numel(idx)) = light_idx(idx);
    def1_grid(f,1:numel(idx))  = has_pulse_atp(idx);
    def2_grid(f,1:numel(idx))  = has_pulse_notch(idx);
end

fprintf('\n=== per-fly trial counts, by light condition and ejection definition ===\n');
fprintf('%-20s %-10s %-10s %-10s %-10s\n','fly','closed','dark','def1','def2');
for f = 1:n_flies
    fprintf('%-20s %-10d %-10d %-10d %-10d\n', fly_short_label(fly_list{f}), ...
        sum(light_grid(f,:)==1), sum(light_grid(f,:)==2), sum(def1_grid(f,:)), sum(def2_grid(f,:)));
end

figure(6); clf
set(gcf,'Name','per-fly overview: light condition + ejection labels','Position',[100,100,900,50+30*n_flies])
ax = axes; hold(ax,'on')
h_im = imagesc(ax,light_grid);
set(h_im,'AlphaData',~isnan(light_grid)) % blank (no-trial) cells show through as white, not a color
colormap(ax,[1,0.85,0.2; 0.15,0.15,0.15]) % 1 = closed loop (yellow), 2 = dark (near-black)
clim(ax,[1,2])
set(ax,'Color',[1,1,1])

[fr1,fc1] = find(def1_grid & ~def2_grid);
[fr2,fc2] = find(def2_grid);
scatter(ax,fc1,fr1,80, 'k',            '*','LineWidth',1.2)
scatter(ax,fc2,fr2,160,[0.85,0.6,0.05],'*','LineWidth',1.5)

yticks(ax,1:n_flies); yticklabels(ax,arrayfun(@(f) fly_short_label(fly_list{f}),1:n_flies,'UniformOutput',false))
xlabel(ax,'trial (chronological within fly)')
axis(ax,'tight')
title(ax,'light condition per trial (color) and ejection labels (stars), per fly')

legend_h = [patch(ax,nan,nan,[1,0.85,0.2]),patch(ax,nan,nan,[0.15,0.15,0.15]), ...
            scatter(ax,nan,nan,80,'k','*'),scatter(ax,nan,nan,160,[0.85,0.6,0.05],'*')];
legend(ax,legend_h,{'closed loop','dark','definition 1 ejection','definition 2 ejection'},'Location','eastoutside')

%% figure: atp and epg calcium traces aligned to stim onset, for every trial with a detected ejection (definition 1)
% same layout as lpsp_p2x2_walking_script.m's own sanity-check figure: one
% panel per trial, atp (red, left axis) and gcamp (green, right axis)
% mean +/- SEM (plot_sem) across that trial's own stims. Panels that also
% clear the stricter definition 2 (EPG notch) are marked in their title.
pulse_trials = find(has_pulse_atp);
n_col_p = ceil(sqrt(numel(pulse_trials)));
n_row_p = ceil(numel(pulse_trials)/n_col_p);

figure(4); clf
set(gcf,'Name','atp + gcamp aligned to stim onset','Position',[50,50,220*n_col_p,180*n_row_p])
for k = 1:numel(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets  = find(diff([false;stims])==1);
    t_onset = xf(onsets);

    subplot(n_row_p,n_col_p,k); hold on
    if has_pulse_notch(i); notch_tag = ', +notch'; else; notch_tag = ''; end
    title(sprintf('%s (trial %d, n=%d stims%s)',fly_short_label(all_data(i).meta),i,numel(onsets),notch_tag),'FontSize',7,'Interpreter','none')

    atp_sig   = sum(all_data(i).atp.d,1);
    gcamp_sig = sum(all_data(i).im.d,1);

    atp_aligned   = nan(numel(onsets),numel(t_common));
    gcamp_aligned = nan(numel(onsets),numel(t_common));
    for s = 1:numel(onsets)
        xb_rel = xb - t_onset(s);
        atp_aligned(s,:)   = interp1(xb_rel,atp_sig,  t_common,'linear',nan);
        gcamp_aligned(s,:) = interp1(xb_rel,gcamp_sig,t_common,'linear',nan);
    end

    yyaxis left
    h = plot_sem(gca,t_common,atp_aligned); h.FaceColor = 'r';
    ylabel('atp')

    yyaxis right
    h = plot_sem(gca,t_common,gcamp_aligned); h.FaceColor = 'g';
    ylabel('gcamp')

    plot([0,0],ylim,':k')
    xlim([-peri_win,peri_win])
end
sgtitle(sprintf('average atp (red, left axis) and gcamp (green, right axis) fluorescence around each stim -- %d/%d trials with a detected ejection (definition 1); "+notch" = also clears definition 2',numel(pulse_trials),n_scored))

%% bump mobility (bump path length / heading path length) per trial
% same convention as lpsp_p2x2_walking_script.m's own section 2: detect
% walking bouts from the fly's own rotational speed (r_speed -- NOT cue
% position, since a trial's visual cue can move without the fly actually
% walking), then regress bump path length (from the re-estimated mu) against
% heading path length (integral of |r_speed|) across those bouts, weighted
% by each bout's own duration. mov_ratio(i) = 1 means the bump moves exactly
% as much as the fly turns; this is the per-trial quantity the pre/post
% comparison below pools per fly.
smooth_window      = 60;    % samples (~1s at xf's ~60Hz), gaussian smoothing window
turn_thresh        = .25;   % rad/s, minimum heading speed to call a bout "walking"
max_gap_frames     = .5*60; % frames of non-walking allowed within a bout before splitting it
min_walking_frames = .5*60; % minimum bout length to keep

mov_ratio    = nan(n_trials,1);
bout_mov_mu  = cell(n_trials,1); % per-bout bump path length
bout_mov_cue = cell(n_trials,1); % per-bout heading path length (despite the name -- kept for consistency with lpsp_p2x2_walking_script.m -- this is r_speed-based, not cue position)
bout_dur     = cell(n_trials,1); % duration of each bout (s), used to weight the regression

for i = 1:n_trials
    xf = all_data(i).ft.xf;
    dt = median(diff(xf));

    mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');

    r_speed_smooth = smoothdata(all_data(i).ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    mov_mu  = nan(length(bout_starts),1);
    mov_cue = nan(length(bout_starts),1);
    dur     = nan(length(bout_starts),1);
    for b = 1:length(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(r_speed_smooth(bout_starts(b):bout_ends(b))),'omitnan')*dt;
        dur(b)     = (bout_ends(b)-bout_starts(b)+1)*dt;
    end

    w = sqrt(dur); % weighted least squares via the standard sqrt(weight) rescaling trick
    mov_ratio(i)    = (w.*mov_cue) \ (w.*mov_mu);
    bout_mov_mu{i}  = mov_mu;
    bout_mov_cue{i} = mov_cue;
    bout_dur{i}     = dur;
end

%% bump mobility before vs after perturbation, pooled per fly, by light condition -- both ejection definitions
% same convention as lpsp_p2x2_walking_script.m's own section 4 (minus the
% genotype split that dataset needed -- single group here): for each fly,
% its OWN first detected-ejection trial (chronological order, which
% all_data already has -- see header point 1's natsortfiles) marks the
% split; "pre" = that fly's trials before it, "post" = after it (excluding
% any further detected-ejection trials). Within pre/post, every walking
% bout from that fly's own trials in ONE light condition is pooled into a
% single weighted regression (min_bouts_per_condition floor), rather than
% averaging each trial's own often-underpowered mov_ratio -- same
% reasoning as the reference script. Run once per ejection definition so
% the two can be compared side by side.
min_bouts_per_condition = 0;

ejection_defs   = {has_pulse_atp, has_pulse_notch};
def_labels      = {'definition 1: atp peak only','definition 2: atp peak + EPG notch'};
fly_pre_defs    = cell(1,2);
fly_post_defs   = cell(1,2);

for def_i = 1:2
    [fly_pre_defs{def_i},fly_post_defs{def_i}] = pre_post_bump_mobility( ...
        ejection_defs{def_i}, fly_num, n_flies, bout_mov_mu, bout_mov_cue, bout_dur, is_dark, min_bouts_per_condition);
end

figure(5); clf
set(gcf,'Name','bump mobility before vs after perturbation','Position',[100,100,800,700])
for def_i = 1:2
    fly_pre  = fly_pre_defs{def_i};
    fly_post = fly_post_defs{def_i};
    for li = 1:2
        subplot(2,2,(def_i-1)*2+li); hold on

        valid = ~isnan(fly_pre(:,li)) & ~isnan(fly_post(:,li));
        pre_vals  = fly_pre(valid,li);
        post_vals = fly_post(valid,li);

        plot([ones(sum(valid),1),2*ones(sum(valid),1)]',[pre_vals,post_vals]','Color',[.5,.5,.5,.3])
        scatter(ones(sum(valid),1), pre_vals,'filled','MarkerFaceAlpha',.5)
        scatter(2*ones(sum(valid),1),post_vals,'filled','MarkerFaceAlpha',.5)
        errorbar([1,2],[mean(pre_vals,'omitnan'),mean(post_vals,'omitnan')], ...
                       [std(pre_vals,'omitnan'),std(post_vals,'omitnan')]/sqrt(max(sum(valid),1)),'-ok','LineWidth',1.5)

        xlim([.5,2.5]); xticks([1,2]); xticklabels({'pre','post'})
        ylabel('bump path length / heading path length')
        title(sprintf('%s\n%s (n=%d flies)',def_labels{def_i},light_labels{li},sum(valid)),'FontSize',9)
    end
end
sgtitle(sprintf('bump mobility before vs after perturbation, pooled per fly (min %d bouts/condition), rows = ejection definition, columns = light condition',min_bouts_per_condition))

%% save figures
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\dopamine_ionto_walking';
if ~isfolder(fig_dir)
    mkdir(fig_dir)
end

fig_handles = findobj('Type','figure');
[~,order] = sort(arrayfun(@(f) f.Number, fig_handles));
fig_handles = fig_handles(order);

save_failed = {};
for k = 1:numel(fig_handles)
    fig = fig_handles(k);
    fig_name = get(fig,'Name');
    if isempty(fig_name); fig_name = sprintf('figure_%d',fig.Number); end
    safe_name = regexprep(fig_name,'[^\w\-]+','_');
    out_path  = fullfile(fig_dir,sprintf('fig%02d_%s.pdf',fig.Number,safe_name));
    try
        exportgraphics(fig,out_path)
        fprintf('saved %s\n', out_path);
    catch ME
        save_failed{end+1} = out_path; %#ok<AGROW>
        warning('could not save %s (%s) -- is it open in another program?', out_path, ME.message);
    end
end
if ~isempty(save_failed)
    fprintf('\n%d figure(s) failed to save -- close them in any viewer and re-run to update:\n', numel(save_failed));
    fprintf('  %s\n', save_failed{:});
end

%% functions

function lbl = fly_short_label(fly_path)
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N, or
    % from a full trial meta path -- strips down to the last 2 non-empty
    % path parts either way.
    parts = strsplit(fly_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        lbl = strjoin(parts(max(1,end-1):end),' ');
    else
        lbl = strjoin(parts(max(1,fly_part(1)-1):fly_part(1)),' ');
    end
end

function [fly_pre,fly_post] = pre_post_bump_mobility(has_pulse_def, fly_num, n_flies, bout_mov_mu, bout_mov_cue, bout_dur, is_dark, min_bouts_per_condition)
    % for each fly: find its own trials (fly_num==f, already chronological
    % since all_data is built in natsortfiles order), take the FIRST trial
    % that clears has_pulse_def as the perturbation split point. "pre" is
    % everything strictly before that trial; "post" is EVERY trial from
    % directly after it onward, through the end of that fly's session --
    % including any later trial that itself also clears has_pulse_def (a
    % second/third ejection doesn't stop counting as "post", it's still
    % after the first perturbation) and regardless of light condition (a
    % dark trial following a closed-loop perturbation, or vice versa,
    % still counts as "post" -- the light-condition split below just bins
    % pre/post separately per condition, it doesn't require post to match
    % whatever condition the perturbation trial itself was in). Bouts are
    % pooled per light condition into one weighted regression -- same
    % pooling as lpsp_p2x2_walking_script.m's own section 4, minus its
    % genotype dimension.
    fly_pre  = nan(n_flies,2); fly_post = nan(n_flies,2); % columns: [closed loop, dark]
    for f = 1:n_flies
        trial_idx = find(fly_num==f);
        pulse_pos = find(has_pulse_def(trial_idx));
        if isempty(pulse_pos); continue; end

        pre_idx  = trial_idx(1:pulse_pos(1)-1);
        post_idx = trial_idx(pulse_pos(1)+1:end);

        for li = 1:2 % 1 = closed loop, 2 = dark
            pre_this  = pre_idx(is_dark(pre_idx)   == (li==2));
            post_this = post_idx(is_dark(post_idx) == (li==2));
            if isempty(pre_this) || isempty(post_this); continue; end

            pre_cue  = vertcat(bout_mov_cue{pre_this});  pre_mu  = vertcat(bout_mov_mu{pre_this});  pre_dur  = vertcat(bout_dur{pre_this});
            post_cue = vertcat(bout_mov_cue{post_this}); post_mu = vertcat(bout_mov_mu{post_this}); post_dur = vertcat(bout_dur{post_this});

            if length(pre_cue) >= min_bouts_per_condition
                w = sqrt(pre_dur);
                fly_pre(f,li) = (w.*pre_cue) \ (w.*pre_mu);
            end
            if length(post_cue) >= min_bouts_per_condition
                w = sqrt(post_dur);
                fly_post(f,li) = (w.*post_cue) \ (w.*post_mu);
            end
        end
    end
end

function h = plot_sem(ax,t,x)
    % shades mean(x) +/- sem(x) over rows of x (one row per replicate,
    % columns matching t); face color is set by the caller. Unchanged from
    % lpsp_p2x2_walking_script.m's own plot_sem.
    t = reshape(t,1,[]);
    m1 = mean(x,1,'omitnan');
    s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

    idx = ~isnan(m1);
    m1 = m1(idx);
    s1 = s1(idx);
    t  = t(idx);

    h = patch(ax,[t,fliplr(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end

function s = mu_rho_from_f(f_raw, smooth_win, f0_pct)
    % re-derive dff/zscore/mu/rho from a moving-average-smoothed copy of an
    % already-extracted per-glomerulus raw fluorescence trace (f_raw, rows
    % = glomeruli, columns = time) -- same math as the tail end of
    % process_im, just starting from f_raw instead of a raw movie, so it
    % can be re-run cheaply (e.g. to try a different smoothing window)
    % without re-touching the raw registered movies.
    f_smooth = smoothdata(f_raw,2,'movmean',smooth_win);
    f0       = prctile(f_smooth,f0_pct,2);
    dff      = (f_smooth - f0) ./ f0;
    zs       = zscore(dff,[],2);

    n_centroid = size(f_raw,1)/2; % per hemisphere; process_im always doubles this for the full alpha vector
    alpha      = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    alpha      = repmat(alpha,1,2);

    [x_tmp,y_tmp] = pol2cart(alpha,zs');
    [mu,rho]      = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = mod(mu,2*pi);
    mu(mu > pi) = mu(mu > pi) - 2*pi;

    s.d   = dff;
    s.z   = zs;
    s.mu  = mu;
    s.rho = rho;
end

function cue_smooth = smooth_circular(cue, win)
    % moving-average smoothing for a wrapped angle (e.g. heading) --
    % unwrap first so the average near a +-pi wrap boundary isn't
    % corrupted, then re-wrap back to (-pi,pi], same convention process_ft
    % already uses for its own cue smoothing.
    cue_smooth = unwrap(cue);
    cue_smooth = smoothdata(cue_smooth,1,'movmean',win,'omitnan');
    cue_smooth = mod(cue_smooth,2*pi);
    cue_smooth(cue_smooth > pi) = cue_smooth(cue_smooth > pi) - 2*pi;
end

function fid = trial_fly_id(meta_path)
    % meta = ...\<date folder>\fly N\<trial folder>\registration
    parts = strsplit(meta_path,filesep);
    parts(cellfun(@isempty,parts)) = [];
    fly_part = find(~cellfun(@isempty,regexpi(parts,'^fly\s*\d+','once')));
    if isempty(fly_part)
        fid = meta_path; % shouldn't happen; fall back to a unique id per trial
    else
        fid = strjoin(parts(1:fly_part(1)),filesep); % date folder + "fly N"
    end
end

function s = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type)
    % unchanged from dopamine_ionto_script.m -- see header point 4 for why
    % this dataset's field names already match what this function expects.
    f_speed = ftData_dat.velFor{:};
    f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}' / 192 * 2 * pi - pi;
    cue(abs(gradient(cue)) > 2) = nan;
    cue     = unwrap(cue);
    cue     = smoothdata(cue,1,ft_type,ft_win,'omitnan');
    cue     = mod(cue,2*pi);
    cue(cue > pi) = cue(cue > pi) - 2*pi;

    s.xf      = seconds(ftData_DAQ.trialTime{:});
    if ismember('volClock',ftData_DAQ.Properties.VariableNames)
        s.xb  = seconds(ftData_DAQ.volClock{:});
    end
    s.f_speed = smoothdata(f_speed,1,ft_type,ft_win);
    s.r_speed = smoothdata(r_speed,1,ft_type,ft_win);
    s.cue     = cue;
end

function s = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct)
    % unchanged from dopamine_ionto_script.m -- imgData is expected as
    % Y x X x T (see header point 2: this dataset's img{1}/img{2} are
    % already in that shape, so no sum-over-Z step happens before this is
    % called).
    imgData = smoothdata(imgData,3,im_type{1},im_win{1});
    imgData = imgData - min(imgData,[],'all');

    [y_mask,x_mask] = find(mask);
    min_axis        = min(range(x_mask),range(y_mask));
    mid             = bwskel(mask,'MinBranchLength',min_axis);
    [y_mid,x_mid]   = find(mid);
    ep              = bwmorph(mid,'endpoints');
    [y0,x0]         = find(ep,1); %#ok<ASGLU>

    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);

    xq          = [-min_axis:(length(x_mid)+min_axis)];
    x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
    y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

    idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');
    x_mid       = x_mid(idx);
    y_mid       = y_mid(idx);

    xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];
    centroids   = centroids(2:2:end-1,:);
    [~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1);

    imgData_2d      = reshape(imgData,[],size(imgData,3));
    centroid_log    = false(2*n_centroid,size(imgData_2d,1));
    for i = 1:2*n_centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);
    f0              = prctile(f_cluster,f0_pct,2);
    dff_cluster     = (f_cluster - f0) ./ f0;
    zscore_cluster  = zscore(dff_cluster,[],2);

    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    alpha       = repmat(alpha,1,2);

    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});

    mask_2d = reshape(mask,[],1);
    f_mask  = mean(imgData_2d(mask_2d,:),1);
    f_nmask = mean(imgData_2d(~mask_2d,:),1);

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.d  = dff_cluster;
    s.f  = f_cluster;
    s.alpha = alpha;
    s.f_mask = f_mask;
    s.f_nmask= f_nmask;
end
