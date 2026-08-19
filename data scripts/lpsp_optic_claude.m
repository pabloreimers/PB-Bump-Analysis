%% lpsp_optic_claude
% Loads the LPsP open-loop optic-flow dataset (lpsp_ol_data_20240104.mat)
% and reports how many flies are in it and each fly's average forward
% walking speed.
%
% This dataset is a table (variable "data", 48 rows = 48 trials), not the
% all_data struct array used elsewhere in this codebase (e.g.
% lpsp_kir_claude.m) -- see optic_flow_analysis.m, section "%% save
% variables", which built and saved it. Each row's own "trial_names"
% entry is only the SESSION DATE (e.g. "20231108"), not a full trial
% folder name or fly identifier -- confirmed directly (48/48 entries are
% bare 8-digit dates, several repeated on the same date for multiple
% trials/flies imaged that day).
%
% Fly identity (which trials came from the same physical fly) is not
% recoverable from the data table or from the raw folder tree alone --
% unlike lpsp_kir_redo's raw tree, this dataset's raw sessions
% (Z:\pablo\lpsp ol\analysis\<date>\<date>-N_LPsP_syt7f\) are organized by
% date and trial number only, with no "fly N" subfolder. The only record
% of which trials belong to the same fly is the hand-written trial->fly
% lookup table in optic_flow_analysis.m itself (section "%% save meta
% data for each fly") -- fly_num is reproduced verbatim below from that
% table, since there is no other source for it.
%
% That hand-written table has 49 trial entries, but this dataset has only
% 48 rows: trial 20231121-3_LPsP_syt7f (assigned fly 13 there) does not
% exist under Z:\pablo\lpsp ol\analysis\20231121\ (confirmed directly --
% only trials 1,2,4,5,6,7 are present for that date; trial 3 sits
% unprocessed at the top level of Z:\pablo\lpsp ol\ instead), so it isn't
% part of this data pull and fly 13 never appears here. Also corrected
% here: the hand-written table lists trial 20231206-1/-2 as fly 23, but no
% 20231206 folder exists anywhere under Z:\pablo\lpsp ol\ -- only
% 20231207 does (confirmed directly) -- a one-day typo in that table, not
% a real extra session; fixed to 20231207 below.
%
% Since data.trial_names only stores the date, each row is matched to its
% own trial NUMBER by re-deriving it from the raw folder tree: for every
% date, the analysis subfolders are listed and sorted by ascending trial
% number, and that order is assumed to match the row order optic_flow_
% analysis.m originally produced via dir()'s own (alphabetical, so also
% ascending for these single-digit trial numbers) traversal order -- the
% per-date row/folder COUNTS are checked to match as a sanity check on
% this assumption (they do, exactly, for all 14 dates in this dataset).
% Trial number + date is then looked up against the fly table above to
% get each row's fly_num.

%% load data
data_dir    = fullfile('..','.data');
source_file = 'lpsp_ol_data_20240104.mat';
base_dir    = 'Z:\pablo\lpsp ol\analysis';
assert(isfolder(base_dir), 'cannot find %s -- check drive mapping', base_dir)

tmp  = load(fullfile(data_dir,source_file),'data');
data = tmp.data;
n_trials = height(data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% hand-written trial -> fly lookup, reproduced from optic_flow_analysis.m
% (date, trial number, fly number, rearing food 'g'=german/'m'=molasses --
% the same "%% save meta data for each fly" table, its 3rd column); see
% header for the one date fix (20231206 -> 20231207) and the one dropped
% trial (20231121-3, fly 13, not present in this data pull). Cross-checked
% directly against that table: all 49 (date,trial)->(fly,food) rows here
% match it exactly (accounting for the one date fix above).
fly_table = {
    '20231108', 2, 1,  'm';
    '20231110', 1, 2,  'm';
    '20231110', 2, 2,  'm';
    '20231111', 1, 3,  'm';
    '20231111', 2, 3,  'm';
    '20231112', 1, 4,  'm';
    '20231112', 2, 4,  'm';
    '20231112', 3, 5,  'g';
    '20231112', 4, 5,  'g';
    '20231113', 1, 6,  'g';
    '20231113', 2, 6,  'g';
    '20231113', 3, 7,  'g';
    '20231113', 4, 7,  'g';
    '20231115', 1, 8,  'g';
    '20231115', 2, 8,  'g';
    '20231115', 3, 9,  'g';
    '20231115', 4, 9,  'g';
    '20231115', 5, 10, 'g';
    '20231115', 6, 10, 'g';
    '20231115', 7, 11, 'g';
    '20231121', 1, 12, 'm';
    '20231121', 2, 12, 'm';
    '20231121', 3, 13, 'm'; % not present in this dataset (see header)
    '20231121', 4, 14, 'm';
    '20231121', 5, 14, 'm';
    '20231121', 6, 15, 'm';
    '20231121', 7, 15, 'm';
    '20231122', 1, 16, 'm';
    '20231127', 1, 17, 'm';
    '20231127', 2, 17, 'm';
    '20231127', 4, 18, 'm';
    '20231204', 1, 19, 'm';
    '20231204', 2, 19, 'm';
    '20231204', 3, 20, 'm';
    '20231204', 4, 20, 'm';
    '20231204', 5, 21, 'm';
    '20231204', 6, 21, 'm';
    '20231204', 7, 22, 'm';
    '20231204', 8, 22, 'm';
    '20231207', 1, 23, 'm'; % date corrected from '20231206' (see header)
    '20231207', 2, 23, 'm';
    '20231215', 1, 24, 'g';
    '20231215', 2, 24, 'g';
    '20231215', 3, 25, 'g';
    '20231215', 4, 26, 'g';
    '20231216', 1, 27, 'g';
    '20231216', 2, 28, 'g';
    '20231218', 1, 29, 'g';
    '20231218', 2, 30, 'g';
    };

%% recover each row's trial number from the raw folder tree, then look up its fly number
trial_num = nan(n_trials,1);
fly_num   = nan(n_trials,1);

[dates_present,~,date_group] = unique(data.trial_names);
for d = 1:numel(dates_present)
    this_date = dates_present{d};
    rows      = find(date_group==d);

    d_listing = dir(fullfile(base_dir,this_date));
    d_listing = d_listing([d_listing.isdir]);
    tnums = nan(numel(d_listing),1);
    for k = 1:numel(d_listing)
        m = regexp(d_listing(k).name,['^',this_date,'-(\d+)_'],'tokens','once');
        if ~isempty(m)
            tnums(k) = str2double(m{1});
        end
    end
    tnums = sort(tnums(~isnan(tnums)));

    if numel(tnums) ~= numel(rows)
        warning('%s: %d trial folders on disk but %d rows in data -- fly assignment for this date may be wrong', ...
            this_date, numel(tnums), numel(rows));
        continue
    end

    for k = 1:numel(rows)
        trial_num(rows(k)) = tnums(k);
        match = strcmp(fly_table(:,1),this_date) & cell2mat(fly_table(:,2))==tnums(k);
        if ~any(match)
            warning('no fly_table entry for %s trial %d (row %d)', this_date, tnums(k), rows(k));
            continue
        end
        fly_num(rows(k)) = fly_table{match,3};
    end
end

n_unmatched = sum(isnan(fly_num));
if n_unmatched > 0
    warning('%d/%d trials could not be assigned a fly number', n_unmatched, n_trials)
end

fly_list = unique(fly_num(~isnan(fly_num)));
n_flies  = numel(fly_list);
fprintf('\n%d trials -> %d flies\n', n_trials, n_flies);

%% rearing food per fly ('g'=german, 'm'=molasses), from fly_table's 4th column
% food is a fly-level trait (constant across a fly's own trials in
% fly_table), so each fly's food is just its first matching row's own
% entry -- checked directly, no fly in fly_table has inconsistent food
% across its rows.
fly_food = repmat(' ',n_flies,1);
for f = 1:n_flies
    match = cell2mat(fly_table(:,3))==fly_list(f);
    foods = unique(fly_table(match,4));
    if numel(foods) > 1
        warning('fly %d has inconsistent food labels in fly_table: %s', fly_list(f), strjoin(foods,','));
    end
    fly_food(f) = foods{1};
end
fprintf('  molasses: %d flies, german: %d flies\n', sum(fly_food=='m'), sum(fly_food=='g'));

%% average forward speed per fly
% data.f_speed is each trial's fictrac forward velocity (velFor, signed --
% positive/negative for forward/backward), smoothed, at the fictrac frame
% rate. "forward speed" here follows the same convention already used for
% this exact dataset in optic_flow_analysis.m ("%% compare the overall
% activity of flies in the two groups": f_dist = mean(abs(f_speed)),
% labeled "mean forward speed (mm/s)") -- the magnitude of forward
% velocity, since a raw (signed) mean would let backward bouts cancel out
% forward ones. A fly's trials are pooled (concatenated) before averaging,
% not averaged trial-by-trial then re-averaged, so a fly with 2 trials
% isn't weighted differently from a fly with 1.
fly_speed = nan(n_flies,1);
for f = 1:n_flies
    rows = find(fly_num==fly_list(f));
    pooled = cat(1,data.f_speed{rows});
    fly_speed(f) = mean(abs(pooled),'omitnan');
end

fprintf('\n=== average forward speed per fly (mm/s) ===\n');
for f = 1:n_flies
    fprintf('  fly %2d (n=%d trials): %.3f mm/s\n', fly_list(f), sum(fly_num==fly_list(f)), fly_speed(f));
end

%% figure: average forward speed per fly
figure(1); clf
set(gcf,'Name','average forward speed per fly')
bar(fly_list,fly_speed)
xlabel('fly number')
ylabel('mean |forward speed| (mm/s)')
title(sprintf('lpsp\\_ol: average forward speed per fly (n=%d flies, %d trials)', n_flies, n_trials))
xticks(fly_list)

%% open loop speed at every time point, from the velocity of ft.cue
% cue (data.cue) is the panel/grating position (rad, wrapped to
% [-pi,pi]) -- already unwrapped, heavily gaussian-smoothed, then
% rewrapped when this dataset was built (see ft_calc in
% optic_flow_analysis.m). Its raw instantaneous gradient still carries a
% lot of high-frequency noise: the underlying cuePos comes off a
% 192-step rotary encoder, so single-frame quantization jitter survives
% the smoothing pass and shows up as brief, non-physical velocity spikes
% after differentiating -- confirmed directly: plain
% abs(gradient(unwrap(cue),xf)) reaches speeds up to ~14.5 rad/s, even
% though the open-loop stimulus was only ever driven at discrete
% commanded speeds from 0 to 1.0 rad/s in 0.1 steps (see
% optic_flow_analysis.m's own os_vec = round(0:.1:1,1), and its
% independently-computed c_speed field, which used a heavy movmedian
% filter for exactly this reason).
%
% Same fix applied here: differentiate, then movmedian-smooth over
% ol_speed_smooth_win samples (1000 samples @ ~60Hz fictrac rate =
% ~16.7s, matching the window c_speed itself used) before rounding to
% the nearest 0.1 rad/s. This recovers exactly the expected 0-1.0 rad/s
% range (confirmed: max over all 48 trials is exactly 1.0) and agrees
% with the dataset's own stored c_speed at ~96.5% of timepoints (mean
% across trials, min 95.1%; the rest sit right at speed-step
% transitions, where any smoothing window necessarily blends two
% adjacent commanded speeds -- unavoidable without the original raw,
% un-smoothed cuePos trace, which isn't stored in this dataset).
ol_speed_smooth_win = 1000; % samples (~16.7s @ 60Hz) -- matches c_speed's own filter window

ol_speed = cell(n_trials,1);
for i = 1:n_trials
    raw_speed   = abs(gradient(unwrap(data.cue{i}),data.xf{i}));
    ol_speed{i} = round(smoothdata(raw_speed,'movmedian',ol_speed_smooth_win),1);
end

match_frac = cellfun(@(a,b) mean(a==b), ol_speed, data.c_speed);
fprintf('\nderived open-loop speed matches this dataset''s own c_speed field on %.1f%% of timepoints (mean across trials, min %.1f%%)\n', ...
    100*mean(match_frac), 100*min(match_frac));

%% mean PB activity (im.d convention: raw dF/F, NOT z-scored), aligned to the fictrac timebase
% previously this used a z-scored mean-across-clusters (im.z convention:
% z-score each cluster over time, then average across the 40 clusters --
% confirmed ~0.99 correlated in shape with the raw version per trial, but
% z-scoring rescales every trial/fly to its own unit variance, which can
% distort cross-fly comparisons for flies with naturally weak vs. strong
% bump signals).
%
% Using im.d instead (raw dF/F, unnormalized) means using data.dff_tot
% directly -- it's already exactly this quantity: bump_calc_pb (in
% optic_flow_analysis.m) computes amp_tot = mean(dff_cluster,1) (mean raw
% dF/F across the 40 clusters) at the imaging frame rate, then
% interpolates it onto the fictrac timebase (xf) before storing it, so no
% further interpolation is needed here to line it up with ol_speed/r_speed.
mean_activity = data.dff_tot;

%% figure: mean PB activity vs. open-loop speed, one line per fly + group mean +/- SEM
% each fly's own trials are pooled (concatenated across timepoints)
% before binning by speed, same per-fly-not-per-trial convention used
% throughout this codebase (e.g. lpsp_kir_claude.m) -- a fly with 2
% trials contributes exactly one line/mean, not two. Speed values are
% compared with a small tolerance rather than exact equality, since
% ol_speed and speed_levels can each round to floating-point
% representations of e.g. 0.3 that differ in the last bit.
speed_levels = 0:0.1:1;
speed_tol    = 1e-6;

fly_activity_by_speed = plot_activity_by_speed(1,ol_speed,mean_activity,fly_num,fly_list,fly_food,speed_levels,speed_tol,n_flies,n_trials);

%% figure: same, restricted to timepoints where the fly's own rotation has been below rot_thresh for a sustained preceding window
% isolates visually-driven activity from self-generated turning -- while
% the fly is actively rotating, r_speed itself correlates with bump
% velocity (see lpsp_kir_claude.m's own gain analysis) and would
% confound a pure open-loop-speed tuning curve.
%
% "not turning" here is a SUSTAINED, causal criterion, not an
% instantaneous one -- adopted from optic_flow_analysis.m's own still/flow
% comparison ("%% establish significance": speed_thresh + speed_win +
% smoothdata(r_speed<speed_thresh,'movmean',[speed_win,0])==1). A
% timepoint only counts as "still" if r_speed has been below rot_thresh
% for every one of the preceding speed_win samples, not just at that
% single instant -- smoothdata's backward-only window ([speed_win,0], 0
% samples forward) applied to a 0/1 indicator equals 1 only when every
% sample in that window was below threshold. This is stricter than a
% plain instantaneous r_speed<rot_thresh check: a frame right after a
% fast turn stays excluded until the fly has been calm for the whole
% window, not just for that one frame. r_speed is already an absolute
% value (see ft_calc).
rot_thresh = 0.5;   % rad/s -- "fly is not turning"; user-adjustable
speed_win  = 3*60;  % samples (~3s @ ~60Hz fictrac rate) -- same preceding-window length optic_flow_analysis.m used

still_activity = cell(n_trials,1);
for i = 1:n_trials
    is_still = smoothdata(data.r_speed{i} < rot_thresh,1,'movmean',[speed_win,0]) == 1;
    a = mean_activity{i};
    a(~is_still) = nan;
    still_activity{i} = a;
end

fly_activity_by_speed_still = plot_activity_by_speed(2,ol_speed,still_activity,fly_num,fly_list,fly_food,speed_levels,speed_tol,n_flies,n_trials, ...
    sprintf(', rotation < %.1f rad/s for preceding %.0fs',rot_thresh,speed_win/60));

%% figure: open-loop speed and PB activity (dff_tot) traces, one subplot per fly
% diagnostic sanity-check figure, not a summary statistic: raw ol_speed
% and mean_activity (=data.dff_tot) traces over time, one tile per fly, so
% the two derived/reused signals feeding every figure above can be
% eyeballed directly. Speed (left axis) and dF/F (right axis) are plotted
% on separate y-axes since their scales don't overlap (0-1 rad/s vs. raw
% dF/F). A fly with more than one trial has them concatenated end-to-end
% by elapsed time (not real wall-clock -- trials aren't recorded
% back-to-back) with a dotted vertical line at the trial boundary, same
% convention as lpsp_kir_claude.m's concat_trials_im.
n_cols = ceil(sqrt(n_flies));
n_rows = ceil(n_flies/n_cols);

figure(4); clf
set(gcf,'Name','open-loop speed and PB activity (dff_tot), per fly','Position',[50,50,1600,900])
t_trace = tiledlayout(n_rows,n_cols,'TileSpacing','compact','Padding','compact');

for f = 1:n_flies
    rows_f = find(fly_num==fly_list(f));
    [x_cat,speed_cat,act_cat,bounds] = concat_fly_trace(data.xf(rows_f),ol_speed(rows_f),mean_activity(rows_f));

    ax = nexttile(t_trace); hold(ax,'on')
    yyaxis(ax,'left')
    plot(ax,x_cat,speed_cat,'Color',[0,0.4470,0.7410],'LineWidth',0.5)
    ylim(ax,[0,1.05])
    yyaxis(ax,'right')
    plot(ax,x_cat,act_cat,'Color',[0.8500,0.3250,0.0980],'LineWidth',0.5)

    for b = bounds(1:end-1)
        xline(ax,b,':','Color',[.6,.6,.6]);
    end
    xlim(ax,[x_cat(1),x_cat(end)])
    title(ax,sprintf('fly %d (%s)',fly_list(f),fly_food(f)),'FontSize',7)
    set(ax,'FontSize',6)
end

xlabel(t_trace,'time (s)')
title(t_trace,'open-loop speed (blue, left axis, rad/s) and PB activity dff\_tot (orange, right axis, dF/F), per fly','Interpreter','none')

%% figure: mean dff_tot vs. open-loop speed, one subplot per fly (recreates optic_flow_analysis.m's figure 5, per fly instead of per trial)
% optic_flow_analysis.m's own figure(5) (in its "%% plot" section) is this
% same accumarray-mean-dF/F-by-speed-bin scatter, one subplot per TRIAL,
% restricted to instantaneous stillness (data.r_speed{i}<0.1) with no
% minimum-sample-count check. Reusing fly_activity_by_speed_still
% (already computed above, restricted to the sustained/causal stillness
% criterion and the well-powered-bin check adopted for the
% rotation-restricted figure) gives the same "mean dF/F per speed level"
% quantity, pooled per FLY instead of per trial -- same dotted reference
% line at that fly's own speed=0 value as the original figure.
figure(5); clf
set(gcf,'Name','mean dff_tot vs. open-loop speed, per fly','Position',[50,50,1600,900])
t5 = tiledlayout(n_rows,n_cols,'TileSpacing','compact','Padding','compact');
for f = 1:n_flies
    ax = nexttile(t5); hold(ax,'on')
    y = fly_activity_by_speed_still(f,:);
    scatter(ax,speed_levels,y,15,'filled')
    if ~isnan(y(1))
        plot(ax,speed_levels([1,end]),y(1)*[1,1],'k:')
    end
    title(ax,sprintf('fly %d (%s)',fly_list(f),fly_food(f)),'FontSize',7)
    set(ax,'FontSize',6)
end
xlabel(t5,'open-loop speed (rad/s)')
ylabel(t5,'mean dF/F (dff\_tot)','Interpreter','none')
title(t5,'mean dff\_tot vs. open-loop speed, per fly (dotted line = that fly''s own speed=0 value; rotation-restricted, well-powered bins only)','Interpreter','none')

%% figure: mean dff_tot vs. order of speed presentation, one subplot per fly (recreates optic_flow_analysis.m's figure 6)
% same idea as figure 6 in optic_flow_analysis.m's "%% plot" section:
% tests whether dF/F tracks WHEN a speed was presented (its order within
% the session) rather than the speed itself -- a monotonic drift across
% the session would show up here as a trend even if the true speed-tuning
% curve (figure 5) is flat. "order" is which of that fly's own present
% speed levels appeared FIRST in time, pooling the fly's trials in their
% own concatenated time order (same concatenation as the per-fly trace
% figure above) -- matching the original's [~,ia]=unique(...);
% [~,tmp_order]=sort(ia) logic, just pooled per fly instead of per trial.
figure(6); clf
set(gcf,'Name','mean dff_tot vs. order of open-loop speed presentation, per fly','Position',[50,50,1600,900])
t6 = tiledlayout(n_rows,n_cols,'TileSpacing','compact','Padding','compact');
for f = 1:n_flies
    rows_f = find(fly_num==fly_list(f));
    [~,speed_cat] = concat_fly_trace(data.xf(rows_f),ol_speed(rows_f),mean_activity(rows_f));

    present = find(~isnan(fly_activity_by_speed_still(f,:)));
    first_seen = nan(size(present));
    for k = 1:numel(present)
        first_seen(k) = find(abs(speed_cat-speed_levels(present(k)))<speed_tol,1,'first');
    end
    [~,ord] = sort(first_seen);
    present_ordered = present(ord);
    y = fly_activity_by_speed_still(f,present_ordered);

    ax = nexttile(t6); hold(ax,'on')
    scatter(ax,1:numel(y),y,15,'filled')
    zero_speed_pos = find(speed_levels(present_ordered)==0,1);
    if ~isempty(zero_speed_pos)
        plot(ax,[1,max(numel(y),2)],y(zero_speed_pos)*[1,1],'k:')
    end
    title(ax,sprintf('fly %d (%s)',fly_list(f),fly_food(f)),'FontSize',7)
    set(ax,'FontSize',6)
end
xlabel(t6,'order of open-loop speed presentation (1st, 2nd, ... distinct speed encountered)')
ylabel(t6,'mean dF/F (dff\_tot)','Interpreter','none')
title(t6,'mean dff\_tot vs. order of speed presentation, per fly (dotted line = that fly''s own speed=0 value)','Interpreter','none')

%% figure: mean dff_tot vs. open-loop speed, all flies overlaid, one panel per rearing food
% same fly_activity_by_speed_still dots as figure 5, but instead of one
% subplot per fly, every fly in a food group is overlaid on one shared
% panel -- a thin gray line connects a given fly's own dots across speed
% (breaking automatically at NaN/underpowered bins), with the dots
% themselves colored by food group so individual flies stay visually
% distinct from the group's overall trend. One panel per food condition
% (linked y-axes) rather than mixing both groups' dots on one panel, so
% each food's own scatter isn't obscured by the other's.
food_groups = {'m','g'};
food_labels = {'molasses','german'};
food_colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980];

figure(7); clf
set(gcf,'Name','mean dff_tot vs. open-loop speed, all flies overlaid, by food','Position',[100,100,1200,500])
axs7 = gobjects(1,numel(food_groups));
for gi = 1:numel(food_groups)
    in_grp = find(fly_food==food_groups{gi});
    axs7(gi) = subplot(1,numel(food_groups),gi); hold on
    for f = in_grp(:)'
        plot(speed_levels,fly_activity_by_speed_still(f,:),'-','Color',[.6,.6,.6],'LineWidth',0.5)
    end
    for f = in_grp(:)'
        scatter(speed_levels,fly_activity_by_speed_still(f,:),20,food_colors(gi,:),'filled')
    end
    grp_mean = mean(fly_activity_by_speed_still(in_grp,:),1,'omitnan');
    grp_sem  = std(fly_activity_by_speed_still(in_grp,:),1,'omitnan')./sqrt(sum(~isnan(fly_activity_by_speed_still(in_grp,:)),1));
    errorbar(speed_levels,grp_mean,grp_sem,'-o','Color',food_colors(gi,:)*.5,'LineWidth',2.5,'MarkerFaceColor',food_colors(gi,:)*.5)
    xlabel('open-loop speed (rad/s)')
    ylabel('mean dF/F (dff\_tot)','Interpreter','none')
    title(sprintf('%s (n=%d flies)',food_labels{gi},numel(in_grp)))
end
linkaxes(axs7,'xy')
sgtitle('mean dff\_tot vs. open-loop speed, all flies overlaid (gray=per-fly connection, dark=mean\pmSEM), by rearing food','Interpreter','none')

%% mean dff_tot per fly: open-loop speed = 0 vs. speed > 0 (moving), restricted to timepoints where the fly itself isn't moving
% recreates optic_flow_analysis.m's own dff_still/dff_flow comparison
% ("%% establish significance"), which collapsed all nonzero c_speed
% values into one "flow" bucket vs. a single "still" (c_speed==0)
% baseline -- but pooled per FLY (not per trial) and restricted using the
% same sustained/causal fly-stillness criterion and well-powered-bin
% check already established above (still_activity, min_n_bin), rather
% than the original's own (looser, differently-thresholded) definition.
min_n_bin_ol = 300; % samples -- same well-powered threshold as plot_activity_by_speed

fly_dff_ol0      = nan(n_flies,1); % mean dF/F while OL speed == 0
fly_dff_olmoving = nan(n_flies,1); % mean dF/F while OL speed > 0 (any nonzero speed, pooled)
for f = 1:n_flies
    rows = find(fly_num==fly_list(f));
    pooled_speed = cat(1,ol_speed{rows});
    pooled_act   = cat(1,still_activity{rows});

    is0       = abs(pooled_speed)<speed_tol & ~isnan(pooled_act);
    is_moving = pooled_speed>speed_tol & ~isnan(pooled_act);

    if sum(is0) >= min_n_bin_ol
        fly_dff_ol0(f) = mean(pooled_act(is0));
    end
    if sum(is_moving) >= min_n_bin_ol
        fly_dff_olmoving(f) = mean(pooled_act(is_moving));
    end
end

n_valid = sum(~isnan(fly_dff_ol0) & ~isnan(fly_dff_olmoving));
fprintf('\n%d/%d flies have well-powered OL=0 AND OL>0 bins (rotation-restricted)\n', n_valid, n_flies);

%% figure: mean dff_tot, OL speed=0 vs. moving, per fly, by rearing food
figure(8); clf
set(gcf,'Name','mean dff_tot: OL speed=0 vs. moving, per fly','Position',[100,100,900,500])
axs8 = gobjects(1,numel(food_groups));
for gi = 1:numel(food_groups)
    in_grp = find(fly_food==food_groups{gi});
    axs8(gi) = subplot(1,numel(food_groups),gi); hold on
    for f = in_grp(:)'
        if ~isnan(fly_dff_ol0(f)) && ~isnan(fly_dff_olmoving(f))
            plot([1,2],[fly_dff_ol0(f),fly_dff_olmoving(f)],'-','Color',[.6,.6,.6],'LineWidth',0.5)
        end
    end
    scatter(ones(numel(in_grp),1),fly_dff_ol0(in_grp),25,food_colors(gi,:),'filled')
    scatter(2*ones(numel(in_grp),1),fly_dff_olmoving(in_grp),25,food_colors(gi,:),'filled')

    grp_mean = [mean(fly_dff_ol0(in_grp),'omitnan'),mean(fly_dff_olmoving(in_grp),'omitnan')];
    grp_sem  = [std(fly_dff_ol0(in_grp),'omitnan'),std(fly_dff_olmoving(in_grp),'omitnan')] ./ ...
        sqrt([sum(~isnan(fly_dff_ol0(in_grp))),sum(~isnan(fly_dff_olmoving(in_grp)))]);
    errorbar([1,2],grp_mean,grp_sem,'-ok','LineWidth',2.5,'MarkerFaceColor','k')

    xlim([0.5,2.5]); xticks([1,2]); xticklabels({'OL speed = 0','OL speed > 0 (moving)'})
    ylabel('mean dF/F (dff\_tot)','Interpreter','none')
    title(sprintf('%s (n=%d flies)',food_labels{gi},numel(in_grp)))
end
linkaxes(axs8,'y')
sgtitle('mean dff\_tot per fly: OL speed=0 vs. moving (any nonzero speed), restricted to fly not turning','Interpreter','none')

%% save all figures as PDF
% same save-every-open-figure-by-its-Name approach as lpsp_kir_claude.m
fig_dir = 'C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\ugly_figures\optic flow';
if ~isfolder(fig_dir)
    mkdir(fig_dir)
end

fig_handles = findobj('Type','figure');
[~,order] = sort(arrayfun(@(f) f.Number, fig_handles));
fig_handles = fig_handles(order); % ascending by figure number, not creation/stacking order

save_failed = {};
for k = 1:numel(fig_handles)
    fig = fig_handles(k);
    fig_name = get(fig,'Name');
    if isempty(fig_name)
        fig_name = sprintf('figure_%d',fig.Number);
    end
    safe_name = regexprep(fig_name,'[^\w\-]+','_'); % filesystem-safe filename
    out_path  = fullfile(fig_dir,sprintf('fig%02d_%s.pdf',fig.Number,safe_name));
    try
        exportgraphics(fig,out_path)
        fprintf('saved %s\n', out_path);
    catch ME
        % e.g. the PDF is open in a viewer and locked for writing -- don't
        % let one such failure abort saving every other figure
        save_failed{end+1} = out_path; %#ok<AGROW>
        warning('could not save %s (%s) -- is it open in another program?', out_path, ME.message);
    end
end
if ~isempty(save_failed)
    fprintf('\n%d figure(s) failed to save -- close them in any viewer and re-run to update:\n', numel(save_failed));
    fprintf('  %s\n', save_failed{:});
end

%% functions
function [x_cat,speed_cat,act_cat,bounds] = concat_fly_trace(xf_cell,speed_cell,act_cell)
    % concatenate one or more trials' xf/speed/activity end-to-end along
    % elapsed time (each trial's own xf restarts near 0, so later trials
    % are offset by the running total), tracking each trial's own end time
    % in bounds for drawing boundary lines.
    x_cat = []; speed_cat = []; act_cat = []; bounds = [];
    x_offset = 0;
    for k = 1:numel(xf_cell)
        x = xf_cell{k}(:) + x_offset;
        speed_cat = [speed_cat; speed_cell{k}(:)]; %#ok<AGROW>
        act_cat   = [act_cat; act_cell{k}(:)]; %#ok<AGROW>
        x_cat     = [x_cat; x]; %#ok<AGROW>
        x_offset  = x(end);
        bounds(end+1) = x(end); %#ok<AGROW>
    end
end
function fly_activity_by_speed = plot_activity_by_speed(fig_offset,speed_cell,activity_cell,fly_num,fly_list,fly_food,speed_levels,speed_tol,n_flies,n_trials,title_suffix)
    % individual fly lines are colored by rearing food (fly_food, 'g'/'m'),
    % with one mean+/-SEM line per food group overlaid in the same color --
    % same "individual points/lines in a category's color + darker mean"
    % convention as lpsp_kir_claude.m's groupplot-style figures.
    %
    % a fly-speed bin is only kept if it's well powered: at least
    % min_n_bin samples went into it (checked AFTER excluding NaNs, so the
    % "restricted to rotation" call -- which NaNs out most samples in most
    % bins -- is held to the same standard as the unrestricted one, not
    % just "at least min_n_bin raw timepoints existed regardless of how
    % many survived filtering"). 300 samples (~5s at this dataset's ~60Hz
    % fictrac rate) is the default -- comfortably above fluorescence's own
    % autocorrelation timescale, so a kept bin isn't just a handful of
    % highly-correlated adjacent frames masquerading as independent
    % samples. Underpowered bins are left as NaN (dropped from that fly's
    % line and from the group mean/SEM at that speed, via 'omitnan'), not
    % zeroed or interpolated.
    if nargin < 11
        title_suffix = '';
    end
    min_n_bin = 300; % samples -- ~5s at ~60Hz; "well powered" cutoff, adjustable

    n_speed = numel(speed_levels);
    fly_activity_by_speed = nan(n_flies,n_speed);
    n_dropped = 0;
    for f = 1:n_flies
        rows = find(fly_num==fly_list(f));
        pooled_speed = cat(1,speed_cell{rows});
        pooled_act   = cat(1,activity_cell{rows});
        for s = 1:n_speed
            in_bin = abs(pooled_speed-speed_levels(s))<speed_tol & ~isnan(pooled_act);
            if sum(in_bin) >= min_n_bin
                fly_activity_by_speed(f,s) = mean(pooled_act(in_bin));
            else
                n_dropped = n_dropped+1;
            end
        end
    end
    fprintf('  %d/%d fly-speed bins dropped for being underpowered (<%d samples)%s\n', ...
        n_dropped, n_flies*n_speed, min_n_bin, title_suffix);

    food_groups = {'m','g'};
    food_labels = {'molasses','german'};
    food_colors = [0,0.4470,0.7410; 0.8500,0.3250,0.0980];

    figure(1+fig_offset); clf; hold on
    set(gcf,'Name',sprintf('mean PB activity vs open-loop speed%s',title_suffix))
    h_grp = gobjects(1,numel(food_groups));
    for gi = 1:numel(food_groups)
        in_grp = fly_food==food_groups{gi};
        if ~any(in_grp)
            continue
        end
        plot(speed_levels,fly_activity_by_speed(in_grp,:)','Color',[food_colors(gi,:),.35],'LineWidth',1)
        grp_mean = mean(fly_activity_by_speed(in_grp,:),1,'omitnan');
        grp_sem  = std(fly_activity_by_speed(in_grp,:),1,'omitnan')./sqrt(sum(~isnan(fly_activity_by_speed(in_grp,:)),1));
        h_grp(gi) = errorbar(speed_levels,grp_mean,grp_sem,'-o','Color',food_colors(gi,:),'LineWidth',2.5,'MarkerFaceColor',food_colors(gi,:));
    end
    legend(h_grp,arrayfun(@(gi) sprintf('%s (n=%d)',food_labels{gi},sum(fly_food==food_groups{gi})),1:numel(food_groups),'UniformOutput',false), ...
        'Location','best')
    xlabel('open-loop speed (rad/s)')
    ylabel('mean PB activity (dF/F)')
    title(sprintf('mean PB activity vs. open-loop speed%s, by rearing food (n=%d flies, %d trials)',title_suffix,n_flies,n_trials))
end
