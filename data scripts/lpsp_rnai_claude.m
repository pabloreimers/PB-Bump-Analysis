%% lpsp_rnai_claude
% Loads the joint LPsP-RNAi dataset (lpsp_rnai_joint_nosmooth_*.mat) and
% counts how many flies belong to each of the 5 genotypes in the
% experiment:
%   lpsp>TH-RNAi,     empty>TH-RNAi     (dopamine synthesis knockdown)
%   lpsp>vGlut-RNAi,  empty>vGlut-RNAi  (glutamate release knockdown)
%   lpsp>mCherry-RNAi                  (negative-control transgene, LPsP driver only -- there is no empty>mCherry-RNAi in this dataset)
%
% Genotype is read from all_data.meta, which is a raw-data path of the
% form (confirmed directly against every one of the 559 trials in
% lpsp_rnai_joint_nosmooth_20260522.mat, all 9 path components deep):
%   Z:\pablo\lpsp_rnai\todo\<project>\<date>\fly N\<trial folder>\registration_001
%
% <project> is NOT the genotype -- it's just which imaging session this
% trial was collected in ("th" or "lpsp_vglutrnai"), and each project
% folder mixes its own experimental trials together with their shared
% negative-control trials. The actual genotype lives in <trial folder>
% (e.g. "20250826-1_epg_syt8s_lpsp_thrnai"), except for a block of
% "blinded" trials (dates 20260120-20260212) where <trial folder> is
% anonymized ("..._blind_N") and the genotype is instead written into the
% "fly N" folder segment itself (e.g. "fly 2 _ empty _ vglut") -- so both
% segments are checked together, not just the trial folder.
%
% Within the "th" project, a genotype-keyword search alone is ambiguous:
% some TH-RNAi trials are named explicitly ("..._lpsp_thrnai"/
% "..._empty_thrnai"), but many others are named with NO target keyword
% at all ("..._lpsp_rnai"/"..._empty_rnai"). This matched the convention
% used by the older data scripts/lpsp_rnai_joint.m (loading a th-specific
% source file, it flagged only contains(geno,'mcherry'/'empty'/'vglut')
% and left every remaining trial as the implicit default = TH-RNAi).
% Checked directly here: the "th" project has exactly 4 keyword
% combinations (lpsp+bare, empty+bare, lpsp+thrnai, empty+thrnai) and NO
% mcherry-labeled trials at all, while "lpsp_vglutrnai" has exactly 3
% (empty+vglut, lpsp+vglut, lpsp+mcherry) and NO bare/unlabeled trials --
% so defaulting unlabeled "th"-project trials to TH-RNAi (not mCherry) is
% what makes every trial in the dataset resolve to one of exactly the 5
% genotypes above, with zero unlabeled and zero genotype-inconsistent
% flies.

%% load data
data_dir    = fullfile('data');
source_file = 'lpsp_rnai_joint_nosmooth_20260522.mat';

tmp = load(fullfile(data_dir,source_file),'all_data');
all_data = tmp.all_data(:);
n_trials = numel(all_data);
fprintf('loaded %d trials from %s\n', n_trials, source_file);

%% genotype (driver x target) per trial, from all_data.meta
driver = cell(n_trials,1); % 'lpsp' | 'empty'
target = cell(n_trials,1); % 'th' | 'vglut' | 'mcherry'

for i = 1:n_trials
    [fly_seg,trial_seg] = meta_fly_and_trial_seg(all_data(i).meta);
    combo = lower(regexprep([fly_seg,trial_seg],'[_\s]',''));

    if contains(combo,'lpsp')
        driver{i} = 'lpsp';
    elseif contains(combo,'empty')
        driver{i} = 'empty';
    else
        driver{i} = '';
        warning('trial %d: could not determine driver from meta "%s"', i, all_data(i).meta);
    end

    if contains(combo,'thrnai')
        target{i} = 'th';
    elseif contains(combo,'vglut')
        target{i} = 'vglut';
    elseif contains(combo,'mcherry')
        target{i} = 'mcherry';
    else
        target{i} = 'th'; % bare "..._rnai" in the "th" project = implicit TH-RNAi (see header)
    end
end
genotype = strcat(driver,'>',target);

%% light condition per trial (closed loop vs. dark), from ft.pattern
is_dark = false(n_trials,1);
for i = 1:n_trials
    is_dark(i) = contains(all_data(i).ft.pattern,'background');
end
fprintf('\ntrials: %d closed-loop, %d dark\n', sum(~is_dark), sum(is_dark));

%% fly ID per trial (project + date + "fly N" folder), and a genotype-consistency check
fly_id = cell(n_trials,1);
for i = 1:n_trials
    fly_id{i} = meta_fly_id(all_data(i).meta);
end
[fly_list,~,fly_num] = unique(fly_id);
n_flies = numel(fly_list);

fly_genotype = cell(n_flies,1);
for f = 1:n_flies
    g = unique(genotype(fly_num==f));
    if numel(g) > 1
        warning('fly %s has inconsistent genotype labels across its trials: %s', fly_list{f}, strjoin(g,', '));
    end
    fly_genotype{f} = g{1};
end
fprintf('%d trials -> %d flies\n', n_trials, n_flies);

%% report: flies per genotype
geno_order = {'lpsp>th','empty>th','lpsp>vglut','empty>vglut','lpsp>mcherry'};

fprintf('\n=== flies per genotype ===\n');
for k = 1:numel(geno_order)
    n_g = sum(strcmp(fly_genotype,geno_order{k}));
    fprintf('  %-14s n=%d flies\n', geno_order{k}, n_g);
end

other_geno = setdiff(unique(fly_genotype),geno_order);
for k = 1:numel(other_geno)
    fprintf('  %-14s n=%d flies  (unexpected genotype -- check meta parsing)\n', other_geno{k}, sum(strcmp(fly_genotype,other_geno{k})));
end

fprintf('\ntotal flies: %d\n', n_flies);

%% bump mobility: ratio of bump path length to heading path length, during walking bouts
% same bump-mobility ("path-length gain") calculation as
% data scripts/lpsp_p2x2_walking_script.m's own section 2 ("how much does
% the bump (mu) move relative to how much the fly turns?"): within each
% sustained-turning "walking bout" (|r_speed|, gaussian-smoothed over
% smooth_window samples, above turn_thresh, with gaps up to
% max_gap_frames bridged via imclose and bouts shorter than
% min_walking_frames dropped via bwareaopen), the bump's path length (sum
% |diff(smoothed mu)|) and the fly's own turning path length (integral of
% |smoothed r_speed| dt) are summed per bout, then related by a
% bout-duration-weighted (w=sqrt(dur)), through-origin least-squares fit:
% mov_ratio = bump path length per unit heading path length. Same
% thresholds too (smooth_window=60 samples/~1s, turn_thresh=.25rad/s,
% max_gap_frames=min_walking_frames=30 frames/0.5s). trial_walking_bouts
% below is copied verbatim from lpsp_kir_claude.m, which itself copied
% this bout-detection algorithm from lpsp_p2x2_walking_script.m unchanged.
%
% The one difference from lpsp_p2x2_walking_script.m: that script fits
% mov_ratio once PER TRIAL. Here (as in lpsp_kir_claude.m) a fly's bouts
% are pooled across every trial it has in a given light condition BEFORE
% fitting, so each fly contributes exactly one mov_ratio per light
% condition regardless of how many trials it has in that condition --
% matching lpsp_kir_claude.m's own per-fly aggregation, not
% lpsp_p2x2_walking_script.m's per-trial one.
mobility_smooth_window      = 60;   % samples (~1s at this rig's ~60Hz fictrac rate), gaussian smoothing window
mobility_turn_thresh        = 0.25; % rad/s, minimum heading speed to call a bout "walking"
mobility_max_gap_frames     = 30;   % frames of non-walking allowed within a bout before splitting it (0.5s)
mobility_min_walking_frames = 30;   % minimum bout length to keep (0.5s)

chunk_bout_mu        = cell(n_trials,1);
chunk_bout_speed     = cell(n_trials,1);
chunk_bout_dur       = cell(n_trials,1);
chunk_mu_smooth      = cell(n_trials,1); % full-trial smoothed bump position (xf timebase, unwrapped) -- for the diagnostic overlay below
chunk_heading_smooth = cell(n_trials,1); % full-trial smoothed heading position (xf timebase, unwrapped) -- ditto
chunk_is_walking     = cell(n_trials,1); % full-trial walking-bout membership (xf timebase) -- for shading the diagnostic overlay
for i = 1:n_trials
    [chunk_bout_mu{i},chunk_bout_speed{i},chunk_bout_dur{i},chunk_mu_smooth{i},chunk_heading_smooth{i},chunk_is_walking{i}] = trial_walking_bouts( ...
        all_data(i),mobility_smooth_window,mobility_turn_thresh,mobility_max_gap_frames,mobility_min_walking_frames);
end

%% per-fly, per-light-condition bump mobility ratio
cond_label = {'closed loop','dark'};

group_defs = struct('geno',{},'dark',{},'label',{});
for gi = 1:numel(geno_order)
    for ci = 1:numel(cond_label)
        group_defs(end+1) = struct('geno',geno_order{gi},'dark',ci-1,'label',sprintf('%s (%s)',geno_order{gi},cond_label{ci})); %#ok<SAGROW>
    end
end

% require a minimum number of walking bouts for a given FLY, WITHIN ONE
% LIGHT CONDITION, before that fly contributes its own point to the group
% plot at all. This is a PER-FLY threshold, not a group-level one -- a
% group is no longer included/excluded by summing bouts across every fly
% that shares its genotype x light-condition; each fly is judged only
% against its own bout count. (An earlier version of this section did
% exactly that group-level sum, which let individual flies with only a
% handful of bouts each into the plot as long as the group's pooled total
% cleared the threshold -- that's the bug being fixed here.)
min_bouts_per_fly = 10;

fly_mob_ratio  = []; % one entry per (fly, light-condition) pair with >=min_bouts_per_fly bouts of its OWN
cat_x_mob      = []; % group_defs index for that entry
fly_mob_fly    = []; % fly_num index for that entry
fly_n_bouts    = []; % this fly's own bout count for that (fly, light-condition) pair
group_rows     = cell(numel(group_defs),1);
for gd = 1:numel(group_defs)
    group_rows{gd} = find(strcmp(genotype,group_defs(gd).geno) & is_dark==group_defs(gd).dark);
end

for gd = 1:numel(group_defs)
    rows = group_rows{gd};
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mov_mu    = cat(1,chunk_bout_mu{trial_list});
        mov_speed = cat(1,chunk_bout_speed{trial_list});
        dur       = cat(1,chunk_bout_dur{trial_list});

        if numel(mov_mu) >= min_bouts_per_fly
            w = sqrt(dur);
            ratio = (w.*mov_speed) \ (w.*mov_mu); % weighted, through-origin: bump path length per unit heading path length
            fly_mob_ratio(end+1) = ratio; %#ok<AGROW>
            cat_x_mob(end+1)     = gd; %#ok<AGROW>
            fly_mob_fly(end+1)   = f; %#ok<AGROW>
            fly_n_bouts(end+1)   = numel(mov_mu); %#ok<AGROW>
        end
    end
end

% a group is shown in the plot if AT LEAST ONE fly qualified for it above
% -- not based on any group-level bout sum.
keep_group = ismember(1:numel(group_defs),cat_x_mob)';

fprintf('\n=== flies per genotype x light-condition group with >=%d walking bouts of their own (threshold applies PER FLY) ===\n', min_bouts_per_fly);
for gd = 1:numel(group_defs)
    n_flies_total = numel(unique(fly_num(group_rows{gd})));
    n_flies_pass  = sum(cat_x_mob==gd);
    excluded_tag = '';
    if ~keep_group(gd)
        excluded_tag = '  EXCLUDED (no fly meets threshold)';
    end
    fprintf('  %-28s %d/%d flies qualify%s\n', group_defs(gd).label, n_flies_pass, n_flies_total, excluded_tag);
end

%% figure: bump mobility (path-length ratio), one point per fly, by genotype x light condition
base_colors    = lines(numel(geno_order));
all_cat_colors = zeros(numel(group_defs),3);
for k = 1:numel(geno_order)
    all_cat_colors(2*k-1,:) = base_colors(k,:);
    all_cat_colors(2*k,:)   = base_colors(k,:);
end

% renumber categories to only the kept (above-threshold) groups, so
% groupplot's contiguous 1:n_cat indexing lines up, and drop any fly
% entries belonging to an excluded group from the plotted data
kept_idx                = find(keep_group);
gd_to_plot_idx           = zeros(numel(group_defs),1);
gd_to_plot_idx(kept_idx) = 1:numel(kept_idx);
plot_mask                = keep_group(cat_x_mob);
cat_x_mob_plot           = gd_to_plot_idx(cat_x_mob(plot_mask));
fly_mob_ratio_plot       = fly_mob_ratio(plot_mask);
cat_labels               = {group_defs(kept_idx).label};
cat_colors               = all_cat_colors(kept_idx,:);

figure(1); clf
set(gcf,'Name','bump mobility: path-length ratio per fly','Position',[100,100,max(700,120*numel(cat_labels)),600])
groupplot(cat_x_mob_plot,fly_mob_ratio_plot,cat_labels,cat_colors)
hold on
plot(xlim,[1,1],':k') % reference: bump moves exactly as much as heading
ylabel('bump path length / heading path length (per fly, bout-duration-weighted)')
title(sprintf('bump mobility relative to fly turning, by genotype and light condition (flies with \\geq%d bouts of their own)',min_bouts_per_fly))

%% diagnostic figures: PB activity (im.z) + bump position (im.mu) + heading, for every fly x light-condition with mobility ratio > 3
% a mov_ratio this high means the fitted bump path length is more than 3x
% the fly's own turning path length within its walking bouts -- almost
% certainly a bump-tracking artifact (the fitted mu jittering/spinning
% independent of real bump movement) rather than a real biological signal,
% so these are worth looking at directly. fly_mob_ratio already only
% contains flies that individually cleared min_bouts_per_fly, so no
% further bout-count filtering is needed here. One figure per qualifying
% (fly, light-condition) instance; if that fly has more than one trial in
% this light condition, they're concatenated frame-by-frame
% (concat_trials_im, same helper and convention as lpsp_kir_claude.m's own
% PB-activity-overlay figure).
mobility_diagnostic_thresh = 1.5;

diag_idx = find(fly_mob_ratio > mobility_diagnostic_thresh);
fprintf('\n=== %d fly x light-condition instance(s) with bump mobility ratio > %g ===\n', numel(diag_idx), mobility_diagnostic_thresh);

for d = 1:numel(diag_idx)
    idx  = diag_idx(d);
    gd   = cat_x_mob(idx);
    f    = fly_mob_fly(idx);
    trial_idx = find(fly_num==f & is_dark==group_defs(gd).dark)';

    [x_cat,alpha,z_cat,mu_cat,heading_cat,is_walking_cat,bounds] = concat_trials_im(all_data,chunk_mu_smooth,chunk_heading_smooth,chunk_is_walking,trial_idx);

    figure(200+d); clf
    set(gcf,'Name',sprintf('mobility diagnostic: %s, %s',fly_short_label(fly_list{f}),group_defs(gd).label),'Position',[100,100,1000,600])
    t = tiledlayout(3,1,'TileSpacing','compact');

    nexttile(t,[2,1]); hold on
    imagesc(x_cat,alpha,z_cat)
    set(gca,'YDir','normal','XLim',[x_cat(1),x_cat(end)],'YLim',[alpha(1),alpha(end)])

    % mov_ratio is a regression over ONLY the walking-bout stretches below
    % -- shade them so it's visually obvious whether bump/heading movement
    % elsewhere in the trial (which does NOT feed mov_ratio) is what's
    % driving any mismatch between "how much movement this panel shows"
    % and mov_ratio's own value.
    bout_edges = find(diff([0;is_walking_cat;0])~=0);
    for b = 1:2:numel(bout_edges)-1
        bout_x = x_cat([bout_edges(b),bout_edges(b+1)-1]);
        patch(bout_x([1,2,2,1]),alpha([1,1,end,end]),[1,1,0],'FaceAlpha',.15,'EdgeColor','none')
        xline(bout_x(1),'Color',[.9,.7,0],'LineWidth',1.2)
        xline(bout_x(2),'Color',[.9,.7,0],'LineWidth',1.2)
    end

    % mu/heading are single-hemisphere angles (-pi:pi), but alpha (and so
    % the imagesc's y-axis) spans -pi:~3pi -- both hemispheres stacked.
    % Drawn once as-is (lines up with the 1st hemisphere band) and once
    % shifted by +2*pi (lines up with the 2nd) so the overlay tracks the
    % bump/heading across the whole panel, not just its bottom half.
    mu_plot = mu_cat;
    mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan; % break the line at circular wraps, not connect across them
    plot(x_cat,mu_plot,'w','LineWidth',1)
    plot(x_cat,mu_plot+2*pi,'w','LineWidth',1)

    heading_plot = heading_cat;
    heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;
    plot(x_cat,heading_plot,'Color',[0,1,1],'LineWidth',1)
    plot(x_cat,heading_plot+2*pi,'Color',[0,1,1],'LineWidth',1)

    for b = bounds(1:end-1)
        xline(b,':','Color',[.6,.6,.6]);
    end

    cb = colorbar; ylabel(cb,'z-score')
    xlabel('frame'); ylabel('PB angle (rad)')
    title(sprintf('%s -- %s -- mov ratio=%.2f -- PB activity (im.z), %d-sample-smoothed bump position (white), heading (-ft.cue, cyan), walking bouts (yellow)', ...
        fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio(idx),mobility_smooth_window),'Interpreter','none')

    % the scatter + fit actually behind mov_ratio: one point per walking
    % bout, this fly's own bouts only, same weighted through-origin fit --
    % same diagnostic idea as lpsp_p2x2_walking_script.m's own "show the
    % exact scatter + fit behind mov_ratio for single trials" section.
    bout_mov_mu    = cat(1,chunk_bout_mu{trial_idx});
    bout_mov_speed = cat(1,chunk_bout_speed{trial_idx});
    bout_dur       = cat(1,chunk_bout_dur{trial_idx});

    nexttile(t); hold on
    scatter(bout_mov_speed,bout_mov_mu,20,bout_dur,'filled')
    xl = xlim; xl(1) = 0;
    plot(xl,fly_mob_ratio(idx)*xl,'-k')
    plot(xl,xl,':','Color',[.6,.6,.6])
    xlim(xl)
    cb2 = colorbar; ylabel(cb2,'bout duration (s)')
    xlabel('heading path length (rad)'); ylabel('bump path length (rad)')
    title(sprintf('%d bouts -- fit slope (mov ratio) = %.2f -- dotted line = slope 1',numel(bout_mov_mu),fly_mob_ratio(idx)))

    fprintf('  fly %-30s %-28s mov_ratio=%.2f (%d bouts)\n',fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio(idx),numel(bout_mov_mu));
end

%% re-estimate bump position from smoothed, z-scored raw fluorescence (im.f), as an alternative to the upstream im.mu fit
% im.mu (used above) comes from whatever fit produced the joint dataset --
% if it's noisier than it needs to be, that noise alone could inflate
% mov_ratio (spurious bump "movement" that isn't real). This recomputes
% bump position from scratch, directly off the raw per-glomerulus
% fluorescence (im.f, same 32-glomerulus x n_frames layout as im.z),
% heavily smoothed in time BEFORE z-scoring:
%   1) gaussian-smooth each glomerulus's raw fluorescence over
%      mobility_smooth_window_s seconds (converted to samples using this
%      trial's own imaging frame period -- imaging frame rate isn't
%      constant across trials in this dataset, e.g. 37800 vs 36006
%      fictrac samples for the same 6306 imaging frames)
%   2) z-score each glomerulus's smoothed trace over time (puts every
%      glomerulus on the same scale, matching what im.z represents)
% then bump position/strength are the angle/length of the mean resultant
% vector of the 32 glomeruli's (smoothed+zscored) activity, each treated
% as a vector pointing at that glomerulus's own angular position
% (im.alpha) -- the same population-vector-average formula used
% throughout this codebase's legacy bump_calc_pb (e.g.
% epg_grab_script.m's own alpha/pol2cart/cart2pol sequence), just applied
% to this new smoothed+zscored signal instead of the original dff_cluster.
%
% mobility_smooth_window_s is reused below (bump mobility using the new
% bump-position estimate) as the SAME real-time smoothing window applied
% to the heading (ft.r_speed) trace, so the two path lengths (bump vs.
% heading) that make up mov_ratio are smoothed at a matched timescale
% rather than the bump alone being heavily smoothed -- im.f (imaging
% timebase, xb) and ft.r_speed (fictrac timebase, xf) run at different,
% and not perfectly identical across trials, sample rates, so
% seconds->samples is converted separately for each using that signal's
% own timebase, never a shared raw sample count.
mobility_smooth_window_s = 10; % seconds

chunk_mu_new  = cell(n_trials,1); % bump position, this new estimate (im timebase)
chunk_rho_new = cell(n_trials,1); % bump vector strength, this new estimate
chunk_fz      = cell(n_trials,1); % smoothed + per-glomerulus zscored fluorescence (32 x n_im), for the diagnostic imagesc below
for i = 1:n_trials
    [chunk_mu_new{i},chunk_rho_new{i},chunk_fz{i}] = trial_smoothed_bump_pva(all_data(i),mobility_smooth_window_s);
end

%% bump mobility using the new bump-position estimate
% same walking-bout DETECTION (and therefore the same bout boundaries and
% duration as the original analysis above -- bout detection depends only
% on ft.r_speed, never on mu) but now the heading trace feeding both bout
% detection and the per-bout heading path length is smoothed over the
% SAME mobility_smooth_window_s used to derive mu_new above, instead
% of the original ~1s mobility_smooth_window -- so bump and heading path
% lengths are computed at matched smoothing timescales. mov_speed/dur are
% NOT reused from chunk_bout_speed/chunk_bout_dur here (unlike an earlier
% version of this section) precisely because that smoothing differs now.
chunk_bout_mu_new        = cell(n_trials,1);
chunk_bout_speed_new     = cell(n_trials,1);
chunk_bout_dur_new       = cell(n_trials,1);
chunk_mu_smooth_new      = cell(n_trials,1); % full-trial smoothed bump position (xf timebase, unwrapped) -- for the diagnostic overlay below
chunk_heading_smooth_new = cell(n_trials,1); % full-trial smoothed heading position (xf timebase, unwrapped) -- ditto
chunk_is_walking_new     = cell(n_trials,1); % full-trial walking-bout membership (xf timebase) -- for shading the diagnostic overlay
for i = 1:n_trials
    [chunk_bout_mu_new{i},chunk_bout_speed_new{i},chunk_bout_dur_new{i},chunk_mu_smooth_new{i},chunk_heading_smooth_new{i},chunk_is_walking_new{i}] = trial_walking_bouts_from_mu(all_data(i),chunk_mu_new{i}, ...
        mobility_smooth_window_s,mobility_turn_thresh,mobility_max_gap_frames,mobility_min_walking_frames);
end

% group_rows (which trials belong to which genotype x light-condition
% group) is reused unchanged from above -- that's a genotype/light-only
% grouping, independent of bump or bout method. min_bouts_per_fly is
% reused unchanged too (still 50, still applied PER FLY, not summed
% across a group).
fly_mob_ratio_new = [];
cat_x_mob_new     = [];
fly_mob_fly_new   = [];
fly_n_bouts_new   = []; % this fly's OWN bout count for that (fly, light-condition) pair
for gd = 1:numel(group_defs)
    rows = group_rows{gd};
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mov_mu    = cat(1,chunk_bout_mu_new{trial_list});
        mov_speed = cat(1,chunk_bout_speed_new{trial_list});
        dur       = cat(1,chunk_bout_dur_new{trial_list});

        if numel(mov_mu) >= min_bouts_per_fly
            w = sqrt(dur);
            ratio = (w.*mov_speed) \ (w.*mov_mu);
            fly_mob_ratio_new(end+1) = ratio; %#ok<AGROW>
            cat_x_mob_new(end+1)     = gd; %#ok<AGROW>
            fly_mob_fly_new(end+1)   = f; %#ok<AGROW>
            fly_n_bouts_new(end+1)   = numel(mov_mu); %#ok<AGROW>
        end
    end
end

% a group is shown if AT LEAST ONE fly qualified for it above -- not
% based on any group-level bout sum.
keep_group_new = ismember(1:numel(group_defs),cat_x_mob_new)';

fprintf('\n=== flies per genotype x light-condition group with >=%d walking bouts of their own, NEW estimate ===\n', min_bouts_per_fly);
for gd = 1:numel(group_defs)
    n_flies_total = numel(unique(fly_num(group_rows{gd})));
    n_flies_pass  = sum(cat_x_mob_new==gd);
    excluded_tag = '';
    if ~keep_group_new(gd)
        excluded_tag = '  EXCLUDED (no fly meets threshold)';
    end
    fprintf('  %-28s %d/%d flies qualify%s\n', group_defs(gd).label, n_flies_pass, n_flies_total, excluded_tag);
end

%% figure: bump mobility (new estimate), one point per fly, by genotype x light condition
kept_idx_new                 = find(keep_group_new);
gd_to_plot_idx_new            = zeros(numel(group_defs),1);
gd_to_plot_idx_new(kept_idx_new) = 1:numel(kept_idx_new);
plot_mask_new                = keep_group_new(cat_x_mob_new);
cat_x_mob_new_plot            = gd_to_plot_idx_new(cat_x_mob_new(plot_mask_new));
fly_mob_ratio_new_plot        = fly_mob_ratio_new(plot_mask_new);
cat_labels_new                = {group_defs(kept_idx_new).label};
cat_colors_new                = all_cat_colors(kept_idx_new,:);

figure(2); clf
set(gcf,'Name','bump mobility (new mu estimate): path-length ratio per fly','Position',[100,100,max(700,120*numel(cat_labels_new)),600])
groupplot(cat_x_mob_new_plot,fly_mob_ratio_new_plot,cat_labels_new,cat_colors_new)
hold on
plot(xlim,[1,1],':k')
ylabel('bump path length / heading path length (per fly, bout-duration-weighted)')
title(sprintf('bump mobility (bump position re-estimated from %ds-smoothed, zscored im.f; heading smoothed over the same %ds), by genotype and light condition (flies with \\geq%d bouts of their own)', ...
    mobility_smooth_window_s,mobility_smooth_window_s,min_bouts_per_fly))

%% diagnostic figures with the new estimate: smoothed+zscored im.f + new bump position + heading, for every fly x light-condition with new mobility ratio > 3
% fly_mob_ratio_new already only contains flies that individually cleared
% min_bouts_per_fly, so no further bout-count filtering is needed here.
diag_idx_new = find(fly_mob_ratio_new > .9);
fprintf('\n=== %d fly x light-condition instance(s) with NEW bump mobility ratio > %g ===\n', numel(diag_idx_new), mobility_diagnostic_thresh);

for d = 1:numel(diag_idx_new)
    idx  = diag_idx_new(d);
    gd   = cat_x_mob_new(idx);
    f    = fly_mob_fly_new(idx);
    trial_idx = find(fly_num==f & is_dark==group_defs(gd).dark)';

    [x_cat,alpha,fz_cat,mu_cat,heading_cat,is_walking_cat,bounds] = concat_trials_im_new(all_data,chunk_fz,chunk_mu_smooth_new,chunk_heading_smooth_new,chunk_is_walking_new,trial_idx);

    figure(300+d); clf
    set(gcf,'Name',sprintf('mobility diagnostic (new mu): %s, %s',fly_short_label(fly_list{f}),group_defs(gd).label),'Position',[100,100,1000,600])
    t = tiledlayout(3,1,'TileSpacing','compact');

    nexttile(t,[2,1]); hold on
    imagesc(x_cat,alpha,fz_cat)
    set(gca,'YDir','normal','XLim',[x_cat(1),x_cat(end)],'YLim',[alpha(1),alpha(end)])

    % mov_ratio is a regression over ONLY the walking-bout stretches below
    % -- shade them so it's visually obvious whether bump/heading movement
    % elsewhere in the trial (which does NOT feed mov_ratio) is what's
    % driving any mismatch between "how much movement this panel shows"
    % and mov_ratio's own value.
    bout_edges = find(diff([0;is_walking_cat;0])~=0);
    for b = 1:2:numel(bout_edges)-1
        bout_x = x_cat([bout_edges(b),bout_edges(b+1)-1]);
        patch(bout_x([1,2,2,1]),alpha([1,1,end,end]),[1,1,0],'FaceAlpha',.15,'EdgeColor','none')
        xline(bout_x(1),'Color',[.9,.7,0],'LineWidth',1.2)
        xline(bout_x(2),'Color',[.9,.7,0],'LineWidth',1.2)
    end

    % see the comment on the equivalent overlay in the first diagnostic
    % loop above: mu/heading are -pi:pi, alpha is -pi:~3pi (both
    % hemispheres), so each trace is drawn twice, once per hemisphere band
    mu_plot = mu_cat;
    mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan;
    plot(x_cat,mu_plot,'w','LineWidth',1)
    plot(x_cat,mu_plot+2*pi,'w','LineWidth',1)

    heading_plot = heading_cat;
    heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;
    plot(x_cat,heading_plot,'Color',[0,1,1],'LineWidth',1)
    plot(x_cat,heading_plot+2*pi,'Color',[0,1,1],'LineWidth',1)

    for b = bounds(1:end-1)
        xline(b,':','Color',[.6,.6,.6]);
    end

    cb = colorbar; ylabel(cb,'z-score')
    xlabel('frame'); ylabel('PB angle (rad)')
    title(sprintf('%s -- %s -- mov ratio=%.2f -- %ds-smoothed zscored im.f, %ds-smoothed bump position (white), %ds-smoothed heading (-ft.cue, cyan), walking bouts (yellow)', ...
        fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio_new(idx),mobility_smooth_window_s,mobility_smooth_window_s,mobility_smooth_window_s),'Interpreter','none')

    % the scatter + fit actually behind mov_ratio: one point per walking
    % bout, this fly's own bouts only, same weighted through-origin fit.
    bout_mov_mu    = cat(1,chunk_bout_mu_new{trial_idx});
    bout_mov_speed = cat(1,chunk_bout_speed_new{trial_idx});
    bout_dur       = cat(1,chunk_bout_dur_new{trial_idx});

    nexttile(t); hold on
    scatter(bout_mov_speed,bout_mov_mu,20,bout_dur,'filled')
    xl = xlim; xl(1) = 0;
    plot(xl,fly_mob_ratio_new(idx)*xl,'-k')
    plot(xl,xl,':','Color',[.6,.6,.6])
    xlim(xl)
    cb2 = colorbar; ylabel(cb2,'bout duration (s)')
    xlabel('heading path length (rad)'); ylabel('bump path length (rad)')
    title(sprintf('%d bouts -- fit slope (mov ratio) = %.2f -- dotted line = slope 1',numel(bout_mov_mu),fly_mob_ratio_new(idx)))

    fprintf('  fly %-30s %-28s mov_ratio=%.2f (%d bouts)\n',fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio_new(idx),numel(bout_mov_mu));
end

%% walking bouts redefined by TOTAL movement (forward + rotational), not rotation alone
% the bout detector above only looks at |r_speed| -- so a fly walking
% briskly straight forward while barely turning can drop BELOW turn_thresh
% and get treated as "not walking", fragmenting one continuous bout of
% forward locomotion into several. This section instead builds a single
% "how much is the fly moving, period" signal from BOTH ft.f_speed
% (forward/backward walking speed) and ft.r_speed (rotational/turning
% speed), and thresholds THAT for bout detection.
%
% Why this isn't a completely well-posed problem, and what I did about it:
% f_speed and r_speed are physically different quantities -- f_speed is a
% TRANSLATIONAL speed (this dataset's units, confirmed directly: signed,
% ~50% negative, i.e. forward/backward, with typical |f_speed| around
% 0.1-0.7 and occasional much larger spikes), while r_speed is an ANGULAR
% velocity (rad/s, signed left/right, typical |r_speed| around 0.02-0.4).
% There's no first-principles "correct" way to add a translational speed
% to an angular velocity -- any combination requires picking a conversion
% factor, which is really an implicit choice of "reference radius" (the
% distance from some pivot at which r_speed's rotation would produce a
% matching tangential speed). Three ways to pick that factor, roughly in
% order of how principled vs. how simple they are:
%   1) PHYSICAL: multiply r_speed by an assumed radius in mm (e.g. the
%      fly's own body length, ~2-2.5mm, or half that if you think of
%      rotation as pivoting near the body center) to get an actual
%      mm/s-equivalent tangential speed, directly addable to f_speed.
%      Most physically interpretable, but "which radius" is itself a
%      modeling choice with no single right answer for what should count
%      as equivalent to translation.
%   2) STATISTICAL: normalize each channel by its own spread (e.g. divide
%      by its SD, or by its own 95th percentile) so both become
%      dimensionless "how many typical-magnitudes is this" values, then
%      combine. Sidesteps the unit problem entirely, but the result
%      depends on what population you normalize against (this trial? this
%      fly? the whole dataset?) and loses direct physical meaning.
%   3) PRECEDENT: data scripts/lpsp_rnai_joint.m (this same repo, an older
%      script) already defines a "was this fly walking" flag as
%      f_speed>0.1 OR abs(r_speed)*2>0.1 -- i.e. r_speed scaled by a
%      factor of 2 before comparing on the same scale as f_speed. Checked
%      directly against this dataset: the ratio of each channel's own 95th
%      percentile |f_speed|/|r_speed| across a few sampled trials ranges
%      from about 1.2x to 4.2x, centered roughly around 2-3x -- so that
%      prior factor of 2 isn't arbitrary, it's in the right ballpark
%      empirically too.
% I went with option 3 (r_speed_weight below), both because it's already
% this codebase's own precedent and because it roughly matches what a
% quick empirical check of this dataset's own speed distributions
% suggests. It is very much a free, arguable parameter -- if the resulting
% bouts still look wrong in the diagnostic figures below, this is the
% first number to change.
%
% The two channels are combined via their EUCLIDEAN magnitude
% (sqrt(f_speed^2 + (r_speed_weight*r_speed)^2)), not a plain sum --
% treating (f_speed, r_speed_weight*r_speed) as if they were two
% orthogonal velocity components and taking the resulting vector length.
% A plain sum would double-count a fly that's both walking forward AND
% turning hard at the same instant; the Euclidean form doesn't inflate
% past whichever single channel would already justify calling this
% "moving" on its own (it's always >= max of the two).
%
% mov_speed itself (the heading path length that mov_ratio is actually
% computed from) is UNCHANGED -- still purely from r_speed, exactly as
% above. Only which stretches of the trial count as "in a bout" changes
% here; what's measured within a bout (bump angle vs. heading angle) does
% not, so mov_ratio keeps its original "radians of bump per radian of
% heading turn" meaning.
mobility_r_speed_weight = 2;    % converts rad/s of r_speed onto roughly f_speed's own scale before combining -- see precedent/rationale above
mobility_move_thresh    = mobility_turn_thresh; % same starting value as the rotation-only threshold above (0.25) -- now applied to the COMBINED signal, so it means something different; retune by eye against the diagnostic figures below if bouts still look wrong

chunk_bout_mu_move        = cell(n_trials,1);
chunk_bout_speed_move     = cell(n_trials,1);
chunk_bout_dur_move       = cell(n_trials,1);
chunk_mu_smooth_move      = cell(n_trials,1);
chunk_heading_smooth_move = cell(n_trials,1);
chunk_is_walking_move     = cell(n_trials,1);
for i = 1:n_trials
    [chunk_bout_mu_move{i},chunk_bout_speed_move{i},chunk_bout_dur_move{i},chunk_mu_smooth_move{i},chunk_heading_smooth_move{i},chunk_is_walking_move{i}] = trial_walking_bouts_totalmovement(all_data(i),chunk_mu_new{i}, ...
        mobility_smooth_window_s,mobility_move_thresh,mobility_r_speed_weight,mobility_max_gap_frames,mobility_min_walking_frames);
end

% group_rows/min_bouts_per_fly reused unchanged -- see the equivalent
% comment in the new-mu-estimate section above.
fly_mob_ratio_move = [];
cat_x_mob_move     = [];
fly_mob_fly_move   = [];
fly_n_bouts_move   = [];
for gd = 1:numel(group_defs)
    rows = group_rows{gd};
    these_flies = unique(fly_num(rows));
    for f = these_flies'
        trial_list = rows(fly_num(rows)==f);
        mov_mu    = cat(1,chunk_bout_mu_move{trial_list});
        mov_speed = cat(1,chunk_bout_speed_move{trial_list});
        dur       = cat(1,chunk_bout_dur_move{trial_list});

        if numel(mov_mu) >= min_bouts_per_fly
            w = sqrt(dur);
            ratio = (w.*mov_speed) \ (w.*mov_mu);
            fly_mob_ratio_move(end+1) = ratio; %#ok<AGROW>
            cat_x_mob_move(end+1)     = gd; %#ok<AGROW>
            fly_mob_fly_move(end+1)   = f; %#ok<AGROW>
            fly_n_bouts_move(end+1)   = numel(mov_mu); %#ok<AGROW>
        end
    end
end

% a group is shown if AT LEAST ONE fly qualified for it above -- not
% based on any group-level bout sum.
keep_group_move = ismember(1:numel(group_defs),cat_x_mob_move)';

fprintf('\n=== flies per genotype x light-condition group with >=%d walking bouts of their own, TOTAL-MOVEMENT bout definition ===\n', min_bouts_per_fly);
for gd = 1:numel(group_defs)
    n_flies_total = numel(unique(fly_num(group_rows{gd})));
    n_flies_pass  = sum(cat_x_mob_move==gd);
    excluded_tag = '';
    if ~keep_group_move(gd)
        excluded_tag = '  EXCLUDED (no fly meets threshold)';
    end
    fprintf('  %-28s %d/%d flies qualify%s\n', group_defs(gd).label, n_flies_pass, n_flies_total, excluded_tag);
end

%% figure: bump mobility (total-movement bout definition), one point per fly, by genotype x light condition
kept_idx_move                  = find(keep_group_move);
gd_to_plot_idx_move             = zeros(numel(group_defs),1);
gd_to_plot_idx_move(kept_idx_move) = 1:numel(kept_idx_move);
plot_mask_move                 = keep_group_move(cat_x_mob_move);
cat_x_mob_move_plot             = gd_to_plot_idx_move(cat_x_mob_move(plot_mask_move));
fly_mob_ratio_move_plot         = fly_mob_ratio_move(plot_mask_move);
cat_labels_move                 = {group_defs(kept_idx_move).label};
cat_colors_move                 = all_cat_colors(kept_idx_move,:);

figure(3); clf
set(gcf,'Name','bump mobility (total-movement bout definition): path-length ratio per fly','Position',[100,100,max(700,120*numel(cat_labels_move)),600])
groupplot(cat_x_mob_move_plot,fly_mob_ratio_move_plot,cat_labels_move,cat_colors_move)
hold on
plot(xlim,[1,1],':k')
ylabel('bump path length / heading path length (per fly, bout-duration-weighted)')
title(sprintf('bump mobility, bouts defined by total movement (f\\_speed + %g*r\\_speed, thresh=%.2f), by genotype and light condition (flies with \\geq%d bouts of their own)', ...
    mobility_r_speed_weight,mobility_move_thresh,min_bouts_per_fly))

%% diagnostic figures with the total-movement bout definition, for every fly x light-condition with mobility ratio > 900
% fly_mob_ratio_move already only contains flies that individually
% cleared min_bouts_per_fly, so no further bout-count filtering is needed
% here.
diag_idx_move = find(fly_mob_ratio_move > 900);
fprintf('\n=== %d fly x light-condition instance(s), total-movement bouts, with mobility ratio > 1 ===\n', numel(diag_idx_move));

for d = 1:numel(diag_idx_move)
    idx  = diag_idx_move(d);
    gd   = cat_x_mob_move(idx);
    f    = fly_mob_fly_move(idx);
    trial_idx = find(fly_num==f & is_dark==group_defs(gd).dark)';

    [x_cat,alpha,fz_cat,mu_cat,heading_cat,is_walking_cat,bounds] = concat_trials_im_new(all_data,chunk_fz,chunk_mu_smooth_move,chunk_heading_smooth_move,chunk_is_walking_move,trial_idx);

    figure(400+d); clf
    set(gcf,'Name',sprintf('mobility diagnostic (total-movement bouts): %s, %s',fly_short_label(fly_list{f}),group_defs(gd).label),'Position',[100,100,1000,600])
    t = tiledlayout(3,1,'TileSpacing','compact');

    nexttile(t,[2,1]); hold on
    imagesc(x_cat,alpha,fz_cat)
    set(gca,'YDir','normal','XLim',[x_cat(1),x_cat(end)],'YLim',[alpha(1),alpha(end)])

    bout_edges = find(diff([0;is_walking_cat;0])~=0);
    for b = 1:2:numel(bout_edges)-1
        bout_x = x_cat([bout_edges(b),bout_edges(b+1)-1]);
        patch(bout_x([1,2,2,1]),alpha([1,1,end,end]),[1,1,0],'FaceAlpha',.15,'EdgeColor','none')
        xline(bout_x(1),'Color',[.9,.7,0],'LineWidth',1.2)
        xline(bout_x(2),'Color',[.9,.7,0],'LineWidth',1.2)
    end

    mu_plot = mu_cat;
    mu_plot(find(abs(diff(mu_plot))>pi)+1) = nan;
    plot(x_cat,mu_plot,'w','LineWidth',1)
    plot(x_cat,mu_plot+2*pi,'w','LineWidth',1)

    heading_plot = heading_cat;
    heading_plot(find(abs(diff(heading_plot))>pi)+1) = nan;
    plot(x_cat,heading_plot,'Color',[0,1,1],'LineWidth',1)
    plot(x_cat,heading_plot+2*pi,'Color',[0,1,1],'LineWidth',1)

    for b = bounds(1:end-1)
        xline(b,':','Color',[.6,.6,.6]);
    end

    cb = colorbar; ylabel(cb,'z-score')
    xlabel('frame'); ylabel('PB angle (rad)')
    title(sprintf('%s -- %s -- mov ratio=%.2f -- total-movement walking bouts (yellow), bump position (white), heading (-ft.cue, cyan)', ...
        fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio_move(idx)),'Interpreter','none')

    bout_mov_mu    = cat(1,chunk_bout_mu_move{trial_idx});
    bout_mov_speed = cat(1,chunk_bout_speed_move{trial_idx});
    bout_dur       = cat(1,chunk_bout_dur_move{trial_idx});

    nexttile(t); hold on
    scatter(bout_mov_speed,bout_mov_mu,20,bout_dur,'filled')
    xl = xlim; xl(1) = 0;
    plot(xl,fly_mob_ratio_move(idx)*xl,'-k')
    plot(xl,xl,':','Color',[.6,.6,.6])
    xlim(xl)
    cb2 = colorbar; ylabel(cb2,'bout duration (s)')
    xlabel('heading path length (rad)'); ylabel('bump path length (rad)')
    title(sprintf('%d bouts -- fit slope (mov ratio) = %.2f -- dotted line = slope 1',numel(bout_mov_mu),fly_mob_ratio_move(idx)))

    fprintf('  fly %-30s %-28s mov_ratio=%.2f (%d bouts)\n',fly_short_label(fly_list{f}),group_defs(gd).label,fly_mob_ratio_move(idx),numel(bout_mov_mu));
end

%% functions
function [fly_seg,trial_seg] = meta_fly_and_trial_seg(meta_path)
    % meta = ...\<project>\<date>\fly N\<trial folder>\registration_NNN
    % returns the "fly N" folder segment and the trial folder segment,
    % since genotype can be written into either one (see header comment).
    %
    % Split on a literal backslash, not filesep: all_data.meta is always a
    % Windows-style path string (this raw data was collected on Windows),
    % regardless of what platform this script itself runs on -- splitting
    % by filesep silently no-ops (leaving one giant "part") on any
    % non-Windows machine, e.g. a Mac, where filesep=='/'.
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        trial_seg = parts{end-1};
        fly_seg   = parts{end-2};
    else
        trial_seg = parts{end};
        fly_seg   = parts{end-1};
    end
end

function fid = meta_fly_id(meta_path)
    % unique fly identifier = project + date + "fly N" folder, i.e.
    % everything up through (and including) the "fly N" segment.
    % Split on a literal backslash, not filesep -- see comment in
    % meta_fly_and_trial_seg above.
    parts = strsplit(meta_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    if ~isempty(regexpi(parts{end},'^registration_\d+$','once')) && numel(parts) > 1
        fly_part_end = numel(parts)-2;
    else
        fly_part_end = numel(parts)-1;
    end
    fid = strjoin(parts(1:fly_part_end),'\');
end

function groupplot(cat_x, values, cat_labels, colors)
    % jittered per-point scatter + mean +/- SEM errorbar per category.
    % adapted from lpsp_kir_claude.m / lpsp_compartments_claude_script.m's
    % own groupplot: that version drew each category's "n=" as a separate
    % floating text() object anchored to a single shared y-coordinate
    % (the axis' own bottom edge) -- fine for the 2-4 categories those
    % scripts show, but with up to 10 genotype x light-condition
    % categories here, adjacent categories sit close enough together
    % (default xlim spacing of 1, +/-0.15 jitter) that a floating text
    % object anchored at a fixed y for every category can visually
    % collide with a neighboring category's own jittered points or its
    % own "n=" text, reading as if labels were duplicated/garbled. n is
    % appended directly into each category's OWN xtick label instead --
    % MATLAB lays out tick labels itself (spacing them), so they can't
    % overlap the plotted data or each other the way a manually positioned
    % text() can. Each label is kept on a SINGLE line (no embedded '\n')
    % -- xticklabels(), given a cell array of strings that themselves
    % contain a literal newline, splits each one on that newline into
    % SEPARATE list entries rather than rendering a 2-line label per tick
    % (confirmed directly: this is what produced the earlier "labels
    % interleaved between group name and n value" bug -- the flattened,
    % doubled-length list was assigned one entry per tick, so alternating
    % ticks got just the group name or just the "(n=X)" half).
    hold on
    n_cat = numel(cat_labels);
    labels_with_n = cell(1,n_cat);
    for c = 1:n_cat
        y = values(cat_x==c);
        y = y(~isnan(y));
        labels_with_n{c} = sprintf('%s (n=%d)',cat_labels{c},numel(y));
        if isempty(y)
            continue
        end
        jitter = (rand(size(y))-.5)*.3;
        scatter(c+jitter,y,20,colors(c,:),'filled','MarkerFaceAlpha',.3)
        errorbar(c,mean(y),std(y)/sqrt(numel(y)),'o','Color',colors(c,:)*.6, ...
            'MarkerFaceColor',colors(c,:)*.6,'LineWidth',2,'MarkerSize',7)
    end
    xticks(1:n_cat); xticklabels(labels_with_n); xtickangle(15)
    xlim([0.5,n_cat+0.5])
end

function [mov_mu,mov_speed,dur,mu_smooth,heading_smooth,is_walking] = trial_walking_bouts(trial,smooth_window,turn_thresh,max_gap_frames,min_walking_frames)
    % detects "walking bouts" (sustained turning) in one trial and returns
    % this trial's own per-bout bump path length (mov_mu) and heading path
    % length (mov_speed, from r_speed), plus each bout's duration -- same
    % method as lpsp_p2x2_walking_script.m's own bump-mobility section,
    % copied verbatim from lpsp_kir_claude.m. Also returns the full-trial
    % smoothed bump position (mu_smooth) and a smoothed heading POSITION
    % (heading_smooth, from -ft.cue) -- both on the fictrac timebase (xf),
    % UNWRAPPED -- purely so a diagnostic figure can display exactly (for
    % mu_smooth) or consistently (for heading_smooth, a position-domain
    % companion to the velocity-based r_speed_smooth actually used in
    % mov_speed) what this function's own smoothing did, instead of
    % plotting raw/unsmoothed traces. Both are unwrapped BEFORE smoothing
    % (a circular variable, gaussian-smoothed while still wrapped, would
    % smear across its false -pi/pi discontinuities) -- callers must wrap
    % back to -pi:pi themselves for display, after smoothing/interpolating.
    %
    % is_walking (full-trial logical, xf timebase) is also returned so a
    % diagnostic figure can shade which stretches of the trial actually
    % fed mov_mu/mov_speed -- mov_ratio is a regression over ONLY these
    % bout windows, not the whole trial.
    xf   = trial.ft.xf;
    dt   = median(diff(xf));
    n_im = numel(trial.im.mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu = interp1(xb,unwrap(trial.im.mu(:)),xf,'linear','extrap');

    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window); % force column, unlike cue, r_speed's native orientation isn't guaranteed
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    heading_smooth = smoothdata(unwrap(-trial.ft.cue(:)),'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    n_bouts   = numel(bout_starts);
    mov_mu    = nan(n_bouts,1);
    mov_speed = nan(n_bouts,1);
    dur       = nan(n_bouts,1);
    for b = 1:n_bouts
        rng          = bout_starts(b):bout_ends(b);
        mov_mu(b)    = sum(abs(diff(mu_smooth(rng))),'omitnan');
        mov_speed(b) = sum(abs(r_speed_smooth(rng)),'omitnan')*dt; % path length = integral of |speed| dt
        dur(b)       = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end

function [x_cat,alpha,z_cat,mu_cat,heading_cat,is_walking_cat,bounds] = concat_trials_im(all_data,chunk_mu_smooth,chunk_heading_smooth,chunk_is_walking,trial_idx)
    % concatenate the im.z/bump-position/heading traces of one or more
    % trials (for one fly, one light condition) end-to-end along the
    % imaging-frame axis, so a fly with more than one trial in a condition
    % still gets a single panel. concatenation is by frame count, not real
    % elapsed time, since trials aren't necessarily the same
    % duration/frame rate -- this is just a continuous strip to look at,
    % not an aligned timebase. Copied from lpsp_kir_claude.m's own
    % concat_trials_im, with three fixes:
    %  1) im.alpha's 32 entries repeat the SAME -pi:pi angular sequence
    %     twice (one full pass per PB hemisphere -- confirmed directly
    %     against this dataset), so unwrap() is required to turn it into
    %     the continuous -pi:~3pi range that matches how every legacy
    %     script in this repo plots im.z (always
    %     imagesc(...,unwrap(all_data(i).im.alpha),...), never raw alpha)
    %     -- using raw alpha here would make the 2nd hemisphere's wedges
    %     plot on top of the 1st's instead of stacked above them.
    %  2) mu/heading are no longer read raw off all_data -- they're pulled
    %     from chunk_mu_smooth/chunk_heading_smooth (trial_walking_bouts's
    %     own smoothed, unwrapped, fictrac-timebase outputs), interpolated
    %     onto this trial's own image timebase (xb) and only THEN wrapped
    %     to -pi:pi for display -- so the overlay shows the actual
    %     smoothed traces the mobility calculation used, not raw/unsmoothed
    %     ones (smoothing must happen before the -pi:pi wrap, never after,
    %     or it would smear across false discontinuities).
    %  3) also returns is_walking_cat (chunk_is_walking, nearest-neighbor
    %     interpolated onto xb) -- which stretches of the panel are inside
    %     a walking bout, i.e. actually fed mov_mu/mov_speed, since
    %     mov_ratio is a regression over ONLY those bouts, not the whole
    %     trial shown here.
    alpha = unwrap(all_data(trial_idx(1)).im.alpha(:));
    z_cat = []; mu_cat = []; heading_cat = []; is_walking_cat = []; x_cat = []; bounds = [];
    x_offset = 0;
    for k = 1:numel(trial_idx)
        i    = trial_idx(k);
        z    = all_data(i).im.z;
        n_im = size(z,2);
        x    = (0:n_im-1)' + x_offset;

        xf = all_data(i).ft.xf;
        xb = linspace(xf(1),xf(end),n_im)';

        mu         = interp1(xf,chunk_mu_smooth{i},xb,'linear','extrap');
        heading    = interp1(xf,chunk_heading_smooth{i},xb,'linear','extrap');
        mu         = wrap_to_pi(mu);
        heading    = wrap_to_pi(heading);
        is_walking = interp1(xf,double(chunk_is_walking{i}),xb,'nearest','extrap') > 0;

        z_cat          = [z_cat, z]; %#ok<AGROW>
        mu_cat         = [mu_cat; mu]; %#ok<AGROW>
        heading_cat    = [heading_cat; heading]; %#ok<AGROW>
        is_walking_cat = [is_walking_cat; is_walking]; %#ok<AGROW>
        x_cat          = [x_cat; x]; %#ok<AGROW>
        x_offset = x(end) + 1;
        bounds(end+1) = x(end); %#ok<AGROW>
    end
end

function x = wrap_to_pi(x)
    % wraps an unwrapped (possibly large-magnitude, cumulative) angle back
    % into -pi:pi for display -- only ever apply this AFTER any smoothing
    % or interpolation of a circular variable, never before (interpolating
    % or gaussian-smoothing a pre-wrapped angle smears across its false
    % -pi/pi jumps).
    x = mod(x,2*pi);
    x(x > pi) = x(x > pi) - 2*pi;
end

function lbl = fly_short_label(fly_path)
    % "<date> fly N" from a fly_id path like ...\<date folder>\fly N.
    % Split on a literal backslash, not filesep -- see comment in
    % meta_fly_and_trial_seg above.
    parts = strsplit(fly_path,{'\','/'});
    parts(cellfun(@isempty,parts)) = [];
    lbl = strjoin(parts(max(1,end-1):end),' ');
end

function [mu_new,rho_new,f_z] = trial_smoothed_bump_pva(trial,smooth_window_s)
    % re-estimates bump position (mu_new) and vector strength (rho_new)
    % directly from raw per-glomerulus fluorescence (im.f, 32 x n_im, same
    % layout as im.z), instead of trusting the upstream im.mu/im.rho fit:
    % each glomerulus's raw trace is gaussian-smoothed over
    % smooth_window_s SECONDS (converted to samples using this trial's own
    % imaging frame period, which isn't constant across trials in this
    % dataset), then z-scored over time (puts every glomerulus on the same
    % scale). mu_new/rho_new are then the angle/length of the mean
    % resultant vector of the 32 (smoothed+zscored) glomeruli, each
    % treated as a vector pointing at its own angular position (im.alpha)
    % -- the same population-vector-average formula this codebase's
    % legacy bump_calc_pb (e.g. epg_grab_script.m) uses on dff_cluster,
    % just applied to this new smoothed+zscored signal instead.
    n_im  = size(trial.im.f,2);
    xf    = trial.ft.xf;
    dt_im = (xf(end)-xf(1)) / (n_im-1);
    win_samples = max(1,round(smooth_window_s/dt_im));

    f_smooth = smoothdata(trial.im.f,2,'gaussian',win_samples);
    f_z      = (f_smooth - mean(f_smooth,2)) ./ std(f_smooth,0,2);

    alpha_row = trial.im.alpha(:)';
    [x_tmp,y_tmp]    = pol2cart(alpha_row,f_z');
    [mu_new,rho_new] = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
end

function [mov_mu,mov_speed,dur,mu_smooth,heading_smooth,is_walking] = trial_walking_bouts_from_mu(trial,mu,smooth_window_s,turn_thresh,max_gap_frames,min_walking_frames)
    % same walking-bout detection and per-bout path lengths as
    % trial_walking_bouts, but on a caller-supplied bump-position trace
    % (mu, on this trial's own im timebase) instead of trial.im.mu, AND
    % with the heading (ft.r_speed) trace smoothed over the SAME real-time
    % window as the bump-position re-estimate upstream (mu is typically
    % already heavily pre-smoothed by whatever produced it, e.g.
    % trial_smoothed_bump_pva) -- so mov_mu and mov_speed reflect matched
    % smoothing timescales, rather than a heavily-smoothed bump against a
    % lightly-smoothed heading. smooth_window_s is given in SECONDS, not
    % samples, and converted to samples using THIS trial's own fictrac
    % sample period (dt = median(diff(xf))) -- necessary because im.f
    % (the imaging timebase mu was originally derived from) and ft.r_speed
    % (fictrac timebase) are not sampled at the same rate, so a single
    % shared sample count would not represent the same real-time window
    % for both signals. mu itself is interpolated onto xf below before
    % smoothing, so at that point both traces share xf's timebase and a
    % single dt-derived sample count correctly applies to both.
    %
    % Also returns the full-trial smoothed bump position (mu_smooth) and a
    % smoothed heading POSITION (heading_smooth, from -ft.cue) -- both on
    % xf, UNWRAPPED -- purely so a diagnostic figure can display exactly
    % (for mu_smooth) or consistently (for heading_smooth, a
    % position-domain companion to the velocity-based r_speed_smooth
    % actually used in mov_speed) what this function's own smoothing did.
    % Both are unwrapped BEFORE smoothing -- a circular variable,
    % gaussian-smoothed while still wrapped, would smear across its false
    % -pi/pi discontinuities -- callers must wrap back to -pi:pi
    % themselves for display, after smoothing/interpolating.
    %
    % is_walking (full-trial logical, xf timebase) is also returned so a
    % diagnostic figure can shade which stretches of the trial actually
    % fed mov_mu/mov_speed -- mov_ratio is a regression over ONLY these
    % bout windows, not the whole trial, so a trace that looks dominated
    % by bump movement outside of them can visually mislead about what
    % mov_ratio itself reflects.
    xf   = trial.ft.xf;
    dt   = median(diff(xf));
    smooth_window = max(1,round(smooth_window_s/dt));
    n_im = numel(mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu_f = interp1(xb,unwrap(mu(:)),xf,'linear','extrap');

    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth      = smoothdata(mu_f,'gaussian',smooth_window);
    heading_smooth = smoothdata(unwrap(-trial.ft.cue(:)),'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    n_bouts   = numel(bout_starts);
    mov_mu    = nan(n_bouts,1);
    mov_speed = nan(n_bouts,1);
    dur       = nan(n_bouts,1);
    for b = 1:n_bouts
        rng          = bout_starts(b):bout_ends(b);
        mov_mu(b)    = sum(abs(diff(mu_smooth(rng))),'omitnan');
        mov_speed(b) = sum(abs(r_speed_smooth(rng)),'omitnan')*dt;
        dur(b)       = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end

function [x_cat,alpha,fz_cat,mu_cat,heading_cat,is_walking_cat,bounds] = concat_trials_im_new(all_data,chunk_fz,chunk_mu_smooth,chunk_heading_smooth,chunk_is_walking,trial_idx)
    % same concatenation as concat_trials_im (one fly, one light
    % condition, joined end-to-end along the imaging-frame axis by frame
    % count), but pulling the smoothed+zscored fluorescence from the
    % precomputed chunk_fz cell array (trial_smoothed_bump_pva's output)
    % instead of all_data(i).im.z. alpha is unwrapped for the same reason
    % as concat_trials_im above: im.alpha repeats the same -pi:pi sequence
    % once per PB hemisphere, so unwrap() is needed to get the continuous
    % -pi:~3pi range this codebase's imagesc panels expect.
    %
    % mu/heading come from chunk_mu_smooth/chunk_heading_smooth --
    % trial_walking_bouts_from_mu's own smoothed, unwrapped, fictrac-
    % timebase outputs (NOT chunk_mu_new, which is only the
    % fluorescence-level-smoothed PVA estimate, before the additional
    % bout-level smoothing mov_mu is actually computed from) --
    % interpolated onto this trial's image timebase (xb) and only THEN
    % wrapped to -pi:pi for display, so the overlay matches exactly what
    % the mobility calculation used. is_walking_cat (chunk_is_walking,
    % nearest-neighbor interpolated onto xb) marks which stretches of the
    % panel actually fed mov_mu/mov_speed -- mov_ratio is a regression over
    % ONLY those bouts, not the whole trial shown here.
    alpha = unwrap(all_data(trial_idx(1)).im.alpha(:));
    fz_cat = []; mu_cat = []; heading_cat = []; is_walking_cat = []; x_cat = []; bounds = [];
    x_offset = 0;
    for k = 1:numel(trial_idx)
        i    = trial_idx(k);
        fz   = chunk_fz{i};
        n_im = size(fz,2);
        x    = (0:n_im-1)' + x_offset;

        xf = all_data(i).ft.xf;
        xb = linspace(xf(1),xf(end),n_im)';

        mu         = interp1(xf,chunk_mu_smooth{i},xb,'linear','extrap');
        heading    = interp1(xf,chunk_heading_smooth{i},xb,'linear','extrap');
        mu         = wrap_to_pi(mu);
        heading    = wrap_to_pi(heading);
        is_walking = interp1(xf,double(chunk_is_walking{i}),xb,'nearest','extrap') > 0;

        fz_cat         = [fz_cat, fz]; %#ok<AGROW>
        mu_cat         = [mu_cat; mu]; %#ok<AGROW>
        heading_cat    = [heading_cat; heading]; %#ok<AGROW>
        is_walking_cat = [is_walking_cat; is_walking]; %#ok<AGROW>
        x_cat          = [x_cat; x]; %#ok<AGROW>
        x_offset = x(end) + 1;
        bounds(end+1) = x(end); %#ok<AGROW>
    end
end

function total_move = trial_total_movement(trial,smooth_window,r_speed_weight)
    % combines forward speed (ft.f_speed, signed +/- for forward/backward)
    % and rotational speed (ft.r_speed, rad/s, signed +/- for left/right)
    % into one non-negative "how much is the fly moving" magnitude, on the
    % fictrac timebase (xf). Each channel is smoothed WHILE STILL SIGNED
    % (so e.g. a fly oscillating rapidly forward/backward or left/right
    % averages toward ~0 net movement, not toward a large "moving" value)
    % and only then rectified -- exactly the same order of operations this
    % script already uses for r_speed_smooth/fly_speed in
    % trial_walking_bouts. r_speed_weight converts r_speed onto roughly
    % f_speed's own scale before combining -- see the header comment on
    % this script's "walking bouts redefined by TOTAL movement" section
    % for why that's an inherently approximate, tunable choice, not a
    % first-principles conversion. The two are combined by Euclidean
    % magnitude (not a plain sum) -- see that same header comment.
    f_speed_smooth = smoothdata(trial.ft.f_speed(:),'gaussian',smooth_window);
    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    total_move = sqrt(f_speed_smooth.^2 + (r_speed_weight*r_speed_smooth).^2);
end

function [mov_mu,mov_speed,dur,mu_smooth,heading_smooth,is_walking] = trial_walking_bouts_totalmovement(trial,mu,smooth_window_s,move_thresh,r_speed_weight,max_gap_frames,min_walking_frames)
    % same as trial_walking_bouts_from_mu (bump path length from a
    % caller-supplied mu, heading path length from r_speed, both smoothed
    % over smooth_window_s converted to samples via this trial's own
    % fictrac sample period), EXCEPT bout membership (is_walking) is
    % decided by trial_total_movement (forward + rotational speed
    % combined) instead of |r_speed_smooth| alone -- so a fly walking
    % forward while barely turning still counts as "walking" and doesn't
    % get fragmented into separate bouts. mov_speed/mov_mu themselves are
    % UNCHANGED (still pure heading-angle / bump-angle path length) --
    % only which time windows qualify as a bout is different here.
    xf   = trial.ft.xf;
    dt   = median(diff(xf));
    smooth_window = max(1,round(smooth_window_s/dt));
    n_im = numel(mu);
    xb   = linspace(xf(1),xf(end),n_im)';

    mu_f = interp1(xb,unwrap(mu(:)),xf,'linear','extrap');

    r_speed_smooth = smoothdata(trial.ft.r_speed(:),'gaussian',smooth_window);
    mu_smooth      = smoothdata(mu_f,'gaussian',smooth_window);
    heading_smooth = smoothdata(unwrap(-trial.ft.cue(:)),'gaussian',smooth_window);
    total_move     = trial_total_movement(trial,smooth_window,r_speed_weight);

    is_walking = total_move > move_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    n_bouts   = numel(bout_starts);
    mov_mu    = nan(n_bouts,1);
    mov_speed = nan(n_bouts,1);
    dur       = nan(n_bouts,1);
    for b = 1:n_bouts
        rng          = bout_starts(b):bout_ends(b);
        mov_mu(b)    = sum(abs(diff(mu_smooth(rng))),'omitnan');
        mov_speed(b) = sum(abs(r_speed_smooth(rng)),'omitnan')*dt;
        dur(b)       = (bout_ends(b)-bout_starts(b)+1)*dt;
    end
end
