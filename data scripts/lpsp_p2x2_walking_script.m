%% load data
%clear all
close all

load('.data/lpsp_p2x2_walking_20260728.mat') %loads all_data

%% subset to trials collected on or after 7/8
trial_date = cellfun(@(x)(str2double(regexp(x,'\d{8}','match','once'))),{all_data.meta});
all_data = all_data(trial_date >= 20260708);

is_empty = arrayfun(@(x)(contains(x.meta,'_empty_')),all_data);
is_lpsp  = arrayfun(@(x)(contains(x.meta,'_lpsp_')),all_data);
group_idx = is_lpsp + 1; %1 = empty, 2 = lpsp
group_labels = {'empty','lpsp'};

n = length(all_data);
fprintf('%i trials on or after 7/8 (%i empty, %i lpsp)\n',n,sum(is_empty),sum(is_lpsp))

%% 1) how accurately does mu track the fly's actual heading?
rho_thresh = .1; %minimum bump vector strength (im.rho) to trust a mu estimate

mu_err  = nan(n,1); %mean absolute circular error between mu and heading, in radians
mu_corr = nan(n,1); %circular-circular correlation between mu and heading

for i = 1:n
    xb  = all_data(i).ft.xb;
    xf  = all_data(i).ft.xf;
    cue = -all_data(i).ft.cue; %cue is stored with the opposite sign convention of mu (see lpsp_p2x2_script_minimal.m)

    mu  = interp1(xb,unwrap(all_data(i).im.mu),xf,'linear','extrap'); %upsample mu from the imaging clock (xb) onto the faster behavior clock (xf)
    rho = interp1(xb,all_data(i).im.rho,xf,'linear','extrap');

    idx = rho > rho_thresh;

    mu_err(i)  = mean(abs(circ_dist(cue(idx),mu(idx))),'omitnan');
    mu_corr(i) = circ_corrcc(mod(cue(idx),2*pi),mod(mu(idx),2*pi));
end

figure(1); clf

subplot(1,2,1); hold on
swarmchart(group_idx,mu_err,'filled','MarkerFaceAlpha',.5)
for g = 1:2
    errorbar(g,mean(mu_err(group_idx==g),'omitnan'),std(mu_err(group_idx==g),'omitnan')/sqrt(sum(group_idx==g)),'ok')
end
xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
ylabel('mean |heading error| (rad)')
title('bump accuracy')

subplot(1,2,2); hold on
swarmchart(group_idx,mu_corr,'filled','MarkerFaceAlpha',.5)
for g = 1:2
    errorbar(g,mean(mu_corr(group_idx==g),'omitnan'),std(mu_corr(group_idx==g),'omitnan')/sqrt(sum(group_idx==g)),'ok')
end
xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
ylabel('circular correlation (mu vs heading)')
title('bump-heading correlation')

sgtitle('mu vs fly heading, trials from 7/8 onward')

%% 2) how much does the bump (mu) move relative to how much the fly turns?
%path-length gain metric, adapted from gain_metric_sandbox.m
smooth_window      = 60;  %samples (~1s at xf's ~60Hz), gaussian smoothing window
turn_thresh        = .25; %rad/s, minimum heading speed to call a bout "walking"
max_gap_frames     = .5*60; %frames of non-walking allowed within a bout before splitting it
min_walking_frames = .5*60; %minimum bout length to keep

mov_ratio = nan(n,1); %ratio of bump path length to heading path length, per trial

for i = 1:n
    xf = all_data(i).ft.xf;
    dt = median(diff(xf));

    mu  = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');
    cue = unwrap(all_data(i).ft.cue);

    mu_smooth  = smoothdata(mu,'gaussian',smooth_window);
    cue_smooth = smoothdata(-cue,'gaussian',smooth_window);
    fly_speed  = [abs(diff(cue_smooth));nan]/dt;

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    mov_mu  = nan(length(bout_starts),1);
    mov_cue = nan(length(bout_starts),1);
    for b = 1:length(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(diff(cue_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
    end

    mov_ratio(i) = mov_cue \ mov_mu; %bump path length per unit heading path length, across all walking bouts in the trial
end

figure(2); clf; hold on
swarmchart(group_idx,mov_ratio,'filled','MarkerFaceAlpha',.5)
for g = 1:2
    errorbar(g,mean(mov_ratio(group_idx==g),'omitnan'),std(mov_ratio(group_idx==g),'omitnan')/sqrt(sum(group_idx==g)),'ok')
end
plot(xlim,[1,1],':k')
xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
ylabel('bump path length / heading path length')
title('bump mobility relative to fly turning, trials from 7/8 onward')

%% 3) label each trial: which fly it belongs to, its order within that fly, lighting condition, and whether it's the perturbation trial
is_dark = arrayfun(@(x)(contains(x.ft.pattern,'background')),all_data); %closed-loop cue vs dark (background pattern)

fly_id       = cell(n,1);
trial_num    = nan(n,1);
atp_activity = nan(n,1); %peak red-channel (atp) signal in this trial - not yet thresholded

for i = 1:n
    parts   = strsplit(all_data(i).meta,'\');
    fly_pos = find(startsWith(parts,'fly '),1);
    fly_id{i} = strjoin(parts(1:fly_pos),'\');
    trial_num(i) = str2double(regexp(parts{fly_pos+1},'-(\d+)_','tokens','once'));

    if isfield(all_data(i),'atp') && isfield(all_data(i).atp,'d') && ~isempty(all_data(i).atp.d)
        atp_activity(i) = max(smoothdata(sum(all_data(i).atp.d,1),'movmean',20));
    end
end

[fly_list,~,fly_ix] = unique(fly_id);
n_flies = length(fly_list);

%% flag pulse trials per fly, relative to that fly's own atp baseline (an absolute threshold doesn't generalize across flies with different expression/baseline levels)
pulse_outlier_factor = 1.5; %trials whose atp_activity is this many scaled-MADs above their fly's own median are flagged as pulse trials
min_trials            = 3; %need at least this many trials (either lighting) for a fly to judge outliers robustly

has_pulse = false(n,1); %true if this trial is a detected perturbation trial

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    valid_idx = trial_idx(~isnan(atp_activity(trial_idx))'); %baseline is built from both closed-loop and dark trials, so a fly's normal dark-trial fluorescence doesn't get mistaken for a pulse

    if length(valid_idx) < min_trials; continue; end

    is_out  = isoutlier(atp_activity(valid_idx),'median','ThresholdFactor',pulse_outlier_factor)';
    cl_mask = ~is_dark(valid_idx); %but only closed-loop trials are ever actually flagged as pulses

    has_pulse(valid_idx(is_out & cl_mask)) = true;
end

%% sanity-check the detector against the expected layout (most flies: 9 trials, pulses at trials 6 and 8)
n_pulse_per_fly   = nan(n_flies,1);
pulse_pos_per_fly = cell(n_flies,1);
n_trials_per_fly  = nan(n_flies,1);
fly_date          = cell(n_flies,1);
fly_num           = nan(n_flies,1);

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    [~,order] = sort(trial_num(trial_idx));
    trial_idx = trial_idx(order);

    n_trials_per_fly(f)  = length(trial_idx);
    pulse_pos_per_fly{f} = find(has_pulse(trial_idx))';
    n_pulse_per_fly(f)   = length(pulse_pos_per_fly{f});

    fly_date{f} = regexp(fly_list{f},'\d{8}','match','once');
    fly_num(f)  = str2double(regexp(fly_list{f},'fly (\d+)','tokens','once'));
end

fprintf('\n%-10s %-6s %-8s %-8s %s\n','date','fly','n_trials','n_pulses','pulse trial(s)')
for f = 1:n_flies
    fprintf('%-10s %-6i %-8i %-8i %s\n',fly_date{f},fly_num(f),n_trials_per_fly(f),n_pulse_per_fly(f),mat2str(pulse_pos_per_fly{f}))
end

fprintf('\n%i of %i flies had exactly 2 detected pulse trials\n',sum(n_pulse_per_fly==2),n_flies)
is_typical = n_trials_per_fly==8;
matches_expected = cellfun(@(p)(isequal(p,[5,7])),pulse_pos_per_fly);
fprintf('of the %i flies with 9 trials, %i had pulses detected at exactly trials 6 and 8\n',sum(is_typical),sum(is_typical & matches_expected))

%% 4) for each fly, average bump mobility before vs after its perturbation trial, split by lighting condition
lighting_labels = {'closed loop','dark'};

fly_group = nan(n_flies,1);            %1 = empty, 2 = lpsp
fly_pre   = nan(n_flies,2);            %columns: [closed loop, dark]
fly_post  = nan(n_flies,2);

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    [~,order] = sort(trial_num(trial_idx));
    trial_idx = trial_idx(order); %chronological order of this fly's trials

    fly_group(f) = group_idx(trial_idx(1));

    pulse_pos = find(has_pulse(trial_idx));
    if isempty(pulse_pos); continue; end %no detected perturbation for this fly - can't split pre/post

    pre_idx  = trial_idx(1:pulse_pos(1)-1);
    post_idx = trial_idx(pulse_pos(end)+1:end);

    for li = 1:2 %1 = closed loop, 2 = dark
        pre_this  = pre_idx(is_dark(pre_idx)   == (li==2));
        post_this = post_idx(is_dark(post_idx) == (li==2));
        if isempty(pre_this) || isempty(post_this); continue; end

        fly_pre(f,li)  = mean(mov_ratio(pre_this),'omitnan');
        fly_post(f,li) = mean(mov_ratio(post_this),'omitnan');
    end
end

fprintf('%i of %i flies had a detectable perturbation trial\n',sum(any(~isnan(fly_pre) | ~isnan(fly_post),2)),n_flies)

%% 5) plot pre vs post bump mobility, split by genotype and lighting condition
figure(3); clf
for g = 1:2
    for li = 1:2
        subplot(2,2,(li-1)*2+g); hold on

        valid = fly_group==g & ~isnan(fly_pre(:,li)) & ~isnan(fly_post(:,li));
        pre_vals  = fly_pre(valid,li);
        post_vals = fly_post(valid,li);

        plot([ones(sum(valid),1),2*ones(sum(valid),1)]',[pre_vals,post_vals]','Color',[.5,.5,.5,.3])
        scatter(ones(sum(valid),1), pre_vals,'filled','MarkerFaceAlpha',.5)
        scatter(2*ones(sum(valid),1),post_vals,'filled','MarkerFaceAlpha',.5)
        errorbar([1,2],[mean(pre_vals,'omitnan'),mean(post_vals,'omitnan')],...
                       [std(pre_vals,'omitnan'),std(post_vals,'omitnan')]/sqrt(sum(valid)),'-ok','LineWidth',1.5)

        xlim([.5,2.5]); xticks([1,2]); xticklabels({'pre','post'})
        ylabel('bump path length / heading path length')
        title(sprintf('%s, %s (n=%i flies)',group_labels{g},lighting_labels{li},sum(valid)))
    end
end
sgtitle('bump mobility before vs after perturbation')

%% 6) summarize the pre -> post change across the four groups
delta = fly_post - fly_pre; %columns: [closed loop, dark]

figure(4); clf
for li = 1:2
    subplot(1,2,li); hold on
    swarmchart(fly_group,delta(:,li),'filled','MarkerFaceAlpha',.5)
    for g = 1:2
        valid = fly_group==g & ~isnan(delta(:,li));
        errorbar(g,mean(delta(valid,li),'omitnan'),std(delta(valid,li),'omitnan')/sqrt(sum(valid)),'ok')
    end
    plot(xlim,[0,0],':k')
    xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
    ylabel('\Delta bump mobility (post - pre)')
    title(lighting_labels{li})
end
sgtitle('change in bump mobility after perturbation')
