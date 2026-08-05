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

mu_err  = nan(n,1); %circular variance of the mu-heading offset - i.e. how much that offset drifts, not its (arbitrary) absolute size
mu_corr = nan(n,1); %circular-circular correlation between mu and heading

for i = 1:n
    xb  = all_data(i).ft.xb;
    xf  = all_data(i).ft.xf;
    cue = -all_data(i).ft.cue; %cue is stored with the opposite sign convention of mu (see lpsp_p2x2_script_minimal.m)

    mu  = interp1(xb,unwrap(all_data(i).im.mu),xf,'linear','extrap'); %upsample mu from the imaging clock (xb) onto the faster behavior clock (xf)
    rho = interp1(xb,all_data(i).im.rho,xf,'linear','extrap');

    idx = rho > rho_thresh;

    mu_err(i)  = circ_var(circ_dist(cue(idx),mu(idx))); %the offset between mu and heading is arbitrary (depends on how mu's zero-point happens to be defined) - what matters is how much it drifts, i.e. its circular variance, not its mean magnitude
    mu_corr(i) = circ_corrcc(mod(cue(idx),2*pi),mod(mu(idx),2*pi));
end

figure(1); clf

subplot(1,2,1); hold on
h_err = swarmchart(group_idx,mu_err,'filled','MarkerFaceAlpha',.5,'XJitterWidth',.3);
for g = 1:2
    errorbar(g,mean(mu_err(group_idx==g),'omitnan'),std(mu_err(group_idx==g),'omitnan')/sqrt(sum(group_idx==g)),'ok')
end
xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
ylabel('circular variance of mu-heading offset')
title('bump accuracy (offset drift)')

subplot(1,2,2); hold on
h_corr = swarmchart(group_idx,mu_corr,'filled','MarkerFaceAlpha',.5,'XJitterWidth',.3);
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

mov_ratio     = nan(n,1); %ratio of bump path length to heading path length, per trial
frac_walking  = nan(n,1); %fraction of this trial the fly spent walking
mean_rho      = nan(n,1); %mean bump vector strength (confidence) over this trial

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

    mov_ratio(i)    = mov_cue \ mov_mu; %bump path length per unit heading path length, across all walking bouts in the trial
    frac_walking(i) = mean(is_walking,'omitnan');
    mean_rho(i)     = mean(all_data(i).im.rho,'omitnan');
end

figure(2); clf; hold on
h_mov = swarmchart(group_idx,mov_ratio,'filled','MarkerFaceAlpha',.5,'XJitterWidth',.3);
for g = 1:2
    errorbar(g,mean(mov_ratio(group_idx==g),'omitnan'),std(mov_ratio(group_idx==g),'omitnan')/sqrt(sum(group_idx==g)),'ok')
end
plot(xlim,[1,1],':k')
xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
ylabel('bump path length / heading path length')
title('bump mobility relative to fly turning, trials from 7/8 onward')

%% 3) label each trial: which fly it belongs to, its order within that fly, lighting condition, and whether it's the perturbation trial
is_dark = arrayfun(@(x)(contains(x.ft.pattern,'background')),all_data); %closed-loop cue vs dark (background pattern)

%a trial is labeled a perturbation trial if its stim-triggered average atp trace shows a defined peak just after
%the stim's rising edge (ft.stims, the ground-truth ejection marker) - rather than comparing raw activity across trials
peri_win    = 5;      %seconds before/after each stim to check for a response
t_common    = linspace(-peri_win,peri_win,101); %common relative-time axis stims get interpolated onto
peak_win    = [0,3];  %seconds after stim onset to look for the atp peak
base_win    = [-peri_win,0]; %seconds before stim onset used as this trial's own baseline
peak_factor = 3;      %the post-stim peak must clear this many baseline-SDs above the pre-stim baseline to count as "a defined peak"

fly_id    = cell(n,1);
trial_num = nan(n,1);
has_pulse = false(n,1); %true if this trial shows a stim-locked atp peak
atp_peri  = cell(n,1);  %this trial's stim-averaged atp trace (on t_common) - kept around for later plotting/QC

for i = 1:n
    parts   = strsplit(all_data(i).meta,'\');
    fly_pos = find(startsWith(parts,'fly '),1);
    fly_id{i} = strjoin(parts(1:fly_pos),'\');
    trial_num(i) = str2double(regexp(parts{fly_pos+1},'-(\d+)_','tokens','once'));

    if ~(isfield(all_data(i),'atp') && isfield(all_data(i).atp,'d') && ~isempty(all_data(i).atp.d)); continue; end

    stims  = logical(all_data(i).ft.stims(:));
    onsets = find(diff([false;stims])==1); %rising edge of each stim
    if isempty(onsets); continue; end %no ejections delivered this trial - can't be a perturbation trial

    xb      = all_data(i).ft.xb;
    xf      = all_data(i).ft.xf;
    t_onset = xf(onsets);
    atp_sig = sum(all_data(i).atp.d,1);

    atp_aligned = nan(length(onsets),length(t_common));
    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s); %imaging clock, re-centered on this stim's onset
        atp_aligned(s,:) = interp1(xb_rel,atp_sig,t_common,'linear',nan);
    end

    atp_peri{i} = mean(atp_aligned,1,'omitnan');

    base_idx = t_common >= base_win(1) & t_common <  base_win(2);
    peak_idx = t_common >  peak_win(1) & t_common <= peak_win(2);

    base_mean = median(atp_peri{i}(base_idx),'omitnan'); %median, not mean/max - a single noisy sample shouldn't drive the call
    base_std  = std(atp_peri{i}(base_idx),'omitnan');
    peak_amp  = median(atp_peri{i}(peak_idx),'omitnan');

    has_pulse(i) = (peak_amp - base_mean) > peak_factor*base_std;
end

[fly_list,~,fly_ix] = unique(fly_id);
n_flies = length(fly_list);

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

%% sanity check: ejections should only ever happen during closed-loop trials
dark_pulse_idx = find(has_pulse & is_dark(:));
if isempty(dark_pulse_idx)
    fprintf('sanity check passed: all %i detected perturbation trials are closed-loop trials\n',sum(has_pulse))
else
    fprintf('WARNING: %i detected perturbation trial(s) are DARK trials - check these:\n',length(dark_pulse_idx))
    for i = dark_pulse_idx'
        f = fly_ix(i);
        fprintf('  %s fly %i, trial %i\n',fly_date{f},fly_num(f),trial_num(i))
    end
end

%% debug: dump ground-truth stims vs. the detector's decision, trial by trial, for one fly
%set debug_date/debug_num to whichever fly looks suspicious and re-run this cell
debug_date = '20260728';
debug_num  = 12;

f = find(strcmp(fly_date,debug_date) & fly_num==debug_num);
if isempty(f)
    fprintf('\nno fly matches %s fly %i\n',debug_date,debug_num)
else
    trial_idx = find(fly_ix == f);
    [~,order] = sort(trial_num(trial_idx));
    trial_idx = trial_idx(order); %pos 1..end, chronological, matching the "pulse trial(s)" column above

    fprintf('\ndebugging %s fly %i (%i trials)\n',debug_date,debug_num,length(trial_idx))
    fprintf('%-4s %-7s %-8s %-10s %-9s %-9s %-9s %-9s %-7s %s\n',...
        'pos','trial#','n_stims','has_pulse','peak_amp','base_mn','base_sd','thresh','z','trial folder (raw metadata)')

    debug_rows = ceil(sqrt(length(trial_idx)));
    debug_cols = ceil(length(trial_idx)/debug_rows);
    figure(11); clf

    for pos = 1:length(trial_idx)
        i = trial_idx(pos);

        stims  = logical(all_data(i).ft.stims(:)); %ground truth: does this trial's raw DAQ record any ejection at all?
        onsets = find(diff([false;stims])==1);

        subplot(debug_rows,debug_cols,pos); hold on

        if ~isempty(atp_peri{i})
            base_idx  = t_common >= base_win(1) & t_common <  base_win(2);
            peak_idx  = t_common >  peak_win(1) & t_common <= peak_win(2);
            base_mean = median(atp_peri{i}(base_idx),'omitnan');
            base_std  = std(atp_peri{i}(base_idx),'omitnan');
            peak_amp  = median(atp_peri{i}(peak_idx),'omitnan');
            z         = (peak_amp-base_mean)/base_std;
            thresh    = base_mean + peak_factor*base_std;

            yl = [min(atp_peri{i}),max([atp_peri{i},thresh])];
            yl = yl + [-1,1]*.1*max(diff(yl),eps);

            patch(base_win([1,2,2,1]),yl([1,1,2,2]),'k','FaceAlpha',.08,'EdgeColor','none') %pre-stim baseline window
            patch(peak_win([1,2,2,1]),yl([1,1,2,2]),'g','FaceAlpha',.08,'EdgeColor','none') %post-stim peak window

            plot(t_common,atp_peri{i},'r','LineWidth',1.5)
            plot([-peri_win,peri_win],[base_mean,base_mean],':k')
            plot([-peri_win,peri_win],[thresh,thresh],':r')
            plot([0,0],yl,':k')

            xlim([-peri_win,peri_win]); ylim(yl)
        else
            base_mean = nan; base_std = nan; peak_amp = nan; z = nan; thresh = nan;
            text(0,0,'no stims this trial','HorizontalAlignment','center')
            axis off
        end

        parts        = strsplit(all_data(i).meta,'\');
        trial_folder = parts{find(startsWith(parts,'fly '),1)+1};

        title(sprintf('pos %i, trial %i, has\\_pulse=%i',pos,trial_num(i),has_pulse(i)))

        fprintf('%-4i %-7i %-8i %-10i %-9.3f %-9.3f %-9.3f %-9.3f %-7.2f %s\n',...
            pos,trial_num(i),length(onsets),has_pulse(i),peak_amp,base_mean,base_std,thresh,z,trial_folder)
    end
    sgtitle(sprintf('%s fly %i: peri-stim atp average per trial (baseline=gray, peak window=green, threshold=dashed red)',debug_date,debug_num))
end

%% 4) for each fly, average bump mobility before vs after its perturbation trial, split by lighting condition
lighting_labels = {'closed loop','dark'};

fly_group = nan(n_flies,1);            %1 = empty, 2 = lpsp
fly_pre   = nan(n_flies,2);            %columns: [closed loop, dark]
fly_post  = nan(n_flies,2);

fly_pre_idx  = cell(n_flies,2);        %the actual trial indices (into all_data) behind fly_pre/fly_post,
fly_post_idx = cell(n_flies,2);        %kept so we can go back and plot specific example trials later

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    [~,order] = sort(trial_num(trial_idx));
    trial_idx = trial_idx(order); %chronological order of this fly's trials

    fly_group(f) = group_idx(trial_idx(1));

    pulse_pos = find(has_pulse(trial_idx));
    if isempty(pulse_pos); continue; end %no detected perturbation for this fly - can't split pre/post

    pre_idx  = trial_idx(1:pulse_pos(1)-1);
    post_idx = trial_idx(pulse_pos(1)+1:end);
    post_idx = post_idx(~has_pulse(post_idx)); %exclude any other detected perturbation trials (e.g. a second stim trial), but keep everything else - including the closed-loop trial between two stim trials

    for li = 1:2 %1 = closed loop, 2 = dark
        pre_this  = pre_idx(is_dark(pre_idx)   == (li==2));
        post_this = post_idx(is_dark(post_idx) == (li==2));
        if isempty(pre_this) || isempty(post_this); continue; end

        fly_pre(f,li)  = mean(mov_ratio(pre_this),'omitnan');
        fly_post(f,li) = mean(mov_ratio(post_this),'omitnan');

        fly_pre_idx{f,li}  = pre_this;
        fly_post_idx{f,li} = post_this;
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
h_delta = gobjects(1,2);
for li = 1:2
    subplot(1,2,li); hold on
    h_delta(li) = scatter(fly_group,delta(:,li),'filled','MarkerFaceAlpha',.5);
    for g = 1:2
        valid = fly_group==g & ~isnan(delta(:,li));
        errorbar(g+.1,mean(delta(valid,li),'omitnan'),std(delta(valid,li),'omitnan')/sqrt(sum(valid)),'ok')
    end
    plot(xlim,[0,0],':k')
    xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
    ylabel('\Delta bump mobility (post - pre)')
    title(lighting_labels{li})
end
sgtitle('change in bump mobility after perturbation')

%% 7) bootstrap test: is the pre -> post change significantly different from zero in each group?
n_boot = 10000;

boot_p  = nan(2,2);   %rows: genotype, columns: lighting
boot_ci = nan(2,2,2); %rows: genotype, columns: lighting, page: [lo,hi]

fprintf('\nbootstrap test of mean(post - pre) ~= 0:\n')
for g = 1:2
    for li = 1:2
        vals = delta(fly_group==g & ~isnan(delta(:,li)),li);
        if isempty(vals); continue; end

        boot_means = nan(n_boot,1);
        for b = 1:n_boot
            boot_means(b) = mean(vals(randi(length(vals),length(vals),1)));
        end

        boot_p(g,li)    = 2*min(mean(boot_means>=0),mean(boot_means<=0)); %two-sided bootstrap p-value
        boot_ci(g,li,:) = prctile(boot_means,[2.5,97.5]);

        fprintf('%-6s %-12s mean = %+.3f, 95%% CI = [%+.3f, %+.3f], p = %.4f (n=%i flies)\n',...
            group_labels{g},lighting_labels{li},mean(vals),boot_ci(g,li,1),boot_ci(g,li,2),boot_p(g,li),length(vals))
    end
end

figure(4)
for li = 1:2
    subplot(1,2,li)
    for g = 1:2
        if isnan(boot_p(g,li)); continue; end
        text(g,max(ylim)*.9,sprintf('p = %.3f',boot_p(g,li)),'HorizontalAlignment','center')
    end
end

%% 8) pick one representative example fly per genotype: has a pre/post delta in both lighting conditions,
%walks a lot, and has confident (high rho) bump estimates
fly_walk_frac = nan(n_flies,1); %fraction of time this fly spent walking, averaged across its trials
fly_rho       = nan(n_flies,1); %mean bump vector strength across this fly's trials

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    fly_walk_frac(f) = mean(frac_walking(trial_idx),'omitnan');
    fly_rho(f)       = mean(mean_rho(trial_idx),'omitnan');
end

example_fly = nan(2,1); %index into fly_list/fly_ix, one per genotype

for g = 1:2
    candidates = find(fly_group==g & all(~isnan(delta),2)); %needs a valid delta in both lighting conditions

    is_good_walker = fly_walk_frac(candidates) >= median(fly_walk_frac(candidates),'omitnan');
    is_clean_bump  = fly_rho(candidates)       >= median(fly_rho(candidates),'omitnan');
    good_candidates = candidates(is_good_walker & is_clean_bump);
    if isempty(good_candidates); good_candidates = candidates; end %fall back if nothing clears both bars
    if isempty(good_candidates); continue; end

    mean_delta = mean(delta(fly_group==g,:),1,'omitnan'); %group-average [closed loop, dark] delta
    dist = sqrt(sum((delta(good_candidates,:) - mean_delta).^2,2));
    [~,best] = min(dist);
    example_fly(g) = good_candidates(best);

    fprintf('%s example fly: %s fly %i (walk frac = %.2f, rho = %.2f, delta = [%+.3f, %+.3f])\n',...
        group_labels{g},fly_date{example_fly(g)},fly_num(example_fly(g)),...
        fly_walk_frac(example_fly(g)),fly_rho(example_fly(g)),delta(example_fly(g),1),delta(example_fly(g),2))
end

%% 9) plot example pre/post traces for the representative fly in each group
for g = 1:2
    f = example_fly(g);
    if isnan(f); continue; end

    figure(4+g); clf
    for li = 1:2
        for phase = 1:2
            if phase == 1
                cand = fly_pre_idx{f,li};
                phase_label = 'pre';
            else
                cand = fly_post_idx{f,li};
                phase_label = 'post';
            end
            if isempty(cand); continue; end

            [~,best] = max(mean_rho(cand)); %show the cleanest trial among the candidates for this cell
            i = cand(best);

            subplot(2,2,(li-1)*2+phase); hold on
            imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
            a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'c','linewidth',1.5); a.YData(abs(diff(a.YData))>pi) = nan;
            a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w','linewidth',1.5); a.YData(abs(diff(a.YData))>pi) = nan;
            axis tight
            yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})
            xlabel('time (s)')
            title(sprintf('%s, %s (trial %i, mov ratio = %.2f)',lighting_labels{li},phase_label,trial_num(i),mov_ratio(i)))
        end
    end
    sgtitle(sprintf('%s example fly: %s fly %i',group_labels{g},fly_date{f},fly_num(f)))
end

%% 10) highlight the two example flies (one per genotype) on every earlier trial- and fly-level plot
example_color = {[1,0,0],[0,.4,1]}; %ring color per genotype's example fly: {empty, lpsp}

for g = 1:2
    f = example_fly(g);
    if isnan(f); continue; end
    c = example_color{g};
    trial_idx = find(fly_ix == f); %all of this fly's trials, for the trial-level plots

    figure(1)
    subplot(1,2,1); scatter(h_err.XData(trial_idx), h_err.YData(trial_idx), 60,c,'linewidth',1)
    subplot(1,2,2); scatter(h_corr.XData(trial_idx),h_corr.YData(trial_idx),60,c,'linewidth',1)

    figure(2)
    scatter(h_mov.XData(trial_idx),h_mov.YData(trial_idx),60,c,'linewidth',1)

    figure(3)
    for li = 1:2
        subplot(2,2,(li-1)*2+g); hold on
        plot([1,2],[fly_pre(f,li),fly_post(f,li)],'-o','Color',c,'LineWidth',.1,'MarkerFaceColor',c,'MarkerSize',5)
    end

    figure(4)
    for li = 1:2
        subplot(1,2,li); hold on
        scatter(h_delta(li).XData(f),h_delta(li).YData(f),60,c,'linewidth',1)
    end
end

note_str = 'red ring/marker = empty example fly, blue = lpsp example fly';
figure(1); sgtitle({'mu vs fly heading, trials from 7/8 onward',note_str})
figure(2); title({'bump mobility relative to fly turning, trials from 7/8 onward',note_str})
figure(3); sgtitle({'bump mobility before vs after perturbation',note_str})
figure(4); sgtitle({'change in bump mobility after perturbation',note_str})

%% 11) compare basic walking statistics across groups: forward speed, rotational speed, rotational displacement
speed_f = nan(n,1); %mean forward speed over the trial (mm/s)
speed_r = nan(n,1); %mean rotational (yaw) speed over the trial (rad/s)
disp_r  = nan(n,1); %total rotational displacement over the trial (heading path length, rad)

for i = 1:n
    f_speed = abs(all_data(i).ft.f_speed);
    speed_f(i) = mean(f_speed(~isoutlier(f_speed)),'omitnan'); %drop outlier frames (e.g. ball-tracking glitches) before averaging
    speed_r(i) = mean(abs(all_data(i).ft.r_speed),'omitnan');
    disp_r(i)  = sum(abs(diff(unwrap(all_data(i).ft.cue))),'omitnan');
end

speed_f(42) = nan; %hard-coded, not sure why this trial has a monotonically increasing forward speed but it does

fly_speed_f = nan(n_flies,2); %columns: [closed loop, dark]
fly_speed_r = nan(n_flies,2);
fly_disp_r  = nan(n_flies,2);

for f = 1:n_flies
    trial_idx = find(fly_ix == f);
    for li = 1:2
        this_idx = trial_idx(is_dark(trial_idx) == (li==2));
        if isempty(this_idx); continue; end

        fly_speed_f(f,li) = mean(speed_f(this_idx),'omitnan');
        fly_speed_r(f,li) = mean(speed_r(this_idx),'omitnan');
        fly_disp_r(f,li)  = mean(disp_r(this_idx),'omitnan');
    end
end

walk_metric        = {fly_speed_f,fly_speed_r,fly_disp_r};
walk_metric_labels = {'mean forward speed (mm/s)','mean rotational speed (rad/s)','rotational displacement (rad/trial)'};

figure(7); clf
for m = 1:3
    for li = 1:2
        subplot(3,2,(m-1)*2+li); hold on
        scatter(fly_group,walk_metric{m}(:,li),'filled','MarkerFaceAlpha',.5)
        for g = 1:2
            valid = fly_group==g & ~isnan(walk_metric{m}(:,li));
            errorbar(g+.1,mean(walk_metric{m}(valid,li),'omitnan'),std(walk_metric{m}(valid,li),'omitnan')/sqrt(sum(valid)),'ok')
        end
        xticks(1:2); xticklabels(group_labels); xlim([.5,2.5])
        ylabel(walk_metric_labels{m})
        title(lighting_labels{li})
    end
end
sgtitle('walking statistics by genotype and lighting condition')

%% identify the fly with the highest mean forward speed (sanity check for outliers)
[max_speed,max_idx] = max(fly_speed_f(:));
[f,li] = ind2sub(size(fly_speed_f),max_idx);
fprintf('\nhighest mean forward speed: %s fly %i, %s, %.2f mm/s\n',fly_date{f},fly_num(f),lighting_labels{li},max_speed)

%% 12) break the bump accuracy (offset drift) metric into pre/post stim and closed-loop/dark trials, as in section 5
%reuses the same per-fly pre/post trial groupings computed in section 4 (fly_pre_idx/fly_post_idx), just averaging mu_err instead of mov_ratio
fly_err_pre  = nan(n_flies,2); %columns: [closed loop, dark]
fly_err_post = nan(n_flies,2);

for f = 1:n_flies
    for li = 1:2
        pre_this  = fly_pre_idx{f,li};
        post_this = fly_post_idx{f,li};
        if isempty(pre_this) || isempty(post_this); continue; end

        fly_err_pre(f,li)  = mean(mu_err(pre_this),'omitnan');
        fly_err_post(f,li) = mean(mu_err(post_this),'omitnan');
    end
end

figure(8); clf
for g = 1:2
    for li = 1:2
        subplot(2,2,(li-1)*2+g); hold on

        valid = fly_group==g & ~isnan(fly_err_pre(:,li)) & ~isnan(fly_err_post(:,li));
        pre_vals  = fly_err_pre(valid,li);
        post_vals = fly_err_post(valid,li);

        plot([ones(sum(valid),1),2*ones(sum(valid),1)]',[pre_vals,post_vals]','Color',[.5,.5,.5,.3])
        scatter(ones(sum(valid),1), pre_vals,'filled','MarkerFaceAlpha',.5)
        scatter(2*ones(sum(valid),1),post_vals,'filled','MarkerFaceAlpha',.5)
        errorbar([1.1,2.1],[mean(pre_vals,'omitnan'),mean(post_vals,'omitnan')],...
                       [std(pre_vals,'omitnan'),std(post_vals,'omitnan')]/sqrt(sum(valid)),'-ok','LineWidth',1.5)

        xlim([.5,2.5]); xticks([1,2]); xticklabels({'pre','post'})
        ylabel('circular variance of mu-heading offset')
        title(sprintf('%s, %s (n=%i flies)',group_labels{g},lighting_labels{li},sum(valid)))
    end
end
sgtitle('bump accuracy (offset drift) before vs after perturbation')

%% 13) sanity check: average atp and gcamp fluorescence around each stim, for every detected perturbation trial
%ft.stims is the ground-truth ejection marker (logical, on the ft.xf clock) - this overlaps every individual stim
%within a trial (aligned to its onset) and averages the atp (atp.d) and gcamp (im.d) signal in a +/-5s window around it
peri_win = 5; %seconds before/after each stim to show
t_common = linspace(-peri_win,peri_win,101); %common relative-time axis stims get interpolated onto

pulse_trials = find(has_pulse);
n_rows = ceil(sqrt(length(pulse_trials)));
n_cols = ceil(length(pulse_trials)/n_rows);

figure(9); clf
for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets  = find(diff([false;stims])==1); %rising edge of each stim
    t_onset = xf(onsets);

    subplot(n_rows,n_cols,k); hold on
    f = fly_ix(i);
    title(sprintf('%s fly %i, trial %i (n=%i stims)',fly_date{f},fly_num(f),trial_num(i),length(onsets)))

    if isempty(onsets)
        axis off
        continue
    end

    atp_sig   = sum(all_data(i).atp.d,1);
    gcamp_sig = sum(all_data(i).im.d,1);

    atp_aligned   = nan(length(onsets),length(t_common));
    gcamp_aligned = nan(length(onsets),length(t_common));
    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s); %imaging clock, re-centered on this stim's onset
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
sgtitle('average atp (red, left axis) and gcamp (green, right axis) fluorescence around each stim')

%% 14) show the decision behind each detected perturbation trial: peri-stim average atp trace vs. its own baseline + threshold
%reuses atp_peri (computed once in section 3) rather than recomputing the alignment - this is exactly what fed has_pulse
%recomputes pulse_trials/n_rows/n_cols here (rather than reusing section 13's) so this figure can't go stale relative to
%has_pulse if section 13 isn't re-run in the same pass - e.g. after re-running section 3 with different detection settings
pulse_trials = find(has_pulse);
n_rows = ceil(sqrt(length(pulse_trials)));
n_cols = ceil(length(pulse_trials)/n_rows);

figure(10); clf
for k = 1:length(pulse_trials)
    i = pulse_trials(k);
    f = fly_ix(i);

    base_idx = t_common >= base_win(1) & t_common <  base_win(2);
    peak_idx = t_common >  peak_win(1) & t_common <= peak_win(2);
    base_mean = median(atp_peri{i}(base_idx),'omitnan');
    base_std  = std(atp_peri{i}(base_idx),'omitnan');
    thresh    = base_mean + peak_factor*base_std;

    subplot(n_rows,n_cols,k); hold on

    yl = [min(atp_peri{i}),max([atp_peri{i},thresh])];
    yl = yl + [-1,1]*.1*max(diff(yl),eps);

    patch(base_win([1,2,2,1]),yl([1,1,2,2]),'k','FaceAlpha',.08,'EdgeColor','none') %pre-stim baseline window
    patch(peak_win([1,2,2,1]),yl([1,1,2,2]),'g','FaceAlpha',.08,'EdgeColor','none') %post-stim peak window

    plot(t_common,atp_peri{i},'r','LineWidth',1.5)
    plot([-peri_win,peri_win],[base_mean,base_mean],':k')
    plot([-peri_win,peri_win],[thresh,thresh],':r')
    plot([0,0],yl,':k')

    xlim([-peri_win,peri_win]); ylim(yl)
    title(sprintf('%s fly %i, trial %i',fly_date{f},fly_num(f),trial_num(i)))
end
sgtitle('perturbation-trial detection: peri-stim atp average (red), baseline window (gray), peak window (green), threshold (dashed red)')

%% local functions
function h = plot_sem(ax,t,x)
%shades mean(x) +/- sem(x) over rows of x (one row per replicate, columns matching t); face color is set by the caller
t = reshape(t,1,[]);
m1 = mean(x,1,'omitnan');
s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

idx = ~isnan(m1);
m1 = m1(idx);
s1 = s1(idx);
t  = t(idx);

h = patch(ax,[t,fliplr(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end
