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
bout_mov_mu   = cell(n,1); %per-bout bump path length, this trial's raw data behind mov_ratio's regression
bout_mov_cue  = cell(n,1); %per-bout heading path length, same
bout_dur      = cell(n,1); %duration of each bout (s) - used to weight the regression (longer bouts = more reliable = more weight)

for i = 1:n
    xf = all_data(i).ft.xf;
    dt = median(diff(xf));

    mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');

    %use r_speed (the ball's actual rotation) rather than cue position - a handful of trials change the visual
    %cue experimentally while the fly isn't walking, which would otherwise show up as spurious "heading" movement
    r_speed_smooth = smoothdata(all_data(i).ft.r_speed(:),'gaussian',smooth_window); %force column - unlike cue, r_speed's native orientation isn't guaranteed by an upstream transpose
    mu_smooth      = smoothdata(mu,'gaussian',smooth_window);
    fly_speed      = abs(r_speed_smooth); %already a speed - no differentiation needed, unlike mu/cue position

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    mov_mu  = nan(length(bout_starts),1);
    mov_cue = nan(length(bout_starts),1); %despite the name (kept for consistency downstream), this is now the fly's r_speed-based path length, not cue position
    dur     = nan(length(bout_starts),1);
    for b = 1:length(bout_starts)
        mov_mu(b)  = sum(abs(diff(mu_smooth(bout_starts(b):bout_ends(b)))),'omitnan');
        mov_cue(b) = sum(abs(r_speed_smooth(bout_starts(b):bout_ends(b))),'omitnan')*dt; %path length = integral of |speed| dt
        dur(b)     = (bout_ends(b)-bout_starts(b)+1)*dt;
    end

    w = sqrt(dur); %weighted least squares via the standard sqrt(weight) rescaling trick: (sqrt(w).*x)\(sqrt(w).*y)
    mov_ratio(i)    = (w.*mov_cue) \ (w.*mov_mu); %bump path length per unit heading path length, longer bouts weighted more heavily
    frac_walking(i) = mean(is_walking,'omitnan');
    mean_rho(i)     = mean(all_data(i).im.rho,'omitnan');
    bout_mov_mu{i}  = mov_mu;
    bout_mov_cue{i} = mov_cue;
    bout_dur{i}     = dur;
end

%% histogram of bout counts per trial - use this to pick an empirical minimum-bout threshold below
n_bouts_per_trial = cellfun(@length,bout_mov_cue);

figure(14); clf
histogram(n_bouts_per_trial,'BinMethod','integers')
xlabel('number of walking bouts in a trial')
ylabel('number of trials')
title('distribution of bout counts per trial')

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

%% how well-powered is the bump-mobility regression? show the exact scatter + fit behind mov_ratio for single trials
%mov_ratio(i) = (sqrt(dur).*cue)\(sqrt(dur).*mu) - this shows that same per-trial weighted regression directly, not pooled
n_example_trials = 6;
[~,sort_ix] = sort(n_bouts_per_trial);
example_ranks  = round(linspace(1,n,n_example_trials)); %span the full range of bout counts, from least- to best-powered
example_trials = sort_ix(example_ranks);

figure(12); clf
for k = 1:n_example_trials
    i = example_trials(k);
    x = bout_mov_cue{i};
    y = bout_mov_mu{i};
    dur = bout_dur{i};
    slope = mov_ratio(i); %the exact value this trial contributes to the rest of the analysis

    subplot(2,3,k); hold on
    scatter(x,y,10+40*dur/max(dur),'filled','MarkerFaceAlpha',.5) %marker area scales with bout duration (the regression weight)
    xl = xlim; xl(1) = 0;
    plot(xl,xl*slope,'r','LineWidth',1.5)
    plot(xl,xl,':k') %reference: bump moves exactly as much as heading (slope = 1)
    xlim(xl)

    xlabel('heading path length per bout (rad)')
    ylabel('bump path length per bout (rad)')
    f = fly_ix(i);
    title(sprintf('%s fly %i, trial %i (n=%i bouts, slope=%.2f)',fly_date{f},fly_num(f),trial_num(i),n_bouts_per_trial(i),slope))
end
sgtitle('bump-mobility regression for single trials - spanning least- to best-powered')

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
min_bouts_per_condition = 15; %minimum pooled bout count (per fly x lighting x pre/post cell) to trust its mov_ratio estimate - set from the bout-count histogram above

fly_group = nan(n_flies,1);            %1 = empty, 2 = lpsp
fly_pre   = nan(n_flies,2);            %columns: [closed loop, dark]
fly_post  = nan(n_flies,2);

fly_pre_idx  = cell(n_flies,2);        %the actual trial indices (into all_data) behind fly_pre/fly_post,
fly_post_idx = cell(n_flies,2);        %kept so we can go back and plot specific example trials later

n_excluded_low_bouts = 0; %how many fly x lighting x pre/post cells got dropped for having too few pooled bouts

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

        %pool every walking bout across all trials in this condition and fit one regression, rather than
        %averaging each trial's own (often underpowered) slope - this is the same fit mov_ratio(i) uses, just
        %on the combined bout set, so each fly x lighting x pre/post cell gets as much data as it actually has
        pre_cue  = vertcat(bout_mov_cue{pre_this});
        pre_mu   = vertcat(bout_mov_mu{pre_this});
        pre_dur  = vertcat(bout_dur{pre_this});
        post_cue = vertcat(bout_mov_cue{post_this});
        post_mu  = vertcat(bout_mov_mu{post_this});
        post_dur = vertcat(bout_dur{post_this});

        if length(pre_cue) >= min_bouts_per_condition
            w = sqrt(pre_dur);
            fly_pre(f,li) = (w.*pre_cue) \ (w.*pre_mu);
        else
            n_excluded_low_bouts = n_excluded_low_bouts + ~isempty(pre_cue);
        end
        if length(post_cue) >= min_bouts_per_condition
            w = sqrt(post_dur);
            fly_post(f,li) = (w.*post_cue) \ (w.*post_mu);
        else
            n_excluded_low_bouts = n_excluded_low_bouts + ~isempty(post_cue);
        end

        fly_pre_idx{f,li}  = pre_this;
        fly_post_idx{f,li} = post_this;
    end
end

fprintf('excluded %i fly x lighting x pre/post cell(s) for having fewer than %i pooled bouts\n',n_excluded_low_bouts,min_bouts_per_condition)

fprintf('%i of %i flies had a detectable perturbation trial\n',sum(any(~isnan(fly_pre) | ~isnan(fly_post),2)),n_flies)

%% show the pooled scatter + regression behind a few example fly x lighting x pre/post pools
%same idea as the single-trial figure above, but now for the pooled fly_pre/fly_post cells themselves.
%marker size scales with each bout's duration, to visualize the weighting applied to the regression.
%only cells that passed the min_bouts_per_condition threshold (i.e. actually fed fly_pre/fly_post) are shown here
pool_list = {}; %columns: {cue, mu, dur, slope, label}

for f = 1:n_flies
    for li = 1:2
        if ~isnan(fly_pre(f,li))
            c = vertcat(bout_mov_cue{fly_pre_idx{f,li}});
            m = vertcat(bout_mov_mu{fly_pre_idx{f,li}});
            dur = vertcat(bout_dur{fly_pre_idx{f,li}});
            pool_list(end+1,:) = {c,m,dur,fly_pre(f,li),sprintf('%s fly %i, %s, pre',fly_date{f},fly_num(f),lighting_labels{li})}; %#ok<AGROW>
        end
        if ~isnan(fly_post(f,li))
            c = vertcat(bout_mov_cue{fly_post_idx{f,li}});
            m = vertcat(bout_mov_mu{fly_post_idx{f,li}});
            dur = vertcat(bout_dur{fly_post_idx{f,li}});
            pool_list(end+1,:) = {c,m,dur,fly_post(f,li),sprintf('%s fly %i, %s, post',fly_date{f},fly_num(f),lighting_labels{li})}; %#ok<AGROW>
        end
    end
end

n_bouts_per_pool = cellfun(@length,pool_list(:,1));

n_example_pools = 6;
[~,sort_ix] = sort(n_bouts_per_pool);
example_ranks = round(linspace(1,size(pool_list,1),n_example_pools)); %span the full range of pool sizes, least- to best-powered
example_pools = sort_ix(example_ranks);

figure(13); clf
for k = 1:n_example_pools
    p = example_pools(k);
    x = pool_list{p,1};
    y = pool_list{p,2};
    dur = pool_list{p,3};
    slope = pool_list{p,4};

    subplot(2,3,k); hold on
    scatter(x,y,10+40*dur/max(dur),'filled','MarkerFaceAlpha',.5) %marker area scales with bout duration (the regression weight)
    xl = xlim; xl(1) = 0;
    plot(xl,xl*slope,'r','LineWidth',1.5)
    plot(xl,xl,':k') %reference: bump moves exactly as much as heading (slope = 1)
    xlim(xl)

    xlabel('heading path length per bout (rad)')
    ylabel('bump path length per bout (rad)')
    title(sprintf('%s (n=%i bouts, slope=%.2f)',pool_list{p,5},n_bouts_per_pool(p),slope))
end
sgtitle('bump-mobility regression, pooled per fly x lighting x pre/post condition - spanning least- to best-powered')

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
linkaxes(get(gcf,"Children"),'y')

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

    median_delta = median(delta(fly_group==g,:),1,'omitnan'); %group-median [closed loop, dark] delta
    dist = sqrt(sum((delta(good_candidates,:) - median_delta).^2,2));
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
                fly_val = fly_pre(f,li); %the pooled, across-trial value shown in figures 3/4 for this cell - not the same thing as this one trial's own mov_ratio
            else
                cand = fly_post_idx{f,li};
                phase_label = 'post';
                fly_val = fly_post(f,li);
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
            title(sprintf('%s, %s (trial %i)\ntrial mov ratio=%.2f, pooled (n=%i trials) mov ratio=%.2f',...
                lighting_labels{li},phase_label,trial_num(i),mov_ratio(i),length(cand),fly_val))
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

%% 15) sensitivity check: how does bump mobility change with different levels of mu smoothing?
%this is exploratory - it does NOT change mov_ratio/fly_pre/fly_post used elsewhere in the script, it just re-runs the
%same trial-level pipeline (bout detection + weighted regression) under different mu-smoothing choices, including a
%version where the smoothing kernel is weighted by rho (bump vector strength / decoding confidence) instead of
%treating every timepoint as equally trustworthy - low-rho samples contribute less to the local average
smoothing_configs = {
    20,  false, 'window=20, unweighted';
    60,  false, 'window=60 (current default), unweighted';
    120, false, 'window=120, unweighted';
    240, false, 'window=240, unweighted';
    60,  true,  'window=60, rho-weighted';
    120, true,  'window=120, rho-weighted';
};

n_configs = size(smoothing_configs,1);
mov_ratio_sweep = nan(n,n_configs);

for c = 1:n_configs
    fprintf('running smoothing sweep %i of %i: %s\n',c,n_configs,smoothing_configs{c,3})
    mov_ratio_sweep(:,c) = bump_mobility_sweep(all_data,smoothing_configs{c,1},turn_thresh,max_gap_frames,min_walking_frames,smoothing_configs{c,2});
end

empty_idx    = group_idx==1; %restrict to empty flies only - no need to bring lpsp into this smoothing check
lighting_idx = is_dark(:)+1; %1 = closed loop, 2 = dark

figure(15); clf
for c = 1:n_configs
    subplot(2,3,c); hold on
    swarmchart(lighting_idx(empty_idx),mov_ratio_sweep(empty_idx,c),'filled','MarkerFaceAlpha',.5,'XJitterWidth',.3)
    for li = 1:2
        idx = empty_idx' & is_dark(:)==(li==2);
        errorbar(li,mean(mov_ratio_sweep(idx,c),'omitnan'),std(mov_ratio_sweep(idx,c),'omitnan')/sqrt(sum(idx)),'ok')
    end
    plot(xlim,[1,1],':k')
    xticks(1:2); xticklabels(lighting_labels); xlim([.5,2.5])
    ylabel('bump path length / heading path length')
    title(smoothing_configs{c,3})
end
sgtitle('bump mobility (empty flies only) across different mu-smoothing choices')

fprintf('\n%-35s %-12s %-12s\n','smoothing config','closed loop','dark')
for c = 1:n_configs
    fprintf('%-35s %-12.3f %-12.3f\n',smoothing_configs{c,3},...
        mean(mov_ratio_sweep(empty_idx' & ~is_dark(:),c),'omitnan'),mean(mov_ratio_sweep(empty_idx' & is_dark(:),c),'omitnan'))
end

%% 16) directional bump mobility: does the bump move differently during net-positive vs net-negative turns?
%note: this doesn't touch mov_ratio/fly_pre/fly_post computed earlier - every existing metric uses abs(), so it's
%blind to a global sign flip and none of it needed to change. this section only builds new, signed quantities.

%% 16a) determine which side of the PB was stimulated, per fly
%convention matches the codebase's existing right_idx logic (e.g. lpsp_p2x2_script_minimal.m / gain_metric_sandbox.m):
%centroids 1:16 brighter than 17:32 during the post-stim peak window = "right". we don't independently verify this
%matches true anatomical left/right - what matters is that it's applied consistently across every fly.
fly_flip_sign = nan(n_flies,1); %+1 = right-stimulated (reference, no flip), -1 = left-stimulated (flip to align)

for f = 1:n_flies
    trial_idx     = find(fly_ix == f);
    pulse_trials_f = trial_idx(has_pulse(trial_idx));
    if isempty(pulse_trials_f); continue; end %can't determine stim side without a detected pulse trial

    hemi1 = nan(length(pulse_trials_f),1); %mean peak-window atp signal, centroids 1:16
    hemi2 = nan(length(pulse_trials_f),1); %mean peak-window atp signal, centroids 17:32

    for k = 1:length(pulse_trials_f)
        i = pulse_trials_f(k);
        stims   = logical(all_data(i).ft.stims(:));
        onsets  = find(diff([false;stims])==1);
        xb      = all_data(i).ft.xb;
        t_onset = all_data(i).ft.xf(onsets);

        h1 = nan(length(onsets),1);
        h2 = nan(length(onsets),1);
        for s = 1:length(onsets)
            in_peak = (xb-t_onset(s)) > peak_win(1) & (xb-t_onset(s)) <= peak_win(2);
            h1(s) = mean(all_data(i).atp.d(1:16,in_peak),'all');
            h2(s) = mean(all_data(i).atp.d(17:32,in_peak),'all');
        end
        hemi1(k) = mean(h1,'omitnan');
        hemi2(k) = mean(h2,'omitnan');
    end

    fly_flip_sign(f) = 2*(mean(hemi1,'omitnan') > mean(hemi2,'omitnan')) - 1; %+1 if hemi1 (1:16) brighter, else -1
end

fprintf('\nstim side: %i right (flip=+1), %i left (flip=-1), %i undetermined (no detected pulse trial)\n',...
    sum(fly_flip_sign==1),sum(fly_flip_sign==-1),sum(isnan(fly_flip_sign)))

%% 16b) classify each walking bout as a net-positive or net-negative turn, in the aligned (flipped) reference frame
%sign convention assumed here: positive raw r_speed = counterclockwise (FicTrac's usual velYaw convention) -
%if this rig's convention is the opposite, just swap the "CW"/"CCW" labels below, the analysis itself is unaffected
bout_turn_sign = cell(n,1);

for i = 1:n
    xf = all_data(i).ft.xf;
    dt = median(diff(xf));

    r_speed_smooth = smoothdata(all_data(i).ft.r_speed(:),'gaussian',smooth_window);
    fly_speed = abs(r_speed_smooth);

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames+1,1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false;is_walking;false]);
    bout_starts = find(d==1);
    bout_ends   = find(d==-1)-1;

    flip = fly_flip_sign(fly_ix(i)); %nan if this fly's stim side is undetermined

    turn_sign = nan(length(bout_starts),1);
    for b = 1:length(bout_starts)
        turn_sign(b) = sign(sum(r_speed_smooth(bout_starts(b):bout_ends(b)),'omitnan')*flip);
    end
    bout_turn_sign{i} = turn_sign; %same length/order as bout_mov_mu{i}/bout_mov_cue{i}/bout_dur{i} from section 2
end

%% 16c) pool bouts by fly x lighting x pre/post x turn direction, and fit mobility separately for each direction
dir_labels = {'CCW','CW'}; %dirn=1 -> turn_sign>0, dirn=2 -> turn_sign<0 (see the sign-convention note above)

fly_pre_dir  = nan(n_flies,2,2); %dims: [fly, lighting, direction]
fly_post_dir = nan(n_flies,2,2);

for f = 1:n_flies
    if isnan(fly_flip_sign(f)); continue; end %can't align this fly without a known stim side

    for li = 1:2
        for phase = 1:2
            if phase==1; cand = fly_pre_idx{f,li}; else; cand = fly_post_idx{f,li}; end
            if isempty(cand); continue; end

            sign_all = vertcat(bout_turn_sign{cand});
            cue_all  = vertcat(bout_mov_cue{cand});
            mu_all   = vertcat(bout_mov_mu{cand});
            dur_all  = vertcat(bout_dur{cand});

            for dirn = 1:2
                this_sign = 3-2*dirn; %dirn=1 -> +1 (CCW), dirn=2 -> -1 (CW)
                idx = sign_all==this_sign;
                if sum(idx) < min_bouts_per_condition; continue; end

                w = sqrt(dur_all(idx));
                slope = (w.*cue_all(idx)) \ (w.*mu_all(idx));

                if phase==1; fly_pre_dir(f,li,dirn) = slope; else; fly_post_dir(f,li,dirn) = slope; end
            end
        end
    end
end

%% 16d) plot paired pre/post directional mobility, by lighting and turn direction, colored by genotype
geno_color = {[1,0,0],[0,.4,1]}; %empty = red, lpsp = blue

figure(16); clf
for dirn = 1:2
    for li = 1:2
        subplot(2,2,(dirn-1)*2+li); hold on
        for g = 1:2
            valid = fly_group==g & ~isnan(fly_pre_dir(:,li,dirn)) & ~isnan(fly_post_dir(:,li,dirn));
            c = geno_color{g};
            pre_vals  = fly_pre_dir(valid,li,dirn);
            post_vals = fly_post_dir(valid,li,dirn);

            plot([ones(sum(valid),1),2*ones(sum(valid),1)]',[pre_vals,post_vals]','Color',[c,.3])
            scatter(ones(sum(valid),1), pre_vals,'filled','MarkerFaceColor',c,'MarkerFaceAlpha',.5)
            scatter(2*ones(sum(valid),1),post_vals,'filled','MarkerFaceColor',c,'MarkerFaceAlpha',.5)
        end
        xlim([.5,2.5]); xticks([1,2]); xticklabels({'pre','post'})
        ylabel('bump path length / heading path length')
        title(sprintf('%s turns, %s',dir_labels{dirn},lighting_labels{li}))
    end
end
sgtitle({'directional bump mobility before vs after perturbation (red=empty, blue=lpsp)',...
    sprintf('all flies aligned as if stimulated on the same (reference) side of the PB - %i of %i flies flipped',...
        sum(fly_flip_sign==-1),sum(~isnan(fly_flip_sign)))})

%% 17) pool EPG calcium (gcamp) around each stim, empty vs lpsp, split by stimulated vs non-stimulated hemisphere
%same peri-stim alignment as figure 9, but now: (1) split the 32 centroids into stim-side/non-stim-side halves using
%each fly's stim-side call from section 16a, and (2) baseline-subtract each stim's trace to 0 at the rising edge (t=0)
pulse_trials = find(has_pulse); %recomputed fresh so this can't go stale
[~,zero_idx] = min(abs(t_common)); %index of t=0 within t_common

gcamp_pool_stim    = {[],[]}; %row per individual stim; column 1 = empty, column 2 = lpsp
gcamp_pool_nonstim = {[],[]};

for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets = find(diff([false;stims])==1);
    if isempty(onsets); continue; end
    t_onset = xf(onsets);

    if fly_flip_sign(fly_ix(i)) == 1
        stim_rows = 1:16; nonstim_rows = 17:32;
    else
        stim_rows = 17:32; nonstim_rows = 1:16;
    end
    gcamp_stim_sig    = sum(all_data(i).im.d(stim_rows,:),1);
    gcamp_nonstim_sig = sum(all_data(i).im.d(nonstim_rows,:),1);

    stim_aligned    = nan(length(onsets),length(t_common));
    nonstim_aligned = nan(length(onsets),length(t_common));
    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s);
        stim_aligned(s,:)    = interp1(xb_rel,gcamp_stim_sig,   t_common,'linear',nan);
        nonstim_aligned(s,:) = interp1(xb_rel,gcamp_nonstim_sig,t_common,'linear',nan);
    end

    stim_aligned    = stim_aligned    - stim_aligned(:,zero_idx);    %baseline-subtract: each stim's own trace starts at 0
    nonstim_aligned = nonstim_aligned - nonstim_aligned(:,zero_idx);

    g = group_idx(i);
    gcamp_pool_stim{g}    = [gcamp_pool_stim{g};stim_aligned];
    gcamp_pool_nonstim{g} = [gcamp_pool_nonstim{g};nonstim_aligned];
end

figure(17); clf
panel_data   = {gcamp_pool_stim,gcamp_pool_nonstim};
panel_labels = {'stimulated hemisphere','non-stimulated hemisphere'};

for p = 1:2
    subplot(1,2,p); hold on
    h = gobjects(1,2);
    for g = 1:2
        if isempty(panel_data{p}{g}); continue; end
        h(g) = plot_sem(gca,t_common,panel_data{p}{g});
        h(g).FaceColor = geno_color{g};
    end
    plot([0,0],ylim,':k')
    plot(xlim,[0,0],':k')
    xlim([-peri_win,peri_win])
    xlabel('time from stim onset (s)')
    ylabel('EPG calcium (summed gcamp dF/F, baseline-subtracted at t=0)')
    legend(h,group_labels,'AutoUpdate','off','Location','best')
    title(sprintf('%s (%i empty stims, %i lpsp stims)',panel_labels{p},size(panel_data{p}{1},1),size(panel_data{p}{2},1)))
end
sgtitle('EPG calcium around stim onset, by hemisphere')

%% 18) average glomerulus heatmap around stims, empty vs lpsp - same full-PB flip as above to align left-stim to right-stim
%flipping the entire 32-row centroid stack (not just the summed signal) puts rows 1:16 as the stimulated hemisphere
%and 17:32 as the non-stimulated hemisphere consistently for every fly, regardless of which side was really stimulated
heatmap_pool = {{},{}}; %cell per genotype; each entry is one stim's 32 x length(t_common) baseline-subtracted heatmap

for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets = find(diff([false;stims])==1);
    if isempty(onsets); continue; end
    t_onset = xf(onsets);

    im_d = all_data(i).im.d; %32 x n_frames
    if fly_flip_sign(fly_ix(i)) == -1
        im_d = flipud(im_d); %align this left-stimulated fly to the right-stimulated reference frame
    end

    g = group_idx(i);
    for s = 1:length(onsets)
        xb_rel  = xb - t_onset(s);
        aligned = interp1(xb_rel,im_d',t_common,'linear',nan)'; %32 x length(t_common)
        aligned = aligned - aligned(:,zero_idx); %baseline-subtract each glomerulus row to 0 at t=0

        heatmap_pool{g}{end+1} = aligned; %#ok<AGROW>
    end
end

mean_heatmap = cell(1,2);
for g = 1:2
    if ~isempty(heatmap_pool{g})
        mean_heatmap{g} = mean(cat(3,heatmap_pool{g}{:}),3,'omitnan');
    end
end

clim_max = max(abs([mean_heatmap{1}(:);mean_heatmap{2}(:)]),[],'omitnan');
n_col = 256;
diverge_cmap = [linspace(0,1,n_col/2)',linspace(0,1,n_col/2)',ones(n_col/2,1);... %blue -> white
                ones(n_col/2,1),linspace(1,0,n_col/2)',linspace(1,0,n_col/2)'];   %white -> red

figure(18); clf
for g = 1:2
    subplot(1,2,g); hold on
    if isempty(mean_heatmap{g}); continue; end

    imagesc(t_common,1:32,mean_heatmap{g})
    plot([-peri_win,peri_win],[16.5,16.5],'k:','LineWidth',1) %hemisphere boundary
    plot([0,0],[.5,32.5],'k:','LineWidth',1) %stim onset
    colormap(gca,diverge_cmap)
    caxis([-clim_max,clim_max])
    axis tight
    yticks([8.5,24.5]); yticklabels({'stimulated hemisphere','non-stim hemisphere'})
    xlabel('time from stim onset (s)')
    c = colorbar; c.Label.String = 'dF/F (baseline-subtracted)';
    title(sprintf('%s (n=%i stims)',group_labels{g},length(heatmap_pool{g})))
end
sgtitle('average glomerulus dF/F around stim onset, aligned to a common stim side')

%% 19) average bump position (mu) trace around stims, empty vs lpsp - same flip and baseline-subtraction as above
%mu is unwrapped first (so the sign flip and later baseline subtraction operate on a continuous trace, not a
%wrapped one), then negated for left-stimulated flies, then anchored to 0 at the stim onset for each individual stim
mu_pool = {[],[]}; %row per individual stim; column 1 = empty, column 2 = lpsp

for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets = find(diff([false;stims])==1);
    if isempty(onsets); continue; end
    t_onset = xf(onsets);

    mu_sig = unwrap(all_data(i).im.mu);
    if fly_flip_sign(fly_ix(i)) == -1
        mu_sig = -mu_sig; %align this left-stimulated fly to the right-stimulated reference frame
    end

    mu_aligned = nan(length(onsets),length(t_common));
    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s);
        mu_aligned(s,:) = interp1(xb_rel,mu_sig,t_common,'linear',nan);
    end

    mu_aligned = mu_aligned - mu_aligned(:,zero_idx); %baseline-subtract: each stim's mu starts at 0 at t=0

    g = group_idx(i);
    mu_pool{g} = [mu_pool{g};mu_aligned];
end

figure(19); clf; hold on
h = gobjects(1,2);
for g = 1:2
    if isempty(mu_pool{g}); continue; end
    h(g) = plot_sem(gca,t_common,mu_pool{g});
    h(g).FaceColor = geno_color{g};
end
plot([0,0],ylim,':k')
plot(xlim,[0,0],':k')
xlim([-peri_win,peri_win])
xlabel('time from stim onset (s)')
ylabel('bump position (rad, relative to stim onset)')
legend(h,group_labels,'AutoUpdate','off','Location','best')
title(sprintf('bump position around stim onset - %i empty stims, %i lpsp stims',size(mu_pool{1},1),size(mu_pool{2},1)))

%% 20) same bump-position analysis, but keep left- and right-stimulated flies separate instead of flipping/pooling them
%a check on whether the flip in section 19 is doing something sensible: left- and right-stim flies should look like
%mirror images of each other here, since no flip is applied - mu is raw/unflipped in this section
mu_pool_side = {{[],[]},{[],[]}}; %{side}{genotype}: side 1 = left-stimulated, side 2 = right-stimulated

for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets = find(diff([false;stims])==1);
    if isempty(onsets); continue; end
    t_onset = xf(onsets);

    mu_sig = unwrap(all_data(i).im.mu); %raw, unflipped - that's the whole point of this comparison

    mu_aligned = nan(length(onsets),length(t_common));
    for s = 1:length(onsets)
        xb_rel = xb - t_onset(s);
        mu_aligned(s,:) = interp1(xb_rel,mu_sig,t_common,'linear',nan);
    end

    mu_aligned = mu_aligned - mu_aligned(:,zero_idx); %baseline-subtract: each stim's mu starts at 0 at t=0

    side = 1 + (fly_flip_sign(fly_ix(i))==1); %1 = left-stimulated, 2 = right-stimulated
    g    = group_idx(i);
    mu_pool_side{side}{g} = [mu_pool_side{side}{g};mu_aligned];
end

side_labels = {'left-stimulated flies','right-stimulated flies'};

figure(20); clf
for side = 1:2
    subplot(1,2,side); hold on
    h = gobjects(1,2);
    for g = 1:2
        if isempty(mu_pool_side{side}{g}); continue; end
        h(g) = plot_sem(gca,t_common,mu_pool_side{side}{g});
        h(g).FaceColor = geno_color{g};
    end
    plot([0,0],ylim,':k')
    plot(xlim,[0,0],':k')
    xlim([-peri_win,peri_win])
    xlabel('time from stim onset (s)')
    ylabel('bump position (rad, relative to stim onset)')
    legend(h,group_labels,'AutoUpdate','off','Location','best')
    title(sprintf('%s (%i empty stims, %i lpsp stims)',side_labels{side},size(mu_pool_side{side}{1},1),size(mu_pool_side{side}{2},1)))
end
sgtitle('bump position around stim onset, left- vs right-stimulated flies kept separate (no PB flip applied)')

%% 21) same average glomerulus heatmap as figure 18, but keep left- and right-stimulated flies separate (no PB flip)
%within a single stim-side bucket, flies are already mutually aligned (they all had the same physical side
%stimulated), so unlike figure 18 no row-flip is needed here - only the y-tick labels differ between the two rows,
%since which physical centroid block (1:16 vs 17:32) counts as "stimulated" swaps between left- and right-stim flies
heatmap_pool_side = {{{},{}},{{},{}}}; %{side}{genotype}: side 1 = left-stimulated, side 2 = right-stimulated

for k = 1:length(pulse_trials)
    i = pulse_trials(k);

    xb    = all_data(i).ft.xb;
    xf    = all_data(i).ft.xf;
    stims = logical(all_data(i).ft.stims(:));

    onsets = find(diff([false;stims])==1);
    if isempty(onsets); continue; end
    t_onset = xf(onsets);

    im_d = all_data(i).im.d; %raw, unflipped - that's the whole point of this comparison

    side = 1 + (fly_flip_sign(fly_ix(i))==1); %1 = left-stimulated, 2 = right-stimulated
    g    = group_idx(i);

    for s = 1:length(onsets)
        xb_rel  = xb - t_onset(s);
        aligned = interp1(xb_rel,im_d',t_common,'linear',nan)'; %32 x length(t_common)
        aligned = aligned - aligned(:,zero_idx); %baseline-subtract each glomerulus row to 0 at t=0

        heatmap_pool_side{side}{g}{end+1} = aligned; %#ok<AGROW>
    end
end

mean_heatmap_side = cell(2,2);
for side = 1:2
    for g = 1:2
        if ~isempty(heatmap_pool_side{side}{g})
            mean_heatmap_side{side,g} = mean(cat(3,heatmap_pool_side{side}{g}{:}),3,'omitnan');
        end
    end
end

clim_vals = [];
for side = 1:2
    for g = 1:2
        if ~isempty(mean_heatmap_side{side,g})
            clim_vals = [clim_vals;mean_heatmap_side{side,g}(:)]; %#ok<AGROW>
        end
    end
end
clim_max_side = max(abs(clim_vals),[],'omitnan');

hemi_labels = {{'non-stim hemisphere','stimulated hemisphere'},{'stimulated hemisphere','non-stim hemisphere'}}; %{side}{1:16 label, 17:32 label}

figure(21); clf
for side = 1:2
    for g = 1:2
        subplot(2,2,(side-1)*2+g); hold on
        if isempty(mean_heatmap_side{side,g}); continue; end

        imagesc(t_common,1:32,mean_heatmap_side{side,g})
        plot([-peri_win,peri_win],[16.5,16.5],'k:','LineWidth',1) %hemisphere boundary
        plot([0,0],[.5,32.5],'k:','LineWidth',1) %stim onset
        colormap(gca,diverge_cmap)
        caxis([-clim_max_side,clim_max_side])
        axis tight
        yticks([8.5,24.5]); yticklabels(hemi_labels{side})
        xlabel('time from stim onset (s)')
        c = colorbar; c.Label.String = 'dF/F (baseline-subtracted)';
        title(sprintf('%s, %s (n=%i stims)',side_labels{side},group_labels{g},length(heatmap_pool_side{side}{g})))
    end
end
sgtitle('average glomerulus dF/F around stim onset, left- vs right-stimulated flies kept separate (no PB flip applied)')

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

function mov_ratio_out = bump_mobility_sweep(all_data,smooth_window,turn_thresh,max_gap_frames,min_walking_frames,weight_by_rho)
    %re-runs the trial-level bump-mobility pipeline (section 2) under an arbitrary mu-smoothing choice, for the
    %smoothing sensitivity check - kept separate from section 2 itself so exploring this doesn't disturb the
    %bout_mov_mu/bout_mov_cue/bout_dur values that the rest of the script's fly-level analysis depends on
    n_trials = length(all_data);
    mov_ratio_out = nan(n_trials,1);

    for i = 1:n_trials
        xf = all_data(i).ft.xf;
        xb = all_data(i).ft.xb;
        dt = median(diff(xf));

        mu  = interp1(xb,unwrap(all_data(i).im.mu),xf,'linear','extrap');
        rho = interp1(xb,all_data(i).im.rho,xf,'linear','extrap');

        if weight_by_rho
            mu_smooth = weighted_gauss_smooth(mu,rho,smooth_window);
        else
            mu_smooth = smoothdata(mu,'gaussian',smooth_window);
        end

        r_speed_smooth = smoothdata(all_data(i).ft.r_speed(:),'gaussian',smooth_window);
        fly_speed = abs(r_speed_smooth);

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

        w = sqrt(dur);
        mov_ratio_out(i) = (w.*mov_cue) \ (w.*mov_mu);
    end
end

function mu_smooth = weighted_gauss_smooth(mu,rho,window)
    %gaussian-kernel smoothing weighted by rho (bump vector strength), so low-confidence decoded headings
    %contribute less to the local average - a weighted analogue of smoothdata(mu,'gaussian',window)
    sigma = window/5; %matches smoothdata's convention: window = 5*sigma
    half  = ceil(3*sigma);
    x = -half:half;
    kernel = exp(-x.^2/(2*sigma^2));
    kernel = kernel/sum(kernel);

    mu  = mu(:);
    w   = rho(:);
    num = conv(mu.*w,kernel,'same');
    den = conv(w,kernel,'same');
    mu_smooth = num./den;
end
