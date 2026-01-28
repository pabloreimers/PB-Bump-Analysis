%% load in both datas

th_filename = uigetfile('.data/*th*.mat');
th_data     = load(['.data/',th_filename]);
vg_filename = uigetfile('.data/*vg*.mat');
vg_data     = load(['.data/',vg_filename]);

all_data    = cat(2,th_data.all_data,vg_data.all_data);

%% create meta data indexes
trial_num   = zeros(length(all_data),1);
dark_idx    = false(length(all_data),1);
empty_idx   = false(length(all_data),1);
mcherry_idx = false(length(all_data),1);
vglut_idx   = false(length(all_data),1);
walk_idx    = false(length(all_data),1);
cue_idx     = false(length(all_data),1);
rho_idx     = false(length(all_data),1);
fly_num     = nan(length(all_data),1);
last_id = '';
fly_counter = 0;

for_thresh = .1; %what is the minimum walking speed
for_length = .2; %what percentage of the trial does the animal have to be walking that speed
cue_thresh = 5;  %how much does the cue have to move to report that it was moving during the trial
rho_thresh = 0.5; %what does the average pva rho have to be to label this brain as having a "bump" during the trial

for i = 1:length(all_data)
    meta_parts  = split(all_data(i).meta,'\');
    geno        = [meta_parts{end-2:end-1}];
    fly_id      = [meta_parts{end-4:end-2}];

    if ~strcmp(fly_id,last_id)
        counter = 0;       
        fly_counter = fly_counter+1;
    end

    counter = counter+1;
    trial_num(i) = counter;
    fly_num(i) = fly_counter;
    last_id = fly_id;
    
    if (sum(all_data(i).ft.f_speed>for_thresh) > length(all_data(i).ft.f_speed)*for_length) 
        walk_idx(i) = true;
    end
    
    if contains(all_data(i).ft.pattern,'background')
        dark_idx(i) = true;
    end

    if contains(geno,'mcherry')
        mcherry_idx(i) = true;
    end

    if contains(geno,'empty')
        empty_idx(i) = true;
    end

    if contains(geno,'vglut')
        vglut_idx(i) = true;
    end

    if sum(abs(diff(unwrap(all_data(i).ft.cue)))) > cue_thresh
        cue_idx = true;
    end

    if mean(all_data(i).im.rho) > rho_thresh
        rho_idx = true;
    end
end


%% extract the integrative gain of each trial, group identical trials
hv_thresh = .4; %what is the heading variability of specific integrative gain window have to be to be counted (low variability heading traces can have any gain and it'll work, becase the fly isn't rotating)
v_thresh  = .1; %what is the maximum loss function value (circvar of circdist) to be counted as a reasonable estimate of the gain (the optimization "worked")

g = {}; %create a cell array for each trial extracting the fit gains which pass the selection criteria
for i = 1:length(all_data)
    g{i} = all_data(i).gain.g(all_data(i).gain.hv > hv_thresh & all_data(i).gain.v < v_thresh);
end
g = reshape(g,[],1);

inc_idx = walk_idx & cue_idx & rho_idx; %create an inclusion index. had to walk, cue had to be working, bump had to be detectable

group_idx = [fly_num,empty_idx,mcherry_idx,vglut_idx,dark_idx];

g = g(inc_idx); %remove all trials to exclude
group_idx = group_idx(inc_idx,:);

[unique_groups,~,ic] = unique(group_idx,'rows');
g_grouped = nan(length(unique_groups),1);
v_grouped = nan(length(unique_groups),1);

for i = 1:length(unique_groups)
    g_grouped(i) = mean(vertcat(g{ic==i}),'omitnan');
    v_grouped(i) = var(vertcat(g{ic==i}),'omitnan');
end 

group_ind   = 1+sum(unique_groups(:,2:5) .* [1,2,4,8],2);
group_labels= { 'LPsP\newlineTH-RNAi\newlineCL\newline',...
                'Empty\newlineTH-RNAi\newlineCL\newline',...
                'LPsP\newlinemCherry-RNAi\newlineCL\newline',...
                'Empty\newlinemCherry-RNAi\newlineCL\newline',...
                'LPsP\newlinevGlut-RNAi\newlineCL\newline',...
                'Empty\newlinevGlut-RNAi\newlineCL\newline',...
                'x',...
                'x',...
                'LPsP\newlineTH-RNAi\newlineDark\newline',...
                'Empty\newlineTH-RNAi\newlineDark\newline',...
                'LPsP\newlinemCherry-RNAi\newlineDark\newline',...
                'Empty\newlinemCherry-RNAi\newlineDark\newline',...
                'LPsP\newlinevGlut-RNAi\newlineDark\newline',...
                'Empty\newlinevGlut-RNAi\newlineDark\newline',...
                'x',...
                'x'};

figure(1); clf

for i = unique(group_ind)'
    subplot(2,1,1); hold on
    scatter(i*ones(sum(group_ind==i),1),g_grouped(group_ind==i),'k','filled','MarkerFaceAlpha',.1);
    errorbar(i+.1,mean(g_grouped(group_ind==i),'omitnan'),std(g_grouped(group_ind==i),'omitnan')/sqrt(sum(group_ind==i)),'or')
    group_labels{i} = [group_labels{i},sprintf('(n = %i)',sum(group_ind==i))];
    
    subplot(2,1,2); hold on
    scatter(i*ones(sum(group_ind==i),1),v_grouped(group_ind==i),'k','filled','MarkerFaceAlpha',.1);
    errorbar(i+.1,mean(v_grouped(group_ind==i),'omitnan'),std(v_grouped(group_ind==i),'omitnan')/sqrt(sum(group_ind==i)),'or')
end
subplot(2,1,1); ylabel('mean gain'); plot(xlim,.8*[1,1],':k'); xticks(1:length(group_labels)); xticklabels(group_labels)
subplot(2,1,2); ylabel('mean variance'); xticks(1:length(group_labels)); xticklabels(group_labels)
linkaxes(get(gcf,"Children"),'x')

%% extract the instantaneous gain for each trial and plot as above
vel_thresh  = .2; %minimum rotational speed of fly for a time point to be counted
bump_thresh = 10; %maximum rotational speed of bump for a time point to be counted
rho_thresh  = .2; %minimum rho of the pva for the time point to be counted
vel_max     = 10; %maximum rotational speed of the fly for a time point to be counted
lag         = 10; %the lag to apply to the fly speed to align it to the bump speed (in fictrac frames)
n_smooth    = 5;

inc_idx = walk_idx & cue_idx & rho_idx; %create an inclusion index. had to walk, cue had to be working, bump had to be detectable
group_idx = [fly_num,empty_idx,mcherry_idx,vglut_idx,dark_idx];

vel_cell = cell(length(all_data),1); %store metrics for each trial in a cell array

for i = 1:length(all_data)
    
    bump_vel = interp1(all_data(i).ft.xb,gradient(unwrap(all_data(i).im.mu)) / mean(diff(all_data(i).ft.xb)),all_data(i).ft.xf);
    fly_vel =  all_data(i).ft.r_speed;
    rho      = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);

    rho      = rho(lag+1:end);
    bump_vel = bump_vel(1+lag:end);
    fly_vel  = fly_vel(1:end-lag);

    for j = 1:n_smooth
        bump_vel = smoothdata(bump_vel,1,'gaussian',1);
        fly_vel  = smoothdata(fly_vel,1,'gaussian',300);
    end

    vel_cell{i} = [fly_vel,bump_vel,rho];
end

vel_cell = vel_cell(inc_idx); %remove all trials to exclude
group_idx = group_idx(inc_idx,:);

[unique_groups,~,ic] = unique(group_idx,'rows');
g_grouped = nan(length(unique_groups),1);
group_ind   = 1+sum(unique_groups(:,2:5) .* [1,2,4,8],2);
group_labels= { 'LPsP\newlineTH-RNAi\newlineCL\newline',...
                'Empty\newlineTH-RNAi\newlineCL\newline',...
                'LPsP\newlinemCherry-RNAi\newlineCL\newline',...
                'Empty\newlinemCherry-RNAi\newlineCL\newline',...
                'LPsP\newlinevGlut-RNAi\newlineCL\newline',...
                'Empty\newlinevGlut-RNAi\newlineCL\newline',...
                'x',...
                'x',...
                'LPsP\newlineTH-RNAi\newlineDark\newline',...
                'Empty\newlineTH-RNAi\newlineDark\newline',...
                'LPsP\newlinemCherry-RNAi\newlineDark\newline',...
                'Empty\newlinemCherry-RNAi\newlineDark\newline',...
                'LPsP\newlinevGlut-RNAi\newlineDark\newline',...
                'Empty\newlinevGlut-RNAi\newlineDark\newline',...
                'x',...
                'x'};

for i = unique(group_ind)' %initialize a iigure for each animal group
    figure(i); clf
    t(i) = tiledlayout("flow");
    title(t(i),group_labels{i})
end

for i = 1:length(unique_groups) %go through each unique id and calculate the gain

    tmp         = vertcat(vel_cell{ic==i});
    fly_vel     = tmp(:,1);
    bump_vel    = tmp(:,2);
    rho         = tmp(:,3);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh & abs(fly_vel) < vel_max;
    b = [ones(sum(idx),1),fly_vel(idx)] \ bump_vel(idx);
    g_grouped(i) = b(2);
    
    figure(group_ind(i))
    nexttile; hold on
    scatter(fly_vel(idx),bump_vel(idx),5,'filled','k','MarkerFaceAlpha',.1)
    plot([-2,2],b(1)+g_grouped(i)*[-2,2],'r')
    text(max(xlim),min(ylim),sprintf('gain: %.2f',g_grouped(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
end

%%
figure(max(group_ind)+1); clf; hold on
scatter(group_ind,g_grouped,'k','filled','MarkerFaceAlpha',.1);
for i = unique(group_ind)'
    group_labels{i} = [group_labels{i},sprintf('(n = %i)',sum(group_ind==i))];
    errorbar(i+.1,mean(g_grouped(group_ind==i)),std(g_grouped(group_ind==i))/sqrt(sum(group_ind==i)),'or');
end
scatter(group_ind,g_grouped,'k','filled','MarkerFaceAlpha',.1);
ylabel('instantaneous gain'); plot(xlim,.8*[1,1],':k'); xticks(1:length(group_labels)); xticklabels(group_labels)

%% extract the offset variability
hv_thresh = .4; %what is the heading variability of specific integrative gain window have to be to be counted (low variability heading traces can have any gain and it'll work, becase the fly isn't rotating)
v_thresh  = .1; %what is the maximum loss function value (circvar of circdist) to be counted as a reasonable estimate of the gain (the optimization "worked")

o = {}; %create a cell array for each trial extracting the fit gains which pass the selection criteria
h = {};
for i = 1:length(all_data)
    tmp_mu = all_data(i).im.mu;
    tmp_mu(all_data(i).im.rho < rho_thresh) = nan;
    tmp_mu = unwrap(tmp_mu);
    tmp = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,tmp_mu,all_data(i).ft.xf));
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    tmp = tmp(abs(all_data(i).ft.r_speed)>r_thresh);
    
    o{i} = tmp(~isnan(tmp));
    h{i} = all_data(i).ft.cue;
end
o = reshape(o,[],1);

inc_idx = walk_idx & cue_idx & rho_idx; %create an inclusion index. had to walk, cue had to be working, bump had to be detectable

group_idx = [fly_num,empty_idx,mcherry_idx,vglut_idx,dark_idx];

o = o(inc_idx); %remove all trials to exclude
h = h(inc_idx); %remove all trials to exclude
group_idx = group_idx(inc_idx,:);

[unique_groups,~,ic] = unique(group_idx,'rows');
o_grouped = nan(length(unique_groups),1);
h_grouped = nan(length(unique_groups),1);

for i = 1:length(unique_groups)
    o_grouped(i) = circ_var(vertcat(o{ic==i}));
    h_grouped(i) = circ_var(vertcat(h{ic==i}));
end

group_ind   = 1+sum(unique_groups(:,2:5) .* [1,2,4,8],2);
group_labels= { 'LPsP\newlineTH-RNAi\newlineCL\newline',...
                'Empty\newlineTH-RNAi\newlineCL\newline',...
                'LPsP\newlinemCherry-RNAi\newlineCL\newline',...
                'Empty\newlinemCherry-RNAi\newlineCL\newline',...
                'LPsP\newlinevGlut-RNAi\newlineCL\newline',...
                'Empty\newlinevGlut-RNAi\newlineCL\newline',...
                'x',...
                'x',...
                'LPsP\newlineTH-RNAi\newlineDark\newline',...
                'Empty\newlineTH-RNAi\newlineDark\newline',...
                'LPsP\newlinemCherry-RNAi\newlineDark\newline',...
                'Empty\newlinemCherry-RNAi\newlineDark\newline',...
                'LPsP\newlinevGlut-RNAi\newlineDark\newline',...
                'Empty\newlinevGlut-RNAi\newlineDark\newline',...
                'x',...
                'x'};

figure(1); clf

for i = unique(group_ind)'
    subplot(2,1,1); hold on
    scatter(i*ones(sum(group_ind==i),1),o_grouped(group_ind==i),'k','filled','MarkerFaceAlpha',.1);
    errorbar(i+.1,mean(o_grouped(group_ind==i),'omitnan'),std(o_grouped(group_ind==i),'omitnan')/sqrt(sum(group_ind==i)),'or')
    group_labels{i} = [group_labels{i},sprintf('(n = %i)',sum(group_ind==i))];
    
    subplot(2,1,2); hold on
    scatter(i*ones(sum(group_ind==i),1),h_grouped(group_ind==i),'k','filled','MarkerFaceAlpha',.1);
    errorbar(i+.1,mean(h_grouped(group_ind==i),'omitnan'),std(v_grouped(group_ind==i),'omitnan')/sqrt(sum(group_ind==i)),'or')
end
subplot(2,1,1); ylabel('offset variance'); plot(xlim,.8*[1,1],':k'); xticks(1:length(group_labels)); xticklabels(group_labels)
subplot(2,1,2); ylabel('heading variance'); xticks(1:length(group_labels)); xticklabels(group_labels)
linkaxes(get(gcf,"Children"),'x')
