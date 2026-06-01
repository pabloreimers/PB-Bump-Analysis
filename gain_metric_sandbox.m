%% scan through masks
fig = figure(1);
i = 1;

while true
    clf
    subplot(2,2,1)
    imagesc(all_data(i).im.mask)
    title(i)

    subplot(2,1,2);
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
    hold on
    if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
    a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c,'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w','linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
    title(all_data(i).meta,'Interpreter','none')
    xlabel('time (s)')
    set(gca,'CLim',[-2,3])

    w = waitforbuttonpress;
    if w == 1 % Check if it was a keyboard press
        key = get(fig, 'CurrentKey'); % Returns 'leftarrow', 'rightarrow', etc.
    end

    if strcmp(key,'leftarrow')
        i = max(i-1,1);
    end

    if strcmp(key,'rightarrow')
        i = min(i+1,length(all_data));
    end
end

%% take path length of walking bouts over the course of trials, then fit the slope of fly path length and bump path length for each walking bout

i = 25;

% compare the movement ratios for each fly
smooth_window = 60;
turn_thresh = .25;
n_frames    = 1*60; %how many seconds after a turn start to keep or analysis
max_gap_frames= .5*60; %how many frames within a walking bout can be 0 before we separate the bout in 2?
min_walking_frames= .5*60;

mu_cell = cell(length(all_data),1);
cue_cell= cell(length(all_data),1);

for i = 1:length(all_data)
    dt = median(diff(all_data(i).ft.xf));

    mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
    cue= unwrap(all_data(i).ft.cue);

    mu_smooth = smoothdata(mu,"gaussian",smooth_window);
    cue_smooth= smoothdata(-cue,"gaussian",smooth_window);
    fly_speed = [abs(diff(cue_smooth));nan]/dt;

    is_walking =  fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames + 1, 1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false; is_walking; false]);
    boutStarts = find(d == 1);
    boutEnds   = find(d == -1) - 1;

    % turning_idx = [abs(diff(cue_smooth));nan]/dt > turn_thresh;
    % turn_starts = find(diff(turning_idx) == 1);
    % turn_inds   = turn_starts + (0:n_frames);
    % turn_idx    = false(length(all_data(i).ft.xf));
    % turn_idx(turn_inds) = true;

    mov_mu = nan(length(boutStarts),1);
    mov_cue = nan(length(boutStarts),1);

    for b = 1:length(boutStarts)

        mov_mu(b)  = sum(abs(diff(mu_smooth(boutStarts(b):boutEnds(b)))),'omitnan');
        mov_cue(b) = sum(abs(diff(cue_smooth(boutStarts(b):boutEnds(b)))),'omitnan');
    end

    mu_cell{i} = mov_mu;
    cue_cell{i} = mov_cue;

    % figure(2); clf; 
    % subplot(2,1,1); hold on
    % 
    % plot(mu_smooth)
    % plot(cue_smooth)
    % plot(turning_idx)
    % 
    % subplot(2,1,2)
    % scatter(mov_cue,mov_mu)
    % title(mov_ratio(i))
end

i = 25
% troubleshoot fig
figure(3); clf

a1 = subplot(3,1,1);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c,'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
idx = round(all_data(i).ft.cue,4) == -.2945;

a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w','linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta,'Interpreter','none')
xlabel('time (s)')
subplot(3,1,2); hold on
plot(all_data(i).ft.xf,fly_speed)
plot(all_data(i).ft.xf,is_walking)
scatter(all_data(i).ft.xf(boutStarts),0,'r*')
scatter(all_data(i).ft.xf(boutEnds),0,'b*')
linkaxes(get(gcf,'Children'),'x')

subplot(3,2,5);
scatter(mov_cue,mov_mu)
title(mov_ratio(i))


%%
% compare the movement ratios for each fly
lag = 7;
smooth_window = 60;
turn_thresh = .25;
max_gap_frames= .5*60; %how many frames within a walking bout can be 0 before we separate the bout in 2?
min_walking_frames= .5*60;
rho_thresh = .2;

mu_cell = cell(length(all_data),1);
cue_cell= cell(length(all_data),1);

for i = 1:length(all_data)
    dt = median(diff(all_data(i).ft.xf));

    mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
    cue= smoothdata(unwrap(all_data(i).ft.cue),'gaussian',30);
    rho= interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);

    mu_lag = mu(1:end-lag);
    cue_lag = cue(1+lag:end);
    rho_lag = rho(1:end-lag);

    good = rho_lag > rho_thresh;

    mu_interp = nan(length(mu_lag),1);
    mu_interp(good) = mu_lag(good);
    mu_interp = fillmissing(mu_interp,'linear');

    % Find runs of bad samples
    d = diff([0; ~good; 0]);
    start_idx = find(d == 1);
    end_idx = find(d == -1) - 1;

    % Remove interpolated regions that are too long
    for k = 1:length(start_idx)

        gap_len = end_idx(k) - start_idx(k) + 1;

        if gap_len > max_gap_frames
            mu_interp(start_idx(k):end_idx(k)) = nan;
        end
    end

    mu_smooth = smoothdata(mu_interp,"gaussian",smooth_window);
    cue_smooth= smoothdata(-cue_lag,"gaussian",smooth_window);

    

    fly_speed = [abs(diff(cue_smooth));nan]/dt;

    is_walking = fly_speed > turn_thresh;
    is_walking = imclose(is_walking,ones(max_gap_frames + 1, 1));
    is_walking = bwareaopen(is_walking,min_walking_frames);

    d = diff([false; is_walking; false]);
    boutStarts = find(d == 1);
    boutEnds   = find(d == -1) - 1;

    mov_mu = nan(length(boutStarts),1);
    mov_cue = nan(length(boutStarts),1);

    for b = 1:length(boutStarts)

        mov_mu(b)  = sum(abs(diff(mu_smooth(boutStarts(b):boutEnds(b)))),'omitnan');
        mov_cue(b) = sum(abs(diff(cue_smooth(boutStarts(b):boutEnds(b)))),'omitnan');
    end

    mu_cell{i} = mov_mu;
    cue_cell{i} = mov_cue;

end


inc_idx = walk_idx & cue_idx & rho_idx; %create an inclusion index. had to walk, cue had to be working, bump had to be detectable

group_idx = [fly_num,empty_idx,vglut_idx,mcherry_idx,dark_idx];

cue_cell = cue_cell(inc_idx); %remove all trials to exclude
mu_cell  = mu_cell(inc_idx);
group_idx = group_idx(inc_idx,:);

[unique_groups,~,ic] = unique(group_idx,'rows');

mov_ratio = nan(length(unique_groups),1);
mov_corr  = nan(length(unique_groups),1);

for i = 1:length(unique_groups)
    if size(vertcat(cue_cell{ic==i}),1) > 20
        mov_ratio(i) = vertcat(cue_cell{ic==i}) \ vertcat(mu_cell{ic==i});
        mov_corr(i)  = corr(vertcat(cue_cell{ic==i}), vertcat(mu_cell{ic==i}))
    end
end 

group_ind   = 1+sum(unique_groups(:,2:5) .* [1,2,4,8],2);
group_labels= { 'LPsP\newlineTH-RNAi\newlineCL\newline',...
    'Empty\newlineTH-RNAi\newlineCL\newline',...
    'LPsP\newlinevGlut-RNAi\newlineCL\newline',...
    'Empty\newlinevGlut-RNAi\newlineCL\newline',...
    'LPsP\newlinemCherry-RNAi\newlineCL\newline',...
    'Empty\newlinemCherry-RNAi\newlineCL\newline',...
    'x',...
    'x',...
    'LPsP\newlineTH-RNAi\newlineDark\newline',...
    'Empty\newlineTH-RNAi\newlineDark\newline',...
    'LPsP\newlinevGut-RNAi\newlineDark\newline',...
    'Empty\newlinevGlut-RNAi\newlineDark\newline',...
    'LPsP\newlinemCherry-RNAi\newlineDark\newline',...
    'Empty\newlinemCherry-RNAi\newlineDark\newline',...
    'x',...
    'x'};

figure(1); clf

for i = unique(group_ind)'
    hold on
    scatter(i*ones(sum(group_ind==i),1),mov_ratio(group_ind==i),'k','filled','MarkerFaceAlpha',.1);
    errorbar(i+.1,mean(mov_ratio(group_ind==i),'omitnan'),std(mov_ratio(group_ind==i),'omitnan')/sqrt(sum(group_ind==i)),'or')
    group_labels{i} = [group_labels{i},sprintf('(n = %i)',sum(~isnan(mov_ratio(group_ind==i))))];

end
ylabel('path length ratio'); plot(xlim,[1,1],':k'); xticks(1:length(group_labels)); xticklabels(group_labels)
linkaxes(get(gcf,"Children"),'x')


%% compare velocity scatter plots, look at gain, lag, corrcoeff

i = 4;

vel_thresh = .2;
bump_max = 10;
rho_thresh = .2;
vel_max = 10;
lag = 7;

g_mat = nan(length(all_data),30);
corr_mat = nan(length(all_data),30);



for i = 1:length(all_data)

    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = gradient(-all_data(i).ft.cue)/fr; %all_data(i).ft.r_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);


    for lag = 1:30
    
    fly_lag  = fly_vel(1:end-lag);
    bump_lag = bump_vel(lag+1:end);
    rho_lag      = rho(lag+1:end);

    idx = abs(fly_lag) > vel_thresh & abs(bump_lag) < bump_max & rho_lag > rho_thresh & abs(fly_lag) < vel_max;
    
    
    if sum(idx) > 60
    g_mat(i,lag) = fly_lag(idx) \ bump_lag(idx);
    corr_mat(i,lag) = corr(fly_lag(idx),bump_lag(idx)); 
    end
 
    end

end

%%
group_idx = empty_idx + 2*vglut_idx + 4*mcherry_idx;
inc_idx = walk_idx & cue_idx & rho_idx & ~dark_idx;

figure(2); clf
subplot(2,2,1)
imagesc(g_mat(inc_idx,:)); colorbar
set(gca,'Clim',[0,2])

subplot(2,2,3)
swarmchart(group_idx(inc_idx),g_mat(inc_idx,7))

subplot(2,2,2)
imagesc(corr_mat(inc_idx,:)); colorbar

subplot(2,2,4)
swarmchart(group_idx(inc_idx),corr_mat(inc_idx,7))


%%
figure(4);clf
cc_mat = nan(length(all_data),30);

for i = 1:length(all_data)
    mu  = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
    cue = -unwrap(all_data(i).ft.cue);
    
    for lag = 1:30
        mu_lag = mu(1:end-lag);
        cue_lag= cue(1+lag:end);
        
        idx = ~isnan(mu_lag) & ~isnan(cue_lag);
        cc_mat(i, lag) = circ_corrcc(mu_lag(idx),cue_lag(idx));
    end
end

imagesc(cc_mat(inc_idx,:))
