%% take path length of walking bouts over the course of trials, then fit the slope of fly path length and bump path length for each walking bout

i = 25;

% compare the movement ratios for each fly
smooth_window = 60;
turn_thresh = .25;
n_frames    = 1*60; %how many seconds after a turn start to keep or analysis
max_gap_frames= .5*60; %how many frames within a walking bout can be 0 before we separate the bout in 2?
min_walking_frames= .5*60;

mov_ratio = nan(length(all_data),1);

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


mov_ratio(i) = mov_cue \ mov_mu;


%
figure(2); clf
group_idx = empty_idx +2*vglut_idx + 4*mcherry_idx;
swarmchart(group_idx(walk_idx & rho_idx & cue_idx & ~dark_idx),mov_ratio(walk_idx & rho_idx & cue_idx & ~dark_idx),...
    'filled','markerfacealpha',.5)
xlim([-.5,4.5])


%% troubleshoot fig
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