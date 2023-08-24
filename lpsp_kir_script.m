%% set experiment params
transfer_flag = false;
local_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\analysis\'; %define the local directory to clear up
remote_dir           = '\\files.med.harvard.edu\Neurobio\wilsonlab\pablo\lpsp_kir\analysis'; %define the remote directory to make sure things are backed up to

f0_pct              = 5;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 10;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 1;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = 0.1;
eb_flag             = true;
%% 
if transfer_flag

date_folders        = cellstr(ls(remote_dir)); %store all the subfolders, which include the dates
date_folders(1:2) = []; %clean up this variable

for i = date_folders'
    tmp_remote= [remote_dir,'\',i{1}];

    remote_flies= cellstr(ls(tmp_remote));
    remote_flies(1:2,:) = [];

    for f = remote_flies'
    
    tmp_remote= [remote_dir,'\',i{1},'\',f{1}];

    remote_trials= cellstr(ls(tmp_remote));
    remote_trials(1:2,:) = [];
    
    for j = remote_trials'
        
        tmp_remote= [remote_dir,'\',i{1},'\',f{1},'\',j{1}];
        tmp_local = [local_dir,'\',i{1},'\',f{1},'\',j{1}];
        fprintf('current folder: %s',[i{1},'\',f{1},'\',j{1}])
        if ~(contains(j{1},'cl') | contains(j{1},'dark'))
            fprintf('\n')
            continue
        end
        mkdir(tmp_local)
        try
        copyfile([tmp_remote,'\registration_001\imagingData*'],tmp_local)
        copyfile([tmp_remote,'\*ficTracData_DAQ*'],tmp_local)
        catch
            fprintf('***ERROR with imaging, ***')
        end
        try
        copyfile([tmp_remote,'\mask*'],tmp_local)
        catch
            fprintf('*** No Mask ***')
        end
        fprintf('\n')
    end
    end
end
end

%% analyze each trial and store appropriate values
dff_tot = {}; %preallocate dff cell
dff_peak= {};
dff_cluster = {};
f_cluster = {};
mu      = {};
rho     = {};
f_speed = {};
r_speed = {};
intHD   = {};
cue     = {};
r_vel   = {};
xf      = {};
cl_idx  = {};
exp_idx = {};
fly_num = {};
exp_date= {};
trial   = {};

i = 0;

dates = cellstr(ls(local_dir));
dates(1:2) = [];
for d = dates'
    d = d{1};
    flies = cellstr(ls([local_dir,'\',d]));
    flies(1:2) = [];
    for f = flies'
        f = f{1};
        trials = cellstr(ls([local_dir,'\',d,'\',f]));
        trials(1:2) = [];
        tmp = cellfun(@(x)(contains(x,'cl')),trials) | cellfun(@(x)(contains(x,'dark')),trials);
        trials = trials(tmp);
        for t = trials'
            t = t{1};
            fprintf('date: %s  fly: %s  trial: %s',d,f,t)
            
            try
            tmp_dir = [local_dir,'\',d,'\',f,'\',t,'\'];
            tmp = ls([tmp_dir,'mask*']);
            load([tmp_dir,tmp])
            tmp = ls([tmp_dir,'imagingData*']);
            load([tmp_dir,tmp])
            tmp = ls([tmp_dir,'*ficTracData*']);
            load([tmp_dir,tmp])
            catch
                fprintf(' *** ERROR ***\n')
                continue
            end
            i = i+1;
            fly_num{i} = f(end);
            exp_date{i} = d;
            trial{i} = t;
            cl_idx{i} = contains(t,'cl');
            exp_idx{i} = contains(t,'LPsP');
            [f_speed{i},r_speed{i},intHD{i},cue{i},r_vel{i}] = ft_calc(ftData_DAQ,n_smooth,f_smooth);
            xf{i}  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed{i})]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
            if isduration(xf{i})
                xf{i} = seconds(xf{i});
            end
            [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i},f_cluster{i}] = bump_calc_eb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
            fprintf('\n')
        end
    end
end

%% store everything in a table for easy access
full_data = table(dff_tot, dff_peak, dff_cluster, f_cluster, mu, rho, f_speed,r_speed,intHD,...
    cue, r_vel, xf,cl_idx,exp_idx,fly_num,exp_date,trial);

%% show the raw traces sorted by percentage of time a vector is shown
rows = 6;
cols = 2;

[~,idx] = sort(cellfun(@(x) sum(x > rho_thresh) / length(x),full_data.rho),'ascend');

f2 = 0;
counter = 13;
for i = idx
    counter = counter+1;

    if counter > 12
        f2 = f2+1;
        figure(f2); clf
        tiledlayout(rows,cols)
        counter = 1;
    end

    nexttile
    colormap(parula)
    imagesc([1,length(full_data.mu{i})],[-pi,pi],full_data.dff_cluster{i})
    hold on
    a = plot(full_data.mu{i},'w');
    a.YData(abs(diff(a.YData))>pi) = nan;
    a.YData(full_data.rho{i} < rho_thresh) = nan;
    xticks([])
    ylabel(['trial num: ', num2str(i)])
    xlabel(sum(full_data.rho{i} > rho_thresh)/length(full_data.rho{i}))
    title(full_data.trial{i}(1:8))
end

%% show the raw traces for each group
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_names = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

rows = 6;
cols = 2;

f2 = 0;
for f = 1:length(unique(group_idx))
counter = 13;
for i = find(group_idx == f)
    counter = counter+1;

    if counter > 12
        f2 = f2+1;
        figure(f2); clf
        set(gcf,'Name',group_names{f})
        tiledlayout(rows,cols)
        counter = 1;
    end

    nexttile
    colormap(parula)
    imagesc([1,length(full_data.mu{i})],[-pi,pi],full_data.dff_cluster{i})
    hold on
    a = plot(-full_data.cue{i},'m','Linewidth',2);
    a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(full_data.mu{i},'w');
    a.YData(abs(diff(a.YData))>pi) = nan;
    a.YData(full_data.rho{i} < .2) = nan;
    xticks([])
    %xlabel({['offset var: ',num2str(offset_var(i,1))],['circ corr: ' num2str(circ_corr(i,1))],['unwrap corr: ', num2str(unwrap_corr(i,1))]})
    title(full_data.trial{i}(1:8))
    ylabel(['trial num: ',num2str(i)])
end
end

%% show mean image for each group
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_names = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

rows = 3;
cols = 4;


f2 = 0;
counter = 13;

dates = cellstr(ls(local_dir));
dates(1:2) = [];
for d = dates'
    d = d{1};
    flies = cellstr(ls([local_dir,'\',d]));
    flies(1:2) = [];
    for f = flies'
        f = f{1};
        trials = cellstr(ls([local_dir,'\',d,'\',f]));
        trials(1:2) = [];
        tmp = cellfun(@(x)(contains(x,'cl')),trials) | cellfun(@(x)(contains(x,'dark')),trials);
        trials = trials(tmp);
        for t = trials'
            t = t{1};
            fprintf('date: %s  fly: %s  trial: %s',d,f,t)
            tmp_dir = [local_dir,'\',d,'\',f,'\',t,'\'];
            try
            tmp = ls([tmp_dir,'imagingData*']);
            load([tmp_dir,tmp])
            catch
                fprintf(' *** ERROR ***\n')
                continue
            end


            counter = counter+1;
        
            if counter > 12
                f2 = f2+1;
                figure(f2); clf
                tiledlayout(rows,cols)
                counter = 1;
            end
        
            nexttile
            imagesc(mean(sum(regProduct,3),4))
            colormap(bone)
            axis equal tight
            xticks([])
            yticks([])
            title([d,' ', f])
            fprintf('\n')
            drawnow
        end
    end
end


%% show the vel_corrs for each group
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_names = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

vel_corr = nan(size(exp_idx,2),20);

for j =1:20
lag = j;
for i = 1:length(vel_corr)
    mu_tmp = mu{i}(1+lag:end);
    rho_tmp = rho{i}(1+lag:end);
    vel_tmp = r_vel{i}(1:end-lag);
    
    mu_tmp = smoothdata(unwrap(mu_tmp),1,'gaussian',f_smooth);
    mu_tmp = mod(mu_tmp+pi,2*pi)-pi;

    tmp = ~isnan(mu_tmp) & (rho_tmp > rho_thresh);
    bump_vel = [diff(mu_tmp(tmp));0]./fr;
    fly_vel  = vel_tmp(tmp);
    tmp = abs(fly_vel) > vel_min & abs(bump_vel) < 10;
    vel_corr(i,j) = corr(bump_vel(tmp),fly_vel(tmp));
end
end


%% plot the results
for f = 1:length(unique(group_idx))
figure(f); clf
set(gcf,'Name',group_names{f})
counter = 0;
rows = floor(sqrt(sum(group_idx==f)));
cols = ceil(sum(group_idx==f)/rows);
for i = find(group_idx == f)
    counter = counter+1;
    
    %[~,lag] = max(vel_corr(i,:),[],2);
    lag = lag_idx(group_idx(i))
    mu_tmp = mu{i}(1+lag:end);
    rho_tmp = rho{i}(1+lag:end);
    vel_tmp = r_vel{i}(1:end-lag);

    mu_vel = [diff(mu_tmp);0]./fr;
    
    idx = abs(vel_tmp) > vel_min & ~isnan(mu_vel) & abs(mu_vel) < 10;
    
    g = [zeros(length(vel_tmp(idx)),1), vel_tmp(idx)] \ mu_vel(idx);
    

    subplot(rows,cols,counter)
    scatter(vel_tmp(idx),mu_vel(idx),2,'filled','Markerfacealpha',0.1,'markerfacecolor',scatter_color)
    hold on
    axis equal
    x = xlim;
    y = ylim;
    plot(x,[0,0],':','color',scatter_color)
    plot([0,0],y,':','color',scatter_color)
    plot(x,x*g(2)+g(1),'r' )
    xlim(x); ylim(y)
    text(x(2),y(2),{['corr: ' num2str(vel_corr(i,lag))],['gain: ', num2str(g(2))]},'HorizontalAlignment','right','VerticalAlignment','top','color',scatter_color)
    if counter>(rows-1)*cols
        xlabel('fly vel (rad/s)')
    end
    if mod(counter,cols) == 1
        ylabel('bump vel (rad/s)')
    end
end

if dark_flag
    set(gcf,'color','none')
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
    try
        set(tmp(i),'color','none')
        set(tmp(i),'xcolor','w')
        set(tmp(i),'ycolor','w')
    end
    end
   
end
end




%% calculate offset variability across groups (time consuming)
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
lag_size = 5; %in frames
max_lag = 1000; %in ms

lag_vec  = 0:lag_size*fr*1000:max_lag;
offset_var = nan(size(exp_idx,2),length(lag_vec));
circ_corr  = nan(size(exp_idx,2),length(lag_vec));
unwrap_corr  = nan(size(exp_idx,2),length(lag_vec));
win_size = 60; %in seconds


for j = 1:length(lag_vec)
lag = (j-1)*lag_size;
for i = 1:length(offset_var)
    mu_tmp = full_data.mu{i}(lag+1:end);
    rho_tmp = full_data.rho{i}(lag+1:end);
    cue_tmp = -full_data.cue{i}(1:end-lag);
    x_tmp  = full_data.xf{i}(1:end-lag);
    r_speed = full_data.r_speed{i}(1:end-lag);

    ov = nan(ceil(x_tmp(end)-win_size),1);
    cc = nan(ceil(x_tmp(end)-win_size),1);
    uc = nan(ceil(x_tmp(end)-win_size),1);
    

    for k = 1:ceil(x_tmp(end))-win_size
        tmp = ~isnan(mu_tmp) & r_speed > 0.2 & (rho_tmp > rho_thresh) & [abs(diff(cue_tmp));0] > 0 & x_tmp' > k & x_tmp' < k+win_size;
        if sum(tmp) < 0.2*sum(x_tmp' > k & x_tmp' < k+win_size)
            continue
        end
        
        tmp1 = mu_tmp;
        tmp1(~tmp) = nan;
        tmp1 = mod(smoothdata(unwrap(tmp1),1,'gaussian',f_smooth)+pi,2*pi)-pi;
        
        tmp3 = unwrap(tmp1);
        tmp4 = unwrap(cue_tmp);
        uc(k)= corr(tmp3(~isnan(tmp1)),tmp4(~isnan(tmp1)));
        
        tmp2 = cue_tmp(~isnan(tmp1));
        tmp1 = tmp1(~isnan(tmp1));

        cc(k) = circ_corrcc(tmp1,tmp2);
        ov(k) = circ_var(circ_dist(tmp1,tmp2));
        fprintf("j: %.2f i: %.2f k: %.2f\n",j/length(lag_vec),i/length(offset_var),k/ceil(x_tmp(end)))
    end
    offset_var(i,j) = mean(ov,'omitnan');
    circ_corr(i,j) = mean(cc,'omitnan');
    unwrap_corr(i,j) = mean(uc,'omitnan');
end
end

%% plot the results
c = [  1,  .25,  0;...
       0,  .25,  1;...
       1, .75,.75;...
     .5, .5, 1];

keep_idx = cellfun(@(x)(sum(x>rho_thresh)/length(x) > 0.5),full_data.rho);

lag_idx = nan(4,1);
[~,lag_idx(1)] = min(mean(offset_var(group_idx == 1,:),1));
[~,lag_idx(2)] = min(mean(offset_var(group_idx == 2,:),1));
[~,lag_idx(3)] = min(mean(offset_var(group_idx == 3,:),1));
[~,lag_idx(4)] = min(mean(offset_var(group_idx == 4,:),1));

p1 = ranksum(offset_var(group_idx==1 & keep_idx,lag_idx(1)),offset_var(group_idx==2 & keep_idx,lag_idx(2)),'Tail','left'); %test that the correlation for empty is greater than for LPsP
p2 = ranksum(offset_var(group_idx==3 & keep_idx,lag_idx(3)),offset_var(group_idx==4 & keep_idx,lag_idx(4)),'Tail','left'); %test that the correlation for empty is greater than for LPsP

o_empty = offset_var(group_idx==1 & keep_idx,lag_idx(1));
o_lpsp  = offset_var(group_idx==2 & keep_idx,lag_idx(2));
idx_empty = randi(length(g_empty),length(g_empty),N);
idx_lpsp  = randi(length(g_lpsp),length(g_lpsp),N);
p1 = sum(mean(o_empty(idx_empty),1) - mean(o_lpsp(idx_lpsp),1) < 0)/ N;
g_empty = offset_var(group_idx==3 & keep_idx,lag_idx(3));
g_lpsp  = offset_var(group_idx==4 & keep_idx,lag_idx(4));
idx_empty = randi(length(g_empty),length(g_empty),N);
idx_lpsp  = randi(length(o_lpsp),length(o_lpsp),N);
p2 = sum(mean(g_empty(idx_empty),1) - mean(g_lpsp(idx_lpsp),1) < 0)/ N;

figure(1); clf

for i = 1:4
    subplot(2,3,1); hold on
    m = mean(offset_var(group_idx==i & keep_idx,lag_idx(i)),1);
    s = std(offset_var(group_idx==i  & keep_idx,lag_idx(i)),1)/sqrt(sum(group_idx==i & keep_idx));
    scatter(group_idx(group_idx==i  & keep_idx),offset_var(group_idx==i  & keep_idx,lag_idx(i)),100,'filled',scatter_color,'MarkerFaceAlpha',.5)
    errorbar(i+.2,m,s,'Color',c(i,:),'LineWidth',2)

    subplot(2,3,4); hold on
    m = mean(offset_var(group_idx==i & keep_idx,:),1);
    s = std(offset_var(group_idx==i  & keep_idx,:),1)/sqrt(sum(group_idx==i & keep_idx));
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);
    
    subplot(2,3,2); hold on
    m = mean(circ_corr(group_idx==i & keep_idx,lag_idx(i)),1);
    s = std(circ_corr(group_idx==i  & keep_idx,lag_idx(i)),1)/sqrt(sum(group_idx==i & keep_idx));
    scatter(group_idx(group_idx==i  & keep_idx),circ_corr(group_idx==i  & keep_idx,lag_idx(i)),100,'filled',scatter_color,'MarkerFaceAlpha',.5)
    errorbar(i+.2,m,s,'Color',c(i,:),'LineWidth',2)

    subplot(2,3,5); hold on
    m = mean(circ_corr(group_idx==i & keep_idx,:),1);
    s = std(circ_corr(group_idx==i  & keep_idx,:),1)/sqrt(sum(group_idx==i & keep_idx));
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);
    
    subplot(2,3,3); hold on
    m = mean(unwrap_corr(group_idx==i & keep_idx,lag_idx(i)),1);
    s = std(unwrap_corr(group_idx==i  & keep_idx,lag_idx(i)),1)/sqrt(sum(group_idx==i & keep_idx));
    scatter(group_idx(group_idx==i  & keep_idx),unwrap_corr(group_idx==i  & keep_idx,lag_idx(i)),100,'filled',scatter_color,'MarkerFaceAlpha',.5)
    errorbar(i+.2,m,s,'Color',c(i,:),'LineWidth',2)

    subplot(2,3,6); hold on
    m = mean(unwrap_corr(group_idx==i & keep_idx,:),1);
    s = std(unwrap_corr(group_idx==i  & keep_idx,:),1)/sqrt(sum(group_idx==i & keep_idx));
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);


end
subplot(2,3,1)
title('Offset Variability','color',scatter_color)
%ylabel({'Offset Variability',sprintf('%d s rolling average',win_size)})
xticks([1,2,3,4]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([.5,4.5])

subplot(2,3,2)
title('Circular Correlation','color',scatter_color)
%ylabel({'Circular Correlation',sprintf('%d s rolling average',win_size)})
xticks([1,2,3,4]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([.5,4.5])

subplot(2,3,3)
title('Unwrapped Correlation','color',scatter_color)
%ylabel({'Unwrapped Correlation',sprintf('%d s rolling average',win_size)})
xticks([1,2,3,4]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([.5,4.5])

% tmp = ylim;
% plot([1,2],tmp(2)*[1,1].*.9,'k','Linewidth',3)
% text(1.5,tmp(2).*.9,sprintf('p = %.3f',p1),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')
% plot([3,4],tmp(2)*[1,1].*.9,'k','Linewidth',3)
% text(3.5,tmp(2).*.9,sprintf('p = %.3f',p2),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')

% % set(gcf,'color','none')
% % set(gca,'color','none','xcolor','w','ycolor','w')
fontsize(gcf,30,'pixels')

subplot(2,3,4)
ylabel(sprintf('%d s rolling average',win_size))
xlabel('lag (ms)')
axis tight

subplot(2,3,5)
xlabel('lag (ms)')
axis tight

subplot(2,3,6)
xlabel('lag (ms)')
axis tight

if dark_flag
    set(gcf,'color','none')
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
    try
        set(tmp(i),'color','none')
        set(tmp(i),'xcolor','w')
        set(tmp(i),'ycolor','w')
    end
    end
    legend(group_labels,'Autoupdate','off','Location','Northeast','textcolor','w', 'edgecolor','w')
end
%% compare pearson correlation of velocities
dark_flag = true;
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_labels = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};
c = [  1,  .25,  0;...
       0,  .25,  1;...
       1, .75,.75;...
     .5, .5, 1];

fr = mean(diff(full_data.xf{1}));
vel_corr = nan(size(exp_idx,2),20);
vel_slope = nan(size(exp_idx,2),20);
vel_bias = nan(size(exp_idx,2),20);
lag_vec  = (0:29)*fr*1000;
if dark_flag
    scatter_color = 'w';
else
    scatter_color = 'k';
end

for j =1:30
lag = j;
for i = 1:length(vel_corr)
    mu_tmp = full_data.mu{i}(lag:end);
    rho_tmp = full_data.rho{i}(lag:end);
    fly_vel = full_data.r_vel{i}(1:end-lag+1);

    bump_vel = [diff(mu_tmp);0]./fr;
    
    tmp = ~isnan(bump_vel) & (rho_tmp > rho_thresh) & abs(fly_vel) > vel_min & abs(bump_vel) < vel_thresh & [abs(diff(full_data.cue{i}(1:end-lag+1)));0] > 0;
    
    vel_corr(i,j) = corr(bump_vel(tmp),fly_vel(tmp));
    b = [ones(sum(tmp),1) bump_vel(tmp)] \ fly_vel(tmp);
    vel_slope(i,j) = b(2);
    vel_bias(i,j) = b(1);
end
end

figure(13); clf
subplot(2,2,1)
hold on
[~,opt_lag] = max(vel_corr,[],2);
% for i = 1:length(colormap)
%     scatter(nan,1,[],c(i,:),'filled')
% end
%legend(group_labels,'Autoupdate','off')
scatter(1:length(opt_lag),lag_vec(opt_lag),100,c(group_idx,:),'filled');
ylabel('optimal lag (ms)')
xlabel('trial number')


keep_idx = cellfun(@(x)(sum(x>rho_thresh)/length(x) > 0.15),full_data.rho);

lag_idx = nan(4,1);
[~,lag_idx(1)] = max(mean(vel_corr(group_idx == 1 & keep_idx,:),1));
[~,lag_idx(2)] = max(mean(vel_corr(group_idx == 2 & keep_idx,:),1));
[~,lag_idx(3)] = max(mean(vel_corr(group_idx == 3 & keep_idx,:),1));
[~,lag_idx(4)] = max(mean(vel_corr(group_idx == 4 & keep_idx,:),1));

p1 = ranksum(vel_corr(group_idx==1 & keep_idx,lag_idx(1)),vel_corr(group_idx==2 & keep_idx,lag_idx(2)),'Tail','right'); %test that the correlation for empty is greater than for LPsP
p2 = ranksum(vel_corr(group_idx==3 & keep_idx,lag_idx(3)),vel_corr(group_idx==4 & keep_idx,lag_idx(4)),'Tail','right'); %test that the correlation for empty is greater than for LPsP


subplot(2,2,2); cla
hold on
plot([1,2],tmp(2)*[1,1].*.9,scatter_color,'Linewidth',3)
text(1.5,tmp(2).*.9,sprintf('p = %.3f',p1),'color',scatter_color,'HorizontalAlignment','center','VerticalAlignment','bottom')
plot([3,4],tmp(2)*[1,1].*.9,scatter_color,'Linewidth',3)
text(3.5,tmp(2).*.9,sprintf('p = %.3f',p2),'color',scatter_color,'HorizontalAlignment','center','VerticalAlignment','bottom')
plot([.5,4.5],[0,0],':','color',scatter_color)
ylabel({'Pearson CorrCoeff:', 'Bump and Fly vel'})
xticks([1,2,3,4]); xticklabels(group_labels)
xlim([.5,4.5])
fontsize(gcf,30,'pixels')

bump_pct = cell2mat(cellfun(@(x)(sum(x>(0:.05:1))/length(x)),full_data.rho,'UniformOutput',false)');
[~,p1] = kstest2(mean(bump_pct(group_idx==1 & keep_idx,:),1),mean(bump_pct(group_idx==2 & keep_idx,:),1)); %test that the correlation for empty is greater than for LPsP
[~,p2] = kstest2(mean(bump_pct(group_idx==3 & keep_idx,:),1),mean(bump_pct(group_idx==4 & keep_idx,:),1)); %test that the correlation for empty is greater than for LPsP


subplot(2,2,3); cla; hold on
for i = 1:4
    subplot(2,2,3)
    m = mean(vel_corr(group_idx==i&keep_idx,:),1);
    s = std(vel_corr(group_idx==i&keep_idx,:),1)/sqrt(sum(group_idx==i&keep_idx));
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);

    subplot(2,2,2)
    scatter(group_idx(keep_idx&group_idx==i),vel_corr(keep_idx&group_idx==i,lag_idx(i)),100,scatter_color,'filled','MarkerFaceAlpha',.5)
    errorbar(i+.2,m(lag_idx(i)),s(lag_idx(i)),'Color',c(i,:),'Linewidth',2)
    
    subplot(2,2,4); hold on
    m = mean(bump_pct(group_idx==i & keep_idx,:),1);
    s = std(bump_pct(group_idx==i & keep_idx,:),1)/sqrt(sum(group_idx==i));
    patch([0:.05:1,fliplr(0:.05:1)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(0:.05:1,m,'Color',c(i,:),'Linewidth',2);
    plot(0:.05:1,bump_pct(group_idx==i & keep_idx,:)','Color',[c(i,:),.1],'Linewidth',1,'HandleVisibility','off');
end
subplot(2,2,3);
xlabel('lag (ms)')
ylabel({'Mean Pearson CorrCoeff';'(+- s.e.m.)'})
fontsize(gcf,30,'pixels')
axis tight
subplot(2,2,2)
ylim([min(ylim),1])
subplot(2,2,4)
legend(group_labels,'Autoupdate','off','Location','Northeast')
xlabel('\rho_{threshold}')
ylabel('\rho > \rho_{threshold} (%)')

if dark_flag
    set(gcf,'color','none')
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
    try
        set(tmp(i),'color','none')
        set(tmp(i),'xcolor','w')
        set(tmp(i),'ycolor','w')
    end
    end
    legend(group_labels,'Autoupdate','off','Location','Northeast','textcolor','w', 'edgecolor','w')
end

%% compare gains
dark_flag = true;
N=1e4;
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_labels = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};
c = [  1,  .25,  0;...
       0,  .25,  1;...
       1, .75,.75;...
     .5, .5, 1];

fr = mean(diff(full_data.xf{1}));
vel_corr = nan(size(exp_idx,2),20);
vel_slope = nan(size(exp_idx,2),20);
vel_bias = nan(size(exp_idx,2),20);
lag_vec  = (0:29)*fr*1000;

for j =1:30
lag = j;
for i = 1:length(vel_corr)
    mu_tmp = full_data.mu{i}(lag:end);
    rho_tmp = full_data.rho{i}(lag:end);
    fly_vel = full_data.r_vel{i}(1:end-lag+1);

    bump_vel = [diff(mu_tmp);0]./fr;
    
    tmp = ~isnan(bump_vel) & (rho_tmp > rho_thresh) & abs(fly_vel) > vel_min & abs(bump_vel) < 10 & [abs(diff(full_data.cue{i}(1:end-lag+1)));0] > 0;
    
    vel_corr(i,j) = corr(bump_vel(tmp),fly_vel(tmp));
    b = [zeros(sum(tmp),1),fly_vel(tmp)] \ bump_vel(tmp);
    vel_slope(i,j) = b(2);
    vel_bias(i,j) = b(1);
end
end

figure(14); clf

keep_idx = cellfun(@(x)(sum(x>rho_thresh)/length(x) > 0.4),full_data.rho);

lag_idx = nan(4,1);
[~,lag_idx(1)] = max(mean(vel_corr(group_idx == 1 & keep_idx,:),1));
[~,lag_idx(2)] = max(mean(vel_corr(group_idx == 2 & keep_idx,:),1));
[~,lag_idx(3)] = max(mean(vel_corr(group_idx == 3 & keep_idx,:),1));
[~,lag_idx(4)] = max(mean(vel_corr(group_idx == 4 & keep_idx,:),1));

p1 = ranksum(vel_slope(group_idx==1 & keep_idx,lag_idx(1)),vel_slope(group_idx==2 & keep_idx,lag_idx(2)),'Tail','right'); %test that the correlation for empty is greater than for LPsP
p2 = ranksum(vel_slope(group_idx==3 & keep_idx,lag_idx(3)),vel_slope(group_idx==4 & keep_idx,lag_idx(4)),'Tail','right'); %test that the correlation for empty is greater than for LPsP

g_empty = vel_slope(group_idx==1 & keep_idx,lag_idx(1));
g_lpsp  = vel_slope(group_idx==2 & keep_idx,lag_idx(2));
idx_empty = randi(length(g_empty),length(g_empty),N);
idx_lpsp  = randi(length(g_lpsp),length(g_lpsp),N);
p1 = sum(mean(g_empty(idx_empty),1) - mean(g_lpsp(idx_lpsp),1) < 0)/ N;
g_empty = vel_slope(group_idx==3 & keep_idx,lag_idx(3));
g_lpsp  = vel_slope(group_idx==4 & keep_idx,lag_idx(4));
idx_empty = randi(length(g_empty),length(g_empty),N);
idx_lpsp  = randi(length(g_lpsp),length(g_lpsp),N);
p2 = sum(mean(g_empty(idx_empty),1) - mean(g_lpsp(idx_lpsp),1) < 0)/ N;

subplot(2,2,1); cla
title({'Gain (slope):', 'Bump and Fly vel'},'color',scatter_color)
hold on
xticks([1,2,3,4]); xticklabels(group_labels)
xlim([.5,4.5])
fontsize(gcf,30,'pixels')

subplot(2,2,3); cla; hold on
for i = 1:4
    %slope
    m = mean(vel_slope(group_idx==i&keep_idx,:),1);
    s = std(vel_slope(group_idx==i&keep_idx,:),1)/sqrt(sum(group_idx==i&keep_idx));
    subplot(2,2,1); hold on
    scatter(group_idx(keep_idx&group_idx==i),vel_slope(keep_idx&group_idx==i,lag_idx(i)),100,scatter_color,'filled','MarkerFaceAlpha',.5)
    errorbar(i+.2,m(lag_idx(i)),s(lag_idx(i)),'Color',c(i,:),'Linewidth',2)
    
    subplot(2,2,3); hold on
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);
    
    %bias
    m = mean(vel_bias(group_idx==i&keep_idx,:),1);
    s = std(vel_bias(group_idx==i&keep_idx,:),1)/sqrt(sum(group_idx==i&keep_idx));
    subplot(2,2,2); hold on
    scatter(group_idx(keep_idx&group_idx==i),vel_bias(keep_idx&group_idx==i,lag_idx(i)),100,scatter_color,'filled','MarkerFaceAlpha',.5)
    errorbar(i+.2,m(lag_idx(i)),s(lag_idx(i)),'Color',c(i,:),'Linewidth',2)
    
    subplot(2,2,4); hold on
    patch([lag_vec,fliplr(lag_vec)],[m+s,fliplr(m-s)],c(i,:),'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off')
    plot(lag_vec,m,'Color',c(i,:),'Linewidth',2);

end

subplot(2,2,1)
tmp = ylim;
plot([1,2],tmp(2)*[1,1].*.9,scatter_color,'Linewidth',3)
text(1.5,tmp(2).*.9,sprintf('p = %.3f',p1),'color',scatter_color,'HorizontalAlignment','center','VerticalAlignment','bottom')
plot([3,4],tmp(2)*[1,1].*.9,scatter_color,'Linewidth',3)
text(3.5,tmp(2).*.9,sprintf('p = %.3f',p2),'color',scatter_color,'HorizontalAlignment','center','VerticalAlignment','bottom')
plot([.5,4.5],[0,0],':','color',scatter_color)

subplot(2,2,3);
xlabel('lag (ms)')
ylabel('(+- s.e.m.)')
axis tight

subplot(2,2,2)
title({'Bias:','Bump and Fly'},'color',scatter_color)
plot([.5,4.5],[0,0],':','color',scatter_color)
xticks([1,2,3,4]); xticklabels(group_labels)
xlim([.5,4.5])

subplot(2,2,4)
legend(group_labels,'Autoupdate','off','Location','Northeast')

fontsize(gcf,30,'pixels')

if dark_flag
    set(gcf,'color','none')
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
    try
        set(tmp(i),'color','none')
        set(tmp(i),'xcolor','w')
        set(tmp(i),'ycolor','w')
    end
    end
    legend(group_labels,'Autoupdate','off','Location','Northeast','textcolor','w', 'edgecolor','w')
end

%% show the polar histograms of fly and bump heading
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_names = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

rows = 4;
cols = 6;

f2 = 0;
for f = 1:length(unique(group_idx))
counter = 25;
for i = find(group_idx == f)
    counter = counter+1;

    if counter > 24
        f2 = f2+1;
        figure(f2); clf
        set(gcf,'Name',group_names{f})
        tiledlayout(rows,cols)
        counter = 1;
    end

    nexttile
    polarhistogram(full_data.mu{i})
    hold on
    polarhistogram(full_data.cue{i})
    thetaticklabels([])
    rticklabels([])
    title({[full_data.trial{i}(1:8)], ['Trial: ', num2str(i)]})

end
end

%% Compare the shape of bump and heading histograms, looking for hotspots
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
rho_thresh = .1;
rho_flag = 1;

% for no thresholding
H_mu  = cellfun(@(x)(entropy(round(x,2))),full_data.mu);
H_cue = cellfun(@(x)(entropy(round(x,2))),full_data.cue);

% for threshold
if rho_flag
H_mu = cellfun(@(x,y)(entropy(round(x(y>rho_thresh),2))),full_data.mu,full_data.rho);
H_cue = cellfun(@(x,y)(entropy(round(x(y>rho_thresh),2))),full_data.cue,full_data.rho);
end


figure(5); clf
subplot(1,2,1)
scatter(group_idx-.15,H_mu,100,'filled','c','MarkerFaceAlpha',.75)
hold on
scatter(group_idx+.15,H_cue,100,'filled','r','MarkerFaceAlpha',.75)
for i = 1:length(H_mu)
    plot([group_idx(i)-.15,group_idx(i)+.15]',[H_mu(i),H_cue(i)]','Color',[0.5,0.5,0.5,.75])
end
ylabel('Entropy (H)')
legend('mu','cue','Location','best')
xticks([1,2,3,4]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([.5,4.5])

subplot(1,2,2)
scatter(group_idx,H_mu-H_cue,100,'filled','m','MarkerFaceAlpha',.75)
hold on
plot([-.5,4.5],[0,0],':','color',scatter_color)
for i = 1:4
    m = mean(H_mu(group_idx==i)-H_cue(group_idx==i));
    s = std(H_mu(group_idx==i)-H_cue(group_idx==i))/sqrt(sum(group_idx==i));
    errorbar(i+.1,m,s,'Linewidth',2,'Color',scatter_color)
end

h_empty = H_mu(group_idx == 1 & keep_idx) - H_cue(group_idx==1 & keep_idx);
h_lpsp  = H_mu(group_idx == 2 & keep_idx) - H_cue(group_idx == 2 & keep_idx);
idx_empty = randi(length(h_empty),length(h_empty),N);
idx_lpsp  = randi(length(h_lpsp),length(h_lpsp),N);
p = sum(mean(h_empty(idx_empty),1) - mean(h_lpsp(idx_lpsp),1) > 0)/ N;
text(1.5,max(ylim),sprintf('p = %.3f',p),'HorizontalAlignment','Center','VerticalAlignment','top','Color',scatter_color)

h_empty = H_mu(group_idx == 3 & keep_idx) - H_cue(group_idx == 3 & keep_idx);
h_lpsp  = H_mu(group_idx == 4 & keep_idx) - H_cue(group_idx == 4 & keep_idx);
idx_empty = randi(length(h_empty),length(h_empty),N);
idx_lpsp  = randi(length(h_lpsp),length(h_lpsp),N);
p = sum(mean(h_empty(idx_empty),1) - mean(h_lpsp(idx_lpsp),1) > 0)/ N;
text(3.5,max(ylim),sprintf('p = %.3f',p),'HorizontalAlignment','Center','VerticalAlignment','top','Color',scatter_color)

ylabel('H_{mu}-H_{cue}')
xticks([1,2,3,4]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([.5,4.5])

fontsize(gcf,30,'pixels')

if dark_flag
    set(gcf,'color','none')
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
    try
        set(tmp(i),'color','none')
        set(tmp(i),'xcolor','w')
        set(tmp(i),'ycolor','w')
    end
    end
end

%% find hotspots in fluoresnce, not just bump
summed_f = cellfun(@(x)(sum(x./prctile(x,f0_pct,2),2,'omitnan')),full_data.f_cluster,'uniformoutput',false);
H = cellfun(@(x)(entropy(zscore(x))),summed_f);
figure(20)
scatter(group_idx,H)
xlim([.5,4.5])
xticks([1:4])
xticklabels(group_labels)
ylabel('H_{fluorescence}')

h_empty = H(group_idx == 1 & keep_idx);
h_lpsp  = H(group_idx == 2 & keep_idx);
idx_empty = randi(length(h_empty),length(h_empty),N);
idx_lpsp  = randi(length(h_lpsp),length(h_lpsp),N);

p = sum(mean(h_empty(idx_empty),1) - mean(h_lpsp(idx_lpsp),1) > 0)/ N
%% Functions
function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function [amp_tot,amp_peak, mu, rho,dff_cluster,f_cluster] = bump_calc_eb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
regProduct = smoothdata(regProduct(:,:,1:end,:),4,'gaussian',b_smooth);

[tmp_y,tmp_x] = find(bwmorph(mask,'shrink','inf'));
if length(tmp_y) == 1
    [y_mask,x_mask] = find(mask);
    figure(1); clf; imagesc(mask);
    tmp = drawellipse('Center',[tmp_x,tmp_y],'SemiAxes',[range(x_mask)/6,range(y_mask)/6]);
    mask = logical(mask - createMask(tmp));
end

%extract midline from mask
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
mid             = bwmorph(mask,'remove');
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);
hold on; scatter(x0,y0,'b','filled')

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size

if centroids(1,1) < centroids(end,1)
    centroids = flipud(centroids);
end

%assign each pixel to a centroid
[~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)

imgData = squeeze(sum(regProduct,3));
tmp     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
background  = tmp(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
%flash_idx   = med_pix > 0.05;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
%dff_cluster = zscore(dff_cluster,[],2);
zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
tmp = dff_cluster;
tmp(:,~flash_idx) = zscore_cluster;
zscore_cluster = tmp;
dff_cluster = zscore_cluster;

% n_planes = size(regProduct,3);
% filt_mu = linspace(max(centroids(:,1)),min(centroids(:,1)),n_planes);
% filt_sig= range(centroids(:,1))/(n_centroid/2);
% filt_x  = [1:size(regProduct,1)]';
% f_cluster = zeros(n_centroid,size(regProduct,4));
% 
% for p = 1:n_planes
% tmp = double(squeeze(regProduct(:,:,p,:))) .* normpdf(filt_x,filt_mu(p),filt_sig);
% 
% 
% imgData_2d      = reshape(double(tmp),[],size(regProduct,4));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
% centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
% for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
%     centroid_log(i, sub2ind(size(regProduct,1,2,4),y_mask(idx==i),x_mask(idx ==i))) = true;
% end
% 
% tmp       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
% f_cluster = f_cluster + tmp;
% end
% 
% f_cluster(:,flash_idx) = nan;
% f0              = prctile(f_cluster,f0_pct,2);            %find the baseline fluorescence in each cluster
% dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
% zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
% tmp = dff_cluster;
% tmp(:,~flash_idx) = zscore_cluster;
% dff_cluster = tmp;
alpha       = linspace(-pi,pi,n_centroid);

for i = 1:n_smooth
    dff_cluster = smoothdata(dff_cluster,2,'gaussian',b_smooth);
end

tmp = smoothdata(repmat(dff_cluster,3,1),1,'gaussian',3);
tmp = tmp(n_centroid+1:end-n_centroid,:);
[x_tmp,y_tmp]   = pol2cart(alpha,tmp');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
for i = 1:n_smooth
mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
rho    = smoothdata(rho,1,'gaussian',b_smooth);
end

mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
amp_tot = mean(dff_cluster,1)';
[~,i] = mink(abs(alpha-mu),2,2); %find indexes corresponding to peak of each time point
i2 = i' + size(dff_cluster,1)*[0:size(dff_cluster,2)-1]; %find the linear index into the peak of each column (time point) value. this was clever :)
amp_peak = mean(dff_cluster(i2),1)'; %extract peak amplitude

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
amp_tot = smoothdata(amp_tot,1,'gaussian',b_smooth); 
amp_peak = smoothdata(amp_peak,1,'gaussian',b_smooth);
end

total_t = max(xf);
xb  = linspace(0,total_t,size(regProduct,4))';

mu          = interp1(xb,unwrap(mu),xf)';
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp_tot     = interp1(xb,amp_tot,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';
end

function [f_speed,r_speed,intHD,cue,r_vel] = ft_calc(ftData_DAQ,n_smooth,f_smooth)

f_speed = ftData_DAQ.velFor{:};                       %store each speed
r_speed = ftData_DAQ.velYaw{:};
intHD   = ftData_DAQ.intHD{:};
cue     = ftData_DAQ.cuePos{:}';


f_speed = f_speed;                                      %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi
cue     = unwrap(cue / 192 * 2*pi - pi);
r_vel  = ftData_DAQ.velYaw{:};

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'gaussian',f_smooth); 
r_speed = smoothdata(r_speed,1,'gaussian',f_smooth);
intHD   = smoothdata(intHD,  1,'gaussian',f_smooth);
cue     = smoothdata(cue,    1,'gaussian',f_smooth);
r_vel   = smoothdata(r_vel,  1,'gaussian',f_smooth);
end

intHD = mod(intHD,2*pi);                                %rewrap heading data, and put between -pi and pi.
intHD(intHD > pi) = intHD(intHD > pi) - 2*pi;
cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;


% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
%    r_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
end