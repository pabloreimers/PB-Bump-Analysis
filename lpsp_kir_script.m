%% set experiment params
transfer_flag = false;
local_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\analysis\'; %define the local directory to clear up
local_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\lpsp_kir\analysis\';
remote_dir           = '\\files.med.harvard.edu\Neurobio\wilsonlab\pablo\lpsp_kir\analysis'; %define the remote directory to make sure things are backed up to

f0_pct              = 7;                                                    %set the percentile that for the baseline fluorescence
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
rho_thresh          = 2e-1;
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
mu      = {};
rho     = {};
f_speed = {};
r_speed = {};
intHD   = {};
cue     = {};
r_vel   = {};
x_pos   = {};
y_pos   = {};
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
            [f_speed{i},r_speed{i},intHD{i},cue{i},r_vel{i},x_pos{i},y_pos{i}] = ft_calc(ftData_DAQ,n_smooth,f_smooth);
            xf{i}  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed{i})]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
            if isduration(xf{i})
                xf{i} = seconds(xf{i});
            end
            [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i}] = bump_calc_eb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
            fprintf('\n')
        end
    end
end

%% store everything in a table for easy access
full_data = table(dff_tot, dff_peak, dff_cluster, mu, rho, f_speed,r_speed,intHD,...
    cue, r_vel, x_pos, y_pos, xf,cl_idx,exp_idx,fly_num,exp_date,trial);

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
    %xlabel(['vel corr: ',num2str(vel_corr(i,1))])
    title(full_data.trial{i}(1:8))
    ylabel(['trial num: ',num2str(i)])
end
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
    hold on
    tmp = circ_dist(-full_data.cue{i},full_data.mu{i});
    a = plot(tmp,'m','Linewidth',2);
    a.YData(abs(diff(a.YData))>pi) = nan;
    a.YData(full_data.rho{i} < rho_thresh) = nan;

    a = plot(-full_data.cue{i},'k','Linewidth',.5);
    a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(full_data.mu{i},'b','Linewidth',.5);
    a.YData(abs(diff(a.YData))>pi) = nan;
    a.YData(full_data.rho{i} < 0) = nan;
    xticks([])
    xlabel(['offset var: ',num2str(circ_var(tmp(~isnan(tmp)&full_data.rho{i} > rho_thresh)))])
    title(full_data.trial{i}(1:8))
    ylabel(['trial num: ',num2str(i)])
end
end

%% show the cartesian trajectory of eachfly by group
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
    hold on
    scatter(x_pos{i},y_pos{i},[],1:length(x_pos{i}),'.'); axis equal tight
    xlabel('x pos (mm)'); ylabel('y pos (mm)')
    xticks([])
    %xlabel(['offset var: ',num2str(circ_var(tmp(~isnan(tmp)&full_data.rho{i} > rho_thresh)))])
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
            colormap(parula)
            axis equal tight
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
    tmp = abs(fly_vel) > vel_thresh & abs(bump_vel) < 10;
    vel_corr(i,j) = corr(bump_vel(tmp),fly_vel(tmp));
end
end

for f = 1:length(unique(group_idx))
figure(f); clf
set(gcf,'Name',group_names{f})
counter = 0;
rows = floor(sqrt(sum(group_idx==f)));
cols = ceil(sum(group_idx==f)/rows);
for i = find(group_idx == f)
    counter = counter+1;
    
    [~,lag] = max(vel_corr,[],2);
    mu_tmp = mu{i}(1+lag:end);
    rho_tmp = rho{i}(1+lag:end);
    vel_tmp = r_vel{i}(1:end-lag);
    
    idx = abs(vel_tmp) > vel_thresh & ~isnan(mu_vel) ;
    
    g = mu_vel(idx) \ [ones(length(r_vel{i}(idx)),1),r_vel{i}(idx)];
    
    subplot(rows,cols,counter)
    scatter(r_vel{i}(idx),mu_vel(idx),2,'filled','Markerfacealpha',0.1)
    hold on
    x = xlim;
    y = ylim;
    plot(x,[0,0],':k')
    plot([0,0],y,':k')
    plot(x,x*g(2)+g(1),'r' )
    xlim(x); ylim(y)
    xlabel(vel_corr(i,lag))
end
end

%% show the polar histograms of fly and bump heading
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx + 1;
group_names = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

for f = 1:length(unique(group_idx))
figure(f); clf
set(gcf,'Name',group_names{f})
counter = 0;
rows = floor(sqrt(sum(group_idx==f)));
cols = ceil(sum(group_idx==f)/rows);
for i = find(group_idx == f)
    counter = counter+1;
    subplot(rows,cols,counter)
    polarhistogram(mu{i})
    hold on
    polarhistogram(-cue{i})
    % a.YData(abs(diff(a.YData))>pi) = nan;
    % a = plot(mu{i},'w');
    % a.YData(abs(diff(a.YData))>pi) = nan;
    % a.YData(rho{i} < prctile(rho{i},50)) = nan;
    % xticks([])
    % xlabel(offset_var(i))
    title(trial{i})
end
end

%% compare offset variability across groups
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx;

offset_var = nan(size(exp_idx,2),20);
win_size = 120; %in seconds

for j =1:20
lag = j;
for i = 1:length(offset_var)
    mu_tmp = full_data.mu{i}(1+lag:end);
    rho_tmp = full_data.rho{i}(1+lag:end);
    cue_tmp = -full_data.cue{i}(1:end-lag);    

    tmp = ~isnan(mu_tmp) & (rho_tmp > rho_thresh) & [abs(diff(cue_tmp));0] > 0;
    tmp = mod(smoothdata(unwrap(circ_dist(mu_tmp(tmp),cue_tmp(tmp))),1,'gaussian',f_smooth)+pi,2*pi) - pi;
    offset_var(i,j) = circ_var(tmp(~isnan(tmp)));
end
end

keep_idx = cellfun(@(x)(sum(x>.2)/length(x) > .1),full_data.rho) & cellfun(@(x)(sum(x>0)/length(x))>.6,full_data.f_speed);

%offset_var = min(offset_var,[],2).*ones(size(offset_var));

[~,lag] = max(mean(offset_var(group_idx == 0 | group_idx == 2,:),1));
[h1,p1] = ttest2(offset_var(group_idx==0 & keep_idx),offset_var(group_idx==1 & keep_idx),'Tail','left','Vartype','unequal'); %test that the correlation for empty is greater than for LPsP
[h2,p2] = ttest2(offset_var(group_idx==2 & keep_idx),offset_var(group_idx==3 & keep_idx),'Tail','left','Vartype','unequal'); %test that the correlation for empty is greater than for LPsP

figure(1); clf
scatter(group_idx,offset_var(:,1),100,'filled','k','MarkerFaceAlpha',.5)
hold on
tmp = ylim;
plot([0,1],tmp(2)*[1,1].*.9,'k','Linewidth',3)
text(.5,tmp(2).*.9,sprintf('p = %.3f',p1),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')
plot([2,3],tmp(2)*[1,1].*.9,'k','Linewidth',3)
text(2.5,tmp(2).*.9,sprintf('p = %.3f',p2),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')
plot([-5,3.5],[0,0],':k')
ylabel('Offset Variability')
xticks([0,1,2,3]); xticklabels({'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'})
xlim([-.5,3.5])
% set(gcf,'color','none')
% set(gca,'color','none','xcolor','w','ycolor','w')
fontsize(gcf,30,'pixels')

%% compare pearson correlation of velocities
exp_idx = cellfun(@(x)contains(x,'LPsP'),full_data.trial);
cl_idx = cellfun(@(x)contains(x,'cl'),full_data.trial);
group_idx = (cl_idx*2) + exp_idx;
group_labels = {'Dark, empty','Dark, LPsP','CL, empty','CL, LPsP'};

vel_corr = nan(size(exp_idx,2),20);

for j =1:20
lag = j;
for i = 1:length(vel_corr)
    mu_tmp = mu{i}(1+lag:end);
    rho_tmp = rho{i}(1+lag:end);
    fly_vel = r_vel{i}(1:end-lag);

    bump_vel = [diff(mu_tmp);0]./fr;
    
    tmp = ~isnan(bump_vel) & (rho_tmp > rho_thresh) & abs(fly_vel) < vel_thresh & abs(fly_vel) > vel_min & abs(bump_vel) < 10 & [abs(diff(cue{i}(1:end-lag)));0] > 0;
    
    vel_corr(i,j) = corr(bump_vel(tmp),fly_vel(tmp));
end
end

figure(13); clf
hold on
[~,lag] = max(vel_corr,[],2);
c = [1,.5,0;0,.5,1;.75,0,.75;0,1,0];
colormap(c)
for i = 1:length(colormap)
    scatter(nan,1,[],c(i,:),'filled')
end
legend(group_labels,'Autoupdate','off')
a = scatter(1:length(lag),lag*.0167*1000,100,group_idx+1,'filled');
ylabel('optimal lag (ms)')
xlabel('trial number')


keep_idx = cellfun(@(x)(sum(x>rho_thresh)/length(x) > rho_thresh),rho);

vel_corr = max(vel_corr,[],2).*ones(size(vel_corr));

[~,lag] = max(mean(vel_corr(group_idx == 0 | group_idx == 2,:),1));
[h1,p1] = ttest2(vel_corr(group_idx==0 & keep_idx,lag),vel_corr(group_idx==1 & keep_idx,lag),'Tail','right','Vartype','unequal'); %test that the correlation for empty is greater than for LPsP
[h2,p2] = ttest2(vel_corr(group_idx==2 & keep_idx,lag),vel_corr(group_idx==3 & keep_idx,lag),'Tail','right','Vartype','unequal'); %test that the correlation for empty is greater than for LPsP


figure(4); clf
scatter(group_idx(keep_idx),vel_corr(keep_idx,lag),100,'filled','k','MarkerFaceAlpha',.5)
hold on
tmp = ylim;
plot([0,1],tmp(2)*[1,1].*.9,'k','Linewidth',3)
text(.5,tmp(2).*.9,sprintf('p = %.3f',p1),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')
plot([2,3],tmp(2)*[1,1].*.9,'k','Linewidth',3)
text(2.5,tmp(2).*.9,sprintf('p = %.3f',p2),'color','k','HorizontalAlignment','center','VerticalAlignment','bottom')
plot([-5,3.5],[0,0],':k')
ylabel('Pearson CorrCoeff: Bump and Fly vel')
xticks([0,1,2,3]); xticklabels(group_labels)
xlim([-.5,3.5])
% set(gcf,'color','none')
% set(gca,'color','none','xcolor','w','ycolor','w')
fontsize(gcf,30,'pixels')



function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function [amp_tot,amp_peak, mu, rho,dff_cluster] = bump_calc_eb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
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

%find the mean activity in each cluster
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
background  = imgData(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
%flash_idx   = med_pix > 0.05;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

n_planes = size(regProduct,3);
filt_mu = linspace(max(centroids(:,1)),min(centroids(:,1)),n_planes);
filt_sig= range(centroids(:,1))/(n_centroid/2);
filt_x  = [1:size(regProduct,1)]';
f_cluster = zeros(n_centroid,size(regProduct,4));

for p = 1:n_planes
tmp = double(squeeze(regProduct(:,:,p,:))) .* normpdf(filt_x,filt_mu(p),filt_sig);


imgData_2d      = reshape(double(tmp),[],size(regProduct,4));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(regProduct,1,2,4),y_mask(idx==i),x_mask(idx ==i))) = true;
end

tmp       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f_cluster = f_cluster + tmp;
end

f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);            %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
tmp = dff_cluster;
tmp(:,~flash_idx) = zscore_cluster;
dff_cluster = tmp;
alpha       = linspace(-pi,pi,n_centroid);

for i = 1:n_smooth
    dff_cluster = smoothdata(dff_cluster,2,'gaussian',b_smooth);
end

[x_tmp,y_tmp]   = pol2cart(alpha,dff_cluster');
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

function [f_speed,r_speed,intHD,cue,r_vel,x_pos,y_pos] = ft_calc(ftData_DAQ,n_smooth,f_smooth)

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

f_vel = ftData_DAQ.velFor{:};
s_vel = ftData_DAQ.velSide{:};
theta = ftData_DAQ.cueAngle{:} / 180 * pi;

x_pos = cumsum(f_vel.*cos(theta') + s_vel.*sin(theta'))/60;
y_pos = cumsum(f_vel.*sin(theta') + s_vel.*cos(theta'))/60;


% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
%    r_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
end