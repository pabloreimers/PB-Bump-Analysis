%% set script params
clear all
close all
mask_flag           = true;
mask_overwrite      = false;
movie_flag          = false;                                                %set whether to play movies or just run straight through the script
pause_time          = 0;                                                    %set the pause time if playing movies
data_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %set main data directory for ease of selecting files
%data_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\2p data\';
f0_pct              = 7;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 12;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 1;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = .1;
eb_flag             = false;
pb_flag             = true;

%% select folder of images
base_dir = uigetdir(data_dir);
all_files = dir([base_dir,'\**\*imagingData*']);

%% [NEW METHOD] open each trial, draw mask, and save
if mask_flag
for i = 4:length(all_files)
    fprintf('%s\n',all_files(i).folder)
    if ~isempty(dir([all_files(i).folder,'\*imagingData*'])) && (isempty(dir([all_files(i).folder,'\*mask*']))  || mask_overwrite)
    

    load([all_files(i).folder,'\',all_files(i).name])
    if size(regProduct,1) > size(regProduct, 2)
        regProduct  = permute(regProduct, [2 1 3 4]);
        save([tmp(j).folder,'\',tmp(j).name],"regProduct",'-v7.3')
    end


    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    imgData2    = imgData;
    top_int     = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
    bot_int     = prctile(imgData2,5,'all');
    imgData2    = max(min(imgData2,top_int),bot_int) - bot_int;
    imgData2    = 256*imgData2/max(imgData2,[],'all');
       
    
    figure(1); clf                                          % clear the current figure
    tmp2 = mean(imgData,3);
    tmp2 = 256*(tmp2 - min(tmp2(:)))/(max(tmp2(:) - min(tmp2(:))));
    image(tmp2); colormap(bone); drawnow;

    mask = roipoly(uint8(mean(imgData2,3)));  

    save([all_files(i).folder,'\mask.mat'],'mask')
    end
end
end

%% store the fictrac, total amplitude, and peak amplitude for each trial
all_trials = dir([base_dir,'\**\*imagingData*']);
n       = length(all_trials);
dff_tot = cell(n,1); %preallocate dff cell
dff_peak= cell(n,1);
dff_cluster = cell(n,1);
mu      = cell(n,1);
rho     = cell(n,1);
f_speed = cell(n,1);
r_speed = cell(n,1);
intHD   = cell(n,1);
cue     = cell(n,1);
r_vel   = cell(n,1);
c_speed = cell(n,1);
xf      = cell(n,1);
trial_names = cell(n,1);


figure(2); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

    date_dir = dir(base_dir);
    date_dir(1:2) = [];

for i = 1:n
    
    if ~isempty(dir([all_trials(i).folder,'\',all_trials(i).name])) && (~isempty(dir([all_trials(i).folder,'\','*mask*'])))
    
    fprintf('%s\n',all_trials(i).folder)

    load([all_trials(i).folder,'\',all_trials(i).name]) % load imaging data
    tmp = dir([all_trials(i).folder,'\*mask*']); %load mask
    load([tmp.folder,'\',tmp.name])
    tmp = dir([erase(all_trials(i).folder,'registration_001'),'\*ficTracData_DAQ*']); %load fictract
    load([tmp.folder,'\',tmp.name])
    
    
    [f_speed{i},r_speed{i},intHD{i},cue{i},r_vel{i},c_speed{i}] = ft_calc(ftData_DAQ,n_smooth,f_smooth);
    xf{i}  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed{i})]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
    if isduration(xf{i})
        xf{i} = seconds(xf{i});
    end

    if pb_flag
        [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i}] = bump_calc_pb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
    end

    nexttile
    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    tmp = mean(imgData,3);
    tmp = 256*(tmp - min(tmp(:)))/(max(tmp(:) - min(tmp(:))));
    image(tmp); colormap('bone'); axis equal tight; title(trial_names{i},'color','w')
    set(gca,'color','none')
    end
end
set(gcf,'color','none')

%% save variables
tmp = cellfun(@(x)(strsplit(x,'\')),{all_trials.folder}','UniformOutput',false);
trial_names = cellfun(@(x)(x{5}),tmp,'UniformOutput',false);
data = table(dff_tot,dff_peak,dff_cluster,mu,rho,f_speed,r_speed,intHD,cue,r_vel,c_speed,xf,trial_names);
save(['lpsp_ol_data_',num2str(yyyymmdd(datetime("today")))],"data")
%% save meta data for each fly

meta = ...
   {{'20231108-2_LPsP_syt7f'}, 1, 'm';...
    {'20231110-1_LPsP_syt7f'}, 2, 'm';...
    {'20231110-2_LPsP_syt7f'}, 2, 'm';...
    {'20231111-1_LPsP_syt7f'}, 3, 'm';...
    {'20231111-2_LPsP_syt7f'}, 3, 'm';...
    {'20231112-1_LPsP_syt7f'}, 4, 'm';...
    {'20231112-2_LPsP_syt7f'}, 4, 'm';...
    {'20231112-3_LPsP_syt7f'}, 5, 'g';...
    {'20231112-4_LPsP_syt7f'}, 5, 'g';...
    {'20231113-1_LPsP_syt7f'}, 6, 'g';...
    {'20231113-2_LPsP_syt7f'}, 6, 'g';...
    {'20231113-3_LPsP_syt7f'}, 7, 'g';...
    {'20231113-4_LPsP_syt7f'}, 7, 'g';...
    {'20231115-1_LPsP_syt7f'}, 8, 'g';...
    {'20231115-2_LPsP_syt7f'}, 8, 'g';...
    {'20231115-3_LPsP_syt7f'}, 9, 'g';...
    {'20231115-4_LPsP_syt7f'}, 9, 'g';...
    {'20231115-5_LPsP_syt7f'}, 10,'g';...
    {'20231115-6_LPsP_syt7f'}, 10,'g';...
    {'20231115-7_LPsP_syt7f'}, 11,'g';...
    {'20231121-1_LPsP_syt7f'}, 12,'m';...
    {'20231121-2_LPsP_syt7f'}, 12,'m';...
    {'20231121-3_LPsP_syt7f'}, 13,'m';...
    {'20231121-4_LPsP_syt7f'}, 14,'m';...
    {'20231121-5_LPsP_syt7f'}, 14,'m';...
    {'20231121-6_LPsP_syt7f'}, 15,'m';...
    {'20231121-7_LPsP_syt7f'}, 15,'m';...
    {'20231122-1_LPsP_syt7f'}, 16,'m';...
    {'20231127-1_LPsP_syt7f'}, 17,'m';...
    {'20231127-2_LPsP_syt7f'}, 17,'m';...
    {'20231127-4_LPsP_syt7f'}, 18,'m';...
    {'20231204-1_LPsP_syt7f'}, 19,'m';...
    {'20231204-2_LPsP_syt7f'}, 19,'m';...
    {'20231204-3_LPsP_syt7f'}, 20,'m';...
    {'20231204-4_LPsP_syt7f'}, 20,'m';...
    {'20231204-5_LPsP_syt7f'}, 21,'m';...
    {'20231204-6_LPsP_syt7f'}, 21,'m';...
    {'20231204-7_LPsP_syt7f'}, 22,'m';...
    {'20231204-8_LPsP_syt7f'}, 22,'m';...
    {'20231206-1_LPsP_syt7f'}, 23,'m';...
    {'20231206-2_LPsP_syt7f'}, 23,'m';...
    {'20231215-1_LPsP_syt7f'}, 24,'g';...
    {'20231215-2_LPsP_syt7f'}, 24,'g';...
    {'20231215-3_LPsP_syt7f'}, 25,'g';...
    {'20231215-4_LPsP_syt7f'}, 26,'g'};%;...;

fly_num = cell2mat(meta(:,2));
fly_food= cellfun(@(x)(strcmp(x,'m')),meta(:,3));

%%
n = size(data,1);

figure(3); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

figure(4); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

figure(5); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

figure(6); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

for i = 1:size(data,1)
still_idx = data.r_speed{i} < 1e-1; %mark when the fly isnt moving
trans_idx = zeros(length(still_idx),1); 

switch fly_food(i)
    case 1
        c = 'w';
    case 0
        c = 'm';
end

%trans_idx = smoothdata(~[diff(data.c_speed{i})==0;0],'movmean',60)>0; %mark all time points around transitions, where it's hard to attribute a signal to a specific cue speed

figure(3)
nexttile
scatter(data.c_speed{i}(still_idx& ~trans_idx),data.dff_tot{i}(still_idx& ~trans_idx),5,'filled','MarkerFaceAlpha',.1)
xlabel('OL speed'); ylabel('dF/F')
title(data.trial_names{i},'color',c)
%text(max(c_speed{i}(still_idx)),min(dff_tot{i}(still_idx)),'while fly is stopped','VerticalAlignment','bottom','HorizontalAlignment','right')

figure(4)
nexttile
scatter(data.r_speed{i}(~trans_idx),data.dff_tot{i}(~trans_idx),5,data.c_speed{i}(~trans_idx),'filled','markerfacealpha',.1)
xlabel('fly speed'); ylabel('dF/F')
title(data.trial_names{i},'color',c)

% tmp = colorbar;
% ylabel(tmp,'OL speed')


tmp1 = data.c_speed{i}*10+1;
tmp3 = accumarray(tmp1(still_idx & ~trans_idx),data.dff_tot{i}(still_idx& ~trans_idx),[],@mean,nan); %,[size(tmp1),size(tmp2)],@mean);

[tmp_speeds,ia,ic] = unique(tmp1(still_idx&~trans_idx));
[~,tmp_order] = sort(ia);
tmp_speeds = (tmp_speeds-1)/10;
tmp3(isnan(tmp3)) = [];

figure(5)
nexttile
scatter(tmp_speeds,tmp3)
hold on
plot(xlim,tmp3(tmp_speeds==0)*[1,1],'k:')
ylabel('mean dF/F'); xlabel('OL speed')
title(data.trial_names{i},'color',c)

figure(6)
nexttile
scatter(tmp_order,tmp3)
hold on
plot(xlim,tmp3(tmp_speeds==0)*[1,1],'k:')
ylabel('mean dF/F'); xlabel('OL order')
title(data.trial_names{i},'color',c)
end



%% establish signifance

%for each fly, find the dff while it is still for with and without optic
%flow

dff_still = cellfun(@(x,y,z) mean(x(y==0 & z<0.1)),...
                    data.dff_tot,data.c_speed,data.r_speed);
dff_flow = cellfun(@(x,y,z) mean(x(y>0 & z<0.1)),...
                    data.dff_tot,data.c_speed,data.r_speed);

figure(10); clf
subplot(2,2,1); hold on
plot([dff_still(fly_food),dff_flow(fly_food)]','-ow','linewidth',2)
plot([dff_still(~fly_food),dff_flow(~fly_food)]','-om','linewidth',2)
xlim([.5,2.5])
ylabel('mean dF/F'); xticks([1,2]); xticklabels({'still','flow'})
title('all trials','color','w')

subplot(2,2,2); hold on
scatter(fly_num(fly_food),dff_flow(fly_food)-dff_still(fly_food),'w','linewidth',2);
scatter(fly_num(~fly_food),dff_flow(~fly_food)-dff_still(~fly_food),'m','linewidth',2);
hold on
plot(xlim,[0,0],':w','linewidth',2)
xlabel('fly'); ylabel('dF/F difference (flow-still)')
axis tight 
y = ylim; ylim([min(-.1,y(1)),y(2)])
title('all trials','color','w')

subplot(2,2,3); hold on
plot([accumarray(fly_num(fly_food),dff_still(fly_food),[],@mean,nan),accumarray(fly_num(fly_food),dff_flow(fly_food),[],@mean,nan)]','-ow','linewidth',2)
plot([accumarray(fly_num(~fly_food),dff_still(~fly_food),[],@mean,nan),accumarray(fly_num(~fly_food),dff_flow(~fly_food),[],@mean,nan)]','-om','linewidth',2)
xlim([.5,2.5])
ylabel('mean dF/F'); xticks([1,2]); xticklabels({'still','flow'})
title('fly average','color','w')

subplot(2,2,4); hold on
dff_diff = accumarray(fly_num,dff_flow-dff_still,[],@mean,nan);
food_ind = logical(accumarray(fly_num,fly_food,[],@mean,nan));
num_vec = 1:max(fly_num);
scatter(num_vec(food_ind),dff_diff(food_ind),'w','linewidth',2);
scatter(num_vec(~food_ind),dff_diff(~food_ind),'m','linewidth',2);
hold on
plot(xlim,[0,0],':w','linewidth',2)
xlabel('fly'); ylabel('dF/F difference (flow-still)')
axis tight
y = ylim; ylim([min(-.1,y(1)),y(2)])
title('fly average','color','w')

fontsize(gcf,20,'pixels')
set(gcf,'color','none')
set(get(gcf,'Children'),'color','none','YColor','w','XColor','w')
pos = get(gca,'Position');
leg = legend('molasses','german','textcolor','w','edgecolor','w');
leg.Position(1:2) = [0.5,0.5] - 0.5*leg.Position(3:4);
set(gca,'Position',pos)

%% bin data by optic flow speed, and plot dff vs fly speed
fs_vec = 0:.1:max(cellfun(@max,data.r_speed)); %create a vector of fly rspeeds that we can interpolate everything on to for easy averaging
os_vec = unique(data.c_speed{1});

dff_mat = nan(n,length(fs_vec),length(os_vec));

for i = 1:length(os_vec)
    dff_mat(:,:,i) = cell2mat(cellfun(@(dff,os,fs) ...
                   (interp1(fs(os==os_vec(i)),dff(os==os_vec(i)),fs_vec)),...
                   data.dff_tot,data.c_speed,data.r_speed,'UniformOutput',  false));
    % dff_mat(:,:,i) = cell2mat(cellfun(@(dff,os,fs) ...
    %            (accumarray(@mean,round(fs(os==os_vec(i)),1)*10+1,dff(os==os_vec(i)))),...
    %            data.dff_tot,data.c_speed,data.r_speed,'UniformOutput',  false));
end

figure(13); clf
plot(fs_vec,squeeze(mean(dff_mat(fly_food,:,:),1,'omitnan'))')
colororder(cool(size(dff_mat,3)))
xlabel('fly speed'); ylabel('dff')
legend(num2str(os_vec))

figure(14); clf
for i = 1:n
    subplot(7,7,i);
    plot(squeeze(dff_mat(i,:,:)))
    colororder(cool(size(dff_mat,3)))
end
linkaxes(get(gcf,'Children'),'x')
%% plot dff averaged over all flies aligned to optic flow
n = size(data,1);
delay_vec = nan(n,1);
df_aligned = nan(n,size(data.c_speed{1},1));
cs_aligned = nan(n,size(data.c_speed{1},1));
si_aligned = false(n,size(data.c_speed{1},1));

for i = 1:n
    delay_vec(i) = finddelay(data.c_speed{1},data.c_speed{i});
    still_idx = abs(data.r_speed{i}) < .1;

    if delay_vec(i) > 0
       cs_aligned(i,1:end-delay_vec(i)+1) = data.c_speed{i}(delay_vec(i):end);
       df_aligned(i,1:end-delay_vec(i)+1) = data.dff_tot{i}(delay_vec(i):end);
       si_aligned(i,1:end-delay_vec(i)+1) = still_idx(delay_vec(i):end);
    else
       cs_aligned(i,-delay_vec(i)+1:end) = data.c_speed{i}(1:end-delay_vec(i):end);
       df_aligned(i,-delay_vec(i)+1:end) = data.dff_tot{i}(1:end-delay_vec(i):end);
       si_aligned(i,-delay_vec(i)+1:end) = still_idx(1:end-delay_vec(i):end);
    end

end

figure(12); clf

subplot(4,2,1); hold on; title('Molasses Food','color','w')
plot(data.xf{1},mean(cs_aligned(fly_food,:),1)','color',[0,.5,1])
ylabel('OL Speed')

subplot(4,2,[3,5,7]); hold on
x = data.xf{1};
y = mean(df_aligned(fly_food,:),1);
s = std(df_aligned(fly_food,:),1)./sqrt(sum(~isnan(df_aligned(fly_food,:)),1));
idx = ~isnan(x.*y.*s);
x = x(idx);
y = y(idx);
s = s(idx);
c = data.c_speed{1}(idx)' > 0;

patch([x,fliplr(x)],[y+s,fliplr(y-s)],.5*[1,1,1])
patch([x,fliplr(x)],[c,zeros(size(c))],[0,.5,1],'FaceAlpha',.5)
plot(x,y,'w')

ylabel('mean dff')
xlabel('time (s)')


subplot(4,2,2); hold on; title('German Food','color','w')
plot(data.xf{1},mean(cs_aligned(~fly_food,:),1)','color',[0,.5,1])

subplot(4,2,[4,6,8]); hold on
x = data.xf{1};
y = mean(df_aligned(~fly_food,:),1);
s = std(df_aligned(~fly_food,:),1)./sqrt(sum(~isnan(df_aligned(fly_food,:)),1));
idx = ~isnan(x.*y.*s);
x = x(idx);
y = y(idx);
s = s(idx);
c = data.c_speed{1}(idx)' > 0;

patch([x,fliplr(x)],[y+s,fliplr(y-s)],.5*[1,1,1])
patch([x,fliplr(x)],[c,zeros(size(c))],[0,.5,1],'FaceAlpha',.5)
plot(x,y,'w')

xlabel('time (s)')

linkaxes(get(gcf,'Children'),'x')
set(get(gcf,'Children'),'color','none','ycolor','w','xcolor','w')
set(get(gcf,'Children'),'color','none','ycolor','w','xcolor','w')
set(gcf,'color','none')

figure(13); clf
df_aligned(~si_aligned) = nan;

subplot(4,2,1); hold on; title('Molasses Food','color','w')
plot(data.xf{1},mean(cs_aligned(fly_food,:),1)','color',[0,.5,1])
ylabel('OL Speed')

subplot(4,2,[3,5,7]); hold on
x = data.xf{1};
y = mean(df_aligned(fly_food,:),1,'omitnan');
s = std(df_aligned(fly_food,:),1,'omitnan')./sqrt(sum(~isnan(df_aligned(fly_food,:)),1,'omitnan'));
idx = ~isnan(x.*y.*s);
x = x(idx);
y = y(idx);
s = s(idx);
c = data.c_speed{1}(idx)' > 0;

patch([x,fliplr(x)],[y+s,fliplr(y-s)],.5*[1,1,1])
patch([x,fliplr(x)],[c,zeros(size(c))],[0,.5,1],'FaceAlpha',.5)
plot(x,y,'w')

ylabel('mean dff')
xlabel('time (s)')


subplot(4,2,2); hold on; title('German Food','color','w')
plot(data.xf{1},mean(cs_aligned(~fly_food,:),1)','color',[0,.5,1])

subplot(4,2,[4,6,8]); hold on
x = data.xf{1};
y = mean(df_aligned(~fly_food,:),1,'omitnan');
s = std(df_aligned(~fly_food,:),1,'omitnan')./sqrt(sum(~isnan(df_aligned(fly_food,:)),1,'omitnan'));
idx = ~isnan(x.*y.*s);
x = x(idx);
y = y(idx);
s = s(idx);
c = data.c_speed{1}(idx)' > 0;

patch([x,fliplr(x)],[y+s,fliplr(y-s)],.5*[1,1,1])
patch([x,fliplr(x)],[c,zeros(size(c))],[0,.5,1],'FaceAlpha',.5)
plot(x,y,'w')

xlabel('time (s)')

linkaxes(get(gcf,'Children'),'x')
set(get(gcf,'Children'),'color','none','ycolor','w','xcolor','w')
set(get(gcf,'Children'),'color','none','ycolor','w','xcolor','w')
set(gcf,'color','none')

%% functions

function [amp_tot,amp_peak, mu, rho,dff_cluster] = bump_calc_pb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
tmp = bwmorph(mask,'thicken',20);
background  = imgData(~reshape(tmp,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = zeros(size(med_pix)); %med_pix > median(med_pix) + 2*std(med_pix);
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = imgData - prctile(background,2,'all');
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length
min_axis        = min(range(x_mask),range(y_mask));
mid             = bwskel(mask,'MinBranchLength',min_axis);  %find the midline as the skeleton, shaving out all sub branches that are smaller than the minimum axis length
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!

xq          = [-min_axis:(length(x_mid)+min_axis)];                             %extend the midline so that it reaches the border of the mask. extrapolate as many points as the minimum axis length
x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');                 %keep only the points that exist within the mask
x_mid       = x_mid(idx);
y_mid       = y_mid(idx);

xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';                        %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                     %take every other so that we dont start at the edges, and all are same size
%assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)

imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
%f_cluster       = f_cluster  - median(bacground,1);
f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;

alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);

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

function [f_speed,r_speed,intHD,cue,r_vel,c_speed] = ft_calc(ftData_DAQ,n_smooth,f_smooth)

f_speed = ftData_DAQ.velFor{:};                       %store each speed
r_speed = ftData_DAQ.velYaw{:};
intHD   = ftData_DAQ.intHD{:};
cue     = ftData_DAQ.cuePos{:}';

f_speed = f_speed;                                      %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi
cue     = unwrap(cue / 192 * 2*pi - pi);

c_speed = movsum([diff(unwrap(ftData_DAQ.cuePos{:}')/192*2*pi);0],60);
c_speed = smoothdata(c_speed,'movmedian',1e3);



for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'gaussian',f_smooth); 
r_speed = smoothdata(r_speed,1,'gaussian',f_smooth);
intHD   = smoothdata(intHD,  1,'gaussian',f_smooth);
cue     = smoothdata(cue,    1,'gaussian',f_smooth);
end

c_speed = round(abs(c_speed),1);


intHD = mod(intHD,2*pi);                                %rewrap heading data, and put between -pi and pi.
intHD(intHD > pi) = intHD(intHD > pi) - 2*pi;
cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;

r_vel  = ftData_DAQ.velYaw{:};
% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
%    r_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
end