%% Set script params
transfer_flag   = false;
mask_flag       = false;
local_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\2p data\';

for_thresh  = .5; %mean forward velocity to be included
time_thresh = .2; %amount of time above threshold forward velociy
ratio_thresh = 0; %ratio of forward to backward walking to be included

%% Get directory from user
base_dir        = uigetdir([],'Select the folder of dated folders to include in analysis');      %ask user for data folder
date_list       = dir(base_dir);   %create list of subfolders
date_list(1:2)  = [];              %delete not folders
trial_list      = {};

n = 0;
for i = 1:size(date_list,1)         %enter each day, and create a handle to the directory for each experiment
    
    tmp_dir      = dir([base_dir,'\',date_list(i).name]);
    tmp_dir(1:2) = [];
    for j = 1:size(tmp_dir,1)
        n = n+1;
        trial_list{n,1}   = [base_dir,'\',date_list(i).name,'\',tmp_dir(j).name];
    end
end

%% Assess locomotion of each trial, and update inclusion logical
n           = size(trial_list,1);
epg_idx     = cellfun(@(x)(contains(x(20:end),'EPG')),trial_list(:,1)); %create a handle to either epg or lpsp imaging
lpsp_idx    = cellfun(@(x)(contains(x(20:end),'LPsP')),trial_list(:,1));
dark_idx    = cellfun(@(x)(contains(x(20:end),'dark')),trial_list(:,1));
f_vel       = cell(n,1);

hold on
for i = 1:n
    fprintf('%.2f%%\n',i/n*100)
    name = ls([trial_list{i},'\*ficTracData_DAQ*']);
    if isempty(name)
        continue
    end
    load([trial_list{i,1},'\',name(1,:)])
    trial_list{i,2} = size(ftData_DAQ,1);

    f_vel{i} = ftData_DAQ(1,:).velFor{:};
    for k = 1:10
        f_vel{i} = smoothdata(f_vel{i},1,'gaussian',50);
    end    
end

%% Plot fictrac analysis to set threshold
figure(1); clf
ax(1) = subplot(3,1,1); hold on; ylabel('Mean Forward')
ax(2) = subplot(3,1,2); hold on; ylabel('Forward/Backward')
ax(3) = subplot(3,1,3); hold on; ylabel('Time > .5mm/s')

scatter(ax(1),nan,nan,'filled','r')
scatter(ax(1),nan,nan,'filled','b')
legend(ax(1),{'EPG','LPsP'},'AutoUpdate','off')

for i = 1:n
    if epg_idx(i)
        c = 'b';
    elseif lpsp_idx(i)
        c = 'r';
    end
    scatter(ax(1),i,mean(f_vel{i}),'filled',c)
    scatter(ax(2),i,sum(f_vel{i} > 0)/sum(f_vel{i} <0),'filled',c)
    scatter(ax(3),i,sum(f_vel{i} > for_thresh)/length(f_vel{i}),'filled',c)
end


%% Set logicals
idx_fun = @(x) sum(x > 0) / sum(x < 0) > ratio_thresh & ... %passes forward to backward velocity ratio
               sum(x > for_thresh)/length(x) > time_thresh; %passes percent time going forward thresh
ft_idx = cellfun(@(x) idx_fun(x),f_vel);                     %define which trials pass the locomotion thresholds


%% download the appropriate files to workplace
if transfer_flag
dest_dir = 'C:\Users\preim\Documents\Wilson Lab\data\to analyze';
tmp_idx = find(ft_idx);
time_vec = nan(size(tmp_idx));

for j = 25:length(tmp_idx)
    tic
    fprintf('%.2f%%....',j/length(tmp_idx)*100)
    i = tmp_idx(j);
    
    try
    tmp = regexp(trial_list{i},'\','split');
    tmp_dir = [dest_dir,'\',tmp{end}];
    if ~isfolder(tmp_dir)
    mkdir(tmp_dir)
    
    name = ls([trial_list{i},'\*ficTracData_DAQ*']);
    copyfile([trial_list{i,1},'\',name(1,:)],tmp_dir)
    folder_name = ls([trial_list{i},'\*registration*']);
    name = ls([trial_list{i},'\*registration*\*imagingData*']);
    copyfile([trial_list{i,1},'\',folder_name,'\',name(1,:)],tmp_dir)
    name = ls([trial_list{i},'\*mask*']);
    if ~isempty(name)
    copyfile([trial_list{i,1},'\',name(1,:)],tmp_dir)
    end
    end

    time_vec(j) = toc;
    catch
        fprintf(['failure for ',num2str(j)])
    end
    fprintf('Estimated TIme Remaining: %.2f minutes\n', nanmean(time_vec) * (length(tmp_idx) - j) / 60)
end
end


%% Go draw and save a mask for each image\
if mask_flag
local_dir = uigetdir(local_dir,'Select local directory of images');
files_list = dir(local_dir);
files_list(1:2) = [];

for i = 1:size(files_list,1)
    disp(i)
    if isempty(ls([local_dir,'\',files_list(i).name,'\*mask*']))
    name = ls([local_dir,'\',files_list(i).name,'\*imagingData*']);
    load([local_dir,'\',files_list(i).name,'\',name])

    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    imgData2    = imgData;
    top_int     = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
    bot_int     = prctile(imgData2,5,'all');
    imgData2    = max(min(imgData2,top_int),bot_int) - bot_int;
    imgData2    = 256*imgData2/max(imgData2,[],'all');
    %imgData2    = smoothdata(imgData2,3,'movmean',avg_win);                      %smooth the data, again for viewing purposes (should this go before the clipping)            
    
    figure(1); clf                                          % clear the current figure
    mask = roipoly(uint8(mean(imgData2,3)));               %this create a black and white mask (logical) of the PB, drawn on a maxZ projection of the image
    figure(2)
    imagesc(subplot(2,1,1),mask); colormap(bone); xticks([]); yticks([])

    save([local_dir,'\',files_list(i).name,'\mask.mat'],'mask')
    end
end
end

%% Create list of files in local directory
local_dir = uigetdir(local_dir,'Select local directory of images to analyze');
files_list = dir(local_dir);
files_list = struct2cell(files_list)';
files_list(1:2,:) = [];

files_info = cellfun(@(x)(regexp(x,'_','split')),files_list(:,1),'UniformOutput',false);
    
n           = size(files_list,1);
epg_idx     = cellfun(@(x)(contains(x,'EPG')),files_list(:,1)); %create a handle to either epg or lpsp imaging
lpsp_idx    = cellfun(@(x)(contains(x,'LPsP')),files_list(:,1));
dark_idx    = cellfun(@(x)(contains(x,'dark')),files_list(:,1));

T = nan(size(files_list,1),9);
D = {size(files_list,1),3};

for i = 1:n
    i/n
    [T(i,1),T(i,2),T(i,3),T(i,4),T(i,5),T(i,6),T(i,7),T(i,8),T(i,9),D{i,1},D{i,2},D{i,3}] = single_analysis([local_dir,'\',files_list{i,1}], true);
end


T = array2table(T,'VariableNames',{'vel_rho','vel_pval','pva_corr','pva_pval','lag','f_r2','r_r2','j_r2', 'num_flash'});

%% Label each trial by fly
[~,~,fly_id] = unique(str2num([cell2mat(cellfun(@(x)x{1}(1:8),files_info,'UniformOutput',false)),...
                                         cellfun(@(x)x{end},files_info),...
                                         num2str(cellfun(@(x)contains(x{2},'LPsP'),files_info))]));
