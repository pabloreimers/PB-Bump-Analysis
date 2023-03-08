%% script params
transfer_flag = false;
f0_pct              = 7;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 10;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 10;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 10;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;   
med_thresh          = 0.05;
pva_thresh          = 0;
%% 
base_dir = 'Z:\pablo\LPsP_CL\to analyze';
dest_dir = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\epg grabda testing\';
folders_list = cellstr(ls(base_dir));
epg_idx  = cellfun(@(x)(contains(x,'EPG')),folders_list);
folders_list = folders_list(epg_idx);

if transfer_flag
for i = 1:size(folders_list,1)
    fprintf('On %i of %i\n',i,size(folders_list,1))
    tmp_dir = [dest_dir,'\',folders_list{i}];
    if ~isfolder(tmp_dir)
    mkdir(tmp_dir)
    end
    
    name = ls([base_dir,'\',folders_list{i},'\*ficTracData_DAQ*']);
    copyfile([base_dir,'\',folders_list{i},'\',name(1,:)],tmp_dir)
    name = ls([base_dir,'\',folders_list{i},'\*imagingData*']);
    copyfile([base_dir,'\',folders_list{i},'\',folder_name,'\',name(1,:)],tmp_dir)
    name = ls([base_dir,'\',folders_list{i},'\*mask*']);
    if ~isempty(name)
    copyfile([base_dir,'\',folders_list{i},'\',name(1,:)],tmp_dir)
    end
end
end
%% analyze each trial
vel_corr = nan(size(folders_list,1),1);
disp_corr = nan(size(folders_list,1),1);
for i = 1:size(folders_list,1)
    if i == 39
        continue
    end
    fprintf('On %i of %i\n',i,size(folders_list,1))
    [~,~,~,~,~,~,~,~,~,~,~,~,vel_corr(i),disp_corr(i)] = single_analysis_da([dest_dir,'\',folders_list{i}],true);
end