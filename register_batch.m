%% find all files
base_dir = uigetdir(); %('Z:\pablo\lpsp_p2x2\todo\');
all_files = dir([base_dir,'\**\*imagingData*.mat']);
all_files = natsortfiles(all_files);

%% check which files have a saved registered image

tic
for i = 1:length(all_files)
    fprintf('processing: %s ',all_files(i).folder)
    
    clear regProduct img
    if isempty(dir([all_files(i).folder,'\*imgData_smooth_reg*.mat']))
        
        
        tmp = dir([all_files(i).folder,'\*imagingData*.mat']);
        load([tmp.folder,'\',tmp.name])
        
        if exist('img','var')
            regProduct = img{1};
        end

        imgData_smooth_reg = normcorre_regProduct(smoothdata(squeeze(sum(regProduct,3)),3,'movmean',5),false); %save motion correction of the image summed across all planes
        
        save([all_files(i).folder,'\imgData_smooth_reg.mat'],'imgData_smooth_reg','-v7.3')
    end
    fprintf('ETR: %.2f mins\n',toc/i * (length(all_files)-i) / 60 )
end