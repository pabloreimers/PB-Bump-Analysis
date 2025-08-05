
%%
all_files = dir('Z:\pablo\hackathon_gain_change\**\registered_movie_z_proj_ch1.tif');
all_files = natsortfiles(all_files);

for i = 1:length(all_files)
    fprintf('noah2pablo %s\n',all_files(i).folder)
    if ~isfile([all_files(i).folder,'\imagingData.mat'])
        imgData = tiffreadVolume([all_files(i).folder,'\',all_files(i).name]);
        imgData = rot90(imgData,2);
        save([all_files(i).folder,'\imagingData.mat'],'imgData')
    end
end

%%
all_files = dir('Z:\pablo\lpsp_cschrimson_redo\**\imagingData_reg*.mat');
all_files = natsortfiles(all_files);

for i = 1:length(all_files)
    fprintf('pablo2pablo %s\n',all_files(i).folder)
    if ~isfile([all_files(i).folder,'\imagingData.mat'])
        load([all_files(i).folder,'\',all_files(i).name]);
        imgData = squeeze(sum(regProduct,3));
        save([all_files(i).folder,'\imagingData.mat'],'imgData')
    end
end