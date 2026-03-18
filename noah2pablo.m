
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
all_files = dir('Z:\pablo\lpsp_lateralized\**\imagingData*.mat');
all_files = natsortfiles(all_files);
% idx = cellfun(@(x)(contains(x,'exclusions') | contains(x,'open loop')),{all_files.folder}');
% all_files(idx) = [];

for i = 1:length(all_files)
    fprintf('pablo2pablo %s\n',all_files(i).folder)
    try
    if ~isfile([all_files(i).folder,'\imagingData.mat']) || ~isfile([all_files(i).folder,'\imgData_reg.mat'])
        load([all_files(i).folder,'\',all_files(i).name]);
        % img{1} = squeeze(sum(img{1},3));
        % try img{2} = squeeze(sum(img{2},3)); end
        % save([all_files(i).folder,'\imagingData.mat'],'img','-v7.3')
        imgData = squeeze(sum(regProduct,3));
        save([all_files(i).folder,'\imagingData.mat'],'imgData')
        imgData = normcorre_regProduct(imgData,false);
        save([all_files(i).folder,'\imgData_reg.mat'],'imgData')
    end
    catch
        fprintf('Error!!!!\n')
    end
end

%%
base_dir = 'Z:\pablo\epg_dlight\**\postreg_ch0'; %uigetdir(); %
all_files = dir([base_dir,'\**\registered_movie_ch0.tif']);
all_files = natsortfiles(all_files);

for i = 1:length(all_files)
    fprintf('%i / %i\n',i,length(all_files))

    curr_file = [all_files(i).folder,'\',all_files(i).name];

    info = imfinfo(curr_file);
    numSlices = numel(info);

    imgData = imread(curr_file, 1);
    imgData = repmat(imgData,1,1,numSlices);

    for k = 2:numSlices
        imgData(:,:,k) = imread(curr_file, k);
    end

    save([all_files(i).folder,'\imgData_denoised.mat'],'imgData')
end

%%

for i = 2:length(all_files)
    movefile([all_files(i).folder,'\imgData_denoised.mat'],[fileparts(fileparts(all_files(i).folder)),'\registration_001\'])
end