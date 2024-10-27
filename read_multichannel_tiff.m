
import ScanImageTiffReader.ScanImageTiffReader

early_abort_flag = false;

base_dir = uigetdir();
all_files = dir([base_dir,'\**\*.tif']);


for i = 1:length(all_files)
    filepath = all_files(i).folder;
    filename = all_files(i).name;
    fprintf('checking %s\n',filepath);
    
    if ~(isfolder([filepath,'\registration']) || isfolder([filepath,'\registration_001'])) || isempty(dir([filepath,'\registration*\*imagingData*.mat']))
        
        mkdir([filepath,'\registration'])
            
        %[filename,filepath] = uigetfile('*.tif');
        file = [filepath,'\',filename];
        % Get image data
        imgData = ScanImageTiffReader(file);
        
        % Get metadata
        siMetadata = scanimageMetadata(file);
        SI = parse_scanimage_metadata(siMetadata);
        
        % get metadata
        metadata = imgData.metadata();
        
        % scanimage object SI
        siMetadata = scanimageMetadata(file);
        siMetadata = parse_scanimage_metadata(siMetadata);
        SI = siMetadata.SI;
        
        % determine number of channels in data
        numChannels = length(SI.hChannels.channelSave);
        
        % get imaging data and permute so x is rows, y is columns
        imgDataSeries = permute(imgData.data(),[2,1,3]);
        
        dims = [size(imgDataSeries,1),...
                size(imgDataSeries,2),...
                SI.hStackManager.numFramesPerVolumeWithFlyback];
        
        clear img
        for i = 1:numChannels
            img{i} = imgDataSeries(:,:,i:numChannels:end);

            if early_abort_flag
                img{i} = img{i}(:,:,1:end-mod(size(img{i},3),dims(3)));
            end

            img{i} = reshape(img{i}(:,:,1:end),dims(1),dims(2),dims(3),[]);
            img{i}(:,:,SI.hStackManager.numFramesPerVolume+1:end,:) = [];
        end
        
    
        save([filepath,'\registration\imagingData_trial001.mat'],'img','-v7.3')

        imgData_reg = normcorre_regProduct(squeeze(sum(img{1},3)),false); %save motion correction of the image summed across all planes
        
        save([filepath,'\registration\imgData_reg.mat'],'imgData_reg','-v7.3')

        imgData_smooth_reg = normcorre_regProduct(smoothdata(squeeze(sum(img{1},3)),3,'movmean',5),false); %save motion correction of the image summed across all planes
        
        save([filepath,'\registration\imgData_smooth_reg.mat'],'imgData_smooth_reg','-v7.3')
    end
end