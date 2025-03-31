base_dir = uigetdir();
all_files = dir([base_dir,'\**\imagingData*.mat']);

bands = (1:50); 

tic
for j = 1:length(all_files)
    filepath = all_files(j).folder;
    filename = all_files(j).name;
    fprintf('checking %s \n ',filepath);

    load([filepath,'\',filename])

    dims = size(img{1});
    f = img{1};
    
    fhat = fft(f,dims(2),2);
    fhat(:,bands,:,:) = 0;
    ffilt = abs(ifft(fhat,dims(2),2));
    
    ffilt = reshape(ffilt,dims(1),dims(2),[]);
    for i = 1:length(ffilt)
        ffilt(:,:,i) = medfilt2(ffilt(:,:,i));
    end
    imgData = uint16(reshape(ffilt,dims));
    
    save([filepath,'\imgData_denoise.mat'],'imgData','-v7.3')

    imgData = normcorre_regProduct(squeeze(sum(imgData,3)),false);
    
    save([filepath,'\imgData_denoise_reg.mat'],'imgData','-v7.3')

    fprintf('\nETR: %.2f mins\n',(toc/j) * (length(all_files) - j) / 60)
end
    