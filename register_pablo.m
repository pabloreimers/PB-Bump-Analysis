function regProduct = register_pablo(imgData)
n_frames = size(imgData,4);
n_planes = size(imgData,3);
sum_movie = squeeze(sum(imgData,3)); %sum over planes
ref_img   = squeeze(sum(sum_movie,3)); %sum over time
[optimizer, metric] = imregconfig('monomodal');
regProduct = nan(size(imgData));


for i = 1:n_frames
    fprintf('registering %i of %i\n',i,n_frames)
    tform = imregtform(sum_movie(:,:,i),ref_img,"translation",optimizer,metric);
    regProduct(:,:,:,i) = imwarp(imgData,tform)
%     for j = 1:n_planes
%         regProduct(:,:,j,i) = imwarp(imgData(:,:,j,i),tform);
%     end
end