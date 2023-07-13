function regProduct = register_pablo(imgData)
n_frames = size(imgData,4);
n_planes = size(imgData,3);
sum_movie = squeeze(sum(imgData,3)); %sum over planes
ref_img   = squeeze(sum(sum_movie,3)); %sum over time
[optimizer, metric] = imregconfig('monomodal');
regProduct = nan(size(imgData));
time_vec = nan(n_frames,1);


for i = 1:n_frames
    tic
    fprintf('registering %i of %i ',i,n_frames)
    tform = imregtform(sum_movie(:,:,i),ref_img,"translation",optimizer,metric);
    regProduct(:,:,:,i) = imwarp(imgData(:,:,:,i),tform);
%     for j = 1:n_planes
%         regProduct(:,:,j,i) = imwarp(imgData(:,:,j,i),tform);
%     end
    time_vec(i) = toc;
    fprintf('time remaining %i minutes\n', round(mean(time_vec,'omitnan')*(n_frames-i)/60))
end