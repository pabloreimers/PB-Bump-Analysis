%% clear all
clear all
close all

%% find path to all relevant files
base_dir = ('Z:pablo\stacks\lpsp_kir_redo\');
all_files = dir([base_dir,'\**\*imagingData*.mat']);

%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        imgData = mean(regProduct,[3,4]);
        figure(1); clf; imagesc(imgData); colormap(bone); drawnow;
        tmp = drawellipse('Center',[size(imgData,2)/2,size(imgData,1)/2],'SemiAxes',[size(imgData,2)/4,size(imgData,1)/4]);
        input('') %move on once user presses enter, after adjusting the ellipse
        mask = createMask(tmp);
        save([fileparts(all_files(i).folder),'\mask.mat'],'mask')
    end
end

%% process and store all values
ft_type= 'gaussian'; %the type of smoothing for fictrac data
ft_win = 60; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'gaussian','gaussian'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {5,5};

n_centroid = 16;
f0_pct = 7;

all_data = struct();

tic
for i = 1:length(all_files)
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{7})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])

    all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
    all_data(i).im = process_im_3d(regProduct, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% plot all
vel_thresh = .2;
bump_thresh = 10;
rho_thresh = .2;
lag = 7;

cc       = nan(length(all_data),1);
gains    = nan(length(all_data),1);
lpsp_idx = false(length(all_data),1);

figure(1); clf
for i = 1:length(all_data)
    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = all_data(i).ft.r_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);

    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;
    subplot(9,9,i); hold on
    scatter(fly_vel(idx),bump_vel(idx),'filled','markerfacealpha',.1)
    y = ylim; x = xlim;
    plot(x,[0,0],':k'); 
    plot([0,0],y,':k');
    x = xlim;
    y = ylim;
    b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); %fit the slope of fly vel and bump vel with an arbitrary offset
    r = corr(fly_vel(idx),bump_vel(idx));
    plot([0,0],y,'k:')
    plot(x,[0,0],':k')
    plot(x,x*b(2) + b(1),'r')
    axis tight equal    
    text(x(2),y(1),sprintf('gain: %.2f\nr: %.2f',b(2),r),'HorizontalAlignment','right','VerticalAlignment','bottom')
    xlabel('fly vel (rad/s)'); ylabel('bump vel (rad/s)')

        cc(i) = r;
    gains(i) = b(2);
    lpsp_idx(i) = contains(all_data(i).meta,'_lpsp_','IgnoreCase',true);
    tmp = strsplit(all_data(i).meta,'\');
    if lpsp_idx(i)
        title(strcat(tmp(5),tmp(6)),'Color','r')
    else
        title(strcat(tmp(5),tmp(6)),'Color','k')
    end

end

%% test significance (bootstrap to ask about a mean difference, we don't know how variances compare)
N = 1e4; %number of bootstrap reps


cc_lpsp = cc(lpsp_idx);
cc_empty= cc(~lpsp_idx);

i_lpsp = randi(length(cc_lpsp),N);
i_empty = randi(length(cc_empty),N);

p = sum(mean(cc_lpsp(i_lpsp),1) - mean(cc_empty(i_empty),1) > 0) / N

%% Functions

function s = process_ft(ftData_DAQ, ft_win, ft_type)

    f_speed = ftData_DAQ.velFor{:};                       %store each speed
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}';
    cue     = smoothdata(unwrap(cue / 192 *2 * pi - pi),1,ft_type,ft_win);
    cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
    cue(cue > pi) = cue(cue > pi) - 2*pi;
    
    s.xf      = seconds(ftData_DAQ.trialTime{:});
    s.f_speed = smoothdata(f_speed,1,ft_type,ft_win); 
    s.r_speed = smoothdata(r_speed,1,ft_type,ft_win);
    s.cue     = cue;
    
end

function s = process_im(regProduct, im_win, im_type, mask, n_centroid, f0_pct)
    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
    mid             = bwmorph(mask,'remove');
    [y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
    ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
    [y0,x0]         = find(ep,1);
    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
    xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
    centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
    if centroids(1,1) < centroids(end,1)
        centroids = flipud(centroids);
    end
    [~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)

    imgData = squeeze(sum(regProduct,3)); %sum across all z slices
    imgData = imgData - prctile(imgData,1,'all'); %baseline subtract

    imgData = smoothdata(imgData,3,im_type{1},im_win{1});

    imgData_2d = reshape(imgData,[],size(imgData,3)); %reshape the data to be all pixels x all time points
    centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
    end
    f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;

    s.mu = mu;
    s.rho= smoothdata(rho,1,im_type{2},im_win{2});
    s.z  = zscore_cluster;
    s.alpha = alpha;
    s.imgData = imgData;
end

function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function s = process_im_3d(regProduct, im_win, im_type, mask, n_centroid, f0_pct)
    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
    mid             = bwmorph(mask,'remove');
    [y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
    xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
    centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
    if centroids(1,1) < centroids(end,1)
        centroids = flipud(centroids);
    end

    centroids3 = [centroids,(1-centroids(:,1)/max(centroids(:,1)))*size(regProduct,3)*4];


    %[~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)
    y_mask3 = repmat(y_mask,size(regProduct,3),1);
    x_mask3 = repmat(x_mask,size(regProduct,3),1);
    z_mask3 = reshape(ones(size(y_mask,1),size(regProduct,3)).*[1:size(regProduct,3)],size(y_mask3,1),1);
    [~,idx] = pdist2(centroids3,[y_mask3,x_mask3,z_mask3.*4],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)
    
    imgData = int16(smoothdata(regProduct,4,im_type{1},im_win{1}));
    imgData = imgData - min(imgData,[],'all');
    
    imgData_2d = reshape(imgData,[],size(imgData,4)); %reshape the data to be all pixels x all time points
    centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask3(idx==i),x_mask3(idx ==i),z_mask3(idx ==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});
    
    imgData = squeeze(sum(imgData,3));
    imgData = 256*(imgData-min(imgData,[],'all'))/(max(imgData,[],'all')-min(imgData,[],'all'));

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.alpha = alpha;
    %s.imgData = imgData;
end