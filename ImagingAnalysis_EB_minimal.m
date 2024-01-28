%% clear all variables
clear all
close all

%% load in data
[filename,filepath] = uigetfile('.mat','Select Registered Movie');
[filename2,filepath2] = uigetfile(filepath,'Select FicTrac Data');
load([filepath,'\',filename])
load([filepath2,'\',filename2])
load([filepath2,'\mask.mat'])

%% define script params
ft_type= 'gaussian'; %the type of smoothing for fictrac data
ft_win = 60; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'gaussian','gaussian'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {5,5};

n_centroid = 16;
f0_pct = 7;

% Process Both
ft_data   = process_ft(ftData_DAQ, ft_win, ft_type);
im_data   = process_im_3d(regProduct, im_win, im_type, mask, n_centroid, f0_pct);

%% Plot results
bump_thresh = 10;
fly_thresh  = .2;
rho_thresh  = .2;
lag = 7;
gain = .7;
fr = mean(diff(ft_data.xf));
figure(1); clf; 

%convince yourself that the fictrac smoothing is honest. cumsum vel and cue
%pos should have a steady offset
subplot(3,1,1); hold on
a = plot(ft_data.xf,-ft_data.cue); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(ft_data.xf,mod(cumsum(ft_data.r_speed*gain*fr),2*pi)-pi); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(ft_data.xf,circ_dist(-ft_data.cue,mod(cumsum(ft_data.r_speed*gain*fr),2*pi)-pi)); a.YData(abs(diff(a.YData))>pi) = nan;
legend('cue','cumsum(vel)','diff','Location','northeastoutside')
axis tight

%convince yourself that the mu estimation is honest.
subplot(3,1,2); hold on
a = plot(ft_data.xf,-ft_data.cue); a.YData(abs(diff(a.YData))>pi) = nan;
xb = linspace(ft_data.xf(1),ft_data.xf(end),size(im_data.mu,1));
a = scatter(xb,im_data.mu,im_data.rho*5,'filled'); a.YData(abs(diff(a.YData))>pi) = nan;
legend('cue','mu','location','northeastoutside')
axis tight

%look at scatter
subplot(3,2,5); hold on
tmp = interp1(xb,unwrap(im_data.mu),ft_data.xf);
bump_vel = gradient(tmp)/fr;
bump_vel = bump_vel(lag+1:end);
fly_vel  = ft_data.r_speed(1:end-lag);
rho = interp1(xb,im_data.rho,ft_data.xf);
rho = rho(lag+1:end);

idx = bump_vel < bump_thresh & abs(fly_vel) > fly_thresh &  rho > rho_thresh;
scatter(fly_vel(idx),bump_vel(idx),'filled','MarkerFaceAlpha',5e-2)
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

text(x(2),y(1),sprintf('gain: %.2f\nr: %.2f',b(2),r),'HorizontalAlignment','left','VerticalAlignment','bottom')
xlabel('fly vel (rad/s)'); ylabel('bump vel (rad/s)')

%% play movie
pause_time = 5e-2;
figure(2); clf

subplot(2,2,1)
a = image(im_data.imgData(:,:,1)); xticks([]); yticks([]);
b = polaraxes('Position',a.Parent.Position,'color','none');
b = polarscatter(b,im_data.mu(1),1,100,'m','filled');
hold on
h = polarscatter(-ft_data.cue(1),1,100,'k','filled');
set(gca,'color','none','RTick',[],'ThetaAxisUnits','radians','ThetaLim',[-pi,pi],'RLim',[0,1],'ThetaColor','m')

subplot(2,2,2)
c = plot(linspace(-pi,pi,size(im_data.z,1)),im_data.z(:,1));
hold on
f = scatter(pi,0,'m','filled');
g = scatter(pi,0,'k','filled');
ylim([min(im_data.z,[],'all'),max(im_data.z,[],'all')])
legend('dff (z)','mu','cue')

subplot(2,1,2)
xb = linspace(min(ft_data.xf),max(ft_data.xf),size(im_data.imgData,3));
cue = interp1(ft_data.xf,ft_data.cue,xb);
d = plot(xb,nan(size(cue)),'k');
hold on
e = plot(xb,nan(size(cue)),'m');
xlim([0,max(xb)])

for i = 1:size(im_data.imgData,3)
    a.CData = im_data.imgData(:,:,i);
    b.ThetaData = im_data.mu(i);
    c.YData = im_data.z(:,i);
    f.XData = im_data.mu(i);
    g.XData = -cue(i);
    h.ThetaData = -cue(i);
    if abs(diff(cue(i:i+1))) < pi
    d.YData(i) = -cue(i);
    end
    if abs(diff(im_data.mu(i:i+1))) < pi
    e.YData(i) = im_data.mu(i);
    end
    pause(pause_time)
end
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
    s.imgData = imgData;
end