%% reset the workspace and set script parameters
clear all
close all
movie_flag          = false;                                                %set whether to play movies or just run straight through the script
pause_time          = 0;                                                    %set the pause time if playing movies
data_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %set main data directory for ease of selecting files
%data_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\2p data\';
f0_pct              = 15;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 5;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 5;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = 5e-3;
%% ask user for image data
[filename,filepath] = uigetfile('.mat','Select Registered Movie',data_dir);
load([filepath,'\',filename])

%% plot the movie of all planes separately
pause_time  = 0;                         %pause this many seconds between each frame
n_plane = size(regProduct,3);
n_frame = size(regProduct,4);
row    = ceil(sqrt(n_plane));
col     = ceil(sqrt(n_plane));

figure(1); clf
for i = 1:n_plane
    h(i) = imagesc(subplot(row,col,i),regProduct(:,:,i,1));
    axis equal tight
    xticks([])
    yticks([])
    title(['plane: ',num2str(i)])
    colormap(bone)
end

if movie_flag
for f = 1:size(regProduct,4)
    for i = 1:size(regProduct,3)
        h(i).CData = regProduct(:,:,i,f);
    end
    pause(pause_time)
end
end

%% Make compress in z, and scale between 0 and 256
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)

imgData2 = imgData;
top_int = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
bot_int = prctile(imgData2,5,'all');
imgData2 = max(min(imgData2,top_int),bot_int) - bot_int;
imgData2 = 256*imgData2/max(imgData2,[],'all');
imgData2 = smoothdata(imgData2,3,'movmean',avg_win);                      %smooth the data, again for viewing purposes (should this go before the clipping)            

%% Play the movie
pause_time  = 0;                         %pause this many seconds between each frame

figure(1); clf                              %clear the current figure
h           = image(imgData(:,:,1));        %initialize an image object where we will update the color values in each frame 
axis equal tight                            %make all pixels square, and keep axis tight for clean look
colormap(bone)                              %set to favorite colormap. I like black and white
colorbar

if movie_flag
for i = 1:size(imgData,3)                   %loop through each frame and update the ColorData with the values for the current frame
    h.CData = imgData(:,:,i);
    pause(pause_time)
end
end

%% plot the pb and draw bounding box over the whole thing
figure(1); clf                                          % clear the current figure
mask = roipoly(uint8(mean(imgData2,3)));               %this create a black and white mask (logical) of the PB, drawn on a maxZ projection of the image
figure(2)
imagesc(subplot(2,1,1),mask); colormap(bone); xticks([]); yticks([])

%% extract midline axis, and subsample for centroid locations (using graph search)
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length
min_axis        = min(range(x_mask),range(y_mask));
mid             = bwskel(mask,'MinBranchLength',min_axis);  %find the midline as the skeleton, shaving out all sub branches that are smaller than the minimum axis length
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);
figure(2); imagesc(subplot(2,1,2),mid); colormap(bone); xticks([]); yticks([])
hold on; scatter(x0,y0,'b','filled')

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!

xq          = [-min_axis:(length(x_mid)+min_axis)];                             %extend the midline so that it reaches the border of the mask. extrapolate as many points as the minimum axis length
x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');                 %keep only the points that exist within the mask
x_mid       = x_mid(idx);
y_mid       = y_mid(idx);

xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';                        %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
cmap        = repmat(make_colors(n_centroid),2,1); %repmat(cbrewer2('set1',n_centroid),2,1);                     %set a colormap to plot each centroid in a different color, and which repeats per hemisphere (note that if a hemisphere has more clusters than colorbrewer can generate, it will repeat colors within each hemisphere).

figure(1); clf
imagesc(max(imgData,[],3))                                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%% assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)
figure(1); clf
imagesc(max(imgData,[],3))                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_mask(idx == i),y_mask(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%% Find the mean activity in each group over time
imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
zscore_cluster  = zscore(f_cluster,[],2);
%dff_cluster     = smoothdata(dff_cluster,2,'gaussian',b_smooth);
figure(3); clf
imagesc(subplot(2,2,1),centroid_log); colormap('bone')
xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
imagesc(subplot(2,2,2),imgData_2d); colormap('bone')
xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
imagesc(subplot(2,1,2),dff_cluster); colormap('bone'); pos = get(gca,'position'); colorbar; set(gca,'Position',pos)
xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')
hold on
scatter(ones(1,2*n_centroid),1:2*n_centroid,100,cmap,'filled','square')

%% Estimate von mises parameters (both sides)
% n_frames    = size(dff_cluster,2);
% alpha       = linspace(-pi,pi,n_centroid);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
% muL         = nan(n_frames,1);                      %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
% muR         = nan(n_frames,1);                      %mu is thetahat
% kappaL      = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
% kappaR      = nan(n_frames,1);
% ampL        = nan(n_frames,1);                      %this will be used to scale the bump amplitude from a generic von mises pdf (integral 1) to something that matches the data
% ampR        = nan(n_frames,1);
% r2L         = nan(n_frames,1);                      %rsquared values, telling how much of the variance in the distribution is explained by our fit (as compared to a flat line which is the mean)
% r2R         = nan(n_frames,1);
% 
% for i = 1:n_frames                                                                      %for each frame
%     [muL(i), kappaL(i)] = circ_vmpar(alpha,dff_cluster(1:n_centroid,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
%     [muR(i), kappaR(i)] = circ_vmpar(alpha,dff_cluster((1:n_centroid) + n_centroid,i)); %do this for each hemisphere separately.
%     if any(isnan(circ_vmpdf(alpha,muL(i),kappaL(i)))) || any(isnan(circ_vmpdf(alpha,muR(i),kappaR(i)))) %if the distribution returns nans (calculated something too sharp, e.g.)
%         muL(i) = nan;
%         muR(i) = nan;
%         kappaL(i) = nan;
%         kappaR(i) = nan;
%         continue
%     end
%     [ampL(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muL(i),kappaL(i)) - dff_cluster(1:n_centroid,i)).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
%     sst                 = sum( (ampL(i)*circ_vmpdf(alpha,muL(i),kappaL(i)) - mean(dff_cluster(1:n_centroid,i)) ).^2);        %find the squared error with the mean
%     r2L(i)              = 1 - ssr/sst;                                                                                       %calculate r2 (which is 1- ratio of our fit/baseline fit).
%     [ampR(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muR(i),kappaR(i)) - dff_cluster((1:n_centroid) + n_centroid,i)).^2),1);
%     sst                 = sum( (ampR(i)*circ_vmpdf(alpha,muR(i),kappaR(i)) - mean(dff_cluster((1:n_centroid)+n_centroid,i)) ).^2);
%     r2R(i)              = 1 - ssr/sst;
% end
% 
% idxL = r2L < 0;
% idxR = r2R < 0;
% 
% %muL(idxL) = 0;
% kappaL(idxL) = 0;
% %muR(idxR) = 0;
% kappaR(idxR) = 0;

%% Estimate von mises (Mel Style)
% n_frames    = size(dff_cluster,2);
% alpha       = linspace(-pi,pi,n_centroid);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
% adj_rs      = nan(n_frames,1);
% mu          = nan(n_frames,1);
% amp         = nan(n_frames,1);
% width       = nan(n_frames,1);
% 
% fo = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[0,0,0,-pi],... % [a,c,k,u]
%     'Upper',[inf,inf,inf,pi],...
%     'StartPoint',[1,0,1,0]);
% ft = fittype('a*exp(k*cos(x-u))+c','options',fo);
% 
% for i = 1:n_frames
%     [f, gof] = fit([alpha,alpha]',zscore_cluster(:,i),ft,...
%         'MaxIter',20000,'MaxFunEvals',20000);
%     adj_rs(i) = gof.adjrsquare;
%     mu(i)        = f.u;
%     amp(i)       = f.a * ( exp(f.k) - exp(-f.k) );
%     width(i)     = 2 * abs( acos( 1/f.k * log( 1/2 *( exp(f.k) + exp(-f.k) ))));
%     fprintf('frame %i / %i\n',i,n_frame)
% end

%% estimate von mises for joint left and right
% n_frames    = size(dff_cluster,2);
% alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
% mu          = nan(n_frames,1);                      %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
% kappa       = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
% amp         = nan(n_frames,1);                      %this will be used to scale the bump amplitude from a generic von mises pdf (integral 1) to something that matches the data
% r2          = nan(n_frames,1);                      %rsquared values, telling how much of the variance in the distribution is explained by our fit (as compared to a flat line which is the mean)
% 
% idx = var(dff_cluster,[],2) > var_thresh;
% 
% for i = 1:n_frames                                                                      %for each frame
%     [mu(i), kappa(i)] = circ_vmpar(alpha(idx),dff_cluster(idx,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
%     if any(isnan(circ_vmpdf(alpha,mu(i),kappa(i))))%if the distribution returns nans (calculated something too sharp, e.g.)
%         mu(i) = nan;
%         kappa(i) = nan;
%         continue
%     end
%     [amp(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha(idx),mu(i),kappa(i)) - dff_cluster(idx,i)).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
%     sst                 = sum( ( amp(i)*circ_vmpdf(alpha(idx),mu(i),kappa(i))   - mean(dff_cluster(idx,i))).^2);        %find the squared error with the mean
%     r2(i)              = 1 - ssr/sst;
% end

%% estimate von mises for joint left and right (zscore data, fit c)
% n_frames    = size(dff_cluster,2);
% alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
% mu          = nan(n_frames,1);                      %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
% kappa       = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
% amp         = nan(n_frames,1);                      %this will be used to scale the bump amplitude from a generic von mises pdf (integral 1) to something that matches the data
% c           = nan(n_frames,1);
% r2          = nan(n_frames,1);                      %rsquared values, telling how much of the variance in the distribution is explained by our fit (as compared to a flat line which is the mean)
% 
% for i = 1:n_frames                                                                      %for each frame
%     [mu(i), kappa(i)] = circ_vmpar(alpha,zscore_cluster(:,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
%     if any(isnan(circ_vmpdf(alpha,mu(i),kappa(i))))%if the distribution returns nans (calculated something too sharp, e.g.)
%         mu(i) = nan;
%         kappa(i) = nan;
%         continue
%     end
%     [tmp,ssr]       = fminsearch(@(a) sum((a(1)*circ_vmpdf(alpha,mu(i),kappa(i))+a(2) - zscore_cluster(:,i)).^2),[1,0]);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
%     amp(i)          = tmp(1);
%     c(i)            = tmp(2);
%     sst             = sum( ( amp(i)*circ_vmpdf(alpha,mu(i),kappa(i))+c(i)   - mean(zscore_cluster(:,i))).^2);        %find the squared error with the mean
%     r2(i)              = 1 - ssr/sst;
% end
% 
% c(isoutlier(c,'percentiles',[1,99])) = nan;
% amp(isoutlier(amp,'percentiles',[1,99])) = nan;
% 
% for i = 1:n_smooth
% mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
% amp    = smoothdata(amp,1,'gaussian',b_smooth);
% kappa  = smoothdata(kappa,1,'gaussian',b_smooth);
% end
% 
% mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
% mu(mu > pi) = mu(mu > pi) - 2*pi;

% %% estimate von mises for joint left and right (dff data, c = mean)
% tmp_data    = dff_cluster;
% n_frames    = size(dff_cluster,2);
% alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
% mu          = nan(n_frames,1);                      %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
% kappa       = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
% amp         = nan(n_frames,1);                      %this will be used to scale the bump amplitude from a generic von mises pdf (integral 1) to something that matches the data
% c           = nan(n_frames,1);
% r2          = nan(n_frames,1);                      %rsquared values, telling how much of the variance in the distribution is explained by our fit (as compared to a flat line which is the mean)
% 
% for i = 1:n_frames                                                                      %for each frame
%     c(i) = mean(tmp_data(:,1));
%     [mu(i), kappa(i)] = circ_vmpar(alpha,tmp_data(:,i) - c(i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
%     if any(isnan(circ_vmpdf(alpha,mu(i),kappa(i))))%if the distribution returns nans (calculated something too sharp, e.g.)
%         mu(i) = nan;
%         kappa(i) = nan;
%         continue
%     end
%     [amp(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,mu(i),kappa(i)) - (tmp_data(:,i)-c(i))).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
%     sst             = sum( ( amp(i)*circ_vmpdf(alpha,mu(i),kappa(i))+c(i)   - mean(tmp_data(:,i))).^2);        %find the squared error with the mean
%     r2(i)              = 1 - ssr/sst;
% end
% 
% c(isoutlier(c,'percentiles',[1,99])) = nan;
% amp(isoutlier(amp,'percentiles',[1,99])) = nan;
% 
% for i = 1:n_smooth
% mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
% amp    = smoothdata(amp,1,'gaussian',b_smooth);
% kappa  = smoothdata(kappa,1,'gaussian',b_smooth);
% end
% 
% mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
% mu(mu > pi) = mu(mu > pi) - 2*pi;

%% Find bump as PVA
tmp_data    = dff_cluster;
alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);

[x_tmp,y_tmp]   = pol2cart(alpha,tmp_data');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
for i = 1:n_smooth
mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
rho    = smoothdata(rho,1,'gaussian',b_smooth);
end

mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
amp = mean(dff_cluster,1)';
[~,i] = mink(abs(alpha-mu),2,2); %find indexes corresponding to peak of each time point
i2 = i' + size(dff_cluster,1)*[0:size(dff_cluster,2)-1]; %find the linear index into the peak of each column (time point) value. this was clever :)
amp_peak = mean(dff_cluster(i2),1)'; %extract peak amplitude

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
amp = smoothdata(amp,1,'gaussian',b_smooth); 
amp_peak = smoothdata(amp_peak,1,'gaussian',b_smooth);
end

%% plot! movie
% c1 = [1,0.5,0];                                                             %define the colors for the left and right bump
% c2 = [0,0.5,1];
% n_frames = size(dff_cluster,2);
% n_ticks  = 4;
% f_data = zscore_cluster;
% 
% figure(4); clf
% subplot(2,2,3)
% hold on
% clear h v p
% h(1) = plot(1:(2*n_centroid),nan(2*n_centroid,1),'k.-');                    %initialize a plot to show the average activity at each centroid over time
% v(1)= plot(1:2*n_centroid,nan(2*n_centroid,1),'Color',c1);
% pos = get(gca,'Position');
% legend('Activity','Fit','autoupdate','off','Location','northoutside')
% tick_vec  = 1:2*n_centroid;
% label_vec = [alpha,alpha];
% idx       = linspace(0,2*n_centroid,n_ticks+1);
% idx(1)    = [];
% set(gca,'Position',pos, 'YLim',[min(f_data,[],'all'),max(f_data,[],'all')],...
%    'XLim',[0,2*n_centroid],'XTick',tick_vec(idx),'XTickLabels',label_vec(idx)/pi)
% xlabel('Cluster (\pi)')
% ylabel('Activity (a.u.)')
% 
% subplot(2,2,1)                                                                      %initialize another plot to show the movie as the extracted data plays
% h(2) = image(imgData2(:,:,1));
% title('processed F')
% hold on
% [y,x] = find(bwmorph(mask,'remove'));                                               %overlay the outline of our mask
% [y,x] = graph_sort(y,x);
% plot(x,y,':','Color',[0.75,0.75,0.75])
% axis equal tight
% colormap(bone)
% xticks([])
% yticks([])
% 
% bump_params = {amp;mu;kappa};                                                    %store all of the bump statistics in a cell array for easy repetitive access in a loop
% bump_params = cellfun(@(x)(smoothdata(x,1,'gaussian',b_smooth)),bump_params,'UniformOutput',false); %smooth all of the bump parameters.
% bp_names    = {'Amplitude','Position (\mu)','Concentration (\kappa)'};
% bp_h        = cell(3,1);
% 
% for i = 1:size(bump_params,1)                                                       %initialize plots for each bump statistic over each frame
%     set(subplot(3,2,2*i),'NextPlot','add','xlim',[0,n_frames],...
%         'ylim',[min([bump_params{i,:}],[],'all'),max([bump_params{i,:}],[],'all')])
%     ylabel(subplot(3,2,2*i),bp_names{i})
%     bp_h{i,1} = scatter(subplot(3,2,2*i),1:n_frames,nan(n_frames,1),'filled','MarkerFaceColor',c1);
% end
% xlabel('frame')
% set(subplot(3,2,4),'YTick',[-pi,0,pi],'YTickLabels',{'-\pi','0','\pi'},'YLim',[-pi,pi])
% title(subplot(3,2,2),'Bump Parameters')
% 
% if movie_flag
% for i = 1:n_frames                                       %at each frame, update:
%     h(1).YData = f_data(:,i);                                    %the average intensity per centroids
%     h(2).CData = imgData2(:,:,i);                                        %the displayed image from the movie
% 
%     v(1).YData = amp(i)*circ_vmpdf(alpha,mu(i),kappa(i)) + c(i);            %generate the fit von mises distribution, over the range of thetas, where it is centered at this time point, and how concentrated. scale it, too
%     
%     for j = 1:size(bump_params,1)                                       %update the current distribution statistics, too
%         bp_h{j,1}.YData(i) = bump_params{j,1}(i);
%     end
%     
%     if i > 1                                                            %because circular, if the instantaneous position jump is more than pi radians, just blank the point (because the shortest path is actually out of frame)
%     if abs(bump_params{2,1}(i) - bump_params{2,1}(i-1)) > pi || r2(i) < vm_thresh
%         bp_h{2,1}.YData(i) = nan;
%     end
%     end
%     
%     pause(pause_time)
% end
% else
%     for j = 1:size(bump_params,1)                                       %update the current distribution statistics, too
%         bp_h{j,1}.YData = bump_params{j,1};
%         bp_h{j,1}.YData(r2 < vm_thresh) = nan;
%     end
% end

%% load fictrac data
[filename2,filepath2] = uigetfile(filepath,'Select FicTrac Data');
load([filepath2,'\',filename2])

%% Process fictrac data
if ~exist('ftData_DAQ','var')
    ftData_DAQ= ftData;
    ftData_DAQ.velFor = ftData.fwSpeed;
    ftData_DAQ.velYaw = ftData.yawSpeed;
end

% f_speed = ftData.fwSpeed{:};                       %store each speed
% r_speed = ftData.yawSpeed{:};
% intHD   = ftData.intHD{:};

f_speed = ftData_DAQ.velFor{:};                       %store each speed
r_speed = ftData_DAQ.velYaw{:};
intHD   = ftData_DAQ.intHD{:};
cue     = ftData_DAQ.cuePos{:}';


f_speed = f_speed;                                      %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi
cue     = unwrap(cue / 192 * 2*pi - pi);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'gaussian',f_smooth); 
r_speed = smoothdata(r_speed,1,'gaussian',f_smooth);
intHD   = smoothdata(intHD,  1,'gaussian',f_smooth);
cue     = smoothdata(cue,    1,'gaussian',f_smooth);
end

intHD = mod(intHD,2*pi);                                %rewrap heading data, and put between -pi and pi.
intHD(intHD > pi) = intHD(intHD > pi) - 2*pi;
cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;

xf  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed)]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
if isduration(xf)
    xf = seconds(xf);
end
total_t = max(xf);
xb  = linspace(0,total_t,size(imgData,3))';
fr  = mean(diff(xb));

% f_speed = interp1(xf,f_speed,xb,[],'extrap');                       %interpolate the fictrac trace to the timepoints of the image (bump)
% r_speed = interp1(xf,r_speed,xb,[],'extrap');
% intHD   = interp1(xf,intHD,xb,[],'extrap');
% cue     = interp1(xf,cue,xb,[],'extrap');

mu          = interp1(xb,unwrap(mu),xf)';
mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp         = interp1(xb,amp,xf)';
amp_peak         = interp1(xb,amp_peak,xf)';

%% Plot bump params over fictrac data
% figure(5); clf
% 
% h(1) = subplot(3,1,1);
% tmp = -intHD; tmp(abs(diff(tmp)) > pi) = nan;
% plot(xb,-tmp,'k'); ylabel('Heading (rad)'); ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'});
% yyaxis right; hold on
% tmp = mu; tmp(abs(diff(tmp)) > pi) = nan; tmp(r2 < vm_thresh) = nan;
% scatter(xb,tmp,'.','MarkerEdgeColor',c1)
% ylabel('Bump Position (\mu)'); yticks([])
% axis tight
% 
% h(2) = subplot(3,1,2);
% plot(xb,f_speed,'k'); ylabel('Forward Speed (mm/s)');
% yyaxis right; hold on
% tmp = amp; tmp(r2 < vm_thresh) = nan;
% scatter(xb,tmp,'.','MarkerEdgeColor',c1)
% ylabel('Bump Amplitude'); yticks([]); axis tight
% 
% h(3) = subplot(3,1,3);
% plot(xb,r_speed,'k'); ylabel('Rotational Speed (rad/s)');
% yyaxis right; hold on
% tmp = amp; tmp(r2 < vm_thresh) = nan;
% scatter(xb,tmp,'.','MarkerEdgeColor',c1)
% ylabel('Bump Amplitude'); yticks([]); axis tight
% xlabel('time (s)')
% 
% linkaxes(h,'x')

%% linear model fictrac params with dF/F
% if cluster_flag
%    dff     = smoothdata(mean(dff_cluster,1)',1,'gaussian',b_smooth);
% else
%     mask_1d = reshape(mask,[],1);                           % reshape the mask into a single vector
%     f       = reshape(imgData,size(mask_1d,1),[]);          % reshape the movie into pixels x frames
%     f       = f(mask_1d,:);                                 % extract only pixels within the mask
%     f0      = prctile(f,f0_pct,2);                    % find the baseline fluorescence in each pixel over time
%     dff     = (f - f0) ./ f0;                               % calculate pixelwise df/f
%     dff     = mean(dff,1)';
%     dff     = smoothdata(dff,1,'gaussian',b_smooth);        % smooth dff
% end
% amp     = dff;

c1 = [1,0.5,0];                                                             %define the colors for the left and right bump
c2 = [0,0.5,1];
fr      = mean(diff(xf)); %find the frame rate of the data

[f_corr,f_lag]  = xcorr(f_speed,amp,ceil(max_lag/fr),'coeff');      % find the max correlation, in steps
[r_corr,r_lag]  = xcorr(r_speed,amp,ceil(max_lag/fr),'coeff');
[~,l_idx]       = max(f_corr);
f_lag           = f_lag(l_idx);
[~,l_idx]       = max(r_corr);
r_lag           = r_lag(l_idx);

f_idx           = xf - f_lag*fr > 0 & xf - f_lag*fr < max(xf); %find the values that are not within the lag (in time)
r_idx           = xf - r_lag*fr > 0 & xf - r_lag*fr < max(xf);

f_speed(~f_idx) = 0;                                            %set those to zero before you shift things
r_speed(~r_idx) = 0;
 
f_speed_lag     = circshift(f_speed,-f_lag);                    %shift the traces to be aligned.
r_speed_lag     = circshift(r_speed,-r_lag);

amp(isnan(amp)) = 0;

[f_fit,f_gof]   = fit(f_speed_lag,amp,'poly1');
[r_fit,r_gof]   = fit(r_speed_lag,amp,'poly1');
[m_fit,m_gof]   = fit([r_speed_lag,f_speed_lag],amp,'poly11');

figure(6);clf
subplot(2,1,1)
plot(xf,f_speed_lag,'k'); ylabel('Forward Velocity (mm/s)')
yyaxis right; plot(xf, amp,'Color',c1); ylabel('Average \DeltaF/F'); ax = gca; ax.YAxis(2).Color = c1;
axis tight
% y = ylim; text(xb(end),y(2),sprintf('r^2: %.2f\n2p lag: %.1fs\nball: %.1fs\n2p: %.1fs',...
%                 f_gof.adjrsquare,f_lag*fr,f_smooth*fr,b_smooth*fr),...
%                 'HorizontalAlignment','right','VerticalAlignment','top')
y = ylim; text(xf(end),y(2),sprintf('r^2: %.2f\n',...
                f_gof.adjrsquare),...
                'HorizontalAlignment','right','VerticalAlignment','top')

subplot(2,1,2)
plot(xf,r_speed_lag,'k'); ylabel('Rotational Speed (rad/s)')
yyaxis right; plot(xf, amp,'Color',c1); ylabel('Average \DeltaF/F'); ax = gca; ax.YAxis(2).Color = c1;
xlabel('time (s)')
axis tight
% y = ylim; text(xb(end),y(2),sprintf('r^2: %.2f\n2p lag: %.1fs\njoint r^2:%.2f',...
%                                 r_gof.adjrsquare,r_lag*fr,m_gof.adjrsquare),...
%                 'HorizontalAlignment','right','VerticalAlignment','top')
y = ylim; text(xf(end),y(2),sprintf('r^2: %.2f\njoint r^2: %.2f',...
                                r_gof.adjrsquare,m_gof.adjrsquare),...
                'HorizontalAlignment','right','VerticalAlignment','top')

%% Plot heading with fluorescence data
tmp_data = dff_cluster;

figure(7); clf; clear h
h(1) = subplot(3,1,1);
imagesc(xb,1:2*n_centroid,tmp_data); xticks([]); ylabel('cluster'); title('\DeltaF/F0'); pos = get(gca,'Position'); colorbar; set(gca,'Position',pos); colormap(bone)
h(2) = subplot(3,1,2);
tmp = tmp_data ./ sum(tmp_data,1);
top_val = prctile(tmp(:),95);
bot_val = prctile(tmp(:),5);
tmp( tmp > top_val) = top_val;
tmp( tmp < bot_val) = bot_val;
imagesc(xb,1:2*n_centroid,tmp); xticks([]); ylabel('cluster'); title('\DeltaF/(F0*sum(F_t))'); pos = get(gca,'Position'); colorbar; set(gca,'Position',pos); colormap(bone)
h(3) = subplot(3,1,3); 
tmp = cue; tmp(abs(diff(tmp)) > pi) = nan; plot(xf,tmp);
axis tight; xlabel('time(s)'); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'}); ylim([-pi,pi]); ylabel('heading')

linkaxes(h,'x')

%% Scatter plot of bump to fly velocity
% bump_vel = [diff(mu);0] / fr; %divide by the frame rate (seconds per frame). since the mu is in units of radians, this should give rad per seconds
% fly_vel  = ftData_DAQ.velYaw{:};
% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
% fly_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
% fly_vel = interp1(xf,fly_vel,xb,[],'extrap');
% vel_idx = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min; %ignore outlier bump speeds with arbitrary threshold
% 
% 
% figure(9); clf
% scatter(bump_vel(r2 > vm_thresh & vel_idx),fly_vel(r2 > vm_thresh & vel_idx),'filled')
% xlabel('Bump Vel (rad/s)'); ylabel('Fly Rot Vel (rad/s)');
% hold on
% [rho,pval] = corr(bump_vel(r2 > vm_thresh & vel_idx),fly_vel(r2 > vm_thresh & vel_idx));
% x = xlim;
% y = ylim;
% text(x(2),y(2),sprintf('$m = %.2f$\n$p < 10^{%i}$',rho,ceil(log10(pval))),'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
% plot([0,0],y,'k:')
% plot(x,[0,0],':k')


%% plot fly heading and bump heading, from PVA
%lag = ceil(fmincon(@(x)(circ_corrcc(cue(1:end-ceil(x)),-mu((ceil(x)+1):end))),20,[],[],[],[],0));
lag = 20;
rho_idx = rho>rho_thresh;
[pva_corr,pva_pval] = circ_corrcc(cue(1:end-ceil(lag)),-mu((ceil(lag)+1):end));

figure(7); clf
subplot(3,1,1)
a = plot(xf,cue,'k','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'}); ylim([-pi,pi])
ylabel('Fly Heading (cue)')

subplot(3,1,2)
scatter(xf(rho_idx),-mu(rho_idx),[],rho(rho_idx),'.')
colormap('bone')
yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'}); ylim([-pi,pi])
pos = get(gca,'Position');
colorbar
set(gca,'Position',pos)
ylabel('PVA')

subplot(3,1,3)
tmp = circ_dist(cue,-mu);
a = plot(xf(rho_idx),tmp(rho_idx),'k','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
ylabel('Offset')
yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})

linkaxes(get(gcf,'Children'),'x')
axis tight; ylim([-pi,pi])
xlabel('time (s)')

%% plot scatter of fly velocity to bump velocity
bump_vel = [diff(unwrap(mu));0] / fr; %divide by the frame rate (seconds per frame). since the mu is in units of radians, this should give rad per seconds
fly_vel  = ftData_DAQ.velYaw{:};
bump_vel = bump_vel(lag+1:end);
fly_vel  = fly_vel(1:end-lag);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
fly_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
end
vel_idx = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min; %ignore outlier bump speeds with arbitrary threshold

figure(9); clf
scatter(bump_vel(vel_idx),fly_vel(vel_idx),5,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.1)
xlabel('Bump Vel (rad/s)'); ylabel('Fly Rot Vel (rad/s)');
hold on
[vel_rho,vel_pval] = corr(bump_vel(vel_idx),fly_vel(vel_idx));
x = xlim;
y = ylim;
text(x(2),y(2),sprintf('$r = %.2f$\n$p < 10^{%i}$',vel_rho,ceil(log10(vel_pval))),'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
plot([0,0],y,'k:')
plot(x,[0,0],':k')

lag = lag*fr;

%% plot bump gains
fly_vel  = ftData_DAQ.velYaw{:};
for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
fly_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
end

cue_gain = median([diff(unwrap(-cue))/fr;0] ./ fly_vel);
mu_gain  = median([diff(unwrap(mu))/fr;0] ./ fly_vel); 

figure(8); clf
hold on
plot(xf,unwrap(cue),'k','linewidth',2)
plot(xf, -unwrap(mu),'b','linewidth',2)
legend(sprintf('cue, gain = %.2f',cue_gain),sprintf('mu, gain = %.2f',mu_gain))
ylabel('Accumulated Rotation')
%% plot heading over fluorescence data (ovie

t_span = 1:100;
figure(10); clf
subplot(3,1,1)
imagesc(dff_cluster(:,t_span) ./ sum(dff_cluster(:,t_span),1))
subplot(3,1,2:3)
hold on
hh(1) = plot(alpha,nan(2*n_centroid,1),'k.-');  
hh(2) = scatter(0,0,'r','filled')
if movie_flag
for i = t_span
    hh(1).YData = dff_cluster(:,i);
    hh(2).XData = intHD(i);
    drawnow
    pause(.1)
end
end
%% 
tmp = unwrap(mu) - median(circ_dist(mu,-cue)); %- mu(1);
tmp   = mod(tmp,2*pi);                                %rewrap heading data, and put between -pi and pi.
tmp(tmp > pi) = tmp(tmp > pi) - 2*pi;
figure(12); clf
a = plot(tmp,'k','LineWidth',2);
a.YData(abs(diff(a.YData)) > pi) = nan;
hold on
tmp = unwrap(cue); % - cue(1);
tmp   = mod(tmp,2*pi);                                %rewrap heading data, and put between -pi and pi.
tmp(tmp > pi) = tmp(tmp > pi) - 2*pi;
a = plot(-tmp,'b','LineWidth',2);
a.YData(abs(diff(a.YData)) > pi) = nan;
axis tight
ylabel('Azimuth (radians)')
yticks([-pi,0,pi])
ylim([-pi,pi])
yticklabels({'-\pi','0','\pi'})
a = gca;
a.FontSize= 20;
xticks([])
xticklabels([])
legend('PVA','Cue')

% %% load fly movie and interpolate at framerate (time consuming)
% [filename2,filepath2] = uigetfile(filepath,'Select fly movie');
% ft_vid = VideoReader([filepath2,'\',filename2]);
% 
% n = ft_vid.NumFrames;
% vid_lum =  zeros(round(n/2),1);
% for i = 1:round(n/2)
%     vid_lum(i) = sum(read(ft_vid,i),'all');
% end
% vid_start = find(vid_lum > mean(vid_lum),1);
% 
% xf = linspace(0,max(xb),n-vid_start);
% amp_int = interp1(xb,amp,xf);
% r_speed_int = interp1(xb,r_speed,xf);
% img_int = imgData;
% top_int = prctile(img_int,98,'all');                                    %clip the extremes and renormalize for viewing
% bot_int = prctile(img_int,5,'all');
% img_int = max(min(img_int,top_int),bot_int) - bot_int;
% img_int = 255*img_int/max(img_int,[],'all');
% img_int = smoothdata(img_int,3,'movmean',avg_win);                      %smooth the data, again for viewing purposes (should this go before the clipping)            
% img_int = permute(interp1(xb,permute(img_int,[3,1,2]),xf),[2,3,1]);     %interpolate for smooth plotting with fictrac movie
% 

% %% show fly, 2p, and amplitude as a movie
% clear h
% figure(7); clf
% subplot(2,2,1); h(1) = image(read(ft_vid,1+vid_start)); xticks([]); yticks([]); axis equal tight
% subplot(2,2,2); h(2) = image(img_int(:,:,i)); xticks([]); yticks([]); axis equal tight
% hold on
% [y,x] = find(bwmorph(mask,'remove'));                                               %overlay the outline of our mask
% [y,x] = graph_sort(y,x);
% plot(x,y,':','Color',[0.75,0.75,0.75])
% colormap(bone)
% 
% subplot(2,1,2); h(3) = plot(xf,nan(size(r_speed_int)),'k','Linewidth',2); xlabel('time (s)'); ylabel('Rotational Speed'); xlim([0,max(xf)]); ylim([min(r_speed),max(r_speed)])
% yyaxis right; h(4)= plot(xf,nan(size(amp_int)),'Color',c1,'Linewidth',2); xlabel('time (s)'); ylabel('\DeltaF/F'); xlim([0,max(xf)]); ylim([min(amp_int),max(amp_int)])
% ax = gca; ax.YAxis(2).Color = c1;
% for i = 1:(n-vid_start)
%     h(1).CData = read(ft_vid,i+vid_start);
%     h(2).CData = img_int(:,:,i);
%     h(3).YData(i) = r_speed_int(i);
%     h(4).YData(i) = amp_int(i);
%     pause(pause_time)
% end