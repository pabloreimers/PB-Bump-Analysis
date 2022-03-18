%% reset the workspace and set script parameters
clear all
close all
movie_flag          = false;                                                %set whether to play movies or just run straight through the script
pause_time          = 0;                                                    %set the pause time if playing movies
data_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %set main data directory for ease of selecting files
f0_pct              = 10;                                                   %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 10;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 10;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds

%% ask user for image data
[filename,filepath] = uigetfile('.mat','Select Registered Movie',data_dir);
load([filepath,'\',filename])

%% plot the movie of all planes separately
pause_time  = 0;                         %pause this many seconds between each frame
n_plane = size(regProduct,3);
n_frame = size(regProduct,4);
row    = ceil(sqrt(n_plane));
col     = floor(sqrt(n_plane));

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

%% Play the movie
pause_time  = 0;                         %pause this many seconds between each frame

figure(1); clf                              %clear the current figure
h           = image(imgData(:,:,1));        %initialize an image object where we will update the color values in each frame 
axis equal tight                            %make all pixels square, and keep axis tight for clean look
colormap(bone)                              %set to favorite colormap. I like black and white

if movie_flag
for i = 1:size(imgData,3)                   %loop through each frame and update the ColorData with the values for the current frame
    h.CData = imgData(:,:,i);
    pause(pause_time)
end
end

%% plot the pb and draw bounding box over the whole thing
figure(1); clf                                          % clear the current figure
mask = roipoly(uint8(max(imgData,[],3)));               %this create a black and white mask (logical) of the PB, drawn on a maxZ projection of the image
figure(2)
imagesc(subplot(2,1,1),mask); colormap(bone); xticks([]); yticks([])

%% extract midline axis, and subsample for centroid locations (using graph search)
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length
mid             = bwskel(mask,'MinBranchLength',min(range(x_mask),range(y_mask)));  %find the midline as the skeleton, shaving out all sub branches that are smaller than the minimum axis length
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);
figure(2); imagesc(subplot(2,1,2),mid); colormap(bone); xticks([]); yticks([])
hold on; scatter(x0,y0,'b','filled')

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!

xq          = linspace(1,length(y_mid),n_centroid*2)';                           %set query points for interpolation (the number of centroids we want)
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
cmap        = repmat(cbrewer2('set1',n_centroid),2,1);                     %set a colormap to plot each centroid in a different color, and which repeats per hemisphere (note that if a hemisphere has more clusters than colorbrewer can generate, it will repeat colors within each hemisphere).

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
figure(3); clf
imagesc(subplot(2,2,1),centroid_log); colormap('bone')
xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
imagesc(subplot(2,2,2),imgData_2d); colormap('bone')
xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
imagesc(subplot(2,1,2),dff_cluster); colormap('bone')
xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')

%% Estimate von mises parameters
n_frames    = size(dff_cluster,2);
alpha       = linspace(-pi,pi,n_centroid);          %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
muL         = nan(n_frames,1);                      %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
muR         = nan(n_frames,1);                      %mu is thetahat
kappaL      = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
kappaR      = nan(n_frames,1);
ampL        = nan(n_frames,1);                      %this will be used to scale the bump amplitude from a generic von mises pdf (integral 1) to something that matches the data
ampR        = nan(n_frames,1);
r2L         = nan(n_frames,1);                      %rsquared values, telling how much of the variance in the distribution is explained by our fit (as compared to a flat line which is the mean)
r2R         = nan(n_frames,1);

for i = 1:n_frames                                                                      %for each frame
    [muL(i), kappaL(i)] = circ_vmpar(alpha,dff_cluster(1:n_centroid,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
    [muR(i), kappaR(i)] = circ_vmpar(alpha,dff_cluster((1:n_centroid) + n_centroid,i)); %do this for each hemisphere separately.
    if any(isnan(circ_vmpdf(alpha,muL(i),kappaL(i)))) || any(isnan(circ_vmpdf(alpha,muR(i),kappaR(i)))) %if the distribution returns nans (calculated something too sharp, e.g.)
        muL(i) = nan;
        muR(i) = nan;
        kappaL(i) = nan;
        kappaR(i) = nan;
        continue
    end
    [ampL(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muL(i),kappaL(i)) - dff_cluster(1:n_centroid,i)).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
    sst                 = sum( (ampL(i)*circ_vmpdf(alpha,muL(i),kappaL(i)) - mean(dff_cluster(1:n_centroid,i)) ).^2);        %find the squared error with the mean
    r2L(i)              = 1 - ssr/sst;                                                                                       %calculate r2 (which is 1- ratio of our fit/baseline fit).
    [ampR(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muR(i),kappaR(i)) - dff_cluster((1:n_centroid) + n_centroid,i)).^2),1);
    sst                 = sum( (ampR(i)*circ_vmpdf(alpha,muR(i),kappaR(i)) - mean(dff_cluster((1:n_centroid)+n_centroid,i)) ).^2);
    r2R(i)              = 1 - ssr/sst;
end

%% plot! movie
c1 = [1,0.5,0];                                                             %define the colors for the left and right bump
c2 = [0,0.5,1];
n_frames = size(dff_cluster,2);
n_ticks  = 4;

figure(4); clf
subplot(2,2,3)
hold on
clear h v p
h(1) = plot(1:(2*n_centroid),nan(2*n_centroid,1),'k.-');                    %initialize a plot to show the average activity at each centroid over time
v(1)= plot(1:n_centroid,nan(n_centroid,1),'Color',c1);
v(2)= plot((1:n_centroid) + n_centroid, nan(n_centroid,1),'Color',c2);
pos = get(gca,'Position');
legend('Activity','Left PB fit','Right PB fit','autoupdate','off','Location','northoutside')
tick_vec  = 1:2*n_centroid;
label_vec = [alpha,alpha];
idx       = linspace(0,2*n_centroid,n_ticks+1);
idx(1)    = [];
set(gca,'Position',pos, 'YLim',[min(dff_cluster,[],'all'),max(dff_cluster,[],'all')],...
    'XLim',[0,2*n_centroid],'XTick',tick_vec(idx),'XTickLabels',label_vec(idx)/pi)
xlabel('Cluster (\pi)')
ylabel('Activity (a.u.)')

subplot(2,2,1)                                                                      %initialize another plot to show the movie as the extracted data plays
h(2) = image(imgData(:,:,i));
title('raw F')
hold on
[y,x] = find(bwmorph(mask,'remove'));                                               %overlay the outline of our mask
[y,x] = graph_sort(y,x);
plot(x,y,':','Color',[0.75,0.75,0.75])
axis equal tight
colormap(bone)
xticks([])
yticks([])

bump_params = {ampL,ampR;muL,muR;kappaL,kappaR};                                                    %store all of the bump statistics in a cell array for easy repetitive access in a loop
bump_params = cellfun(@(x)(smoothdata(x,1,'gaussian',b_smooth)),bump_params,'UniformOutput',false); %smooth all of the bump parameters.
bp_names    = {'Amplitude','Position (\mu)','Concentration (\kappa)'};
bp_h        = cell(3,2);

for i = 1:size(bump_params,1)                                                       %initialize plots for each bump statistic over each frame
    set(subplot(3,2,2*i),'NextPlot','add','xlim',[0,n_frames],...
        'ylim',[min([bump_params{i,:}],[],'all'),max([bump_params{i,:}],[],'all')])
    ylabel(subplot(3,2,2*i),bp_names{i})
    bp_h{i,1} = plot(subplot(3,2,2*i),1:n_frames,nan(n_frames,1),'Color',c1);
    bp_h{i,2} = plot(subplot(3,2,2*i),1:n_frames,nan(n_frames,1),'Color',c2);
end
xlabel('frame')
set(subplot(3,2,4),'YTick',[-pi,0,pi],'YTickLabels',{'-\pi','0','\pi'},'YLim',[-pi,pi])
title(subplot(3,2,2),'Bump Parameters')

for i = 1:n_frames                                       %at each frame, update:
    h(1).YData = dff_cluster(:,i);                                    %the average intensity per centroids
    h(2).CData = imgData(:,:,i);                                        %the displayed image from the movie

    v(1).YData = ampL(i)*circ_vmpdf(alpha,muL(i),kappaL(i));            %generate the fit von mises distribution, over the range of thetas, where it is centered at this time point, and how concentrated. scale it, too
    v(2).YData = ampR(i)*circ_vmpdf(alpha,muR(i),kappaR(i));
    
    for j = 1:size(bump_params,1)                                       %update the current distribution statistics, too
        bp_h{j,1}.YData(i) = bump_params{j,1}(i);
        bp_h{j,2}.YData(i) = bump_params{j,2}(i);
    end
    
    if i > 1                                                            %because circular, if the instantaneous position jump is more than pi radians, just blank the point (because the shortest path is actually out of frame)
    if abs(bump_params{2,1}(i) - bump_params{2,1}(i-1)) > pi
        bp_h{2,1}.YData(i) = nan;
    end
    if abs(bump_params{2,2}(i) - bump_params{2,2}(i-1)) > pi
        bp_h{2,2}.YData(i) = nan;
    end
    end
    
    pause(pause_time)
end

%% load fictrac data
[filename2,filepath2] = uigetfile(filepath);
load([filepath2,'\',filename2])

%% Process fictrac data
f_speed = ftData.fwSpeed{:};                            %store each speed
r_speed = ftData.yawSpeed{:};
intHD   = ftData.intHD{:};

f_speed = abs(f_speed);                                 %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'gaussian',f_smooth); 
r_speed = smoothdata(r_speed,1,'gaussian',f_smooth);
intHD   = smoothdata(intHD, 1, 'gaussian',f_smooth);
end

intHD = mod(intHD,2*pi);                                %rewrap heading data, and put between -pi and pi.
intHD(intHD > pi) = intHD(intHD > pi) - 2*pi;
total_t = max(ftData.frameTimes{:});                    %find how long the trial was

xf  = ftData.frameTimes{:};                             %create vectors with timestamps for each trace
xb  = linspace(0,total_t,size(imgData,3))';
fr  = mean(diff(xb));

f_speed = interp1(xf,f_speed,xb);                       %interpolate the fictrac trace to the timepoints of the image (bump)
r_speed = interp1(xf,r_speed,xb);
intHD   = interp1(xf,intHD,xb);

%% Plot bump params over fictrac data
muL     = smoothdata(muL,1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
muR     = smoothdata(muR,1,'gaussian',b_smooth);
ampL    = smoothdata(ampL,1,'gaussian',b_smooth);
ampR    = smoothdata(ampR,1,'gaussian',b_smooth);
kappaL  = smoothdata(kappaL,1,'gaussian',b_smooth);
kappaR  = smoothdata(kappaR,1,'gaussian',b_smooth);

figure(5); clf

h(1) = subplot(3,1,1);
tmp = intHD; tmp(abs(diff(tmp)) > pi) = nan;
plot(xb,tmp,'k'); ylabel('Heading (rad)'); ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'});
yyaxis right; hold on
plot(xb,muL,'-','Color',c1)
plot(xb,muR,'-','Color',c2)
ylabel('Bump Position (\mu)'); yticks([])
axis tight

h(2) = subplot(3,1,2);
plot(xb,f_speed,'k'); ylabel('Abs Forward Speed (mm/s)');
yyaxis right; hold on
plot(xb, ampL,'-','Color',c1)
plot(xb, ampR,'-','Color',c2)
ylabel('Bump Amplitude'); yticks([]); axis tight

h(3) = subplot(3,1,3);
plot(xb,r_speed,'k'); ylabel('Rotational Speed (mm/s)');
yyaxis right; hold on
plot(xb, ampL,'-','Color',c1)
plot(xb, ampR,'-','Color',c2)
ylabel('Bump Amplitude'); yticks([]); axis tight
xlabel('time (s)')

linkaxes(h,'x')

%% linear model fictrac params with dF/F
mask_1d = reshape(mask,[],1);                           % reshape the mask into a single vector
f       = reshape(imgData,size(mask_1d,1),[]);          % reshape the movie into pixels x frames
f       = f(mask_1d,:);                                 % extract only pixels within the mask
f0      = prctile(f,f0_pct,'dim',2);                    % find the baseline fluorescence in each pixel over time
dff     = (f - f0) ./ f0;                               % calculate pixelwise df/f
dff     = mean(dff,1)';
dff     = smoothdata(dff,1,'gaussian',b_smooth);        % smooth dff

dff     = smoothdata(mean(dff_cluster,1)',1,'gaussian',b_smooth);
amp     = dff;

%amp     = smoothdata(mean([ampL,ampR],2),1,'gaussian',b_smooth);

[f_corr,f_lag]  = xcorr(f_speed,amp,ceil(max_lag/fr));
[r_corr,r_lag]  = xcorr(r_speed,amp,ceil(max_lag/fr));
[~,l_idx]       = max(f_corr);
f_lag           = f_lag(l_idx);
[~,l_idx]       = max(r_corr);
r_lag           = r_lag(l_idx);

f_idx           = xb - f_lag*fr > 0 & xb - f_lag*fr < max(xb);
r_idx           = xb - r_lag*fr > 0 & xb - r_lag*fr < max(xb);

f_speed(~f_idx) = 0;
r_speed(~r_idx) = 0;

f_speed_lag     = circshift(f_speed,round(-f_lag*fr));
r_speed_lag     = circshift(r_speed,round(-r_lag*fr));

amp(isnan(amp)) = 0;

[f_fit,f_gof]   = fit(f_speed_lag,amp,'poly1');
[r_fit,r_gof]   = fit(r_speed_lag,amp,'poly1');
[m_fit,m_gof]   = fit([r_speed_lag,f_speed_lag],amp,'poly11');

figure(6);clf
subplot(2,1,1)
plot(xb,f_speed,'k'); ylabel('Abs Forward Speed (mm/s)')
yyaxis right; plot(xb, amp,'Color',c1); ylabel('Average dF/F'); ax = gca; ax.YAxis(2).Color = c1;
axis tight
y = ylim; text(xb(end),y(2),sprintf('r^2: %.2f\n2p lag: %.1fms\nball: %.1fs\n2p: %.1fs',...
                f_gof.adjrsquare,f_lag*fr*1e3,f_smooth*fr,b_smooth*fr),...
                'HorizontalAlignment','right','VerticalAlignment','top')

subplot(2,1,2)
plot(xb,r_speed,'k'); ylabel('Rotational Speed (mm/s)')
yyaxis right; plot(xb, amp,'Color',c1); ylabel('Average dF/F'); ax = gca; ax.YAxis(2).Color = c1;
xlabel('time (s)')
axis tight
y = ylim; text(xb(end),y(2),sprintf('r^2: %.2f\n2p lag: %.1fms\njoint r^2:%.2f',...
                                r_gof.adjrsquare,r_lag*fr*1e3,m_gof.adjrsquare),...
                'HorizontalAlignment','right','VerticalAlignment','top')
