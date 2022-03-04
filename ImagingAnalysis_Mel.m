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

for f = 1:size(regProduct,4)
    for i = 1:size(regProduct,3)
        h(i).CData = regProduct(:,:,i,f);
    end
    pause(pause_time)
end

%% Play the movie
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
pause_time  = 0.05;                         %pause this many seconds between each frame

figure(1); clf                              %clear the current figure
h           = imagesc(imgData(:,:,1));      %initialize an image object where we will update the color values in each frame 
axis equal tight                            %make all pixels square, and keep axis tight for clean look
colormap(bone)                              %set to favorite colormap. I like black and white

for i = 1:size(imgData,3)                   %loop through each frame and update the ColorData with the values for the current frame
    h.CData = imgData(:,:,i);
    pause(pause_time)
end

%% plot the pb and draw bounding box over the whole thing
figure(1); clf              % clear the current figure
imagesc(max(imgData,[],3))  %(max intensity over time of summed intensity over planes)
colormap(bone)
axis equal tight
mask = roipoly;             %this create a black and white mask (logical) of the PB
figure(2)
imagesc(subplot(2,1,1),mask); colormap(bone); xticks([]); yticks([])
%% extract midline axis, and subsample for centroid locations
n_centroid = 9;                                 %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
mid   = bwmorph(mask,'thin',inf);               %find the cartesian coordinates of the midline
[y,x] = find(mid);
ep    = bwmorph(mid,'endpoints');               %find the endpoints of the midline
[y0,x0] = find(ep,1);

r_com = sum([x,y]) / numel(x);                  %now find the center of mass of all the points, and sort by the angle to center of mass to order them appropriately
[~,I] = sort(atan2(y-r_com(2), x-r_com(1)),'descend');
x = x(I);
y = y(I);
k = find(x==x0 & y==y0);                        %point are sorted in order, but we want to start at the end points
x = circshift(x,-(k-1),1);                       %so find how far offset the starting point is in the array, and circshift it backwards.
y = circshift(y,-(k-1),1);

idx = round(linspace(1,length(y),n_centroid*2));  %linearly subsample the midline into the desired number of centroids.
centroids = [y(idx),x(idx)];
cmap = repmat(cbrewer2('set1',n_centroid),2,1);   %set a colormap to plot each centroid in a different color, and which repeats per hemisphere

figure(1); clf
imagesc(max(imgData,[],3))                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
plot(x,y,'w')                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled') %show the centroids in each of their colors
axis equal tight
figure(2); imagesc(subplot(2,1,2),mid); colormap(bone); xticks([]); yticks([])
hold on; scatter(r_com(1),r_com(2),'r','filled'); scatter(x0,y0,'b','filled')

%% assign each pixel to a centroid
[y_coor,x_coor] = find(mask);                   %find coordinates of every pixel in the mask
[dists,idx] = pdist2(centroids,[y_coor,x_coor],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask
figure(1); hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_coor(idx == i),y_coor(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end

%% Find the mean activity in each group over time
f0_pct = 10;

imgData_2d = reshape(imgData,[],size(imgData,3));               %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log = false(2*n_centroid,size(imgData_2d,1));          %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                          %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_coor(idx==i),x_coor(idx ==i))) = true;
end

avg_intensity = centroid_log * imgData_2d ./ sum(centroid_log,2); %the summed activity in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f0            = prctile(avg_intensity,f0_pct,2);
avg_intensity = (avg_intensity - f0) ./ f0;
figure(1); clf
imagesc(subplot(2,2,1),centroid_log)
xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
imagesc(subplot(2,2,2),imgData_2d)
xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
imagesc(subplot(2,1,2),avg_intensity)
xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')
%% Estimate von mises parameters
alpha = linspace(0,2*pi,n_centroid);                        %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
muL = nan(size(avg_intensity,2),1);                         %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
muR = nan(size(avg_intensity,2),1);                         %mu is thetahat
kappaL = nan(size(avg_intensity,2),1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
kappaR = nan(size(avg_intensity,2),1);
ampL   = nan(size(avg_intensity,2),1);
ampR   = nan(size(avg_intensity,2),1);

for i = 1:size(avg_intensity,2)                                                             %for each frame
    [muL(i), kappaL(i)] = circ_vmpar(alpha,avg_intensity(1:n_centroid,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
    [muR(i), kappaR(i)] = circ_vmpar(alpha,avg_intensity((1:n_centroid) + n_centroid,i)); %do this for each hemisphere separately.
    ampL(i)             = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muL(i),kappaL(i)) - avg_intensity(1:n_centroid,i)).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
    ampR(i)             = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muR(i),kappaL(i)) - avg_intensity((1:n_centroid) + n_centroid,i)).^2),1);
    if muL(i) < 0                                                                           %circ_stats does things from -pi:pi, but for convenience I will keep everything 0:2pi
        muL(i) = muL(i) + 2*pi;
    end
    if muR(i) < 0
        muR(i) = muR(i) + 2*pi;
    end
end

%% plot!
c1 = [1,0.5,0];                                                             %define the colors for the left and right bump
c2 = [0,0.5,1];
n_frames = size(avg_intensity,2);
pause_time = 0;                                                           %set the pause time

figure(1); clf
subplot(2,2,3)
hold on
h(1) = plot([alpha,alpha+2*pi],nan(2*n_centroid,1),'k.-');                    %initialize a plot to show the average activity at each centroid over time
v(1) = plot(alpha,nan(n_centroid,1),'Color',c1);
v(2) = plot(alpha + 2*pi, nan(n_centroid,1),'Color',c2);
pos = get(gca,'Position');
legend('Activity','Left PB fit','Right PB fit','autoupdate','off','Location','northoutside')
set(gca,'Position',pos, 'YLim',[min(avg_intensity,[],'all'),max(avg_intensity,[],'all')],...
    'XLim',[0,4*pi],'XTick',[alpha,2*pi+alpha(2:end)],'XTickLabels',[alpha,alpha(2:end)]/pi)
xlabel('Cluster (\pi)')
ylabel('Activity (a.u.)')

p(1) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal','Color',c1);      %initialize bump parameters on this, shown as little error bars with mu and sigma (1-kappa)
p(2) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal','Color',c2);

subplot(2,2,1)                                                                      %initialize another plot to show the movie as the extracted data plays
h(2) = imagesc(imgData(:,:,i));
axis equal tight
colormap(bone)
xticks([])
yticks([])

bump_params = {ampL,ampR;muL,muR;kappaL,kappaR};                                    %store all of the bump statistics in a cell array for easy repetitive access in a loop
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
set(subplot(3,2,4),'YTick',[0,pi,2*pi],'YTickLabels',{'0','\pi','2\pi'},'YLim',[0,2*pi])
title(subplot(3,2,2),'Bump Parameters')

for i = 1:n_frames                                       %at each frame, update:
    h(1).YData = avg_intensity(:,i);                                    %the average intensity per centroids
    h(2).CData = imgData(:,:,i);                                        %the displayed image from the movie
 
    p(1).XData = muL(i);                                                %the mean of the left distribution
    p(1).XNegativeDelta = 1 - kappaL(i);                                %the spread of the left distribution
    p(1).XPositiveDelta = 1 - kappaL(i);
    p(2).XData = muR(i) + 2*pi;                                         %same for right, just plotted +2pi
    p(2).XNegativeDelta = 1 - kappaR(i);
    p(2).XPositiveDelta = 1 - kappaR(i);

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