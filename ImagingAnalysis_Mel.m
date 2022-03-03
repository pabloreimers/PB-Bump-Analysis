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

%% extract midline axis, and subsample for centroid locations
n_centroid = 9;                                 %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
smooth_factor = 50;
mid  = bwmorph(mask,'thin',inf);                %thin the mask into a single midline. this is done by thinning infinite times, until mask cannot be thinned any more
[y,x] = find(mid);                              %find the cartesian coordinates of the midline
x = smooth(x,smooth_factor);                    %smooth the coordinates of the midnight. this is set empirically to 15, CHECK THIS. Basically, it makes the midnight monotonic, and not jagged.
y = smooth(y,smooth_factor);
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

%% assign each pixel to a centroid
[y_coor,x_coor] = find(mask);                   %find coordinates of every pixel in the mask
[dists,idx] = pdist2(centroids,[y_coor,x_coor],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask
figure(1); hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_coor(idx == i),y_coor(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end

%% Find the mean activity in each group over time
imgData_2d = reshape(imgData,[],size(imgData,3));               %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log = false(2*n_centroid,size(imgData_2d,1));          %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                          %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_coor(idx==i),x_coor(idx ==i))) = true;
end

avg_intensity = centroid_log * imgData_2d ./ sum(centroid_log,2); %the summed activity in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame

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
pause_time = 0.2;                                                       %set the pause time

figure(1); clf
subplot(2,1,1)
hold on
h(1) = plot([alpha,alpha+2*pi],nan(2*n_centroid,1));                    %initialize a plot to show the average activity at each centroid over time
ylim([0,max(avg_intensity,[],'all')])
xlim([0,4*pi])
xticks([0:pi:4*pi])
xticklabels(0:4)
xlabel('Position (\pi)')
ylabel('Activity')

p(1) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal');      %initialize bump parameters on this, shown as little error bars with mu and sigma (1-kappa)
p(2) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal');

subplot(2,1,2)                                                          %initialize another plot to show the movie as the extracted data plays
h(2) = imagesc(imgData(:,:,i));
axis equal tight
colormap(bone)

for i = 1:size(avg_intensity,2)                                         %at each frame, update:
    h(1).YData = avg_intensity(:,i);                                    %the average intensity per centroids
    h(2).CData = imgData(:,:,i);                                        %the displayed image from the movie
 
    p(1).XData = muL(i);                                                %the mean of the left distribution
    p(1).XNegativeDelta = 1 - kappaL(i);                                %the spread of the left distribution
    p(1).XPositiveDelta = 1 - kappaL(i);
    p(1).YNegativeDelta = ampL(i) / max(ampL);
    p(2).XData = muR(i) + 2*pi;                                         %same for right, just plotted +2pi
    p(2).XNegativeDelta = 1 - kappaR(i);
    p(2).XPositiveDelta = 1 - kappaR(i);
    p(2).YNegativeDelta = ampR(i) / max(ampR);
    
    pause(pause_time)
end

%% the same as above, except overlaying hemispheres
pause_time = 0.1;
c1 = [1,0.5,0];
c2 = [0,0.5,1];

figure(1); clf
subplot(2,1,1)
hold on
h1(1) = plot(alpha,nan(n_centroid,1),'Color',c1);                           %initialize a plot to show intensity at each centroid (represented by an angle, alpha)
h1(2) = plot(alpha,nan(n_centroid,1),'Color',c2);                           %create a unique handle for the left and right PB
ylim([0,max(avg_intensity,[],'all')])
xlim([0,2*pi])
xticks([0:pi:2*pi])
xticklabels(0:2)
xlabel('Position (\pi)')
ylabel('Activity')
legend('Left PB','Right PB','autoupdate','off')
p(1) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal','Color',c1); %same for the bump parameters, displayed as errorbars shows middle and spread of estimated von mises
p(2) = errorbar(nan,max(avg_intensity,[],'all'),nan,'horizontal','Color',c2);

subplot(2,1,2)                                                              %plot the movie in real time, too
h(2) = imagesc(imgData(:,:,i));
axis equal tight
colormap(bone)

for i = 1:size(avg_intensity,2)                                             %for each frame, update:
    h1(1).YData = avg_intensity(1:n_centroid,i);                            %the intensity in the left PB
    h1(2).YData = avg_intensity((1:n_centroid) + n_centroid,i);             %the intrnsity in the right PB
    h(2).CData = imgData(:,:,i);                                            %the shown movie

    p(1).XData = muL(i);                                                    %the bump parameters
    p(1).XNegativeDelta = 1 - kappaL(i);
    p(1).XPositiveDelta = 1 - kappaL(i);
    p(2).XData = muR(i);
    p(2).XNegativeDelta = 1 - kappaR(i);
    p(2).XPositiveDelta = 1 - kappaR(i);
    
    pause(pause_time)
end