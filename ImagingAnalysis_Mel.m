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
pause_time  = 0;                         %pause this many seconds between each frame

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

%% extract midline axis, and subsample for centroid locations (using graph search)
n_centroid = 20;                                 %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
[x,y] = find(mask);
mid   = bwskel(mask,'MinBranchLength',min(range(x),range(y)));  %find the midline as the skeleton, shaving out all sub branches
[y,x] = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep    = bwmorph(mid,'endpoints');               %find the endpoints of the midline
[y0,x0] = find(ep,1);
figure(2); imagesc(subplot(2,1,2),mid); colormap(bone); xticks([]); yticks([])
hold on; scatter(x0,y0,'b','filled')

I = knnsearch([x,y],[x,y],'k',3);   %find the 3 nearest neighbors to each point in the midline (including self)
G = graph;                          %initialize a graph object
s = 1:size(x);                      %set starting nodes as the index to each point
G = addedge(G,s,I(:,2));            %add an edge connecting each point to its nearest neighbor
G = addedge(G,s,I(:,3));            %add an edge connecting each points to its second nearest neighbor
k = find(x==x0 & y==y0);            %find the index of our starting point
v = bfsearch(G,k);                  %discover all points in the graph, starting at our starting point. return the index of discover (such that you traverse the shortest path)
x = x(v);                           %reorder the coordinates in order of discovery
y = y(v);

xq = linspace(1,length(y),n_centroid*2)';           %set query points for interpolation
centroids = [interp1(1:length(y),y,xq),interp1(1:length(x),x,xq)]; %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
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
n_frames = size(avg_intensity,2);
alpha = linspace(-pi,pi,n_centroid);                        %create a vector map to represent each centroid as an angle, where the ends represent the same angle (0, 2pi)
muL = nan(n_frames,1);                         %initialize vectors to store the estimate mean and kappa of each hemisphere, if modelled as a von mises distribution
muR = nan(n_frames,1);                         %mu is thetahat
kappaL = nan(n_frames,1);                      %kappa is a measure of concentration. it is bounded [0,1] and is inversely proporitional to variance.
kappaR = nan(n_frames,1);
ampL   = nan(n_frames,1);
ampR   = nan(n_frames,1);
r2L    = nan(n_frames,1);
r2R    = nan(n_frames,1);

for i = 1:size(avg_intensity,2)                                                             %for each frame
    [muL(i), kappaL(i)] = circ_vmpar(alpha,avg_intensity(1:n_centroid,i));                %fit von mises paramters, where each angle is represented by alpha and the strength (or counts) of each angle are represented by the average intensity
    [muR(i), kappaR(i)] = circ_vmpar(alpha,avg_intensity((1:n_centroid) + n_centroid,i)); %do this for each hemisphere separately.
    if any(isnan(circ_vmpdf(alpha,muL(i),kappaL(i)))) || any(isnan(circ_vmpdf(alpha,muR(i),kappaR(i))))
        muL(i) = nan;
        muR(i) = nan;
        kappaL(i) = nan;
        kappaR(i) = nan;
        continue
    end
    [ampL(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muL(i),kappaL(i)) - avg_intensity(1:n_centroid,i)).^2),1);   %fit the amplitude of the von mises by minimizing the sum of the squared difference with scaled VM and data at current frame
    sst                 = sum( (ampL(i)*circ_vmpdf(alpha,muL(i),kappaL(i)) - mean(avg_intensity(1:n_centroid,i)) ).^2);
    r2L(i)              = 1 - ssr/sst;
    [ampR(i),ssr]       = fminsearch(@(a) sum((a*circ_vmpdf(alpha,muR(i),kappaR(i)) - avg_intensity((1:n_centroid) + n_centroid,i)).^2),1);
    sst                 = sum( (ampR(i)*circ_vmpdf(alpha,muR(i),kappaR(i)) - mean(avg_intensity((1:n_centroid)+n_centroid,i)) ).^2);
    r2R(i)              = 1 - ssr/sst;
end

%% plot! movie
c1 = [1,0.5,0];                                                             %define the colors for the left and right bump
c2 = [0,0.5,1];
n_frames = size(avg_intensity,2);
n_ticks  = 4;
pause_time = .5;                                                           %set the pause time

figure(1); clf
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
set(gca,'Position',pos, 'YLim',[min(avg_intensity,[],'all'),max(avg_intensity,[],'all')],...
    'XLim',[0,2*n_centroid],'XTick',tick_vec(idx),'XTickLabels',label_vec(idx)/pi)
xlabel('Cluster (\pi)')
ylabel('Activity (a.u.)')

subplot(2,2,1)                                                                      %initialize another plot to show the movie as the extracted data plays
h(2) = imagesc(imgData(:,:,i));
hold on
[y,x] = find(bwmorph(mask,'remove'));
[y,x] = graph_sort(y,x);
plot(x,y,'w')
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
set(subplot(3,2,4),'YTick',[-pi,0,pi],'YTickLabels',{'-\pi','0','\pi'},'YLim',[-pi,pi])
title(subplot(3,2,2),'Bump Parameters')

for i = 1:n_frames                                       %at each frame, update:
    h(1).YData = avg_intensity(:,i);                                    %the average intensity per centroids
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

%% plot fictract parameters
% if nanmedian(r2L) < 0
%     ampL = mean(avg_intensity(1:n_centroid,:),1);
%     muL  = zeros(size(muL));
% end
% if nanmedian(r2R) < 0
%     ampR = mean(avg_intensity((1:n_centroid)+n_centroid,:),1);
%     muR  = zeros(size(muR));
% end
%smooth_window = 5e-3; %smoothing window, in seconds
total_t = max(ftData.frameTimes{1});
b_smooth = 1; %(n_frames / total_t) * 5e-3;
f_smooth = 20; %(length(ftData.frameTimes{1}) / total_t) * 100e-3;
figure(3); clf
h(1) = subplot(3,1,1);
tmp = smoothdata(ftData.intHD{1},1,'movmean',f_smooth);
tmp(tmp>pi) = tmp(tmp>pi) - 2*pi;
idx = abs(diff(tmp)) > pi;
tmp(idx) = nan;
plot(ftData.frameTimes{1},tmp,'k')
ylabel('heading')
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylim([-pi,pi])
yyaxis right
plot(linspace(0,total_t,n_frames),smoothdata(muL,1,'movmean',b_smooth),'-','Color',c1)
hold on
plot(linspace(0,total_t,n_frames),smoothdata(muR,1,'movmean',b_smooth),'-','Color',c2)
ylabel('Position')
yticks([])

h(2) = subplot(3,1,2);
plot(ftData.frameTimes{1},smoothdata(abs(ftData.yawSpeed{1}),1,'movmean',f_smooth),'k')
ylabel('Rotational speed')
yyaxis right
plot(linspace(0,total_t,n_frames),smoothdata(ampL,1,'movmean',b_smooth),'-','Color',c1)
hold on
plot(linspace(0,total_t,n_frames),smoothdata(ampR,1,'movmean',b_smooth),'-','Color',c2)
ylabel('Amplitude')

h(3) = subplot(3,1,3);
plot(ftData.frameTimes{1},smoothdata(ftData.fwSpeed{1},1,'movmean',f_smooth),'k')
ylabel('Forward speed')
yyaxis right
plot(linspace(0,total_t,n_frames),smoothdata(ampL,1,'movmean',b_smooth),'-','Color',c1)
hold on
plot(linspace(0,total_t,n_frames),smoothdata(ampR,1,'movmean',b_smooth),'-','Color',c2)
ylabel('Amplitude')
xlabel('time (s)')
yyaxis left
y = ylim;
text(total_t,y(2),sprintf('bump smooth: %.2fs \nfictrac smooth: %.2fs',...
     b_smooth*(total_t/n_frames),f_smooth*(total_t/length(ftData.frameTimes{1}))),...
     'VerticalAlignment','top','HorizontalAlignment','right')
 ylim(y)

linkaxes(h,'x')
xlim([0,total_t])