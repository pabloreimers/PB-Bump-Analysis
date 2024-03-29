%% reset the workspace and set script parameters
clear all
close all
movie_flag          = false;                                                %set whether to play movies or just run straight through the script
rsn_flag            = false;
pause_time          = 0;                                                    %set the pause time if playing movies
data_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %set main data directory for ease of selecting files
%data_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\2p data\';
f0_pct              = 5;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = .25*10*5;                                              %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter. gaussian sigma is 1/5 the window length, so this is sigma*fr*5, in seconds
f_smooth            = .5*60;                                              %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 1;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = 0;%5e-3;
z_flag = true;
%% ask user for image data
[filename,filepath] = uigetfile('.mat','Select Registered Movie',data_dir);
load([filepath,'\',filename])

%% plot the movie of all planes separately
%regProduct = imgData;
n_plane = size(regProduct,3);
n_frame = size(regProduct,4);
row     = ceil(sqrt(n_plane));
col     = ceil(sqrt(n_plane));

figure(1); clf
for i = 1:n_plane
    h(i) = imagesc(subplot(row,col,i),mean(regProduct(:,:,i,:),4));
    axis equal tight
    xticks([])
    yticks([])
    title(['plane: ',num2str(i)])
    colormap(bone)
end

if movie_flag
for f = 1650:size(regProduct,4)
    for i = 1:size(regProduct,3)
        h(i).CData = mean(regProduct(:,:,i,f:f),4);
    end
    pause(pause_time)
end
end

%% remove scan noise from image
regProduct = smoothdata(regProduct,4,'movmean',5);

if rsn_flag
rm_bands = 5;
dims = size(regProduct);

im_2d = double(reshape(permute(regProduct,[2,1,3,4]),size(regProduct,2),[]));
n = size(im_2d,1);
b = mean(im_2d,1);
x = im_2d - b;
fhat = fft(x,[],1);
fhat(rm_bands:(n-rm_bands+2),:) = 0;
im_rsn = ifft(fhat,[],1) + b;

regProduct = permute(reshape(im_rsn,dims(2),dims(1),dims(3),dims(4)),[2,1,3,4]);
end

%% Make compress in z, and scale between 0 and 256
%regProduct = smoothdata(regProduct(:,:,1:end,:),4,'gaussian',b_smooth);

imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
h0 = prctile(imgData(:),1);
h1 = prctile(imgData(:),99);
imgData     = 256*(imgData - h0)/(h1 - h0); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)

imgData2 = imgData;
top_int = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
bot_int = prctile(imgData2,1,'all');
imgData2 = max(min(imgData2,top_int),bot_int) - bot_int;
imgData2 = 256*imgData2/max(imgData2,[],'all');
%imgData2 = smoothdata(imgData2,3,'movmean',avg_win);                      %smooth the data, again for viewing purposes (should this go before the clipping)            

%% Play the movie
figure(1); clf                              %clear the current figure
h           = image(imgData(:,:,1));        %initialize an image object where we will update the color values in each frame 
axis equal tight                            %make all pixels square, and keep axis tight for clean look
colormap(bone)                              %set to favorite colormap. I like black and white
colorbar

if movie_flag
for i = 1650:size(imgData,3)                   %loop through each frame and update the ColorData with the values for the current frame
    h.CData = imgData(:,:,i);
    pause(pause_time)
end
end

%% plot the eb and draw bounding box over the whole thing
figure(1); clf                                          % clear the current figure
tmp = mean(imgData,3);
tmp = 256*(tmp - min(tmp(:)))/(max(tmp(:) - min(tmp(:))));
image(tmp)
tmp = drawellipse('Center',[size(tmp,2)/2,size(tmp,1)/2],'SemiAxes',[size(tmp,2)/4,size(tmp,1)/4]);
input('') %move on once user presses enter, after adjusting the ellipse
mask = createMask(tmp);
tmp.SemiAxes = tmp.SemiAxes / 3;
mask = logical(mask - createMask(tmp));
%mask = roipoly(uint8(tmp));               %this create a black and white mask (logical) of the PB, drawn on a maxZ projection of the image
figure(2)
imagesc(subplot(2,1,1),mask); colormap(bone); xticks([]); yticks([])

%% extract midline axis, and subsample for centroid locations (using graph search)
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
mid             = bwmorph(mask,'remove');
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);
figure(2); imagesc(subplot(2,1,2),mid); colormap(bone); xticks([]); yticks([])
hold on; scatter(x0,y0,'b','filled')

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
cmap        = make_colors(n_centroid); %repmat(cbrewer2('set1',n_centroid),2,1);                     %set a colormap to plot each centroid in a different color, and which repeats per hemisphere (note that if a hemisphere has more clusters than colorbrewer can generate, it will repeat colors within each hemisphere).

if centroids(1,1) < centroids(end,1)
    centroids = flipud(centroids);
end

figure(1); clf
imagesc(mean(imgData,3))                                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%% assign each pixel to a centroid
[~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)
figure(1); clf
imagesc(mean(imgData,3))                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
for i = 1:n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_mask(idx == i),y_mask(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%% Find the mean activity in each group over time
tmp     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
tmp     = reshape(tmp,[],size(tmp,3)); %reshape imgData to be allpixels x frames
background  = tmp(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > inf;
%flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
%flash_idx   = med_pix > 0.05;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));


%tmp = smoothdata(imgData,3,'movmean',5);
imgData_2d      = reshape(tmp,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
%dff_cluster = zscore(dff_cluster,[],2);
zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
tmp = dff_cluster;
tmp(:,~flash_idx) = zscore_cluster;
zscore_cluster = tmp;
dff_cluster = zscore_cluster;
%dff_cluster     = smoothdata(dff_cluster,2,'gaussian',b_smooth);
figure(3); clf
imagesc(subplot(2,2,1),centroid_log); colormap('bone')
xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
imagesc(subplot(2,2,2),imgData_2d); colormap('bone')
xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
imagesc(subplot(2,1,2),dff_cluster); colormap('bone'); pos = get(gca,'position'); colorbar; set(gca,'Position',pos)
xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')
hold on
scatter(ones(n_centroid,1),1:n_centroid,[],cmap,'filled')

 %% Find the mean activity in each group over time, 3 dimensional
% tmp     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
% tmp     = reshape(tmp,[],size(tmp,3)); %reshape imgData to be allpixels x frames
% background  = tmp(~reshape(mask,[],1),:); %extract points that are outside of the mask)
% med_pix     = sum(background,1);
% med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
% flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
% %flash_idx   = med_pix > 0.05;
% flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));
% 
% 
% n_planes = size(regProduct,3);
% %y_scaled = (centroids(:,1) - min(centroids(:,1)))/ (max(centroids(:,1)) - min(centroids(:,1)));
% filt_mu = linspace(max(centroids(:,1)),min(centroids(:,1)),n_planes);
% filt_sig= range(centroids(:,1));
% filt_x  = [1:size(regProduct,1)]';
% %p_scaled = linspace(1,0,n_planes);
% %y_scaled = repmat(linspace(0,1,size(regProduct,1))',1,size(regProduct,2));
% f_cluster = zeros(n_centroid,n_frame);
% 
% for p = 1:n_planes
% %z_mult = abs(y_scaled - p_scaled(p));   %multiply teh values in each cluster by a multiplier which depends on the plane we are in. dorsal (bottom) is passed through at early planes, ventral passed through at late planes
% %tmp = double(squeeze(regProduct(:,:,p,:))) .* (1-abs(y_scaled - p_scaled(p)));
% tmp = double(squeeze(regProduct(:,:,p,:))) .* normpdf(filt_x,filt_mu(p),filt_sig);
% 
% 
% imgData_2d      = reshape(double(tmp),[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
% centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
% for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
%     centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
% end
% 
% tmp       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
% f_cluster = f_cluster + tmp;
% end
% 
% f_cluster(:,flash_idx) = nan;
% f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
% dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
% zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
% tmp = dff_cluster;
% tmp(:,~flash_idx) = zscore_cluster;
% zscore_cluster = tmp;
% if z_flag
% dff_cluster = zscore_cluster;
% end
% %dff_cluster     = smoothdata(dff_cluster,2,'gaussian',b_smooth);
% figure(3); clf
% imagesc(subplot(2,2,1),centroid_log); colormap('bone')
% xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
% imagesc(subplot(2,2,2),imgData_2d); colormap('bone')
% xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
% imagesc(subplot(2,1,2),dff_cluster); colormap('bone'); pos = get(gca,'position'); colorbar; set(gca,'Position',pos)
% xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')
% hold on
% scatter(ones(n_centroid,1),1:n_centroid,[],cmap,'filled')

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

%% Find bump as PVA
tmp_data    = dff_cluster;
alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
alpha       = linspace(-pi,pi,n_centroid);

[x_tmp,y_tmp]   = pol2cart(alpha,tmp_data');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
% for i = 1:n_smooth
% mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
% rho    = smoothdata(rho,1,'gaussian',b_smooth);
% end

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

% if mode(diff(cue)) == 0
%     cue = -cumsum(ftData_DAQ.velYaw{:} / 6);
% end

f_speed = f_speed;                                      %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi
cue     = unwrap(cue / 192 * 2*pi - pi);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'rlowess',f_smooth); 
r_speed = smoothdata(r_speed,1,'rlowess',f_smooth);
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
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp         = interp1(xb,amp,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';


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

[f_corr,f_lag]  = xcorr(f_speed,amp,ceil(max_lag*fr),'coeff');      % find the max correlation, in steps
[r_corr,r_lag]  = xcorr(r_speed,amp,ceil(max_lag*fr),'coeff');
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
lag = 15;
rho_idx = rho>rho_thresh;
[pva_corr,pva_pval] = circ_corrcc(cue(1:end-ceil(lag)),-mu((ceil(lag)+1):end));

figure(11); clf
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
gains = nan(20,20);
rhos = nan(20,20);

for ii = 1
    for jj = 3
lag = jj*5;
bump_vel = [diff(unwrap(mu));0] / fr; %divide by the frame rate (seconds per frame). since the mu is in units of radians, this should give rad per seconds
fly_vel  = ftData_DAQ.velYaw{:};
bump_vel = bump_vel(lag+1:end);
fly_vel  = fly_vel(1:end-lag);
rho_thresh = 0.1; %prctile(rho,10);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
fly_vel = smoothdata(fly_vel,1,'rlowess',ii*10);
bump_vel= smoothdata(bump_vel,1,'rlowess',ii*10);
end
vel_idx = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min & rho(lag+1:end) > rho_thresh; %ignore outlier bump speeds with arbitrary threshold

figure(9); clf
scatter(fly_vel(vel_idx),bump_vel(vel_idx),5,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.1)
xlabel('Fly Vel (rad/s)'); ylabel('Bump Vel (rad/s)');
hold on
[vel_rho,vel_pval] = corr(bump_vel(vel_idx),fly_vel(vel_idx));
x = xlim;
y = ylim;
b = [ones(sum(vel_idx),1),fly_vel(vel_idx)]\bump_vel(vel_idx); %fit the slope of fly vel and bump vel with an arbitrary offset
text(x(2),y(2),sprintf('$r = %.2f$\n$gain = %.2f$\n$lag = %.0f ms$',vel_rho,b(2),lag*fr*1000),'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
plot([0,0],y,'k:')
plot(x,[0,0],':k')
plot(fly_vel(vel_idx),fly_vel(vel_idx)*b(2) + b(1),'r')
axis tight

gains(ii,jj) = b(2);
rhos(ii,jj) = vel_rho;
%drawnow
%title(ii)
%pause(.2)
fprintf('ii:%i, jj:%i\n',ii,jj)
    end
end

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
plot(xf,unwrap(cue) - cue(1),'k','linewidth',2)
plot(xf, -unwrap(mu) + mu(1),'b','linewidth',2)
legend(sprintf('cue, gain = %.2f',cue_gain),sprintf('mu, gain = %.2f',mu_gain))
ylabel('Accumulated Rotation')
%% plot heading over fluorescence data (ovie
% 
% t_span = 1:100;
% figure(10); clf
% subplot(3,1,1)
% imagesc(dff_cluster(:,t_span) ./ sum(dff_cluster(:,t_span),1))
% subplot(3,1,2:3)
% hold on
% hh(1) = plot(alpha,nan(2*n_centroid,1),'k.-');  
% hh(2) = scatter(0,0,'r','filled')
% if movie_flag
% for i = t_span
%     hh(1).YData = dff_cluster(:,i);
%     hh(2).XData = intHD(i);
%     drawnow
%     pause(.1)
% end
% end


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

%% plot fly in cartesian coordinates

f_vel = ftData_DAQ.velFor{:};
s_vel = ftData_DAQ.velSide{:};
theta = ftData_DAQ.cueAngle{:} / 180 * pi;

x_pos = cumsum(f_vel.*cos(theta') + s_vel.*sin(theta'))*fr;
y_pos = cumsum(f_vel.*sin(theta') + s_vel.*cos(theta'))*fr;

figure(15)
subplot(2,1,1)
scatter(x_pos,y_pos,[],1:length(x_pos),'.'); axis equal
xlabel('x pos (mm)'); ylabel('y pos (mm)')
axis tight
subplot(2,2,3)
scatter(xf,theta,[],1:length(theta),'.'); ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'}); ylabel('heading (rad)'); xlabel('time (s)'); axis tight
subplot(2,2,4)
scatter(xf,f_vel,[],1:length(f_vel),'.'); ylabel('forward vel (mm/s)'); xlabel('time (s)'); axis tight
%% Play the movie
pause_time  = 0;                         %pause this many seconds between each frame

figure(14); clf                              %clear the current figure
set(gcf,'color','none')
subplot(3,2,[3,5])
set(gca,'color','none','ycolor','w','xcolor','w')
h           = image(imgData(:,:,1));        %initialize an image object where we will update the color values in each frame 
axis equal tight                            %make all pixels square, and keep axis tight for clean look
colormap(bone)                              %set to favorite colormap. I like black and white
hold on
[y_out,x_out] = find(bwmorph(mask,'remove'));
[y_out,x_out] = graph_sort(y_out,x_out);
plot(x_out,y_out,'w:','linewidth',2)
scatter(x_mask,y_mask,5,cmap(idx,:),'filled','MarkerFaceAlpha',.05);
xticks([])
yticks([])
subplot(3,2,[4,6])
h2 = polarscatter(0,1,1e3,'m','filled');
set(gca,'ThetaAxisUnits','radians','ThetaTick',(0:1/4:2)*pi)
rlim([0,1]); rticks([])
%ylim([min(fly_vel),max(fly_vel)])
title('Fly Heading (rad)','color','w')
set(gca,'color','none','thetacolor','w')
[~,tmp] = sort(idx,'ascend');
% subplot(2,3,5)
% h6 = scatter(x_mask,y_mask,5,cmap(idx,:),'filled');
% set(gca,'color','none','YDir','Reverse')
subplot(3,1,1)
h3 = imagesc(xb,alpha,dff_cluster);
h3.CData = nan(size(dff_cluster));
hold on
h5 = plot(xb,nan(1,length(xb)),'Color',[.75,0,.75],'Linewidth',2);
h4 = plot(xb,nan(1,length(xb)),'Color','w','Linewidth',2);
cue_tmp = interp1(xf,cue,xb);
mu_tmp = interp1(xf,mu,xb);
%vel_tmp = interp1(xf,fly_vel,xb);
if isnan(mu(1))
    mu(1) = 0;
end
cue_tmp = mod(unwrap(cue_tmp) - cue(1) - mu(1),2*pi);
cue_tmp(cue_tmp>pi) = cue_tmp(cue_tmp>pi) - 2*pi;
xlabel('time (s)')
ylabel('azimuth (rad)')
legend('heading','pva','color','none','textcolor','w')
yticks([-pi,0,pi])
yticklabels({'-\pi','0','\pi'})
set(gca,'color','none','xcolor','w','ycolor','w')

fontsize(gcf,20,'pixels')
a_values = imgData_2d(reshape(mask,[],1),:);
a_values = a_values - prctile(a_values(:),10);
a_values = (a_values / prctile(a_values(:),99))*256;

if movie_flag
    %v = VideoWriter('sample.mp4');
    %v.FrameRate = 60;
    %open(v)
for i = 1:size(imgData,3)                   %loop through each frame and update the ColorData with the values for the current frame
    h.CData = imgData(:,:,i);
    h2.ThetaData = -cue_tmp(i);
    h3.CData(:,i) = dff_cluster(:,i);
    % h6.AlphaData = a_values(:,i);
    % h6.MarkerFaceAlpha = 'flat';
    if i>1 && abs(mu_tmp(i) - mu_tmp(i-1)) < pi
        h4.YData(i) = mu_tmp(i);
    end
    if i>1 && abs(cue_tmp(i) - cue_tmp(i-1)) < pi
        h5.YData(i) = -cue_tmp(i);
    end
    pause(pause_time)
    %writeVideo(v,getframe(gcf))
end
end

h3.CData = dff_cluster;
h4.YData = mu_tmp;
h4.YData(abs(diff(h4.YData)) > pi) = nan;
h5.YData = -cue_tmp;
h5.YData(abs(diff(h5.YData)) > pi) = nan;
% if movie_flag
%     close(v)
% end
%% functions
function [y_warp,x_warp] = circ_warp(y,x) %create a function to warp coordinates so that the x and y axis are the same range
    x_b = min(x);
    y_b = min(y);
    x_m = range(x);
    y_m = range(y);
    m = max(x_m,y_m);

    x_warp = m*((x-x_b)/x_m) + x_b;
    y_warp = m*((y-y_b)/y_m) + y_b;

end

function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end