close all
clear all

%% Select the movie to work with

[filename,foldername] = uigetfile('*.tif'); % select the registered tif
imgData = tiffreadVolume([foldername,filename]);

%% play the movie
smooth_win = 10;

figure(1); clf
subplot(2,1,1)
a = image(imgData(:,:,1));

for i = 1:1e3
    a.CData = mean(imgData(:,:,i:i+smooth_win),3)*5;
    %b.CData = mean(imgData_raw(:,:,i:i+smooth_win),3);
    pause(1e-2)
    drawnow
end

%% draw mask
figure(1); clf; imagesc(mean(imgData,3)*10); axis equal tight; drawnow;
mask = roipoly();


%% process things
ft_type= 'movmean'; %the type of smoothing for fictrac data
ft_win = 10; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'gaussian','movmean'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {10,1};
n_centroid = 20;
f0_pct = 7;

tmp2 = dir([fileparts(foldername(1:end-1)),'\*ficTracData_DAQ.mat']);
load([tmp2.folder,'\',tmp2.name])
tmp2 = dir([fileparts(foldername(1:end-1)),'\*ficTracData_dat.mat']);
load([tmp2.folder,'\',tmp2.name])

all_data.ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
all_data.im = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct);
all_data.meta = fileparts(foldername(1:end-1));
all_data.ft.stims = ftData_DAQ.stim{1};

%% plot results
rho_thresh = .1;

figure(1); clf
subplot(3,1,1); hold on
plot(all_data.ft.xf,all_data.ft.f_speed)
plot(all_data.ft.xf,abs(all_data.ft.r_speed))
ylabel('speeds'); legend('forward','yaw')

subplot(3,1,2); hold on
a = plot(all_data.ft.xb,all_data.im.mu); a.YData(abs(diff(a.YData))>pi) = nan; 
a.YData(all_data.im.rho<rho_thresh) = nan;
a = plot(all_data.ft.xf,-all_data.ft.cue); a.YData(abs(diff(a.YData))>pi) = nan;
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})
set(gca,'YDir','reverse')
legend('mu','cue')

subplot(3,1,3)
imagesc(all_data.ft.xb,unwrap(all_data.im.alpha),all_data.im.z)

linkaxes(get(gcf,'Children'),'x')
axis tight
hold on
a = scatter(all_data.ft.xf,min(ylim)*ones(size(all_data.ft.xf)),'r','linewidth',5); a.YData(~all_data.ft.stims) = nan;


%% look at gain
lag = 10;

figure(2); clf
mu_tmp = all_data.im.mu;
mu_tmp(all_data.im.rho<rho_thresh) = nan;
mu_tmp = interp1(all_data.ft.xb,unwrap(mu_tmp),all_data.ft.xf,'linear','extrap');
cue_tmp = -unwrap(all_data.ft.cue);

mu_tmp = mu_tmp(1+lag:end);
cue_tmp = cue_tmp(1:end-lag);
xf_tmp = all_data.ft.xf(1:end-lag);

mu_speed = gradient(mu_tmp) / mean(diff(all_data.ft.xf));
cue_speed = gradient(cue_tmp) / mean(diff(all_data.ft.xf));
fly_speed = all_data.ft.r_speed(1+lag:end);

subplot(2,1,1); hold on
plot(xf_tmp,mu_tmp)
plot(xf_tmp,cue_tmp)

subplot(2,2,3); hold on
idx = abs(cue_speed) > .1 & ~isnan(mu_speed);
scatter(cue_speed(idx),mu_speed(idx),'.'); plot(xlim,[0,0],'k:'); plot([0,0],ylim,'k:')
xlabel('cue speed (rad/s)'); ylabel('mu speed (rad/s)')
g = cue_speed(idx) \ mu_speed(idx);
text(max(xlim),max(ylim),sprintf('gain: %.2f\nlag: %dms',g,round(1e3*lag*mean(diff(xf_tmp)))),...
        'HorizontalAlignment','right','VerticalAlignment','top')

subplot(2,2,4); hold on
idx = abs(fly_speed) > .1 & ~isnan(mu_speed);
scatter(fly_speed(idx),mu_speed(idx),'.'); plot(xlim,[0,0],'k:'); plot([0,0],ylim,'k:')
xlabel('fly speed (rad/s)'); ylabel('mu speed (rad/s)')
g = fly_speed(idx) \ mu_speed(idx);
text(max(xlim),max(ylim),sprintf('gain: %.2f\nlag: %dms',g,round(1e3*lag*mean(diff(xf_tmp)))),...
        'HorizontalAlignment','right','VerticalAlignment','top')

%% functions
function s = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct)
    
    imgData = smoothdata(imgData,3,im_type{1},im_win{1});
    imgData = imgData - min(imgData,[],'all');

    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length
min_axis        = min(range(x_mask),range(y_mask));
mid             = bwskel(mask,'MinBranchLength',min_axis);  %find the midline as the skeleton, shaving out all sub branches that are smaller than the minimum axis length
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!

xq          = [-min_axis:(length(x_mid)+min_axis)];                             %extend the midline so that it reaches the border of the mask. extrapolate as many points as the minimum axis length
x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');                 %keep only the points that exist within the mask
x_mid       = x_mid(idx);
y_mid       = y_mid(idx);

xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';                        %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                     %take every other so that we dont start at the edges, and all are same size
%assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)

    imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
    centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    alpha       = repmat(alpha,1,2);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});
    
    % imgData = squeeze(sum(imgData,3));
    % imgData = 256*(imgData-min(imgData,[],'all'))/(max(imgData,[],'all')-min(imgData,[],'all'));

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.d  = dff_cluster;
    s.f  = f_cluster;
    s.alpha = alpha;
    %s.imgData = imgData;
end

function s = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type)

    f_speed = ftData_dat.velFor{:};                       %store each speed
    f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}' / 192 * 2 * pi - pi;
    cue(abs(gradient(cue)) > 2) = nan;
    cue     = unwrap(cue);
    cue     = smoothdata(cue,1,ft_type,ft_win,'omitnan');
    cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
    cue(cue > pi) = cue(cue > pi) - 2*pi;
    
    s.xf      = seconds(ftData_DAQ.trialTime{:});
    if ismember('volClock',ftData_DAQ.Properties.VariableNames);
        s.xb      = seconds(ftData_DAQ.volClock{:});
    end
    s.f_speed = smoothdata(f_speed,1,ft_type,ft_win); 
    s.r_speed = smoothdata(r_speed,1,ft_type,ft_win);
    s.cue     = cue;
    
end