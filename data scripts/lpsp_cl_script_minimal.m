%% clear all
clear all
close all

%% find path to all relevant files
base_dir = ('Z:\pablo\lpsp_cl\to analyze\20*\');
all_files = dir([base_dir,'\**\*imagingData*.mat']);

%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    if ~isfile([all_files(i).folder,'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        imgData = mean(regProduct,[3,4]);
        imgData = imgData - min(imgData,[],'all');
        imgData = imgData/max(imgData,[],'all')*256;
        figure(1); clf; imagesc(imgData); colormap(bone); drawnow;
        mask = roipoly(uint8(imgData));
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

%all_data = struct();

tic
for i = 84:length(all_files)
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end})
    load([all_files(i).folder,'\',all_files(i).name])
    load([all_files(i).folder,'\mask.mat'])
    tmp2 = dir([all_files(i).folder,'\*ficTracData_DAQ*']);
    load([tmp2.folder,'\',tmp2.name])

    all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
    all_data(i).im = process_im(regProduct, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% remove empties
empty_idx = nan(length(all_data),1);

for i = 1:length(all_data)
empty_idx(i) = isempty(all_data(i).ft);
end

all_data = all_data(~empty_idx);

%% linear model dff against r_speed
lag = 10; %apply a 10 frame lag in ft time, so 160ms, between fly speed and fluorescence
vel_thresh = .5;
bump_thresh = 10;
rho_thresh = .1;

lpsp_idx = nan(length(all_data),1);
dark_idx = nan(length(all_data),1);
r_R2 = nan(length(all_data),1);
f_R2 = nan(length(all_data),1);
j_R2 = nan(length(all_data),1);
vel_corr = nan(length(all_data),1);
fly_id = cell(length(all_data),1);


for i = 1:length(all_data)
    try
    lpsp_idx(i) = contains(all_data(i).meta,'_lpsp_','IgnoreCase',true);
    dark_idx(i) = contains(all_data(i).meta,'_dark','IgnoreCase',true);
    tmp = strsplit(all_data(i).meta,'\');
    tmp2 = regexp(tmp{end},'\d*','Match');
    fly_id{i} = [tmp2{1},'_',tmp2{end}];

    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = all_data(i).ft.r_speed;
    for_vel  = all_data(i).ft.f_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);
    amp      = interp1(xb,mean(all_data(i).im.d,1),xf);

    fly_vel  = fly_vel(1:end-lag);
    for_vel  = for_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);
    amp      = amp(lag+1:end);

    mdl = fitlm(abs(fly_vel),amp);
    r_R2(i)    = mdl.Rsquared.Adjusted;
    mdl = fitlm(for_vel,amp);
    f_R2(i)    = mdl.Rsquared.Adjusted;
    mdl = fitlm([abs(fly_vel),for_vel],amp);
    j_R2(i)    = mdl.Rsquared.Adjusted;

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;
    vel_corr(i)    = corr(fly_vel(idx),bump_vel(idx));

    end
end

[~,~,fly_num] = unique(fly_id);

%% plot example fly
for i = 1:length(all_data)
all_data(i).ft.xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),length(all_data(i).im.d));
end

i = 10;
figure(8); clf
subplot(3,1,1)
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue);
a.YData(abs(diff(a.YData))>pi) = nan;
yticks([-pi,0,pi])
yticklabels({'\-pi',0,'\pi'})
yticklabels({'-\pi',0,'\pi'})
ylim([-pi,pi])
set(gca,'Ydir','reverse')

subplot(3,1,2)
tmp = all_data(i).im.d;
tmp(:,sum(tmp,1)>10) = nan;
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp)

subplot(3,1,3)
a = plot(all_data(i).ft.xb,all_data(i).im.mu);
a.YData(abs(diff(a.YData))>pi) = nan;
yticks([-pi,0,pi])
yticklabels({'\-pi',0,'\pi'})
yticklabels({'-\pi',0,'\pi'})
ylim([-pi,pi])
set(gca,'Ydir','reverse')

linkaxes(get(gcf,'Children'),'x')

%% Plot results
group_order = {'EPG > GRAB(DA2m) (CL)','EPG > GRAB(DA2m) (dark)','LPsP > syt7f (CL)','LPsP > syt7f (dark)'};
ind = dark_idx + 2*lpsp_idx + 1;
c = hsv(4);

figure(1); clf
subplot(1,2,1)
X = [f_R2,r_R2,j_R2];
plot(X(ind==3,:)','.-k')
ylabel('R-squared')
xticks(1:3); xticklabels({'forward','rotational','joint'});xlim([.5,3.5])
title({'LPsP > syt7f','Closed Loop'})

subplot(1,2,2)
X = [f_R2,r_R2,j_R2];
plot(X(ind==4,:)','.-k')
ylabel('R-squared')
xticks(1:3); xticklabels({'forward','rotational','joint'});xlim([.5,3.5])
title({'LPsP > syt7f','Dark'})

linkaxes(get(gcf,'Children'),'y')
ylim([0,.5])
fontsize(gcf,20,'pixels')

figure(2)
subplot(2,2,1);
scatter(ind,vel_corr);
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Correlation Coefficent', '(fly vs bump vel)'})
hold on; plot(xlim,[0,0],'k:')

subplot(2,2,2)
scatter(ind,r_R2); hold on
scatter(ind+.1,f_R2)


subplot(2,1,2)
scatter(fly_num,vel_corr,[],c(ind,:))


%% assess how bump amplitude compare to calibration error

grab_idx = contains({all_data.meta},'EPG');

fr = mean(diff(all_data(i).ft.xf));
mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
mu_speed = gradient(mu) / fr;
dff = interp1(all_data(i).ft.xb,max(all_data(i).im.d,[],1),all_data(i).ft.xf);

figure(1); clf

for lag = 10
mu_speed_tmp = mu_speed(1:end-lag);
fly_speed_tmp= all_data(i).ft.r_speed(1+lag:end);
dff_tmp = dff(1:end-lag);

idx = abs(fly_speed_tmp) > .1;
scatter(fly_speed_tmp(idx),mu_speed_tmp(idx),[],dff_tmp(idx),'filled','MarkerFaceAlpha',.2)
xlabel('fly speed (rad/s)'); ylabel('bump speed (rad/s)')
pause(.1)
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


function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function s = process_im(regProduct, im_win, im_type, mask, n_centroid, f0_pct)
    
    imgData = squeeze(sum(smoothdata(regProduct,4,im_type{1},im_win{1}),3));
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

function out = resample_pablo(N,x)
    idx = randi(length(x),length(x),N);
    out = x(idx);
end


function cc = find_cc(lag,fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh)
    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;
    cc = corr(fly_vel(idx),bump_vel(idx));
end