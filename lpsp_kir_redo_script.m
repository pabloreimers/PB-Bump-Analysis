%% set experiment params
calc_flag = true;
local_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\analysis\'; %define the local directory to clear up
remote_dir           = 'Z:\pablo\lpsp_kir\analysis\'; %define the remote directory to make sure things are backed up to

f0_pct              = 5;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = .5*10*5;                                              %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter. gaussian sigma is 1/5 the window length, so this is sigma*fr*5, in seconds
f_smooth            = .5*60*5;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 1;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = 0.1;

%% find paths to data
all_data = dir([remote_dir,'\**\denoised_regProduct.mat']);

%% clear data
n = length(all_data);

if calc_flag && ~exist('dff_tot','var')
dff_tot = cell(n,1); %preallocate all variables if we are storing new
dff_peak= cell(n,1);
dff_cluster = cell(n,1);
f_cluster = cell(n,1);
mu      = cell(n,1);
rho     = cell(n,1);
f_speed = cell(n,1);
r_speed = cell(n,1);
intHD   = cell(n,1);
cue     = cell(n,1);
r_vel   = cell(n,1);
xf      = cell(n,1);
cl_idx  = cell(n,1);
exp_idx = cell(n,1);
fly_num = cell(n,1);
exp_date= cell(n,1);
trial   = cell(n,1);
x_pos   = cell(n,1);
y_pos   = cell(n,1);
end

%% create masks
n = length(all_data);

for i = 1:length(all_data)
    fprintf('%s\n',all_data(i).folder)
    if isempty(dir([all_data(i).folder,'\denoised_mask.mat']))
        load([all_data(i).folder,'\',all_data(i).name]) %load in the denoised regProduct
        figure(1); clf
        imagesc(mean(sum(regProduct,4),3))
        drawnow
        tmp = drawellipse('Center',[size(regProduct,2)/2,size(regProduct,1)/2],'SemiAxes',[size(regProduct,2)/4,size(regProduct,1)/4]);
        input('') %move on once user presses enter, after adjusting the ellipse
        mask = createMask(tmp);
        tmp.SemiAxes = tmp.SemiAxes / 3;
        mask = logical(mask - createMask(tmp));

        save([all_data(i).folder,'\denoised_mask.mat'],'mask')
    else
        load([all_data(i).folder,'\denoised_mask.mat'])
        load([all_data(i).folder,'\',all_data(i).name]) %load in the denoised regProduct
    end

    if calc_flag
        tmp = split(all_data(i).folder,'\');
        fly_num{i} = tmp{6};
        exp_date{i} = tmp{5};
        trial{i} = tmp{7};
        cl_idx{i} = contains(tmp{7},'cl');
        exp_idx{i} = contains(tmp{7},'LPsP');

        tmp = dir([fileparts(all_data(i).folder),'\*fictracData_DAQ.mat']); %load fictrac data
        load([tmp.folder,'\',tmp.name])

        [f_speed{i},r_speed{i},intHD{i},cue{i},r_vel{i},x_pos{i},y_pos{i}] = ft_calc(ftData_DAQ,n_smooth,f_smooth);
        xf{i}  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed{i})]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
        if isduration(xf{i})
            xf{i} = seconds(xf{i});
        end
        [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i},f_cluster{i}] = bump_calc_eb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
        figure(3);
        subplot(3,1,1)
        a = plot(xf{i},-cue{i},'m');
        a.YData(abs(diff(a.YData))>pi) = nan;
        subplot(3,1,2)
        a = plot(xf{i},mu{i},'k');
        a.YData(abs(diff(a.YData))>pi) = nan;
        subplot(3,1,3)
        imagesc(dff_cluster{i})
        drawnow
    end
end

%% save data
data = table(dff_cluster,dff_tot,dff_peak,exp_date,exp_idx,f_cluster,f_speed,fly_num,mu,rho,r_speed,intHD,cue,r_vel,xf, cl_idx, trial,x_pos,y_pos);
save(['lpsp_kir_data_denoise_',num2str(yyyymmdd(datetime("today")))],'data')
%% Functions
function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function [amp_tot,amp_peak, mu, rho,dff_cluster,f_cluster] = bump_calc_eb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
regProduct = smoothdata(regProduct(:,:,1:end,:),4,'gaussian',b_smooth);

[tmp_y,tmp_x] = find(bwmorph(mask,'shrink','inf'));
if length(tmp_y) == 1
    [y_mask,x_mask] = find(mask);
    figure(1); clf; imagesc(mask);
    tmp = drawellipse('Center',[tmp_x,tmp_y],'SemiAxes',[range(x_mask)/6,range(y_mask)/6]);
    mask = logical(mask - createMask(tmp));
end

%extract midline from mask
[y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
mid             = bwmorph(mask,'remove');
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);
hold on; scatter(x0,y0,'b','filled')

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size

if centroids(1,1) < centroids(end,1)
    centroids = flipud(centroids);
end

%assign each pixel to a centroid
[~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)

imgData = squeeze(sum(regProduct,3));
tmp     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
background  = tmp(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
%flash_idx   = med_pix > 0.05;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
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

% n_planes = size(regProduct,3);
% filt_mu = linspace(max(centroids(:,1)),min(centroids(:,1)),n_planes);
% filt_sig= range(centroids(:,1))/(n_centroid/2);
% filt_x  = [1:size(regProduct,1)]';
% f_cluster = zeros(n_centroid,size(regProduct,4));
% 
% for p = 1:n_planes
% tmp = double(squeeze(regProduct(:,:,p,:))) .* normpdf(filt_x,filt_mu(p),filt_sig);
% 
% 
% imgData_2d      = reshape(double(tmp),[],size(regProduct,4));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
% centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
% for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
%     centroid_log(i, sub2ind(size(regProduct,1,2,4),y_mask(idx==i),x_mask(idx ==i))) = true;
% end
% 
% tmp       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
% f_cluster = f_cluster + tmp;
% end
% 
% f_cluster(:,flash_idx) = nan;
% f0              = prctile(f_cluster,f0_pct,2);            %find the baseline fluorescence in each cluster
% dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
% zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
% tmp = dff_cluster;
% tmp(:,~flash_idx) = zscore_cluster;
% dff_cluster = tmp;
alpha       = linspace(-pi,pi,n_centroid);

% for i = 1:n_smooth
%     dff_cluster = smoothdata(dff_cluster,2,'gaussian',b_smooth);
% end

% tmp = smoothdata(repmat(dff_cluster,3,1),1,'gaussian',3);
% tmp = tmp(n_centroid+1:end-n_centroid,:);
tmp = dff_cluster;
[x_tmp,y_tmp]   = pol2cart(alpha,tmp');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
% for i = 1:n_smooth
% mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
% rho    = smoothdata(rho,1,'gaussian',b_smooth);
% end

mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
amp_tot = mean(dff_cluster,1)';
[~,i] = mink(abs(alpha-mu),2,2); %find indexes corresponding to peak of each time point
i2 = i' + size(dff_cluster,1)*[0:size(dff_cluster,2)-1]; %find the linear index into the peak of each column (time point) value. this was clever :)
amp_peak = mean(dff_cluster(i2),1)'; %extract peak amplitude

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
amp_tot = smoothdata(amp_tot,1,'gaussian',b_smooth); 
amp_peak = smoothdata(amp_peak,1,'gaussian',b_smooth);
end

total_t = max(xf);
xb  = linspace(0,total_t,size(regProduct,4))';

mu          = interp1(xb,unwrap(mu),xf)';
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp_tot     = interp1(xb,amp_tot,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';
end

function [f_speed,r_speed,intHD,cue,r_vel,x_pos,y_pos] = ft_calc(ftData_DAQ,n_smooth,f_smooth)

f_speed = ftData_DAQ.velFor{:};                       %store each speed
r_speed = ftData_DAQ.velYaw{:};
intHD   = ftData_DAQ.intHD{:};
cue     = ftData_DAQ.cuePos{:}';


f_speed = f_speed;                                      %turn each velocity into a speed
r_speed = abs(r_speed);
intHD   = unwrap(intHD);                                %unwrap heading to perform circular smoothing. keeps radians continuous, so that smoothing 0 and 2pi doesnt go to 1pi
cue     = unwrap(cue / 192 * 2*pi - pi);
r_vel  = ftData_DAQ.velYaw{:};

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
f_speed = smoothdata(f_speed,1,'gaussian',f_smooth); 
r_speed = smoothdata(r_speed,1,'gaussian',f_smooth);
intHD   = smoothdata(intHD,  1,'gaussian',f_smooth);
cue     = smoothdata(cue,    1,'gaussian',f_smooth);
r_vel   = smoothdata(r_vel,  1,'gaussian',f_smooth);
end

intHD = mod(intHD,2*pi);                                %rewrap heading data, and put between -pi and pi.
intHD(intHD > pi) = intHD(intHD > pi) - 2*pi;
cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;

f_vel = ftData_DAQ.velFor{:};
s_vel = ftData_DAQ.velSide{:};
theta = ftData_DAQ.cueAngle{:} / 180 * pi;

x_pos = cumsum(f_vel.*cos(theta') + s_vel.*sin(theta'))/60;
y_pos = cumsum(f_vel.*sin(theta') + s_vel.*cos(theta'))/60;

% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
%    r_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
end
