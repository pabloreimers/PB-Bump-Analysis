%% set script params
clear all
close all
mask_flag           = true;
mask_overwrite      = true;
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
rho_thresh          = 1e-2;

%% select folder of images
base_dir = uigetdir(data_dir);
trials   = dir(base_dir);
trials(1:2) = [];
trials   = {trials.name}';
tmp = regexp( trials, '-\d+_', 'match', 'once' );
[~,reindex] = sort(-cellfun(@(x)sscanf(x,'%d'),tmp));
trials = trials(reindex);

%% extract gains
tmp   = cellfun(@(x)(regexp(x,'_','split')),trials,'UniformOutput',false);
gains = cellfun(@(x)(str2double(x{end})),tmp);
dark  = cellfun(@(x)(contains(x,'dark')),trials);

%% open each trial, draw mask, and save
if mask_flag
for i = 1:length(trials)
    if isempty(dir([base_dir,'\',trials{i},'\*mask*'])) || mask_overwrite
    tmp = dir([base_dir,'\',trials{i},'\registration*\*imagingData*']);
    load([tmp.folder,'\',tmp.name])
    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    imgData2    = imgData;
    top_int     = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
    bot_int     = prctile(imgData2,5,'all');
    imgData2    = max(min(imgData2,top_int),bot_int) - bot_int;
    imgData2    = 256*imgData2/max(imgData2,[],'all');
       
    
    figure(1); clf                                          % clear the current figure
    mask = roipoly(uint8(mean(imgData2,3))); 
    save([base_dir,'\',trials{i},'\mask.mat'],'mask')
    end
end
end
    
%% store the fictrac, total amplitude, and peak amplitude for each trial
n       = length(trials);
dff_tot = {n,1}; %preallocate dff cell
dff_peak= {n,1};
mu      = {n,1};
rho     = {n,1};
f_speed = {n,1};
r_speed = {n,1};
intHD   = {n,1};
cue     = {n,1};
r_vel   = {n,1};
xf      = {n,1};

for i = 1:length(trials)
    disp(i)

    tmp = dir([base_dir,'\',trials{i},'\registration*\*imagingData*']);
    load([tmp.folder,'\',tmp.name])
    tmp = dir([base_dir,'\',trials{i},'\*ficTracData_DAQ*']);
    load([tmp.folder,'\',tmp.name])
    tmp = dir([base_dir,'\',trials{i},'\*mask*']);
    load([tmp.folder,'\',tmp.name])
    
    [f_speed{i},r_speed{i},intHD{i},cue{i},r_vel{i}] = ft_calc(ftData_DAQ,n_smooth,f_smooth);
    xf{i}  = mean(diff(ftData_DAQ.trialTime{:})) * [1:length(f_speed{i})]; %create a time vector for the fictrac data. have to remake because of the error frames                             %create vectors with timestamps for each trace
    if isduration(xf{i})
        xf{i} = seconds(xf{i});
    end
    [dff_tot{i},dff_peak{i},mu{i},rho{i}] = bump_calc(mask,squeeze(sum(regProduct,3)),n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
    
end

%% verify cue gain

gains_set = unique(gains);

re_cue    = {length(trials),length(gains_set)}; %store the cue position for each trial for each of the tested gains

for i = 1:length(gains_set)
    for j = 1:length(trials)
         tmp = cumsum(r_vel{j} * gains_set(i) * mean(diff(xf{i})));
         for k = 1:10
             tmp = smoothdata(tmp,'gaussian',50);
         end
         re_cue{j,i} = tmp;

    end
end

%% plotunwrapped cues over each other
figure(2); clf
str = {size(gains_set)};
for j = 1:length(trials)
    subplot(length(trials),2,2*j - 1)
    hold on
    for i = 1:length(gains_set)
        plot(xf{j},(unwrap(re_cue{j,i}) - re_cue{j,i}(1)) ,'linewidth',2)
        str{i} = sprintf('gain = %.1f',gains_set(i));
    end
    plot(xf{j},unwrap(mu{j}) - mu{j}(1),'linewidth',2)
    
    pos = get(gca,'Position');
    legend(str{:},'mu','location','northeast outside')
    set(gca,'Position',pos)
    axis tight
    ylabel(sprintf('Accumulated\nRotation\n(radians)'))
    ylabel('$\int\hat{\theta}$','Interpreter','latex')
    if dark(j)
        tmp = 'dark';
    else
        tmp = 'cue';
    end
    y = ylim;
    text(0,max(y),sprintf('trained with gain = %.1f, %s',gains(j),tmp),'VerticalAlignment','top')
    fontsize(gca,.25,'inches')
end


%% running error calculation with each gain
% bin_time = 2; %set the bin time to be one second
% lag      = 1/3;
% errs = {size(re_cue)};
% gain_inf = {size(re_cue)};
% 
% for j = 1:length(trials)
%     fr = mean(diff(xf{j}));
%     bin_size = floor(bin_time / fr);
%     lag_frames = floor(lag / fr);
% 
%     for i = 1:length(gains_set)
%         tot_err = 0;
%         errs{j,i} = nan((length(re_cue{j,i})-bin_size-lag_frames),1);
%         gain_inf{j,i} = nan((length(re_cue{j,i})-bin_size-lag_frames),1);
%         for k = 1:(length(re_cue{j,i})-bin_size-lag_frames)
%             fprintf('j:%i, i:%i, k:%.2f\n',j, i,k/(length(re_cue{j,i})-bin_size-lag_frames))
%             mu_tmp  = unwrap(mu{j}((k:k+bin_size) + lag_frames));
%             mu_tmp  = mu_tmp - mu_tmp(1);
%             cue_tmp = unwrap(re_cue{j,i}(k:k+bin_size));
%             cue_tmp = cue_tmp - cue_tmp(1);
%             errs{j,i}(k) = sum(abs(mu_tmp - cue_tmp));% / sum(r_speed{j}((k:k+bin_size) + lag_frames));% / sum(abs(diff(cue_tmp))); %find the positional error in the current binned window. sum up all the error in the window. divide by the frame rate to keep this a consistent metric.
%             gains_inf{j,i}(k) = ((mu{j}(k + lag_frames) - mu{j}(k+1+lag_frames))) * fr / r_vel{j}(k);
%         end
%         %errs{j,i} = ; %total error is all the summed error divided by however many windows were taken
%     end
% end

%% find gain as a running metric
bin_time = 10; %set the bin time to be one second
lag_init = 1/3;
gain_init= .7;
gain_inf = cell(length(trials),1);
lag_inf  = cell(length(trials),1);
x0 = [0.7, lag_init, 0];
lb = [0,    0, -5];
ub = [2,    1,  5];

for j = 1:length(trials)
    fr = mean(diff(xf{j}));
    bin_size = floor(bin_time / fr);
    %lag_frames = floor(lag / fr);
    num_bins = 3*floor(length(mu{j}) / bin_size);
    gain_inf{j} = nan(num_bins-3,1);

    for k = 1:(num_bins-3)
        start_idx = (k-1)*(bin_size/3) + 1;
        obj_fun = @(x)bump_error(x,mu{j}(start_idx:start_idx+bin_size),r_vel{j}(start_idx:start_idx+bin_size),fr);
        x  = fmincon(obj_fun,x0,[],[],[],[],lb,ub);
        gain_inf{j}(k) = x(1);
        lag_inf{j}(k) = x(2);
    end
end
%% plot inferred gains
figure(2); subplot(2,2,2); cla
yyaxis right
ax = gca;
ax.YAxis(2).Color = 'k';
ax.YAxis(1).Color = 'none';
hold on
str = {};
for i = 1:length(gain_inf)
    if dark(i)
        c = [0,0,0];
    else
        c = [0,.5,0];
    end
    swarmchart(i*ones(length(gain_inf{i})),gain_inf{i},'filled','XJitterWidth',.5,'MarkerFaceColor',c,'MarkerFaceAlpha',.5)
    errorbar(i,mean(gain_inf{i}),std(gain_inf{i})/sqrt(length(gain_inf{i})),'.r')
    str{i} = sprintf('Trial: %d (%.1f)',i,gains(i));
end
xticks([1:length(gain_inf)])
xticklabels(str)
ylabel('inferred gain')
% 
% idx = find(dark);
% subplot(2,2,4); cla
% c = [1,0.5,0; 0,0.5,1; .75, 0, .75];
% yyaxis right

% ax = gca;
% ax.YAxis(2).Color = 'k';
% ax.YAxis(1).Color = 'none';
% ax.ColorOrder = c
% x  = 0:.1:2;
% y  = nan(length(x)-1,length(trials));
% str = {};
% for i = 1:length(trials)
%     y(:,i) = histcounts(gain_inf{i},x)/length(gain_inf{i});
%     str{i} = sprintf('Trial: %d, Trained: %.1f',i,gains(i));
% end
% 
% bar(x(1:end-1),y(:,dark))
% legend(str(dark))
% xlabel('gain'); ylabel('proportion')

%% find pvalues of histograms being different (paired ttest2)
h = nan(length(trials));
p = nan(length(trials));

for i = 1:length(trials)
    for j = 1:length(trials)
        [h(i,j),p(i,j)] = ttest2(gain_inf{i},gain_inf{j},'Vartype','unequal');
    end
end

bonf = (numel(p) - length(trials)) / 2;
figure(2); subplot(2,2,4)
heatmap(p * bonf)
title('paired t-test p-values (bonferroni correction)')
%% plot errors with traces
% figure(2);
% for j = 1:length(trials)
%         subplot(length(trials),2,2*j); cla
%         hold on
%     for i = 1:length(gains_set)
%         plot(errs{j,i},'linewidth',.5); %plot(medfilt1(gains_inf{j,i}),'linewidth',.5)
%     end
% end
%         %heatmap(errs,'XData',gains_set)
% subplot(length(trials),2,2)
% title(sprintf('Rolling Positional Error\nbin = %.1f seconds',bin_time),'Fontsize',20)
% subplot(length(trials),2,2*length(trials))
% xlabel('gains','Fontsize',72*.5)
% 
% %% show cue_gain vs mu_gain for each trial
% c       = [0,0.5,0; 0,0,0];
% 
% figure(8); clf
% hold on
% for i = 1:length(trials)
%     fr      = mean(diff(xf{i})); %find the frame rate of the data
%     vel = smoothdata(r_vel{i},1,'gaussian',f_smooth);
%     idx = abs(vel) > vel_min;
%     
% 
%     tmp = [diff(unwrap(-cue{i}))/fr;0];
% 
%     cue_gain = vel(idx) \ tmp(idx);
% 
%     tmp = [diff(unwrap(mu{i}))/fr;0];
%     mu_gain  = vel(idx) \ tmp(idx);
% 
%     subplot(1,2,1)
%     hold on
%     scatter(gains(i),mu_gain,100,'filled','MarkerFaceColor',c(dark(i)+1,:))
% 
%     subplot(length(trials),2,2*i)
%     hold on
%     a = plot(xf{i},-cue{i},'Color',c(dark(i)+1,:),'linewidth',2);
%     a.YData(abs(diff(a.YData))>pi) = nan;
% 
%     tmp = unwrap(mu{i}) + median(circ_dist(-cue{i},mu{i}));
%     tmp   = mod(tmp,2*pi);                                %rewrap heading data, and put between -pi and pi.
%     tmp(tmp > pi) = tmp(tmp > pi) - 2*pi;
% 
%     a = plot(xf{i},tmp,'Color', [1,.5,0],'linewidth',2);
%     a.YData(abs(diff(a.YData))>pi) = nan;
%     a.YData(rho{i}<rho_thresh) = nan;
%     axis tight
%     yticks([-pi,-0,pi]); yticklabels({'-\pi','0','\pi'});ylim([-pi,pi]);ylabel(sprintf('gain = %.2f',gains(i)))
% end
% 
% subplot(1,2,1)
% xlabel('cue gain')
% ylabel('mu gain')
% axis equal
% x = xlim;
% plot([x(1),x(2)],[x(1),x(2)],':k')
% %% plot distribution of rotational speed and gains
% figure(1); clf
% for i = 1:n
% 
%     idx = r_speed{i} > vel_min;
% subplot(3,2,1)
% hold on
% swarmchart(gains(i)*ones(length(dff_peak{i}),1),r_speed{i},'filled','MarkerFaceAlpha',.2)
% ylabel('r speed')
% subplot(3,2,3)
% hold on
% swarmchart(gains(i)*ones(sum(idx),1),dff_peak{i}(idx),'filled','MarkerFaceAlpha',.2)
% ylabel('dff peak')
% subplot(3,2,5)
% hold on
% swarmchart(gains(i)*ones(sum(idx),1),dff_peak{i}(idx)./r_speed{i}(idx),'filled','MarkerFaceAlpha',.2)
% ylabel('dff peak / r speed')
% xlabel('gain')
% 
% subplot(3,2,2)
% hold on
% swarmchart(gains(i)*ones(length(dff_peak{i}),1),r_speed{i},'filled','MarkerFaceAlpha',.2)
% ylabel('r speed')
% subplot(3,2,4)
% hold on
% swarmchart(gains(i)*ones(sum(idx),1),dff_tot{i}(idx),'filled','MarkerFaceAlpha',.2)
% ylabel('dff tot')
% subplot(3,2,6)
% hold on
% swarmchart(gains(i)*ones(sum(idx),1),dff_tot{i}(idx)./r_speed{i}(idx),'filled','MarkerFaceAlpha',.2)
% ylabel('dff tot / r speed')
% xlabel('gain')
% end
% subplot(3,2,1); plot(xlim,vel_min*[1,1],'k','linewidth',2)
% subplot(3,2,2); plot(xlim,vel_min*[1,1],'k','linewidth',2)
% %% bootstrap a bunch of timepoints and compare dff to rspeed ratio
% N       = 1e4;
% frames  = length(dff_peak{i});
% idx     = randi(frames,N,frames);
% p_peak  = nan(length(trials));
% p_tot   = nan(length(trials));
% p_speed = nan(length(trials));
% 
% for i = 1:length(trials)
%     for j= 1:length(trials)
%         p_peak(i,j) = sum(mean(dff_peak{i}(idx)./r_speed{i}(idx),2) > mean(dff_peak{j}(idx)./r_speed{j}(idx),2))/N;
%         p_tot(i,j)  = sum(mean(dff_tot{i}(idx)./r_speed{i}(idx),2) > mean(dff_tot{j}(idx)./r_speed{j}(idx),2))/N;
%         p_speed(i,j) = sum(mean(r_speed{i}(idx),2) > mean(r_speed{j}(idx),2))/N;
%     end
% end
% 
% %%
% N = 1e2;
% m_peak = nan(N,length(trials));
% c_peak = nan(N,length(trials));
% m_tot  = nan(N,length(trials));
% c_tot  = nan(N,length(trials));
% 
% for i = 1:length(trials)
%     tmp_idx = r_speed{i} > .2;
%     tmp_speed = r_speed{i}(tmp_idx);
%     tmp_peak  = dff_peak{i}(tmp_idx);
%     tmp_tot   = dff_tot{i}(tmp_idx);
% 
%     tmp_idx = randi(length(tmp_speed),N,length(tmp_speed));
% 
%     for j = 1:N
%         sprintf('i: %i, j: %i',i,j)
%         tmp = fitlm(tmp_speed(tmp_idx(j,:)),tmp_peak(tmp_idx(j,:)));
%     
%         m_peak(j,i) = tmp.Coefficients.Estimate(2);
%         c_peak(j,i) = tmp.Coefficients.Estimate(1);
% 
%         tmp = fitlm(tmp_speed(tmp_idx(j,:)),tmp_tot(tmp_idx(j,:)));
%     
%         m_tot(j,i)  = tmp.Coefficients.Estimate(2);
%         c_tot(j,i)  = tmp.Coefficients.Estimate(1);
%     end
% end
% 
% %% 
% c = [1, 0.5, 0; ...
%     0, 0.5, 1; ...
%     .75, 0, 0.75;...
%     .75, .75, 0; ...
%     0, .75, .75; ...
%     0, 0, 0];
% %c = jet(6)
% 
% figure(4); clf
% for i = 1:length(gains)
%     subplot(1,2,1)
%     hold on
%     swarmchart(gains(i)*ones(size(m_tot(:,i))),m_tot(:,i),[],c(i,:),'filled')
%     title('Total vs Speed (slope)')
%     xlabel('gain')
% 
%     subplot(1,2,2)
%     hold on
%     swarmchart(gains(i)*ones(size(m_peak(:,i))),m_peak(:,i),[],c(i,:),'filled')
%     title('Peak vs Speed (slope)')
%     xlabel('gain')
% end
% 
% %% show p values in heatmap
% subplot(2,2,1)
% heatmap(p_speed)
% title('r speed')
% 
% subplot(2,2,3)
% heatmap(p_peak)
% title('dff peak')
% 
% subplot(2,2,4)
% heatmap(p_tot)
% title('dff tot')
% 
% %% plot trajectories
% figure(2); clf
% for i = 1:length(trials)
%     subplot(length(trials),1,i);
%     a = plot(xf,cue{i});
%     a.YData(abs(diff(a.YData))>pi) = nan;
%     yticks([-pi,0,pi]);yticklabels({'-\pi','0','\pi'});
%     axis tight; ylim([-pi,pi])
%     text(max(xf),pi,sprintf('gain = %.1f',gains(i)),'HorizontalAlignment','right','VerticalAlignment','top')
%     ylabel(num2str(i))
% end


%% functions
function mse = bump_error(x,mu,r_vel,fr)
    gain = x(1);
    lag  = x(2);
    offset=x(3);
    lag_frames = max(floor(lag / fr),1);
    cue = cumsum(r_vel * gain * fr);
    for k = 1:10
        cue = smoothdata(cue,'gaussian',50);
    end
    mu_tmp = unwrap(mu(lag_frames:end));
    mu_tmp = mu_tmp - mu_tmp(1);
    cue_tmp= cue(1:end-lag_frames+1) - cue(1);
    mse = mean((mu_tmp - cue_tmp+offset).^2);     

end
function [f_speed,r_speed,intHD,cue,r_vel] = ft_calc(ftData_DAQ,n_smooth,f_smooth)

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

r_vel  = ftData_DAQ.velYaw{:};
% for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
%    r_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
% end
end

function [amp_tot,amp_peak, mu, rho] = bump_calc(mask, imgData, n_centroid, f0_pct, n_smooth, b_smooth, xf)
%extract midline from mask
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
centroids   = centroids(2:2:end-1,:);

%assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)

%find the mean activity in each cluster
imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f0              = prctile(f_cluster,f0_pct,2);            %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.

alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);

[x_tmp,y_tmp]   = pol2cart(alpha,dff_cluster');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
for i = 1:n_smooth
mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
rho    = smoothdata(rho,1,'gaussian',b_smooth);
end

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
xb  = linspace(0,total_t,size(imgData,3))';

mu          = interp1(xb,unwrap(mu),xf)';
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp_tot     = interp1(xb,amp_tot,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';
end