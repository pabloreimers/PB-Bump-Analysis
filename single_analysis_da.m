function [vel_rho,vel_pval,pva_corr,pva_pval,lag,f_r2,r_r2,j_r2,num_flash,cue,mu,rho,vel_corr,disp_corr] = single_analysis_da(base_dir,images_flag)
%% script params
f0_pct              = 7;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 10;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 10;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 10;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;   
med_thresh          = 0.05;
pva_thresh          = 0;

%% load the relevant files for analysis from the input dir
name = ls([base_dir,'\*imagingData*']);
load([base_dir,'\',name],'regProduct');
name = ls([base_dir,'\*ficTrac*']);
load([base_dir,'\',name],'ftData_DAQ');
name = ls([base_dir,'\*mask*']);
load([base_dir,'\',name],'mask');

% find the flashes as points outside mask
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
background  = imgData(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = median(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + std(med_pix);
flash_idx   = med_pix > med_thresh;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',20));


figure(5); clf
plot(med_pix,'linewidth',2)
hold on
%plot([0,length(med_pix)],(median(med_pix) + std(med_pix))*[1,1],'linewidth',2)
plot([0,length(med_pix)], med_thresh*[1,1],'linewidth',2)
legend('Median Pixel Val, dF/F','Inclusion Threshold')
num_flash = sum(flash_idx);

%process the data to make image figures
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = 256*(imgData - min(imgData(:,:,~flash_idx),[],'all'))/(max(imgData(:,:,~flash_idx),[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
imgData2    = imgData;
top_int     = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
bot_int     = prctile(imgData2,5,'all');
imgData2    = max(min(imgData2,top_int),bot_int) - bot_int;
imgData2    = 256*imgData2/max(imgData2,[],'all');
imgData2    = smoothdata(imgData2,3,'movmean',avg_win);                      %smooth the data, again for viewing purposes (should this go before the clipping)            

figure(1); clf                                          % clear the current figure
image(uint8(mean(imgData2,3))); 
colormap('bone')
axis equal tight
xticks([]); yticks([])


tic
%% extract midline from mask
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

figure(2); clf
imagesc(max(imgData,[],3))                                      %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)
figure(2); clf
image(uint8(mean(imgData2,3)));                    %plot the image again with max intensity over time to show the whole pb
colormap(bone)
hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_mask(idx == i),y_mask(idx == i),'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
axis equal tight

%find the mean activity in each cluster
imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
end

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f0              = prctile(f_cluster(:,~flash_idx),f0_pct,2);            %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.

figure(3); clf
imagesc(subplot(2,2,1),centroid_log); colormap('bone')
xlabel('all pixels'); ylabel('cluster'); title('Centroid Logical')
imagesc(subplot(2,2,2),imgData_2d); colormap('bone')
xlabel('frames'); ylabel('all pixels'); title('Pixel Intensity (2D)')
imagesc(subplot(2,1,2),dff_cluster); colormap('bone'); pos = get(gca,'position'); colorbar; set(gca,'Position',pos)
xlabel('frame'); ylabel('cluster'); title('Grouped Intensity')

%% find pva estimate
alpha       = repmat(linspace(-pi,pi,n_centroid),1,2);

[x_tmp,y_tmp]   = pol2cart(alpha,dff_cluster');
[mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
amp             = mean(dff_cluster,1)';
amp_peak        = max(dff_cluster,[],1); %extract max amplitude at each time point

for i = 1:n_smooth
mu     = smoothdata(unwrap(mu),1,'gaussian',b_smooth); %smooth all bump parameters. this was done before in a cell array called for plotting. just do it do the actual variables now for ease of calling
rho    = smoothdata(rho,1,'gaussian',b_smooth);
amp    = smoothdata(amp,1,'gaussian',b_smooth);
amp_peak    = smoothdata(amp_peak,1,'gaussian',b_smooth);
end

mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;


%% process fictrac params
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
flash_idx   = logical(interp1(xb,double(flash_idx),xf));

%% find and plot locomotion correlations
fr      = mean(diff(xf)); %find the frame rate of the data

[f_corr,f_lag]  = xcorr(f_speed(~flash_idx),amp(~flash_idx),ceil(max_lag/fr),'coeff');      % find the max correlation, in steps
[r_corr,r_lag]  = xcorr(r_speed(~flash_idx),amp(~flash_idx),ceil(max_lag/fr),'coeff');      % find the max correlation, in steps
[~,l_idx]       = max(f_corr);
f_lag           = f_lag(l_idx);
[~,l_idx]       = max(r_corr);
r_lag           = r_lag(l_idx);

f_idx           = xb - f_lag*fr > 0 & xb - f_lag*fr < max(xb); %find the values that are not within the lag (in time)
r_idx           = xb - r_lag*fr > 0 & xb - r_lag*fr < max(xb);

f_speed(~f_idx) = 0;                                            %set those to zero before you shift things
r_speed(~r_idx) = 0;
 
f_speed_lag     = circshift(f_speed,-f_lag);                    %shift the traces to be aligned.
r_speed_lag     = circshift(r_speed,-r_lag);

amp(isnan(amp)) = 0;

[f_fit,f_gof]   = fit(f_speed_lag(~flash_idx),amp(~flash_idx),'poly1');
[r_fit,r_gof]   = fit(r_speed_lag(~flash_idx),amp(~flash_idx),'poly1');
[m_fit,m_gof]   = fit([r_speed_lag(~flash_idx),f_speed_lag(~flash_idx)],amp(~flash_idx),'poly11');


c1 = [1,0.5,0];
figure(6);clf
subplot(2,1,1)
plot(xf,f_speed_lag,'k','Linewidth',2); ylabel('Forward Velocity (mm/s)')
yyaxis right; a = plot(xf, amp,'Color',c1,'linewidth',2); ylabel('Average \DeltaF/F'); ax = gca; ax.YAxis(2).Color = c1;
a.YData(flash_idx) = nan;
axis tight
y = ylim; text(xf(end),y(2),sprintf('r^2: %.2f\n',...
                f_gof.adjrsquare),...
                'HorizontalAlignment','right','VerticalAlignment','top')

subplot(2,1,2)
plot(xf,r_speed_lag,'k','linewidth',2); ylabel('Rotational Speed (rad/s)')
yyaxis right; a = plot(xf, amp,'Color',c1,'linewidth',2); ylabel('Average \DeltaF/F'); ax = gca; ax.YAxis(2).Color = c1;
a.YData(flash_idx) = nan;
xlabel('time (s)')
axis tight
y = ylim; text(xf(end),y(2),sprintf('r^2: %.2f\njoint r^2: %.2f',...
                                r_gof.adjrsquare,m_gof.adjrsquare),...
                'HorizontalAlignment','right','VerticalAlignment','top')

f_r2 = f_gof.adjrsquare;
r_r2 = r_gof.adjrsquare;
j_r2 = m_gof.adjrsquare;
%% plot fly heading and bump heading, from PVA
%lag = ceil(fmincon(@(x)(circ_corrcc(cue(1:end-ceil(x)),-mu((ceil(x)+1):end))),20,[],[],[],[],0));
lag = 20;
pva_idx = rho > pva_thresh;
cue_lag = cue(1:end-ceil(lag));
mu_lag  = mu((ceil(lag)+1):end);
[pva_corr,pva_pval] = circ_corrcc(cue_lag(pva_idx((ceil(lag)+1):end)),-mu_lag(pva_idx(ceil(lag)+1:end)));

figure(7); clf
subplot(3,1,1)
a = plot(xf,cue,'k','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'}); ylim([-pi,pi])
ylabel('Fly Heading (cue)')

subplot(3,1,2)
scatter(xf(~flash_idx & pva_idx'),-mu(~flash_idx & pva_idx'),[],rho(~flash_idx & pva_idx'),'.')
colormap('bone')
yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'}); ylim([-pi,pi])
pos = get(gca,'Position');
colorbar
set(gca,'Position',pos)
ylabel('PVA')

subplot(3,1,3)
a = plot(xf,circ_dist(cue,-mu),'k','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
a.YData(flash_idx | ~pva_idx') = nan;
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
flash_lag= flash_idx(lag+1:end);
pva_lag  = pva_idx(lag+1:end);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
fly_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
end
vel_idx = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min; %ignore outlier bump speeds with arbitrary threshold

figure(9); clf
scatter(bump_vel(vel_idx & ~flash_lag' & pva_lag),fly_vel(vel_idx & ~flash_lag' & pva_lag),5,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.1)
xlabel('Bump Vel (rad/s)'); ylabel('Fly Rot Vel (rad/s)');
hold on
[vel_rho,vel_pval] = corr(bump_vel(vel_idx & ~flash_lag' & pva_lag),fly_vel(vel_idx & ~flash_lag' & pva_lag));
x = xlim;
y = ylim;
text(x(2),y(2),sprintf('$r = %.2f$\n$p < 10^{%i}$',vel_rho,ceil(log10(vel_pval))),'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
plot([0,0],y,'k:')
plot(x,[0,0],':k')

lag = lag*fr;

%% show fluorescence over heading
tmp_data = dff_cluster;
tmp_idx = med_pix > median(med_pix) + 2*std(med_pix);
tmp_data(:,tmp_idx) = nan;

figure(4); clf
subplot(3,1,1)
imagesc(xb,1:2*n_centroid,tmp_data); xticks([]); ylabel('cluster'); title('\DeltaF/F0'); pos = get(gca,'Position'); colorbar; set(gca,'Position',pos); colormap(bone)
h(2) = subplot(3,1,2);
tmp = tmp_data ./ (sum(tmp_data,1));
top_val = prctile(tmp(:),95);
bot_val = prctile(tmp(:),5);
tmp( tmp > top_val) = top_val;
tmp( tmp < bot_val) = bot_val;
imagesc(xb,1:2*n_centroid,tmp); xticks([]); ylabel('cluster'); title('\DeltaF/(F0*sum(F_t))'); pos = get(gca,'Position'); colorbar; set(gca,'Position',pos); colormap(bone)
h(3) = subplot(3,1,3); 
tmp = cue; tmp(abs(diff(tmp)) > pi) = nan; plot(xf,tmp);
axis tight; xlabel('time(s)'); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'}); ylim([-pi,pi]); ylabel('heading')

linkaxes(h,'x')

%% Plot peak fluorescence against bump disp vs fly disp
lag = 20;
bump_vel = [diff(unwrap(mu));0] / fr; %divide by the frame rate (seconds per frame). since the mu is in units of radians, this should give rad per seconds
fly_vel  = ftData_DAQ.velYaw{:};
bump_vel = bump_vel(lag+1:end);
fly_vel  = fly_vel(1:end-lag);
tmp = amp_peak(lag+1:end);

for i = 1:n_smooth                                      %smooth fictrac data n times, since fictrac is super noisy.
fly_vel = smoothdata(fly_vel,1,'gaussian',f_smooth);
end
vel_idx = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min; %ignore outlier bump speeds with arbitrary threshold

figure(13); clf
scatter(bump_vel(vel_idx) ./ fly_vel(vel_idx),tmp(vel_idx),[],fly_vel(vel_idx),'filled','MarkerFaceAlpha',0.2)
xlabel('bump vel / fly vel')
ylabel('peak \DeltaF/F')
tmp = corrcoef(bump_vel(idx)./fly_vel(idx),tmp(idx));
vel_corr = tmp(2);
colormap(parula)

figure(14); clf
win = ceil(1/fr);
fly_disp  = abs(circ_dist(cue(1:end-win),cue(win+1:end)));
fly_disp  = fly_disp(lag+1:end);
bump_disp = abs(circ_dist(mu(1:end-win),mu(win+1:end)));
bump_disp = bump_disp(1:end-lag);

%tmp = max(interp1(xb,smoothdata(dff_cluster,1,'gaussian',n_centroid/2)',xf)',[],1);
tmp = amp_peak;
tmp = tmp(win+1:end);
tmp = tmp(1:end-lag);
idx = fly_disp > 0.05;
scatter(log(fly_disp(idx)./bump_disp(idx)),tmp(idx),[],'filled','MarkerFaceAlpha',.1)
xlabel('$ln(\frac{\Delta_{bump}}{\Delta_{fly}}) (1s)$','Interpreter','latex','Fontsize',50)
ylabel('$max(\frac{\Delta F}{F0})$','Interpreter','latex','Fontsize',50)
tmp = corrcoef(fly_disp(idx)./bump_disp(idx),tmp(idx));
disp_corr = tmp(2);


%% save all figures
if images_flag
% saveas(figure(1),[base_dir,'\mean_image.png'])
% saveas(figure(2),[base_dir,'\ROI.png'])
% saveas(figure(3),[base_dir,'\dff_calc.png'])
% saveas(figure(4),[base_dir,'\fluorescence.png'])
% saveas(figure(5),[base_dir,'\flash.png'])
% saveas(figure(6),[base_dir,'\locomotion.png'])
% saveas(figure(7),[base_dir,'\pva.png'])
% saveas(figure(9),[base_dir,'\velocities.png'])
saveas(figure(13),[base_dir,'\velocity_gain.png'])
saveas(figure(14),[base_dir,'\displacement_gain.png'])
end
toc
end

