%% set script params
clear all
close all
mask_flag           = true;
mask_overwrite      = false;
movie_flag          = false;                                                %set whether to play movies or just run straight through the script
pause_time          = 0;                                                    %set the pause time if playing movies
data_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %set main data directory for ease of selecting files
%data_dir            = 'C:\Users\preim\Documents\Wilson Lab\data\2p data\';
f0_pct              = 7;                                                    %set the percentile that for the baseline fluorescence
n_centroid          = 20;                                                   %this is how many centroids per hemisphere. centroids = bins = glomeruli, basically
b_smooth            = 12;                                                   %define how many frames to smooth 2p data. both for bump parameters, and for fluorescence. gaussian filter.
f_smooth            = 50;                                                   %set how many frames to smooth for fictrac. gaussian, and repeated n times because very noisy
n_smooth            = 1;                                                   %set how many times to perform gaussian smoothing on fictrac
max_lag             = 1e3;                                                  %max lag to search for optimal cross correlation between dF/F and kinematics, in seconds
cluster_flag        = 0;                                                    %find dff by cluster or for whole pb (do we see a bump)
avg_win             = 5;                                                    %when showing raw 2p data, set window averaging size
vm_thresh           = 0;
var_thresh          = 0;
vel_thresh          = 10;                                                   %exclude points in bump to fly vel correlation that are faster than 10rad/s
vel_min             = 1e-1;                                                 %exclude points in bump to fly vel correlation where fly is slower than .01rad/s (effectively just fictrac noise)
rho_thresh          = 2e-4;
eb_flag             = true;
pb_flag             = false;
tmp_labels = {'CL = 0.7','Dark (0.7)','Ramp = 0.7 -> 0.3','Dark (0.3)','Ramp 0.3 -> 0.8','Dark (0.8)'};

%% select folder of images
base_dir = uigetdir(data_dir);
trials   = dir(base_dir);
trials(1:2) = [];
trials   = {trials.name}';
tmp = regexp( trials, '-\d+_', 'match', 'once' );
[~,reindex] = sort(-cellfun(@(x)sscanf(x,'%d'),tmp));
trials = trials(reindex);

%% open each trial, draw mask, and save
if mask_flag
for i = 1:length(trials)
    if isempty(dir([base_dir,'\',trials{i},'\*mask*'])) || mask_overwrite
    tmp = dir([base_dir,'\',trials{i},'\registration*\*imagingData*']);
    tmp = dir([base_dir,'\',trials{i},'\*imagingData*']);
    load([tmp.folder,'\',tmp.name])
    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    imgData2    = imgData;
    top_int     = prctile(imgData2,99,'all');                                    %clip the extremes and renormalize for viewing
    bot_int     = prctile(imgData2,5,'all');
    imgData2    = max(min(imgData2,top_int),bot_int) - bot_int;
    imgData2    = 256*imgData2/max(imgData2,[],'all');
       
    
    figure(1); clf                                          % clear the current figure
    tmp = mean(imgData,3);
    tmp = 256*(tmp - min(tmp(:)))/(max(tmp(:) - min(tmp(:))));
    image(tmp); colormap(bone); drawnow;
    if eb_flag
        tmp = drawellipse('Center',[size(tmp,2)/2,size(tmp,1)/2],'SemiAxes',[size(tmp,2)/4,size(tmp,1)/4]);
        input('') %move on once user presses enter, after adjusting the ellipse
        mask = createMask(tmp);
    end
    if pb_flag
        mask = roipoly(uint8(mean(imgData2,3)));  
    end
    save([base_dir,'\',trials{i},'\mask.mat'],'mask')
    end
end
end
    
%% store the fictrac, total amplitude, and peak amplitude for each trial
n       = length(trials);
dff_tot = {n,1}; %preallocate dff cell
dff_peak= {n,1};
dff_cluster = {n,1};
mu      = {n,1};
rho     = {n,1};
f_speed = {n,1};
r_speed = {n,1};
intHD   = {n,1};
cue     = {n,1};
r_vel   = {n,1};
xf      = {n,1};


figure(2); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
tiledlayout(rows,cols)

for i = 1:length(trials)
    disp(i)

    tmp = dir([base_dir,'\',trials{i},'\registration*\*imagingData*']);
    %tmp = dir([base_dir,'\',trials{i},'\*imagingData*']);
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
    if eb_flag
    [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i}] = bump_calc_eb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
    end
    if pb_flag
        [dff_tot{i},dff_peak{i},mu{i},rho{i},dff_cluster{i}] = bump_calc_pb(mask,regProduct,n_centroid,f0_pct,n_smooth,b_smooth,xf{i});
    end

    nexttile
    imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
    imgData     = 256*(imgData - min(imgData,[],'all'))/(max(imgData,[],'all') - min(imgData,[],'all')); %linearly rescale the scanimage data so that it can be shown as an image (without scaling the image for each plane, which can change baseline brightness in the visuals)
    tmp = mean(imgData,3);
    tmp = 256*(tmp - min(tmp(:)))/(max(tmp(:) - min(tmp(:))));
    image(tmp); colormap('bone'); axis equal tight; title(tmp_labels{i},'color','w')
    set(gca,'color','none')
end
set(gcf,'color','none')

%% find the correlation between fly vel and bump vel in each trial
lag = 10;
figure(9); clf
rows = floor(sqrt(n));
cols = ceil(n/rows);
for i = 1:n
    rho_thresh = prctile(rho{i},20);
    fr      = mean(diff(xf{i})); %find the frame rate of the data
    bump_vel = [diff(unwrap(mu{i}));0] / fr;
    bump_vel = bump_vel(lag+1:end);
    fly_vel  = r_vel{i}(1:end-lag);
    vel_idx  = abs(fly_vel) < vel_thresh & abs(bump_vel) < vel_thresh & abs(fly_vel) > vel_min & rho{i}(lag+1:end) > rho_thresh; %ignore outlier bump speeds with arbitrary threshold
    subplot(rows,cols,i)
    scatter(fly_vel(vel_idx),bump_vel(vel_idx),5,'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.1)
    xlabel('Fly Vel (rad/s)'); ylabel('Bump Vel (rad/s)');
    hold on
    [vel_rho,vel_pval] = corr(bump_vel(vel_idx),fly_vel(vel_idx));
    x = xlim;
    y = ylim;
    b = bump_vel(vel_idx) \ [ones(sum(vel_idx),1),fly_vel(vel_idx)]; %fit the slope of fly vel and bump vel with an arbitrary offset
    text(x(2),y(2),sprintf('$r = %.2f$\n$gain = %.2f$\n$lag = %.0f ms$',vel_rho,b(2),lag*fr*1000),'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','top')
    plot([0,0],y,'k:')
    plot(x,[0,0],':k')
    plot(fly_vel(vel_idx),fly_vel(vel_idx)*b(2) + b(1),'r')
    axis tight
end
linkaxes(get(gcf,'Children'))

%% plot the dff and heading trace of each trial
figure(10); clf; set(gcf,'color','none')
%tmp_labels = {'0.5','(0.5)','0.5 -> 1','(1)','1 -> 0.5','(0.5)'};
for i = 1:n
    subplot(n,1,i)
    yyaxis left
    set(gca,'xcolor','w','ycolor','w','color','none')

    if eb_flag
        alpha = linspace(-pi,pi,size(dff_cluster{i},1));
        ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'})
        hold on
    else
        alpha = linspace(-pi,3*pi,size(dff_cluster{i},1));
        ylim([-pi,3*pi]); yticks([-pi,pi,3*pi]); yticklabels({'0','2\pi','4\pi'})
        hold on
    end

    
    
    imagesc(xf{i},alpha,dff_cluster{i})
    hold on
    a = plot(xf{i},mu{i},'linewidth',2,'Color',[1,1,1]);
    a.YData(abs(diff(a.YData)) > pi) = nan;
    a.YData(rho{i}<rho_thresh) = nan;
    tmp_dist = circ_dist(mu{i},-cue{i});
    tmp_offset = circ_mean(tmp_dist(~isnan(tmp_dist)));
    tmp_cue = mod(unwrap(-cue{i}+pi) + tmp_offset,2*pi) - pi;
    
    b = plot(xf{i},tmp_cue,'linewidth',2,'Color',[.75,0,.75]);
    b.YData(abs(diff(b.YData)) > pi) = nan;
    if pb_flag
        b = plot(xf{i},-cue{i}+2*pi,'-','linewidth',2,'Color',[.75,0,.75]);
        b.YData(abs(diff(b.YData)) > pi) = nan;
    end
    
    ylabel('\DeltaF/F','Rotation',0)
    axis tight
    colormap(bone)
    yyaxis right
    set(gca,'ycolor','[.75,0,.75],','color','none')
    ylabel(sprintf('heading\n%s',tmp_labels{i}),'Rotation',0)
    yticks([])
end
xlabel('time (s)')
fontsize(gcf,20,'pixels')

%% evaluate DA peak as gain reporter
lag = 450; %in ms
lag = round((lag/1000)/fr);
%tmp_labels = {'Closed Loop','Dark','Closed Loop','Dark'};
figure(20); clf
for i = 1:n
    r_speed_lag     = r_speed{i}(1:end-lag);
    dff_peak_lag    = dff_peak{i}(lag+1:end);
    mu_lag          = mu{i}(lag+1:end);
    rho_lag         = rho{i}(lag+1:end);


    subplot(n,1,i)
    hold on
    yyaxis left
    plot(xf{i}(1:end-lag),dff_peak_lag,'Linewidth',2,'Color',[1,.5,0])
    set(gca,'ycolor',[1,.5,0])
    ylabel('peak dF/F')
    yyaxis right
    plot(xf{i}(1:end-lag),r_speed_lag,'LineWidth',2,'Color',[0.1,0.7,1])
    ylabel('r speed')
    set(gca,'color','none','ycolor',[0.1,.7,1],'xcolor','w')
    text(max(xlim),max(ylim),tmp_labels{i},'HorizontalAlignment',   'right','VerticalAlignment','top','color','white')
end

fontsize(gcf,20,'pixels')
set(gcf,'color','none')


figure(21); clf
for i = 1:n
    r_speed_lag     = r_speed{i}(1:end-lag);
    dff_peak_lag    = dff_peak{i}(lag+1:end);
    mu_lag          = mu{i}(lag+1:end);
    rho_lag         = rho{i}(lag+1:end);

    subplot(3,2,i)
    idx = r_speed_lag > 1e-1 & rho_lag > prctile(rho_lag,10);
    tmp = [abs(diff(unwrap(mu_lag)));0]./fr;
    tmp(tmp > prctile(tmp,80)) = prctile(tmp,80);
    scatter(r_speed_lag(idx),dff_peak_lag(idx),10,tmp(idx),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    colormap(cool)
    ylabel('dFF peak')
    xlabel('r speed')
    title(tmp_labels{i},'color','w')
    pos = get(gca,'Position');
    a = colorbar('color','w');
    ylabel(a,'bump speed')
    %set(gca,'Position',pos)
    set(gca,'color','none','ycolor','w','xcolor','w')
    set(gcf,'color','none')
end
fontsize(gcf,20,'pixels')

c = cool(2);

figure(22); clf
for i = 1:n
    subplot(3,2,i); cla
    hold on
    scatter(nan,nan,10,c(1,:),'filled','MarkerEdgeColor','none')
    scatter(nan,nan,10,c(2,:),'filled','MarkerEdgeColor','none')
    legend('slow','fast','textcolor','w','color','none','autoupdate','off')

    r_speed_lag     = r_speed{i}(1:end-lag);
    dff_peak_lag    = dff_peak{i}(lag+1:end);
    mu_lag          = mu{i}(lag+1:end);
    rho_lag         = rho{i}(lag+1:end);

    tmp = [abs(diff(unwrap(mu_lag)));0]./fr;
    tot_idx = r_speed_lag > 1e-1 & rho_lag > prctile(rho_lag,10);
    slow_idx = tmp < prctile(tmp,40);
    idx = tot_idx & slow_idx;


    scatter(r_speed_lag(idx),dff_peak_lag(idx),10,c(1,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    idx = tot_idx & ~slow_idx;
    scatter(r_speed_lag(idx),dff_peak_lag(idx),10,c(2,:),'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    ylabel('dFF peak')
    xlabel('r speed')
    title(tmp_labels{i},'color','w')
    %set(gca,'Position',pos)
    set(gca,'color','none','ycolor','w','xcolor','w')
    set(gcf,'color','none')
end

fontsize(gcf,20,'pixels')

%% fit model for dFF from rspeed and bump speed
lag = 0; %in ms
lag = round((lag/1000)/fr);

figure(24); clf; set(gcf,'color','none')
r_mat = nan(n,3,40);

for lag = 1:40
for i = 1:n
    r_speed_lag     = r_speed{i}(1:end-lag);
    dff_peak_lag    = dff_peak{i}(lag+1:end);
    mu_lag          = mu{i}(lag+1:end);
    rho_lag         = rho{i}(lag+1:end);
    mu_speed_lag    = [abs(diff(unwrap(mu_lag)));0]./fr;
    
    tot_idx = r_speed_lag > 0 & rho_lag > prctile(rho_lag,20) & ~isnan(mu_speed_lag) & mu_speed_lag < 6;
    
    [f,gof] = fit(r_speed_lag(tot_idx),dff_peak_lag(tot_idx),'poly1');
    r_mat(i,1,lag) = gof.adjrsquare;

    subplot(n,2,(2*i) - 1)
    scatter(r_speed_lag(tot_idx),dff_peak_lag(tot_idx),10,[0.75,0,0.75],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    set(gca,'color','none','xcolor','w','ycolor','w'); xlabel('r speed'); ylabel('dff peak')

    [f,gof] = fit([r_speed_lag(tot_idx),mu_speed_lag(tot_idx)],dff_peak_lag(tot_idx),'poly11');
    r_mat(i,2,lag) = gof.adjrsquare;
    subplot(n,2,(2*i))
    scatter3(r_speed_lag(tot_idx),mu_speed_lag(tot_idx),dff_peak_lag(tot_idx),10,[0.75,0,0.75],'filled','MarkerEdgeColor','none','MarkerFaceAlpha',.3)
    set(gca,'color','none','xcolor','w','ycolor','w','zcolor','w'); xlabel('r speed'); ylabel('bump speed'); zlabel('dff_peak')

    r_mat(i,3,lag) = corr(r_speed_lag(tot_idx),mu_speed_lag(tot_idx));
end
lag
end
%% fit the gain in a non-instantaneous way
% reps = 1;
% win_size = 300; %in seconds for piecewise optimization
% win_res = 0.5; %in seconds
% params_cell = {n,1};
% fval_cell = {n,1};
% 
% figure(11); clf
% a = scatter([1,2,3],[0,0,0],'filled');
% ylim([0,1])
% xticks([1,2,3]); xticklabels({'pieces','reps','trials'})
% ylabel('% done')
% title('Optimization Process')
% for i = 1:n
% params = nan(floor(max(xf{i})),2,reps);
% fval   = nan(floor(max(xf{i})),1,reps);
% for j = 1:reps
% fr      = mean(diff(xf{i})); %find the frame rate of the data
% 
% for t = 1:(max(xf{i})-win_size-20*fr)
% idx = xf{i} > t & xf{i} < t+win_size;
% 
% theta_tmp = cumsum(r_vel{i}(idx)*fr);
% 
% % xq = [t:win_res:t+win_size]';
% % mu_tmp = interp1(xf{i}(idx),mu{i}(idx),xq,[],'extrap');
% % theta_tmp = interp1(xf{i}(idx),theta_tmp,xq,[],'extrap');
% % rho_tmp = interp1(xf{i}(idx),rho{i}(idx),xq,[],'extrap');
% mu_tmp = mu{i}(idx);
% rho_tmp = rho{i}(idx);
% 
% obj_fun = @(x)(loss_fun(x,theta_tmp,mu_tmp,rho_tmp));
% %x = [gain,bias,lag,offset]
% lb = [0, 0];
% ub = [2, 20];
% %x0 = [.7, 10];
% x0 = rand(1,2).*(ub-lb) + lb;
% opts = optimoptions('fmincon','FiniteDifferenceStepSize',[sqrt(eps),.2]);
% [params(t,:,j),fval(t,1,j)] = fmincon(obj_fun,x0,[],[],[],[],lb,ub,[],opts);
% a.YData(1) = t/max(xf{i});
% drawnow
% end
% params_cell{i} = params;
% fval_cell{i} = fval;
% a.YData(2) = j/reps;
% drawnow
% end
% a.YData(3) = i/n;
% drawnow
% end





%% plot the fit values
figure(12)
rows = ceil(sqrt(n));
cols = ceil(n/rows);
for i = 1:n
    subplot(rows,cols,i)
    colormap(copper)
    scatter(params(:,i,1),params(:,i,2)*fr*1000,20,fval(:,i),'filled')
    title(tmp_labels{i})
    ylabel('lag (ms)'); xlabel('gain');
end
fontsize(gcf,20,'pixels')
shg

%% show the best fits
figure(13); clf
for i = 1:n
[~,idx] = min(fval(:,i));
x = squeeze(params(idx,i,:));
lag = ceil(x(2));

subplot(n,1,i)
hold on
a = plot(xf{i}(1:end-lag),mu{i}(1:end-lag),'linewidth',2,'Color',0*[1,1,1]);
a.YData(abs(diff(a.YData)) > pi) = nan;
tmp = rho{i}(1:end-lag);
a.YData(tmp < prctile(rho{i},20)) = nan;



theta_hat = theta_calc(x,r_vel{i}*fr);
b = plot(xf{i}(1:end-lag),theta_hat(lag+1:end),'linewidth',2,'Color',[.75,0,.75]);
b.YData(abs(diff(b.YData)) > pi) = nan;

axis tight
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})

ylabel(tmp_labels{i},'Rotation',0)
text(max(xlim),min(ylim),sprintf('g: %.2f\nlag: %.2f (ms)',x.*[1;fr*1000]),'HorizontalAlignment','left','VerticalAlignment','bottom')
end

xlabel('time (s)')
subplot(n,1,1)
pos =get(gca,'Position');
legend('$(\int{(\dot{\theta})*g})$','$\hat{\theta}$','Interpreter','latex','Location','northoutside')
set(gca,'Position',pos)
fontsize(gcf,20,'pixels')





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

function [amp_tot,amp_peak, mu, rho,dff_cluster] = bump_calc_eb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
[tmp_y,tmp_x] = find(bwmorph(mask,'shrink','inf'));
if length(tmp_y) == 1
    [y_mask,x_mask] = find(mask);
    figure(1); imagesc(mask);
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

%find the mean activity in each cluster
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
background  = imgData(~reshape(mask,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
%flash_idx   = med_pix > 0.05;
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

n_planes = size(regProduct,3);
filt_mu = linspace(max(centroids(:,1)),min(centroids(:,1)),n_planes);
filt_sig= range(centroids(:,1));
filt_x  = [1:size(regProduct,1)]';
f_cluster = zeros(n_centroid,size(regProduct,4));

for p = 1:n_planes
tmp = double(squeeze(regProduct(:,:,p,:))) .* normpdf(filt_x,filt_mu(p),filt_sig);


imgData_2d      = reshape(double(tmp),[],size(regProduct,4));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
    centroid_log(i, sub2ind(size(regProduct,1,2,4),y_mask(idx==i),x_mask(idx ==i))) = true;
end

tmp       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
f_cluster = f_cluster + tmp;
end

f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);            %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.

zscore_cluster  = zscore(dff_cluster(:,~flash_idx),[],2);
tmp = dff_cluster;
tmp(:,~flash_idx) = zscore_cluster;
zscore_cluster = tmp;
dff_cluster = zscore_cluster;

alpha       = linspace(-pi,pi,n_centroid);

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
xb  = linspace(0,total_t,size(regProduct,4))';

mu          = interp1(xb,unwrap(mu),xf)';
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp_tot     = interp1(xb,amp_tot,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';
end

function [amp_tot,amp_peak, mu, rho,dff_cluster] = bump_calc_pb(mask, regProduct, n_centroid, f0_pct, n_smooth, b_smooth, xf)
imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = reshape(imgData,[],size(imgData,3)); %reshape imgData to be allpixels x frames
tmp = bwmorph(mask,'thicken',20);
background  = imgData(~reshape(tmp,[],1),:); %extract points that are outside of the mask)
med_pix     = sum(background,1);
med_pix     = (med_pix - prctile(med_pix,10)) / prctile(med_pix,10);
flash_idx   = med_pix > median(med_pix) + 2*std(med_pix);
flash_idx   = logical(smoothdata(flash_idx,'gaussian',5));

imgData     = squeeze(sum(regProduct,3));   %regProduct is X x Y x Plane x Frames matrix. sum fluorescence across all planes, and squeeze into an X x Y x frames matrix
imgData     = imgData - prctile(background,2,'all');
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

f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
%f_cluster       = f_cluster  - median(bacground,1);
f_cluster(:,flash_idx) = nan;
f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
dff_cluster     = (f_cluster - f0) ./ f0;

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
xb  = linspace(0,total_t,size(regProduct,4))';

mu          = interp1(xb,unwrap(mu),xf)';
mu          = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
mu(mu > pi) = mu(mu > pi) - 2*pi;
rho         = interp1(xb,rho,xf)';
amp_tot     = interp1(xb,amp_tot,xf)';
amp_peak    = interp1(xb,amp_peak,xf)';
end

function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function loss = loss_fun(x,theta,mu,rho)
    gain = x(1);
    %bias = x(2) / 1000;
    lag = ceil(x(2));
    %offset = x(4);
    theta_hat = theta*gain; %subtract off walking bias from the observed rotational velocity, and then multiply this by the gain and integrate
    %theta_hat = theta_hat+offset; %subtract off the arbitrary offset
    theta_hat = mod(theta_hat,2*pi); %wrap this around -pi to pi
    theta_hat(theta_hat > pi) = theta_hat(theta_hat > pi) - 2*pi;

    tmp1 = theta_hat(1:end-lag); %take the accumulated rotation at time t
    tmp2 = mu(lag+1:end); %compare to the estimated heading at time t+lag
    rho = rho(lag+1:end); %apply the same shift to the vector strengths
    idx = rho > prctile(rho,10); %create an index to only take strong estimates
    tmp1(~idx) = nan;
    tmp2(~idx) = nan;    

    tmp  = circ_dist(tmp1,tmp2);
    loss = circ_var(tmp(~isnan(tmp))); %find the circular distance between the estimated heading (shifted backwards) and the integrated rotational velocity (shifted forwards)
end

function theta_hat = theta_calc(x,r_vel)
    gain = x(1);
    %bias = x(2) / 1000;
    lag = ceil(x(2));
    %offset = x(4);
    theta_hat = cumsum((r_vel)*gain); %subtract off walking bias from the observed rotational velocity, and then multiply this by the gain and integrate
    %theta_hat = theta_hat+offset; %subtract off the arbitrary offset
    theta_hat = mod(theta_hat,2*pi); %wrap this around -pi to pi
    theta_hat(theta_hat > pi) = theta_hat(theta_hat > pi) - 2*pi;
end