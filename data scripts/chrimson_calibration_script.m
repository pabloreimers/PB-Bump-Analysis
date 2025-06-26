%% load in data
base_dir = uigetdir(); %('Z:\pablo\chrimson_calibration\to do\'); %uigetdir(); %
all_files = dir([base_dir,'\**\*imagingData*.mat']);
all_files = natsortfiles(all_files);
%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        
        if ~exist('regProduct','var')
            regProduct = img{1};
        end

        imgData = squeeze(sum(regProduct,3));


        top_pct = prctile(imgData,98,'all');
        bot_pct = prctile(imgData,5,'all');
        
        imgData(imgData>top_pct) = top_pct;
        imgData(imgData<bot_pct) = bot_pct;        

        figure(1); clf; imagesc(mean(imgData,3)); colormap(bone); axis equal tight; drawnow;
        mask = roipoly();
        save([fileparts(all_files(i).folder),'\mask.mat'],'mask')
    end
end

%% process and store all values
ft_type= 'movmean'; %the type of smoothing for fictrac data
ft_win = 10; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'movmean','movmean'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {1,1};
n_centroid = 16;
f0_pct = 7;

all_data = struct();

tic
for i = 1:length(all_files)
    clear img regProduct 

    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    load([tmp2.folder,'\',tmp2.name])

    if ~exist('regProduct','var')
        regProduct = img{1};
    end

    imgData = squeeze(sum(regProduct,3));

    all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
    all_data(i).im = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    if ~ismember('xb',fieldnames(all_data(i).ft))
        xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),size(all_data(i).im.d,2));
    end

    try
     all_data(i).ft.stims = ftData_DAQ.stim{1};
    end
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% show
tmp_str = '20250221\fly 1';
trial_num = 2;



tmp_ind = find(cellfun(@(x)(contains(x,tmp_str)),{all_data.meta}'));
i = tmp_ind(trial_num);

idx = logical(interp1(all_data(i).ft.xf,double(~all_data(i).ft.stims),all_data(i).ft.xb,'linear','extrap'));
h0 = prctile(all_data(i).im.f_nmask(idx),5);

figure(1); clf
subplot(3,1,1)
plot(all_data(i).ft.xb,(all_data(i).im.f_mask - all_data(i).im.f_nmask)/h0); ylabel('dFF')
%plot(all_data(i).ft.xb,mean(all_data(i).im.d,1))
 ylabel('(f_{mask} - f_{nmask})/f0_{nmask}')
title(all_data(i).meta)
subplot(3,1,2)
plot(all_data(i).ft.xf,all_data(i).ft.stims/10); ylabel('stim light')
subplot(3,1,3)
plot(all_data(i).ft.xf,abs(all_data(i).ft.r_speed)); ylabel('rot speed (rad/s)')

linkaxes(get(gcf,'Children'),'x')

%% show baseline traces for each fly

figure(2); clf; 

subplot(3,1,1); hold on; ylabel('dff'); ylabel('$\frac{F_{mask} - F_{nonmask}}{F_{5,nonmask}}$','Interpreter','Latex','Rotation',0,'fontsize',20); axis tight
subplot(3,2,3); hold on; ylabel('$\frac{F_{5,mask}}{F_{5,nonmask}}$','Interpreter','Latex','Rotation',0,'fontsize',30); xticks([0,1]); xticklabels({'+','CsChrimson'}); title('baseline Ca++')
subplot(3,2,4); hold on; ylabel('$\frac{F_{95,mask}}{F_{5,nonmask}}$','Interpreter','Latex','Rotation',0,'fontsize',30); xticks([0,1]); xticklabels({'+','CsChrimson'}); title('peak Ca++ (non-stim)')
subplot(3,1,3); hold on; ylabel('$\frac{F_{5,mask}}{F_{5,nonmask}}$','Interpreter','Latex','Rotation',0,'fontsize',30)
for i = 1:3

    tmp_idx = contains(all_data(i).meta,'cschrimson');
    chill_idx = contains(all_data(i).meta,'chilled');
    idx = logical(interp1(all_data(i).ft.xf,double(~all_data(i).ft.stims),all_data(i).ft.xb,'linear','extrap'));

    h0 = prctile(all_data(i).im.f_nmask(idx),5);
    h1 = prctile(all_data(i).im.f_mask(idx),5);
    h2 = prctile(all_data(i).im.f_mask(idx),90);

    if  (h2 - h1) / h0 < .1 %| ~tmp_idx %| any(all_data(i).ft.stims)
        continue
    end

    subplot(3,1,1)
    plot(all_data(i).ft.xb,(all_data(i).im.f_mask - all_data(i).im.f_nmask)/h0,'Color',[tmp_idx,.5*chill_idx,1-tmp_idx])

    subplot(3,2,3)
    scatter3(tmp_idx+(rand(1)-.5)/10,(h1)/h0,i,100,'MarkerEdgeColor',[tmp_idx,.5*chill_idx,1-tmp_idx])

    subplot(3,2,4)
    scatter(tmp_idx+(rand(1)-.5)/10,(h2)/h0,100,'MarkerEdgeColor',[tmp_idx,.5*chill_idx,1-tmp_idx])

    subplot(3,1,3)
    scatter(i,(h1)/h0,100,'MarkerEdgeColor',[tmp_idx,.5*chill_idx,1-tmp_idx])
end

%% show effect of light stimulation on both groups

figure(3); clf;
subplot(2,1,1); hold on; axis tight; ylabel('dF/F (+)')
subplot(2,1,2); hold on; axis tight; ylabel('dF/F (CsChrimson)')
for i = 1:length(all_data)
    if any(all_data(i).ft.stims)
        tmp_idx = contains(all_data(i).meta,'cschrimson');
        subplot(2,1,1+tmp_idx)
        plot(all_data(i).ft.xb,mean(all_data(i).im.d,1),'Color',[tmp_idx,.5,1-tmp_idx])
    end
end

chrim_idx = contains({all_data.meta},'cschrimson');
figure(4); clf;  hold on
tmp = find(~chrim_idx);
for i = 1:length(tmp)
    if any(all_data(tmp(i)).ft.stims)
    plot(all_data(tmp(i)).ft.xb,mean(all_data(tmp(i)).im.d,1) + 2*i,'w')
    a = plot(all_data(tmp(i)).ft.xf,0*all_data(tmp(i)).ft.xf + 2*i,'Color',[1,.6,.6],'linewidth',1.5); a.YData(~all_data(tmp(i)).ft.stims)=nan;
    end
end
title('mean dff traces for chrimson flies')
xlabel('time (s)')
set(gca,'Color','none','ycolor','w','xcolor','w')
set(gcf,'Color','none','InvertHardCopy','off')

%% show average dff during stim and off stim for both groups
in_stim = nan(length(all_data),1);
out_stim= nan(length(all_data),1);
stim_idx = false(length(all_data),1);
chrim_idx = contains({all_data.meta},'cschrimson')';

for i = 1:length(all_data)

    idx = logical(interp1(all_data(i).ft.xf,double(~all_data(i).ft.stims),all_data(i).ft.xb,'linear','extrap'));

    f0 = prctile(all_data(i).im.f_mask(idx),5);
    

    stim_idx(i) = any(all_data(i).ft.stims);
    in_stim(i) = (mean(all_data(i).im.f_mask(~idx)) - f0) / f0;
    out_stim(i)= (mean(all_data(i).im.f_mask(idx)) - f0) / f0;
end

figure(1); clf; hold on
plot([1,2],[out_stim(chrim_idx&stim_idx),in_stim(chrim_idx&stim_idx)],'.-k')
plot([3,4],[out_stim(~chrim_idx&stim_idx),in_stim(~chrim_idx&stim_idx)],'.-k')


xlim([.5,4.5])
ylabel('average dff')
xticks([])
text(.5,0,'cschrimson','HorizontalAlignment','right','VerticalAlignment','top')
text(.5,-.1,'light','HorizontalAlignment','right','VerticalAlignment','top')
text([1,2,3,4],[0,0,0,0],{'+','+','-','-'},'HorizontalAlignment','right','VerticalAlignment','top')
text([1,2,3,4],[0,0,0,0]-.1,{'-','+','-','+'},'HorizontalAlignment','right','VerticalAlignment','top')

%% look at behavioral effects of chrimson vs control
r_thresh = .1;
f_thresh = .1;
f_max    = 20;

figure(4); clf;  hold on
tmp = find(chrim_idx & stim_idx);
for i = 1:length(tmp)
    a = plot(all_data(tmp(i)).ft.xf,abs(all_data(tmp(i)).ft.r_speed) + 4*i,'k'); a.YData(a.YData > 5+4*i) = nan;
    a = plot(all_data(tmp(i)).ft.xf,0*all_data(tmp(i)).ft.xf + 4*i,'Color',[1,.6,.6],'linewidth',1.5); a.YData(~all_data(tmp(i)).ft.stims)=nan;
end
title('mean f speed traces for chrimson flies')
xlabel('time (s)')


in_stim_r = nan(length(all_data),1);
out_stim_r= nan(length(all_data),1);
in_stim_f = nan(length(all_data),1);
out_stim_f= nan(length(all_data),1);
stim_idx = false(length(all_data),1);
chrim_idx = contains({all_data.meta},'cschrimson')';

for i = 1:length(all_data)    

    stim_idx(i) = any(all_data(i).ft.stims);
    mov_idx = (abs(all_data(i).ft.r_speed) > r_thresh | all_data(i).ft.f_speed > f_thresh) & abs(all_data(i).ft.f_speed) < f_max;
    in_stim_r(i) = mean(abs(all_data(i).ft.r_speed(logical(all_data(i).ft.stims) & mov_idx)));
    out_stim_r(i)= mean(abs(all_data(i).ft.r_speed(~logical(all_data(i).ft.stims) & mov_idx)));
    in_stim_f(i) = mean(all_data(i).ft.f_speed(logical(all_data(i).ft.stims) & mov_idx));
    out_stim_f(i)= mean(all_data(i).ft.f_speed(~logical(all_data(i).ft.stims) & mov_idx));

end

figure(5); clf; subplot(2,1,1); hold on
plot([1,2],[out_stim_r(chrim_idx&stim_idx),in_stim_r(chrim_idx&stim_idx)],'.-k')
plot([3,4],[out_stim_r(~chrim_idx&stim_idx),in_stim_r(~chrim_idx&stim_idx)],'.-k')
ylabel('average r speed')
xlim([.5,4.5])
xticks([])

subplot(2,1,2); hold on
plot([1,2],[out_stim_f(chrim_idx&stim_idx),in_stim_f(chrim_idx&stim_idx)],'.-k')
plot([3,4],[out_stim_f(~chrim_idx&stim_idx),in_stim_f(~chrim_idx&stim_idx)],'.-k')

ylim([-1,max(ylim)])
b = min(ylim);
xlim([.5,4.5])
ylabel('average f speed')
xticks([])
text(.5,b,'cschrimson','HorizontalAlignment','right','VerticalAlignment','top')
text(.5,b-range(ylim)/10,'light','HorizontalAlignment','right','VerticalAlignment','top')
text([1,2,3,4],b+[0,0,0,0],{'+','+','-','-'},'HorizontalAlignment','right','VerticalAlignment','top')
text([1,2,3,4],b+[0,0,0,0]-range(ylim)/10,{'-','+','-','+'},'HorizontalAlignment','right','VerticalAlignment','top')


%% plot bump things for an example fly
i = 3;
lag = 20;

figure(6); clf
subplot(2,1,1); 
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.d)
hold on
a = plot(all_data(i).ft.xf,max(ylim)*ones(length(all_data(i).ft.xf),1),'r','linewidth',2);
a.YData(~all_data(i).ft.stims) = nan;
title(all_data(i).meta)

subplot(2,1,2); hold on
plot(all_data(i).ft.xf,cumsum(all_data(i).ft.r_speed/60),'m')
plot(all_data(i).ft.xb,unwrap(all_data(i).im.mu),'k')
a = plot(all_data(i).ft.xf,max(ylim)*ones(length(all_data(i).ft.xf),1),'r','linewidth',2);
a.YData(~all_data(i).ft.stims) = nan;
legend('fly','bump')
ylabel('unwrapped position')
set(gca,'YDir','reverse')

linkaxes(get(gcf,'Children'),'x')
axis tight
% look at correlation between bump and cue 



figure(7); clf

m_speed = gradient(interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf,'linear','extrap')) * 60;
r_speed = all_data(i).ft.r_speed;
m_speed = m_speed(1+lag:end);
r_speed = r_speed(1:end-lag);
idx = abs(r_speed) > r_thresh;
stims = all_data(i).ft.stims(1:end-lag);

subplot(2,1,1); hold on
plot(all_data(i).ft.xf(1:end-lag),r_speed,'m')
plot(all_data(i).ft.xf(1:end-lag),m_speed,'k')
a = plot(all_data(i).ft.xf,min(ylim)*ones(length(all_data(i).ft.xf),1),'r','linewidth',2);
a.YData(~stims) = nan;
ylabel('rotational speed (rad/s)')
xlabel('time (s)')
legend('fly','bump')
title(all_data(i).meta)

subplot(2,2,3)
hold on
scatter(r_speed(idx),m_speed(idx),'.')
%scatter(r_speed(idx&stims),m_speed(idx&stims),'.r')
c = corr(r_speed(idx),m_speed(idx));
g = polyfit(r_speed(idx),m_speed(idx),1);
b = g(2);
g = g(1);
plot([-2,2],g(1)*[-2,2],'m','linewidth',2)
plot([-2,2],[-2,2],':k','linewidth',2)
text(max(xlim),0,sprintf('corr: %.2f\ngain: %.2f\nbias: %.2f\nlag: %dms',c,g,b,round(lag/60*1e3)),'horizontalalignment','left','VerticalAlignment','middle')
xlabel('fly speed (rad/s)')
ylabel('bump speed (rad/s)')
axis tight

%% calculate bump things for all flies

g_vec = nan(length(all_data),1);
b_vec = nan(length(all_data),1);
c_vec = nan(length(all_data),1);

chrim_idx = contains({all_data.meta},'cschrimson')';
chill_idx = contains({all_data.meta},'chill')';
stims_idx = false(size(chrim_idx));

for i = 1:length(all_data)
    stims_idx(i) = any(all_data(i).ft.stims);

    try
    m_speed = gradient(interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf,'linear','extrap')) * 60;
    r_speed = all_data(i).ft.r_speed;
    m_speed = m_speed(1+lag:end);
    r_speed = r_speed(1:end-lag);
    idx = abs(r_speed) > r_thresh;
    stims = all_data(i).ft.stims(1:end-lag);

    g = polyfit(r_speed(idx),m_speed(idx),1);
    g_vec(i) = g(1);
    b_vec(i) = g(2);
    [c_vec(i)] = corr(r_speed(idx),m_speed(idx));
    catch
        fprintf('fly %i failed\n',i)
    end
end

figure(8); clf
group_idx = stim_idx + 2*chrim_idx;
num_group = length(unique(group_idx));

hold on
scatter(group_idx,b_vec,'r')
scatter(group_idx(chill_idx),b_vec(chill_idx),'b')
plot(xlim,[0,0],':k')
legend('warm','chilled')
ylabel('bias')
xlim([-.5,3.5])
xticks([])
b = min(ylim);
text(-.5,b,'cschrimson','HorizontalAlignment','right','VerticalAlignment','top')
text(-.5,b-range(ylim)/10,'light','HorizontalAlignment','right','VerticalAlignment','top')
text(unique(group_idx),b+zeros(num_group,1),{'-','-','+','+'},'HorizontalAlignment','right','VerticalAlignment','top')
text(unique(group_idx),b+zeros(num_group,1)-range(ylim)/10,{'-','+','-','+'},'HorizontalAlignment','right','VerticalAlignment','top')

%% Functions

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


function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

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

    mask_2d = reshape(mask,[],1);
    f_mask = mean(imgData_2d(mask_2d,:),1);
    f_nmask= mean(imgData_2d(~mask_2d,:),1);

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.d  = dff_cluster;
    s.f  = f_cluster;
    s.alpha = alpha;
    s.f_mask = f_mask;
    s.f_nmask= f_nmask;

    %s.imgData = imgData;
end

function out = resample_pablo(N,x)
    idx = randi(length(x),length(x),N);
    out = x(idx);
end


function [cc,gain,bias] = find_cc(lag,fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh)
    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;

    cc  = corr(fly_vel(idx),bump_vel(idx));
    b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); 
    gain = b(2);
    bias = b(1);
end

function [gain,bias] = find_gain(fly_vel,mu,fr)

    int_pos = cumsum(fly_vel * fr);

    tmp1 = mu(~isnan(mu));
    tmp2 = int_pos(~isnan(mu));
    obj_fun = @(x) circ_var(circ_dist(tmp1,tmp2*x(1) + x(2)));
    x0 = [.7,0];
    lb = [-10,-10];
    ub = [10,10];
   
    x = fmincon(obj_fun,x0,[],[],[],[],lb,ub);
    gain = x(1);
    bias = x(2);
end

function h = plot_sem(ax,t,x)

m1 = mean(x,1,'omitnan');
s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

idx = ~isnan(m1);
m1 = m1(idx);
s1 = s1(idx);
t  = t(idx);


h = patch(ax,[t;flipud(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end

function rgb_image = mat2rgb(v,map)
    minv = min(v(:));
    maxv = max(v(:));
    ncol = size(map,1);
    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
    rgb_image = ind2rgb(s,map);
end