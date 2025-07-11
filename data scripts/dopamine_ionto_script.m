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

%all_data = struct();

tic
for i = length(all_data):length(all_files)
    clear img regProduct 

    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    load([tmp2.folder,'\',tmp2.name])

    all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
    all_data(i).im = process_im(squeeze(sum(img{1},3)), im_win, im_type, mask, n_centroid, f0_pct);
    if length(img) > 1; all_data(i).atp = process_im(squeeze(sum(img{2},3)), im_win, im_type, mask, n_centroid, f0_pct); end
    all_data(i).meta = all_files(i).folder;

    % if ~ismember('xb',fieldnames(all_data(i).ft))
    %     xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),size(all_data(i).im.d,2));
    % end

    try
     all_data(i).ft.stims = ftData_DAQ.stim{1};
    end
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%%
i  = 147;

tmp = zeros([size(all_data(i).im.d),3]);
tmp(:,:,1) = all_data(i).atp.d*5;
tmp(:,:,2) = all_data(i).im.d*5;

figure(1); clf
subplot(4,1,1); hold on
plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1)/max(sum(all_data(i).atp.d,1)),'r')
plot(all_data(i).ft.xb,sum(all_data(i).im.d,1)/max(sum(all_data(i).im.d,1)),'g')
plot(xlim,mean(sum(all_data(i).im.d,1)/max(sum(all_data(i).im.d,1)))*[1,1],':k')
title(all_data(i).meta)

subplot(4,1,2);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp)

subplot(4,1,3); hold on
%imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.d)
plot(all_data(i).ft.xb,unwrap(all_data(i).im.mu))
plot(all_data(i).ft.xf,-unwrap(all_data(i).ft.cue))

subplot(4,1,4)
plot(all_data(i).ft.xf,all_data(i).ft.stims)

linkaxes(get(gcf,'Children'),'x')
axis tight

%%
rows = ceil(sqrt(length(all_data)));
cols = ceil(length(all_data)/rows);
figure(3); clf
for i = 1:length(all_data)

    tmp = zeros([size(all_data(i).im.d),3]);
    tmp(:,:,1) = all_data(i).atp.d/2;
    tmp(:,:,2) = all_data(i).im.d*2;

    a1 = subplot(rows,cols,i);
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp)
    pos = get(gca,'Position');
    pos(2) = pos(2)+pos(4);
    pos(4) = .03;
    a2 = axes('Position',pos,'color','none');
    hold on
    plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1)/max(sum(all_data(i).atp.d,1)),'r')
    plot(all_data(i).ft.xb,sum(all_data(i).im.d,1)/max(sum(all_data(i).im.d,1)),'g')
    xticks([])

    linkaxes([a1,a2],'x')
    axis tight
end

pos=  get(gca,'Position');
legend('dopamine','gcamp','location','NortheastOutside')
set(gca,'Position',pos)

%% create an averaged trace
win = [-10,60];
dt = .1;
t = win(1):dt:win(2);

num_pulse = 0;
for i = 1:length(all_data)
    num_pulse = num_pulse+sum(diff(all_data(i).ft.stims)>0);
end

pulses = nan(num_pulse,length(t));
mus = nan(num_pulse,length(t));
cues = nan(num_pulse,length(t));
pulse_length = nan(num_pulse,1);
pulse_left = nan(num_pulse,1);

counter = 0;
for i = 1:length(all_data)
    
    
    stim_times = all_data(i).ft.xf(find(diff(all_data(i).ft.stims)>0)); %find the times at which pulse started
    stim_lengths = all_data(i).ft.xf(find(floor(diff(all_data(i).ft.stims))<-1)) - all_data(i).ft.xf(find(diff(all_data(i).ft.stims)>0));
    stim_left = sum(all_data(i).atp.d(1:end/2,:),[1,2]) < sum(all_data(i).atp.d(end/2:end,:),[1,2]);

    tmp_t   = all_data(i).ft.xb(1):dt:all_data(i).ft.xb(end);
    tmp_trace = interp1(all_data(i).ft.xb,sum(all_data(i).im.d,1),tmp_t);
    tmp_mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),tmp_t);
    tmp_cue = interp1(all_data(i).ft.xf,unwrap(all_data(i).ft.cue),tmp_t);

    for j = 1:length(stim_times)
        counter = counter+1;
        idx = tmp_t > stim_times(j)+win(1) & tmp_t < stim_times(j)+win(2); %extract stim times
        try idx(tmp_t > stim_times(j+1)) = false; end
        pulses(counter,1:sum(idx)) = tmp_trace(idx);
        tmp = tmp_mu(idx);
        mus(counter,1:sum(idx)) = tmp - tmp(round(abs(win(1))/dt));
        tmp = tmp_cue(idx);
        cues(counter,1:sum(idx)) = tmp - tmp(round(abs(win(1))/dt));
        pulse_length(counter) = stim_lengths(j);
        pulse_left(counter) = stim_left;
    end
end

idx = all(isnan(pulses),1);
pulses(:,idx) = [];
mus(:,idx) = [];
cues(:,idx) = [];
t(idx) = [];

figure(4); clf
subplot(3,1,1:2)
imagesc(t,1:num_pulse,pulses)
ylabel('individual pulses')

subplot(3,1,3); hold on
m = mean(pulses(pulse_length>1,:),1,'omitnan');
s = std(pulses(pulse_length>1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)] - m(1),[1,0,0])

m = mean(pulses(pulse_length<1,:),1,'omitnan');
s = std(pulses(pulse_length<1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)] - m(1),[1,.8,.8])

ylabel('mean + s.e.m.')
xlabel('time from stim start (s)')
legend('2s','.5s','autoupdate','off')
linkaxes(get(gcf,'Children'),'x')
axis tight
plot(xlim,[0,0],':k')

figure(5); clf
subplot(1,2,1); hold on; set(gca,'YDir','reverse')
m = mean(mus(pulse_length>1 & pulse_left==1,:),1,'omitnan');
s = std(mus(pulse_length>1 & pulse_left==1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],[1,.5,.5])

m = mean(mus(pulse_length>1 & pulse_left==0,:),1,'omitnan');
s = std(mus(pulse_length>1 & pulse_left==0,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],[.5,.5,1])


m = mean(cues(pulse_length<1 & pulse_left==1,:),1,'omitnan');
s = std(cues(pulse_length<1 & pulse_left==1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],'c')

m = mean(cues(pulse_length<1 & pulse_left==0,:),1,'omitnan');
s = std(cues(pulse_length<1 & pulse_left==0,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],'m')


plot(xlim,[0,0],':k')
title('2s eject')
ylabel('unwrapped bump position')
xlabel('time since stim start')

legend('left (mu)','right (mu)','left (cue)','right (cue)')

subplot(1,2,2); hold on; set(gca,'YDir','reverse')
m = mean(mus(pulse_length<1 & pulse_left==1,:),1,'omitnan');
s = std(mus(pulse_length<1 & pulse_left==1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],[1,.5,.5])

m = mean(mus(pulse_length<1 & pulse_left==0,:),1,'omitnan');
s = std(mus(pulse_length<1 & pulse_left==0,:),1,'omitnan') ./ sqrt(sum(~isnan(pulses),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],[.5,.5,1])


plot(xlim,[0,0],':k')
title('.5s eject')
fontsize(gcf,20,'pixels')

%% find the effect on fluorescence
dff_aligned = nan(32,900,length(all_data));
stim_length = nan(length(all_data),1);

figure(4); clf

for i = 1:length(all_data)
    start_idx = find(diff(all_data(i).ft.stims)>0,1);
    end_idx = find(diff(all_data(i).ft.stims)<0,1);
    stim_length(i) = all_data(i).ft.xf(end_idx) - all_data(i).ft.xf(start_idx);
end

for i = 1:length(all_data)
[~,ind] = max(sum(all_data(i).atp.d,2)); %find the peak stimulated glomerulus
tmp = circshift(all_data(i).im.d,size(all_data(i).im.d,1)/2 - ind,1);

%start_block = all_data(i).ft.xb < 60;
%end_block = all_data(i).ft.xb > all_data(i).ft.xb(end) - 60;
dff_aligned(:,1:299,i) = tmp(:,1:299);
dff_aligned(:,301:900,i) = tmp(:,end-599:end);
end

subplot(3,1,1); imagesc(mean(dff_aligned,3)); title({'mean aligned fluorescence','all trials'})
subplot(3,1,2); imagesc(mean(dff_aligned(:,:,stim_length<1),3)); title('0.5s'); ylabel('aligned glomeruli')
subplot(3,1,3); imagesc(mean(dff_aligned(:,:,stim_length>1),3)); title('2s'); xlabel({'frames','first ~30 seconds, last ~60 seconds'})


%% extract first trials for each fly
first_idx = false(length(all_data),1);
tmp_str = [];
for i = 1:length(all_data)
    if ~strcmp(tmp_str,all_data(i).meta(1:45))
        tmp_str = all_data(i).meta(1:45);
        first_idx(i) = true;
    end
end
    

%% show things for walking flies
i = 9;
figure(1); clf
a1 = subplot(2,1,1);
tmp = zeros([size(all_data(i).im.d),3]);
tmp(:,:,1) = all_data(i).atp.d*3;
tmp(:,:,2) = all_data(i).im.d;

imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp)
hold on
a = plot(all_data(i).ft.xb,all_data(i).im.mu); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue); a.YData(abs(diff(a.YData))>pi) = nan;

pos = get(gca,'Position');
pos(2) = pos(2)+pos(4);
pos(4) = .05;
a2 = axes('Position',pos,'color','none');
hold on
plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1)/max(sum(all_data(i).atp.d,1)),'r')
plot(all_data(i).ft.xb,sum(all_data(i).im.d,1)/max(sum(all_data(i).im.d,1)),'g')
xticks([])



a3 = subplot(2,1,2); hold on
a = plot(all_data(i).ft.xb,unwrap(all_data(i).im.mu) - all_data(i).im.mu(1)); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-unwrap(all_data(i).ft.cue) + all_data(i).ft.cue(1)); a.YData(abs(diff(a.YData))>pi) = nan;
a = scatter(all_data(i).ft.xf,zeros(length(all_data(i).ft.stims),1)); a.YData(~all_data(i).ft.stims) = nan;
set(gca,'YDir','reverse')
legend('mu','cue')

linkaxes([a1,a2,a3],'x')
axis tight
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