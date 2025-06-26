%% 
close all
clear all

%% load in data
base_dir = uigetdir(); %('Z:\pablo\gain_change\to do\'); %uigetdir(); %
all_files = dir([base_dir,'\**\*imagingData*.mat']);
all_files = natsortfiles(all_files);
%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])

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
im_win = {5,1};
n_centroid = 16;
f0_pct = 7;
r_thresh = .1;
rho_thresh = .1;

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
    % tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    % load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\csv\trialSettings.csv']);
    tmp2 = readtable([tmp2.folder,'\',tmp2.name]);

    imgData = squeeze(sum(regProduct,3));

    all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
    all_data(i).im = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    all_data(i).ft.pattern = tmp2.patternPath{1};

    try
     all_data(i).ft.stims = ftData_DAQ.stim{1};
    end
    
    
    dr = all_data(i).ft.r_speed;
    dr(abs(dr)>r_thresh) = dr(abs(dr)>r_thresh) - mean(dr(abs(dr)>r_thresh));
    all_data(i).ft.heading = cumsum(dr)/60;
    

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% calculate integrative gain
win_sec = 20;
lag = 20;
fr = 60;
win_frames = win_sec*fr;
r_thresh = .1;
rho_thresh = .1;
g = cell(length(all_data),1);
v = cell(length(all_data),1);

tic
for i = 1:length(all_data)
    fprintf('processing: %i ',i)

    %extract the best unwrapped estimate of mu by nan-ing low confidence
    %values and unwrapping and re-interp onto fictrac timescale
    m = all_data(i).im.mu;
    m(all_data(i).im.rho<rho_thresh) = nan;
    m = interp1(all_data(i).ft.xb,unwrap(m),all_data(i).ft.xf);
    m = smoothdata(m,1,"gaussian",60);
    amp = interp1(all_data(i).ft.xb,sum(all_data(i).im.d,1),all_data(i).ft.xf);

    %extract the fly's heading (no gain applied) and apply all lags
    h = all_data(i).ft.heading;

    m = m(lag+1:end);
    h = h(1:end-lag);


    % for each second
    g_tmp = nan(floor((length(m) - win_frames)/fr),1);
    v_tmp = nan(floor((length(m) - win_frames)/fr),1);
    
    for j = 1:length(g_tmp)
        f = j*fr;
        h_tmp = h(f:f+win_frames);
        if circ_var(h_tmp) > .1
        m_tmp = m(f:f+win_frames) - m(f);
        fun = @(x)(circ_var(circ_dist(m_tmp,h_tmp*x),[], [], [],'omitnan')); %find the gain and bias that best fits bump position to fly position over a window
        tmp = inf;
        for k = 0:.05:5 %evaluate the loss function at many possible gains, and save gain that gives the minimum value
            if fun(k) < tmp
                g_tmp(j) = k;
                tmp = fun(k);
            end
        end
        v_tmp(j) = tmp;
        end
    end

    g{i} = g_tmp;
    v{i} = v_tmp;
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_data)-i) / 60 / 60)
end

%% create figure to show example
i = 2;
binedges = 0:.05:5;

figure(1); clf
a1 = subplot(3,1,1);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)
xlabel('time (s)')

a2 = subplot(3,1,2);
plot(g{i}); hold on; plot(xlim,[.8,.8],':k')
ylabel('integrative gain')

linkaxes([a1,a2],'x')

subplot(3,2,5); hold on
h = histogram(g{i},'BinEdges',binedges,'FaceAlpha',.1);
plot(binedges(1:end-1),h.Values,'k','Linewidth',2)
xlabel('integrative gain')
ylabel('counts')

subplot(3,2,6); 
scatter(g{i},v{i},'filled','MarkerFaceAlpha',.5)
xlabel('integrative gain')
ylabel('func value')

%% show each
tmp_str = '20250611\fly 1';
tmp_ind = find(cellfun(@(x)(contains(x,tmp_str)),{all_data.meta}'));

rows = length(tmp_ind);

figure(1); clf
for ii = 1:rows
    a1= subplot(rows,1,ii);
    i = tmp_ind(ii);

    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
    hold on
    if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
    a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c); a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
   % a = plot(all_data(i).ft.xf,mod(all_data(i).ft.heading,2*pi) - pi,'g'); a.YData(abs(diff(a.YData))>pi) = nan;

    % pos = get(gca,'Position');
    % pos(2) = pos(2)+pos(4);
    % pos(4) = .05;
    % a2 = axes('Position',pos,'color','none');
    % hold on
    % scatter(1:length(all_data(i).gain.g),all_data(i).gain.g,'.')
    %ylim([-1,5])
    %linkaxes([a1,a2],'x')
end

subplot(rows,1,1)
title(all_data(tmp_ind(1)).meta)

%%
win_sec = 10;
lag = 10;
fr = 60;
win_frames = win_sec*fr;
r_thresh = .1;
rho_thresh = .1;

figure(2); clf
for ii = 1:rows
    
    i = tmp_ind(ii);
    dm = gradient(interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf)) * fr;
    dr = all_data(i).ft.r_speed;
    rho = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);
    
    dm(abs(dr) < r_thresh & rho < rho_thresh) = nan;

    dm = dm(lag+1:end);
    dr = dr(1:end-lag);

    g = nan(size(dm));

    for j = 1:length(dm)-win_frames
        dr_tmp = dr(j:j+win_frames);
        dm_tmp = dm(j:j+win_frames);
        idx = ~isnan(dm_tmp);
        g(j) = dr_tmp(idx) \ dm_tmp(idx);
    end
        
        
     subplot(rows,1,ii)
     plot(g)
     hold on
     plot(xlim,[.8,.8],':k')
     title(mean(g,'omitnan'))
end

figure(3); clf
for ii = 1:rows
    subplot(rows,1,ii);
    i = tmp_ind(ii);

    hold on
    plot(all_data(i).ft.xf,-unwrap(all_data(i).ft.cue),'m');
    plot(all_data(i).ft.xb,unwrap(all_data(i).im.mu),'k');
    axis tight
    plot(xlim,repmat(min(ylim):2*pi:max(ylim),2,1)','Color',.8*[1,1,1])
    %set(gca,'YGrid',')
end

figure(4); clf
for ii = 1:rows
    
    i = tmp_ind(ii);
    m = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
    dr = all_data(i).ft.r_speed;
    rho = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);
    
    m(abs(dr) < r_thresh & rho < rho_thresh) = nan;

    m = m(lag+1:end);
    dr = dr(1:end-lag);

    g = nan(size(m));

    for j = 1:length(dm)-win_frames
        dr_tmp = dr(j:j+win_frames);
        m_tmp = m(j:j+win_frames);
        idx = ~isnan(m_tmp);
        fun = @(x)(var(m_tmp - cumsum(x(1)*dr_tmp)/fr,'omitnan')); %find the gain and bias that best fits bump position to fly position over a window
    

        g(j) = fminsearch(fun,1);
    end
        
        
     subplot(rows,1,ii)
     plot(g)
     hold on
     plot(xlim,[.8,.8],':k')
     title(mean(g,'omitnan'))
end

%% plot peak fluorescence as a function of fly speed and bump speed
i = 16;
lag = 10;
figure(3); clf
dm = gradient(interp1(all_data(i).ft.xb,smoothdata(unwrap(all_data(i).im.mu),3,'movmean',3),all_data(i).ft.xf)) * fr;
dr = all_data(i).ft.r_speed;
dc = -gradient(unwrap(all_data(i).ft.cue)) * fr;
amp= interp1(all_data(i).ft.xb,sum(all_data(i).im.d,1),all_data(i).ft.xf);
rho= interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);

dm = dm(lag+1:end);
dr = dr(1:end-lag);
dc = dc(1:end-lag);
amp = amp(lag+1:end);
rho = rho(lag+1:end);

idx = rho > rho_thresh & abs(dr) > r_thresh;
subplot(2,2,1)
scatter(abs(dr(idx)),abs(dm(idx)),[],amp(idx),'filled','MarkerFaceAlpha',.1)
xlabel('fly speed')
ylabel('bump speed')
axis equal

subplot(2,2,2)
scatter(abs(dr(idx)),abs(dc(idx)),[],amp(idx),'filled','MarkerFaceAlpha',.1)
xlabel('fly speed')
ylabel('cue speed')
axis equal

t = 1:round(max(all_data(i).ft.xf));
g = nan(size(t));
for j = 1:length(t)
    idx = all_data(i).ft.xf > j-1 & all_data(i).ft.xf < j;
    
    g(j) = all_data(i).ft.r_speed(idx) \ -gradient(all_data(i).ft.cue(idx))*60;
end

subplot(2,1,2); hold on
plot(all_data(i).ft.xb,max(all_data(i).im.d,[],1))
scatter(t,g,'.')
ylim([0,1.5])
yyaxis right; plot(all_data(i).ft.xf,all_data(i).ft.r_speed)

%% cross-correlate max fluorescence with compass gain
i = 18;
rho_thresh = .2;
dt = 1/60;
t = 0:dt:max(all_data(i).ft.xb);
fr_f = mean(diff(all_data(i).ft.xf));
fr_b = mean(diff(all_data(i).ft.xb));
figure(3); clf; tiledlayout(4,2)

m = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),t,'linear','extrap')';
%m = smoothdata(m,1,'movmean',.5/dt);
dm = gradient(m); % find the displacement of pva between imaging frames 

r = interp1(all_data(i).ft.xf,cumsum(all_data(i).ft.r_speed*fr_f),t,'linear','extrap')';
%r = smoothdata(r,1,'movmean',.5/dt);
dr = gradient(r); % find the displacement of the fly heading between imaging frames

rho = interp1(all_data(i).ft.xb,all_data(i).im.rho,t,'linear','extrap');
m(rho<rho_thresh) = nan;
dm(rho<rho_thresh) = 0; %if we have a poor estimate of the bump, then set the displacement to 0.

amp = interp1(all_data(i).ft.xb,max(all_data(i).im.z,[],1),t,'linear','extrap')';

%plot the image with heading data
a1= nexttile([1,2]);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)

%plot the calculated bump and fly position
a2 = nexttile([1,2]); hold on
scatter(t,m,'.')
scatter(t,r,'.')
xlabel('time (s)'); ylabel('unwrapped position (radians)')
legend('bump','fly')

linkaxes([a1,a2],'x')
axis tight

%cross correlate the bump and fly displacement
nexttile; hold on
[r,lags] = xcorr(dm,dr,ceil(5/dt),'normalized');
stem(lags*dt,r)
plot([0,0],ylim,'k')
xlabel('lag (s)'); ylabel('corr bump disp x fly disp')

[~,idx] = max(r); %find the optimal lag
lag = lags(idx); 
dr = dr(1:end-lag);
dm = dm(1+lag:end);
dg = abs(dm);
amp = amp(1+lag:end);
idx = abs(dr) > .01 & abs(dm) > .01;

nexttile
scatter(dr(idx),dm(idx),[],amp(idx),'filled','MarkerFaceAlpha',.4)
xlabel('fly displacement'); ylabel('bump displacement')
text(max(xlim),max(ylim),sprintf('gain: %.2f',dr(idx) \ dm(idx) ),'HorizontalAlignment','right','VerticalAlignment','top')

% cross correlate bump amp with bump disp per fly disp, with appropriate
% lag applied
nexttile; hold on
[r,lags] = xcorr(dg,amp,ceil(5/dt),'normalized');
stem(lags*dt,r)
plot([0,0],ylim,'k')
xlabel('lag (s)'); ylabel('corr amp x bump disp')

nexttile
[~,idx] = max(r); %find the optimal lag
lag = lags(idx);
dr = dr(1:end+lag);
dg = dg(1:end+lag);
amp = amp(1-lag:end);
scatter(abs(dr),amp,[],dg,'filled','MarkerFaceAlpha',.4)
xlabel('abs(fly displacement)'); ylabel('max amp')
colormap(parula)
%% Functions

function s = process_ft(ftData_DAQ, ft_win, ft_type)

    %f_speed = ftData_dat.velFor{:};                       %store each speed
    %f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    f_speed = ftData_DAQ.velFor{:};
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