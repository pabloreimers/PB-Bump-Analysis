%% clear all
clear all
close all

%% find path to all relevant files
base_dir = ('Z:\pablo\stacks\lpsp_kir_redo\');
all_files = dir([base_dir,'\**\*denoised*.mat']);

%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        imgData = mean(regProduct,[3,4]);
        figure(1); clf; imagesc(imgData); colormap(bone); drawnow;
        tmp = drawellipse('Center',[size(imgData,2)/2,size(imgData,1)/2],'SemiAxes',[size(imgData,2)/4,size(imgData,1)/4]);
        input('') %move on once user presses enter, after adjusting the ellipse
        mask = createMask(tmp);
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

all_data = struct();

tic
for i = 1:length(all_files)
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{7})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])

    all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
    all_data(i).im = process_im_3d(regProduct, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% plot all
dark_mode = true;

vel_thresh = .2;
bump_thresh = 10;
rho_thresh = .2;
vel_max = 5;
lag = 7;

cc       = nan(length(all_data),1);
gains    = nan(length(all_data),1);
lpsp_idx = false(length(all_data),1);
dark_idx = false(length(all_data),1);
fly_id   = cell(length(all_data),1);


for i = 1:length(all_data)
    lpsp_idx(i) = contains(all_data(i).meta,'_lpsp_','IgnoreCase',true);
    
    tmp = strsplit(all_data(i).meta,'\');
    if lpsp_idx(i)
        title(strcat(tmp(5),tmp(6)),'Color','r')
    else
        title(strcat(tmp(5),tmp(6)),'Color','k')
    end

    tmp2 = regexp(tmp{7},'\d*','Match');
    dark_idx(i) = mod(str2double(tmp2{2}),2) == 0;
    fly_id{i} = [tmp{5},'_',tmp{6}];
end
   
group_order = {'empty (CL)','LPsP (CL)','empty (dark)','LPsP (dark)'};
ind = 2*dark_idx + lpsp_idx + 1;
t = cell(size(unique(ind)));
for i = unique(ind)'
figure(i); clf
cols = ceil(sqrt(sum(ind==i)));
rows = ceil(sum(ind==i)/cols);
t{i} = tiledlayout(rows,cols);
end
    
if dark_mode
    c = 'w';
else
    c = 'k';
end

for i = 1:length(all_data)
    figure(ind(i)); nexttile; hold on
    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = [0;-diff(all_data(i).ft.cue)]/fr; %all_data(i).ft.r_speed;
    bump_vel = [0;diff(interp1(xb,unwrap(all_data(i).im.mu),xf))]/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);

    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh & abs(fly_vel) < vel_max;
    scatter(fly_vel(idx),bump_vel(idx),5,c,'filled','markerfacealpha',.1)
    axis equal
    y = ylim; x = xlim;
    plot(x,[0,0],':','Color',c); 
    plot([0,0],y,':','Color',c);
    b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); %fit the slope of fly vel and bump vel with an arbitrary offset
    r = corr(fly_vel(idx),bump_vel(idx));
    plot([0,0],y,':','Color', c)
    plot(x,[0,0],':','Color', c)
    plot(x,x*b(2) + b(1),'r')
    %text(x(2),y(1),sprintf('gain: %.2f\nr: %.2f',b(2),r),'HorizontalAlignment','right','VerticalAlignment','bottom','color',c)
    text(x(2),y(1),sprintf('gain: %.2f',b(2)),'HorizontalAlignment','right','VerticalAlignment','bottom','color',c)
    xlim(x); ylim(y);

    cc(i) = r;
    gains(i) = b(2);

    set(gca,'xcolor',c,'ycolor',c,'color','none')
end

[~,~,fly_num] = unique(fly_id);

for i = unique(ind)'
    %fontsize(figure(i),20,'pixels')
    title(t{i},group_order(i),'fontsize',10,'color',c)
    %xlabel(t{i},'fly vel (rad/s)','fontsize',30,'color',c); ylabel(t{i},'bump vel (rad/s)','fontsize',30,'color',c)
    if dark_mode
        set(figure(i),'color','none')
    end
end


%% confirm that heading traces look reasonable
group_order = {'empty (CL)','LPsP (CL)','empty (dark)','LPsP (dark)'};
ind = 2*dark_idx + lpsp_idx + 1;
t = cell(size(unique(ind)));

for i = unique(ind)'
figure(i); clf
cols = 2;
rows = ceil(sum(ind == i)/cols);
t{i} = tiledlayout(rows,cols);
end

for i = 1:length(all_data)
    cue = all_data(i).ft.cue;
    xb = linspace(min(all_data(i).ft.xf),max(all_data(i).ft.xf),length(all_data(i).im.mu));
    
    figure(ind(i)); nexttile; hold on
    imagesc(xb,all_data(i).im.alpha,all_data(i).im.z); colormap(parula)
    a = plot(xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'m'); a.YData(abs(diff(a.YData))>pi) = nan;
    axis tight
    xticks([]); yticks([])
    text(max(xf),0,sprintf('gain: %.2f\nR: %.2f',gains(i),cc(i)),'VerticalAlignment','middle','HorizontalAlignment','left','color','w')
    ylabel(fly_num(i))
    set(gca,'color','none','ycolor','w','xcolor','w')
end

for i = unique(ind)'
    figure(i)
    xticks('auto'); xlabel('time (s)')
    fontsize(gcf,20,'pixels')
    title(t{i},group_order(i),'fontsize',40,'color','w')
    set(gcf,'color','none')
end

%% compare correlation coefficients and gains
dark_mode = true;

cmap = [0,1,1;
        1,0,0];
group_order = {'empty (CL)','LPsP (CL)','empty (dark)','LPsP (dark)'};
ind = 2*dark_idx + lpsp_idx + 1;
figure(5); clf

if dark_mode
    c = 'w';
else
    c = 'k';
end

subplot(2,2,1)
scatter(ind, cc,[],cmap(lpsp_idx+1,:),'filled','MarkerFaceAlpha',.5); hold on
m = accumarray(ind,cc,[],@mean);
s = accumarray(ind,cc,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'o','Color',c,'linewidth',2)
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Correlation Coefficent', '(fly vs bump vel)'})
set(gca,'color','none','ycolor',c,'xcolor',c)
ylim([-.1,1.2])
subplot(2,2,2)
scatter(ind, gains,[],cmap(lpsp_idx+1,:),'filled','MarkerFaceAlpha',.5); hold on
m = accumarray(ind,gains,[],@mean);
s = accumarray(ind,gains,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'o','Color',c,'linewidth',2)
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Gain', '(fly vs bump vel)'})
set(gca,'color','none','ycolor',c,'xcolor',c)
ylim([0,1.5])
if dark_mode; set(gcf,'color','none','InvertHardcopy','off'); end
%% test significance (bootstrap to ask about a mean difference, we don't know how variances compare)
N = 1e4;

%do for corr coeff
lpsp_cl  = mean(resample_pablo(N,cc(~dark_idx & lpsp_idx)),1);
empty_cl = mean(resample_pablo(N,cc(~dark_idx & ~lpsp_idx)),1);
p_cl     = sum( lpsp_cl - empty_cl > 0) / N;

lpsp_dark  = mean(resample_pablo(N,cc(dark_idx & lpsp_idx)),1);
empty_dark = mean(resample_pablo(N,cc(dark_idx & ~lpsp_idx)),1);
p_dark     = sum( lpsp_dark - empty_dark > 0) / N;

figure(5);subplot(2,2,3);cla
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     'w.','CData',[1,1,1]);%[0 0.4470 0.7410]);
set(gca,'color','none','ycolor','w','xcolor','w')

xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Correlation Coefficient (R)', 'Resampled Means'})
hold on
y = ylim;
plot([1,2],[y(2),y(2)],'w'); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
plot([3,4],[y(2),y(2)],'w'); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
set(gca,'color','none','ycolor','w','xcolor','w')


%repeat for gains
lpsp_cl  = mean(resample_pablo(N,gains(~dark_idx & lpsp_idx)),1);
empty_cl = mean(resample_pablo(N,gains(~dark_idx & ~lpsp_idx)),1);
p_cl     = sum( lpsp_cl - empty_cl > 0) / N;

lpsp_dark  = mean(resample_pablo(N,gains(dark_idx & lpsp_idx)),1);
empty_dark = mean(resample_pablo(N,gains(dark_idx & ~lpsp_idx)),1);
p_dark     = sum( lpsp_dark - empty_dark > 0) / N;

figure(5);subplot(2,2,4); cla; hold on
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     '.','CData',[1,1,1]);%[0 0.4470 0.7410]);


xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Gain', 'Resampled Means'})
hold on
y = ylim;
plot([1,2],[y(2),y(2)],'w'); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
plot([3,4],[y(2),y(2)],'w'); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
set(gca,'color','none','ycolor','w','xcolor','w')
fontsize(gcf,20,'pixels')
%% check the interleaving
c = [min(hsv(2)+.6,1);...
    0.75*hsv(2)];

figure(8); clf
subplot(2,1,1); hold on
for i = unique(ind)'
    scatter(fly_num(ind==i),cc(ind==i),[],c(i,:),'filled');
end
ylabel({'Correlation Coefficient','(bump vs fly vel)'})
legend(group_order,'Location','NortheastOutside','color','none','textcolor','w')
set(gca,'ycolor','w','xcolor','w','color','none')

subplot(2,1,2); hold on
for i = unique(ind)'
    scatter(fly_num(ind==i),gains(ind==i),[],c(i,:),'filled');
end
xlabel('fly num')
ylabel({'Gain','(bump vs fly vel)'})
legend(group_order,'Location','NortheastOutside','color','none','textcolor','w')
set(gca,'ycolor','w','xcolor','w','color','none')

fontsize(gcf,20,'pixels')
set(gcf,'color','none')

%% compare walking stats
dark_mode = true;
if dark_mode
    c = 'w';
else
    c = 'k';
end

mean_f = nan(length(all_data),1);
mean_r = nan(length(all_data),1);
mean_abs = nan(length(all_data),1);
mean_rho = nan(length(all_data),1);

for i = 1:length(all_data)
    mean_f(i) = mean(all_data(i).ft.f_speed);
    mean_r(i) = mean(all_data(i).ft.r_speed);
    mean_abs(i) = mean(abs(all_data(i).ft.r_speed));
    mean_rho(i) = mean(all_data(i).im.rho);
end

figure(9); clf
subplot(2,3,1); hold on
scatter(ind,mean_f,[],c); xticks(unique(ind)); xticklabels(group_order); ylabel('mean F vel (mm/s)')
m = accumarray(ind,mean_f,[],@mean);
s = accumarray(ind,mean_f,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xlim([.5,4.5])

lpsp_cl  = mean(resample_pablo(N,mean_f(~dark_idx & lpsp_idx)),1);
empty_cl = mean(resample_pablo(N,mean_f(~dark_idx & ~lpsp_idx)),1);
lpsp_dark  = mean(resample_pablo(N,mean_f(dark_idx & lpsp_idx)),1);
empty_dark = mean(resample_pablo(N,mean_f(dark_idx & ~lpsp_idx)),1);

p_cl     = sum( lpsp_cl - empty_cl > 0) / N;
p_dark     = sum( lpsp_dark - empty_dark > 0) / N;

subplot(2,3,4);cla; hold on
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     '.',c);
xticks(unique(ind)); xticklabels(group_order); ylabel({'mean F speed (mm/s)','Resampled means'})
y = ylim;
plot([1,2],[y(2),y(2)],c); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)
plot([3,4],[y(2),y(2)],c); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)

subplot(2,3,2); hold on
scatter(ind,mean_r,[],c); xticks(unique(ind)); xticklabels(group_order); ylabel('mean R vel (rad/s)')
m = accumarray(ind,mean_r,[],@mean);
s = accumarray(ind,mean_r,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xlim([.5,4.5])

subplot(2,3,3); hold on
scatter(ind,mean_abs,[],c); xticks(unique(ind)); xticklabels(group_order); ylabel('mean R speed (rad/s)')
m = accumarray(ind,mean_abs,[],@mean);
s = accumarray(ind,mean_abs,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xlim([.5,4.5])

subplot(2,3,6); hold on
scatter(ind,mean_rho,[],c); xticks(unique(ind)); xticklabels(group_order); ylabel('mean vector strength')
m = accumarray(ind,mean_rho,[],@mean);
s = accumarray(ind,mean_rho,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xlim([.5,4.5])

fontsize(gcf,20,'pixels')

if dark_mode
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
        set(tmp(i),'color','none','xcolor','w','ycolor','w')
    end
    set(gcf,'color','none')
    set(gcf,'InvertHardCopy','Off');
else
    set(gcf,'color','white')
end

%% compare optimal lags
opt_lag = nan(length(all_data),1);
opt_cc  = nan(length(all_data),1);

vel_thresh = .2;
bump_thresh = 10;
rho_thresh = .2;

for i = 1:length(all_data)
    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = all_data(i).ft.r_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);

    obj_fun = @(lag) -find_cc(lag,fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh);

    opt_lag(i) = ga(obj_fun,1,[],[],[],[],0,40,[],1);
    opt_cc(i)  = find_cc(opt_lag(i),fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh);
end

%% Plot results
dark_mode = true;
if dark_mode
    c = 'w';
else
    c = 'k';
end


figure(11); clf
subplot(2,2,1); hold on
scatter(ind,opt_lag*fr*1000,c,'filled','markerfacealpha',abs(dark_mode-.3)); ylabel('optimal lag (ms)'); xticks(unique(ind)); xticklabels(group_order); xlim([.5,4.5])
m = accumarray(ind,opt_lag,[],@mean);
s = accumarray(ind,opt_lag,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m*fr*1000,s*fr*1000,'ro')

subplot(2,2,2); hold on
scatter(ind,opt_cc,c); ylabel('optimal corr coeff'); xticks(unique(ind)); xticklabels(group_order); xlim([.5,4.5])
m = accumarray(ind,opt_cc,[],@mean);
s = accumarray(ind,opt_cc,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')


lpsp_cl  = mean(resample_pablo(N,opt_lag(~dark_idx & lpsp_idx)),1);
empty_cl = mean(resample_pablo(N,opt_lag(~dark_idx & ~lpsp_idx)),1);
lpsp_dark  = mean(resample_pablo(N,opt_lag(dark_idx & lpsp_idx)),1);
empty_dark = mean(resample_pablo(N,opt_lag(dark_idx & ~lpsp_idx)),1);

p_cl     = sum( lpsp_cl - empty_cl > 0) / N;
p_dark     = sum( lpsp_dark - empty_dark > 0) / N;

subplot(2,2,3);cla; hold on
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     '.',c);
xticks(unique(ind)); xticklabels(group_order); ylabel({'optimal lag','Resampled means'})
y = ylim;
plot([1,2],[y(2),y(2)],c); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)
plot([3,4],[y(2),y(2)],c); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)

lpsp_cl  = mean(resample_pablo(N,opt_cc(~dark_idx & lpsp_idx)),1);
empty_cl = mean(resample_pablo(N,opt_cc(~dark_idx & ~lpsp_idx)),1);
lpsp_dark  = mean(resample_pablo(N,opt_cc(dark_idx & lpsp_idx)),1);
empty_dark = mean(resample_pablo(N,opt_cc(dark_idx & ~lpsp_idx)),1);

p_cl     = sum( lpsp_cl - empty_cl > 0) / N;
p_dark     = sum( lpsp_dark - empty_dark > 0) / N;

subplot(2,2,4);cla; hold on
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     '.',c);
xticks(unique(ind)); xticklabels(group_order); ylabel({'optimal corr coeff','Resampled means'})
y = ylim;
plot([1,2],[y(2),y(2)],c); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)
plot([3,4],[y(2),y(2)],c); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)

fontsize(gcf,20,'pixels')

if dark_mode
    tmp = get(gcf,'Children');
    for i = 1:length(tmp)
        set(tmp(i),'color','none','xcolor','w','ycolor','w')
    end
    set(gcf,'color','none')
    set(gcf,'InvertHardCopy','Off');
else
    set(gcf,'color','white')
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

function s = process_im(regProduct, im_win, im_type, mask, n_centroid, f0_pct)
    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
    mid             = bwmorph(mask,'remove');
    [y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
    ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
    [y0,x0]         = find(ep,1);
    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
    xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
    centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
    if centroids(1,1) < centroids(end,1)
        centroids = flipud(centroids);
    end
    [~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)

    imgData = squeeze(sum(regProduct,3)); %sum across all z slices
    imgData = imgData - prctile(imgData,1,'all'); %baseline subtract

    imgData = smoothdata(imgData,3,im_type{1},im_win{1});

    imgData_2d = reshape(imgData,[],size(imgData,3)); %reshape the data to be all pixels x all time points
    centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
    end
    f_cluster       = centroid_log * imgData_2d ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;

    s.mu = mu;
    s.rho= smoothdata(rho,1,im_type{2},im_win{2});
    s.z  = zscore_cluster;
    s.alpha = alpha;
    s.imgData = imgData;
end

function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function s = process_im_3d(regProduct, im_win, im_type, mask, n_centroid, f0_pct)

    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length                                           %find the cartesian coordinates of points in the mask, just to find the minimum axis length
    mid             = bwmorph(mask,'remove');
    [y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!
    xq          = linspace(1,length(y_mid),2*(n_centroid) + 1)';                         %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
    centroids   = centroids(2:2:end-1,:);                                                 %take every other so that we dont start at the edges, and all are same size
    if centroids(1,1) < centroids(end,1)
        centroids = flipud(centroids);
    end

    centroids3 = [centroids,(1-centroids(:,1)/max(centroids(:,1)))*size(regProduct,3)*4];


    %[~,idx] = pdist2(whiten(centroids),[whiten(y_mask),whiten(x_mask)],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)
    y_mask3 = repmat(y_mask,size(regProduct,3),1);
    x_mask3 = repmat(x_mask,size(regProduct,3),1);
    z_mask3 = reshape(ones(size(y_mask,1),size(regProduct,3)).*[1:size(regProduct,3)],size(y_mask3,1),1);
    [~,idx] = pdist2(centroids3,[y_mask3,x_mask3,z_mask3.*4],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean distance on whitened positions (so major and minor axes normalized to be same)
    
    imgData = int16(smoothdata(regProduct,4,im_type{1},im_win{1}));
    imgData = imgData - min(imgData,[],'all');
    
    imgData_2d = reshape(imgData,[],size(imgData,4)); %reshape the data to be all pixels x all time points
    centroid_log    = false(n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask3(idx==i),x_mask3(idx ==i),z_mask3(idx ==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});
    
    imgData = squeeze(sum(imgData,3));
    imgData = 256*(imgData-min(imgData,[],'all'))/(max(imgData,[],'all')-min(imgData,[],'all'));

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
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