%% clear all
% clear all
% close all

%% find path to all relevant files
base_dir = uigetdir(); %('Z:\pablo\lpsp_p2x2\todo\');
all_files = dir([base_dir,'\**\*imagingData*.mat']);
all_files = natsortfiles(all_files);

%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
       
        ch1 = squeeze(mean(img{1},3));
        ch2 = squeeze(mean(img{2},3));
        
        figure(1); clf; imagesc(mean(ch1,3)); axis equal tight; colormap(parula); drawnow;
        mask = roipoly();
        save([fileparts(all_files(i).folder),'\mask.mat'],'mask')
    end
end

%% process and store all values
ft_type= 'gaussian'; %the type of smoothing for fictrac data
ft_win = 30; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'movmean','gaussian'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {3,5};

n_centroid = 16;
f0_pct = 7;
%all_data = struct();

tic
for i = 5:length(all_files)
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    load([tmp2.folder,'\',tmp2.name])

    regProduct = img{1};

    all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
    all_data(i).im = process_im(img{1}, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    if ~ismember('xb',fieldnames(all_data(i).ft))
        xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),size(all_data(i).im.d,2));
    end

    all_data(i).atp = process_im(img{2}, im_win, im_type, mask, n_centroid, f0_pct);

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% plot
i = 3;
alpha = unwrap(all_data(i).im.alpha);
fr = mean(diff(all_data(i).ft.xb));

[~,xloc] = findpeaks(sum(all_data(i).atp.d,1),'MinPeakProminence',9,'MinPeakDistance',30/fr);
[~,yloc] = max(sum(all_data(i).atp.d,2));

figure(1); clf
subplot(2,1,1); 
imagesc(all_data(i).ft.xb,alpha,all_data(i).atp.d)
colormap(gca,[linspace(0,1,256)',zeros(256,2)])
hold on
scatter(all_data(i).ft.xb(xloc),alpha(yloc),'r*');
axis tight
title(all_data(i).meta)

subplot(2,1,2);
imagesc(all_data(i).ft.xb,alpha,all_data(i).im.d)
colormap(gca,'parula')
hold on
scatter(all_data(i).ft.xb(xloc),alpha(yloc),'r*');
%a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi)=nan;
%a = plot(all_data(i).ft.xb,all_data(i).im.mu+2*pi,'w'); a.YData(abs(diff(a.YData))>pi)=nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'m','linewidth',2); a.YData(abs(diff(a.YData))>pi)=nan;
axis tight
xlabel('time (s)')
pos = get(gca,'Position');
a3 = axes('Position',[pos(1),pos(2)+pos(4),pos(3),.03],'color','none'); 
plot(all_data(i).ft.xb,sum(all_data(i).atp.f,1)/max(sum(all_data(i).atp.f,1)),'r','linewidth',2)
xticks([]); yticks([])

linkaxes(get(gcf,'Children'),'x')
axis tight

fontsize(gcf,20,'pixels')

    ax = get(gcf,'Children');
    if dark_mode
        set(gcf,'color','none','InvertHardcopy','off')
        
        for i = 1:length(ax)
            if contains(class(ax(i)),'Axes')
                ax(i).Title.Color = 'w';
                set(ax(i),'Color','none','ycolor','w','xcolor','w')
            else
                ax(i).Color = 'w';
            end
        end
    end

%% extract mu aligned pulses
win_start = -30;
win_end = 60;

c_pulses = {};
m_pulses = {};
a_pulses = {};
d_pulses = {};
z_pulses = {};
f_pulses = {};
dark_idx = {};
exp_idx  = {};
right_idx= {};
first_idx= {};
fr = 10;
tmp_x   = linspace(0,300,300*fr);    
tmp_win = floor(win_start*fr):ceil(win_end*fr); %this is the additive index to the frames to extract for a given pulse window
tmp_t   = tmp_win/10;

counter = 0;
for i = 1:length(all_data)
    tmp_f   = interp1(all_data(i).ft.xb,all_data(i).im.f',tmp_x,'linear','extrap')';%interpolate everything into the same framerate
    tmp_d   = interp1(all_data(i).ft.xb,all_data(i).im.d',tmp_x,'linear','extrap')';
    tmp_z   = interp1(all_data(i).ft.xb,all_data(i).im.z',tmp_x,'linear','extrap')';
    tmp_a   = interp1(all_data(i).ft.xb,all_data(i).atp.d',tmp_x,'linear','extrap')';
    
    tmp_m   = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),tmp_x,'linear','extrap')'; tmp_m = mod(tmp_m,2*pi); tmp_m(tmp_m>pi) = tmp_m(tmp_m>pi) - 2*pi;
    tmp_cue = interp1(all_data(i).ft.xf,unwrap(all_data(i).ft.cue),tmp_x,'linear','extrap')';tmp_cue = mod(tmp_cue,2*pi); tmp_cue(tmp_cue>pi) = tmp_cue(tmp_cue>pi) - 2*pi; 
   
    tmp_atp = sum(tmp_a,2);
    tmp_right = sum(tmp_atp(1:end/2,:),'all') > sum(tmp_atp(end/2:end,:),'all');
    [~,loc] = findpeaks(smoothdata(max(tmp_a,[],1),'movmean',5),'MinPeakProminence',5,'MinPeakDistance',50*fr);
    %fr = mean(diff(all_data(i).ft.xb)); %find the new framerate (amount of time per frame)
    %tmp_win = floor(win_start/fr):ceil(win_end/fr); %this is the additive index to the frames to extract for a given pulse window
    %findpeaks(smoothdata(max(tmp_a,[],1),'movmean',5),'MinPeakProminence',1.5,'MinPeakDistance',50/fr);

    if i==1 || ~strcmp(all_data(i).meta(1:49),last_fly)
        last_fly = all_data(i).meta(1:49);
        first_log = true;
    end

    if numel(loc) < 2
        continue
    end

    for j = loc

        n = length(all_data(i).im.alpha);
        pad1 = nan(n,sum(j+tmp_win < 1)); %create padding for images
        pad2 = nan(n,sum(j+tmp_win > length(all_data(i).ft.xb)));
        tmp_idx = j+tmp_win((size(pad1,2)+1):(end-size(pad2,2))); %extract all of the appropriate frames. if we exceed the window in either direction, just fill it with the nearest calculated gain. a bit hacky
        
        counter = counter+1;
        d_pulses{counter} = [pad1,tmp_d(:,tmp_idx),pad2];
        z_pulses{counter} = [pad1,tmp_z(:,tmp_idx),pad2];       
        f_pulses{counter} = [pad1,tmp_f(:,tmp_idx),pad2];       
        a_pulses{counter} = [pad1,tmp_a(:,tmp_idx),pad2];
        m_pulses{counter} = [pad1(1,:),tmp_m(tmp_idx)',pad2(1,:)];
        c_pulses{counter} = [pad1(1,:),tmp_cue(tmp_idx)',pad2(1,:)];
        dark_idx{counter} = contains(all_data(i).meta,'dark');
        exp_idx{counter} = ~contains(all_data(i).meta,'control');
        right_idx{counter} = tmp_right;
        first_idx{counter} = first_log;
        
        if first_log
            first_log = false;
        end

    end
end

dark_idx = logical(cell2mat(dark_idx));
exp_idx = logical(cell2mat(exp_idx));
right_idx = logical(cell2mat(right_idx));
first_idx = logical(cell2mat(first_idx));

alpha = all_data(1).im.alpha;
[~,ind] = min(abs(tmp_win*fr)); %find the index to the frame of the atp peak
ind = round(ind- .5/fr):ind;

atp_loc = cellfun(@(x)(find(sum(x,2,'omitnan')==max(sum(x,2,'omitnan')),1,'first')),a_pulses);
atp_loc(atp_loc > n/2) = atp_loc(atp_loc>n/2) - n/2;
mu_loc  = cellfun(@(x)(find(abs(mean(x(ind),'omitnan')-alpha)==min(abs(mean(x(ind),'omitnan')-alpha)),1,'first')),m_pulses);

overlap_idx = abs(atp_loc - mu_loc) < 4;

% plot results
group_idx = exp_idx + 2*right_idx;
group_labels = {'con (left)','exp (left)','con (right)','exp (right)'};
for i = unique(group_idx)
    figure(i+1); clf
    tmp_ind = find(group_idx==i);

    rows = ceil(sqrt(length(tmp_ind)));
    cols = ceil(length(tmp_ind)/rows);

    for j = 1:length(tmp_ind)
        subplot(rows,cols,j)
        tmp_win = linspace(win_start,win_end,size(d_pulses{tmp_ind(j)},2));
    
        imagesc(tmp_win,unwrap(alpha),d_pulses{tmp_ind(j)})
        hold on
        %a = plot(tmp_win,m_pulses{tmp_ind(j)},'k');  a.YData(abs(diff(a.YData))>pi) = nan;
        a = plot(tmp_win,-c_pulses{tmp_ind(j)},'m','linewidth',2);  a.YData(abs(diff(a.YData))>pi) = nan;
    end
    sgtitle(group_labels{i+1})

    ax = get(gcf,'Children');
    if dark_mode
        set(gcf,'color','none','InvertHardcopy','off')
        
        for i = 1:length(ax)
            if contains(class(ax(i)),'Axes')
                ax(i).Title.Color = 'w';
                set(ax(i),'Color','none','ycolor','w','xcolor','w')
            else
                ax(i).Color = 'w';
            end
        end
    end
end



%% plot sweeps
dark_mode = true;
if dark_mode
    c = 'w';
else
    c = 'k';
end


figure(5); clf
subplot(1,2,1); hold on
a = plot_sem(gca,tmp_t',cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(exp_idx & right_idx)','UniformOutput',false))); a.FaceColor = 'r';
a = plot_sem(gca,tmp_t',cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(exp_idx & ~right_idx)','UniformOutput',false))); a.FaceColor = 'c';
%a = plot_sem(gca,tmp_t',-cell2mat(cellfun(@(x)(unwrap(x-x(1))),c_pulses(exp_idx)','UniformOutput',false))); a.FaceColor = 'g';


plot([win_start,win_end],[0,0],':','Color',c)
scatter(-.5,0,100,'r*')
title('experimental')
xlabel('time post stim (s)')
ylabel('unwrapped bump position (rad)')
axis tight
set(gca,'YDir','reverse')

subplot(1,2,2); hold on
a = plot_sem(gca,tmp_t',cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(~exp_idx & right_idx)','UniformOutput',false))); a.FaceColor = 'r';
a = plot_sem(gca,tmp_t',cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(~exp_idx & ~right_idx)','UniformOutput',false))); a.FaceColor = 'c';
%a = plot_sem(gca,tmp_t',-cell2mat(cellfun(@(x)(unwrap(x-x(1))),c_pulses(~exp_idx)','UniformOutput',false))); a.FaceColor = 'g';
plot([win_start,win_end],[0,0],':','Color',c)
scatter(-.5,0,100,'r*')
title('control')
xlabel('time post stim (s)')
ylabel('unwrapped bump position (rad)')
h = legend('right','left','Location','Southeast','Color',(1-dark_mode)*[1,1,1],'TextColor',c)
axis tight
set(gca,'YDir','reverse')

linkaxes(get(gcf,'Children'))
fontsize(gcf,20,'pixels')

 ax = get(gcf,'Children');
    if dark_mode
        set(gcf,'color','none','InvertHardcopy','off')

        for i = 1:length(ax)
            if contains(class(ax(i)),'Axes')
                ax(i).Title.Color = 'w';
                set(ax(i),'Color','none','ycolor','w','xcolor','w')
            else
                ax(i).Color = 'w';
            end
        end
    end

%% assess changes to fluorescence
group_idx = exp_idx + 2*right_idx + 4*overlap_idx;
group_labels = {'con (left)','exp (left)','con (right)','exp (right)'};
group_labels = [group_labels,group_labels];

d_pulses_aligned = cell(size(d_pulses));
z_pulses_aligned = cell(size(z_pulses));
f_pulses_aligned = cell(size(f_pulses));
a_pulses_aligned = cell(size(a_pulses));

for i = 1:length(d_pulses)
    if i == 158
        a=1;
    end
    [~,k] = max(sum(a_pulses{i},2,'omitnan'));

    tmp1 = circshift(d_pulses{i}(1:n/2,:),-k+n/4,1);
    tmp2 = circshift(d_pulses{i}(n/2+1:end,:),-k+n/4,1);
    d_pulses_aligned{i} = [tmp1;tmp2];

    tmp1 = circshift(z_pulses{i}(1:n/2,:),-k+n/4,1);
    tmp2 = circshift(z_pulses{i}(n/2+1:end,:),-k+n/4,1);
    z_pulses_aligned{i} = [tmp1;tmp2];

    tmp1 = circshift(f_pulses{i}(1:n/2,:),-k+n/4,1);
    tmp2 = circshift(f_pulses{i}(n/2+1:end,:),-k+n/4,1);
    f_pulses_aligned{i} = [tmp1;tmp2];

    tmp1 = circshift(a_pulses{i}(1:n/2,:),-k+n/4,1);
    tmp2 = circshift(a_pulses{i}(n/2+1:end,:),-k+n/4,1);
    a_pulses_aligned{i} = [tmp1;tmp2];
end

figure(6); clf
for i = unique(group_idx)
    subplot(2,4,i+1)
    imagesc(tmp_t+.5,1:n,mean(cell2mat(permute(d_pulses_aligned(group_idx==i),[3,1,2])),3,'omitnan'))
    hold on
    if contains(group_labels{i+1},'left'); scatter(0,3*n/4,'r*')
    else  ;  scatter(0,n/4,'r*'); end
    plot(xlim,[n/2,n/2]+.5,'w','linewidth',2)
    colormap(parula)
    title([group_labels(i+1),sprintf('n=%d',sum(group_idx==i))]);
    ylabel({'zscore dff','aligned'})
end
linkaxes(get(gcf,'Children'))

% figure(3); clf
% for i = unique(group_idx)
%     tmp_mat = cell2mat(permute(f_pulses_aligned(group_idx==i),[3,1,2]));
%     tmp_mat = mean(tmp_mat(:,tmp_t<-10,:),2,'omitnan') - mean(tmp_mat(:,tmp_t>30,:),2,'omitnan');
%     tmp_mat = permute(tmp_mat,[3,1,2]);
% 
%     subplot(2,1,floor(i/2)+1)
%     a = plot_sem(gca,(1:size(tmp_mat,2))',tmp_mat); if mod(i,2) == 0; a.FaceColor = 'b'; end
%     title(group_labels(i+1));
%     ylabel({'zscore dff','aligned'})
% end
% linkaxes(get(gcf,'Children'))

%% assess the bump mobility before and after ATP
uni_flies = unique(cellfun(@(x)(x(1:50)),{all_data.meta},'UniformOutput',false));
num_flies = length(uni_flies);

gains = cell(num_flies,1);
con_idx = false(num_flies,1);

figure(10);clf

for i = 1:num_flies
    idx = contains({all_data.meta},uni_flies{i});
    tmp_gains = nan(sum(idx),1);
    ind = find(idx);
    con_idx(i) = contains(all_data(ind(1)).meta,'control');

    for j = 1:length(ind)
        mu  = unwrap(all_data(ind(j)).im.mu);
        cue = -unwrap(all_data(ind(j)).ft.cue);
        tmp_gains(j) = (mu(end) - mu(1)) / (cue(end) - cue(1));
        subplot(6,4,i); hold on
        if con_idx(i)
            plot(all_data(ind(j)).ft.xb,mu,'r')
        else
            plot(all_data(ind(j)).ft.xb,mu,'b')
        end
        plot(all_data(ind(j)).ft.xf,cue,'k')
    end
    
    
    gains{i} = tmp_gains;
end

figure(11); 
subplot(1,2,1); hold on
cellfun(@(x)(plot(x,'r')),gains(con_idx))
subplot(1,2,2); hold on
cellfun(@(x)(plot(x,'b')),gains(~con_idx))

%% look at fluoresnce in each trial
figure(12); clf
for i = 1:num_flies
    idx = contains({all_data.meta},uni_flies{i});
    tmp_gains = nan(sum(idx),1);
    ind = find(idx);
    con_idx(i) = contains(all_data(ind(1)).meta,'control');

    for j = 1:length(ind)
        subplot(num_flies,length(ind),(i-1)*(length(ind)) + j); hold on
        imagesc(all_data(ind(j)).ft.xb,all_data(ind(j)).im.alpha,all_data(ind(j)).im.d)
        a = plot(all_data(ind(j)).ft.xb,all_data(ind(j)).im.mu,'m'); a.YData(abs(diff(a.YData))>pi) = nan;
        axis tight
        xticks([]); yticks([]);
        ylabel(round(gains{i}(j),2))
    end
end


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