%% 
close all
clear all

%% load in data
base_dir = 'Z:\pablo\lpsp_p2x2_redo\'; %uigetdir(); %
all_files = dir([base_dir,'\**\*imagingData.mat']);
all_files = natsortfiles(all_files);

idx = cellfun(@(x)(contains(x,'exclusions') | contains(x,'open loop')),{all_files.folder}');
all_files(idx) = [];
%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])

        imgData = img{1};
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
for i = 32:length(all_files)
    if i<length(all_data) && strcmp(all_files(i).folder,all_data(i).meta); continue; end %if we've already processed a file, move on
    if i < length(all_data); all_data(i+1:end+1) = all_data(i:end); end %if the current file to process is missing from the data struct, insert it to the middle by shifting all_data down 1 and rewriting all_data(i)
    % clear img regProduct 
    
    all_data(i).meta = all_files(i).folder;

    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
   
        
    load([tmp2.folder,'\',tmp2.name])
   
    tmp2 = dir([fileparts(all_files(i).folder),'\csv\trialSettings.csv']);
    tmp2 = readtable([tmp2.folder,'\',tmp2.name]);
    
    all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
    all_data(i).im = process_im(img{1}, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).atp = process_im(img{2}, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).ft.stims = ftData_DAQ.stim{1};
    all_data(i).ft.pattern = tmp2.patternPath{1};

    
    
    dr = all_data(i).ft.r_speed;
    %dr(abs(dr)>r_thresh) = dr(abs(dr)>r_thresh) - mean(dr(abs(dr)>r_thresh));
    all_data(i).ft.heading = cumsum(dr)/60;
    

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% go through all trials and add fictrac metadata

idx = false(length(all_data),1);
for i = 1:length(all_data)
    try
    tmp2 = dir([fileparts(all_data(i).meta),'\csv\trialSettings.csv']);
    tmp2 = readtable([tmp2.folder,'\',tmp2.name]);
    all_data(i).ft.pattern = tmp2.patternPath{1};
     dr = all_data(i).ft.r_speed;
    all_data(i).ft.heading = cumsum(dr)/60;
    catch
        idx(i) = true;
    end
end
%% calculate integrative gain
win_sec = 20;
lag = 20;
fr = 60;
win_frames = win_sec*fr;
win = [-10,10];
r_thresh = .1;
rho_thresh = .1;
g = cell(length(all_data),1);
v = cell(length(all_data),1);

tic
for i = 1:length(all_data)
    if isempty(all_data(i).ft); continue; end
    fprintf('processing: %i ',i)

    %extract the best unwrapped estimate of mu by nan-ing low confidence
    %values and unwrapping and re-interp onto fictrac timescale
    m = all_data(i).im.mu;
    m(all_data(i).im.rho<rho_thresh) = nan;
    n_frames = 1:length(m);
    m = interp1(all_data(i).ft.xb(n_frames),unwrap(m),all_data(i).ft.xf);
    m = smoothdata(m,1,"gaussian",60);
    amp = interp1(all_data(i).ft.xb(n_frames),sum(all_data(i).im.d,1),all_data(i).ft.xf);

    %extract the fly's heading (no gain applied) and apply all lags
    h = reshape(all_data(i).ft.heading,[],1);
    m = reshape(m,[],1);

    m = m(lag+1:end);
    h = h(1:end-lag);
    xf= all_data(i).ft.xf(1:end-lag);


    % for each second
    % g_tmp = nan(floor((length(m) - win_frames)/fr),1);
    % v_tmp = nan(floor((length(m) - win_frames)/fr),1);
    
    xt = 0:ceil(max(all_data(i).ft.xf));
    g_tmp = nan(length(xt),1);
    v_tmp = nan(length(xt),1);

    for j = 1:length(g_tmp)
        % f = j*fr;
        % h_tmp = h(f:f+win_frames);
        idx = xf > j+win(1) & xf < j+win(2);
        h_tmp = h(idx);
        m_tmp = m(idx);
        if circ_var(h_tmp) > .1
        % m_tmp = m(f:f+win_frames) - m(f);
        m_tmp = m_tmp - m_tmp(1);
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

    all_data(i).gain.g = g_tmp;
    all_data(i).gain.v = v_tmp;
    all_data(i).gain.xt = xt;
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_data)-i) / 60 / 60)
end

%% plot heading traces
idx = find(cellfun(@(x)(contains(x,'20260128\fly 3')),{all_data.meta})); %,6,'last');
dark_mode = true;
figure(2); clf
c1 = [zeros(256,1),linspace(0,1,256)',zeros(256,1)];
c2 = c1(:,[2,1,3]);
bin_edges = -pi:.1:pi;

for i = 1:length(idx)
    a2 = subplot(length(idx),1,i);
    imagesc(all_data(idx(i)).ft.xb,unwrap(all_data(idx(i)).im.alpha),all_data(idx(i)).atp.d,'AlphaData',1);
    colormap(a2,c2)
    yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})

    a1 = axes('Position',get(gca,  'Position')); 
    imagesc(all_data(idx(i)).ft.xb,unwrap(all_data(idx(i)).im.alpha),all_data(idx(i)).im.z,'AlphaData',1);
    colormap(a1,c1)
    yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})
    xticks([])
    % 
    set(gca,'color','none')
    hold on
    [~,ind] = max(sum(all_data(idx(i)).atp.d,2));
    tmp_alpha = unwrap(all_data(idx(i)).im.alpha);
    scatter(0,tmp_alpha(ind),'r*')

    if contains(all_data(idx(i)).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
    a = plot(all_data(idx(i)).ft.xf,-all_data(idx(i)).ft.cue,c);
    a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(all_data(idx(i)).ft.xb,all_data(idx(i)).im.mu,'w');
    a.YData(abs(diff(a.YData))>pi) = nan;
    xticks([]);yticks([])

    pos = get(gca,'Position');
    pos1 = [pos(1),pos(2)+pos(4),pos(3),.03];
    a3 = axes('Position',pos1,'color','none'); 
    plot(all_data(idx(i)).ft.xb,2*sum(all_data(idx(i)).atp.f,1)/max(sum(all_data(idx(i)).atp.f,1)),'r','linewidth',2)   
    hold on
    % plot(all_data(idx(i)).ft.xf,abs(all_data(idx(i)).ft.r_speed)/max(abs(all_data(idx(i)).ft.r_speed)))
    % plot(all_data(idx(i)).ft.xf,abs(all_data(idx(i)).ft.f_speed)/max(abs(all_data(idx(i)).ft.f_speed)))
%    plot(all_data(idx(i)).gain.xt,all_data(idx(i)).gain.g)
    yticks([]);xticks([])
    axis tight
    set(gca,'Color','none')
    linkaxes([a1,a2,a3],'x')
    linkaxes([a1,a2],'y')

    pos2 = [pos(1)+pos(3)+.01,pos(2),.05,pos(4)];
    ax = axes('Position',pos2,'Color','none','XAxisLocation','top');
    offset = circ_dist(-all_data(idx(i)).ft.cue,interp1(all_data(idx(i)).ft.xb,unwrap(all_data(idx(i)).im.mu),all_data(idx(i)).ft.xf));
    histogram(offset,bin_edges,'Orientation','horizontal','edgeColor','none','FaceColor',c)
    box(ax,'off')
    ax.YAxisLocation =  'right'; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};

    if i == length(idx)
        pos = get(a3,'Position');
        legend(a3,'atp','r speed','f speed','Location','NorthwestOutside')
        set(a3,'Position',pos)

        pos = get(a1,'Position');
        legend(a1,'heading','pva','Location','SouthwestOutside')
        set(a1,'Position',pos)
    end

end


ax = get(gcf,'Children');
if dark_mode
    set(gcf,'color','none','InvertHardcopy','off')
    
    for i = 1:length(ax)
        if contains(class(ax(i)),'Axes')
            ax(i).Title.Color = 'w';
            set(ax(i),'Color','none','ycolor','w','xcolor','w')
        else
            ax(i).TextColor = 'w';
        end
    end
end

%% Show gains for CL, Dark, Stimmed, prestimmed.
bin_edges = 0:.1:5;
stim_decay = 20;

stim_idx = false(length(all_data),1);
post_idx = false(length(all_data),1);
dark_idx = false(length(all_data),1);
empty_idx= false(length(all_data),1);

last_str = '';

for i = 1:length(all_data)
    dark_idx(i) = contains(all_data(i).ft.pattern,'background');
    empty_idx(i) = contains(all_data(i).meta,'control');
    tmp_ind = strfind(all_data(i).meta,'fly');
    fly_str = all_data(i).meta(tmp_ind-9:tmp_ind+5);
    
    if strcmp(fly_str,last_str) && post_idx(i-1) %if it's the same fly and the last trial was post stim, this trial is post stim
        post_idx(i) = true;
    end
    last_str = fly_str; % now make this fly the last fly 
    

    tmp = all_data(i).gain.g; %hold gain temporarily

    [~,loc] = findpeaks(sum(all_data(i).atp.d,1),'MinPeakProminence',2); %find the atp ejections

    if ~isempty(loc)   % if there are any stims in this trial
        stim_idx(i) = true; %set stims to be true
        post_idx(i) = true; %set the post stim idx to be true
        stim_time = all_data(i).ft.xb(loc); % remove all gains that fall within the stim window
        tmp_idx = any(all_data(i).gain.xt>stim_time & all_data(i).gain.xt<stim_time+stim_decay,1);
        tmp(tmp_idx) = [];
    end

    g{i} = tmp(~isnan(tmp) & tmp ~= 0);
end

figure(11); clf
t = tiledlayout(2,2);

group_idx = post_idx + 2*dark_idx;
group_names = {'pre stim, CL','post stim CL','pre stim Dark','post stim Dark'};
for i = 0:3
nexttile; hold on
histogram(cell2mat(g(empty_idx & group_idx==i)),'Normalization','probability','BinEdges',bin_edges,'FaceColor',[0,.5,1])
histogram(cell2mat(g(~empty_idx & group_idx==i)),'Normalization','probability','BinEdges',bin_edges,'FaceColor',[1,.5,0])
title(group_names{i+1})
end
xlabel(t,'Integrative Gain')
ylabel(t,'Probability')
title(t,sprintf('Gain by condition, excluding %is post stim',stim_decay))

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

function h = plotsem(t,x,c)
    m = mean(x,1,"omitnan");
    s = std(x,[],1,'omitnan') ./ sqrt(sum(~isnan(x),1));
    t = reshape(t,1,[]);
    
    h = patch([t,fliplr(t)],[m+s,m-s],c);
end