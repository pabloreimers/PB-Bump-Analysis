%% 
%close all
%clear all

%% load in data
base_dir = 'Z:\pablo\hackathon_gain_change\'; %uigetdir(); %
all_files = dir([base_dir,'\**\*imagingData.mat']);
all_files = natsortfiles(all_files);
%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        
        % if exist('regProduct','var')
        % imgData = squeeze(sum(regProduct,3));
        % end
        % clear regProduct

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
r_thresh = .1;
rho_thresh = .1;

%all_data = struct();

tic
for i = 1:length(all_files)
    if i<length(all_data) && strcmp(all_files(i).folder,all_data(i).meta); continue; end %if we've already processed a file, move on
    if i < length(all_data); all_data(i+1:end+1) = all_data(i:end); end %if the current file to process is missing from the data struct, insert it to the middle by shifting all_data down 1 and rewriting all_data(i)
    % clear img regProduct 
    
    all_data(i).meta = all_files(i).folder;

    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    if isempty(tmp2) %if this is a Noah Trial
        tmp2 = dir([fileparts(all_files(i).folder),'\*pablo_ft.mat']); 
        try load([tmp2.folder,'\',tmp2.name]); catch continue; end
        all_data(i).ft = ft;
        all_data(i).ft.xb = all_data(i).ft.xb(1:length(all_data(i).im.mu));
        field_names = fieldnames(all_data(i).ft);
        for j = 1:length(field_names)
            try all_data(i).ft.(field_names{j}) = reshape(all_data(i).ft.(field_names{j}),[],1); end
        end
        all_data(i).ft.pattern = 'berg4';
        dr = -all_data(i).ft.r_speed;
        
    else %ten this is a Pablo Trial
        load([tmp2.folder,'\',tmp2.name])
        all_data(i).ft = process_ft(ftData_DAQ, ft_win, ft_type);
        tmp2 = dir([fileparts(all_files(i).folder),'\csv\trialSettings.csv']);
        tmp2 = readtable([tmp2.folder,'\',tmp2.name]);
        all_data(i).ft.pattern = tmp2.patternPath{1};
        dr = all_data(i).ft.r_speed;
        % imgData = squeeze(sum(regProduct,3));
    end
    
    
    
    
    % tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    % load([tmp2.folder,'\',tmp2.name])
   
    
    all_data(i).im = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct);
    

    

    try
     all_data(i).ft.stims = ftData_DAQ.stim{1};
    end
    
    
    
    %dr(abs(dr)>r_thresh) = dr(abs(dr)>r_thresh) - mean(dr(abs(dr)>r_thresh));
    all_data(i).ft.heading = cumsum(dr)/60;
    

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
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
for i = 33:length(all_data)
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

%% calculate sliding gain
win = [-10,10];
lag = 5;
fr = 60;
r_thresh = .1;
rho_thresh = .1;
g = cell(length(all_data),1);
v = cell(length(all_data),1);

tic
for i = 33:length(all_data)
    if isempty(all_data(i).ft); continue; end
    fprintf('processing: %i ',i)

    %extract the best unwrapped estimate of mu by nan-ing low confidence
    %values and unwrapping and re-interp onto fictrac timescale
    m = all_data(i).im.mu;
    %m(all_data(i).im.rho<rho_thresh) = nan;
    n_frames = 1:length(m);
    m = interp1(all_data(i).ft.xb(n_frames),unwrap(m),all_data(i).ft.xf);
    
    dm = gradient(m) * 60;
    dr = all_data(i).ft.r_speed;
    df = all_data(i).ft.f_speed;

    xf = all_data(i).ft.xf;
    xb = all_data(i).ft.xb;

    g_tmp = nan(length(xb),1);

    dm = dm(1:end-lag);
    dr = dr(lag:end);
    df = df(lag:end);
    xf = xf(1:end-lag);

    for t = 1:length(xb)
        idx = (xf > xb(t)+win(1)) & (xf < xb(t)+win(2));
        dr_tmp = dr(idx);
        df_tmp = df(idx);
        dm_tmp = dm(idx);

        % if all(all_data(i).ft.cue_brightness(idx)==0)
        %     a=1;
        % end
        %speed_corr = (dr_tmp) \ (df_tmp);
        %dr_tmp = dr_tmp - df_tmp'*speed_corr;
        

        % [r,lags]  = xcorr(dr_tmp,dm_tmp,30);
        % lag = lags(r==max(r));
        % 
        % if lag > 0
        %     dr_tmp = dr_tmp(lag:end);
        %     df_tmp = df_tmp(lag:end);
        %     dm_tmp = dm_tmp(1:end-lag+1);
        % else
        %     dr_tmp = dr_tmp(1:end+lag);
        %     df_tmp = df_tmp(1:end+lag);
        %     dm_tmp = dm_tmp(1-lag:end);
        % end

        a = [dr_tmp,df_tmp] \ dm_tmp;
        if ~isempty(dm_tmp) 
            %g_tmp(t) =  dr_tmp(~isnan(dm_tmp)) \ dm_tmp(~isnan(dm_tmp));
            g_tmp(t) = a(1);
        end
    end
    
    all_data(i).gain.inst_g = g_tmp;
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_data)-i) / 60 / 60)
end

%% calculate sliding variance (casting)
win = [-5,5];
lag = 20;
fr = 60;
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
    n_frames = 1:length(m);
    m = interp1(all_data(i).ft.xb(n_frames),unwrap(m),all_data(i).ft.xf);
    h = reshape(all_data(i).ft.heading,[],1);
    c = all_data(i).ft.cue;
    
    xf = all_data(i).ft.xf;
    xb = all_data(i).ft.xb;

    mv = nan(length(xb),1);
    hv = nan(length(xb),1);
    cv = nan(length(xb),1);

    for t = 1:length(xb)
        idx = (xf > xb(t)+win(1)) & (xf < xb(t)+win(2));
        h_tmp = h(idx);
        c_tmp = c(idx);
        m_tmp = m(idx);
        [r,lags]  = xcorr(h_tmp,h_tmp,30);
        lag = lags(r==max(r));
        
        if lag > 0
            h_tmp = h_tmp(lag:end);
            c_tmp = c_tmp(lag:end);
            m_tmp = m_tmp(1:end-lag+1);
        else
            h_tmp = h_tmp(1:end+lag);
            c_tmp = c_tmp(1:end+lag);
            m_tmp = m_tmp(1-lag:end);
        end
        if ~isempty(m_tmp)
            hv(t) =  circ_var(h_tmp);
            cv(t) =  circ_var(c_tmp);
            mv(t) =  circ_var(m_tmp);
        end
    
    end
    
    all_data(i).gain.hv = hv;
    all_data(i).gain.cv = cv;
    all_data(i).gain.mv = mv;
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_data)-i) / 60 / 60)
end

%% create figure to show example
i = 40;
binedges = 0:.05:5;

figure(1); clf
a1 = subplot(3,1,1);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c); a.YData(abs(diff(a.YData))>pi) = nan;
idx = round(all_data(i).ft.cue,4) == -.2945;
h = mod(all_data(i).ft.heading,2*pi) - pi;
a = plot(all_data(i).ft.xf,h,'r'); a.YData([abs(diff(a.YData))'>pi ; 0] | ~idx) = nan;

a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)
xlabel('time (s)')

a2 = subplot(6,1,3); hold on
a=plot(all_data(i).ft.xf,circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf))); a.YData(abs(diff(a.YData))>pi) =nan;
% try patch([all_data(i).ft.xf,flipud(all_data(i).ft.xf)],...
%         [(max(ylim)-min(ylim))*(1-all_data(i).ft.cue_brightness/max(all_data(i).ft.cue_brightness))+min(ylim)),...
%         min(ylim)*ones(size(all_data(i).ft.xf))]); end

try tmp = 2*pi*(all_data(i).ft.cue_brightness==0)-pi;tmp(1) = -pi; tmp(end) = -pi; patch(all_data(i).ft.xf,tmp,'k','FaceAlpha',.1,'EdgeColor','none');end
ylabel('offset')

a3 = subplot(6,1,4); hold on
%plot(all_data(i).ft.xf,all_data(i).ft.f_speed)
plot(all_data(i).gain.xt,all_data(i).gain.g); hold on; plot(xlim,[.8,.8],':k'); plot(xlim,[1.6,1.6],':k')
%plot(all_data(i).ft.xb,all_data(i).gain.hv); hold on; plot(xlim,[.8,.8],':k'); plot(xlim,[1.6,1.6],':k')
%plot(all_data(i).ft.xb,all_data(i).gain.mv); hold on; plot(xlim,[.8,.8],':k'); plot(xlim,[1.6,1.6],':k')
plot(all_data(i).ft.xb,all_data(i).gain.inst_g)
try plot(all_data(i).ft.xf,all_data(i).ft.gain); catch; plot(all_data(i).ft.xf(1:length(all_data(i).gain.vr)),all_data(i).gain.vr); end
%try plot(all_data(i).ft.xf,all_data(i).ft.cue_brightness/150); end
%plot(all_data(i).ft.xf,all_data(i).ft.gain)
%a=plot(all_data(i).ft.xf,circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,all_data(i).im.mu,all_data(i).ft.xf))); a.YData(abs(diff(a.YData))>pi) =nan;
ylabel('integrative gain')


linkaxes([a1,a2,a3],'x')
axis tight
ylim([0,5])

subplot(3,2,5); hold on
h = histogram(all_data(i).gain.inst_g,'BinEdges',binedges,'FaceAlpha',.5,'Normalization','probability');
h = histogram(all_data(i).gain.g,'BinEdges',binedges,'FaceAlpha',.5,'Normalization','probability');
%plot(binedges(1:end-1),h.Values,'k','Linewidth',2)
%h = histogram(all_data(i).gain.g(end/2:end),'BinEdges',binedges,'FaceAlpha',.5);
%plot(binedges(1:end-1),h.Values,'k','Linewidth',2)

xlabel('integrative gain')
ylabel('counts')
legend('linear comb','integrative')
%legend('early trial','late trial')

subplot(3,2,6); 
scatter(all_data(i).gain.g,all_data(i).gain.v,'filled','MarkerFaceAlpha',.5)
xlabel('integrative gain')
ylabel('func value')

%% work on pca dark decoder single fly
i = 40;
dark_idx = all_data(i).ft.cue_brightness' ==0 & all_data(i).ft.f_speed'>.5;

figure(2); clf
scatter(all_data(i).ft.f_speed(dark_idx)',all_data(i).ft.r_speed(dark_idx),'.')

%% show early and late casting variance
figure(4); clf; 
subplot(1,2,1); hold on
for i = find(test_idx)'
    tmp_r = all_data(i).gain.mv ./ all_data(i).gain.hv;
    tmp_r(all_data(i).gain.hv < .1) = nan;
    plot([mean(tmp_r(1:end/2),'omitnan'),mean(tmp_r(end/2:end),'omitnan')],'o-k')
end
xlim([0,3])
plot(xlim,[1,1],":k",'linewidth',4)
xticks([1,2]); xticklabels({'early trial','late trial'})
title('mean(bump variance / heading variance)')

subplot(1,2,2); hold on
for i = find(test_idx)'
    tmp_r = all_data(i).gain.mv ./ all_data(i).gain.cv;
    tmp_r(all_data(i).gain.cv < .1) = nan;
    plot([mean(tmp_r(1:end/2),'omitnan'),mean(tmp_r(end/2:end),'omitnan')],'o-k')
end
xlim([0,3])
plot(xlim,[1,1],":k",'linewidth',4)
xticks([1,2]); xticklabels({'early trial','late trial'})
title('mean(bump variance / cue variance)')
fontsize(gcf,40,'pixels')
%%
figure(2); clf; hold on

m = all_data(i).im.mu;
m(all_data(i).im.rho<rho_thresh) = nan;
m = interp1(all_data(i).ft.xb,unwrap(m),all_data(i).ft.xf);
dm = gradient(m) * 60;
dr = -all_data(i).ft.r_speed;

dm = dm - mean(dm,'omitnan');
dr = dr - mean(dr,'omitnan');
idx1 = abs(dr) > r_thresh & all_data(i).ft.cue_brightness == 0 & all_data(i).ft.xf < median(all_data(i).ft.xf) & ~isnan(dm);
idx2 = abs(dr) > r_thresh & all_data(i).ft.cue_brightness == 0 & all_data(i).ft.xf > median(all_data(i).ft.xf) & ~isnan(dm);
scatter(dr(idx1),dm(idx1),[],'r','filled','MarkerFaceAlpha',.8)
scatter(dr(idx2),dm(idx2),[],'b','filled','MarkerFaceAlpha',.8)
xlabel('fly speed')
ylabel('bump speed')
legend(sprintf('dark 1 (gain: %.2f)',dr(idx1)'\dm(idx1)'),sprintf('dark 2 (gain: %.2f)',dr(idx2)'\dm(idx2)'))
fontsize(gcf,20,'pixels')

%% look at bump offsets
edges = -pi:.1:pi;
rows = ceil(sqrt(length(all_data)));
cols = ceil(length(all_data)/rows);

figure(1); clf
t = tiledlayout(rows,cols);
pop_off = nan(length(all_data),2);
for i = 1:length(all_data)
    nexttile
    if contains(all_data(i).meta,'ramp'); c = 'r'; else; c = 'b'; end
    offset = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,all_data(i).im.mu,all_data(i).ft.xf));
    h = histogram(offset,'FaceColor',c,'EdgeColor','none','BinEdges',edges,'Normalization','probability');

    [x,y] = pol2cart(h.BinEdges(1:end-1),h.Values);
    [pop_off(i,1),pop_off(i,2)] = cart2pol(mean(x),mean(y));
    %histogram(all_data(i).ft.cue,'FaceColor',c,'EdgeColor','none')
    
    %histogram(all_data(i).ft.f_speed,'BinEdges',-2:.5:25,'EdgeColor','none')
    yticks([])
    xticks([])
    hold on
    %plot([10,10],ylim,'r')
    plot([0,0],ylim,':k')
end

xlabel(t,'Cue-Bump Offset')
ylabel(t,'Counts')
set(gcf,'color','w')
fontsize(gcf,40,'pixels')
%fontsize(gca,20,'pixels')
xticks([-pi,0,pi]); xticklabels({'-\pi',0,'\pi'})

figure(2); clf
polarscatter(pop_off(berg_idx,1),pop_off(berg_idx,2),400,'rx')
hold on
polarscatter(pop_off(~berg_idx,1),pop_off(~berg_idx,2),400,'bo')
title('Cue-Bump Offset Vector Strength (pva)')
legend('berg1','berg4')
fontsize(gcf,40,'pixels')
%% save the VR gain for each trial

for i = 34:length(all_data)
    i
    t = 1:round(max(all_data(i).ft.xf));
    g = nan(size(t));
    r = all_data(i).ft.r_speed(1:end-2);
    c = -gradient(all_data(i).ft.cue(3:end))*60;
    c(abs(c)>50) = nan;
    for j = 1:length(t)
        idx = all_data(i).ft.xf(1:end-2) > j-1 & all_data(i).ft.xf(1:end-2) < j+5 & abs(r) > .5 & abs(c) > .5;
        
        
        if sum(idx) == 0; g(j) = nan;else; g(j) = r(idx) \ c(idx); end
    end
    g = interp1(t,g,all_data(i).ft.xf);
    all_data(i).gain.vr = g;
end

%% look at gain quality as a function of step vs ramp
step_idx = false(length(all_data),1);
ramp_idx = false(length(all_data),1);

for i= 1:length(all_data)
    try g = all_data(i).ft.gain; catch g = all_data(i).gain.vr; end
    step_idx(i) = ((sum(g < .85) / length(g)) > .4) & ((sum(g > 1) / length(g)) > .4); %during ramps, the time spent below .85 gain should be less than 40% of the trial, and time spent above .85 should be more than 40% of the trial  
    ramp_idx(i) = ((sum(g > .85) / length(g)) > .4) & ~step_idx(i);
end

figure(3); clf
for i = 1:length(all_data)
    subplot(2,1,1); cla
    hold on;
    try g = all_data(i).ft.gain; catch g = all_data(i).gain.vr; end
    plot(linspace(0,1,length(all_data(i).gain.g)),all_data(i).gain.g)
    plot(linspace(0,1,length(g)),g)
    

    subplot(2,1,2); hold on 
    scatter(step_idx(i) + 2*ramp_idx(i),var(all_data(i).gain.g(end/2:end),'omitnan') )
end
xlim([-1,2])
%% Look at casting as a function of gain

rows = ceil(sqrt(length(all_data)));
cols = ceil(length(all_data)/rows);


figure(1); clf
t = tiledlayout(rows,cols);
for i = 1:length(all_data)
    nexttile
    try
        histogram(all_data(i).ft.gain)
    catch
        histogram(all_data(i).gain.vr)
    end
end
xlabel(t,'Compass Gain')

figure(2); clf
t = tiledlayout(rows,cols);
for i = 1:length(all_data)
    nexttile
    try
        g = interp1(linspace(0,1,length(all_data(i).gain.g)),all_data(i).gain.g,linspace(0,1,length(all_data(i).ft.gain)));
        scatter(all_data(i).ft.gain,g,'.')
    catch
        g = interp1(linspace(0,1,length(all_data(i).gain.g)),all_data(i).gain.g,linspace(0,1,length(all_data(i).gain.vr)));
        scatter(all_data(i).gain.vr,g,'.')
    end
end
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