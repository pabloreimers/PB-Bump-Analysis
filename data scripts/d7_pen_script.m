base_dir = ('Z:\pablo\d7_pen\20260715\'); %uigetdir(); %
all_files = dir([base_dir,'\**\*imagingData*.mat']);
all_files = natsortfiles(all_files);

%% create masks

for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        
        imgData = squeeze(sum(regProduct,[3,4]));
        top_pct = prctile(imgData,98,'all');
        bot_pct = prctile(imgData,5,'all');

        imgData(imgData>top_pct) = top_pct;
        imgData(imgData<bot_pct) = bot_pct;        


        % figure(1); clf; imagesc(imgData); colormap(bone); axis equal tight; drawnow;
        % mask = roipoly();

        mask = imgData > prctile(imgData(:),60);
        tmp = regionprops(mask);
        tmp = sort([tmp.Area],'descend');
        if (tmp(1) / tmp(2)) > 1.5
            mask = bwareafilt(mask,1);
        else
            se = strel('line',10,0);
            mask = imdilate(bwareafilt(mask,2),se);
        end

        imagesc(mask); drawnow

        save([fileparts(all_files(i).folder),'\mask.mat'],'mask')
    end
end

%% store traces

all_data = struct();

tic
for i = 1:length(all_files)
    
    %read out which is being processed
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    
    %store time vectors
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])

    all_data(i).ft.stims = logical(ftData_DAQ.stim{1});
    all_data(i).ft.xb    = seconds(ftData_DAQ.volClock{:});
    all_data(i).ft.xf    = seconds(ftData_DAQ.trialTime{:});

    %load in the movie and store as summed across z, offset subtract, interpolate to DAQ timing
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask.mat'])
    imgData = squeeze(sum(regProduct,3));
    imgData = imgData - prctile(imgData,1,'all');
    imgData = permute(imgData,[3,1,2]);
    imgData = interp1(all_data(i).ft.xb,imgData,all_data(i).ft.xf,'linear','extrap');
    imgData = permute(imgData,[2,3,1]);

    %process
    all_data(i).im.in_stim  = mean(imgData(:,:,all_data(i).ft.stims),3); %save stills from in and out of stim light
    all_data(i).im.out_stim = mean(imgData(:,:,~all_data(i).ft.stims),3);

    imgData_2D              = reshape(imgData,[],length(all_data(i).ft.stims)); %convert to 2D, all_pix x frames
    all_data(i).im.trace    = mean(imgData_2D(mask(:),:),1) - mean(imgData_2D(~mask(:),:),1);
    
    %store meta
    all_data(i).meta = all_files(i).folder;
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end
    
%% store logical indices
intensity_ind = mod(0:(length(all_data)-1),3) + 1;
pen2_log = contains({all_data.meta}','pen2');


%% plot still images
max_pix = 1e3;
n = 256;
half = floor(n/2);

blue = [linspace(0,1,half)', linspace(0,1,half)', ones(half,1)];
red  = [ones(half,1), linspace(1,0,half)', linspace(1,0,half)'];

cmap = [blue; red];

figure(1); clf
for i = 1:length(all_data)
    subplot(length(all_data),3,3*(i-1) + 1)
    imagesc(all_data(i).im.out_stim)
    clim(gca,[0,max_pix])
    axis equal tight

    if pen2_log(i)
        ylabel({'PEN2',['Int Lvl: ',num2str(intensity_ind(i))]}, 'Rotation',0);
    else
        ylabel({'PEN1',['Int Lvl: ',num2str(intensity_ind(i))]}, 'Rotation',0);
    end
    

    subplot(length(all_data),3,3*(i-1) + 2)
    imagesc(all_data(i).im.in_stim)
    clim(gca,[0,max_pix])
    axis equal tight
    
    subplot(length(all_data),3,3*(i-1) + 3)
    imagesc(all_data(i).im.in_stim - all_data(i).im.out_stim)
    colormap(gca,cmap)
    clim(gca,[-max_pix max_pix]/2)
    axis equal tight
end

tmp = get(gcf,'Children');
for i = 1:length(tmp)
    tmp(i).XTick = [];
    tmp(i).YTick = [];
end

subplot(length(all_data),3,1)
title('Out Stim')

subplot(length(all_data),3,2)
title('In Stim')

subplot(length(all_data),3,3)
title('In - Out')

%% show traces
figure(2); clf

% extract stim starts
win_edges = [-2,5] * 60; %frames to extract sampled at 60 frames per second
win_frames = win_edges(1):win_edges(2);

pulses = {};

for i = 1:length(all_data)
    stim_starts = find(diff(all_data(i).ft.stims) > 0);
    stim_wins   = stim_starts + win_frames;
    tmp_trace   = all_data(i).im.trace - median(all_data(i).im.trace);
    pulses{i}   = tmp_trace(stim_wins);
end

% plot each intensity level
for int = 1:3
    subplot(3,3,0+int)
    hold on
    for i = find(intensity_ind == int)
        plot(all_data(i).ft.xf,all_data(i).im.trace - median(all_data(i).im.trace))
    end
    
    subplot(3,3,3+int); hold on
    plot(win_frames / 60, cell2mat(pulses(intensity_ind == int & ~pen2_log')'),'Color',[1,0,0,.1])
    plot(win_frames / 60, cell2mat(pulses(intensity_ind == int & pen2_log')'),'Color',[0,0,0,.1])
    
    subplot(3,3,6+int); hold on
    h = plot_sem(gca,win_frames / 60, cell2mat(pulses(intensity_ind == int & ~pen2_log')')); h.FaceColor = 'r'; h.EdgeColor = 'none';
    h = plot_sem(gca,win_frames / 60, cell2mat(pulses(intensity_ind == int & pen2_log')')); h.FaceColor  = 'k'; h.EdgeColor = 'none';
    plot(win_frames / 60,zeros(size(win_frames)),'k:')
end
legend('PEN1','PEN2','color','none')

tmp = get(gcf,'Children');
linkaxes(tmp,'y')
set(tmp,'Color','none')

for i = 1:3
    subplot(3,3,i)
    tmp_y = ylim;
    patch([all_data(i).ft.xf;flipud(all_data(i).ft.xf)],...
          [all_data(i).ft.stims * 10;...
           zeros(size(all_data(i).ft.xf))] + tmp_y(1),...
           'r','FaceAlpha',.5,'EdgeColor','none')
    xlim([min(all_data(i).ft.xf),max(all_data(i).ft.xf)])
end
%% functions
function h = plot_sem(ax,t,x)
t = reshape(t,1,[]);
m1 = mean(x,1,'omitnan');
s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

idx = ~isnan(m1);
m1 = m1(idx);
s1 = s1(idx);
t  = t(idx);


h = patch(ax,[t,fliplr(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end