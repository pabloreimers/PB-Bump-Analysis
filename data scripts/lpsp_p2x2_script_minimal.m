
%% clear all
clear all
close all

%% find path to all relevant files
base_dir = uigetdir(); %('Z:\pablo\lpsp_p2x2\todo\');
all_files = dir([base_dir,'\**\*imagingData*.mat']);

%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    if ~isfile([fileparts(all_files(i).folder),'\mask.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        


ch1 = squeeze(mean(img{1},3));
ch2 = squeeze(mean(img{2},3));

t1 = prctile(ch1(:),90);
%t2 = prctile(ch2(:),90);
b1 = prctile(ch1(:),5);
%b2 = prctile(ch2(:),5);
figure(1)
%image(256*(mean(ch1,3)-b1)/(t1-b1) - 256*(mean(ch2,3)-b2)/(t2-b2)); axis equal tight


        %imgData = 256*(mean(ch1,3)-b1)/(t1-b1) - 256*(mean(ch2,3)-b2)/(t2-b2);
        
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
for i = 1:length(all_files)
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    %load([all_files(i).folder,'\',all_files(i).name])
    %load([fileparts(all_files(i).folder),'\mask.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    load([tmp2.folder,'\',tmp2.name])

    %regProduct = img{1};

    all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
    %all_data(i).im = process_im(img{1}, im_win, im_type, mask, n_centroid, f0_pct);
    %all_data(i).meta = all_files(i).folder;

    if ~ismember('xb',fieldnames(all_data(i).ft))
        xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),size(all_data(i).im.d,2));
    end

    %all_data(i).atp = process_im(img{2}, im_win, im_type, mask, n_centroid, f0_pct);

    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end




%% plot heading traces
idx = 1:6;

figure(2); clf
c1 = [zeros(256,1),linspace(0,1,256)',zeros(256,1)];
c2 = c1(:,[2,1,3]);

for i = idx
    a2 = subplot(length(idx),1,i);
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).atp.d,'AlphaData',1);
    colormap(a2,c2)
    yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})

    a1 = axes('Position',get(gca,  'Position')); 
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z,'AlphaData',1);
    colormap(a1,c1)
    yticks([-pi,0,pi]); yticklabels({'-\pi','0','\pi'})
    xticks([])

    set(gca,'color','none')
    hold on
    a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'m');
    a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w');
    a.YData(abs(diff(a.YData))>pi) = nan;
    xticks([]);yticks([])

    pos = get(gca,'Position');
    a3 = axes('Position',[pos(1),pos(2)+pos(4),pos(3),.03],'color','none'); 
    plot(all_data(i).ft.xb,sum(all_data(i).atp.f,1)/max(sum(all_data(i).atp.f,1)))   
    hold on
    plot(all_data(i).ft.xf,abs(all_data(i).ft.r_speed)/max(abs(all_data(i).ft.r_speed)))
    plot(all_data(i).ft.xf,abs(all_data(i).ft.f_speed)/max(abs(all_data(i).ft.f_speed)))
    yticks([]);xticks([])
    axis tight
    set(gca,'Color','none')
    linkaxes([a1,a2,a3],'x')
    linkaxes([a1,a2],'y')

    if i == length(idx)
        pos = get(gca,'Position');
        legend('atp','r speed','f speed','Location','NortheastOutside')
        set(gca,'Position',pos)
        
        pos = get(a1,'Position');
        legend(a1,'heading','pva','Location','SoutheastOutside')
        set(a1,'Position',pos)
    end

end


%% extract atp pulses

% %% show gain as a function of atp
% vel_thresh = .2;
% bump_thresh = 10;
% rho_thresh = .1;
% lag = 8;
% win_size = 5; %in seconds
% 
% figure(3); clf
% for i = 1:length(all_data)
% 
% 
% 
%     tmp_r   = all_data(i).ft.r_speed(1:end-lag);
%     tmp_mu  = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
%     tmp_mu  = tmp_mu(lag+1:end);
%     fr      = mean(diff(all_data(i).ft.xf));
%     win     = ceil(win_size/fr); %in frames
%     n = length(tmp_r) - win;
%     tmp1 = nan(n,1);
%     tmp2 = nan(n,1);
% 
%     for j = 1:n
%         [tmp1(j),tmp2(j)] = find_gain(tmp_r(j:j+win),tmp_mu(j:j+win),fr);
%     end
%     gains{i} = tmp1;
%     bias{i} = tmp2;
%     subplot(length(all_data),1,i); hold on
%     plot(tmp1)
%     plot(tmp_mu)
% end

%% show bump amp movie
figure(4); clf
n = 4;

subplot(4,1,1); hold on
h1 = scatter(0,0,'m','filled');
h2 = scatter(0,0,'k','filled');
xlim([-pi,pi]);

subplot(4,1,[2:4]); hold on
h3 = plot(all_data(n).atp.d(:,1),'r');
h4 = plot(all_data(n).im.d(:,1),'g');
legend('atp','dff')
ylim([0,2])

cue     = -interp1(all_data(n).ft.xf,all_data(n).ft.cue,all_data(n).ft.xb);

for i = 700:length(all_data(n).ft.xb)
    h1.XData = cue(i);
    h2.XData = all_data(n).im.mu(i);
    h3.YData = all_data(n).atp.d(:,i);
    h4.YData = all_data(n).im.d(:,i);
     
    pause(1e-2)
end

%% set up channels for movie (timeconsuming)
fly_num = 53;

ch1 = squeeze(mean(img{1},3));
ch2 = squeeze(mean(img{2},3));

ch1 = smoothdata(ch1,3,'movmean',5);
ch2 = smoothdata(ch2,3,'movmean',5);

t1 = prctile(ch1(:),99.9);
t2 = prctile(ch2(:),99.9);
b1 = prctile(ch1(:),5);
b2 = prctile(ch2(:),5);

% ch1 = smoothdata((ch1 - b1) / (t1 - b1),3,'movmean',5);
% ch2 = smoothdata((ch2 - b2) / (t2 - b2),3,'movmean',5);
ch1 = (ch1 - b1) / (t1 - b1);
ch2 = (ch2 - b2) / (t2 - b2);

ch1 = permute(ch1,[3,1,2]);
ch2 = permute(ch2,[3,1,2]);
ch1 = interp1(all_data(fly_num).ft.xb,ch1,all_data(fly_num).ft.xf);
ch2 = interp1(all_data(fly_num).ft.xb,ch2,all_data(fly_num).ft.xf);
ch1 = permute(ch1,[2,3,1]);
ch2 = permute(ch2,[2,3,1]);

vid = read(vidObj);
ft_vid_sum = squeeze(sum(vid,[1,2,3]));
start_idx = find(ft_vid_sum > mean(ft_vid_sum),1,'first');
end_idx = find(ft_vid_sum > mean(ft_vid_sum),1,'last');
vid = vid(:,:,:,start_idx:start_idx+length(all_data(1).ft.xf));

%xv = interp1(1:length(all_data(1).ft.xb),all_data(1).ft.xb,1:size(vid,4));

%%
pause_time = 0;
mu_tmp = interp1(all_data(fly_num).ft.xb,all_data(fly_num).im.mu,all_data(fly_num).ft.xf);

figure(8); clf
set(gcf,'Color','none')

cr = [linspace(0,1,256)',zeros(256,2)];
cg = cr(:,[2,1,3]);

h1 = subplot(3,2,1); 
a1 = image(cat(3,ch2(:,:,1),ch1(:,:,1),zeros(100,256)));
axis equal tight


h2 = subplot(2,2,2); 
a2 = image(vid(:,:,:,1));
axis equal tight

h3 = subplot(3,2,3);
a3 = image(ch1(:,:,1));
colormap(h3,cg)
axis equal tight

h4 = subplot(3,2,5); 
a4 = image(ch1(:,:,1) - ch2(:,:,1));
colormap(h4,cr);
axis equal tight

h5 = subplot(2,2,4);
%a5 = polarscatter(mu_tmp(1),1,200,'filled','w'); hold on
a6 = polarscatter(all_data(fly_num).ft.cue(1),1,200,'filled','m'); 
set(gca,'RLim',[0,1],'Color','none','ThetaColor','w','ThetaAxisUnits','Radians');

writerObj = VideoWriter('example_p2x2.avi');
writerObj.FrameRate = 60;
open(writerObj);

for i = 1:length(all_data(fly_num).ft.xf)
a1.CData = cat(3,ch2(:,:,i),ch1(:,:,i),zeros(100,256));
a2.CData = vid(:,:,:,i);
a3.CData  = 256*(ch1(:,:,i));
a4.CData  = 256*(ch2(:,:,i));
%a5.ThetaData = mu_tmp(i);
a6.ThetaData = -all_data(fly_num).ft.cue(i);
pause(pause_time)
drawnow
writeVideo(writerObj,getframe(gcf));
fprintf('%.2f\n',i/length(all_data(fly_num).ft.xf))
end
close(writerObj);

%% show population effects
win_start = -10;
win_end   = 20; %in seconds
d_pulses = {};
m_pulses = {};
f_pulses = {};
r_pulses = {};
a_pulses = {};
dark_idx = {};

figure(10)
counter = 1;
for i = 1:length(all_data)
    [~,loc] = findpeaks(sum(all_data(i).atp.d,1),'MinPeakProminence',9);
    tmp_win = floor(win_start/mean(diff(all_data(i).ft.xb))):ceil(win_end/mean(diff(all_data(i).ft.xb)));
    tmp_win2= floor(win_start/mean(diff(all_data(i).ft.xf))):ceil(win_end/mean(diff(all_data(i).ft.xf)));

    for j = loc
        [~,k] = min(abs(all_data(i).ft.xf - all_data(i).ft.xb(j)));
        d_pulses{counter} = all_data(i).im.d(:,j+tmp_win);
        m_pulses{counter} = all_data(i).im.mu(j+tmp_win);
        c_pulses{counter} = all_data(i).ft.cue(k+tmp_win2);
        f_pulses{counter} = all_data(i).ft.f_speed(k+tmp_win2);
        r_pulses{counter} = all_data(i).ft.r_speed(k+tmp_win2);
        a_pulses{counter} = all_data(i).atp.d(:,j+tmp_win);
        dark_idx{counter} = contains(all_data(i).meta,'dark');
        
        counter = counter+1;
    end
end

right_idx = cellfun(@(x)(sum(x(1:size(x,1)/2,:),'all') > sum(x(size(x,1)/2:end,:),'all')),a_pulses);
dark_idx = logical(cell2mat(dark_idx));
t = all_data(1).ft.xb(1:length(tmp_win)) + win_start;
t2 = all_data(1).ft.xf(1:length(tmp_win2)) + win_start;

group_idx = right_idx + dark_idx*2;
group_labels = {'left cl','right cl','left dark','right dark'};

for k = unique(group_idx)

tmp = find(group_idx == k);
cols = 4;
rows = ceil(length(tmp)/cols);

figure(k+1); clf
sgtitle(group_labels{k+1})
for j = 1:length(tmp)
    i = tmp(j);
    subplot(rows,cols,j)
    imagesc(t,unwrap(all_data(1).im.alpha),d_pulses{i})
    [~,pos] = max(sum(a_pulses{i},2));
    y = unwrap(all_data(1).im.alpha);
   
    hold on
    scatter(0,y(pos),'*r')
    a = plot(t,m_pulses{i},'w'); a.YData(abs(diff(a.YData))>pi) = nan;
    a = plot(t2,-c_pulses{i},'m'); a.YData(abs(diff(a.YData))>pi) = nan;
end

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

%% look at LTD on non-atp sweeps

d_sweeps = {};
m_sweeps = {};
f_sweeps = {};
r_sweeps = {};
a_sweeps = {};
l_sweeps = {};

counter = 1;
for i = 1:length(all_data)
    [~,loc] = findpeaks(smoothdata(diff(unwrap(all_data(i).ft.cue)),1,'movmedian',500),'MinPeakDistance',1000,'MinPeakHeight',.015);
    tmp_win = floor(win_start/mean(diff(all_data(i).ft.xb))):ceil(win_end/mean(diff(all_data(i).ft.xb)));
    tmp_win2= floor(win_start/mean(diff(all_data(i).ft.xf))):ceil(win_end/mean(diff(all_data(i).ft.xf)));

    for k = loc'
        [~,j] = min(abs(all_data(i).ft.xf(k) - all_data(i).ft.xb));
        d_sweeps{counter} = all_data(i).im.d(:,j+tmp_win);
        m_sweeps{counter} = all_data(i).im.mu(j+tmp_win);
        c_sweeps{counter} = all_data(i).ft.cue(k+tmp_win2);
        f_sweeps{counter} = all_data(i).ft.f_speed(k+tmp_win2);
        r_sweeps{counter} = all_data(i).ft.r_speed(k+tmp_win2);
        a_sweeps{counter} = all_data(i).atp.d(:,j+tmp_win);
        l_sweeps{counter} = all_data(i).im.f(:,j+tmp_win);
        
        counter = counter+1;
    end
end

%% population average
figure(6); clf
t = all_data(1).ft.xb(1:length(tmp_win)) + win_start;
t2 = all_data(1).ft.xf(1:length(tmp_win2)) + win_start;

m = cell(length(unique(group_idx)),1);
c = cell(length(unique(group_idx)),1);
for i = unique(group_idx)
m{i+1} = cell2mat(m_pulses(group_idx==i));
m{i+1} = unwrap(m{i+1}-m{i+1}(1,:))';
c{i+1} = -cell2mat(c_pulses(group_idx==i));
c{i+1} = unwrap(c{i+1}-c{i+1}(1,:))';
end

subplot(2,2,1); hold on
title('CL')
plot(t,m{1},'r')
plot(t,m{2},'c')
axis tight
plot([0,0],ylim,'k:')
plot(xlim,[0,0],'k:')
ylabel('unwrapped pva')
set(gca,'YDir','reverse')

subplot(2,2,3); hold on
h = plot_sem(gca,t,m{1}); h.FaceColor = 'r';
h = plot_sem(gca,t,m{2}); h.FaceColor = 'c';
h = plot_sem(gca,t2,c{1}); h.FaceColor = 'm';
h = plot_sem(gca,t2,c{2}); h.FaceColor = 'g';

legend('mu (right)','mu (left)','cue (right)','cue (left)','autoupdate','off','Location','Northwest')
axis tight
plot([0,0],ylim,'k:')
plot(xlim,[0,0],'k:')
xlabel('time (s)')
ylabel('unwrapped pva')
set(gca,'YDir','reverse')

subplot(2,2,2); hold on
title('Dark')
plot(t,m{3},'r')
plot(t,m{4},'c')
axis tight
plot([0,0],ylim,'k:')
plot(xlim,[0,0],'k:')
set(gca,'YDir','reverse')

subplot(2,2,4); hold on
x = cell(length(unique(group_idx)),1);
for i = unique(group_idx)
x{i+1} = cell2mat(m_pulses(group_idx==i));
x{i+1} = unwrap(x{i+1}-x{i+1}(1,:))';
end

h = plot_sem(gca,t,x{3}); h.FaceColor = 'r';
h = plot_sem(gca,t,x{4}); h.FaceColor = 'c';
h = plot_sem(gca,t2,c{3}); h.FaceColor = 'm';
h = plot_sem(gca,t2,c{4}); h.FaceColor = 'g';

legend('mu (right)','mu (left)','cue (right)','cue (left)','autoupdate','off','Location','Northwest')
axis tight
plot([0,0],ylim,'k:')
plot(xlim,[0,0],'k:')
xlabel('time (s)')
set(gca,'YDir','reverse')

fontsize(gcf,20,'pixels')

ax = get(gcf,'Children');
if dark_mode
    set(gcf,'color','none','InvertHardcopy','off')
    
    for i = 1:length(ax)
        if contains(class(ax(i)),'Axes')
            ax(i).Title.Color = 'w';
            set(ax(i),'Color','none','ycolor','w','xcolor','w')
        else
            ax(i).Color = 'none';
            ax(i).TextColor = 'w';
        end
    end
end


%% look at locomotor effects
dark_mode = true;

figure(10); clf
for i = unique(group_idx)
    subplot(2,2,i+1)
    hold on
    plot(t2,cell2mat(f_pulses(group_idx==i)),'Color',[[0,0,0]+dark_mode,.3])
    plot_sem(gca,t2,cell2mat(f_pulses(group_idx==i))')
    axis tight
    title(group_labels(i+1))
    plot([0,0],ylim,':','Color',[0,0,0]+dark_mode)
    ylabel('mm/s'); ylim([-3,10])
end
sgtitle('F Speed')
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

%% assess changes to dff aftr pulse
dark_mode = true;
% cB = [linspace(1,0,256)',linspace(1,0,256)',ones(256,1)];
% cR = cB(:,[3,1,2]);
cG = [zeros(256,1),linspace(0,1,256)',zeros(256,1)];
cR = cG(:,[2,1,3]);

group_idx = right_idx + dark_idx*2;
group_labels = {'left cl','right cl','left dark','right dark'};

figure(11); clf; hold on
for j = unique(group_idx)
    tmpD = nan([size(d_pulses{1}),sum(group_idx==j)]);
    tmpA = nan([size(a_pulses{1}),sum(group_idx==j)]);
    subplot(2,2,j+1)
    counter = 1;
for i = find(group_idx==j)
[~,k] = max(sum(a_pulses{i},2));
tmp1 = circshift(d_pulses{i},-k+ceil(.5*length(all_data(1).im.alpha)),1);
tmp2 = circshift(a_pulses{i},-k+ceil(.5*length(all_data(1).im.alpha)),1);
idx = t<10;
tmpD(:,:,counter) = tmp1;
tmpA(:,:,counter) = tmp2;

counter = counter+1;
end


im1 = mat2rgb(mean(tmpD,3),cG);
im2 = mat2rgb(mean(tmpA,3),cR);

image((im1+im2)/2)

image(t,[],cat(3,1.5*mean(tmpA,3)./max(mean(tmpA,3),[],'all'),...
                     mean(tmpD,3)./max(mean(tmpD,3),[],'all'),...
                     zeros(size(mean(tmpD,3)))))
text(max(xlim),max(ylim),'ATP','Color','r','HorizontalAlignment','right','VerticalAlignment','bottom')
text(max(xlim),max(ylim)-2,'dF/F','Color','g','HorizontalAlignment','right','VerticalAlignment','bottom')
ylabel('aligned glomeruli')
xlabel('time (s)')
title(group_labels{j+1})
end
fontsize(gcf,20,'pixels')

ax = get(gcf,'Children');
if dark_mode
    set(gcf,'color','none','InvertHardcopy','off')
    set(ax,'Color','none','ycolor','w','xcolor','w')
    for i = 1:length(ax)
        ax(i).Title.Color = 'w';
    end
end
%% assess changes to dff aftr pulse in the sweeps
dark_mode = true;
% cB = [linspace(1,0,256)',linspace(1,0,256)',ones(256,1)];
% cR = cB(:,[3,1,2]);
cG = [zeros(256,1),linspace(0,1,256)',zeros(256,1)];
cR = cG(:,[2,1,3]);

group_idx = right_idx + dark_idx*2;
group_labels = {'left cl','right cl','left dark','right dark'};

figure(12); clf; hold on
tmpL = nan([size(l_sweeps{1}),length(l_sweeps)]);
tmpL = nan([size(d_sweeps{1}),length(d_sweeps)]);
tmpA = nan([size(a_sweeps{1}),length(a_sweeps)]);
pulse_idx = false(length(d_sweeps),1);

counter = 1;
for i = 1:length(f_sweeps)
    if any(a_sweeps{i}>1,'all')
        [~,k] = max(sum(a_sweeps{i},2));
        pulse_idx(i) = true;
    end
tmp1 = circshift(d_sweeps{i},-k+ceil(.5*length(all_data(1).im.alpha)),1);
tmp2 = circshift(a_sweeps{i},-k+ceil(.5*length(all_data(1).im.alpha)),1);
tmp3 = circshift(l_sweeps{i},-k+ceil(.5*length(all_data(1).im.alpha)),1);


tmpD(:,:,counter) = tmp1;
tmpA(:,:,counter) = tmp2;
tmpL(:,:,counter) = tmp3;

counter = counter+1;
end


% im1 = mat2rgb(mean(tmpL,3),cG);
% im2 = mat2rgb(mean(tmpA,3),cR);
% 
% image((im1+im2)/2)

image(t,[],cat(3,.5*mean(tmpA(:,:,~pulse_idx),3)./max(mean(tmpA(:,:,~pulse_idx),3),[],'all'),...
                     mean(tmpD(:,:,~pulse_idx),3)./max(mean(tmpD(:,:,~pulse_idx),3),[],'all'),...
                     zeros(size(mean(tmpL,3)))))
text(max(xlim),max(ylim),'ATP','Color','r','HorizontalAlignment','right','VerticalAlignment','bottom')
text(max(xlim),max(ylim)-2,'Raw F','Color','g','HorizontalAlignment','right','VerticalAlignment','bottom')
ylabel('aligned glomeruli')
xlabel('time (s)')
title(group_labels{j+1})
axis tight

fontsize(gcf,20,'pixels')

ax = get(gcf,'Children');
if dark_mode
    set(gcf,'color','none','InvertHardcopy','off')
    set(ax,'Color','none','ycolor','w','xcolor','w')
    for i = 1:length(ax)
        ax(i).Title.Color = 'w';
    end
end

%% compare max fluorescence per glomerulus before and after pulse
max_first = nan(length(all_data(i).im.alpha),length(all_data));
max_last  = nan(length(all_data(i).im.alpha),length(all_data));
atp_check  = nan(length(all_data(i).im.alpha),length(all_data));

for i = 1:length(all_data)
    [~,k] = max(sum(all_data(i).atp.d,2));
    tmp_idx = max(all_data(i).atp.d,[],1)-smoothdata(max(all_data(i).atp.d,[],1),2,'movmean',1000) < .5;
    first_idx = 1:length(all_data(i).im.d) < length(all_data(i).im.d)/2;

    max_first(:,i) = circshift(max(all_data(i).im.d(:,tmp_idx&first_idx),[],2),-k+ceil(.5*length(all_data(i).im.alpha)));
    max_last(:,i) = circshift(max(all_data(i).im.d(:,tmp_idx&~first_idx),[],2),-k+ceil(.5*length(all_data(i).im.alpha)));
    atp_check(:,i) = circshift(max(all_data(i).atp.d,[],2),-k+ceil(.5*length(all_data(i).im.alpha)));
end

plot_sem(gca,(1:32)',max_first'-max_last')
%% make movie of pulses
pause_time = 1e-1;
figure(12); clf; hold on

k = 65;
h1 = plot(unwrap(all_data(1).im.alpha),d_pulses{1}(:,1),'g');
h2 = plot(unwrap(all_data(1).im.alpha),a_pulses{1}(:,1),'r');
h3 = scatter(0,0,'m','filled');
ylim([-.1,2])

tmp_c = interp1(t2,c_pulses{k},t);
for i = 1:length(d_pulses{1})
    h1.YData = d_pulses{k}(:,i);
    h2.YData = a_pulses{k}(:,i);
    h3.XData = tmp_c(i);
    drawnow
    pause(pause_time)
end

%% look at gain
gains = nan(length(all_data),1);
pre_idx = false(length(all_data),1);

for i = 1:length(all_data)
    fr = mean(diff(all_data(i).ft.xf));
    mu_tmp = interp1(all_data(i).ft.xb,all_data(i).im.mu,all_data(i).ft.xf,'linear','extrap');
    rho_tmp = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf,'linear','extrap');
    atp_tmp = interp1(all_data(i).ft.xb,sum(all_data(i).atp.d,1),all_data(i).ft.xf,'linear','extrap');
    cue_tmp = all_data(i).ft.cue;
    m_speed = gradient(mu_tmp)/fr;

    [~,loc] = findpeaks(atp_tmp,'MinPeakProminence',9);
    idx = true(size(m_speed));
    idx(loc) = false; idx = smoothdata(idx,'movmean',ceil([20,3]/fr)); idx(idx<1) = false;
    idx(abs(all_data(i).ft.r_speed) < .1) = false;
    idx(abs(m_speed)>10) = false;
    idx = logical(idx);


    pre_idx(i) = isempty(loc);
    gains(i) = m_speed(idx) \ all_data(i).ft.r_speed(idx);
end
    
%% plot all
vel_thresh = .2;
bump_thresh = 10;
rho_thresh = .1;
lag = 8;

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

    tmp2 = regexp(tmp{end-1},'\d*','Match');
    dark_idx(i) = mod(str2double(tmp2{2}),2) == 0;
    fly_id{i} = [tmp{4},'_',tmp{5}];
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
    

for i = 1:length(all_data)
    figure(ind(i)); nexttile; hold on
    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = all_data(i).ft.r_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);

    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;
    scatter(fly_vel(idx),bump_vel(idx),'w','filled','markerfacealpha',.1)
    axis equal
    y = ylim; x = xlim;
    plot(x,[0,0],':w'); 
    plot([0,0],y,':w');
    b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); %fit the slope of fly vel and bump vel with an arbitrary offset
    r = corr(fly_vel(idx),bump_vel(idx));
    plot([0,0],y,'w:')
    plot(x,[0,0],'w:')
    plot(x,x*b(2) + b(1),'r')
    text(x(2),y(1),sprintf('gain: %.2f\nr: %.2f',b(2),r),'HorizontalAlignment','right','VerticalAlignment','bottom','color','w')
    xlim(x); ylim(y);

    cc(i) = r;
    gains(i) = b(2);

    set(gca,'xcolor','w','ycolor','w','color','none')
end

[~,~,fly_num] = unique(fly_id);

for i = unique(ind)'
    fontsize(figure(i),20,'pixels')
    title(t{i},group_order(i),'fontsize',40,'color','w')
    xlabel(t{i},'fly vel (rad/s)','fontsize',30,'color','w'); ylabel(t{i},'bump vel (rad/s)','fontsize',30,'color','w')
    set(gcf,'color','none','InvertHardcopy','off')
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
    %xb = linspace(min(all_data(i).ft.xf),max(all_data(i).ft.xf),length(all_data(i).im.mu));
    xb = all_data(i).ft.xb;

    figure(ind(i)); nexttile; hold on
    imagesc(xb,unwrap(all_data(i).im.alpha),all_data(i).im.z); colormap(parula)
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
group_order = {'empty (CL)','LPsP (CL)','empty (dark)','LPsP (dark)'};
ind = 2*dark_idx + lpsp_idx + 1;
figure(5); clf

subplot(2,2,1)
scatter(ind, cc,'w'); hold on
m = accumarray(ind,cc,[],@mean);
s = accumarray(ind,cc,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Correlation Coefficent', '(fly vs bump vel)'})
set(gca,'color','none','ycolor','w','xcolor','w')

subplot(2,2,2)
scatter(ind, gains,'w'); hold on
m = accumarray(ind,gains,[],@mean);
s = accumarray(ind,gains,[],@(x)(std(x)/sqrt(length(x))));
errorbar(unique(ind)+.1,m,s,'ro')
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Gain', '(fly vs bump vel)'})
set(gca,'color','none','ycolor','w','xcolor','w')
set(gcf,'color','none','InvertHardCopy','Off')
%% test significance (bootstrap to ask about a mean difference, we don't know how variances compare)
N = 1e4;
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
xticks(1:4); xticklabels(group_order); xlim([.5,4.5])
ylabel({'Correlation Coefficient (R)', 'Resampled Means'})
hold on
y = ylim;
plot([1,2],[y(2),y(2)],'w'); text(1.5,y(2),sprintf('p = %.3f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
plot([3,4],[y(2),y(2)],'w'); text(3.5,y(2),sprintf('p = %.3f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color','w')
set(gca,'color','none','ycolor','w','xcolor','w')

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
set(gcf,'color','none','InvertHardcopy','off')

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
opt_gain= nan(length(all_data),1);

for i = 1:length(all_data)
    xf = all_data(i).ft.xf;
    xb = linspace(min(xf),max(xf),size(all_data(i).im.mu,1));
    fr = mean(diff(xf));

    fly_vel  = all_data(i).ft.r_speed;
    bump_vel = gradient(interp1(xb,unwrap(all_data(i).im.mu),xf))/fr;
    rho      = interp1(xb,all_data(i).im.rho,xf);

    obj_fun = @(lag) -find_cc(lag,fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh);

    opt_lag(i) = ga(obj_fun,1,[],[],[],[],0,40,[],1);
    [opt_cc(i),opt_gain(i)]  = find_cc(opt_lag(i),fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh);
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

p_cl     = sum( lpsp_cl - empty_cl < 0) / N;
p_dark     = sum( lpsp_dark - empty_dark < 0) / N;

subplot(2,2,3);cla; hold on
a = swarmchart(ones(N,4) .* [1:4],[...
     empty_cl', ...
     lpsp_cl', ...
     empty_dark', ...
     lpsp_dark'], ...
     '.',c);
xticks(unique(ind)); xticklabels(group_order); ylabel({'optimal lag','Resampled means'})
y = ylim;
plot([1,2],[y(2),y(2)],c); text(1.5,y(2),sprintf('p = %.6f',p_cl),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)
plot([3,4],[y(2),y(2)],c); text(3.5,y(2),sprintf('p = %.6f',p_dark),'HorizontalAlignment','center','VerticalAlignment','bottom','color',c)

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

%% check fluorsence and bump width
n = length(all_data);

figure(10); clf
for i = 1:n
    %subplot(2,2,ind(i)); hold on
    scatter(ind(i),mean(all_data(i).im.f,'all')); hold on
end
% for i = unique(ind)'
%     title(subplot(2,2,i),group_order(i))
% end
linkaxes(get(gcf,'Children'),'y')


%% Functions

function s = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type)

    f_speed = ftData_dat.velFor{:};                       %store each speed
    f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}';
    cue     = smoothdata(unwrap(cue / 192 *2 * pi - pi),1,ft_type,ft_win);
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

m1 = mean(x,1);
s1 = std(x,1)./sqrt(size(x,1));

h = patch(ax,[t;flipud(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end

function rgb_image = mat2rgb(v,map)
    minv = min(v(:));
    maxv = max(v(:));
    ncol = size(map,1);
    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
    rgb_image = ind2rgb(s,map);
end