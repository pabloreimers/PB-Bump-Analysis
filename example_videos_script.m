
vid_name = 'example_EPG_GRAB(DA2m).avi';
base_dir = 'Z:\pablo\lpsp_cl\to analyze\to do\20221212\20221212-1_EPG_GRAB(DA2m)_cl\';
%base_dir = 'Z:\pablo\lpsp_cl_redo\20240305\fly 1\20240305-1_epg_syt8m\';
%base_dir = 'Z:\pablo\lpsp_cl\to analyze\to do\20221208\20221208-1_LPsP_syt7f_cl\';

ft_vid = dir([base_dir,'FicTracData\fictrac-raw*avi']);
tmp = dir([base_dir,'*ficTracData_DAQ*']);
load([tmp.folder,'\',tmp.name]);
tmp = dir([base_dir,'registration*\*imagingData*']);
load([tmp.folder,'\',tmp.name]);

%%
xf = seconds(ftData_DAQ.trialTime{1});
xb = linspace(0,xf(end),size(regProduct,4));
ch1 = squeeze(sum(regProduct,3));

ch1 = smoothdata(ch1,3,'gaussian',10);

t1 = prctile(ch1(:),99.9);
b1 = prctile(ch1(:),5);
%b2 = prctile(ch2(:),5);

ch1 = 256*(ch1 - b1) / (t1 - b1);

ch1 = permute(ch1,[3,1,2]);
ch1 = interp1(xb,ch1,xf);
ch1 = permute(ch1,[2,3,1]);
%%
vidObj = VideoReader([ft_vid.folder,'\',ft_vid.name]);
vid = read(vidObj);
ft_vid_sum = squeeze(sum(vid(1:150,:,1,:),[1,2,3]));
start_idx = find(ft_vid_sum > mean(ft_vid_sum),1,'first');
vid = vid(:,:,:,start_idx:start_idx+length(xf));

%%
cue     = ftData_DAQ.cuePos{:}';
cue     = smoothdata(unwrap(cue / 192 *2 * pi - pi),1,'movmean',30);
cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;

%%
pause_time = 0;

figure(8); clf
set(gcf,'Color','none')

h1 = subplot(2,2,1); 
a1 = polarscatter(cue(1),1,200,'filled','m'); 
set(gca,'RLim',[0,1],'linewidth',5,'Color','none','ThetaColor','w','ThetaAxisUnits','Radians');

h2 = subplot(2,2,2); 
a2 = image(vid(:,:,:,1));
axis equal tight

h3 = subplot(2,1,2);
a3 = image(ch1(:,:,1));
colormap(h3,bone)
title('EPG > GRAB(DA2m), Closed Loop','Color','w')
xlabel('2022 12 12, Trial 1','Color',[.75,.75,.75])
fontsize(gca,20,'pixels')
axis equal tight

%%
writerObj = VideoWriter(vid_name);
writerObj.FrameRate = 120;
open(writerObj);

n = 5000; %length(xf)

for i = 21000:25000
a1.ThetaData = cue(i);
a2.CData = vid(:,:,:,i);
a3.CData  = ch1(:,:,i);
pause(pause_time)
drawnow
writeVideo(writerObj,getframe(gcf));
fprintf('%.2f\n',i/n)
end
close(writerObj);



%% Example EPG>syt8m

%load in everything
%base_dir = 'Z:\pablo\lpsp_rnai\20251021\fly 2\20251021-5_epg_syt8s_empty_thrnai';
base_dir = 'Z:\pablo\pizza_talks\epg example\20251028-4_epg_syt8m';
base_dir = 'Z:\pablo\lpsp_cl_redo\20240229\fly 1\20240229-1_epg_syt8m\';
base_dir = 'C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\20240229-1_epg_syt8m';
tmp = dir([base_dir,'\registration_001\imagingData*.mat']);
load([tmp.folder,'\',tmp.name])
tmp = dir([base_dir,'\**\*ficTracData_DAQ.mat']);
load([tmp.folder,'\',tmp.name])
tmp = dir([base_dir,'\FicTracData\*fictrac-raw*.avi']);
vidObj = VideoReader([tmp.folder, '\',tmp.name]);
load('C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\lpsp_cl_redo_data_20240306.mat')

%% sum across z
imgReg = squeeze(sum(regProduct,3));

%% register the image
imgReg = nan(size(imgData));
[optimizer,metric] = imregconfig("monomodal");
template = mean(imgData,3);
for i = 1:size(imgData,3)
    imgReg(:,:,i) = imregister(imgData(:,:,i),template,"translation",optimizer,metric);
end
%% rescale image
h0 = prctile(imgReg(:),20);
h1 = prctile(imgReg(:),99.9);

imgReg_scaled = rot90((imgReg - h0)/(h1-h0)*256,2);

%% load in movie and crop to just the imaging part
vid = read(vidObj);
vid = squeeze(vid(:,:,1,:));
ft_vid_sum = squeeze(sum(vid(1:150,:,:),[1,2]));
start_idx = find(ft_vid_sum > mean(ft_vid_sum),1,'first');
end_idx   = find(ft_vid_sum > mean(ft_vid_sum),1,'last');
vid = vid(:,:,start_idx:end_idx);

%% interpret movie and image to be the same number of frames
xf = seconds(ftData_DAQ.trialTime{1});
%xb = seconds(ftData_DAQ.volClock{1});
xb = linspace(xf(1),xf(end),size(imgReg_scaled,3));
xm = linspace(xb(1),xb(end),size(vid,3));

imgReg_scaled = permute(imgReg_scaled,[3,1,2]);
imgReg_scaled = interp1(xb,imgReg_scaled,xm);
imgReg_scaled = permute(imgReg_scaled,[2,3,1]);

cue     = ftData_DAQ.cuePos{:}';
cue     = interp1(xf,unwrap(cue),xm);
cue     = smoothdata(unwrap(cue / 192 *2 * pi - pi),1,'movmean',30);
cue     = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
cue(cue > pi) = cue(cue > pi) - 2*pi;

%% 

vid_name = 'LPsP_syt8m_cl_color';
figure(1); clf; set(gcf,'Color','none')
cmap1 = [linspace(0,.85,255)',linspace(0,.4,255)',linspace(0,.8,255)'];
cmap2 = [linspace(0,1,255)',linspace(0,1,255)',linspace(0,1,255)'];

subplot(3,1,1)
h1 = image(imgReg_scaled(:,:,1));
colormap(gca,cmap1)
axis equal tight; xticks([]); yticks([])
subplot(3,1,2:3)
h2 = imagesc(vid(:,:,1));
colormap(gca,cmap2)
axis equal tight; xticks([]); yticks([])
pos = get(gca,'Position');
pos = [pos(1)-.1,pos(2)-.1,pos(3)+.2,pos(4)+.2];
polaraxes('Position',pos,'Color','none','rlim',[0,1],'RTick',[],'ThetaTick',[]); hold on
%h3 = polarscatter(cue(i),1,500,'filled','sc');
%pos = [.25,pos(2)+.1,.5,pos(4)];
% axes('Position',pos,'Color','none','XColor','none','YColor','none'); hold on; xlim([-1,1]); ylim([-1,1])
% h3 = scatter(cos(cue(i)),sin(cue(i)),500,'filled','sc');


% writerObj = VideoWriter(vid_name);
% writerObj.FrameRate = 120;
% open(writerObj);

for i = 1.7e4:2e4 %size(imgReg_scaled,3)-10
    h1.CData = mean(imgReg_scaled(:,:,i:i+10),3);
    h2.CData = vid(:,:,i+10);
    xlabel(subplot(3,1,2:3),sprintf('t: %i',round(xm(i))),'Color','w')
    %h3.ThetaData = cue(i);
    % h3.XData = cos(cue(i+10));
    % h3.YData = sin(cue(i+10));
    drawnow
    %pause(1e-2)
    % writeVideo(writerObj,getframe(gcf));
end
% close(writerObj);
