
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
base_dir = 'Z:\pablo\pizza_talks\epg example\20251028-2_epg_syt8m';
%base_dir = 'Z:\pablo\lpsp_cl_redo\20240229\fly 1\20240229-1_epg_syt8m\';
%base_dir = 'C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\20240229-1_epg_syt8m';
tmp = dir([base_dir,'\registration_001\imagingData.mat']);
load([tmp.folder,'\',tmp.name])
tmp = dir([base_dir,'\**\*ficTracData_DAQ.mat']);
load([tmp.folder,'\',tmp.name])
tmp = dir([base_dir,'\FicTracData\*fictrac-raw*.avi']);
vidObj = VideoReader([tmp.folder, '\',tmp.name]);
load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\pizza_talks_examples.mat')

%% sum across z
% imgReg = squeeze(sum(regProduct,3));
imgReg = imgData;
%% register the image
imgReg = nan(size(imgData));
[optimizer,metric] = imregconfig("monomodal");
template = mean(imgData,3);
for i = 1:size(imgData,3)
    imgReg(:,:,i) = imregister(imgData(:,:,i),template,"translation",optimizer,metric);
end
%% rescale image
h0 = prctile(imgReg(:),20);
h1 = prctile(imgReg(:),99.5);

imgReg_scaled = rot90((imgReg - h0)/(h1-h0)*256,2);

%% load in movie and crop to just the imaging part
vid = read(vidObj);
vid = squeeze(vid(:,:,1,:));
ft_vid_sum = squeeze(sum(vid(1:150,:,:),[1,2]));
start_idx = find(ft_vid_sum > mean(ft_vid_sum),1,'first');
end_idx   = find(ft_vid_sum > mean(ft_vid_sum),1,'last');
vid = vid(1:256,:,start_idx:end_idx);

%% interpret movie and image to be the same number of frames
xf = seconds(ftData_DAQ.trialTime{1});
xb = seconds(ftData_DAQ.volClock{1});
%xb = linspace(xf(1),xf(end),size(imgReg_scaled,3));
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
clear h
vid_name = 'EPG_syt8s_cl_color_big';
t_min = 20;
t_max = 80;

figure(1); clf; set(gcf,'Color','none')
cmap1 = [linspace(0,0,255)',linspace(0,1,255)',linspace(0,1,255)'];
cmap2 = [linspace(0,1,255)',linspace(0,1,255)',linspace(0,1,255)'];

subplot(3,3,1)
h1 = image(imgReg_scaled(:,:,1));
colormap(gca,cmap1)
axis equal tight; xticks([]); yticks([])
pos1 = get(gca,'Position');

subplot(3,3,4)
h2 = image(vid(1:256,:,1));
colormap(gca,cmap2)
axis equal tight; xticks([]); yticks([])


subplot(3,3,2:3)
h(1) = imagesc(all_data(2).ft.xb,unwrap(all_data(2).im.alpha),all_data(2).im.z);
pos = get(gca,'Position');
a = colorbar; set(a,'Color', 'w','ycolor','w'); ylabel(a,{'z-scored','\DeltaF/F'}, 'Rotation',0)
yticks([])
set(gca,'Colormap',cmap1,'CLim',clim+[0,-.3]*range(clim),'Position',pos,'xcolor','w','ycolor','w','Fontsize',20)
xlim([t_min,t_max]); xticklabels([])


subplot(3,3,5:6)
h(2) = plot(all_data(2).ft.xf,-all_data(2).ft.cue,'w','linewidth',2); h(2).YData(abs(diff(h(2).YData))>pi) = nan;
set(gca,'Color','none','xcolor','w','ycolor','w','YDir','reverse','Fontsize',30)
xlim([t_min,t_max]); xticklabels([])
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'})

subplot(3,3,8:9); 
h(3) = plot(all_data(2).ft.xf,-all_data(2).ft.cue,'w','linewidth',2); h(3).YData(abs(diff(h(3).YData))>pi) = nan;
hold on
h(4) = plot(all_data(2).ft.xb,all_data(2).im.mu,'Color',[0,1,1],'linewidth',2); h(4).YData(abs(diff(h(4).YData))>pi) = nan;
set(gca,'Color','none','xcolor','w','ycolor','w','YDir','reverse','Fontsize',30)
xlim([t_min,t_max]); xticklabels([])
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'})
xlabel('time (s)')

%xlabel(all_data(2).meta,'color','w')

%fontsize(gcf,40,'pixels')
%%

writerObj = VideoWriter('tmp.avi');
writerObj.FrameRate = 60;
open(writerObj);


i_vec = find(xm > t_min & xm < t_max);
for i = fliplr(i_vec) %size(imgReg_scaled,3)-10
    h1.CData = mean(imgReg_scaled(:,:,i:i+10),3);
    h2.CData = vid(1:256,:,i+10);
    
    t = xm(i);
    h(1).CData(:,h(1).XData > t) = nan;
    h(2).YData(h(2).XData > t) = nan;
    h(3).YData(h(3).XData > t) = nan;
    h(4).YData(h(4).XData > t) = nan;
    

    drawnow
    %pause(1e-2)
    writeVideo(writerObj,getframe(gcf));
end
 close(writerObj);


 v = VideoReader('tmp.avi');
 num_frames = v.NumFrames;
 writerObj = VideoWriter(vid_name);
writerObj.FrameRate = 60;
open(writerObj);

 for i = num_frames:-1:1
    frame = read(v, i); 
    
    % Write the frame to the new video file
    writeVideo(writerObj, frame);
 end
 close(writerObj)