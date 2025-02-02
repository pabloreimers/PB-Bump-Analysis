
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