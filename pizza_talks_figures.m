%% compartmentalization figures
load('C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\lpsp_cl_data_20240206 (1).mat'); %load in the all_data
load('C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\20221212-1_EPG_GRAB(DA2m)_cl_1\imagingData_reg_ch1_trial_001.mat')

keep_idx = false(length(all_data));
for i = 1:length(all_data)
    keep_idx(i) = ~isempty(all_data(i).meta); % check which trials to keep, some are empty
end
all_data = all_data(keep_idx);
%%
for i = 1:length(all_data)
    all_data(i).ft.xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),length(all_data(i).im.mu)); %add in time vector for imaging
end

i = 79; %20221212-1_EPG_GRAB(DA2m)_cl_1'
cmap =  [linspace(0,0,255)',linspace(0,1,255)',linspace(0,.5,255)'];
t_min = 0;
t_max = 100;
% i = 51; %20221212-1_EPG_GRAB(DA2m)_cl_1'
% cmap =  [linspace(0,.84,255)',linspace(0,.4,255)',linspace(0,.8,255)'];
% t_min = 200;
% t_max = 400;
figure(1); clf
subplot(4,1,1)
h(1) = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'Color','w','linewidth',2); h(1).YData(abs(diff(h(1).YData))>pi) = nan;
xticks(t_min:10:t_max); xticklabels([])
yticks([-pi,pi]); yticklabels({'-\pi','\pi'}); ylim([-pi,pi]); set(gca,'YDir','reverse')
ylabel('heading','Rotation',0)

subplot(4,1,2:3)
h(2) = imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z);
xticks([]); yticks([])
pos = get(gca,'Position'); a =colorbar; set(gca,'Position',pos,'Units','normalized','Colormap',cmap,'CLim',clim+[0.2,0]*range(clim),'ycolor','g')
ylabel({'EPG','GRAB(DA2m)','\DeltaF/F','(z-score)'},'Rotation',0)

ax(1)= axes('Position',[pos(1)-.105,pos(2),.1,pos(4)]);
h(4) = imagesc(rot90(mean(imgData,3),3));
xticks([]); yticks([])
set(gca,'Colormap',cmap,'CLim',clim+[0.2,0]*range(clim))

subplot(4,1,4)
h(3) = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color','w','linewidth',2); h(3).YData(abs(diff(h(3).YData))>pi) = nan;
yticks([-pi,pi]); yticklabels({'-\pi','\pi'}); ylim([-pi,pi]); set(gca,'YDir','reverse')
ylabel('bump position','Rotation',0)
xticks(t_min:10:t_max); xticklabels([])
xlabel('time (s)')

hax = get(gcf,'Children');
linkaxes(hax([1,4,5]),'x')
xlim(hax(1),[t_min,t_max])
set(get(gcf,'Children'),'color','none','xcolor','w','ycolor','w')
set(gcf,'Color','none','InvertHardCopy','off')
fontsize(gcf,20,'pixels')

%% compartmentalization movie

writerObj = VideoWriter('tmp.avi');
writerObj.FrameRate = 60;
open(writerObj);
t_vec = t_max:-.1:(t_min+2);
pos_vec = linspace((pos(1)+pos(3)),pos(1),length(t_vec));
for j = 1:length(t_vec)
    t = t_vec(j);
    h(1).YData(h(1).XData > t) = nan;
    h(2).CData(:,h(2).XData > t) = nan;
    h(3).YData(h(3).XData > t) = nan;  
    [~,ind] = min(abs(h(2).XData - t));
    h(4).CData = rot90(mean(imgData(:,:,ind-10:ind),3),3);
    %ax(1).Position(1) = pos_vec(j);
    drawnow
    writeVideo(writerObj,getframe(gcf));
end
 close(writerObj);

 v = VideoReader('tmp.avi');
 num_frames = v.NumFrames;
 writerObj = VideoWriter('EPG_GRAB_heatmap.avi');
writerObj.FrameRate = 60;
open(writerObj);

 for i = num_frames:-1:1
    frame = read(v, i); 
    
    % Write the frame to the new video file
    writeVideo(writerObj, frame);
 end
 close(writerObj)

 %% moving bump video
 figure(1); clf
 i = 79; %20221212-1_EPG_GRAB(DA2m)_cl_1'
 cmap =  [linspace(0,0,255)',linspace(0,1,255)',linspace(0,.5,255)'];
t_min = 0;
t_max = 100;
h = imagesc(rot90(mean(imgData,3),2));
axis equal tight
xticks([]); yticks([])
set(gca,'Colormap',cmap,'CLim',clim+[0.2,-.2]*range(clim))

t_ind = find(all_data(i).ft.xb > t_min & all_data(i).ft.xb < t_max);

 writerObj = VideoWriter('EPG_GRAB_vid.avi');
writerObj.FrameRate = 60;
open(writerObj);
for t = t_ind
    h.CData = rot90(mean(imgData(:,:,t:t+10),3),2);
    drawnow
    writeVideo(writerObj,getframe(gcf));
end
 close(writerObj)

%% atp figures

% load in images
files= {'C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\20240806-15_epg_syt7f_th_p2x2\registration\imagingData_trial001.mat',...
        'C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\20240806-16_epg_syt7f_th_p2x2\registration\imagingData_trial001.mat'};
img_fig = {};

for i = 1:2
load(files{i})

a = max(squeeze(sum(img{2},3)),[],3);
a = 3*a/max(a,[],'all');
a = a .* reshape([1,0,0],1,1,[]);

d = max(squeeze(sum(img{1},3)),[],3);
d = 2*d/max(d,[],'all');
d = d .* reshape([0,.7,.7],1,1,[]);
img_fig{i} = a+d;
end

%% plot everything

cmap = [linspace(0,0,255)',linspace(0,0.7,255)',linspace(0,.7,255)'];
c    = [0,.7,.7;...
        .3,.3,1];


ind = [19,29,20,30];

figure(1); clf
for j = 1:4
    clear h
i = ind(j);
h(1) = subplot(2,8,(j-1)*4 + [2:4]);

imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z);
hold on
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color',c(ceil(j/2),:),'linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'Color',[.5,.5,.5],'linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
pos = get(gca,'Position');
a = colorbar; a.Color = 'w'; ylabel(a,{'\DeltaF/F','z-scored'},'Rotation',0)
set(gca,'Colormap',cmap,'CLim',clim+[0.2,0]*range(clim),'Position',pos)
yticks([]); xticklabels([]); xticks([160:5:180])


pos = [pos(1),pos(2)+pos(4),pos(3),.05];
h(2) = axes('Units','normalized','Position',pos);

plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1,'omitnan'),'r','linewidth',2)
xticks([]); yticks([])

% h(3) = subplot(4,1,4);
% imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).atp.z);

linkaxes(h,'x')
xlim([160,180])
set(h(1),'xcolor','w','ycolor','w')
set(h(2),'xcolor','none','ycolor','none','Color','none')

end

subplot(2,8,1)
image(rot90(img_fig{1},3))
xticks([]);yticks([])
axis equal tight

subplot(2,8,9)
image(rot90(img_fig{2},3))
xticks([]);yticks([])
axis equal tight


set(gcf,'Color','none','InvertHardCopy','off')
fontsize(gcf,20,'pixels')
%% ejection video
 figure(1); clf
 i = 19; %20240806-16_epg_syt7f_th_p2x2
 cmap =  [linspace(0,0,255)',linspace(0,.7,255)',linspace(0,.7,255)'];


a = squeeze(sum(img{2},3));
a = 3*a/max(a,[],'all');
a = a .* reshape([1,0,0],1,1,1,[]);

d = squeeze(sum(img{1},3));
d = 5*d/max(d,[],'all');
d = d .* reshape([0,.7,.7],1,1,1,[]);
imgData = a+d;
imgData = permute(imgData,[1,2,4,3]);

%%
h = image(rot90(mean(imgData,4),2));
axis equal tight
xticks([]); yticks([])

t_min = 20;
t_max = 80;
t_ind = find(all_data(i).ft.xb > t_min & all_data(i).ft.xb < t_max);

writerObj = VideoWriter('ATP_P2X2_example.avi');
writerObj.FrameRate = 20;
open(writerObj);
for t = t_ind'
    h.CData = rot90(mean(imgData(:,:,:,t:t+10),4),2);
    drawnow
    writeVideo(writerObj,getframe(gcf));
end
close(writerObj)


%% rnai figures


i = 8; %lpsp > vglut-rnai
cmap = [linspace(0,1,255)',linspace(0.0,0,255)',linspace(0,1,255)'];
t_min = 20;
t_max = 150;
% 
% i = 245; % lpsp > th-rnai
% cmap = [linspace(0,0,255)',linspace(0.0,0.7,255)',linspace(0,0.7,255)'];
% t_min = 300;
% t_max = 500;

% i = 189; empty > th-rnai
% cmap = [linspace(0,1,255)',linspace(0,0.5,255)',linspace(0,0,255)'];
% t_min = 250;
% t_max = 450;
binedges = 0:.05:5;
dark_mode = true;

figure(1); clf
a1 = subplot(3,1,1);
h = imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
%if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color',cmap(end,:),'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'color',[.5,.5,.5],'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;

title(all_data(i).meta,'color','w')
xticks([t_min:20:t_max]); xticklabels({});
set(gca,'Colormap',cmap)
pos= get(gca,'Position'); colorbar; set(gca,'Position',pos,'CLim',clim+[0.3,-.1]*range(clim))

a2 = subplot(6,1,3); hold on
offset = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
a=plot(all_data(i).ft.xf,offset,'color',[.5,.5,.5],'linewidth',2); a.YData(abs(diff(a.YData))>pi) =nan;
%plot(all_data(i).ft.xf,all_data(i).ft.f_speed)
patch(all_data(i).ft.xf,2*pi*(all_data(i).ft.stims/10)-pi,'r','FaceAlpha',.1,'EdgeColor','none')
ylabel('offset','Rotation',0)
a2.YTick = [-pi,0,pi]; a2.YTickLabels = {'-\pi','0','\pi'}; a2.YLim = [-pi,pi];

xticks([t_min:20:t_max]); xticklabels({});
pos = get(gca,'Position');
pos = [pos(1)+pos(3)+.01,pos(2),.05,pos(4)];
ax = axes('Position',pos,'Color','none','XAxisLocation','top');
histogram(offset,-pi:.1:pi,'Orientation','horizontal','edgeColor','none','FaceColor',[.5,.5,.5],'Normalization','probability')
box(ax,'off')
ax.YAxisLocation =  'right'; ax.YLim = [-pi,pi]; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};
linkaxes([a2,ax],'y')

a3 = subplot(6,1,4); hold on
plot(all_data(i).gain.xt,all_data(i).gain.g,'linewidth',2)
ylabel('gain','Rotation',0);
linkaxes([a1,a2,a3],'x')
xlim([t_min,t_max])
%xlim([min(all_data(i).ft.xb),max(all_data(i).ft.xb)])
ylim([0,5])
plot(xlim,[.8,.8],'w:'); %plot(xlim,[1.6,1.6],':k')

xticks([t_min:20:t_max]); xticklabels({});

subplot(3,2,5); hold on
h = histogram(all_data(i).gain.g,'BinEdges',binedges,'FaceAlpha',.8,'Normalization','probability','EdgeColor','none');

xlabel('gain')
ylabel('counts')

subplot(3,2,6); 
histogram(all_data(i).ft.f_speed,'edgecolor','none')
xlabel('f speed')

set(get(gcf,'Children'),'color','none','xcolor','w','ycolor','w')
set(gcf,'Color','none','InvertHardCopy','off')
fontsize(gcf,20,'pixels')

%% kir figures

i = 45;
cmap = [linspace(0,0,255)',linspace(0.0,0.7,255)',linspace(0,0.7,255)'];
t_min = 0;
t_max = 600;

% i = 189;
% cmap = [linspace(0,1,255)',linspace(0,0.5,255)',linspace(0,0,255)'];
% t_min = 250;
% t_max = 450;
binedges = 0:.05:5;
dark_mode = true;

figure(1); clf
a1 = subplot(3,1,1);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
%if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color',cmap(end,:),'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'color',[.5,.5,.5],'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;

title(all_data(i).meta,'color','w')
xticks([t_min:20:t_max]); xticklabels({});
set(gca,'Colormap',cmap)
pos= get(gca,'Position'); colorbar; set(gca,'Position',pos)

a2 = subplot(6,1,3); hold on
offset = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
a=plot(all_data(i).ft.xf,offset,'color',[.5,.5,.5],'linewidth',2); a.YData(abs(diff(a.YData))>pi) =nan;
ylabel('offset','Rotation',0)
a2.YTick = [-pi,0,pi]; a2.YTickLabels = {'-\pi','0','\pi'}; a2.YLim = [-pi,pi];

xticks([t_min:20:t_max]); xticklabels({});
pos = get(gca,'Position');
pos = [pos(1)+pos(3)+.01,pos(2),.05,pos(4)];
ax = axes('Position',pos,'Color','none','XAxisLocation','top');
histogram(offset,-pi:.1:pi,'Orientation','horizontal','edgeColor','none','FaceColor',[.5,.5,.5],'Normalization','probability')
box(ax,'off')
ax.YAxisLocation =  'right'; ax.YLim = [-pi,pi]; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};
linkaxes([a2,ax],'y')

% a3 = subplot(6,1,4); hold on
% plot(all_data(i).gain.xt,all_data(i).gain.g,'linewidth',2)
% ylabel('gain','Rotation',0);
% linkaxes([a1,a2,a3],'x')
% xlim([t_min,t_max])
% %xlim([min(all_data(i).ft.xb),max(all_data(i).ft.xb)])
% ylim([0,5])
% plot(xlim,[.8,.8],'w:'); %plot(xlim,[1.6,1.6],':k')
% 
% xticks([t_min:20:t_max]); xticklabels({});
% 
% subplot(3,2,5); hold on
% h = histogram(all_data(i).gain.g,'BinEdges',binedges,'FaceAlpha',.8,'Normalization','probability','EdgeColor','none');
% 
% xlabel('gain')
% ylabel('counts')
% 
% subplot(3,2,6); 
% histogram(all_data(i).ft.f_speed,'edgecolor','none')
% xlabel('f speed')

linkaxes([a1,a2])
set(get(gcf,'Children'),'color','none','xcolor','w','ycolor','w')
set(gcf,'Color','none','InvertHardCopy','off')
fontsize(gcf,20,'pixels')