%% analysis pipeline figures
%load('Z:\pablo\epg_gain_change\20250624\fly 4\20250624-7_epg_syt8m_gain_change_8\registration_001\imagingData_reg_ch1_trial_001.mat')
%load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\epg_gain_change_20250630.mat')

b = [linspace(1,0,255)',linspace(1,0,255)',linspace(1,0,255)'];
figure(7); clf
subplot(4,4,1)
imagesc(squeeze(sum(regProduct,[3,4])))
set(gca,'Colormap',b.^2)
xticks([]);yticks([])
axis equal tight

subplot(4,4,2)
imagesc(mask)
xticks([]);yticks([])
axis equal tight
set(gca,'Colormap',b.^2)

subplot(4,4,3) %run this section fom ImagingAnalysis_Mel
imagesc(squeeze(sum(regProduct,[3,4])))                                    %plot the image again with max intensity over time to show the whole pb
hold on
plot(x_mid,y_mid,'w')                                                   %plot the midline in white
scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
set(gca,'Colormap',b.^2)
xticks([]); yticks([])
axis equal tight

subplot(4,4,4)
imagesc(squeeze(sum(regProduct,[3,4]))) 
hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_mask(idx == i),y_mask(idx == i),5,'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end
set(gca,'Colormap',b.^2)
xticks([]); yticks([])
axis equal tight

subplot(4,1,2)
imagesc(all_data(32).ft.xb,unwrap(all_data(32).im.alpha),all_data(32).im.z)
xlim([350,550])
c = clim;
set(gca,'Colormap',b,'CLim',[min(clim)+.1*range(clim),max(clim)-.1*range(clim)])
pos = get(gca,'Position');colorbar;set(gca,'Position',pos)

subplot(4,1,3)
a = plot(all_data(32).ft.xb,all_data(32).im.mu,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([350,550])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('PVA')
set(gca,'YDir','reverse')

subplot(4,1,4)
a = plot(all_data(32).ft.xf,-all_data(32).ft.cue,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([350,550])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('Heading (VR)')
set(gca,'YDir','reverse')

%% Compartmentalization Figures
load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\lpsp_cl_data_20240206.mat')
load('Z:\pablo\lpsp_cl\20221122\20221122-1_EPG_GRABDA(2m)_1\registration_001\imagingData_reg_ch1_trial_001.mat')

i = 38;
b = [linspace(1,0,255)',linspace(1,.7,255)',linspace(1,0,255)'];
figure(7); clf
subplot(4,4,1)
imagesc(squeeze(sum(regProduct,[3,4])))
set(gca,'Colormap',b,'CLim',[min(clim)+.1*range(clim),max(clim)-.7*range(clim)])
xticks([]);yticks([])
axis equal tight

subplot(4,4,2)
imagesc(mask)
xticks([]);yticks([])
axis equal tight
set(gca,'Colormap',b.^2)

% subplot(4,4,3) %run this section fom ImagingAnalysis_Mel
% imagesc(squeeze(sum(regProduct,[3,4])))                                    %plot the image again with max intensity over time to show the whole pb
% hold on
% plot(x_mid,y_mid,'w')                                                   %plot the midline in white
% scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
% set(gca,'Colormap',b.^2)
% xticks([]); yticks([])
% axis equal tight

subplot(4,4,4)
imagesc(squeeze(sum(regProduct,[3,4]))) 
hold on
for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
    scatter(x_mask(idx == i),y_mask(idx == i),5,'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
end
set(gca,'Colormap',b.^2)
xticks([]); yticks([])
axis equal tight

h(1) = subplot(4,1,2);
all_data(i).ft.xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),length(all_data(i).im.mu));
tmp = all_data(i).im.z(:,1:end);
tmp(:,sum(tmp,1)>30) = nan;
imagesc(all_data(i).ft.xb(1:end),unwrap(all_data(i).im.alpha),tmp)
xlim([100,220])
c = clim;
set(gca,'Colormap',b,'CLim',clim+[0,0]*range(clim))
pos = get(gca,'Position');colorbar;set(gca,'Position',pos)

h(2) = subplot(4,1,3);
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([100,220])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('PVA')
set(gca,'YDir','reverse')

h(3) = subplot(4,1,4);
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([100,220])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('Heading (VR)')
set(gca,'YDir','reverse')

linkaxes(h,'x')

%%
%% Compartmentalization Figures
%load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\lpsp_cl_data_20240206.mat')
%load('Z:\pablo\lpsp_cl\20221123\20221123-7_LPsP_syt7f_cl_2\registration_001\imagingData_reg_ch1_trial_001.mat')
i = 39;
b = [linspace(1,.7,255)',linspace(1,0,255)',linspace(1,.7,255)'];
figure(7); clf
subplot(4,4,1)
imagesc(squeeze(sum(regProduct,[3,4])))
set(gca,'Colormap',b,'CLim',clim+[0,0]*range(clim))
xticks([]);yticks([])
axis equal tight

subplot(4,4,2)
imagesc(mask)
xticks([]);yticks([])
axis equal tight
set(gca,'Colormap',b.^2)

% subplot(4,4,3) %run this section fom ImagingAnalysis_Mel
% imagesc(squeeze(sum(regProduct,[3,4])))                                    %plot the image again with max intensity over time to show the whole pb
% hold on
% plot(x_mid,y_mid,'w')                                                   %plot the midline in white
% scatter(centroids(:,2),centroids(:,1),[],cmap,'filled')         %show the centroids in each of their colors
% set(gca,'Colormap',b.^2)
% xticks([]); yticks([])
% axis equal tight
% 
% subplot(4,4,4)
% imagesc(squeeze(sum(regProduct,[3,4]))) 
% hold on
% for i = 1:2*n_centroid                          %overlay each pixel in its indexed color onto the pb image. 
%     scatter(x_mask(idx == i),y_mask(idx == i),5,'filled','MarkerFaceColor',cmap(i,:),'MarkerFaceAlpha',0.2)
% end
% set(gca,'Colormap',b.^2)
% xticks([]); yticks([])
% axis equal tight

h(1) = subplot(4,1,2);
all_data(i).ft.xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),length(all_data(i).im.mu));
tmp = all_data(i).im.z(:,1:end);
tmp(:,sum(tmp,1)>100) = nan;
imagesc(all_data(i).ft.xb(1:end),unwrap(all_data(i).im.alpha),tmp)
xlim([200,500])
c = clim;
set(gca,'Colormap',b,'CLim',clim+[0,0]*range(clim))
pos = get(gca,'Position');colorbar;set(gca,'Position',pos)

h(2) = subplot(4,1,3);
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([200,500])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('PVA')
set(gca,'YDir','reverse')

h(3) = subplot(4,1,4);
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'k');
a.YData(abs(diff(a.YData))>pi) = nan;
xlim([200,500])
ylim([-pi,pi])
yticks([-pi,0,pi])
yticklabels({'-\pi',0,'\pi'})
ylabel('Heading (VR)')
set(gca,'YDir','reverse')

linkaxes(h,'x')

%% atp figures
% load('Z:\pablo\lpsp_p2x2\open loop\drive\20240806\fly 4\20240806-15_epg_syt7f_th_p2x2\registration\imagingData_trial001.mat')
% load('Z:\pablo\lpsp_p2x2\open loop\drive\20240806\fly 4\20240806-16_epg_syt7f_th_p2x2\registration\imagingData_trial001.mat')

i = 30;

figure(4); clf


a = max(squeeze(sum(img{2},3)),[],3);
a = 3*a/max(a,[],'all');
a = 1- a .* reshape([0,.7,.7],1,1,[]);

d = max(squeeze(sum(img{1},3)),[],3);
d = 2*d/max(d,[],'all');
d = 1- d .* reshape([1,.2,.2],1,1,[]);

figure(7); clf
subplot(4,4,1)
image(a+d - 1)
xticks([]);yticks([])
axis equal tight

h(1) = subplot(4,1,2);
b = [linspace(1,0,255)',linspace(1,0.7,255)',linspace(1,.7,255)'];
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z);
hold on
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'c','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'k','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)
set(gca,'Colormap',b,'CLim',clim+[0.2,0]*range(clim))

h(2) = subplot(4,1,3);
plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1),'r','linewidth',2)

h(3) = subplot(4,1,4);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).atp.z);

linkaxes(h,'x')
xlim([160,180])

%% kir figures
ind = [23,31];

figure(5); clf; hold on

for j = 1:2
    subplot(1,2,j); hold on
    i = ind(j)
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
scatter(fly_vel(idx),bump_vel(idx),c,'filled','markerfacealpha',.1)
axis equal
xlim([-6,6]); ylim([-6,6])
y = ylim; x = xlim;
plot(x,[0,0],':','Color',c); 
plot([0,0],y,':','Color',c);
b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); %fit the slope of fly vel and bump vel with an arbitrary offset
r = corr(fly_vel(idx),bump_vel(idx));
plot([0,0],y,':','Color', c)
plot(x,[0,0],':','Color', c)
plot(x,x*b(2) + b(1),'r')
text(x(2),y(1),sprintf('gain: %.2f\nr: %.2f',b(2),r),'HorizontalAlignment','right','VerticalAlignment','bottom','color',c)
xlim(x); ylim(y);
yticks([-5,0,5])
end
linkaxes(get(gcf,'Children'))

%% thrnai figures
i = 0;

figure(1); clf; 
subplot(3,3,1)
c='k';
hold on
errorbar(0*ones(sum(empty_idx & walk_idx & dark_idx==i),1),mean_g(empty_idx & walk_idx & dark_idx==i),sem_g(empty_idx & walk_idx & dark_idx==i),'o','Color',[0,.5,1])
errorbar(1*ones(sum(~empty_idx & walk_idx & dark_idx==i),1),mean_g(~empty_idx & walk_idx & dark_idx==i),sem_g(~empty_idx & walk_idx & dark_idx==i),'o','Color',[1,.5,0])
xticks([0,1]); xticklabels({'Empty','LPsP'}); ylabel('Integrative Gain')
plot(xlim,[.8,.8],':k')
axis padded; set(gca,'Color','none','ycolor',c,'xcolor',c)

legend(sprintf('empty>TH-RNAi (%i)',length(unique(fly_num(walk_idx &empty_idx & dark_idx == i)))),...
   sprintf('lpsp>TH-RNAi (%i)',length(unique(fly_num(walk_idx & ~empty_idx & dark_idx == i)))),...
   'textcolor',c,'Location','best')


subplot(3,3,2); hold on
 scatter(0*ones(sum(empty_idx & walk_idx & dark_idx==i),1),var_o(empty_idx & walk_idx & dark_idx==i),'o','Color',[0,.5,1])
scatter(1*ones(sum(~empty_idx & walk_idx & dark_idx==i),1),var_o(~empty_idx & walk_idx & dark_idx==i),'o','Color',[1,.5,0])
xticks([0,1]); xticklabels({'Empty','LPsP'}); ylabel('Offset Variance (Circular)')
axis padded; set(gca,'Color','none','ycolor',c,'xcolor',c)




ind=[248,189];
b = {[linspace(1,0,255)',linspace(1,0.7,255)',linspace(1,.7,255)'];...
    [linspace(1,1,255)',linspace(1,0.5,255)',linspace(1,0,255)']};
for j = 1:2
    i = ind(j);
a1 = subplot(6,1,3+2*(j-1));

imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
set(gca,'Colormap',b{j},'CLim',clim+[0.2,0]*range(clim))
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'k','linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
idx = round(all_data(i).ft.cue,4) == -.2945;

a = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color',b{j}(end,:),'linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)
xlabel('time (s)')


a2 = subplot(6,1,4+2*(j-1)); hold on
offset = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
a=plot(all_data(i).ft.xf,offset); a.YData(abs(diff(a.YData))>pi) =nan;
 patch(all_data(i).ft.xf,2*pi*(all_data(i).ft.stims/10)-pi,'r','FaceAlpha',.1,'EdgeColor','none')
ylabel('offset')
a2.YTick = [-pi,0,pi]; a2.YTickLabels = {'-\pi','0','\pi'}; a2.YLim = [-pi,pi];
pos = get(gca,'Position');
pos = [pos(1)+pos(3)+.01,pos(2),.05,pos(4)];
ax = axes('Position',pos,'Color','none','XAxisLocation','top');
histogram(offset,-pi:.1:pi,'Orientation','horizontal','edgeColor','none')
box(ax,'off')
ax.YAxisLocation =  'right'; ax.YLim = [-pi,pi]; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};
linkaxes([a1,a2],'x')
xlim(a1,[200,400])
end