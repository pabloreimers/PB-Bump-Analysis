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
%load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\lpsp_cl_data_20240206.mat')
%load('Z:\pablo\lpsp_cl\20221122\20221122-1_EPG_GRABDA(2m)_1\registration_001\imagingData_reg_ch1_trial_001.mat')

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
i = 52;
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