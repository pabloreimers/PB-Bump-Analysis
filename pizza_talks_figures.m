%% compartmentalization figures
load('C:\Users\preim\Documents\GitHub\PB-Bump-Analysis\.data\lpsp_cl_data_20240206 (1).mat'); %load in the all_data
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

figure(1); clf
subplot(4,1,1)
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'Color','w','linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
xticks(0:10:100); xticklabels([])
yticks([-pi,pi]); yticklabels({'-\pi','\pi'}); ylim([-pi,pi]); set(gca,'YDir','reverse')
ylabel('heading','Rotation',0)

subplot(4,1,2:3)
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
xticks([]); 
pos = get(gca,'Position'); colorbar; set(gca,'Position',pos,'Colormap',cmap,'CLim',clim+[.4,0]*range(clim))
ylabel({'z-scored','\DeltaF/F'},'Rotation',0)


subplot(4,1,4)
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'Color','w','linewidth',2); a.YData(abs(diff(a.YData))>pi) = nan;
yticks([-pi,pi]); yticklabels({'-\pi','\pi'}); ylim([-pi,pi]); set(gca,'YDir','reverse')
ylabel('bump position','Rotation',0)
xticks(0:10:100); xticklabels([])
xlabel('time (s)')

linkaxes(get(gcf,'Children'),'x')
xlim([0,100])
set(get(gcf,'Children'),'color','none','xcolor','w','ycolor','w')
set(gcf,'Color','none','InvertHardCopy','off')
fontsize(gcf,20,'pixels')


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
% a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'k','linewidth',2);
% a.YData(abs(diff(a.YData))>pi) = nan;
pos = get(gca,'Position');
a = colorbar; a.Color = 'w'; ylabel(a,'\DeltaF/F','Rotation',0)
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

