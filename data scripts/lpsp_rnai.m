%% create figure to show example
i = 16;
binedges = 0:.05:5;
dark_mode = false;

figure(1); clf
a1 = subplot(3,1,1);
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
hold on
if contains(all_data(i).ft.pattern,'background'); c = 'm'; else; c = 'c'; end
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c); a.YData(abs(diff(a.YData))>pi) = nan;
idx = round(all_data(i).ft.cue,4) == -.2945;

a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w'); a.YData(abs(diff(a.YData))>pi) = nan;
title(all_data(i).meta)
xlabel('time (s)')

a2 = subplot(6,1,3); hold on
offset = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
a=plot(all_data(i).ft.xf,offset); a.YData(abs(diff(a.YData))>pi) =nan;
 patch(all_data(i).ft.xf,2*pi*(all_data(i).ft.stims/10)-pi,'r','FaceAlpha',.1,'EdgeColor','none')
ylabel('offset')
a2.YTick = [-pi,0,pi]; a2.YTickLabels = {'-\pi','0','\pi'}; a2.YLim = [-pi,pi];
pos = get(gca,'Position');
pos = [pos(1)+pos(3)+.01,pos(2),.05,pos(4)];
ax = axes('Position',pos,'Color','none','XAxisLocation','top');
histogram(offset,bin_edges,'Orientation','horizontal','edgeColor','none')
box(ax,'off')
ax.YAxisLocation =  'right'; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};

a3 = subplot(6,1,4); hold on
plot(all_data(i).ft.xb,all_data(i).gain.inst_g)
ylabel('integrative gain')

linkaxes([a1,a2,a3],'x')
xlim([min(all_data(i).ft.xb),max(all_data(i).ft.xb)])
ylim([0,5])
plot(xlim,[.8,.8],':'); %plot(xlim,[1.6,1.6],':k')

subplot(3,2,5); hold on
h = histogram(all_data(i).gain.g,'BinEdges',binedges,'FaceAlpha',.8,'Normalization','probability','EdgeColor','none');

xlabel('gain')
ylabel('counts')
legend('integrative','color','none','textcolor','w')

subplot(3,2,6); 
scatter(all_data(i).gain.g,all_data(i).gain.v,'filled','MarkerFaceAlpha',.5)
xlabel('integrative gain')
ylabel('func value')

%% show historams in the CL and the dark (integrative gain)
g = {};
for i = 1:length(all_data)
    g{i} = all_data(i).gain.g;
end
g = reshape(g,[],1);

figure(3); clf
t = tiledlayout(1,2);
for i = 0:1
    nexttile; hold on
    set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(g(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[0:.1:2],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(g(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[0:.1:2],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
           sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
           'textcolor','w')
    title(sprintf('Dark = %i',i),'Color','w')
end

title(t,'Integrative Gain','color','w')
set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% show historams in the CL and the dark (offset)
o = {};
for i = 1:length(all_data)
    tmp = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    o{i} = tmp(~isnan(tmp));
end
o = reshape(o,[],1);

figure(3); clf
t = tiledlayout(1,2);
for i = 0:1
    nexttile; hold on
    set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(o(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(o(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
           sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
           'textcolor','w')
    title(sprintf('Dark = %i',i),'Color','w')
end

title(t,'Offset','color','w')
set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% show historams in the CL and the dark (bump and cue position occupancy)
m = {};
c = {};
for i = 1:length(all_data)
    tmp = all_data(i).ft.cue;
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    c{i} = tmp(~isnan(tmp));
    
    tmp = all_data(i).im.mu;
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    m{i} = tmp(~isnan(tmp));
end
c = reshape(c,[],1);
m = reshape(m,[],1);

figure(3); clf
t = tiledlayout(1,2);
for i = 0:1
    nexttile; hold on
    set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(m(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(m(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
           sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
           'textcolor','w')
    title(sprintf('Dark = %i',i),'Color','w')
end

title(t,'Cue Position','color','w')
set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% compare entropy of cue heading, bump position, and their substraction
bin_edges = -pi:.1:pi;

ent_cue = nan(length(all_data),1);
ent_mu = nan(length(all_data),1);

for i = 1:length(all_data)
    idx = abs(all_data(i).ft.f_speed) > .01;
    ent_cue(i) = norm_entropy(all_data(i).ft.cue(idx),bin_edges);
    tmp = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf);
    tmp = mod(tmp,2*pi);
    tmp(tmp>pi) = tmp(tmp>pi) - 2*pi;
    ent_mu(i) = norm_entropy(tmp(idx),bin_edges);
end

ent = [ent_cue,ent_mu];

figure(4); clf
t = tiledlayout(1,2);
nexttile; hold on
plot(ent(empty_idx & ~dark_idx,:)','Color',[0,.5,1])
plot(ent(~empty_idx & ~dark_idx,:)','Color',[1,.5,0])

nexttile; hold on
scatter(1*ones(sum(empty_idx),1),ent_mu(empty_idx) ./ ent_cue(empty_idx),[],[0,.5,1])
scatter(2*ones(sum(~empty_idx),1),ent_mu(~empty_idx) ./ ent_cue(~empty_idx),[],[1,.5,0])

t.Children(2).XTick = [1,2]; t.Children(2).XTickLabel = {'Cue','Mu'}; title(t.Children(2), 'Entropy');
t.Children(1).XTick = [1,2]; t.Children(1).XTickLabel = {'Empty','LPsP'}; title(t.Children(1), 'Difference');

%% show scatter of bump and fly vel for every fly
figure(8); clf
rows = ceil(sqrt(length(all_data)));
cols = ceil(length(all_data)/rows);

rows = 2;
cols = ceil(length(all_data)/rows);

t = tiledlayout(rows,cols);
lag = 10;
g = nan(length(all_data),1);

for i = 1:length(all_data)
    nexttile; hold on

    dm = gradient(interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf)) * 60;
    dr = all_data(i).ft.r_speed;

    dm = smoothdata(dm(1+lag:end),1,'movmean',5);
    dr = smoothdata(dr(1:end-lag),1,'movmean',5);
    
    idx = abs(dr) > .1 & ~isnan(dm);
    if empty_idx(i); c = [0,.5,1]; else; c=[1,.5,0];end
    scatter(dr(idx),dm(idx),5,'filled','MarkerFaceColor',c)
    
    b = [ones(sum(idx),1),dr(idx)] \ dm(idx);
    g(i) = b(2);
    plot([-2,2],b(1)+g(i)*[-2,2],'r')
    text(max(xlim),min(ylim),sprintf('gain: %.2f',g(i)),'HorizontalAlignment','right','VerticalAlignment','bottom')
end

xlabel(t,'Fly Speed (rad/s)'); ylabel(t,'Bump Speed (rad/s)')

nexttile
c = [1,.5,0; 0,.5,1];
group_idx = empty_idx + 2*dark_idx;
scatter(group_idx,g,[],c(empty_idx+1,:));ylabel({'instantaneous gain','lag = 160ms'},'Rotation',0); xlim([-.5,3.5]); xticks([0:3]); xticklabels({'LPsP\newlineCL','Empty\newlineCL','LPsP\newlineDark','Empty\newlineDark'}); set(gca,'YAxisLocation','right')