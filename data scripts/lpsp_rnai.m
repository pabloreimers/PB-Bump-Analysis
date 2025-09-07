%% create figure to show example
i = 74;
bin_edges = 0:.05:5;
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
ax.YAxisLocation =  'right'; ax.YLim = [-pi,pi]; ax.YTick = [-pi,0,pi]; ax.YTickLabels = {'-\pi','0','\pi'};

a3 = subplot(6,1,4); hold on
plot(all_data(i).ft.xb,all_data(i).gain.inst_g)
ylabel('integrative gain')

linkaxes([a1,a2,a3],'x')
xlim([min(all_data(i).ft.xb),max(all_data(i).ft.xb)])
ylim([0,5])
plot(xlim,[.8,.8],':'); %plot(xlim,[1.6,1.6],':k')

subplot(3,2,5); hold on
h = histogram(all_data(i).gain.g,'BinEdges',binedges,'FaceAlpha',.8,'Normalization','probability','EdgeColor','none');

xlabel('integrative gain')
ylabel('counts')

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
    %set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(g(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[0:.1:2],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(g(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[0:.1:2],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',length(unique(fly_num(empty_idx & dark_idx == i)))),...
           sprintf('lpsp>TH-RNAi (%i)',length(unique(fly_num(~empty_idx & dark_idx == i)))),...
           'textcolor','k')
    title(sprintf('Dark = %i',i),'color','k')
end

title(t,'Integrative Gain','color','k')
%set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% show historams in the CL and the dark (offset)
o = {};
for i = 1:length(all_data)
    tmp_idx = all_data(i).ft.f_speed > f_thresh;
    tmp = circ_dist(-all_data(i).ft.cue,interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),all_data(i).ft.xf));
    tmp = tmp(tmp_idx);
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
    %set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(o(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(o(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
           sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
           'textcolor','k')
    title(sprintf('Dark = %i',i),'Color','k')
end

title(t,'Offset (mean centered)','color','k')
%set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% show historams in the CL and the dark (bump and cue position occupancy)
m = {};
c = {};
for i = 1:length(all_data)
    tmp = all_data(i).ft.cue;
    tmp_idx = all_data(i).ft.f_speed > f_thresh;
    tmp = tmp(tmp_idx);
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    c{i} = tmp(~isnan(tmp));
    
    tmp = all_data(i).im.mu;
    tmp_idx = interp1(all_data(i).ft.xf,all_data(i).ft.f_speed,all_data(i).ft.xb) > f_thresh;
    tmp = tmp(tmp_idx);
    tmp_mean = atan2(mean(sin(tmp),'omitnan'),mean(cos(tmp),'omitnan'));
    tmp = tmp - tmp_mean;
    tmp(tmp<-pi) = tmp(tmp<-pi) + 2*pi;
    m{i} = tmp(~isnan(tmp));
end
c = reshape(c,[],1);
m = reshape(m,[],1);
%c = [c;m];

figure(3); clf
t = tiledlayout(2,2);
for i = 0:1
    nexttile; hold on
   % set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(c(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(c(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    % legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
    %        sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
    %        'textcolor','k')
    title('cue','Color','k')
end

for i = 0:1
    nexttile; hold on
    %set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(m(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(m(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-pi:.1:pi],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    % legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
    %        sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
    %        'textcolor','k')
    title('bump','color','k')
end
title(t,'Position (mean centered)','color','k')
%set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')

%% compare entropy of cue heading, bump position, and their substraction
rho_thresh = .1;
f_thresh = .5;
tmp_idx = ~dark_idx;
bin_edges = -pi:.1:pi;

ent_cue = nan(length(all_data),1);
ent_mu = nan(length(all_data),1);

for i = 1:length(all_data)
    idx = abs(all_data(i).ft.f_speed) > f_thresh & interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf) > rho_thresh;
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
plot(ent(empty_idx & tmp_idx,:)','Color',[0,.5,1])
plot(ent(~empty_idx & ~tmp_idx,:)','Color',[1,.5,0])
axis padded

nexttile; hold on
scatter(1*ones(sum(empty_idx & tmp_idx),1),ent_mu(empty_idx & tmp_idx) ./ ent_cue(empty_idx & tmp_idx),[],[0,.5,1])
scatter(2*ones(sum(~empty_idx & tmp_idx),1),ent_mu(~empty_idx & tmp_idx) ./ ent_cue(~empty_idx & tmp_idx),[],[1,.5,0])
plot(xlim,[1,1],':k')
axis padded


t.Children(2).XTick = [1,2]; t.Children(2).XTickLabel = {'Cue','Mu'}; title(t.Children(2), 'Entropy');
t.Children(1).XTick = [1,2]; t.Children(1).XTickLabel = {'Empty','LPsP'}; title(t.Children(1), 'Ratio (mu/cue)');

%% show scatter of bump and fly vel for every fly
figure(8); clf
rows = ceil(sqrt(length(all_data)));
cols = ceil(length(all_data)/rows);

% rows = 2;
% cols = ceil(length(all_data)/rows);

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

%% compare behavior of the two groups

f_thresh = .1;
r = nan(length(all_data),1);
f = nan(length(all_data),1);
p = nan(length(all_data),1);


for i = 1:length(all_data)
    tmp_idx = all_data(i).ft.f_speed > f_thresh;

    r(i) = mean(abs(all_data(i).ft.r_speed(tmp_idx)),'omitnan');
    f(i) = mean(abs(all_data(i).ft.f_speed(tmp_idx)),'omitnan');
    p(i) = sum(all_data(i).ft.f_speed > f_thresh) / length(all_data(i).ft.f_speed);
end

figure(6); clf
nexttile
scatter(group_idx,r,[],c(empty_idx+1,:))
title('mean r speed'); axis padded
xlim([-.5,3.5]); xticks([0:3]); xticklabels({'LPsP\newlineCL','Empty\newlineCL','LPsP\newlineDark','Empty\newlineDark'});

nexttile
scatter(group_idx,f,[],c(empty_idx+1,:))
title('mean f speed'); axis padded
xlim([-.5,3.5]); xticks([0:3]); xticklabels({'LPsP\newlineCL','Empty\newlineCL','LPsP\newlineDark','Empty\newlineDark'});

nexttile
scatter(group_idx,p,[],c(empty_idx+1,:))
title('walking time (% of trial)'); axis padded
xlim([-.5,3.5]); xticks([0:3]); xticklabels({'LPsP\newlineCL','Empty\newlineCL','LPsP\newlineDark','Empty\newlineDark'});

%%
%% show historams in the CL and the dark (offset)
f_thresh = -inf;
f = {};
for i = 1:length(all_data)
    tmp_idx = all_data(i).ft.f_speed > f_thresh;
    tmp = all_data(i).ft.f_speed;
    tmp = tmp(tmp_idx);
    f{i} = tmp(~isnan(tmp));
end
f = reshape(f,[],1);

figure(3); clf
t = tiledlayout(1,2);
for i = 0:1
    nexttile; hold on
    %set(gca,'color','none','ycolor','w','xcolor','w')
    tmp = reshape(cell2mat(f(empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-2:.1:10],'Normalization','Probability','FaceColor',[0,.5,1],'FaceAlpha',.8)
    tmp = reshape(cell2mat(f(~empty_idx & dark_idx == i)),1,[]);
    histogram(tmp(~isnan(tmp)),'BinEdges',[-2:.1:10],'Normalization','Probability','FaceColor',[1,.5,0],'FaceAlpha',.8)
    legend(sprintf('empty>TH-RNAi (%i)',sum(empty_idx & dark_idx == i)),...
           sprintf('lpsp>TH-RNAi (%i)',sum(~empty_idx & dark_idx == i)),...
           'textcolor','k')
    title(sprintf('Dark = %i',i),'Color','k')
end

title(t,'Offset (mean centered)','color','k')
%set(gcf,'color','none','InvertHardcopy','off')
fontsize(gcf,20,'pixels')