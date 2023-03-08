%% load data
clear all
load(uigetfile('Select Processed Imaging Data'))

%% assign fly ids
files_info = cellfun(@(x)(regexp(x,'_','split')),files_list(:,1),'UniformOutput',false);
[~,~,fly_id] = unique(str2num([cell2mat(cellfun(@(x)x{1}(1:8),files_info,'UniformOutput',false)),...
                                         cellfun(@(x)x{end},files_info),...
                                         num2str(cellfun(@(x)contains(x{2},'LPsP'),files_info))]))

%% Plot the r2 of linear fits for forward, rotation, and joint velocity fits
tmp_idx = lpsp_idx & include_idx;
tmp = T(tmp_idx,:);
title_str = 'LPsP>syt7f dF/F vs velocity fits';
num_fly = length(unique(fly_id(tmp_idx)));

colormap(hsv(num_fly))
tmp_x = [1,2,3].*[ones(height(tmp),3)];
figure(1); clf
swarmchart(tmp_x(:),tmp{:,end-3:end-1}(:),[],'filled','xjitterwidth',.25,'MarkerFaceColor',.5*[1,1,1])
%swarmchart(tmp_x(:),tmp{:,end-2:end}(:),[],repmat(fly_id(tmp_idx),3,1),'filled','xjitterwidth',.25);%,'MarkerFaceColor',.5*[1,1,1])
xticks([1,2,3])
xticklabels({'Forward','Rotational','Joint'})
ylabel('R^2')
title(title_str)
num_trial = size(tmp_x,1);

x = xlim;
y = ylim;
text(x(2),y(1),sprintf('%i trials from %i flies (both closed loop and dark)',num_trial, num_fly),...
    'HorizontalAlignment','right','VerticalAlignment','bottom')

%% Plot the correlation between bump and fly velocity and bootstrap
N = 1e4;
tmp_idx     = lpsp_idx & include_idx;
title_str   = 'LPsP>syt7f Bump and Fly Velocity Correlation';

figure(2); clf; hold on
[~,~,tmp]   = unique(fly_id(tmp_idx));
d_idx       = ismember(find(tmp_idx),find(dark_idx));
tmp_data    = T.vel_rho(tmp_idx);
scatter(tmp(~d_idx),tmp_data(~d_idx),'filled', 'MarkerFaceAlpha',0.75, 'MarkerFaceColor',[0,0.75,0]);
scatter(tmp(d_idx),tmp_data(d_idx),'filled', 'MarkerFaceAlpha',.75, 'MarkerFaceColor',[.5,.5,.5]);

ylabel('Pearson Correlation Coefficient')
xlabel('Fly #')
title(title_str)
xticks(1:max(tmp))

%bootstrap whether greater than 0
tmp = T.vel_rho(tmp_idx&~dark_idx);
idx = randi(length(tmp),[N,length(tmp)]);
p_L = sum((mean(tmp(idx),2) < 0)) / N;

tmp = T.vel_rho(tmp_idx&dark_idx);
idx = randi(length(tmp),[N,length(tmp)]);
p_D = sum((mean(tmp(idx),2) < 0)) / N;

str1 = sprintf('Cue (p < 10^{%i})',ceil(log10(p_L)));
str2 = sprintf('Dark  (p < 10^{%i})',ceil(log10(p_D)));
legend(str1,str2,'Location','bestoutside','AutoUpdate','off')
pos = get(gca,'Position');
%legend('Cue','Dark','Location','bestoutside')
%set(gca,'Position',pos)
x = xlim;
y = ylim;
plot(x,[0,0],':k')
text(x(2),y(1),sprintf('p-values from mean of %i bootstrap resamples > 0',N),'HorizontalAlignment','right','VerticalAlignment','bottom')

%% Plot the correlation between bump and fly velocity and bootstrap (one point per fly)
N = 1e4;
tmp_idx     = lpsp_idx & include_idx;
title_str   = 'LPsP>syt7f Bump and Fly Velocity Correlation';

figure(3); clf; hold on
[~,~,tmp_num]   = unique(fly_id(tmp_idx));
tmp_dark        = ismember(find(tmp_idx),find(dark_idx));
tmp_rhos        = T.vel_rho(tmp_idx);

tmp_table = table(tmp_num,tmp_dark,tmp_rhos);
tmp_stats = grpstats(tmp_table,{'tmp_num','tmp_dark'},'mean','DataVars','tmp_rhos');

scatter(tmp_stats.tmp_num(~tmp_stats.tmp_dark),tmp_stats.mean_tmp_rhos(~tmp_stats.tmp_dark),'filled', 'MarkerFaceAlpha',0.75, 'MarkerFaceColor',[0,0.75,0]);
scatter(tmp_stats.tmp_num(tmp_stats.tmp_dark),tmp_stats.mean_tmp_rhos(tmp_stats.tmp_dark),'filled', 'MarkerFaceAlpha',.75, 'MarkerFaceColor',[.5,.5,.5]);

ylabel('Pearson Correlation Coefficient')
xlabel('Fly #')
title(title_str)
xticks(1:max(tmp_stats.tmp_num))

%bootstrap whether greater than 0
tmp = tmp_stats.mean_tmp_rhos(~tmp_stats.tmp_dark);
idx = randi(length(tmp),[N,length(tmp)]);
p_L = sum((mean(tmp(idx),2) < 0)) / N;

tmp = tmp_stats.mean_tmp_rhos(tmp_stats.tmp_dark);
idx = randi(length(tmp),[N,length(tmp)]);
p_D = sum((mean(tmp(idx),2) < 0)) / N;

str1 = sprintf('Cue (p < 10^{%i})',ceil(log10(p_L)));
str2 = sprintf('Dark  (p < 10^{%i})',ceil(log10(p_D)));
legend(str1,str2,'Location','bestoutside','AutoUpdate','off')
pos = get(gca,'Position');
%legend('Cue','Dark','Location','bestoutside')
%set(gca,'Position',pos)
x = xlim;
y = ylim;
plot(x,[0,0],':k')
text(x(2),y(1),sprintf('p-values from mean of %i bootstrap resamples > 0',N),'HorizontalAlignment','right','VerticalAlignment','bottom')

figure(4); clf; hold on
for i = unique(tmp_stats.tmp_num)'
    a = plot([0,1],[nan,nan],'-k');
    b = scatter([0,1],[nan,nan],[],[.5,.5,.5;0,.5,0],'filled','MarkerFaceAlpha',0.5);
    if any(tmp_stats.tmp_num == i & tmp_stats.tmp_dark)
    a.YData(1) = tmp_stats.mean_tmp_rhos(tmp_stats.tmp_num == i & tmp_stats.tmp_dark);
    b.YData(1) = tmp_stats.mean_tmp_rhos(tmp_stats.tmp_num == i & tmp_stats.tmp_dark);
    end
    if any(tmp_stats.tmp_num == i & ~tmp_stats.tmp_dark)
    a.YData(2) = tmp_stats.mean_tmp_rhos(tmp_stats.tmp_num == i & ~tmp_stats.tmp_dark); 
    b.YData(2) = tmp_stats.mean_tmp_rhos(tmp_stats.tmp_num == i & ~tmp_stats.tmp_dark);
    end
end
plot([-.5,1.5],[0,0],'k:')
xticks([0,1])
xticklabels({'Dark','Cue'})
xlim([-.5,1.5])
ylabel('Pearson Correlation Coefficient')
text(1.5,0,sprintf('N = %i flies',i),'HorizontalAlignment','right','VerticalAlignment','bottom')
title(title_str)
%%
figure(3)
scatter(1:length(T.r_r2),T.vel_rho)