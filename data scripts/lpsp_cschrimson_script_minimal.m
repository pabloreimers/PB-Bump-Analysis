%% start fresh
 % close all
 % clear all

%% load in data
base_dir = ('Z:\pablo\lpsp_cschrimson\to do\curr'); %uigetdir(); %
all_files = dir([base_dir,'\**\*imgData_smooth_reg*.mat']);
all_files = natsortfiles(all_files);
%% make sure that each file has a mask
for i = 1:length(all_files)
    fprintf('checking mask: %s\n',all_files(i).folder)
    clear img regProduct 

    if ~isfile([fileparts(all_files(i).folder),'\mask_smooth_reg.mat'])
        load([all_files(i).folder,'\',all_files(i).name])
        
        imgData = imgData_smooth_reg;

        top_pct = prctile(imgData,98,'all');
        bot_pct = prctile(imgData,5,'all');
        
        imgData(imgData>top_pct) = top_pct;
        imgData(imgData<bot_pct) = bot_pct;        

        figure(1); clf; imagesc(mean(imgData,3)); colormap(bone); axis equal tight; drawnow;
        mask = roipoly();
        save([fileparts(all_files(i).folder),'\mask_smooth_reg.mat'],'mask')
    end
end

%% process and store all values
ft_type= 'gaussian'; %the type of smoothing for fictrac data
ft_win = 60; %the window over which smoothing of fictrac data occurs. gaussian windows have std = win/5.
im_type= {'gaussian','gaussian'}; %there's two smoothing steps for the im data. one that smooths the summed z-stacks, another that smooths the estimated mu and rho
im_win = {10,10};
n_centroid = 16;
f0_pct = 7;

%all_data = struct();

tic
for i = length(all_data):length(all_files)
    clear img regProduct imgData_reg

    try
    tmp = strsplit(all_files(i).folder,'\');
    fprintf('processing: %s ',tmp{end-1})
    load([all_files(i).folder,'\',all_files(i).name])
    load([fileparts(all_files(i).folder),'\mask_smooth_reg.mat'])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])
    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_dat.mat']);
    load([tmp2.folder,'\',tmp2.name])
    % tmp3 = dir([fileparts(all_files(i).folder),'\FicTracData\*.avi']);
    % vidObj = VideoReader([tmp3.folder,'\',tmp3.name]);
    % vid  = read(vidObj);

    
    imgData = imgData_smooth_reg;

    all_data(i).ft = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type);
    all_data(i).im = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct);
    all_data(i).meta = all_files(i).folder;

    if ~ismember('xb',fieldnames(all_data(i).ft))
        xb = linspace(all_data(i).ft.xf(1),all_data(i).ft.xf(end),size(all_data(i).im.d,2));
    end

    %all_data(i).atp = process_im(img{2}, im_win, im_type, mask, n_centroid, f0_pct);

    % lasing = squeeze(sum(vid(1:end/3,:,:,:),[1,2,3]));
    % ind1 = find(lasing > mean(lasing),1,'first');
    % ind2 = find(lasing > mean(lasing),1,'last');
    % 
    % all_data(i).ft.lasing = squeeze(sum(vid(1:end/5,:,:,:),[1,2,3]));
    % all_data(i).ft.stims = squeeze(sum(vid(1:end/5,1:end/5,:,:),[1,2,3]));
    try
    all_data(i).ft.stims = ftData_DAQ.stim{1};
    end
    end
    fprintf('ETR: %.2f hours\n',toc/i * (length(all_files)-i) / 60 / 60)
end

%% compare gains
tmp_win = 30;
fr = 1;
r_thresh = .1;
rho_thresh = 0;

gain_vel = nan(length(all_data),60);
gain_pos = cell(length(all_data),2);
first_idx = false(length(all_data),1);
dangle_idx= false(length(all_data),1);

last_str = 'blank';

for i = 1:length(all_data)

    tmp_str = all_data(i).meta(1:50);
    
    if ~strcmp(tmp_str,last_str) && circ_var(circ_dist(cumsum(all_data(i).ft.r_speed)/60*.8,-unwrap(all_data(i).ft.cue))) < .5
        first_idx(i) = true;
        last_str = tmp_str;
    end

    if  circ_var(circ_dist(cumsum(all_data(i).ft.r_speed)/60*.8,-unwrap(all_data(i).ft.cue))) > .5
        dangle_idx(i) = true;
    end

    dc = gradient(unwrap(all_data(i).ft.cue)) * 60;
    dm = interp1(all_data(i).ft.xb,gradient(unwrap(all_data(i).im.mu)),all_data(i).ft.xf) * 60;
    dr = all_data(i).ft.r_speed;
    rho = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);

    for lag = 1:size(gain_vel,2)
        tmp_dr = dr(1:end-lag);
        tmp_dm = dm(lag+1:end);
        tmp_rho= rho(lag+1:end);
        idx = abs(tmp_dr) > r_thresh &  tmp_rho > rho_thresh & ~isnan(tmp_dm);
        
        gain_vel(i,lag) = tmp_dr(idx) \ tmp_dm(idx);
    end


    % tmp_mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),1:floor(all_data(i).ft.xb(end))); %interpolate mu for 1 timepoint per second
    % tmp_rs = interp1(all_data(i).ft.xf,all_data(i).ft.r_speed,1:floor(all_data(i).ft.xb(end)));
    % tmp_gain = nan(length(tmp_mu)-tmp_win,2);
    % tmp_fval = nan(length(tmp_mu)-tmp_win,1);
    % tmp_mov  = nan(length(tmp_mu)-tmp_win,1);

    tmp_t = floor(all_data(i).ft.xb(1)):fr:ceil(all_data(i).ft.xb(end));
    tmp1 = unwrap(all_data(i).im.mu);
    tmp1 = interp1(all_data(i).ft.xb,tmp1,tmp_t,'linear','extrap');
    tmp2 = cumsum(all_data(i).ft.r_speed)/60;
    tmp2 = interp1(all_data(i).ft.xf,tmp2,tmp_t,'linear','extrap');

    tmp_gain = nan(length(tmp_t)-tmp_win,2);
    tmp_fval = nan(length(tmp_t)-tmp_win,1);
    tmp_dist = nan(length(tmp_t)-tmp_win,1);
    
    for t = 1:(length(tmp_t) - tmp_win/fr)
        % fun = @(x)(circ_var(circ_dist(...
        %            cumsum(tmp_rs(t:t+tmp_win))*x,tmp_mu(t:t+tmp_win))'));
        % [tmp_gain(t,:),tmp_fval(t)] = fminsearch(fun,1);
        % tmp_mov(i) = sum(abs(tmp_rs(i:i+tmp_win)));
    % 
        
    fun = @(x)(sum((...
           tmp2(t:t+(tmp_win/fr))*x(1)+x(2) - tmp1(t:t+(tmp_win/fr))).^2));
    

    [tmp_gain(t,:),tmp_fval(t)] = fminsearch(fun,[1,0]);
    tmp_dist(t) = sum(abs(diff(tmp2(t:t+tmp_win))));
        
    end
    fprintf('fly: %i\n',i)
    all_data(i).gain.g = tmp_gain(:,1);
    all_data(i).gain.b = tmp_gain(:,2);
    all_data(i).gain.f = tmp_fval;
    all_data(i).gain.d = tmp_dist;

    i / length(all_data)
end

empty_idx = cellfun(@(x)(contains(x,'empty')),{all_data.meta});


%% show dff with flashes
tmp_str = '20250128\fly 1';
trial_num = 2;

tmp_ind = find(cellfun(@(x)(contains(x,tmp_str)),{all_data.meta}'));
i = tmp_ind(trial_num);
%i = 19;
r_thresh = .1;
rho_thresh = 0;

% ind1 = find(all_data(i).ft.lasing > mean(all_data(i).ft.lasing),1,'first');
% ind2 = find(all_data(i).ft.lasing > mean(all_data(i).ft.lasing),1,'last');
% stims = interp1(1:(ind2-ind1+1),all_data(i).ft.stims(ind1:ind2),1:length(all_data(i).ft.xf),'linear','extrap');
%stims = interp1(1:length(all_data(i).ft.stims),stims,1:length(all_data(i).ft.xf),'linear','extrap');
stims = all_data(i).ft.stims/10;

figure(1); clf
subplot(6,1,1); hold on
title(all_data(i).meta)
plot(all_data(i).ft.xf,all_data(i).ft.f_speed)
plot(all_data(i).ft.xf,all_data(i).ft.r_speed)
legend('forward','rotation')
ylabel('speed')

subplot(6,1,2)
plot(all_data(i).ft.xf,stims)
hold on
tmp_gain = all_data(i).gain.g;
tmp_gain(all_data(i).gain.f > 100) = nan;
scatter(1:length(tmp_gain),tmp_gain,'.')
plot(xlim,[.8,.8],':','Color',dark_mode*[1,1,1])
ylim([-5,10])
ylabel('instantaneous gain')

subplot(3,1,2); hold on
tmp = all_data(i).im.z;
%idx = sum(tmp,1) > mean(sum(tmp,1)) + std(sum(tmp,1));
%tmp(:,idx) = nan;
%tmp = tmp ./ sum(tmp,1);
%tmp_idx = sum(all_data(i).im.d,1) > 100;
%tmp(:,tmp_idx) = nan;
imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp)
a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,'m','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;
a = plot(all_data(i).ft.xb,all_data(i).im.mu,'w','linewidth',2);
a.YData(abs(diff(a.YData))>pi) = nan;

linkaxes(get(gcf,'Children'),'x')
axis tight
colormap(parula)

subplot(3,2,5) % instantaneous gain
dc = gradient(unwrap(all_data(i).ft.cue))*60;
dm = interp1(all_data(i).ft.xb,gradient(unwrap(all_data(i).im.mu)),all_data(i).ft.xf) * 60;
m  = interp1(all_data(i).ft.xb,all_data(i).im.mu,all_data(i).ft.xf);
idx = abs(all_data(i).ft.r_speed) > r_thresh & interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf) > rho_thresh & ~isnan(dm);
scatter(all_data(i).ft.r_speed(idx),dm(idx),10,m(idx),'filled','MarkerFaceAlpha',.2)
axis tight
hold on
gl = polyfit(all_data(i).ft.r_speed(idx & all_data(i).ft.r_speed>0),dm(idx & all_data(i).ft.r_speed>0),1);
gr = polyfit(all_data(i).ft.r_speed(idx & all_data(i).ft.r_speed<0), dm(idx & all_data(i).ft.r_speed<0),1);
x = xlim;
plot([0,x(2)],gl(1)*[0,x(2)] + gl(2),'r','linewidth',2)
plot([x(1),0],gr(1)*[x(1),0] + gr(2),'g','linewidth',2)
text(max(xlim),max(ylim),sprintf('gain L: %.2f\n gain R: %.2f',gl(1),gr(1)),'Color',(dark_mode)*[1,1,1])
xlabel('fly speed'); ylabel('bump speed')
title('instantaeous gain')


subplot(3,2,6) % positional gain
tmp_win = 20;
tmp_mu = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),1:floor(all_data(i).ft.xb(end))); %interpolate mu for 1 timepoint per second
tmp_rs = interp1(all_data(i).ft.xf,all_data(i).ft.r_speed,1:floor(all_data(i).ft.xb(end)));
tmp_gain = nan(length(tmp_mu)-tmp_win,2);
tmp_fval = tmp_gain;
tmp_mov  = tmp_gain;
j = i;
for i = 1:length(tmp_gain)
    fun = @(x)(circ_var(circ_dist(...
               cumsum(tmp_rs(i:i+tmp_win))*x,tmp_mu(i:i+tmp_win))'));
% 


    %fun = @(x)(sum(cumsum(tmp_rs(i:i+tmp_win) + x(1))*x(2) - tmp_mu(i:i+tmp_win)).^2);
    [tmp_gain(i,:),tmp_fval(i)] = fminsearch(fun,1);
    tmp_mov(i) = sum(abs(tmp_rs(i:i+tmp_win)));
    fprintf('time: %i\n',i)
end
tmp_gain(tmp_fval>prctile(tmp_fval,50)) = nan;
tmp_gain2 = all_data(j).gain.g;
tmp_gain2(all_data(j).gain.f > median(all_data(j).gain.f)) = nan;
h = histogram(tmp_gain(:,1),'binwidth',.1,'Normalization','probability');
hold on; histogram(tmp_gain2,'binwidth',.1,'Normalization','probability');
title('integrative gain')
legend('calc1','calc2','autoupdate','off')

%scatter(tmp_fval,tmp_gain)
hold on
plot([1,1],ylim,':k')

%xlim([-2,10])

ax = get(gcf,'Children');
if dark_mode
    set(gcf,'color','none','InvertHardcopy','off')
    
    for i = 1:length(ax)
        if contains(class(ax(i)),'Axes')
            ax(i).Title.Color = 'w';
            set(ax(i),'Color','none','ycolor','w','xcolor','w')
        else
            ax(i).Color = 'w';
        end
    end
end

%% show positional gain histogram broken down by type
gain_empty = {};
fval_empty = {};
bias_empty = {};
dist_empty = {};        
gain_chrim = {};
fval_chrim = {};
bias_chrim = {};
dist_chrim = {};


curr_fly = '';
num_empty = 0;
num_chrim = 0;
for i = 1:length(all_data)
    if strcmp(curr_fly,all_data(i).meta(1:50)) %~any(all_data(i).ft.stims)
        continue
    end
    
    curr_fly = all_data(i).meta(1:50);

    if contains(all_data(i).meta,'empty') | contains(all_data(i).meta,'+')
        num_empty = num_empty+1;
        gain_empty = [gain_empty;{all_data(i).gain.g}];
        fval_empty = [fval_empty;{all_data(i).gain.f}];
        bias_empty = [bias_empty;{all_data(i).gain.b}];
        dist_empty = [dist_empty;{all_data(i).gain.d}];

    else
        num_chrim = num_chrim+1;
        gain_chrim = [gain_chrim;{all_data(i).gain.g}];
        fval_chrim = [fval_chrim;{all_data(i).gain.f}];
        bias_chrim = [bias_chrim;{all_data(i).gain.b}];
        dist_chrim = [dist_chrim;{all_data(i).gain.d}];
    end
end

gain_empty = cell2mat(gain_empty);
fval_empty = cell2mat(fval_empty);
bias_empty = cell2mat(bias_empty);
dist_empty = cell2mat(dist_empty);        
gain_chrim = cell2mat(gain_chrim);
fval_chrim = cell2mat(fval_chrim);
bias_chrim = cell2mat(bias_chrim);
dist_chrim = cell2mat(dist_chrim);

%
fval_max = 100;
dist_thresh = 10;
figure(8); clf
subplot(3,3,1); hold on; histogram(bias_empty); xlabel('bias'); ylabel('empty')
subplot(3,3,2); hold on; histogram(dist_empty); plot(dist_thresh*[1,1],ylim,'r'); xlabel('dist'); 
subplot(3,3,3); hold on; histogram(fval_empty); plot(fval_max*[1,1],ylim,'r'); xlabel('fval');

subplot(3,1,2); hold on
histogram(gain_empty(fval_empty<fval_max & dist_empty>dist_thresh),'Normalization','probability','Binwidth',.1)
histogram(gain_chrim(fval_chrim<fval_max & dist_chrim>dist_thresh),'Normalization','probability','Binwidth',.1)
legend(sprintf('empty (%d)',num_empty),sprintf('cschrimson (%d)',num_chrim))
xlabel({'integrative gain','baseline (30s window)'})
ylabel('%')
%% show movie
% i = 19;
% pause_time = 1e-2;
% 
% tmp_cue = interp1(all_data(i).ft.xf,all_data(i).ft.cue,all_data(i).ft.xb);
% tmp_img = smoothdata(imgData_reg,3,'movmean',10);
% tmp_img = 256*(tmp_img - min(tmp_img,[],'all')) / (max(tmp_img,[],'all') - min(tmp_img,[],'all')); 
% 
% figure(8); clf
% subplot(2,1,1)
% h(1) = image(tmp_img(:,:,1));
% subplot(2,1,2); hold on
% h(2) = plot(unwrap(all_data(i).im.alpha),all_data(i).im.z(:,1));
% h(3) = scatter(-tmp_cue(1),0,'r','filled');
% ylim([-1,5])
% 
% for ii = 1:length(all_data(i).ft.xb)
%     h(1).CData = tmp_img(:,:,ii);
%     h(2).YData = all_data(i).im.z(:,ii);
%     h(3).XData = -tmp_cue(ii);
%     title(all_data(i).ft.xb(ii))
%     drawnow
%     pause(pause_time);
% end


%%
%first_idx = second_idx;

figure(8); clf
hold on
c = cool(size(gain_vel,2));
for i = 1:size(gain_vel,2)
    scatter(empty_idx(first_idx)+(.5/size(gain_vel,2))*i,gain_vel(first_idx,i),'.','CData',c(i,:))
end
xlim([-1,2])
ylim([-1,10])
ylabel('gain')
xticks([0,1]); xticklabels({'lpsp','empty'})

figure(9); clf
rows = floor(sqrt(sum(first_idx)));
cols = ceil(sum(first_idx)/rows);
tiledlayout(rows, cols)

lag = 20;
for i = find(first_idx)'
    
    dc = gradient(unwrap(all_data(i).ft.cue))  * 60;
    dm = interp1(all_data(i).ft.xb,gradient(unwrap(all_data(i).im.mu)),all_data(i).ft.xf) * 60;
    dr = all_data(i).ft.r_speed;
    rho = interp1(all_data(i).ft.xb,all_data(i).im.rho,all_data(i).ft.xf);

    tmp_dr = dr(1:end-lag);
    tmp_dm = dm(lag+1:end);
    tmp_rho= rho(lag+1:end);
    idx = abs(tmp_dr) > r_thresh &  tmp_rho > rho_thresh & ~isnan(tmp_dm);
    
    if contains(all_data(i).meta,'empty')
        c = 'k';
    else
        c = 'r';
    end
    nexttile
    scatter(tmp_dr(idx),tmp_dm(idx),c,'.')
    hold on
    plot([-10,10],[-10,10],':k')
    
end
linkaxes(get(get(gcf,'Children'),'Children'))

figure(10); clf
tiledlayout(ceil(sum(first_idx)/2),2)
for i = find(first_idx)'
    nexttile
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),all_data(i).im.z)
    hold on

    if contains(all_data(i).meta,'empty')
        c = 'k';
    else
        c = 'r';
    end
    a = plot(all_data(i).ft.xf,-all_data(i).ft.cue,c,'linewidth',2);
    a.YData(abs(diff(a.YData))>pi) = nan;
    title(all_data(i).meta,'Color','w')
end
if dark_mode
set(gcf,'Color','none','InvertHardCopy','off')
end
figure(11); clf
tiledlayout(ceil(sum(first_idx)/2),2)

[fly_names,~,fly_idx] = unique(cellfun(@(x)(x(1:50)),{all_data.meta},'UniformOutput',false));
num_fly = length(fly_names);

for i = 1:num_fly
    nexttile
    for j = find(fly_idx == i)'
        histogram(gain_pos{j,1}(gain_pos{j,2} < median(gain_pos{j,2})),'binwidth',.1)
        hold on
    end
    plot([1,1],ylim,':k')
    xlim([-1,5])    
    title(all_data(i).meta)
end

%% see how gain changes over trials
[fly_names,~,fly_idx] = unique(cellfun(@(x)(x(1:50)),{all_data.meta},'UniformOutput',false));
num_fly = length(fly_names);

figure(11); clf; hold on
for i = 1:(num_fly)
    plot(empty_idx(fly_idx==i) + .1*(1:sum(fly_idx==i)),gain_vel(fly_idx==i,1)','o-','Color',.8*[1,1,1])
end
xlim([-1,2])
    
plot(empty_idx(fly_idx==i) + .1*(1:sum(fly_idx==i)),gain_vel(fly_idx==i,1)','o-','Color',.8*[1,0,0])


figure(12); clf
rows = ceil(num_fly/2);
rows = 2;
for i = num_fly-3:(num_fly)
    subplot(rows,2,i-(num_fly-3)+1)
    hold on
    for j = find(fly_idx==i)'
        histogram(gain_pos{j}(gain_pos{j}(:,2)<median(gain_pos{j}(:,2)),1),'binwidth',.1)
    end
end
%% work on a positional gain metric
i = 91;
figure(2); clf
subplot(2,1,1)
plot(all_data(i).ft.xb,unwrap(all_data(i).im.mu))
hold on
plot(all_data(i).ft.xf,cumsum(all_data(i).ft.r_speed)/60)
legend('$$\hat{\theta}$$ (bump)','$$\int\dot{\theta}$$ (int vel)','Interpreter','Latex')

subplot(2,1,2)
tmp_win = 10;
tmp_t = floor(all_data(i).ft.xb(1)):ceil(all_data(i).ft.xb(end));
tmp1 = unwrap(all_data(i).im.mu);
tmp1 = interp1(all_data(i).ft.xb,tmp1,tmp_t,'linear','extrap');
tmp2 = cumsum(all_data(i).ft.r_speed)/60;
tmp2 = interp1(all_data(i).ft.xf,tmp2,tmp_t,'linear','extrap');

tmp_gain = nan(length(tmp_t)-tmp_win,2);
tmp_fval = nan(length(tmp_t)-tmp_win,1);

hold on
for t = 1:length(tmp_t)-tmp_win
    fun = @(x)(sum((...
           tmp2(t:t+tmp_win)*x(1)+x(2) - tmp1(t:t+tmp_win)).^2));
    

    [tmp_gain(t,:),tmp_fval(t)] = fminsearch(fun,[1,0]);
    %scatter(t,fun(tmp_gain(t,:)),'.','r')
    scatter(t,tmp_gain(t,1),'.','r')
end

linkaxes(get(gcf,'Children'),'x')

%% plot result
figure(5); clf
subplot(2,2,1)
histogram(gain_pos{1,1}(:,1),'binwidth',.1)
hold on
histogram(gain_pos{2,1}(:,1),'binwidth',.1)
title('LPsP > CsChrimson, Light Dissect','Color','w')
xlim([-1,5])

subplot(2,2,2)
histogram(gain_pos{19,1}(:,1),'binwidth',.1)
hold on
histogram(gain_pos{20,1}(:,1),'binwidth',.1)
title('Empty > CsChrimson, Light Dissect','Color','w')
xlim([-1,5])

subplot(2,2,3)
histogram(gain_pos{46,1}(gain_pos{46,2}<20,1),'binwidth',.1)
hold on
histogram(gain_pos{47,1}(gain_pos{47,2}<20,1),'binwidth',.1)
title('LPsP > CsChrimson, IR Dissect','Color','w')
xlim([-1,5])

subplot(2,2,4)
histogram(gain_pos{25,1}(gain_pos{25,2}<20,1),'binwidth',.1)
hold on
histogram(gain_pos{26,1}(gain_pos{26,2}<20,1),'binwidth',.1)
title('Empty > CsChrimson, IR Dissect','Color','w')
set(get(gcf,'Children'),'xcolor','w','ycolor','w','color','none')
xlim([-1,5])
legend('pre stim','post stim','TextColor','w')

%% show individual movie
i = 16;
pause_time = 1e-1;

%load([all_files(i).folder,'\',all_files(i).name])

tmp_img = imgData_smooth_reg; %smoothdata(imgData_reg,3,'movmean',20);
b1 = min(tmp_img,[],'all');
b2 = prctile(tmp_img(:),99.99); %max(tmp_img,[],'all');
tmp_img = 256*(tmp_img - b1) ./ (b2 - b1);
% 
tmp_img2 = smoothdata(imgData_reg,3,'movmean',5);
b1 = min(tmp_img2,[],'all');
b2 = prctile(tmp_img2(:),99.99); %max(tmp_img,[],'all');
tmp_img2 = 256*(tmp_img2 - b1) ./ (b2 - b1);



stims = interp1(all_data(i).ft.xf,all_data(i).ft.stims,all_data(i).ft.xb);
tmp_cue = interp1(all_data(i).ft.xf,-all_data(i).ft.cue,all_data(i).ft.xb, 'linear','extrap');

figure(13); clf
subplot(2,2,1); h1 = image(tmp_img(:,:,1)); axis equal tight
subplot(2,2,2); h2 = image(tmp_img2(:,:,1)); axis equal tight
subplot(2,2,3); h3 = plot(unwrap(all_data(i).im.alpha),all_data(i).im.z(:,1));ylim([-1,5]);
hold on;        h4 = scatter(0,0,'m','filled');
                h5 = scatter(-pi,0,'r');
for f = 900:size(tmp_img,3)
    h1.CData = tmp_img(:,:,f);
    h2.CData = tmp_img2(:,:,f);
    h3.YData = all_data(i).im.z(:,f);
    h4.XData = tmp_cue(f);
    h5.YData = stims(f);
    xlabel(all_data(i).ft.xb(f))
    pause(pause_time)
end

%% Functions

function s = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type)

    f_speed = ftData_dat.velFor{:};                       %store each speed
    f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}' / 192 * 2 * pi - pi;
    cue(abs(gradient(cue)) > 2) = nan;
    cue     = unwrap(cue);
    cue     = smoothdata(cue,1,ft_type,ft_win,'omitnan');
    cue   = mod(cue,2*pi);                                %rewrap heading data, and put between -pi and pi.
    cue(cue > pi) = cue(cue > pi) - 2*pi;
    
    s.xf      = seconds(ftData_DAQ.trialTime{:});
    if ismember('volClock',ftData_DAQ.Properties.VariableNames);
        s.xb      = seconds(ftData_DAQ.volClock{:});
    end
    s.f_speed = smoothdata(f_speed,1,ft_type,ft_win); 
    s.r_speed = smoothdata(r_speed,1,ft_type,ft_win);
    s.cue     = cue;
    
end


function x_white = whiten(x,dim)
    if nargin < 2
        dim = 1;
    end
    x_white = (x - mean(x,dim)) ./ range(x,dim);
end

function s = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct)
    
    imgData = smoothdata(imgData,3,im_type{1},im_win{1});
    imgData = imgData - min(imgData,[],'all');

    [y_mask,x_mask] = find(mask);                                             %find the cartesian coordinates of points in the mask, just to find the minimum axis length
min_axis        = min(range(x_mask),range(y_mask));
mid             = bwskel(mask,'MinBranchLength',min_axis);  %find the midline as the skeleton, shaving out all sub branches that are smaller than the minimum axis length
[y_mid,x_mid]   = find(mid);                                              %by definition, the skeleton has to be at least as long as the min width of the roi, so shave out subbranches that are shorter than that.
ep              = bwmorph(mid,'endpoints');                               %find the endpoints of the midline
[y0,x0]         = find(ep,1);

[x_mid,y_mid]   = graph_sort(x_mid,y_mid);                                      %align the points of the midline starting at the first pointpoint and going around in a circle. this requires that the midline be continuous!

xq          = [-min_axis:(length(x_mid)+min_axis)];                             %extend the midline so that it reaches the border of the mask. extrapolate as many points as the minimum axis length
x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');                 %keep only the points that exist within the mask
x_mid       = x_mid(idx);
y_mid       = y_mid(idx);

xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';                        %set query points for interpolation (the number of centroids we want). we'll create twice as many points and take every other so that clusters on the edges arent clipped
centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];  %interpolate x and y coordinates, now that they are ordered, into evenly spaced centroids (this allows one to oversample the number of pixels, if desired)
centroids   = centroids(2:2:end-1,:);                                                     %take every other so that we dont start at the edges, and all are same size
%assign each pixel to a centroid
[~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1); %find the index of the centroid that is closest to each pixel in the mask. using euclidean, but maybe chebychev (chessboard)

    imgData_2d      = reshape(imgData,[],size(imgData,3));                  %reshape the data into a 2D pixels with dimensions AllPixels x Frames, where each entry is an intensity
    centroid_log    = false(2*n_centroid,size(imgData_2d,1));               %initialize a logical matrix that is of dimensions Centroids  x AllPixels
    for i = 1:2*n_centroid                                                  %For each centroid, define which pixels are contained in that centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx ==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);     %the summed fluorescence in each group will be [centroids x pixels] * [pixels  x frames], and dividing by the number of pixels in each group gives the average intensity at each frame
    f0              = prctile(f_cluster,f0_pct,2);                          %find the baseline fluorescence in each cluster
    dff_cluster     = (f_cluster - f0) ./ f0;                               %find the dF/F in each cluster. this puts everything on the same scale and eliminates baseline differences.
    zscore_cluster  = zscore(dff_cluster,[],2);
    
    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    alpha       = repmat(alpha,1,2);
    
    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);                                %rewrap heading data, and put between -pi and pi.
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});
    
    % imgData = squeeze(sum(imgData,3));
    % imgData = 256*(imgData-min(imgData,[],'all'))/(max(imgData,[],'all')-min(imgData,[],'all'));

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.d  = dff_cluster;
    s.f  = f_cluster;
    s.alpha = alpha;
    %s.imgData = imgData;
end

function out = resample_pablo(N,x)
    idx = randi(length(x),length(x),N);
    out = x(idx);
end


function [cc,gain,bias] = find_cc(lag,fly_vel,bump_vel,rho,vel_thresh,bump_thresh,rho_thresh)
    fly_vel  = fly_vel(1:end-lag);
    bump_vel = bump_vel(lag+1:end);
    rho      = rho(lag+1:end);

    idx = abs(fly_vel) > vel_thresh & abs(bump_vel) < bump_thresh & rho > rho_thresh;

    cc  = corr(fly_vel(idx),bump_vel(idx));
    b = [ones(sum(idx),1),fly_vel(idx)]\bump_vel(idx); 
    gain = b(2);
    bias = b(1);
end

function [gain,bias] = find_gain(fly_vel,mu,fr)

    int_pos = cumsum(fly_vel * fr);

    tmp1 = mu(~isnan(mu));
    tmp2 = int_pos(~isnan(mu));
    obj_fun = @(x) circ_var(circ_dist(tmp1,tmp2*x(1) + x(2)));
    x0 = [.7,0];
    lb = [-10,-10];
    ub = [10,10];
   
    x = fmincon(obj_fun,x0,[],[],[],[],lb,ub);
    gain = x(1);
    bias = x(2);
end

function h = plot_sem(ax,t,x)

m1 = mean(x,1,'omitnan');
s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

idx = ~isnan(m1);
m1 = m1(idx);
s1 = s1(idx);
t  = t(idx);


h = patch(ax,[t;flipud(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end

function rgb_image = mat2rgb(v,map)
    minv = min(v(:));
    maxv = max(v(:));
    ncol = size(map,1);
    s = round(1+(ncol-1)*(v-minv)/(maxv-minv));
    rgb_image = ind2rgb(s,map);
end