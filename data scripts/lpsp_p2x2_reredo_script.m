%% extract pulses
win_start = 0;
win_end = 15;

c_pulses = {};
m_pulses = {};
a_pulses = {};
d_pulses = {};
z_pulses = {};
f_pulses = {};
dark_idx = {};
exp_idx  = {};
right_idx= {};

counter = 0;
for i = 1:length(all_data)
    mu  = interp1(all_data(i).ft.xb, unwrap(all_data(i).im.mu),all_data(i).ft.xf); %interpolate imaging data to fictrac rate
    atp = interp1(all_data(i).ft.xb, all_data(i).atp.d',all_data(i).ft.xf)';
    dff = interp1(all_data(i).ft.xb, all_data(i).im.d',all_data(i).ft.xf)'; %interpolate imaging data to fictrac rate

    loc = find(diff(all_data(i).ft.stims)>5); %find when the stim happened
    
    sum_atp = sum(all_data(i).atp.f,2);
    for t = 1:length(loc)
        tmp_start = loc(t) + win_start*60;
        tmp_end   = loc(t) + win_end*60;
        counter = counter+1;
        
        c_pulses{counter} = all_data(i).ft.cue(tmp_start:tmp_end);
        m_pulses{counter} = mu(tmp_start:tmp_end);
        d_pulses{counter} = dff(:,tmp_start:tmp_end);
        right_idx{counter} = mean(sum_atp(1:end/2)) > mean(sum_atp(end/2:end));
        exp_idx{counter} = ~contains(all_data(i).meta,'empty');
        dark_idx{counter} = contains(all_data(i).ft.pattern,'background');
    end
end

dark_idx = logical(cell2mat(dark_idx));
exp_idx = logical(cell2mat(exp_idx));
right_idx = logical(cell2mat(right_idx));

%% extract mu aligned pulses
win_start = -5;
win_end = 30;

c_pulses = {};
m_pulses = {};
a_pulses = {};
d_pulses = {};
z_pulses = {};
f_pulses = {};
dark_idx = {};
exp_idx  = {};
right_idx= {};
first_idx= {};
move_idx = {};
stim_length= {};
fr = 10;
tmp_x   = linspace(0,630,630*fr);    
tmp_win = floor(win_start*fr):ceil(win_end*fr); %this is the additive index to the frames to extract for a given pulse window
tmp_t   = tmp_win/10;

counter = 0;
for i = 1:length(all_data)
    tmp_f   = interp1(all_data(i).ft.xb,all_data(i).im.f',tmp_x)';%interpolate everything into the same framerate
    tmp_d   = interp1(all_data(i).ft.xb,all_data(i).im.d',tmp_x)';
    tmp_z   = interp1(all_data(i).ft.xb,all_data(i).im.z',tmp_x)';
    tmp_a   = interp1(all_data(i).ft.xb,all_data(i).atp.d',tmp_x)';
    
    tmp_m   = interp1(all_data(i).ft.xb,unwrap(all_data(i).im.mu),tmp_x)'; tmp_m = mod(tmp_m,2*pi); tmp_m(tmp_m>pi) = tmp_m(tmp_m>pi) - 2*pi;
    %if mean(diff(unwrap(tmp_m)),'omitnan') > 0; tmp_m = -tmp_m; end
    tmp_cue = interp1(all_data(i).ft.xf,unwrap(all_data(i).ft.cue),tmp_x)';tmp_cue = mod(tmp_cue,2*pi); tmp_cue(tmp_cue>pi) = tmp_cue(tmp_cue>pi) - 2*pi; 
   
    tmp_atp = sum(tmp_a,2,'omitnan');
    tmp_right = sum(tmp_atp(1:end/2,:),'all') > sum(tmp_atp(end/2:end,:),'all');
    [~,loc] = findpeaks(smoothdata(max(tmp_a,[],1),'movmean',5),'MinPeakProminence',1,'MinPeakDistance',50*fr);
    %fr = mean(diff(all_data(i).ft.xb)); %find the new framerate (amount of time per frame)
    %tmp_win = floor(win_start/fr):ceil(win_end/fr); %this is the additive index to the frames to extract for a given pulse window
    %findpeaks(smoothdata(max(tmp_a,[],1),'movmean',5),'MinPeakProminence',1.5,'MinPeakDistance',50/fr);
 
    loc = find(diff(all_data(i).ft.stims)>5)';
    if i==1 || ~strcmp(all_data(i).meta(1:40),last_fly)
        last_fly = all_data(i).meta(1:40);
        first_log = true;
    end

    if numel(loc) < 2
        continue
    end


    stim_length_curr = diff(all_data(i).ft.xf(find(abs(diff(all_data(i).ft.stims))>1,2)));
    for j = loc

        n = length(all_data(i).im.alpha);
        
        tmp_idx = j+tmp_win; %extract all of the appropriate frames. if we exceed the window in either direction, just fill it with the nearest calculated gain. a bit hacky
        
        
        pad1 = nan(n,sum(tmp_idx < 1)); %create padding for images
        pad2 = nan(n,sum(tmp_idx > length(tmp_x)));
        tmp_idx = tmp_idx(tmp_idx > 1 & tmp_idx < length(tmp_x)+1);

        counter = counter+1;
        d_pulses{counter} = [pad1,tmp_d(:,tmp_idx),pad2];
        z_pulses{counter} = [pad1,tmp_z(:,tmp_idx),pad2];       
        f_pulses{counter} = [pad1,tmp_f(:,tmp_idx),pad2];       
        a_pulses{counter} = [pad1,tmp_a(:,tmp_idx),pad2];
        m_pulses{counter} = [pad1(1,:),tmp_m(tmp_idx)',pad2(1,:)];
        c_pulses{counter} = [pad1(1,:),tmp_cue(tmp_idx)',pad2(1,:)];
        dark_idx{counter} = contains(all_data(i).meta,'dark') || contains(all_data(i).ft.pattern,'background');
        exp_idx{counter} = ~contains(all_data(i).meta,'empty');
        right_idx{counter} = tmp_right;
        first_idx{counter} = first_log;
        stim_length{counter} = stim_length_curr;
        
        if first_log
            first_log = false;
        end

    end
end

dark_idx = logical(cell2mat(dark_idx));
exp_idx = logical(cell2mat(exp_idx));
right_idx = logical(cell2mat(right_idx));
first_idx = logical(cell2mat(first_idx));
stim_length = cell2mat(stim_length);

alpha = all_data(i).im.alpha;
[~,ind] = min(abs(tmp_win*fr));
ind = round(ind- 2/fr):ind;

%% plot sweeps
tmp_t = win_start:1/60:win_end;

figure(2); clf
%subplot(1,2,1); hold on
hold on
a = plot_sem(gca,tmp_t,cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(~exp_idx & right_idx & ~dark_idx),'UniformOutput',false))'); a.FaceColor = 'r';
a = plot_sem(gca,tmp_t,cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses(~exp_idx & ~right_idx & ~dark_idx),'UniformOutput',false))'); a.FaceColor = 'b';
%a = plot_sem(gca,tmp_t,-cell2mat(cellfun(@(x)(unwrap(x-x(1))),c_pulses(exp_idx),'UniformOutput',false))'); a.FaceColor = 'g';

plot([win_start,win_end],[0,0],':k')
scatter(0,0,100,'r*')
%title('Right Stim','color','w')
%xlabel('time post stim (s)')
%ylabel({'unwrapped', 'bump', 'position', '(\pi rad)'},'Rotation',0)
pos = get(gca,'Position');
legend('LPsP > P2X2 (n=9)','Empty > P2X2 (n=10)','Location','Northeastoutside','textcolor','k','edgecolor','k','Color','none')
axis tight
set(gca, 'YDir','reverse','xcolor','k','ycolor','k','Color','none')
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'})
fontsize(gcf,30,'pixels')


%% functions
function h = plot_sem(ax,t,x)

m1 = mean(x,1,'omitnan');
s1 = std(x,1,'omitnan')./sqrt(sum(~isnan(x),1));

idx = ~isnan(m1);
m1 = m1(idx);
s1 = s1(idx);
t  = t(idx);


h = patch(ax,[t,fliplr(t)],[m1+s1,fliplr(m1-s1)],'r','FaceAlpha',.5);
end