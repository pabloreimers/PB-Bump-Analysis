%% extract pulses
win_start = -2;
win_end = 20;

c_pulses = {};
m_pulses = {};
a_pulses = {};
d_pulses = {};
z_pulses = {};
min_pulses = {};
max_pulses = {};
mean_pulses = {};
dark_idx = {};
exp_idx  = {};
right_idx= {};
fly_num_cell = {};

tmp = cellfun(@(x)(x(1:47)),{all_data.meta}','UniformOutput',false);
[~,~,fly_num] = unique(tmp);

counter = 0;
for i = 1:length(all_data)
    mu  = interp1(all_data(i).ft.xb, unwrap(all_data(i).im.mu),all_data(i).ft.xf); %interpolate imaging data to fictrac rate
    atp = interp1(all_data(i).ft.xb, all_data(i).atp.d',all_data(i).ft.xf)';
    dff = interp1(all_data(i).ft.xb, all_data(i).im.d',all_data(i).ft.xf)'; %interpolate imaging data to fictrac rate

    loc = find(diff(all_data(i).ft.stims)>5); %find when the stim happened
    
    sum_atp = sum(all_data(i).atp.f,2);
    for t = 1:length(loc)
        %if no atp ejection, don't count it
        if (mean(atp(:,loc(j)+120)) - mean(atp(:,loc(j))))/(mean(atp(:,loc(j)))+1) < 2
            continue
        end
        tmp_start = loc(t) + win_start*60;
        tmp_end   = loc(t) + win_end*60;
        counter = counter+1;
        
        c_pulses{counter} = all_data(i).ft.cue(tmp_start:tmp_end);
        m_pulses{counter} = mu(tmp_start:tmp_end);
        d_pulses{counter} = dff(:,tmp_start:tmp_end);
        min_pulses{counter} = min(dff(:,tmp_start:tmp_end),[],1,'omitnan')';
        max_pulses{counter} = max(dff(:,tmp_start:tmp_end),[],1,'omitnan')';
        mean_pulses{counter} = mean(dff(:,tmp_start:tmp_end),1,'omitnan')';
        right_idx{counter} = mean(sum_atp(1:end/2)) > mean(sum_atp(end/2:end));
        exp_idx{counter} = ~contains(all_data(i).meta,'empty');
        dark_idx{counter} = contains(all_data(i).ft.pattern,'background');
        fly_num_cell{counter} = fly_num(i);
    end
end

dark_idx = logical(cell2mat(dark_idx));
exp_idx = logical(cell2mat(exp_idx));
right_idx = logical(cell2mat(right_idx));
fly_num_cell = cell2mat(fly_num_cell);


%% plot sweeps
tmp_t = win_start:1/60:win_end;

tmp_idx   = {right_idx & ~dark_idx, ~right_idx & ~dark_idx} ;

pulse_mat = cell2mat(cellfun(@(x)(unwrap(x-x(1))),m_pulses,'UniformOutput',false))';

figure(2); clf

for i = 1:2
subplot(1,2,i)
hold on
a = plot_sem(gca,tmp_t,pulse_mat(tmp_idx{i} & exp_idx,:)); a.FaceColor = 'r'; a.HandleVisibility = 'off';
plot(tmp_t,mean(pulse_mat(tmp_idx{i} & exp_idx,:),1,'omitnan'),'r','linewidth',2)
a = plot_sem(gca,tmp_t,pulse_mat(tmp_idx{i} & ~exp_idx,:)); a.FaceColor = 'b'; a.HandleVisibility = 'off';
plot(tmp_t,mean(pulse_mat(tmp_idx{i} & ~exp_idx,:),1,'omitnan'),'b','linewidth',2)

plot([win_start,win_end],[0,0],':k')
scatter(0,0,100,'r*')
%title('Minimum Fluorescence')
xlabel('time post stim (s)')
legend(sprintf(['Right stim',' (N = %i)'],length(unique(fly_num_cell(tmp_idx{i} & right_idx)))),...
       sprintf(['Left stim',' (N = %i)'],length(unique(fly_num_cell(tmp_idx{i} & ~right_idx)))),...
        'Location','Northeast','textcolor','k','edgecolor','k','Color','none')
axis tight
set(gca, 'YDir','reverse','xcolor','k','ycolor','k','Color','none')
ylim([-pi,pi]); yticks([-pi,0,pi]); yticklabels({'-\pi',0,'\pi'})
end

ylabel(subplot(1,2,1),{'unwrapped', 'bump', 'position', '(\pi rad)'},'Rotation',0)

%% plot fluorescence

tmp_t = win_start:1/60:win_end;

tmp_idx   = ~dark_idx ;
tmp_label = 'Empty > P2X2';


figure(2); clf
pulse_cell = {min_pulses,max_pulses,mean_pulses};
labels_cell = {'Min','Max','Mean'};
for i = 1:3
pulse_mat = cell2mat(cellfun(@(x)(unwrap(x-x(1))),pulse_cell{i},'UniformOutput',false))';

subplot(2,2,i)
hold on
a = plot_sem(gca,tmp_t,pulse_mat(tmp_idx & exp_idx,:)); a.FaceColor = 'r'; a.HandleVisibility = 'off';
plot(tmp_t,mean(pulse_mat(tmp_idx & exp_idx,:),1,'omitnan'),'r','linewidth',2)
a = plot_sem(gca,tmp_t,pulse_mat(tmp_idx & ~exp_idx,:)); a.FaceColor = 'b'; a.HandleVisibility = 'off';
plot(tmp_t,mean(pulse_mat(tmp_idx & ~exp_idx,:),1,'omitnan'),'b','linewidth',2)
plot(tmp_t,0*tmp_t,'k:')

set(gca,'Color','none')
xlabel('time (s)')
ylabel([labels_cell{i},' Fluoresnce'])
legend('LPsP > P2X2','Empty > P2X2','Color','none')
end

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