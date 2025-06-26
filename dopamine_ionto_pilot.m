figure(1); clf

for i = 1:20
    tmp = zeros([size(all_data(i).im.d),3]);
    tmp(:,:,1) = all_data(i).atp.d*2;
    tmp(:,:,2) = all_data(i).im.d;
    
    a1 = subplot(5,4,i);
    imagesc(all_data(i).ft.xb,unwrap(all_data(i).im.alpha),tmp); xticks([]); yticks([])

    pos = get(gca,'Position');
    pos(2) = pos(2) - .03;
    pos(4) = .03;
    a2 = axes('Position',pos,'color','none');
    hold on
    plot(all_data(i).ft.xb,sum(all_data(i).atp.d,1)/max(sum(all_data(i).atp.d,1)),'r')
    plot(all_data(i).ft.xb,sum(all_data(i).im.d,1)/max(sum(all_data(i).im.d,1)),'g')
    xticks([]); yticks([])

    linkaxes([a1,a2],'x')
    axis tight
end

%% extract each pulse and align

win = [-10,60]; % in seconds
pulse_tmp = nan(ceil((win(2) - win(1)) / mean(diff(all_data(1).ft.xb)))+10,1);
pulse_cell = {};
length_cell = [];
counter = 1;

for i = 1:length(all_data)
    stim_ind = find(diff(all_data(i).ft.stims)>0);
    stim_stop = find(diff(all_data(i).ft.stims)<0);
    stim_length = unique(round(all_data(i).ft.xf(stim_stop))) - unique(round(all_data(i).ft.xf(stim_ind)));
    [~,stim_ind] = min(abs(all_data(i).ft.xf(stim_ind)' - all_data(i).ft.xb));

    fr = mean(diff(all_data(i).ft.xb));
    stim_start = stim_ind + round(win(1)/fr);
    stim_end = stim_ind + round(win(2)/fr);

    tmp = sum(all_data(i).im.d,1) / max(sum(all_data(i).im.d,1));

    for j = 1:length(stim_ind)
        tmp2 = pulse_tmp;
        tmp3 = tmp(stim_start(j):min(stim_end(j),length(tmp)));
        tmp2(1:length(tmp3)) = tmp3;
        pulse_cell{counter} = tmp2;
        length_cell{counter} = stim_length(j);
        counter = counter+1;
    end
end
pulse_cell = cell2mat(pulse_cell)';
pulse_cell(:,all(isnan(pulse_cell),1)) = [];
length_cell = cell2mat(length_cell)';
figure(2); clf; hold on
t = linspace(win(1),win(2),size(pulse_cell,2));
m = mean(pulse_cell(length_cell==1,:),1,'omitnan');
s = std(pulse_cell(length_cell==1,:),1,'omitnan') ./ sqrt(sum(~isnan(pulse_cell),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],'r')

t = linspace(win(1),win(2),size(pulse_cell,2));
m = mean(pulse_cell(length_cell==2,:),1,'omitnan');
s = std(pulse_cell(length_cell==2,:),1,'omitnan') ./ sqrt(sum(~isnan(pulse_cell),1));
patch([t,fliplr(t)],[m+s,fliplr(m-s)],'g')
%plot(t,m)

