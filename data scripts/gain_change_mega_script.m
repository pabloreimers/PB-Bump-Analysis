%% create meta data indexes
trial_num = zeros(length(all_data),1);
dark_idx  = false(length(all_data),1);
empty_idx = false(length(all_data),1);
walk_idx  = false(length(all_data),1);
gain_vr  = nan(length(all_data),1);

last_str = '';
for i = 1:length(all_data)
    tmp_str = all_data(i).meta(1:45);


    if ~strcmp(tmp_str,last_str)
        counter = 0;
        last_str = tmp_str;
    end
    counter = counter+1;
    trial_num(i) = counter;
    last_str = tmp_str;
    
    if sum(all_data(i).ft.f_speed>0) > length(all_data(i).ft.f_speed)/2
        walk_idx(i) = true;
    end
    
    if contains(all_data(i).ft.pattern,'background')
        dark_idx(i) = true;
    end

    if contains(all_data(i).meta,'empty')
        empty_idx(i) = true;
    end

    gain_vr(i) = median(all_data(i).gain.vr,'omitnan');
end

fly_num = [1;cumsum(diff(trial_num)<-1)+2];

%% scatter vr gain against decoded gain

gain_vr = nan(length(all_data),1);
gain_mu = nan(length(all_data),1);

for i = 1:length(all_data)
    gain_vr(i) = round(median(all_data(i).gain.vr,'omitnan'),1);
    gain_mu(i) = median(all_data(i).gain.g(all_data(i).gain.v<.5),'omitnan');
end

%% find flies with really good closed loop trials
good_idx = false(length(all_data),1);
tmp1 = nan(length(all_data),1);
tmp2 = nan(length(all_data),1);


for i = 1:length(all_data)
    tmp1(i) = var(all_data(i).gain.g,'omitnan');
    tmp2(i) =  abs(gain_vr(i)-gain_mu(i));
end

good_idx = tmp1 < .8 & tmp2<.15;

figure(2); clf
nexttile; hold on
scatter(gain_vr(~dark_idx & good_idx),gain_mu(~dark_idx & good_idx),'c','filled')
plot(xlim,xlim,':k')
title('selected good cl trials')
good_flies = unique(fly_num(good_idx));
xlabel('vr gain'); ylabel('compass gain')

%
nexttile; hold on
tmp_idx = ismember(fly_num,good_flies) & trial_num<2 | good_idx;
scatter(gain_vr(~dark_idx & tmp_idx),gain_mu(~dark_idx & tmp_idx),'c','filled')
scatter(gain_vr(dark_idx & tmp_idx),gain_mu(dark_idx & tmp_idx),'m')
plot(xlim,xlim,'k:')
legend('closed loop','dark')
xlabel('vr gain'); ylabel('compass gain')
title('all trials from "good" flies')
