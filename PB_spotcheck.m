load('Z:\pablo\lpsp_p2x2_walking\20260727\fly 3\20260727-17_epg_8m_empty_p2x2\registration\imgData_smooth_reg.mat')
load('Z:\pablo\lpsp_p2x2_walking\20260727\fly 3\20260727-17_epg_8m_empty_p2x2\mask_reg.mat')
load('C:\Users\ReimersPabloAlejandr\Documents\GitHub\LPsP_2p\MelData\PB-Bump-Analysis\.data\lpsp_p2x2_walking_20260728_post0708.mat')
fly_num = 118;

%% interactive version (scrub with slider / arrow keys)
addpath(fullfile(fileparts(mfilename('fullpath')),'claude'))
ft_vid_file = 'Z:\pablo\lpsp_p2x2_walking\20260727\fly 3\20260727-17_epg_8m_empty_p2x2\FicTracData\fictrac-raw-20260727_162955_trial_001.avi';
pb_spotcheck_gui(imgData_smooth_reg, all_data(fly_num), 'VideoFile', ft_vid_file);

%% pick any fly/trial from all_data via dropdowns (loads movie + FicTrac video itself)
addpath(fullfile(fileparts(mfilename('fullpath')),'claude'))
pb_spotcheck_browser(all_data);

%% original animation

figure(1); clf
subplot(2,2,1)
a1 = image(imgData_smooth_reg(:,:,1));
tmp = 200*(imgData_smooth_reg / prctile(imgData_smooth_reg,99,'all'));

subplot(2,1,2)
a2 = imagesc(all_data(fly_num).ft.xb,unwrap(all_data(fly_num).im.alpha),nan(size(all_data(fly_num).im.z)));
hold on
a5 = plot(nan(size(all_data(fly_num).ft.xf)),all_data(fly_num).ft.cue,'c','Linewidth',2);
a5.YData(abs(diff(a5.YData))>pi) = nan;

subplot(2,2,2)
a3 = plot(all_data(fly_num).im.alpha,all_data(fly_num).im.z(:,1));
a3.XData(end/2) = nan;
hold on
a4 = plot(all_data(fly_num).im.mu(1),0,'ro');
ylim([-2,3])
xlim([-pi,pi])

for i = 1:length(imgData_smooth_reg)
    a1.CData = tmp(:,:,i);
    a2.CData(:,i) = all_data(fly_num).im.z(:,i);
    a3.YData = all_data(fly_num).im.z(:,i);
    a4.XData = all_data(fly_num).im.mu(i);
    a5.XData(all_data(fly_num).ft.xf < all_data(fly_num).ft.xb(i)) = all_data(fly_num).ft.xf(all_data(fly_num).ft.xf < all_data(fly_num).ft.xb(i));
    pause(mean(diff(all_data(fly_num).ft.xb))/3);
    drawnow
end