%% load in data
[filename,filepath] = uigetfile('.mat','Select Imaging Data');
[filename2,filepath2] = uigetfile(filepath,'Select FicTrac Data');
load([filepath,'\',filename])
load([filepath2,'\',filename2])

%% Process and find correlations
imgData = squeeze(sum(regProduct,3));

velYaw = smoothdata(ftData.yawSpeed{:},1,'gaussian',30);
%velYaw = smoothdata(ftData_DAQ.velYaw{:},1,'gaussian',30);

xf = linspace(0,10,size(velYaw,1));
xb = linspace(0,10,size(imgData,3));

all_pix = reshape(imgData,[],size(imgData,3));
all_pix = interp1(xb,all_pix',xf)';

right_rho = reshape(corr(all_pix',max(velYaw,0)),size(imgData,1),size(imgData,2));
left_rho = reshape(corr(all_pix',max(-velYaw,0)),size(imgData,1),size(imgData,2));

%% Plot 
im_gain = 30;
tot_rho = (right_rho .* reshape([1,.5,0],1,1,3) * im_gain) + ...
          (left_rho .* reshape([0,.5,1],1,1,3) * im_gain);

figure(1); clf
subplot(2,2,1)
image((right_rho .* reshape([1,.5,0],1,1,3) * im_gain))
title('right correlations','color','w')

subplot(2,2,2)
image(left_rho .* reshape([0,.5,1],1,1,3) * im_gain)
title('left correlations','color','w')

subplot(2,2,3)
image(tot_rho)
title('overlay','color','w')

subplot(2,4,7)
plot(xf,velYaw)
xlabel('time (min)','color','w')
ylabel('Rotational Velocity (rad/s)','color','w')
set(gca,'color','none','ycolor','w','xcolor','w')

subplot(2,4,8)
barh([sum(velYaw<0),sum(velYaw>0)]/60)
xlabel('time','color','w'); yticklabels({'left turns','right turns'})
set(gca,'color','none','ycolor','w','xcolor','w')

set(gcf,'color','none')
annotation(gcf,'textbox',[.4,.4,.1,.1],'String',filepath,'color','w')
