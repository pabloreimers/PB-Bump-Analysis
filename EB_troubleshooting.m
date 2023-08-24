%% ask user for image data
[filename,filepath] = uigetfile('.mat','Select Registered Movie',data_dir);
load([filepath,'\',filename])

%% show EB in 3D
mean_pix = mean(regProduct,4); %take the mean of each pixel over time
dims = size(mean_pix);

x = repmat([1:dims(2)],dims(1),1);
y = repmat([1:dims(1)]',1,dims(2));



figure(1); clf; hold on
for i = 1:dims(3)
    val = reshape(mean_pix(:,:,i),1,[]);
    val = val - prctile(val,5);
    val = val/max(val);

    scatter3(x(:),y(:),i*ones(numel(x),1),[],max(val,0),'filled')
end

%%
figure(20); clf; 
subplot(2,1,1); axis equal; hold on
a = image(imgData(:,:,1)); colormap(bone)
subplot(2,1,2); hold on
b = plot(alpha,dff_cluster(:,1));
c = scatter(0,0,'filled','r');
ylim([-1,3])
xlim([-pi,pi])
for i = 1:length(imgData)
    a.CData = imgData2(:,:,i);
    b.YData = dff_cluster(:,i);
    c.XData = interp1(xf,mu,xb(i));
    drawnow
    input('')
end
