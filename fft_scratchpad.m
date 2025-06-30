%%
bands = (1:50); %define which bands to define as noise, in pixel length


%% try it out on single line

figure(1); clf

imgData = img{1}(:,:,3,:); %extract the single line signal to work with


f = imgData(3,:,1)';

subplot(3,1,1)
plot(f)

n = length(f);
fhat = fft(f,n);
PSD = fhat.*conj(fhat)/n;
freq = (0:n);
L = 1:floor(n/2);

subplot(3,1,2)
plot(freq(L),PSD(L))

fhat(1:bands) = 0;
ffilt = ifft(fhat);

subplot(3,1,3)
plot(f); hold on
plot(abs(ffilt))


%% expand out to image
figure(2); clf

imgData = img{1}(:,:,6,:); %extract the single line signal to work with


f = double(imgData(:,:,1)');
subplot(3,1,1)
imagesc(f')

n = size(f,1);
fhat = fft(f,n);
PSD = fhat.*conj(fhat)/n;
freq = (0:n);
L = 1:floor(n/2);

fhat(1:bands,:) = 0;
ffilt = ifft(fhat,n);
subplot(3,1,2)
imagesc(abs(ffilt)')

subplot(3,1,3)
diff_f = f -abs(ffilt);
plot(diff_f)

%% save a denoised stack
dims = size(img{1});
f = img{1};
f_dn = nan(size(imgData));

fhat = fft(f,dims(2),2);
fhat(:,bands,:,:) = 0;
ffilt = abs(ifft(fhat,dims(2),2));

ffilt = reshape(ffilt,dims(1),dims(2),[]);
for i = 1:length(ffilt)
    ffilt(:,:,i) = medfilt2(ffilt(:,:,i));
end
ffilt = reshape(ffilt,dims);