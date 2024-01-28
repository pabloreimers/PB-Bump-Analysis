stopband = [10 20]; %set emperically for now, stopband frequency indices keep between 2 and half x length . . . hopefully scan noise is fairly constant across recordings

[filename1,filepath1] = uigetfile(".mat",'Select Raw Movie');
[filename2,filepath2] = uigetfile(filepath1,'Select Denoised Movie');

load([filepath1,filename1]);
rp_raw = regProduct;

load([filepath2,filename2]);
rp_dn = regProduct;

%% smooth both

rp_raw_smooth = smoothdata(rp_raw,4,'gaussian',30);
rp_dn_smooth  = smoothdata(rp_dn,4,'gaussian',30);

%% remove scan noise

rp_raw_rsn = cx_fft_filter_1d(rp_raw_smooth,stopband);
rp_dn_rsn  = cx_fft_filter_1d(rp_dn_smooth,stopband);

%% build a remove scna noise function
rm_bands = 5;

im_2d = reshape(permute(rp_dn,[2,1,3,4]),size(rp_dn,2),[]);
n = size(im_2d,1);
b = mean(im_2d,1);
x = im_2d - b;
fhat = fft(x,[],1);
fhat(rm_bands:(n-rm_bands+2),:) = 0;
im_rsn = ifft(fhat,[],1) + b;

rp_dn_rsn = permute(reshape(im_rsn,size(rp_dn,2),size(rp_dn,1),size(rp_dn,3),size(rp_dn,4)),[2,1,3,4]);

%% show the new denoised image

figure(3); clf
for i = 1:size(rp_dn_rsn,4)
    subplot(1,2,1)
    imagesc(sum(rp_dn(:,:,:,i),3))
    subplot(1,2,2)
    imagesc(sum(rp_dn_rsn(:,:,:,i),3))
    pause(1e-1)
end


%% plot all
figure(1); clf

im_cell = {rp_raw,rp_dn;...
           rp_raw_smooth,rp_dn_smooth;...
           rp_raw_rsn,rp_dn_rsn};

for i = 1:6
    h(i) = imagesc(subplot(3,2,i),sum(im_cell{i}(:,:,:,1),3));
end

for j = 1:size(rp_raw,4)
    for i = 1:6
    h(i).CData = sum(im_cell{i}(:,:,:,j),3);
    end
        drawnow
    pause(1e-2)
end


function imout = cx_fft_filter_1d(imin, stopband)

szperm = [size(imin, 2), size(imin, 1), size(imin, 3), size(imin,4)];
imin = permute(imin, [2 1 3 4]);
imin = reshape(imin, size(imin, 1), []);  %collapse z and t because we believe dominant structure is through true time (not volume time)


Ts = 1;                                                             % Sampling Interval
Fs = 1/Ts;                                                          % Sampling Frequency
Fn = Fs/2;                                                          % Nyquist Frequency
L = size(imin,1);                                                   % Length Of ‘data’ Vector
t = 1:L*Ts;
%t = linspace(0, 1, L)*Ts;% Time Vector

Fv = linspace(0, 1, fix(L/2)+1)*Fn;                                 % Frequency Vector (One-Sided FFT)
Iv = 1:length(Fv);                                                 % Index Vector
Ivkf = 1:stopband(1);                                                      % 8First keepfreq FFT Frequencies
Fkfth = Fv(stopband(1));                                                     % Frequency Corresponding To keepfreqth Element
Fkfth2 = Fv(stopband(2));                                                     % Frequency Corresponding To keepfreqth Element


FLen = 40;                                                          % Discrete Filter Order
b_filt = fir1(FLen, [Fkfth/Fn Fkfth2/Fn], 'stop', chebwin(FLen+1,30));                                       % Design FIR Filter


imout = zeros(size(imin), 'single');

for indi = 1:size(imin, 2) %loop over lines, filtering

    data = double(imin(:,indi));
    mnd = mean(data);
    data = data - mnd;
    tmpout = fftfilt(b_filt, data);
    imout(:,indi) = tmpout + mnd;

end


%put back in 4d
imout = reshape(imout, szperm);
imout = permute(imout, [2 1 3 4]);
if ~isa(imout, 'single')
    imout = single(imout);
end

datmin_raw = min(imout(:));
datmax_raw = max(imout(:));
imout = imout - datmin_raw;
% if datmax_raw > 2^16-1
%     "ERROR, CLIPPING REQUIRED, CHANGE OUTPUT TYPE"
%     error
% end
% imout = uint16(imout);

end

function regProduct_rsn = rsn(regProduct)
    
end