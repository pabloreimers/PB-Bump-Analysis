clear all
close all

import ScanImageTiffReader.ScanImageTiffReader
%file = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\20230420-1_EPG_7f_LPsP_kir_cl\20230420-1_EPG_7f_LPsP_kir_104117_trial_001_00001.tif';
[filename,path] = uigetfile('C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\')
file = [path,filename];
%%
import ScanImageTiffReader.ScanImageTiffReader
        
% Get image data
imgSI = ScanImageTiffReader(file);

% scanimage object SI
siMetadata = scanimageMetadata(file);
siMetadata = parse_scanimage_metadata(siMetadata);
SI = siMetadata.SI;

%%
imgData = permute(imgSI.data(),[2,1,3]);

imgData = reshape(imgData,size(imgData,1),...
                          size(imgData,2),...
                          SI.hStackManager.numFramesPerVolumeWithFlyback,...
                          []);
imgData = imgData(:,:,1:SI.hStackManager.numFramesPerVolume,:);


%% try different kinds of smoothing
clf

subplot(2,2,1)
histogram(imgData(:))
xlabel('nothing')

subplot(2,2,2)
tmp = smoothdata(imgData,4,'gaussian',5);
histogram(tmp(:))
xlabel('temporal smooth')

subplot(2,2,3)
tmp = imgaussfilt(imgData,5);
histogram(tmp(:))
xlabel('neighbor smooth')

subplot(2,2,4)
tmp = smoothdata(imgData,4,'gaussian',5);
tmp = imgaussfilt(imgData,1);
histogram(tmp(:))
xlabel('temporal then neighbor')

set(get(gcf,'Children'),'YScale','log')
%% play the movie
figure(2); clf
%tmp = smoothdata(imgData,4,'gaussian',5);
tmp = imgData;
tmp = smoothdata(tmp,4,'gaussian',2);
tmp = 256* tmp / prctile(tmp(:),99);
% tmp(tmp > prctile(tmp(:),99)) = 0;
% tmp(tmp < prctile(tmp(:),1)) = 0;
% tmp = smoothdata(tmp,4,'gaussian',2);


a = cell(size(imgData,3),1);
rows = ceil(sqrt(length(a)));
cols = ceil(length(a)/rows);
for i = 1:length(a)
    subplot(rows,cols,i)
    a{i} = image(tmp(:,:,i,1));
end

for i = 1:size(imgData,4)
    for j = 1:size(imgData,3)
        a{j}.CData = tmp(:,:,j,i);
        drawnow
    end
end