[filename,path] = uigetfile('*.tif');

%%
t = Tiff([path,filename]);
samps = read(t);
samps = double(samps(:));

sig = std(max(samps,0));
mu  = mean(max(samps,0));

figure
histogram(samps,'EdgeColor','none'); set(gca,'YScale','log')
ylabel('counts'); xlabel('pixel value')
text(max(xlim),max(ylim),...
    sprintf('\\mu: %.2f\n\\sigma: %.2f\nQE: %.2f',mu,sig,sig/mu),...
    'HorizontalAlignment','right','VerticalAlignment','top')
 
%%
