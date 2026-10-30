%% fig_dlight_noise_anatomy
% Breaks one raw epg>dLight trial into its three components:
%   1. scan/line noise  -- a ~13 cyc/line oscillation on the analog baseline,
%                          random in phase from line to line and frame to frame
%   2. shot noise       -- individual photons, each a single-pixel pulse of
%                          ~50-1500 counts (median ~335) that stands far above the baseline
%   3. dLight signal    -- the PB arch, only visible as a rate of photon arrivals
%                          after summing many frames
% Key observation: the no-photon baseline is confined to about +/-50 counts
% (the negative tail stops dead at -50), so the line noise lives entirely
% below a threshold that every photon pulse clears. Counting photons
% (raw > thr) therefore removes the line noise in a single frame, with no
% temporal averaging, which pb_remove_scan_noise (10-frame mean) cannot do.

%% 0. load a chunk of the first trial
addpath(fullfile(fileparts(mfilename('fullpath')), '..', '..', 'claude'));
expDir  = 'Z:\pablo\epg_dlight\20260216\fly 1\20260216-1_epg_dLight';
volumes = [1 400];        % ~38 s at 10.5 vol/s; enough to see the PB
thr     = 50;             % counts; baseline max |value| ~50, photons above
plane   = 5;              % plane shown in the single-frame panels
outDir  = fullfile(fileparts(mfilename('fullpath')), '..', 'exports');

[mov, md] = pb_load_raw_tif(expDir, 'volumes', volumes);
m  = single(mov);                         % Y X Z T
fr = squeeze(m(:,:,plane,:));             % Y X T, one plane
ph = single(m > thr);                     % photon events
mk = load(fullfile(expDir,'mask.mat')).mask;   % hand-drawn PB mask from earlier processing

% whole-trial per-plane means (~50 s streaming read; Tiff class would stop at 65535 frames)
wt = pb_read_raw_tif_planes(md.tifPath, md.nplanes, {}, 'planeMeans', true, 'verbose', false);
inout = arrayfun(@(z) mean(wt.planeMean(:,:,z) .* mk, 'all') / nnz(mk) / ...
    (mean(wt.planeMean(:,:,z) .* ~mk, 'all') / nnz(~mk)), md.valid_planes);

%% 1. statistics
vals = double(fr(:));
F  = abs(fft(fr - mean(fr,2), [], 2)).^2;          S  = squeeze(mean(F,[1 3]));
frT = squeeze(ph(:,:,plane,:));
FT = abs(fft(frT - mean(frT,2), [], 2)).^2;        ST = squeeze(mean(FT,[1 3]));
mfr = mean(fr,3);
FM = abs(fft(mfr - mean(mfr,2), [], 2)).^2;        SM = mean(FM,1);
k  = 0:128;
mu = mean(fr,3); va = var(fr,0,3); pfit = polyfit(mu(:), va(:), 1);

%% 2. figure
t0 = 200;                                  % example volume
fig = figure('Color','w','Position',[50 50 1500 900],'Visible','on');
tl = tiledlayout(fig,3,4,'TileSpacing','compact','Padding','compact');
title(tl, sprintf('%s  trial 1 -- raw noise anatomy (plane %d, vols %d-%d)', ...
    '20260216-1 epg>dLight', plane, volumes), 'Interpreter','none');

nexttile; imagesc(fr(:,:,t0)); axis image off; colormap(gca,gray); clim([-50 600]); colorbar
title(sprintf('raw single frame (vol %d)', t0))

nexttile; imagesc(fr(:,:,t0)); axis image off; colormap(gca,gray); clim([-40 40]); colorbar
title('same frame, clim \pm40: line noise')

nexttile; r = double(fr(40,:,t0)); plot(r,'k'); hold on
r2 = r; r2(r2>thr) = NaN; plot(r2,'Color',[0.85 0.33 0.1],'LineWidth',1.5);
yline(thr,'b--'); xlim([1 256]); ylim([-60 700]); box off
xlabel('x (fast-scan px)'); ylabel('counts'); title('one line: photons vs baseline')
legend({'raw','baseline (\le thr)','thr'},'Box','off','Location','northeast')

nrm = @(s) s(1:129) / median(s(40:100));   % flat high-freq floor = 1
nexttile; semilogy(k, nrm(S), 'k'); hold on
semilogy(k, nrm(ST), 'Color',[0 0.45 0.74]);
semilogy(k, nrm(SM), 'Color',[0.47 0.67 0.19]);
xline(13,':'); xlim([1 128]); box off
xlabel('cycles / line'); ylabel('power (norm.)'); title('fast-axis spectrum')
legend({'raw frame','photon-thresholded frame','400-vol mean'},'Box','off')

nexttile; edges = -150:10:2000;
histogram(vals, edges, 'EdgeColor','none','FaceColor','k'); set(gca,'YScale','log'); hold on
histogram(-vals(vals<0), edges, 'DisplayStyle','stairs','EdgeColor',[0.85 0.33 0.1]);
xline(thr,'b--'); box off; xlabel('counts'); ylabel('pixels')
title('pixel histogram'); legend({'raw','mirrored negative side','thr'},'Box','off')

nexttile; imagesc(frT(:,:,t0)); axis image off; colormap(gca,flipud(gray)); colorbar
title(sprintf('photon events, one frame (n=%d)', nnz(frT(:,:,t0))))

nexttile; imagesc(sum(ph(:,:,:,t0),3)); axis image off; colormap(gca,flipud(gray)); colorbar
title('photons, one volume z-summed (~95 ms)')

nexttile; plot(mu(:), va(:), '.', 'Color',[0 0 0 0.1]); hold on
xx = [0 max(mu(:))]; plot(xx, polyval(pfit,xx), 'r'); box off
xlabel('pixel mean (counts)'); ylabel('pixel variance')
title(sprintf('shot noise: var = %.0f x mean', pfit(1)))

nexttile; imagesc(imgaussfilt(mean(sum(ph,3),4),1.5)); axis image off; colormap(gca,hot); colorbar
hold on; contour(mk,[.5 .5],'c'); title('photons z-sum, 400-vol mean (\sigma 1.5)')

nexttile; imagesc(imgaussfilt(sum(wt.planeMean(:,:,md.valid_planes),3),1.5)); axis image off
colormap(gca,hot); colorbar; hold on; contour(mk,[.5 .5],'c')
title(sprintf('raw z-sum, whole trial (%d vols, unregistered)', wt.nVolumes))

nexttile; bar(inout,'k'); yline(1,':'); box off; ylim([0.8 1.5])
xlabel('plane'); ylabel('mean in mask / out'); title('PB contrast, whole trial')

nexttile; pm = squeeze(mean(ph,[1 2 4])); bar(pm*md.px_height*md.px_width,'k'); box off
xlabel('plane'); ylabel('photons / frame'); title('photon rate by plane')

%% 3. export
if ~exist(outDir,'dir'), mkdir(outDir); end
exportgraphics(fig, fullfile(outDir,'fig_dlight_noise_anatomy.png'), 'Resolution', 130);
