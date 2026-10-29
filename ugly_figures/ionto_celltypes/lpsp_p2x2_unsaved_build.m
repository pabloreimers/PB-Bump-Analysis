%% lpsp_p2x2_unsaved_build
% Builds an all_data struct for the two lpsp_p2x2 experiments that were
% recorded and registered on Z: but never saved to any .mat in .data:
%   Z:\pablo\lpsp_p2x2_redo   (Jan-Feb 2026, EPG syt8s, lpsp/empty>P2X2, genotype in trial name)
%   Z:\pablo\lpsp_p2x2_syt8s  (Jun-Jul 2026, EPG syt8s, blinded; empties tagged "fly N _ empty")
% Found by scanning every Z:\pablo\lpsp_p2x2* trial folder against every
% lpsp_p2x2*.mat (2026-10-29): 195 trials, 0 with a saved copy.
%
% Every trial has registration\imagingData.mat with img = {EPG, red} as
% Y x X x T doubles (no Z dim -- same as dopamine_ionto_walking_claude.m
% header point 2), a trial-level mask.mat, ficTracData_DAQ/_dat, and
% csv\trialSettings.csv. Both channels come from the same imagingData.mat
% so they share one registration. process_ft/process_im are copied
% unchanged from dopamine_ionto_walking_claude.m (movmean, window 1, i.e.
% no extra smoothing), and ft.stims = ftData_DAQ.stim{1} is stored, since
% the per-stim ejection analysis needs it.
%
% Runs top to bottom (batch-friendly); a trial that fails to load is
% skipped and listed at the end rather than aborting the build.

%% find trials
roots = {'Z:\pablo\lpsp_p2x2_redo','Z:\pablo\lpsp_p2x2_syt8s'};
all_files = [];
for r = 1:numel(roots)
    assert(isfolder(roots{r}),'cannot find %s -- check drive mapping',roots{r})
    all_files = [all_files; dir(fullfile(roots{r},'**','registration','imagingData.mat'))]; %#ok<AGROW>
end
all_files = natsortfiles(all_files);
n_trials  = numel(all_files);
fprintf('found %d trials\n',n_trials);

%% build
ft_type = 'movmean'; ft_win = 10;
im_type = {'movmean','movmean'}; im_win = {1,1};
n_centroid = 16; f0_pct = 7;

all_data = struct('ft',{},'im',{},'atp',{},'meta',{});
failed   = {};
tic
for i = 1:n_trials
    trial_dir = fileparts(all_files(i).folder);
    fprintf('%d/%d %s ',i,n_trials,trial_dir)
    try
        s = load(fullfile(all_files(i).folder,'imagingData.mat'),'img');
        m = load(fullfile(trial_dir,'mask.mat'),'mask');
        f = dir(fullfile(trial_dir,'*ficTracData_DAQ.mat')); q = load(fullfile(f(1).folder,f(1).name),'ftData_DAQ');
        f = dir(fullfile(trial_dir,'*ficTracData_dat.mat')); w = load(fullfile(f(1).folder,f(1).name),'ftData_dat');
        ts = readtable(fullfile(trial_dir,'csv','trialSettings.csv'));

        k = numel(all_data)+1;
        all_data(k).ft = process_ft(q.ftData_DAQ, w.ftData_dat, ft_win, ft_type);
        all_data(k).ft.pattern = ts.patternPath{1};
        all_data(k).ft.stims   = q.ftData_DAQ.stim{1};
        all_data(k).im   = process_im(s.img{1}, im_win, im_type, m.mask, n_centroid, f0_pct);
        all_data(k).atp  = process_im(s.img{2}, im_win, im_type, m.mask, n_centroid, f0_pct);
        all_data(k).meta = all_files(i).folder;
    catch ME
        failed{end+1} = sprintf('%s (%s)',trial_dir,ME.message); %#ok<SAGROW>
        fprintf('FAILED ')
    end
    fprintf('ETR: %.2f h\n', toc/i*(n_trials-i)/3600)
end

out_file = fullfile('.data','lpsp_p2x2_redo_syt8s_20261029.mat'); % run from the repo root: addpath('ugly_figures/ionto_celltypes'); lpsp_p2x2_unsaved_build
save(out_file,'all_data','-v7.3')
fprintf('saved %d trials to %s\n',numel(all_data),out_file);
fprintf('%d failed:\n',numel(failed)); fprintf('  %s\n',failed{:});

%% functions
function s = process_ft(ftData_DAQ, ftData_dat, ft_win, ft_type)
    % unchanged from dopamine_ionto_walking_claude.m
    f_speed = ftData_dat.velFor{:};
    f_speed = interp1(ftData_dat.trialTime{1},f_speed,seconds(ftData_DAQ.trialTime{1}),'linear','extrap');
    r_speed = ftData_DAQ.velYaw{:};
    cue     = ftData_DAQ.cuePos{:}' / 192 * 2 * pi - pi;
    cue(abs(gradient(cue)) > 2) = nan;
    cue     = unwrap(cue);
    cue     = smoothdata(cue,1,ft_type,ft_win,'omitnan');
    cue     = mod(cue,2*pi);
    cue(cue > pi) = cue(cue > pi) - 2*pi;

    s.xf      = seconds(ftData_DAQ.trialTime{:});
    if ismember('volClock',ftData_DAQ.Properties.VariableNames)
        s.xb  = seconds(ftData_DAQ.volClock{:});
    end
    s.f_speed = smoothdata(f_speed,1,ft_type,ft_win);
    s.r_speed = smoothdata(r_speed,1,ft_type,ft_win);
    s.cue     = cue;
end

function s = process_im(imgData, im_win, im_type, mask, n_centroid, f0_pct)
    % unchanged from dopamine_ionto_walking_claude.m (imgData is Y x X x T)
    imgData = smoothdata(imgData,3,im_type{1},im_win{1});
    imgData = imgData - min(imgData,[],'all');

    [y_mask,x_mask] = find(mask);
    min_axis        = min(range(x_mask),range(y_mask));
    mid             = bwskel(mask,'MinBranchLength',min_axis);
    [y_mid,x_mid]   = find(mid);

    [x_mid,y_mid]   = graph_sort(x_mid,y_mid);

    xq          = [-min_axis:(length(x_mid)+min_axis)];
    x_mid       = round(interp1(1:length(x_mid),x_mid,xq,'linear','extrap'));
    y_mid       = round(interp1(1:length(y_mid),y_mid,xq,'linear','extrap'));

    idx         = ismember([x_mid',y_mid'],[x_mask,y_mask],'rows');
    x_mid       = x_mid(idx);
    y_mid       = y_mid(idx);

    xq          = linspace(1,length(y_mid),2*(n_centroid*2) + 1)';
    centroids   = [interp1(1:length(y_mid),y_mid,xq),interp1(1:length(x_mid),x_mid,xq)];
    centroids   = centroids(2:2:end-1,:);
    [~,idx] = pdist2(centroids,[y_mask,x_mask],'euclidean','smallest',1);

    imgData_2d      = reshape(imgData,[],size(imgData,3));
    centroid_log    = false(2*n_centroid,size(imgData_2d,1));
    for i = 1:2*n_centroid
        centroid_log(i, sub2ind(size(imgData),y_mask(idx==i),x_mask(idx==i))) = true;
    end
    f_cluster       = centroid_log * double(imgData_2d) ./ sum(centroid_log,2);
    f0              = prctile(f_cluster,f0_pct,2);
    dff_cluster     = (f_cluster - f0) ./ f0;
    zscore_cluster  = zscore(dff_cluster,[],2);

    alpha       = linspace(-pi,pi-(2*pi/n_centroid),n_centroid);
    alpha       = repmat(alpha,1,2);

    [x_tmp,y_tmp]   = pol2cart(alpha,zscore_cluster');
    [mu,rho]        = cart2pol(mean(x_tmp,2),mean(y_tmp,2));
    mu = smoothdata(unwrap(mu),1,im_type{2},im_win{2});
    mu = mod(mu,2*pi);
    mu(mu > pi) = mu(mu > pi) - 2*pi;
    rho = smoothdata(rho,1,im_type{2},im_win{2});

    s.mu = mu;
    s.rho= rho;
    s.z  = zscore_cluster;
    s.d  = dff_cluster;
    s.f  = f_cluster;
    s.alpha = alpha;
end
