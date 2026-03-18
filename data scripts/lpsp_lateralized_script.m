%% load in data
base_dir = 'Z:\pablo\lpsp_lateralized\'; %uigetdir(); %
all_files = dir([base_dir,'\**\imgData_reg.mat']);
all_files = natsortfiles(all_files);

right_rho = cell(length(all_files),1);
left_rho  = cell(length(all_files),1);
dark_idx  = false(length(all_files),1);

for i = 1:length(all_files)
    fprintf('%i / %i\n',i,length(all_files))
    try    
    load([all_files(i).folder,'\',all_files(i).name])  

    tmp2 = dir([fileparts(all_files(i).folder),'\*ficTracData_DAQ.mat']);
    load([tmp2.folder,'\',tmp2.name])

    tmp2 = dir([fileparts(all_files(i).folder),'\csv\trialSettings.csv']);
    tmp2 = readtable([tmp2.folder,'\',tmp2.name]);
    dark_idx(i) = contains(tmp2.patternPath{1},'background');

    velYaw = smoothdata(ftData_DAQ.velYaw{:},1,'gaussian',10);

    xf = linspace(0,10,size(velYaw,1));
    xb = linspace(0,10,size(imgData_reg,3));
    
    all_pix = reshape(imgData_reg,[],size(imgData_reg,3));
    all_pix = interp1(xb,all_pix',xf)';
    
    right_rho{i} = reshape(corr(all_pix',max(velYaw,0)),size(imgData_reg,1),size(imgData_reg,2));
    left_rho{i} = reshape(corr(all_pix',max(-velYaw,0)),size(imgData_reg,1),size(imgData_reg,2));
    end
end

%%

[~,~,fly_num] = unique(cellfun(@(x)(x(1:40)),{all_files.folder}','UniformOutput',false));

im_gain = 10;

t = tiledlayout("flow");
for i = 1:length(right_rho)
    tot_rho = (right_rho{i} .* reshape([1,.5,0],1,1,3) * im_gain) + ...
                (left_rho{i} .* reshape([0,.5,1],1,1,3) * im_gain);
    
    nexttile
    image(rot90(tot_rho,2))
    xticks([])
    yticks([])
    title(sprintf('fly num: %i\ndark: %i',fly_num(i),dark_idx(i)),'Color','w')
end
set(gcf,'Color','none','InvertHardcopy','off')