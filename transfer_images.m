local_dir = uigetdir();
remote_dir = uigetdir();
folders   = dir(localdir);
folders(1:2) = [];

for i = 1:length(folders)
    tmp = dir([local_dir,'\',folders(i).name,'\*mask*']);
    for j = 1:length(tmp)
        copyfile([local_dir,'\',folders(i).name,'\',tmp(j).name],[remote_dir,'\',folders(i).name,'\',tmp(j).name])
    end
end