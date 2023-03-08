local_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\'; %define the local directory to clear up
remote_dir           = '\\files.med.harvard.edu\Neurobio\Wilson Lab\pablo\lpsp_kir\'; %define the remote directory to make sure things are backed up to

%% check whether files exist in the remote dir
date_folders        = cellstr(ls(local_dir)); %store all the subfolders, which include the dates
date_folders(1:2) = []; %clean up this variable

for i = date_folders'
    local_trials = cellstr(ls([local_dir,i{1}]));
    local_trials(1:2,:) = [];
    remote_trials= cellstr(ls([remote_dir,i{1}]));
    remote_trials(1:2,:) = [];
    
    for j = local_trials'
        local_contents = cellstr(ls([local_dir,i{1},'\',j{1}]));
        remote_contents = cellstr(ls([remote_dir,i{1},'\',j{1}]));
        
        for k = local_contents'
            if isfolder([local_dir,i{1},'\',j{1},'\',k{1}])
                a=1
            elseif ~exists([remote_dir,i,'\',j,'\',k])
                    copyfile [local_dir,i,'\',j] [remote_dir,i,'\',j]
            end
        end
    end
end

