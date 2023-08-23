local_dir            = 'C:\Users\ReimersPabloAlejandr\Documents\Data\2p data\current\'; %define the local directory to clear up
remote_dir           = '\\files.med.harvard.edu\Neurobio\wilsonlab\pablo\gain_change\'; %define the remote directory to make sure things are backed up to

%% check whether files exist in the remote dir
date_folders        = cellstr(ls(local_dir)); %store all the subfolders, which include the dates
date_folders(1:2) = []; %clean up this variable

for i = date_folders'
    tmp_local = [local_dir,'\',i{1}];
    tmp_remote= [remote_dir,'\',i{1}];
    
    local_flies = cellstr(ls(tmp_local));
    local_flies(1:2,:) = [];
    remote_flies= cellstr(ls(tmp_remote));
    remote_flies(1:2,:) = [];

    for f = local_flies'
    
    tmp_local = [local_dir,'\',i{1},'\',f{1}];
    tmp_remote= [remote_dir,'\',i{1},'\',f{1}];

    local_trials = cellstr(ls(tmp_local));
    local_trials(1:2,:) = [];
    remote_trials= cellstr(ls(tmp_remote));
    remote_trials(1:2,:) = [];

    for j = local_trials'
        tmp_local = [local_dir,'\',i{1},'\',f{1},'\',j{1}];
        tmp_remote= [remote_dir,'\',i{1},'\',f{1},'\',j{1}];

        local_contents = cellstr(ls(tmp_local));
        local_contents(1:2) = [];
        remote_contents = cellstr(ls(tmp_remote));
        remote_contents(1:2) = [];
        
        fprintf(['current folder: ',j{1},'\n'])

        outer_local = tmp_local;
        outer_remote = tmp_remote;
        for k = local_contents'
            tmp_local = [outer_local,'\',k{1}];
            tmp_remote= [outer_remote,'\',k{1}];            


            if isfolder(tmp_local) %if the current content is a folder, check whether it exists in the remote. if yes, delete it, if no, copy and then delete 
                if isfolder(tmp_remote)
                    if ~contains(k{1},'registration')
                        rmdir(tmp_local,'s')
                    end
                else
                    mkdir(tmp_remote)
                    copyfile(tmp_local,tmp_remote)
                end
                
            elseif ~isfile(tmp_remote) %if it's a not a folder, check whether it exists in remote dir. if no, copy it and delete it locally
                copyfile(tmp_local, [remote_dir,i{1},'\',f{1},'\',j{1}])
                if ~contains(k{1},'imagingData') && ~contains(k{1},'ficTracData_DAQ') && ~contains(k{1},'mask')
                    delete(tmp_local)
                end
            else %if it does exist in the remote dir, delete it locally unless it's the imaging file
                if ~contains(k{1},'imagingData') && ~contains(k{1},'ficTracData_DAQ') && ~contains(k{1},'mask')
                    delete(tmp_local)
                end
            end
        end
    end
    end
    fprintf(['current folder: ',j{1},'\n'])
end

