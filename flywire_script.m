%% clear all vars
clear all
close all

%% import the data
data.names      = readtable('flywire_data\names.csv');
data.synapses   = readtable('flywire_data\synapse_coordinates.csv');
data.type       = readtable('flywire_data\consolidated_cell_types.csv');
data.connections = readtable('flywire_data\connections.csv');

%% define which neurons to subselect
names = {'LPsP','EPG','PEN1','PEG'};

root_ids = data.type.root_id(contains(data.type.primary_type,names));

data.names = data.names(any(data.names.root_id == root_ids',2),:);
data.type = data.type(any(data.type.root_id == root_ids',2),:);
data.connections = data.connections(any(data.connections.pre_root_id == root_ids',2) | ...
                                    any(data.connections.post_root_id == root_ids',2),:);

%find the indexes to all non-nan cells in synapses
idx1 = ~isnan(data.synapses.pre_root_id)

%%
last_pre = 'tmp';
last_post= 'tmp';

for i = 1:size(data.synapses,1) %go through the synapses table, and fill in nan values with the corresponding cell id
    if isnan(data.synapses.pre_root_id(i))
        data.synapses.pre_root_id(i) = last_pre;
    else
        last_pre = data.synapses.pre_root_id(i);
    end

    if isnan(data.synapses.post_root_id(i))
        data.synapses.post_root_id(i) = last_post;
    else
        last_post = data.synapses.post_root_id(i);
    end
end

%%
[~,tmp_idx] = unique(data.type.root_id);
T = join(data.connections,data.type(tmp_idx,:),'LeftKeys','pre_root_id','RightKeys','root_id');

%% plot LPsP post-synaptic partners, with nt_type

lpsp_ids = data.type.root_id(contains(data.type.primary_type,'LPsP'));

lpsp_syn = data.connections(data.connections.pre_root_id == lpsp_ids(1),:);



