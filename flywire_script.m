%% clear all vars
clear all
close all

%% import the data
fulldata.names      = readtable('.flywire_data\names.csv');
fulldata.synapses   = readtable('.flywire_data\synapse_coordinates.csv');
fulldata.type       = readtable('.flywire_data\consolidated_cell_types.csv');
fulldata.connections = readtable('.flywire_data\connections.csv');

%% subselect the data for faster performance
names = {'LPsP'};

root_ids = fulldata.type.root_id(contains(fulldata.type.primary_type,names));
partner_ids = [fulldata.connections.pre_root_id(any(fulldata.connections.post_root_id == root_ids',2));...
               fulldata.connections.post_root_id(any(fulldata.connections.pre_root_id == root_ids',2))];
all_ids = [root_ids;partner_ids];

data.names = fulldata.names(any(fulldata.names.root_id == all_ids',2),:);
data.type = fulldata.type(any(fulldata.type.root_id == all_ids',2),:);
data.connections = fulldata.connections(any(fulldata.connections.pre_root_id == all_ids',2) | ...
                                    any(fulldata.connections.post_root_id == all_ids',2),:);

all_ind = find(~isnan(fulldata.synapses.pre_root_id)); %extract index to all new pre synapse cells
pre_ind = find(any(fulldata.synapses.pre_root_id==all_ids',2)); %extract index to presynapses cells we care about
post_ind = find(any(fulldata.synapses.post_root_id==root_ids',2)); %extract index to post-synapses that we care about

start_ind1 = pre_ind; %find where presynapses start for our labelled cells
tmp_ind = all_ind - pre_ind'; %create matrix showing distance to all start_inds from our designated cell starts
tmp_ind(tmp_ind<=0) = inf; %send anything that is backwards to infinity. keep only forward indices
end_ind1 = start_ind1 + min(tmp_ind)'-1;

s = start_ind1;
e = end_ind1;

% start_ind2 = all_ind - post_ind'; %create matrix showing distance to each new presynapse type from our cells of interest as post-synapse
% start_ind2(start_ind2>0) = -inf; %send all indexes after our desired neurons to -inf
% start_ind2 = post_ind + max(start_ind2)'; %the starting cell we care about will be the least negative index distance
% start_ind2 = unique(start_ind2); %okay, now we have all the cells upstream of our cells of interest
% tmp_ind = all_ind - start_ind2'; %repeat above process, but for new cells
% tmp_ind(tmp_ind<=0) = inf; %send anything that is backwards to infinity. keep only forward indices
% end_ind2 = start_ind2 + min(tmp_ind)'-1;
% 
% s = unique([start_ind1;start_ind2]);
% e = unique([end_ind1;end_ind2]);

l = repmat(1:length(s),2,1);
ind=eval(['[' sprintf('s(%d):e(%d) ',l(:)') ']']);

data.synapses = fulldata.synapses(ind,:);
%% go through the synapses table, and fill in nan values with the corresponding cell id
tic
for i = 1:size(data.synapses,1) 
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
    fprintf('%.2f ETR: %.2f hours\n',i/size(data.synapses,1),(toc/i)*(size(data.synapses,1)-i)/60/60)
end

%%
[~,tmp_idx] = unique(data.type.root_id);
T = join(data.connections,data.type(tmp_idx,:),'LeftKeys','pre_root_id','RightKeys','root_id');

%% plot LPsP post-synaptic partners, with nt_type

lpsp_ids = data.type.root_id(contains(data.type.primary_type,'LPsP'));

lpsp_syn = data.connections(data.connections.pre_root_id == lpsp_ids(1),:);