%% clear all vars
clear all
close all

%% import the data
fulldata.names      = readtable('.flywire_data\names.csv');
fulldata.synapses   = readtable('.flywire_data\synapse_coordinates.csv');
fulldata.type       = readtable('.flywire_data\consolidated_cell_types.csv');
fulldata.connections = readtable('.flywire_data\connections.csv');

%% subselect the data for faster performance
names = {'LPsP','EPG','PEN1','PEG'};

root_ids = fulldata.type.root_id(contains(fulldata.type.primary_type,names));
partner_ids = [fulldata.connections.pre_root_id(ismember(fulldata.connections.post_root_id,root_ids));...
               fulldata.connections.post_root_id(ismember(fulldata.connections.pre_root_id,root_ids))];
all_ids = [root_ids;partner_ids];

data.names = fulldata.names(ismember(fulldata.names.root_id,all_ids),:);
data.type = fulldata.type(ismember(fulldata.type.root_id,all_ids),:);
data.connections = fulldata.connections((ismember(fulldata.connections.pre_root_id,all_ids) | ...
                                         ismember(fulldata.connections.post_root_id,all_ids)),:); %weird indexing to use column vectors for speed

all_ind = find(~isnan(fulldata.synapses.pre_root_id)); %extract index to all new pre synapse cells
pre_ind = find(ismember(fulldata.synapses.pre_root_id,all_ids)); %extract index to presynapses cells we care about
post_ind = find(ismember(fulldata.synapses.post_root_id',root_ids)); %extract index to post-synapses that we care about

start_ind1 = pre_ind; %find where presynapses start for our labelled cells
tmp_ind = all_ind - pre_ind'; %create matrix showing distance to all start_inds from our designated cell starts
tmp_ind(tmp_ind<=0) = inf; %send anything that is backwards to infinity. keep only forward indices
end_ind1 = start_ind1 + min(tmp_ind)'-1;

s = start_ind1;
e = end_ind1;

l = repmat(1:length(s),2,1);
ind=eval(['[' sprintf('s(%d):e(%d) ',l(:)') ']']);

data.synapses = fulldata.synapses(ind,:);
%% go through the synapses table, and fill in nan values with the corresponding cell id
pre_ind = [find(~isnan(data.synapses.pre_root_id));size(data.synapses,1)+1];

if any(isnan(data.synapses.pre_root_id))
for i = 1:(length(pre_ind)-1)
    data.synapses.pre_root_id(pre_ind(i):pre_ind(i+1)-1) = data.synapses.pre_root_id(pre_ind(i));
    fprintf('%.2f\n',i/length(pre_ind))
end
end

post_ind = [find(~isnan(data.synapses.post_root_id));size(data.synapses,1)+1];

if any(isnan(data.synapses.post_root_id))
for i = 1:(length(post_ind)-1)
    data.synapses.post_root_id(post_ind(i):post_ind(i+1)-1) = data.synapses.post_root_id(post_ind(i));
    fprintf('%.2f\n',i/length(post_ind))
end
end

%% plot ps196 input to both LPsP

ps196 = data.type.root_id(find(contains(data.type.primary_type,'PS196')));
lpsp = data.type.root_id(find(contains(data.type.primary_type,'GLNO')));

conns = nan(length(ps196),length(lpsp));

figure(1); clf
subplot(2,1,1); hold on
for i = 1:length(ps196)
    for j = 1:length(lpsp)
        conns(i,j) = sum(data.connections.syn_count(...
                                                    data.connections.pre_root_id==ps196(i) & ...
                                                    data.connections.post_root_id==lpsp(j)));
        tmp_idx = data.synapses.pre_root_id==ps196(i) & data.synapses.post_root_id == lpsp(j);
        scatter(data.synapses.x(tmp_idx),data.synapses.y(tmp_idx),'.')
    end
end

subplot(2,1,2); heatmap(conns'); xlabel('PS196'); ylabel('LPsP')

%% show LPsP output in relation to targeted synapses
epg_id = data.type.root_id(contains(data.type.primary_type,'EPG'));
lpsp_id= data.type.root_id(contains(data.type.primary_type,'LPsP'));
pen1_id= data.type.root_id(contains(data.type.primary_type,'PEN1'));
peg_id= data.type.root_id(contains(data.type.primary_type,'PEG'));
d7_id= data.type.root_id(contains(data.type.primary_type,'Delta7'));

figure(4); clf
idx1 = ismember(data.synapses.pre_root_id,lpsp_id) & ismember(data.synapses.post_root_id,epg_id);


h(1) = subplot(2,3,1); hold on; set(gca,'YDir','reverse')
scatter(data.synapses.x(idx1),data.synapses.y(idx1),'filled','k')
idx2 = ismember(data.synapses.pre_root_id,epg_id) & ismember(data.synapses.post_root_id,d7_id);
scatter(data.synapses.x(idx2),data.synapses.y(idx2),'filled')
title('EPG > Delta7')
subplot(2,1,2); hold on
dist = min(pdist2([data.synapses.x(idx1),data.synapses.y(idx1),data.synapses.z(idx1)],...
              [data.synapses.x(idx2),data.synapses.y(idx2),data.synapses.z(idx2)]),[],1); %find distance from each synapse in idx1 to the nearest synapse in idx2
histogram(dist(dist<5000),'binwidth',100,'Normalization','probability')



h(2) = subplot(2,3,2); hold on; set(gca,'YDir','reverse')
scatter(data.synapses.x(idx1),data.synapses.y(idx1),'filled','k')
idx2 = ismember(data.synapses.pre_root_id,epg_id) & ismember(data.synapses.post_root_id,pen1_id);
scatter(data.synapses.x(idx2),data.synapses.y(idx2),'filled')
title('EPG > PEN1')
subplot(2,1,2); hold on
dist = min(pdist2([data.synapses.x(idx1),data.synapses.y(idx1),data.synapses.z(idx1)],...
              [data.synapses.x(idx2),data.synapses.y(idx2),data.synapses.z(idx2)]),[],1); %find distance from each synapse in idx1 to the nearest synapse in idx2
histogram(dist(dist<5000),'binwidth',100,'Normalization','probability')

h(3) = subplot(2,3,3); hold on; set(gca,'YDir','reverse')
scatter(data.synapses.x(idx1),data.synapses.y(idx1),'filled','k')
idx2 = ismember(data.synapses.pre_root_id,epg_id) & ismember(data.synapses.post_root_id,peg_id);
scatter(data.synapses.x(idx2),data.synapses.y(idx2),'filled')
title('EPG > PEG')
subplot(2,1,2); hold on
dist = min(pdist2([data.synapses.x(idx1),data.synapses.y(idx1),data.synapses.z(idx1)],...
              [data.synapses.x(idx2),data.synapses.y(idx2),data.synapses.z(idx2)]),[],1); %find distance from each synapse in idx1 to the nearest synapse in idx2
histogram(dist(dist<5000),'binwidth',100,'Normalization','probability')
legend('EPG > Delta7','EPG > PEN1','EPG > PEG')
title('LPsP > EPG distance to')

linkaxes(h)

%% plot nearest neighbor distance histograms

idx1 = ismember(data.synapses.pre_root_id,lpsp_id) & ismember(data.synapses.post_root_id,epg_id);
idx2 = ismember(data.synapses.pre_root_id,epg_id) & ismember(data.synapses.post_root_id,d7_id);
dist = min(pdist2([data.synapses.x(idx1),data.synapses.y(idx1),data.synapses.z(idx1)],...
              [data.synapses.x(idx2),data.synapses.y(idx2),data.synapses.z(idx2)]),[],2); %find distance from each synapse in idx1 to the nearest synapse in idx2

