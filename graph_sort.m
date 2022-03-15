function [x,y] = graph_sort(x,y,x0,y0)
if nargin < 3
    x0 = x(1);
    y0 = y(1);
end
G = graph(pdist2([x,y],[x,y],'Chebychev') == 1);    %create a graph that connects adjacent points in an image. chebychev distance is chess move distance, aka touching pixels
H = rmedge(G,find(all(ismember(G.Edges.EndNodes,find(degree(G) > 2)),2))); %if a pixel is connected to more than 2 pixels, remove the unnecessary edge

k = find(degree(H) == 1, 1);        %define the starting point as the first node with only one adjacent pixel. THIS IS FRAGILE FOR CIRCLES
v = dfsearch(H,k);                  %discover all points in the graph, starting at our starting point. return the index of discovery (such that you traverse the shortest path)
x = x(v);                           %reorder the coordinates in order of discovery
y = y(v);

%% OLD CODE
% I = knnsearch([x,y],[x,y],'k',8,'Distance','Chebychev');   %find the 3 nearest neighbors to each point (including self) in chebychev distance (chess move distance). basically picks out the two pixels each pixel is touching
% G = graph;                          %initialize a graph object
% s = 1:size(x);                      %set nodes as the index to each point
% G = addedge(G,s,I(:,2));            %add an edge connecting each point to its nearest neighbor
% G = addedge(G,s,I(:,3));            %add an edge connecting each points to its second nearest neighbor
% k = find(x==x0 & y==y0);            %find the index of our starting point