%% fig_TEMPLATE
% Copy this file to fig_<description>.m for a new figure. Run one %% section
% at a time (Ctrl+Enter).

%% load data
clear all
load(uigetfile('Select Processed Imaging Data'))

%% plot
figure('color','w'); clf
% ...

%% export
exportgraphics(gcf, fullfile('..','exports','fig_TEMPLATE.png'), 'Resolution', 300)
