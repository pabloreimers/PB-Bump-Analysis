function [H] = circ_entropy(X, dt)
%% This function takes an array and returns the Shannon Entropy of an array.
% A method for finding out how uniform the distribution of a circular
% variable is

%% Create an empirical probability distribution p for the values of X:

bins = -pi:dt:pi;
%bins(end) = bins(end) + 1; % eliminates the edge case by setting max(X) < bins(end)
X = sort(X);
p = zeros(1,length(bins)-1);
for i = 1:(length(bins)-1)
    p(i) = sum((X>=bins(i)) & (X<bins(i+1)));
end
p = p/sum(p);   % normalize to PDF
p(p==0) = [];   % eliminate empty bins
H = sum(-p.*log2(p)) / log2(length(bins));