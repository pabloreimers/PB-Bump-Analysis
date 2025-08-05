function [H] = entropy(X, varargin)
%% This function takes an array and returns the Shannon Entropy of an array.
% A probability density function is derived from the array.
%
% By default the program selects a number of bins proportional to the square
% root of the number of points in the array, however the Freedman?Diaconis rule ('fd') 
% and Sturge's rule ('sturges') can be selected.
%
%
% Ex: X = randn(1000,1);
%
% H_default = Entropy_Array(X);
% H_fd = Entropy_Array(X,'fd');
%
%
% 3/30/16 - Ben Lucas
%% Create an empirical probability distribution p for the values of X:
%% Select number of bins - This uses a square root method by default,
if isempty(varargin)
    num_bins = ceil(sqrt(length(X)));
else method = varargin{1}
    if strcmpi(method,'fd')
        num_bins = ceil(range(X)./(2*iqr(X)*length(X)^(-1/3)));
    elseif strcmpi(method, 'sturges')
        num_bins = ceil(1 + log2(length(X)));
    elseif strcmpi(method, 'sqrt')
        num_bins = ceil(sqrt(length(X)));
    end
end
bins = linspace(min(X), max(X), num_bins);
bins(end) = bins(end) + 1; % eliminates the edge case by setting max(X) < bins(end)
X = sort(X);
p = zeros(1,length(bins)-1);
for i = 1:(length(bins)-1)
    p(i) = sum((X>=bins(i)) & (X<bins(i+1)));
end
p = p/sum(p);   % normalize to PDF
p(p==0) = [];   % eliminate empty bins
H = sum(-p.*log2(p));