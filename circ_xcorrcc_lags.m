function [r,lags] = circ_xcorrcc_lags(x, y, max_lag)
% CIRC_XCORRCC_LAGS  Cross-correlation of two circular signals up to max_lag.
% Uses circ_corrcc on overlapping samples like xcorr does for linear data.
%
% Inputs:
%   x, y     - angle vectors (radians), same length
%   max_lag  - nonnegative integer, <= numel(x)-1
%
% Outputs:
%   r        - correlation at each lag (column), length 2*max_lag+1
%   lags     - integer lags (column), from -max_lag:max_lag

    x = x(:); y = y(:);
    assert(numel(x)==numel(y), 'x and y must have the same length.');
    N = numel(x);
    if nargin < 3 || isempty(max_lag), max_lag = N-1; end
    max_lag = min(max_lag, N-1);

    lags = (-max_lag:max_lag).';
    r    = nan(numel(lags),1);

    for i = 1:numel(lags)
        k = lags(i);
        if k >= 0
            xi = x(1:N-k);  yi = y(1+k:N);
        else
            kk = -k;        xi = x(1+kk:N); yi = y(1:N-kk);
        end
        % optional NaN handling
        m = ~isnan(xi) & ~isnan(yi);
        if any(m)
            r(i) = circ_corrcc(xi(m), yi(m));  % assumes circ_corrcc returns rho
        end
    end
    lags = -lags; % this sign inversion was neccessary for the output to be not inversely correlated with the output of xcorr.
end
