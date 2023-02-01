function [trc]= compute_trc(v0, lags)
% This script computes the time-resolved correlations between
% lags present in the vector "lags".

% trc is efficiently computed by the product of shifted, z-scored
% timeseries to itself.

% For convenience, the script always returns
% trc values with length 't' irrespective of lags.
% The additional values are set to nan.  For instance, a lag =1,
% generates t-1 unique time-resolved correlations, and 1 nan value.

%% INPUTS

% v0      -> is a timeseries (n x t)

% lags    -> is a vector of integer lags

%% OUTPUTS

% trc     -> is a k x t vector of time-resolved correlations

[n,t]=size(v0);
% zscore timeseries along time
v0= v0- mean(v0,1);
v0= v0./std(v0,[],1);

% initialize time-resolved correlation
trc=zeros(length(lags),t);

for kk=1:length(lags)  % iterate over k lags

    trc(kk,1:lags(kk))= nan;
    trc(kk, lags(kk)+1:end)= 1./(n-1).*sum(v0(:,1:end-lags(kk)).*v0(:,lags(kk)+1:end));

end

