function [trc]= compute_trc(v0, lags)
% % This script directly computes the time-resolved correlations between 
% lags present in kvec
% Alternatively,  we can say it computes kth -diagonal elements of the
% activity correlation matrix (t x t)
% Fastest way is to shift the zscored timeseries by k and multiply (along time) by itself.

%% INPUTS 
%  v0 is a timeseries (n x t)
% kvec is a vector of integer lags 
%% OUTPUTS 
% trc is a k x t vector of time-resolved correlations
% with integer lags in kvec

[n,t]=size(v0);
k=length(lags); 

% zscore timeseries along time 
v0= v0- mean(v0,1); 
v0= v0./std(v0,[],1);

% initialize time-resolved correlation
trc=zeros(length(lags),t); 
for kk=1:length(lags)  % iterate over k lags
     % compute correlations using product of shifted zscorte timeseries to itself
       trc(kk,1:lags(kk))= nan; 
    trc(kk, lags(kk)+1:end)= 1./(n-1).*sum(v0(:,1:end-lags(kk)).*v0(:,lags(kk)+1:end)); 
       %  set last kk values to nan (since these do not exist)
  
end

