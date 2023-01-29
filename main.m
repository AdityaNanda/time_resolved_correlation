
clc;
close all;
clear;
addpath('dfa/'); 
addpath('avalanches/'); 

mat= load('data/ephys_human_seeg.mat');

v0= mat.v0;
time= mat.time;
[n,t]= size(v0);   % regions and timepoints

%% compute trime-resolved correlation
lags=[1:3]; % trc at 3 lags
trc= compute_trc(v0,lags);

% build model timeseries
v1= mean_var_corr(lags,trc, n);
% compute trc in model timeseries
trc1= compute_trc(v1,lags);

figure; 
tiledlayout(1,max(lags))
for jj=1:max(lags)
   nexttile; 
 plot(trc(jj,:), trc1(jj,:), '.');
end

%% mean-var
mv= [mean(v0,1); var(v0,[],1)]; 
v2=  mean_var(mv,n); 

figure; 
tiledlayout(1,2); 
nexttile; 
plot(mean(v2, 1),mean(v0,1)); 
nexttile; 
plot(std(v2,[],1), std(v0,[],1)); 




