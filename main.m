
clc;
close all;
clear;

% human ECoG data
% v0, timeseries of a single subject
% time, time vector

mat= load('data/human_ieeg.mat');
v0= mat.v0;
time= mat.time;
[n,t]= size(v0);   % regions and timepoints

%% time-resolved correlation model

% compute trime-resolved correlation
lags=[1:3]; % trc at 3 lags
trc= compute_trc(v0,lags);

% build model timeseries
v1= mean_var_corr(lags,trc, n);
% compute trc in model timeseries
trc1= compute_trc(v1,lags);

% visual depiction of trc and trc1
figure; 
tiledlayout(1,max(lags))
for jj=1:max(lags)
   nexttile; 
 plot(trc(jj,:), trc1(jj,:), '.');
end

%% time-resolved mean and variance model 

mv= [mean(v0,1); var(v0,[],1)]; 
v2=  mean_var(mv,n); 

figure; 
tiledlayout(1,2); 
nexttile; 
plot(mean(v2, 1),mean(v0,1),'.'); 
nexttile; 
plot(std(v2,[],1), std(v0,[],1),'.'); 




