
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

%% Detrended fluctuation analyisis

intervals = logspace(2,4,100);  % twin= 10 to 10000 timepoints
intervals= floor(intervals); 

beta0= dfa(v0, intervals); 
beta1= dfa(v1, intervals); 
beta2= dfa(v2, intervals); 

%% Avalanche dynamics
powerexp_fit('reset',1e5); 
 [Fn,tau, Fl,alph,...
 shape_cllps, avg_n_given_l, signuz_inv ]= compute_avalanche_dynamics(v0); 

 [Fn1,tau1, Fl1,alph1,...
 shape_cllps1, avg_n_given_l1, signuz_inv1 ]= compute_avalanche_dynamics(v1); 

 [Fn2,tau2, Fl2,alph2,...
 shape_cllps2, avg_n_given_l2, signuz_inv2 ]= compute_avalanche_dynamics(v2); 


close all; 
tiledlayout(1,3); 
nexttile; 
loglog(Fn); 
hold on;
loglog(Fn1); 
loglog(Fn2); 
nexttile; 
loglog(Fl); 
hold on;
loglog(Fl1); 
loglog(Fl2)






