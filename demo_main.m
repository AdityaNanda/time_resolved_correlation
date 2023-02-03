clc; clear; close all;

% load human ECoG data
% v0, timeseries of a single subject
% time, time vector

data = load('data/human_ieeg.mat');
v0 = data.v0;
time = data.time;
[n,t] = size(v0);   % regions and timepoints

%% compute constraints

% compute time-resolved correlation at 3 adjacent timepoints
k = 3;
trc0 = compute_trc(v0, k);

% compute mean-variance
mv0 = [mean(v0); var(v0)];

%% time-resolved correlation model

% build model timeseries
v1 = mean_var_corr(n, trc0, mv0);

% compute trc in model timeseries
trc1 = compute_trc(v1, k);

% plot empirical against model lags
figure;

tiledlayout(1, 3)

nexttile;
plot(mean(v0), mean(v1), '.');
xlabel('empirical time-resolved mean')
ylabel('model time-resolved mean')
axis square; axis tight;

nexttile;
plot(var(v0), var(v1), '.');
xlabel('empirical time-resolved variance')
ylabel('model time-resolved variance')
axis square; axis tight;

nexttile;
hold on;
for jj = max(k):-1:1
    plot(trc0(jj,:), trc1(jj,:), '.');
end
legend("lag " + (max(k):-1:1), 'location', 'northwest')
xlabel('empirical time-resolved correlations')
ylabel('model time-resolved correlations')
axis square; axis tight;

%% time-resolved mean and variance model

% build model timeseries
v2 = mean_var(n, mv0);

% plot empirical against model lags
figure;

nexttile;
plot(mean(v0), mean(v2), '.');
xlabel('empirical time-resolved mean')
ylabel('model time-resolved mean')
axis square; axis tight;

nexttile;
plot(var(v0), var(v2), '.');
xlabel('empirical time-resolved variance')
ylabel('model time-resolved variance')
axis square; axis tight;
