function trc = compute_trc(v0, k)
% COMPUTE_TRC     computes the time-resolved correlations for up to k lags
%
%       trc = compute_trc(v0, k)
%
% INPUTS
%
%       v0          (n x t) matrix
%       k           maximum integer lag
%
% OUTPUTS
%       trc         (k x t) time-resolved correlations up to k lags

[n, t] = size(v0);
% zscore timeseries along time
v0 = v0 - mean(v0, 1);
v0 = v0 ./ std(v0, [], 1);

% initialize time-resolved correlation
trc = zeros(k,t);

for i = 1:k  % iterate over k lags
    trc(i, 1:i) = nan;
    trc(i, i+1:end) = sum(v0(:, 1:end-i) .* v0(:, i+1:end)) ./ (n - 1);
end
