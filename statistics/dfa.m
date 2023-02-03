function [beta1, flucts] = dfa(vt,intervals)
% DFA       function for detrended fluctuation analysis
%
%       [beta1, flucts] = dfa_wrapper(ts, intervals)
%
%   This function computes windowed fluctutations ('flucts') in each window
%   and fits a power-law with an exponenet 'beta1'. The function is
%   optimized for computing the fluctuations for 'n' timeseries in parallel.
%
% Inputs:
%
%   v0              (n x t) amplitude-envelope timeseries (see Methods).
%
%   intervals       (1 x d) vector of d integer time windows
%                           over which fluctuations are computed
%
% Outputs:
%
%   beta1           (n x 1) vector of power-law slope for each 1 x t timeseries
%
%   flucts          (n x length(intervals)) of windowed fluctuations in each window

[n, t] = size(vt);
intervals = reshape(intervals, [], 1);
flucts = zeros(n, length(intervals));

parfor dd = 1:length(intervals)
    twin = intervals(dd);
    nd = floor(t / twin);
    N1 = nd * twin;

    csum = cumsum(vt(:,1:N1), 2);
    yc = (csum - csum(:,end));

    % compute trend to be subtracted
    xm = ([(1:twin)', ones(twin,1)]);
    xpmat = compute_proj_matrix(xm,twin);
    yct = reshape(yc, n, twin, []);

    % multiply xct for each node using pagemtimes
    Yn = reshape(pagemtimes(yct, xpmat), n, N1);
    flucts(:, dd) = sqrt(sum((yc - Yn).^2, 2) / N1);
end

% compute slope in power-law
xx= [log(intervals), ones(size(intervals))];
yy= log(flucts);
temp = pinv(xx)*yy';
beta1 = temp(1, :);

end

function xxp = compute_proj_matrix(x, twin)
% This script computes x * pinv(x)

% Given x =[1:twin; ones(twin,1)]';
% pinv(x)*x is the projection matrix for detrending

xxt(1,1) = sum((1:twin).^2);
xxt(2,2) = twin;
xxt(1,2) = sum(1:twin);
xxt(2,1) = xxt(1,2);
xxp = x*(pinv(xxt)*x');

end