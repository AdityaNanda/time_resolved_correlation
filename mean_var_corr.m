function v1 = mean_var_corr(n, trc, mv)
% MEAN_VAR_CORR     sampling timeseries with time-resolved correlation, variance and mean
%
%       v1 = mean_var_corr(n, trc, mv)
%
% INPUTS
%
%       n           number of nodes in model data
%
%       trc         (k x t) matrix of time-resolved correlations constraints
%                   this matrix can be obtained with the script compute_trc.m
%
%                   height of trc (k) determines the maximum lag for which
%                   correlations are to be preserved. some examples include:
%                       k = 1   % preserve correlations between adjacent timepoints.
%                               % (typical case for ephys results in Nanda et al., 2023).
%
%                       k = 3   % preserve correlations at the next three nearest points.
%
%
%       mv          (2 x t) matrix of mean and variance at each timepoint (optional)
%                   by default v1 is set to have zero mean unit variance at each timepoint
%
% OUTPUTS
%       v1          (n x t) synthetic timeseries with time-resolved
%                       mean, var, and corr constraints.

[k, t] = size(trc);

% check inputs
if ~exist('n', 'var') || isempty(n)
    error('Number of regions must be provided.');
end

% zero-mean and unit variance
if ~exist('mv', 'var') || isempty(mv)
    mv = [zeros(1,t); ones(1,t)];
end

v1 = zeros(n, t);                           % initialize output
seed = rand(n, 1);                          % sample only enough seed columns to start sampling
seed = zscore(seed, [], 1);                 % sample zscored seed vector
v1(:, 1) = seed;                            % the first timepoint is not constrained so use seed

for tt = 2:t                                % iterate over timepoints
    nb_ind = tt - (1:k);                    % index of constrained timepoints
    nb_ind(nb_ind<1) = [];                  % remove timepoints less than 1
    ts = v1(:, nb_ind);                     % activity at previous m timepoints
    c1 = trc(:, tt);                        % correlations to the previous m timepoints
    c1(isnan(c1)) = [];

    % compute b  and A (setting up equations as Ax=b)
    b = [(n-1).*c1; 0];                     % (n-1) is a scaling parameter
    A = [ts'; ones(1,n)];                   % coefficient matrix
    Z = null(A);
    xstar = lsqminnorm(A,b);
    d = sqrt((n-1) - norm(xstar)^2);        % use Pythagoras theorem to compute d
    Nvar = normrnd(0, 1, [size(Z,2), 1]);   % sample normal distribution
    q = Nvar ./ vecnorm(Nvar, 2, 1);        % uniformly sample q
    v1(:, tt) = xstar + d.*Z*q;             % sample timeseries
end

% scale back to user supplied mean and variance
v1 = v1 .* sqrt(mv(2,:)) + mv(1,:);
