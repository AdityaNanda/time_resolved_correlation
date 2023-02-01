function [v1]= mean_var_corr(lags,trc,n, mv)

% MEAN_VAR_CORR - > sampling timeseries with
%                                                   time-resolved correlation
%                                                   (and mean and variance)

%% INPUTS

% lags            - > is a k x 1 vector of integers lags for which
%                       correlations are to be preserved.

%                       If lags = 1 (typical case for all ephys results in paper),
%                       then correlations between adjacent timepoints are preserved.

%                       If lags= [1,2,3], correlations to the 3 nearest timepoints are
%                       preserved and so on

%                       lags should be contiguous integers and start with 1.
%                       For instance, inputs of the form
%                       lags =[1,6, 9] or lags= [2,3,4] are not allowed.

% trc             -> is the set of (k x t) matrix of time-resolved correlations to be satisfied
%                      where k is length(lags)

%                      if lags= 1, then trc must be 1 x t
%                      if lags= [1,2,3], trc must be 3 x t

%                      Use the "compute_trc.m" script to obtain trc

% 'n' is the number of nodes.

%% OPTIONAL INPUTS
% 'mv'           -> (mean-variance) must contain a 2 x t array of mean and variance at each
%                      timepoint
%                      If this is not supplied, v1 is standardized at each timepoint
%                      by setting mean=0 and var=1

%% OUTPUTS
% v1               -> is an n x t surrogate timeseries.
%                       If 'mv' is supplied, v1 has mean mv(1,:) and variance mv(2,:)
%                       v1 satisfies the correlations as prescribed in
%                       'lags' and 'trc'
%                       If mv is not supplied,
%                       v1 is z-scored along time (each timepoint has mean 0 and std 1)

[k,t]= size(trc);

% check inputs
if ~exist('n', 'var')|| isempty(n)
    error('Number of regions must be provided with time-resolved correlation');
end
if  k~=length(lags)
    error('Enough time-resolved corr values are not provided');
end

if ~exist('mv', 'var') || isempty(mv)
    mv= [zeros(1,t); ones(1,t)];  % zero-mean and unit variance
end

nbmin= min(lags);  % closest correlation constraint
seed= rand(n, nbmin);  % sample only enough seed columns to start sampling
seed= zscore(seed,[],1);   % sample zscored seed vector

v1= zeros(n,t); % initialize output
v1(:, 1:nbmin)= seed;  % the first timepoint is not constrained so use seed

for tt=nbmin+1:t    % iterate over timepoints
    nb_ind= tt-lags; % index of constrained timepoints
    nb_ind(nb_ind<1)=[];   % remove timepoints less than 1
    ts= v1(:, nb_ind);  %  activity at previous m timepoints
    c1= trc(:, tt);  %  correlations to the previous m timepoints
    c1(isnan(c1))=[];

    % compute b  and A ( setting up equations as Ax=b)
    b= [(n-1).*c1;0];   % (n-1) is a scaling parameter (see derivation in Methods)
    A= [ts'; ones(1,n)];  % coefficients
    Z=null(A);
    xstar= lsqminnorm(A,b);
    d= sqrt((n-1) - norm(xstar)^2);  % use pythagoras theorem to compute d
    Nvar= normrnd(0,1,[size(Z,2),1]);  % sample normal dist
    q= Nvar./(vecnorm(Nvar,2,1));   % uniformly sample q
    v1(:, tt)=xstar+d.*Z*q;
end
v1= v1.*sqrt(mv(2,:)) + mv(1,:); % scale back to user supplied mean and variance
end
