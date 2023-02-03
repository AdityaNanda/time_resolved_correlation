function v1 = mean_var(n, mv)
% MEAN_VAR          sampling timeseries with time-resolved correlation, variance and mean
%
%       v1 = mean_var_corr(n, mv)
%
% INPUTS
%
%       n           number of nodes in model data
%
%       mv          (2 x t) matrix of mean and variance at each timepoint
%
% OUTPUTS
%       v1          (n x t) synthetic timeseries with time-resolved
%                       mean, var, and corr constraints.

[~, t]= size(mv);

% define equations for analytical computation of nullspace
null_func = @(y) [ ...
    y(2)-y(1)*(n-2)-y(3), ...           % null constraint
    y(3)^2+y(1)^2*(n-2)+y(2)^2-1, ...   % unit norm constraint
    y(3)^2-2*y(1)*y(2)+(n-3)*y(1)^2]';  % orthonormality constraint

% solve equation using fsolve
opts = optimoptions('fsolve', ...
    'display', 'off', 'FunctionTolerance', 1e-30, ...
    'algorithm', 'Levenberg-Marquardt');
x0 = fsolve(null_func, rand(3,1),opts);
alpha = x0(1);
beta = x0(2);
gamma = x0(3);

% sample q
Nvar = normrnd(0, 1, [n-1, t]);         % sample normal dist
q= Nvar ./ vecnorm(Nvar, 2, 1);         % uniformly sample q

% compute the product Zq
% this can be done without writing out the matrix Z for efficiency
Zq = zeros(n,t);
Zq(1  , 1:t) = sum(-gamma.*q,1);
Zq(2:n, 1:t) = repmat(sum(-alpha.*q), n-1,1);
Zq(2:n, 1:t) = Zq(2:end, :) + q.*(beta+alpha);

% compute solution and scale to user-defined mean and variance
d = sqrt(n - 1);
v1 = d.*Zq;
v1 = v1 .* sqrt(mv(2, :)) + mv(1, :);
