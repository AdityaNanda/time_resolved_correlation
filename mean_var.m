function [v1] = mean_var(mv,n)

% MEAN_VAR_SAMPLER   ->  sampling timeseries with
%                                            constrained mean and variance
%                                            along timepoints

% INPUTS
% 'mv'              ->  (mean-variance) must contain a 2 x t array
%                           of mean and variance at each timepoint

% 'n'                 ->  is the number of nodes/regions


% OUTPUTS
% 'v1'               -> is an n x t surrogate timeseries.
%                         v1 has mean mv(1,:) and variance mv(2,:)

[~, t]= size(mv);


% define equations for parameters of nullspace (Methods)
null_func= @(y) [ y(2)-y(1)*(n-2)-y(3),...   % null constraint
    y(3)^2+y(1)^2*(n-2)+y(2)^2-1, ...    % unit norm constraint
    y(3)^2-2*y(1)*y(2)+(n-3)*y(1)^2]';  % orthonormality constraint


% solve equation using fsolve
opt=optimoptions('fsolve','display','off','FunctionTolerance', 1e-30, ...
    'algorithm', 'Levenberg-Marquardt');
x0= fsolve(null_func, rand(3,1),opt);
alpha=x0(1);
beta=x0(2);
gamma=x0(3);


% sample q
Nvar= normrnd(0,1,[n-1,t]);  % sample normal dist
q= Nvar./(vecnorm(Nvar,2,1));   % uniformly sample q


% compute the product Zq
% this can be done without writing out the matrix Z for efficiency
Zq =zeros(n,t);
Zq(1,1:t)= sum(-gamma.*q,1);
Zq(2:n,1:t)= repmat(sum(-alpha.*q), n-1,1);
Zq(2:n,1:t)= Zq(2:end, :)+q.*(beta+alpha);


% compute solution and scale to user-defined mean and variance
% xtilde =0
d=sqrt(n-1);
v1= d.*Zq;
v1= v1.*sqrt(mv(2,:))+mv(1,:);

end

