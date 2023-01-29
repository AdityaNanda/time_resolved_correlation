function [xs, hsol] = single_sum_sampler_3D(A,b,samp,xsig)
% analytical sampling for when 
% A is a row vector, b is a scalar
% coeffs  
% size(A)= [1,t]

%% additionally 
% if b OR xsig is a vector (A must always be a row vector)
% then xs and hsol are 3D vector 
% the third dimension returns solution for each element of b or xsig

% OUTPUTS
% xs is the sampled solution
% hsol is the homogenous solution

if ~exist('samp', 'var') || isempty(samp)
    samp=1; 
end
assert(isvector(A));
n2= length(A);
ffs= @(x0) null_func(x0,n2); 

opt=optimoptions('fsolve','display','off');
x0= fsolve(ffs, rand(3,1),opt); 
ah=x0(1); 
bh=x0(2); 
ch=x0(3); 

ns= length(b); 
if length(xsig)==1
  xisg= repmat(xsig,ns);
else 
  assert(length(xsig)==ns);  
  % is xsig is not of length 1, assert b and xsig have same length
end

Samp=samp*ns;

rda = RandomDirAlgo(n2-1,Samp); 
qvec =zeros(n2,Samp); 

qvec(1,:)= sum(-ch.*rda,1); 
qvec(2:end,:)= repmat(sum(-ah.*rda), n2-1,1); 
qvec(2:end,:)= qvec(2:end, :)+rda.*(bh+ah);

% x0= lsqminnorm(A,double(b)); % do not use mldivide for this-->  lsqminnorm is best (2/19/2020)
%radi = sqrt(xnorm.^2-norm(x0)^2); 
%radi= sqrt(((n2-1).*xsig.^2)+n2*mean(x0).^2-norm(x0).^2); 
% for this case,  x0= b*ones(n2,1)
% thus mean(x0)reduces to b
% and norm x0 reduces to sqrt(n2)*b 
  radi= sqrt(((n2-1).*xsig.^2));
%  radi= sqrt(((n2-1).*xsig.^2));
 hsol= reshape(radi,1,1,ns).*reshape(qvec, [n2, samp,ns]); 
xs= reshape(b, 1,1,ns)+hsol;

end

