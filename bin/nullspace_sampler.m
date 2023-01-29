function [x1, varargout]=...
  nullspace_sampler(A, b, sig1,Z, samp)

% % A is a m x t matrix of m timeseries
%  b is m x 1 vector of products
% The output x1 is a 1 x t vector 
% Ax1=b; 

% sig1 constrains the std. dev of the output x1
% sig1 can also constrain norm by setting it to be negative
% thus, sig1=1 =>  std(x1)=1
% sig1= -1 =>  norm(x1)=1

% Z is the nullspace for matrix A.  size(Z)= t x (t-m)
% If Z is empty or not supplied, "null.m" is used (can be slow)
% samp is the number of samples (default is 1)

% OUTPUTS
% x1 is the sampled timeseries vector that satidfies the constraints A*x1=b
% varargout{1} = Z
% this is the "new nullspace " for the following step
%  null ([A; x1' ]) is obtained by rotating the old nullspace.  Z *Z_q

t= size(A,2);
idx = find(isnan(b));   % if b contain nan, ignore them
b(idx)=[]; A(idx,:)=[];

% compute minimum norm solution
xmn= lsqminnorm(A,b);

% if Z is not supplied, compute Z
if ~exist('Z', 'var') || isempty(Z)
 [q,~]=qr(A'); 
 Z= q(:, size(A,1)+1:end);
end

if ~exist('samp', 'var') || isempty(samp)
 samp=1; 
end

% compute 'd' to satisfy std or norm constraints
mean_xn= mean(xmn);
if sig1>=0 % constrain std
  d= sqrt((t-1)*sig1.^2+(t)*mean_xn^2-norm(xmn).^2); 

elseif sig1<0  % constrain the norm
  norm1=-sig1;
  d=sqrt(norm1^2- norm(xmn)^2);
end

if ~isreal(d)
  warning('std or norm constraints not satisfied. Proceed with caution');  
  x1=reshape(xmn,length(xmn),1);
  varargout{1}=[];
  return;
end

mz= length(Z(1,:));  % dimension of q

Nvar= normrnd(0,1,[mz,samp]);  % sample normal dist
q= Nvar./(vecnorm(Nvar,2,1));   % uniformly sample q

x1= reshape(xmn+d.*Z*q, [], samp); % x1= xmn +dZq

% update nullspace for next step
if nargout==2 && size(Z,2)>1  % if Z is not a single vector, compute new nullspace
    varargout{1}= null_expander(Z,q);
else
  varargout{1}=[];
end

end

