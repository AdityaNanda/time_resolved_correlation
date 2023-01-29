function xxp= compute_proj_matrix(x,twin)
%% This script computes x * pinv(x)

% Given x =[1:twin; ones(twin,1)]';  
% pinv(x)*x is the projection matrix for detrending 

xxt(1,1)= sum((1:twin).^2); 
xxt(2,2)=twin; 
xxt(1,2)= sum(1:twin);
xxt(2,1)=xxt(1,2); 
xxp= x*(pinv(xxt)*x'); 







