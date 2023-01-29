function [beta1, flucts]=dfa(vt,intervals)

%   [beta1, flucts]=dfa_wrapper(ts, intervals)
% computes windowed fluctutations ('flucts') in each window and 
% fits a power-law with an exponenet 'beta1'. The function is heavily
% optimized for computing the fluctuations for 'n' timeseries parallely
% 
%% INPUTS 
% 'v0'                   -> n x t timeseries 
% 'intervals'         -> 1 x d vector of d time windows over which
%                              fluctuations are computed

%% OUTPUTS
% 'beta1'             -> n x 1 vector of power-law slope for each 1 x t
%                             timeseries
% 'flucts'             -> n x length(intervals)  of windowed fluctuations in each window 

[n,t]=size(vt);
intervals= reshape(intervals, [],1);
flucts= zeros(n, length(intervals));

parfor dd = 1:length(intervals)
  twin= intervals(dd); 
nd=floor(t/twin); % number of windows
N1=nd*twin;

csum=cumsum(vt(:,1:N1),2);  
yc= (csum-csum(:,end));
% compute trend to be subtracted
xm= ([[1:twin]', ones(twin,1)]);
xpmat= compute_proj_matrix(xm,twin);  % x'*pinv(x)'
yct= reshape(yc,n,twin,[]);
% multiply xct for each node using pagemtimes
Yn= reshape(pagemtimes(yct,xpmat), n,N1); 
flucts(:,dd)=sqrt(sum((yc-Yn).^2,2)/N1);
end

% compute slope in power-law
xx= [log(intervals),ones(size(intervals)) ]; 
yy= log(flucts); 
 temp=pinv(xx)*yy'; 
beta1=temp(1, :); 

end
