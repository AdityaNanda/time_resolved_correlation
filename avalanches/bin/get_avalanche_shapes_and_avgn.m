function [sh_clps, avg_n, lrange, signuz_inv, varargout]...
    = get_avalanche_shapes_and_avgn(En, lind)

% This function computes the avalanche shape collapse profiles 
% by averaging over all avalanches of the same length 
 % sh_clps is a cell array of collapsed avalanches 
 % lrange is a vector of avalanache lengths 
 % avg_n is a vector of average avalanche size for each avalanche length in
 % lrange
 % thus lrange is the same size as avg_n
 % signuz_inv is the slope for the power law relating lrange and avg_n
 % sh_clps_num is the number of "collapsed " avalanches 
 % Thus if % sh_clps_num(20) = 200 then 200 avalanches were collapsed to 
 % obtain the shape profile in sh_clps{20}
 
if ~nnz(En)
  emp=nan; 
  sh_clps=emp; 
  avg_n=emp; 
  lrange=emp; 
  signuz_inv=emp; 
  varargout{1}=emp; 
  varargout{2}=emp;
  return;
end
[An ,Al, A_sta, A_fin]= get_An_and_Al(En); 

 lrange= sort(unique(Al)); 
avg_n=lrange; 
raw_shapes= cell(length(lrange),1); 
sh_clps=raw_shapes; 
sh_clps_num=lrange; 
for ll=1:length(lrange)
    ix= find(Al==lrange(ll)); 
    avg_n(ll)=  mean(An(ix));  
   
    for jj=1:length(ix) % iterate over all avalanches of same length
    raw_shapes{ll}(jj,:)= En(A_sta(ix(jj)):A_fin(ix(jj)));
    sh_clps_num(ll)= uint16(size(raw_shapes{ll},1)); % number of avalanches being averaged over
    sh_clps{ll}= single(mean(raw_shapes{ll},1)); 
    end
%      close all; plot(sh_clps{ll}'./lrange(ll)); 
%      pause(0.1)
end
if ~exist('lind', 'var') || isempty(lind) || lind>length(lrange)
  lind=length(lrange); 
end

 [signuz_inv, intrcpt, signuz_fit]=power_law(lrange(1:lind),avg_n(1:lind)); 
 
 switch nargout-4
     case 1
 varargout{1}= sh_clps_num; 
     case 2
        varargout{1}= sh_clps_num; 
 varargout{2}= signuz_fit; 
 end
end

