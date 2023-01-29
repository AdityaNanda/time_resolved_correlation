function v1= timesum_sampler(gs,tp, n,samp)

% this script can sample timeseries with 
% constrain global signal and variance 


close all; 
t= length(gs); 
 A= 1./n*ones(1,n); 

 temp =single_sum_sampler_3D(A,0,samp*t,1);
 v1= reshape(temp, [n,t,samp]); 
 v1= v1.*tp; 
 v1=v1+gs; 
end
      