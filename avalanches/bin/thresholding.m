% This script converts ecog potentials 
% into spikes 

function vt= thresholding(v,sd_scale,exc)

% thres= prctile(v',tile);  % nodewise prctile
% This returns a point-process timeseries given a sdancdard deviation 
%  only excursions in the negative potential are considered by default
% if exc is nonzero, the excursions in both directions are considered
if ~exist('exc', 'var') || isempty(exc)
  exc=0; % default value  --> only nLFPs will consider
end

m0= mean(v,2); 
s0= std(v,[],2); 
vt= int16(zeros(size(v))); % zeros  

if exc   % both positive and negative excursions from mean
  thres1= m0+sd_scale.*s0; 
  va= abs(v); 
  vt(va>thres1)=true;  
else
   thres1= m0-sd_scale.*s0; 
  va=v;  % this takes the negative excursions
  vt(va<thres1)=true;  
end


 
end