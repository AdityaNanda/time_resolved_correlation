function [Fn, tau,llh_n, Fl,alph,llh_l]= ...
    get_avalanche_exponents(En)
% En is a single timeseries of avalanche events obtained from 
% makes_en.m or makes_en_3d.m 
% this function first computes starting and terminating 
% indices of each avalanche event in En
% An and Al are computed for each event 
% Fn and Fl are the frequencies corresponding to Fn and Fl  respectively 
% tau and alph are the exponents 
% llh_n and llh_l are the liklihoods 

if ~nnz(En)
  emp= nan; 
  Fn=emp; 
  tau=emp; 
  llh_n=emp; 
  Fl=emp; 
  alph=emp; 
  llh_l=emp; 
  return;
end
  
[An ,Al]= get_An_and_Al(En);

Fn=get_distribution(An);
Fl=get_distribution(Al);

[~, tau, ~,~,~,~,~, llh_n] = ...
    powerexp_fit(Fn, 1, 1);
 [~, alph,~,~,~,~,~,  llh_l] = ...
     powerexp_fit(Fl, 1,1);

end