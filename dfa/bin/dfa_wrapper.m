
function [beta1, flucts]=dfa_wrapper(v0, intervals)

%   [beta1, flucts]=dfa_wrapper(ts, intervals)
% computes windowed fluctutations ('flucts') in each window and 
% fits a power-law with an exponenet 'bets'
% 
%% INPUTS 
% 'v0'                   -> n x t timeseries 
% 'intervals'         -> 1 x d vector of d time windows over which
%                              fluctuations are computed

%% OUTPUTS
% 'beta1'             -> n x 1 vector of power-law slope for each 1 x t
%                             timeseries
% 'flucts'             -> n x length(intervals)  of windowed fluctuations in each window 
                          
intervals= reshape(intervals, [],1);
flucts=dfa(v0,intervals);

xx= [log(intervals),ones(size(intervals)) ]; 
yy= log(flucts); 
 temp=pinv(xx)*yy'; 
beta1=temp(1, :);  % checked

end

