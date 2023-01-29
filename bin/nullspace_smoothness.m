function [v1]= nullspace_smoothness(nb,varargin)

% This function builds surrogate timeseries with 
% constrained time-resolved smoothness between neighbouring timepoints 
% and constrained mean and variance

%% INPUTS

% 'v0' - > empirical timeseries matrix ( size - ->  n x t)
%  n is number of ROIs. t is number of timepoints 

%nbvec - > is a vector of indices for time neighbours,
% correlations to which are constrained.  
%  If nbvec = [1] (typical case for canonical ), 
%           then correlations between adjacent timepoints are preserved.  
% If nbvec= [1,2,3], 
%            preserves correlations to the 3 nearest timepoints and so on
% nbvec should be contiguous.  
% That is, inputs of the form nbvec =[1,6, 9] are not considered 
% 
% nbvec can start with non-adjacent timepoints.  
% nbvec= [6,7,8] is legal 
 %            preserves correlations to 6th , 7th and 8th timepoints (closer timepoints (1,2,3 etc )are 'free' or unconstrained))
 %  f0 is the set of (max(nbvec) x t) matrix of correlations to be satisfied
% that is, if f0 is a 3 x 1000 matrix and (nbvec= [1,2,3])
%  the three adjacent correlations across 1000 timepoints will be
%  constrained.  

%% OUTPUTS 
% v1 is an n x t surrogate timeseries. 
% v1 is z-scored along time (each timepoint has mean 0 and std 1)
% v1 satisfies the correlations as prescribed in nbvec
base= mfilename('fullpath');
addpath([base(1:end-31), 'core_files/']);

check_cont_int= @(nb) all(diff(nb) <=1); 

p= inputParser; 
addRequired(p, 'nb', check_cont_int); 
addOptional(p,'data', []); 
addOptional(p,'similarity', []); 
addOptional(p, 'n', []); 
parse(p,nb,varargin{:}); 

f0= p.Results.similarity; 
v0= p.Results.data; 
nbvec= p.Results.nb; 
n= p.Results.n; 

if ~isempty(v0)
[n,t]= size(v0);
tmean= mean(v0,1); 
tstd= std(v0,[],1); 
% compute correlations using corr_kdiags.m 
f0= zeros(length(nbvec), size(v0,2)); 

for jj=1:length(nbvec)
f0(jj,:)= corr_kdiags(v0,nbvec(jj)); 
end

else  % f0 is provided 
    
    [m,t]= size(f0);
    if isempty(n)
        error('Number of regions must be provided with smoothness');
    end
    if  m~=length(nbvec)
        error('Enough smoothness values are not provided'); 
    end
tmean= zeros(1,t); 
tstd= ones(1,t); 
end
nbmin= min(nbvec);  % closest correlation constraint
seed= rand(n, nbmin);  % sample only enough seed columns to start sampling
seed= zscore(seed,[],1);   % sample zscored seed vector 

v1= zeros(n,t); % initialize output
cont= zeros(n,1);  
v1(:,1:nbmin)= seed;  % the first min(nbvec) timepoints are not constrained so use seed
samp=1;  % number of samples 
xsig=1;  % set std=1 

for tt=nbmin+1:t
    
  nb_ind= tt-nbvec; % index of constrained timepoints
nb_ind(nb_ind<1)=[];   % necessary untill tt= nbmax+1 
        ts= v1(:, nb_ind);  %  activity at previous m timepoints 
      c1= f0(:, tt);  %  correlations to the previous m timepoints
 c1(isnan(c1))=[];
      % compute b  and A ( setting up equations as Ax=b)
  b= [(n-1).*c1;0];   % (n-1 ) is a scaling parameter (see derivation in Methods)
  A= [ts'; ones(1,n)];  % coefficients 

[v1(:,tt)]= timeseries_sampler(A,b,samp,xsig);
  
end

v1= (v1.*tstd)+tmean;

end








