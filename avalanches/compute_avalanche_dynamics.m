function [Fn,tau, Fl,alph,...
 shape_cllps, avg_n_given_l, signuz_inv ]= compute_avalanche_dynamics(v)

% compute avalanche dynamics
% v is a n x t LFP recording 

%% OUTPUTS

% Fn                       -> probability density function of avalanche sizes
% tau                      -> powerlaw exponent fit to Fn
% Fl                        -> probability density func of avlanche duration 
% alph                     -> power-law exponent fit to Fl
% shape_cllps          -> cell array containing averaged shapes of
%                                avalanches
% avg_n_given_l      -> Average avlaanche size given a certain duration
% signuz_inv            -> Power-law exponent fit to avalanche duration and
%                                avg_n_given_l

% convert to point-process
sd_scale=2.5; % std threshold
m0= mean(v,2);  
s0= std(v,[],2);
vt= int16(zeros(size(v))); % zeros
thres1= m0-sd_scale.*s0;  % this takes the negative LFP excursions
va=v;
vt(va<thres1)=true;
% compute avalanche events
En=sum(vt);
En([1:find(~En,1,'first') find(~En,1,'last'):end])=0;       %remove clipped events

if ~nnz(En)
    Fn=nan;
    tau=nan;
    llh_n=nan;
    Fl=nan;
    alph=nan;
    llh_l=nan;
    return;
end
% compute avalanche size and avalanche durations
En_ind=find(En);                                            %indices of event occurence
A_start=En_ind([true, En_ind(2:end)-1~=En_ind(1:end-1)]);     %commencement of avalanches
A_finish=En_ind([En_ind(1:end-1)+1~=En_ind(2:end), true]);     %termination of avalanches

An=zeros(size(A_start)); 
for g=1:length(A_start)
    An(g)=sum(En(A_start(g):A_finish(g)));                       %size of avalanches
    
end
Al=A_finish-A_start+1;        % duration of avalanches

%% compute probability distributions for An and Al
Fn=get_distribution(An);
Fl=get_distribution(Al);

tau = powerexp_fit(Fn, 1,1);
alph= powerexp_fit(Fl, 1,1);

%% compute shape collapse and signuz_inv
lrange= sort(unique(Al)); 
avg_n_given_l=lrange; 
raw_shapes= cell(length(lrange),1); 
shape_cllps=raw_shapes; 
sh_clps_num=lrange; 

for ll=1:length(lrange) % iterate over avalanche durations
    ix= find(Al==lrange(ll));   % all avalanches with length ll
    avg_n_given_l(ll)=  mean(An(ix));   % average size of all avalanches with length ll
    
    for jj=1:length(ix) % iterate over all avalanches of length ll
    raw_shapes{ll}(jj,:)= En(A_start(ix(jj)):A_finish(ix(jj)));
    end
    % number of avalanches being averaged over
    sh_clps_num(ll)= uint16(size(raw_shapes{ll},1)); 
    % compute shape collapse
    shape_cllps{ll}= mean(raw_shapes{ll},1);  
end
%% compute power-law between lrange and avg_n
lind= 12; 
x_soln= polyfit(log(lrange(1:lind)),log(avg_n_given_l(1:lind)),1);
signuz_inv= x_soln(1);  % intercept of powerlaw

 end

function F=get_distribution(X)
[F, ~ ,lim]= histi(int32(X));             %distribution and its range
F=[zeros(1,lim(1)-1) F];                %ensure distribution starts at 1
end


