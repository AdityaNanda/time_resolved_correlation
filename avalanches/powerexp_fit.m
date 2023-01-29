function [alph, d,gamm,xmin, xmax, Px, Cx,...
  llhmax] = powerexp_fit(F, rmin, rmax, lx, tip)
%[d alph gamm xmin xmax P C] = powerexp_fit(F, rmin, rmax, lx, tip)
%
%   Fit power law with exponential cutoff
%   Inputs,
%       F,      input frequency histogram
%       rmin,   range of xmin
%       rmax,   range of xmax
%
%       e.g.    powerexp_fit(F, rmin, rmax, -0.1:0.001:0.1, 1.01:0.01:3);
%
%   Outputs,
%       d,      Kolmogorov-Smirnov statistic
%       alph,   exponent of the power law
%       gamm,   exponent of the exponential
%       xmin,   lower bound on distribution
%       xmax,   upper bound on distribution
%       P,      perfect probability distribution for perfect power-law

% P(alph, gamm) = x^(-alph) * exp(-gamm*x);
% sum( P(alph,gamm,xmin:end) ./ Konst(alph,gamm,xmin) ) = 1;
% sum( P(alph,gamm,xmin:xmax)./(Konst(alph,gamm,xmin)-Konst(alph,gamm,xmax+1)) ) = 1;

% preparation of reference tables
persistent Pag Konst Alph Gamm la lg
if ischar(F) && strcmp(F,'reset')
    mlx=rmin;                                               %this is not rmin in the standard sense
    Alph=1:0.01:5.9;
     Gamm=0; %-0.8:0.002:0.8;
     la = length(Alph);
    lg = length(Gamm);
    Pag = inf(la,lg,mlx);
    Konst = inf(la,lg,mlx);
    for i=1:la
        for j=1:lg
            Pag_ij     = ((1:mlx).^(-Alph(i))) .* exp(-Gamm(j)*(1:mlx));
            idx = find([Pag_ij inf]>Pag_ij(1),1) - 1;
            Pag_ij=Pag_ij(1:idx);
            Pag(i,j,1:idx) = Pag_ij;
            CPag_ij = cumsum(Pag_ij);
            Konst(i,j,1:idx) = CPag_ij(end) - [0 CPag_ij(1:end-1)];
        end
    end
    return
end

C=fliplr(cumsum(fliplr(F)));
if ~exist('lx','var')
    lx=length(F);
end
% lx=lx+1;
if ~exist('tip','var')
    tip='kuiper';
end

% get best alph and gamm for each xmin and xmax pair
% and get D for that pair
% for each xmin, the best alph and gamm maximize:
% Log-likelihood = -n*log(Konst) - alph*sum(ln(F)) - gamm*sum(F)
c = cumsum(log(1:lx).*F(1:lx));
SlogX = c(end) - [0 c(1:end-1)];                            %cumulative sum Ln(F) (F>=x, 1<=x<=rmax)
c = cumsum(   (1:lx).*F(1:lx));
SX    = c(end) - [0 c(1:end-1)];                            %cumulative sum (F) (F>=x, 1<=x<=rmax)

Ai = zeros(rmin,rmax);
Gi = zeros(rmin,rmax);
D  = zeros(rmin,rmax);
for g=1:rmin
    for h=1:rmax
        LLh=-inf(la,lg);
        for i=1:la
            for j=1:lg
                if ~isinf(Konst(i,j,lx-h+1))
                    LLh(i,j) = ...
                        - (    C(g) -     C(lx-h+1)) .* log(Konst(i,j,g)-Konst(i,j,lx-h+1)) ...
                        - (SlogX(g) - SlogX(lx-h+1)) .* Alph(i) ...
                        - (   SX(g) -    SX(lx-h+1)) .* Gamm(j);
                end
            end
        end
        
        [i j]=find(LLh==max(LLh(:)),1);
        Ai(g,h)=i;
        Gi(g,h)=j;
        llhmax= LLh(i,j);  % likelihood of best estimate
        
        CPagAiGi = cumsum(squeeze(Pag(i,j,g:lx-h)./(Konst(i,j,g)-Konst(i,j,lx-h+1)))).';
        Cemp = (C(g:lx-h)-C(lx-h+1))./(C(g)-C(lx-h+1));
        Cthe = (CPagAiGi(end) - [0 CPagAiGi(1:end-1)]);
        assert(isequal(size(Cemp), size(Cthe)), 'Cthe and Cemp must have the same size');
        
        if abs(Cthe(1)-1)<1e-10
            Cthe(1)=1;
        else
            error('The first value of the theoretical cumulative distribution is not 1.')
        end
        d = (Cemp-Cthe)./sqrt(Cthe.*(1-Cthe));
        if ~isreal(d)
            error('d contains an imaginary value.')
        end
        
        switch tip
            case 'ks'
                D(g,h) = max(abs(d));
            case 'kuiper'
                D(g,h) = max(d)+max(-d);
        end
        
    end
end

[imin imax]=find(D==min(D(:)),1);
i=Ai(imin,imax);
j=Gi(imin,imax);
d=D(imin,imax);
alph=Alph(i);
gamm=Gamm(j);
xmin=imin;
xmax=lx-imax;
Px = squeeze(Pag(i,j,1:xmax)./(Konst(i,j,xmin)-Konst(i,j,xmax+1))).';   %best model PDF
Cx = sum(Px) - [0 cumsum(Px)];    % best model CDF

