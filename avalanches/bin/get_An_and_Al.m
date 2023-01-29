function [An ,Al, A_sta, A_fin]= get_An_and_Al(En)

En_ind=find(En);                                            %indices of event occurence
A_sta=En_ind([true, En_ind(2:end)-1~=En_ind(1:end-1)]);     %commencement of avalanches
A_fin=En_ind([En_ind(1:end-1)+1~=En_ind(2:end), true]);     %termination of avalanches

An=zeros(size(A_sta));
for g=1:length(A_sta)
    An(g)=sum(En(A_sta(g):A_fin(g)));                       %size of avalanches
    
end
 Al=A_fin-A_sta+1;        



