function val= null_func(x0,n)
%% see analytical derivation of nullspace (for single constraint)
% methods 

a=x0(1); 
b=x0(2); 
c=x0(3); 

val(1)= b-a*(n-2)-c; 
val(2)= c^2+a^2*(n-2)+b^2-1; 
val(3)= c^2-2*a*b+(n-3)*a^2;
end