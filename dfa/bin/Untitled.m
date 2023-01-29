
clc; 

v0 =rand(10,1e5); 
intervals= floor(logspace(1,4,100)); 

 [beta1, flucts]=dfa(v0, intervals); 

 