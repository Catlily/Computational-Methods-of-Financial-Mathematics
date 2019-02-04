function [S]=geoBrownianMotionTerm(sigma,S0,T)
if T>3
    error('T should be smaller than 3')
end
data=xlsread('Term_Structure.xlsx');
delta_t=3/12;
N=T/delta_t;
end
