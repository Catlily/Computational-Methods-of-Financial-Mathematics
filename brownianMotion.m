function [W]=brownianMotion(mu, sigma, N, T, W0)
%Time vector
delta_t= (0:1:N)'/N; 
delta_t= delta_t*T;
%Running sum of N(0,1/N) variables
SN = [0; cumsum(randn(N,1))]/sqrt(N); 
%Standard Brownian motion
SBM = SN.*sqrt(delta_t);
W = W0+mu*delta_t+sigma*SBM;
plot(W)
end
