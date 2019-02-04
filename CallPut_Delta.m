function [C, Cdelta, P, Pdelta] = CallPut_Delta(S,K,r,sigma,tau,div)
% tau = time to expiry (T-t)
if nargin < 6
div = 0.0;
end
if tau > 0
d1 = (log(S/K) + (r + 0.5*sigma^2)*(tau)*ones(size(S))) / (sigma*sqrt(tau));
d2 = d1 - sigma*sqrt(tau);
N1 = 0.5*(1+erf(d1/sqrt(2))); N2 = 0.5*(1+erf(d2/sqrt(2)));
C = exp(-div*tau) * S.*N1-K*exp(-r*(tau))*N2; Cdelta = exp(-div*tau) * N1;
P = C + K*exp(-r*tau) - exp(-div*tau)*S; Pdelta = Cdelta - exp(-div*tau);
else
C = max(S-K,0); Cdelta = 0.5*(sign(S-K) + 1);
P = max(K-S,0); Pdelta = Cdelta - 1;
end