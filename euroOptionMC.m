function [callPrice, putPrice, callCI, putCI] = euroOptionMC(n,S0,r,sigma,K,T,delta)
S = S0.* exp((r-0.5*sigma.^2)*T + sigma*sqrt(T)*randn(n,1));
%call
callCI = exp(-r*T)*max(K-S,0);
callPrice=mean(callCI)
callstd = sqrt(sum((callCI-callPrice).^2)/n-1);
calltemp = abs(norminv(delta/2))*callstd/sqrt(n);
callconf = [callPrice-calltemp,callPrice + calltemp]
%put
putCI = exp(-r*T)*max(S-K,0);
putPrice=mean(putCI)
putstd = sqrt(sum((putCI-putPrice).^2)/n-1);
puttemp = abs(norminv(delta/2))*putstd/sqrt(n);
putconf = [putPrice-puttemp,putPrice + puttemp]
end