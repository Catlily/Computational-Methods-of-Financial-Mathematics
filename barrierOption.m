function[callPrice,putPrice,callCI,putCI] =barrierOption(n,S0,r,sigma,K,b,T,M,upDown,inOut,delta)
callPayoffs = zeros(n,1);
putPayoffs = zeros(n,1);
deltaT=T/M;
S=zeros(M+1,1);
S(1)=S0;
mu=r;
Z=randn(M+1,1);

if (upDown~='D')&(upDown~='d')&(upDown~='U')&(upDown~='u')
error('Wrong upDown input.')
end

if (inOut~= 'I')&(inOut~='i')&(inOut~='O')&(inOut~='o')
error('Wrong inout input.')
end

for i = 2:M+1
S(i) = S(i-1)*exp((mu-0.5*sigma^2)*deltaT+sigma*sqrt(deltaT)*Z(i));
callPayoffs(i) = max(0,S(M)-K);
putPayoffs(i) = max(0,K-S(M));

hit = 0;
if (max(S)==b)&(upDown=='u'- upDown=='U')
hit = 1;
elseif (min(S)==b)&&(upDown=='d'- upDown=='D')
hit = 1;
end

if (inOut=='i')-(inOut=='I')
if hit == 1
Indicator = 1;
else
Indicator = 0;
end
end

if (inOut =='o')-(inOut =='O')
if hit == 1
Indicator = 0;
else
Indicator = 1;
end
end

callPayoffs(i) = callPayoffs(i)*Indicator;
putPayoffs(i) = putPayoffs(i)*Indicator;
end
callPrice = exp(-r*T) * mean(callPayoffs)
putPrice = exp(-r*T) * mean(putPayoffs)
z_Score = norminv(1-delta/2,0,1);
callStDev = callPayoffs - callPrice;
callStDev = 1/(n-1)* sum(callStDev.^2);
callCI = [callPrice - z_Score*callStDev/sqrt(n),callPrice + z_Score*callStDev/sqrt(n)]
putStDev = putPayoffs - putPrice;
putStDev = 1/(n-1)* sum( putStDev.^2);
putCI = [putPrice - z_Score*putStDev/sqrt(n),putPrice + z_Score*putStDev/sqrt(n)]