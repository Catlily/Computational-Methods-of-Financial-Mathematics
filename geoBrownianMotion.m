function[S]=geoBrownianMotion(mu,sigma,N,T,S0)
deltaT=T/N;
S=zeros(N+1,1);
S(1)=S0;
Z=randn(N+1,1);
for i = 2:N+1
S(i) = S(i-1)*exp((mu-0.5*sigma^2)*deltaT+sigma*sqrt(deltaT)*Z(i));
end
plot(S)
end


