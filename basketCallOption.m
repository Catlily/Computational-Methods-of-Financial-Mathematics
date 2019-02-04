function[optionPrice]=basketCallOption(n,S0,mu,Sigma,c,T,K,r,N)
d=length(mu);
A=zeros(d,d);
v=zeros(d,1);
if Sigma ~= Sigma
 error('Sigma is not symmetric') ;
end
if min(eig(Sigma))<= 0
 error('Sigma is not positive definite');
end
for j = 1:d
 for i = j:d
 v(i) = Sigma(i,j);
 for q = 1:(j-1)
 v(i) = v(i) - A(j,q)*A(i,q);
 end
 A(i,j) = v(i)/sqrt(v(j));
 end
end
deltaT = T/N;
Z = randn(d,N+1);
S = zeros(d,N+1);
S(:,1) = S0;

for j= 1:n
for i = 1:d
 for q = 1:N
 S(i,q+1) = S(i,q)*exp((mu(i)-0.5*Sigma(i,i))*deltaT+sqrt(deltaT)*(A(i,:)*Z(:,q+1)));
 end
end
callpayoffs = zeros(1,N);
for i = 1:K
callpayoffs(1,i)= max(sum(c.*S(:,i)- K),0);
optionPrice(j) = exp(-r*T)*mean(callpayoffs);
end
end
end