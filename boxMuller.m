function[Z] = boxMuller(mu,sigma,N)
u1=rand(1,N); 
u2=rand(1,N); 
val1=sqrt(-2*log(u1)).*cos(2*pi*u2);
val2=sqrt(-2*log(u1)).*sin(2*pi*u2); 
x=mu+sigma*val1;
y=mu+sigma*val2;
Z=x'
end