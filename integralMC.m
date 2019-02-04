function [est, actualError] = integralMC(n)
x=rand(1,n);
y=x.^3;
f=sym('x^3');
integralf=int(f,0,1);
est=sum(y)/n
actualError=integralf-est
end
