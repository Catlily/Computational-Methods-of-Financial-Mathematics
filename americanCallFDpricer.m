function[callPriceMatrix,t,S] =americanCallFDpricer(Smax,sigma,r,K,T,spaceSteps,timeSteps,epsilon)
I=spaceSteps+1;
N=timeSteps+1;
ds=Smax/(I-1);
S=0:ds:Smax;
dt=T/(N-1);
t=0:dt:T;

maxSK=(max(S-K,0))';

j=[2:I-1];
jall=[1:I]'-1;
alphacoeff = 0.25*dt*( sigma^2*(jall.^2) - r*jall);
betacoeff = -0.5*dt*( sigma^2*(jall.^2) + r );
gammacoeff = 0.25*dt*( sigma^2*(jall.^2) + r*jall);

%Boundary conditions
iall=[1:N];
coeffmat1=diag(-alphacoeff(3:I-1),-1)+diag(1-betacoeff(2:I-1))+diag(-gammacoeff(2:I-2),1);
coeffmat2=diag(alphacoeff(3:I-1),-1)+diag(1+betacoeff(2:I-1))+diag(gammacoeff(2:I-2),1);

%matrix devision
c=zeros(I,N);
c(:,N)=maxSK;
C(1,iall)=0;
C(I,iall)=Smax-K*exp(-r*(N-iall)*dt);
cEdge=zeros(I-2,1);

for i=N-1:-1:1
    cEdge(end)=(c(I,i)+c(I,i+1))*gammacoeff(I-1);
    c(j,i)=(coeffmat1\(coeffmat2*c(j,i+1)+cEdge));
end

%SOR I set w=1.5
c=zeros(I,N);
c(:,N)=maxSK;
c(1,iall)=0;
c(I,iall)=(Smax-K)*exp(-r*(N-iall)*dt);
Cchange=zeros(I,1);
cold=zeros(I,1);
rc=zeros(I,1);
cEdge=zeros(I-2,1);
maxcount=50;
w=1.5;
for i=N-1:-1:1
    c(j,i)=c(j,i+1);
    cEdge(end)=(c(I,i+1))*gammacoeff(I-1);
    rc(j)=(coeffmat2*c(j,i+1)+cEdge);
    Cchangenorm=1;
    counter=0;
    while (Cchangenorm>epsilon)&&(counter<maxcount)
        for m=2:I-1
            cold(m)=C(m,i);
            c(m,i)=max(maxSK(m),(c(m,i)+(w.*(1./(1-betacoeff(m))).*(rc(m)+alphacoeff(m).*c((m-1),i)-(1-betacoeff(m)).*c(m,i)+gammacoeff(m).*c(m+1,i)))));
        end
        Cchange(j)=c(j,i)-cold(j);
        Cchangenorm=(norm(Cchange));
        counter=counter+1;
    end
end
callPriceMatrix=c;
surf(t,S,callPriceMatrix);
shading interp;
end

