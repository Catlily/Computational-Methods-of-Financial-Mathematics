function[putPriceMatrix,t,S] =americanPutFDpricer(Smax,sigma,r,K,T,spaceSteps,timeSteps,epsilon)
I=spaceSteps+1;
N=timeSteps+1;
ds=Smax/(I-1);
S=0:ds:Smax
dt=T/(N-1);
t=0:dt:T;

maxKS=(max(K-S,0))';

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
p=zeros(I,N);
p(:,N)=maxKS;
p(I,iall)=0;
p(1,iall)=K;
pEdge=zeros(I-2,1);

for i=N-1:-1:1
    pEdge(1)=(p(1,i)+p(1,i+1))*alphacoeff(2);
    p(j,i)=(coeffmat1\(coeffmat2*p(j,i+1)+pEdge));
end
%SOR I set w=1.5
p=zeros(I,N);
p(:,N)=maxKS;
p(I,iall)=0;
p(1,iall)=K;
Pchange=zeros(I,1);
pold=zeros(I,1);
rp=zeros(I,1);
pEdge=zeros(I-2,1);
w=1.5;
for i=N-1:-1:1
    p(j,i)=p(j,i+1);
    pEdge(1)=(p(1,i+1))*alphacoeff(2);
    rp(j)=(coeffmat2*p(j,i+1)+pEdge);
    Pchangenorm=1;
    counter=0;
    while Pchangenorm>epsilon
        for m=2:I-1
            pold(m)=p(m,i);
            p(m,i)=max(maxKS(m),(p(m,i)+(w.*(1./(1-betacoeff(m))).*(rp(m)+alphacoeff(m).*p((m-1),i)-(1-betacoeff(m)).*p(m,i)+gammacoeff(m).*p(m+1,i)))));
        end
        Pchange(j)=p(j,i)-pold(j);
        Pchangenorm=(norm(Pchange));
        counter=counter+1;
    end
end
putPriceMatrix=p;
surf(t,S,putPriceMatrix);
shading interp;
end

