%European option
%ACall and APut are  matrices that result from the spatial discretization 
function[ACall,APut] = spatialCoeffsBS(r,sigma,spacePoints,xMax,xMin)
%Spacial discretization 
I=spacePoints-1;
dx=(xMax-xMin)/I;
%parameters
xi=xMin:dx:(xMax-2*dx);
alpha=1/2*r*(xMin+dx*xi)/dx+1/2*sigma^2*xi.^2/dx^2;
beta=-1/2*r*xi/dx+1/2*sigma^2*xi.^2/dx^2;
gamma=r+sigma^2*xi.^2/dx^2;
gamma=[r,gamma,r];
%Matrix
main=-diag(gamma);
uper=diag([0,alpha],1);
lower=diag([beta,0],-1);
APut=main+uper+lower
ACall=APut(2:I,2:I)
end