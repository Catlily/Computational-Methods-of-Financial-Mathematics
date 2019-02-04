%European option
%ACall and APut are  matrices that result from the spatial discretization 
function[alpha,beta,gamma] = euro(r,sigma,spacePoints,xMax,xMin)
%Spacial discretization 
I=spacePoints-1;
dx=(xMax-xMin)/I;
%Iteration
xi=xMin:dx:xMax
alpha=1/2*r*(xMin+dx*xi)/dx+1/2*sigma^2*xi.^2/dx^2
beta=-1/2*r*xi/dx+1/2*sigma^2*xi.^2/dx^2
gamma=r+sigma^2*xi.^2/dx^2
end