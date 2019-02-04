function [Am,Ap,B,S] = GridU1D(r,q,T,sigma, Smin,Smax,M,N,theta)
% Uniform grid in original variables
f7 = 1;
dt = T/M;
S = linspace(Smin,Smax,N+1)';
dS = (Smax - Smin) / N;
a = (dt*sigma^2*S(f7+(1:N-1)).^2) ./ (2*dS^2);
b = (dt*(r-q)*S(f7+(1:N-1))) / (2*dS);
d = a - b;
m = - 2*a - dt*r;
u = a + b;
P = spdiags([d m u], 0:2, N-1, N+1);
Am = speye(N-1) - theta *P(:,f7+(1:N-1));
Ap = speye(N-1) + (1-theta)*P(:,f7+(1:N-1));
B = P(:,f7+[0 N]);
end