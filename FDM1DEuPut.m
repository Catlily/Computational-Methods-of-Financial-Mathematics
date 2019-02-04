function P = FDM1DEuPut(S,E,r,T,sigma,q,Smin,Smax,M,N,theta)
% Finite difference method for European put
f7 = 1;
[Am,Ap,B,Svec] = GridU1D(r,q,T,sigma,Smin,Smax,M,N,theta);
V0 = max(E - Svec,0);
% Solve linear system for successive time steps
[L,U] = lu(Am);
for j = M-1:-1:0
V1 = V0;
V0(f7+0) = (E-Svec(f7+0))*exp(-r*(T-j*T/M)); % dt=T/M
b = Ap*V1(f7+(1:N-1)) + theta *B*V0(f7+[0 N])...
+ (1-theta)*B*V1(f7+[0 N]);
V0(f7+(1:N-1)) = U\(L\b);
end

for j = M-1:-1:0
V1 = V0;
V0(f7+0) = (E-Svec(f7+0))*exp(-r*(T-j*T/M)); % dt=T/M
b = Ap*V1(f7+(1:N-1)) + theta *B*V0(f7+[0 N])...
+ (1-theta)*B*V1(f7+[0 N]);
V0(f7+(1:N-1)) = U\(L\b);
end
P = interp1(Svec,V0,S,'spline');