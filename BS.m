%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r=0.2; % Interest rate
sigma=0.25; % Volatility of the underlying
M=1600; % Number of time points
N=160; % Number of share price points
Smax=20; % Maximum share price considered
Smin=0; % Minimum share price considered
T=1.; % Maturation (expiry)of contract
E=10; % Exercise price of the underlying
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt=(T/M); % Time step
ds=(Smax-Smin)/N; % Price step
% Initializing the matrix of the option value
v=zeros(N,M);
% Initial conditions prescribed by the European Call payoff at expiry:V(S,T)=max(S-E,0);
v(1:N,1)=max((Smin+(0:N-1)*ds-E),zeros(1,N))';
% Boundary conditions prescribed by the European Call:
v(1,2:M)=zeros(M-1,1)'; % V(0,t)=0
v(N,2:M)=((N-1)*ds+Smin)-E*exp(-r*(1:M-1)*dt); % V(S,t)=S-Eexp[-r(T-t)] as S-infinite
% Determining the matrix coeficients of the explicit algorithm
aa=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)-r*(1:N-2))';
bb=1-dt*(sigma*sigma*(1:N-2).*(1:N-2)+r)';
cc=0.5*dt*(sigma*sigma*(1:N-2).*(1:N-2)+r*(1:N-2))';
% Implementing the explicit algorithm
for i=2:M,
v(2:N-1,i)=bb.*v(2:N-1,i-1)+cc.*v(3:N,i-1)+aa.*v(1:N-2,i-1);
end
% Reversal of the time components in the matrix as the solution of the BlackScholes
% equation was performed backwards
v=fliplr(v);
% Figure of the value of the option, V(S,t), as a function of S
% at three different times:t=0, T/2 and T (expiry).
figure(1)
plot(Smin+ds*(0:N-1),v(1:N,1)','r-',Smin+ds*(0:N-1),v(1:N,round(M/2))','g-',Smin+ds*(0:N-1),v(1:N,M)','b-');xlabel('S');
ylabel('V(S,t)');
title('European Call Option within the Explicit Method');
% Figure of the Value of the option, V(S,t)
figure(2)
mesh(Smin+ds*(0:N-1),dt*(0:M-1),v(1:N,1:M)')
title('European Call Option value, V(S,t), within the Explicit Method')
xlabel('S')
ylabel('t')