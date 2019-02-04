function [call]= option_price_call_american_fin_diff_implicit(S,X,r,sigma,t,no_s_steps, no_t_steps,q) 
%   S: Spot price of underlying
%   X: Strike Price
%   r: interest rate
%   sigma: volatility
%	t: time to maturity
%	no_s_steps: steps in s
%	no_t_steps: steps in t
%   q: Dividend

if nargin<7, error('Not enough input arguments. At least 7 inputs'); end
% If no dividend entered, then dividend = 0;

if nargin==7
    q=0;
end

sigma_sqr = sigma*sigma;
%M=0; % no_S_steps should be even number..
%if (mod(no_s_steps,2)==1) % not even number..
 %   M=no_s_steps+1; 
 %else 
    M=no_s_steps;
    %end
ds=S/M;
frac=(S-M*ds)/ds;
M=3*M;
dt=t/no_t_steps;
kmax=no_t_steps;
w=1.5; % Relaxation Parameter (Wilmott, p 697)

max_iter=100; % Max Iteration..
acc=0.00000001; % Accuracy

Vold=zeros(M+1,1);
Vnew=zeros(M+1,1);
g=zeros(M+1,1);

% Initial Conditions for Call Payoff.. Similarly Put Payoff also can be set..
for i=1:M+1
    temp=(i-1)*ds;
    if temp > X
        Vold(i) = temp-X;
    end
    g(i)=Vold(i);
end

% Triagonal elements..    
a=zeros(M+1,1);
b=zeros(M+1,1);
c=zeros(M+1,1);

% Triagonal Elements for constant sigma, r, q and dt..
for j=2:M
    a(j) = 0.5*j*dt*(r-q+sigma_sqr*j);
	b(j) = 1.0+dt*(r+sigma_sqr*j*j);
	c(j) = 0.5*j*dt*(q-r+sigma_sqr*j);
end

% Time stepping..
for j=1:kmax-1
    tj=j*dt;
    % Initial Guess
    for i=2:M
        Vnew(i)=Vold(i);
    end
    k=0;
    Vnew(1)=0;
    Vnew(M+1)=ds*M*exp(-q*tj)-X*exp(-r*tj);
    norm=acc+1;
    while (k < max_iter) & (norm > acc)
        norm=0;
        for i=2:M
            y=(Vold(i)+a(i)*Vnew(i+1)+c(i)*Vnew(i-1))/b(i);
            temp=y-Vnew(i);
            norm=norm+temp*temp;
            Vnew(i)=Vnew(i)+w*temp;
            if Vnew(i)<g(i)
                Vnew(i)=g(i);
            end
            if i==M/3
%                [Vnew(i) y temp]
            end
        end
        k=k+1; % Increase the iteration index..
    end
    for i=2:M
        Vold(i)=Vnew(i);
    end
    Vold
    plot(Vold)
end