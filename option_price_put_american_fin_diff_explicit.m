function [put]= option_price_put_american_fin_diff_explicit(S,X,r,sigma,t,no_s_steps, no_t_steps) 
%   S: Spot price of underlying
%   X: Strike Price
%   r: interest rate
%   sigma: volatility
%	t: time to maturity
%	no_s_steps: steps in s
%	no_t_steps: steps in t

sigma_sqr = sigma*sigma;
M=0; % no_S_steps should be even number..
if (mod(no_s_steps,2)==1) % not even number..
    M=no_s_steps+1; 
else 
    M=no_s_steps;
end
delta_S = 2*S/M;
S_values=zeros(M+2,1);
for m=1:M+1
    S_values(m) = m*delta_S;
end
N=no_t_steps;
delta_t = t/N;
    
a=zeros(M+1,1);
b=zeros(M+1,1);
c=zeros(M+1,1);
f_next=zeros(M+2,1);
f=zeros(M+2,1);

r1=1.0/(1.0+r*delta_t);
r2=delta_t/(1.0+r*delta_t);
    
for j=2:M
    a(j) = r2*0.5*j*(-r+sigma_sqr*j);
	b(j) = r1*(1.0-sigma_sqr*j*j*delta_t);
	c(j) = r2*0.5*j*(r+sigma_sqr*j);
end
    
for m=1:M+1
    f_next(m)=max(0.0,X-S_values(m));
end

for tt=N:-1:1
    f(1)=X;
	for m=2:M
        f(m)=a(m)*f_next(m-1)+b(m)*f_next(m)+c(m)*f_next(m+1);
        f(m)=max(f(m),X-S_values(m)); % check for exercise
    end
	f(M+1) = 0;
	for m=1:M+1
        f_next(m) = f(m);
    end
    f_next
    plot(f_next)
end