function[u] = linConGenerator(m,a,c,seed,N)
x(1)=mod(seed*a+c,m);
for i=1:N-1
    x(i+1)=mod(x(i)*a+c,m);
    u(i+1)=x(i+1)/m;
end
disp(x)