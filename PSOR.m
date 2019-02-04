function x1 = PSOR(x1,A,b,payoff,tol,maxit)
% Projected SOR for Ax = b
it = 0; 
converged = 0;
N = length(x1);
while ~converged
x0 = x1;
for i = 1:N-1
x1(i) = 0.5*(b(i)-A(i,:)*x1)/A(i,i) + x1(i);
x1(i) = max(payoff(i), x1(i));
end
converged = all((abs(x1-x0)./(abs(x0)+1)) < tol );
it=it+1; 
if it>maxit
    error('Maxit reached')
end
end