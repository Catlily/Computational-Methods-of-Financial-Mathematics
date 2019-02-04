clear ;
N=10; % number of steps used
a=0; b=1; % endpoints of solution
dx=(b-a)/N; % step size
% save x and y values in vectors so we can plot solution
x(1)=a; y(1)=1; % initial x and y values
for i=1:N
f=y(i); % calculate the function value at (xi,yi)
y(i+1)=y(i)+dx*f;
x(i+1)=x(i)+dx;
end
plot(x,y,'b-')
hold on
exact=exp(x);
plot(x,exact,'r--');
hold off
error=abs(exact(N+1)-y(N+1))
    