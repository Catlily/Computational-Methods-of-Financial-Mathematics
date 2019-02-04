a=0; b=2; % endpoints
y0=1; yN=3; % y values at endpoints
N=100; % number of points
h=(b-a)/N; % x step size
x=a+h:h:b-h; % set up vector of x points (interior points)
% set up the matrix and solve Ay=RHS for y
A=(2*h^2-2)*diag(ones(1,N-1)); % main diagonal elements
A=A+(1-3*h/2)*diag(ones(1,N-2),-1); % one below main diagonal
A=A+(1+3*h/2)*diag(ones(1,N-2),1); % one above main diagonal
RHS=2*h^2*x;
RHS(1)=RHS(1)-(1-3*h/2)*y0;
RHS(N-1)=RHS(N-1)-(1+3*h/2)*yN;
y=A\RHS';
plot(x,y)