function [randValues] = multiVarNormal(mu, Sigma, N)
d=length(mu);
u1=rand(d,N);
u2=rand(d,N);
val1=sqrt(-2*log(u1)).*cos(2*pi*u2);
n = length(Sigma);
L = zeros(n,n);
%Check if the matrix is symmetric
if all( all( Sigma == Sigma' ) ) & min( eig( Sigma ) ) > 0
 for i=1:d
    L(i,i) = sqrt( Sigma(i, i) - L(i, :)*L(i, :)' );
   for j=(i + 1):d
        L(j, i) = ( Sigma(j, i) - L(i, :)*L(j, :)' )/L(i, i)
   end
 end
randValues=mu+L'*val1
else 
    disp('sigma is not symmetric and positive definite')
end
end