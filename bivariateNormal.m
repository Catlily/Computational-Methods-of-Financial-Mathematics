function[X] = bivariateNormal(mu,sigma)
%use eigenvalue to check
e=eig(sigma)
if sigma==sigma' & e(1)>0 & e(2)>0 
C=chol(sigma,'lower')
Z =boxMuller(0,1,2);
X =mu+C'*Z;
else 
    disp('sigma is not symmetric and positive definite')
end
end