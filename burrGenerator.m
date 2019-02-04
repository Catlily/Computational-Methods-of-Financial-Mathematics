function[x] = burrGenerator(N,c,k)
u=rand(N,1);
x=((1-u).^(-1/k)-1).^(1/c);
end