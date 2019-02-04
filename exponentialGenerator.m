function[x] = exponentialGenerator(N,theta)
u=rand(N,1);
x=-theta*log(u);
end