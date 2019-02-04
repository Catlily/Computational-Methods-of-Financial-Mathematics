function X = XGenerator(N)
%g(v)=1,Convenience distribution U(0,1)
%choose c=1/2
c=1/2;
i=1;
      while i <=N
           v = rand(1);
           u = rand(1);
           if u <= 20.*v.*(1-v).^3/c
               X(i)=v;
               i = i + 1;
           end
      end
       
end