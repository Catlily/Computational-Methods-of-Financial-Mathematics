% LU
[L,U]=lu(coeffmat1);
for i=N-1:-1:1
    pEdge(1)=(p(1,i)+p(1,i+1))*alphacoeff(2);
    p(j,i)=(U\(L\(coeffmat2*p(j,i+1)+pEdge)));
end
% matrix devision
for i=N-1:-1:1
    pEdge(1)=(p(1,i)+p(1,i+1))*alphacoeff(2);
    p(j,i)=(coeffmat1\(coeffmat2*p(j,i+1)+pEdge));
end
