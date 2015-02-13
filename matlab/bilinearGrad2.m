function dN = bilinearGrad2(ii, p, dX, X)
sw = [-1, -1;
      -1,  1;
       1, -1;
       1,  1];
size = 1/(X(4,1) - X(1,1))^2;
dN=zeros(1,2);
dN(1,1)=size*sw(ii,1)*sw(ii,2)*dX(2);
dN(1,2)=size*sw(ii,2)*sw(ii,1)*dX(1);
end