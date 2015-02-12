%@param p point in natural coordinate
%@return a row vector of derivatives
function dN = bilinearGrad(ii, p, X)
sw = [-1, -1;
      -1,  1;
       1, -1;
       1,  1];
size = 0.5/(X(4,1) - X(1,1));
dN=zeros(1,2);
dN(1,1)=size*sw(ii,1)*(1+sw(ii,2)*p(2));
dN(1,2)=size*sw(ii,2)*(1+sw(ii,1)*p(1));
end