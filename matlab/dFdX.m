%@param p natural coord
function dF = dFdX(p, dX, x, X)
  dF = zeros(2,2);
  for ii =1:4
    dN = bilinearGrad2(ii, p, dX, X);
    u = x(ii,:) - X(ii,:);
    dF = dF + u'*dN;
  end
end