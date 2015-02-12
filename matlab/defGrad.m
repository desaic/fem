%@param p natural coord
function F = defGrad(p, x, X)
  F = eye(2);
  for ii =1:4
    dN = bilinearGrad(ii,p,X);
    u = x(ii,:) - X(ii,:);
    F = F + u'*dN;
  end
end