function p = natCoord(x, X)
  len = X(4,1) - X(1,1);
  p = (2/len)*(x-X(1,:))-[1 1];
end