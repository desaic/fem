function M = lookupStructures(a, s)
sizea = size(a);
sizem = size(s,2);
nx = sqrt(sizem);
ny = nx;
M = zeros(sizea(1) * nx, sizea(2) * ny);
for ii = 1:sizea(1)
  for jj = 1:sizea(2)
    matIdx = a(ii,jj)+1;
    structure = s(matIdx, :);
    %1 2 3 4 -> 1 3
    %           2 4
    structure = reshape(structure, [nx ny]);
    M( ((ii-1)*ny + 1) : ii*ny, ((jj-1)*nx + 1) : jj*nx) = structure;
  end
end
end