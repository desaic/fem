function [A, M43, U, W] = makePeriodBC(K, d)
l4 = length(d{4});
l3 = length(d{3});
nEdgeDof = (l4 - l3)/2;
nEdgeVert = floor(nEdgeDof/3);
nFaceDof = l3 - nEdgeDof;
i = zeros(l4,1);
j = zeros(l4,1);

for vert = 1:nEdgeVert
  for i1 = 1:3
      dof = 3*(vert-1)+i1;
      idx = 3*(dof-1) + 1;
      i(idx) = 9*(vert-1) + i1;
      i(idx+1) = 9*(vert-1) + i1 + 3;
      i(idx+2) = 9*(vert-1) + i1 + 6;
      j(idx:idx+2) = dof;
  end
end

for i1 = 1:nFaceDof
  idx = 3*nEdgeDof + i1;
  i(idx) = idx;
  j(idx) = nEdgeDof + i1;
end

v = ones(l4, 1);
M43 = sparse(i,j,v,l4,l3);
%M43*ones(l3,1)
%M43*(1:l3)'
K21 = K(d{2}, d{1});
K22 = K(d{2}, d{2});
K23 = K(d{2}, d{3});
K24 = K(d{2}, d{4});
K31 = K(d{3}, d{1});
K33 = K(d{3}, d{3});
K34 = K(d{3}, d{4});
K41 = K(d{4}, d{1});
K44 = K(d{4}, d{4});
A = [K22, K23 + K24*M43;
    (K23+K24*M43)', K33 + K34*M43 + (K34*M43)' + M43'*K44*M43];
U = [K21; K31 + M43'*K41];
W = [K24; K34 + M43'*K44];
end