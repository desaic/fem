function gradN = shapeFunGrad(qq, i)
sw=[-1.0,-1.0,-1.0;
 -1.0,-1.0, 1.0;
 -1.0, 1.0,-1.0;
 -1.0, 1.0, 1.0;
  1.0,-1.0,-1.0;
  1.0,-1.0, 1.0;
  1.0, 1.0,-1.0;
  1.0, 1.0, 1.0];
quadrature=[-0.57735, -0.57735, -0.57735;
 -0.57735, -0.57735,  0.57735;
 -0.57735,  0.57735, -0.57735;
 -0.57735,  0.57735,  0.57735;
  0.57735, -0.57735, -0.57735;
  0.57735, -0.57735,  0.57735;
  0.57735,  0.57735, -0.57735;
  0.57735,  0.57735,  0.57735];
p = quadrature(qq,:);
gradN = zeros(3,1);
gradN(1) = sw(i,1) * (1.0 + sw(i,2) * p(2)) * (1.0 + sw(i,3) * p(3));
gradN(2) = sw(i,2) * (1.0 + sw(i,1) * p(1)) * (1.0 + sw(i,3) * p(3));
gradN(3) = sw(i,3) * (1.0 + sw(i,1) * p(1)) * (1.0 + sw(i,2) * p(2));
gradN = 0.25*gradN;
end