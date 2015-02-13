%@param p  natural coordinate in fine element
function F = embedDefGrad(X, qf, qc, Xf, Xc)
pf = natCoord(X, Xf);
pc = natCoord(X, Xc);
Ff = defGrad(pf, qf+Xf, Xf);
Fc = defGrad(pc, qc+Xc, Xc);
F  = Fc*Ff;
 %compute fine displacement
  pf=natCoord(X,Xf);
  Nf=bilinear(pf);
  uf=zeros(1,2);
  for ii = 1:4
    uf = uf + Nf(ii)*qf(ii,:);
  end
  
for dim = 1:2
  dX = zeros(1,2);
  dX(dim)=1;
  dF=dFdX(pc,dX,qc+Xc,Xc);
  F(:,dim) = F(:,dim) + dF * uf';
end
end