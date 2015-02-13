%@param X reference coordinate
function u=embedDisp(X, qf, qc, Xf, Xc)
  %compute fine displacement
  pf=natCoord(X,Xf);
  Nf=bilinear(pf);
  uf=zeros(1,2);
  for ii = 1:4
    uf = uf + Nf(ii)*qf(ii,:);
  end
  %natural coordinate in coarse element
  pc = natCoord(X,Xc);
  Fc = defGrad(pc,qc+Xc,Xc);
  Nc = bilinear(pc);
  %coarse displacement
  uc = zeros(1,2);
  for ii = 1:4
      uc = uc + Nc(ii)*qc(ii,:);
  end
  u = uc + uf*Fc';
end