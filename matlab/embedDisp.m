%@param p  natural coordinate in fine element
function u=embedDisp(p, qf, qc, Xf, Xc)
  %compute fine displacement
  Nf=bilinear(p);
  uf=zeros(1,2);
  for ii = 1:4
    uf = uf + Nf(ii)*qf(ii,:);
  end
  %point in reference frame
  Xp=zeros(1,2);
  for ii = 1:4
    Xp = Xp + Nf(ii) * Xf(ii,:);
  end
  %natural coordinate in coarse element
  pc = natCoord(Xp,Xc);
  Fc = defGrad(pc,qc+Xc,Xc);
  Nc = bilinear(pc);
  %coarse displacement
  uc = zeros(1,2);
  for ii = 1:4
      uc = uc + Nc(ii)*qc(ii,:);
  end
  u = uc + uf*Fc';
end