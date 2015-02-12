%@param p  natural coordinate in fine element
%@param pf natural coordinates of fine element
%          vertices in coarse element
function u=embedDisp(p, qf, qc, pf, Xf, Xc)
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
  Fc = defGrad(pc,qc,Xc);
  Nc = bilinear(pc);
  %coarse displacement
  uc = zeros(1,2);
  for ii = 1:4
      uc = uc + Nc(ii)*qc(ii,:);
  end
  u = uc + Fc*uf;
end