%@param p  natural coordinate in fine element
%@param pf natural coordinates of fine element
%          vertices in coarse element
function u=embedDisp(p, qf, qc, pf, Xf, Xc)
  %compute fine displacement
  Nf=bilinear(p);
end