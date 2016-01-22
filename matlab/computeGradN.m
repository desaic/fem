function gradN = computeGradN()
gradN = zeros(8,8,3);
for qq= 1:8
  for ii=1:8
    gradN(qq,ii,:) = shapeFunGrad(qq,ii);
  end
end
end