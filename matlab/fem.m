function fem(e,x)
  size = x(e(1,8),1) - x(e(1,1),1);
  gradN = computeGradN();
  gradN = (1.0/size)*gradN;
  K = getStiffness(e,x);
end

function defGrad()

end

function K = getStiffness(e,x)

end

function f = getForce(e,x)

end