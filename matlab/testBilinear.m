%test extrapolation behavior of bilinear weights
N = 20;
v = linspace(-1.5,1.5,N);
npt = N^2;
%displacements
q = zeros(4,2);
q(4, :) = [0.1 0.2];
q(1,:) = [0.2 0.1];
q(2,:) = [0 -0.1];
u = zeros(npt,2);
x = zeros(npt,2);
w = zeros(npt,4);
idx = 1;
square=[0 0;1 0 ; 0 1; 1 1];
for ii = 1:N
    for jj = 1:N
      x(idx, :) = [v(ii) v(jj)];
      w = bilinear(x(idx,:));
     % u(idx,:) = x(idx,:);
      for kk = 1:length(w);
    %    u(idx,:) = u(idx,:) + w(kk) * q(kk,:);
      u(idx,:) = u(idx,:) + w(kk) * (q(kk,:)+square(kk,:));
      end
      idx = idx + 1;
    end
end

scatter(u(:,1),u(:,2));

testp = [0.1 -0.2];
pn = natCoord(testp,square);
dN = zeros(4,2);
for ii = 1:4
  dN(ii,:) = bilinearGrad(ii, pn, square);
end
numDN=zeros(4,2);
h = 0.001;
for dim = 1:2
  testp(dim) = testp(dim) + h;
  pn = natCoord(testp,square);
  Np = bilinear(pn);

  testp(dim) = testp(dim) - 2*h;
  pn = natCoord(testp,square);
  Nm = bilinear(pn);
  
  numDN(:,dim) = (0.5/h)*(Np-Nm);
  
  testp(dim) = testp(dim) + h;
end
pn = natCoord(testp,square);
dN
numDN
q'*dN+eye(2)
F=defGrad(pn,square+q,square)