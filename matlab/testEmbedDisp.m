Xf=[0.5 0.5;
    0.5 1;
    1  0.5;
    1 1];
Xc=[0 0;
    0 1;
    1 0;
    1 1];
qf=zeros(4,2);
qc=zeros(4,2);
qc(1,:) = [0.2 0.1];
qc(2,:) = [0 -0.1];
qc(4, :) = [0.1 0.2];

%qf(1,:) = [0.05 0.1];
%qf(3,:) = [-0.05 -0.05];
%qf(4, :) = [0.1 0.1];

N = 20;
v = linspace(-1, 1, N);
npt = N^2;

u = zeros(npt,2);
p = zeros(npt,2);
idx = 1;
square=[0 0;1 0 ; 0 1; 1 1];
for ii = 1:N
    for jj = 1:N
      p(idx, :) = [v(ii) v(jj)];
      w = bilinear(p(idx,:));
      Xp=w'*Xf;
      u(idx,:) = embedDisp(p(idx,:), qf, qc, Xf, Xc);
      u(idx,:) = u(idx,:) + Xp;
      idx = idx + 1;
    end
end

scatter(u(:,1),u(:,2));
