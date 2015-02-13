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
qc(1,:) = [0.0 0.2];
qc(2,:) = [0.2 0.0];
qc(4, :) = [0.1 -0.2];

qf(1,:) = [0.0 0.05];
qf(2,:) = [-0.05 0.0];
qf(4, :) = [0.1 -0.1];

scale = 0.4;
Xf = scale *Xf;
Xc = scale *Xc;
qf = scale * qf;
qc = scale * qc;

N = 20;
v = scale * linspace(0.5, 1, N);
npt = N^2;

u = zeros(npt,2);
p = zeros(npt,2);
idx = 1;
for ii = 1:N
    for jj = 1:N
      p(idx, :) = [v(ii) v(jj)];
      u(idx,:) = embedDisp(p(idx,:), qf, qc, Xf, Xc);
      u(idx,:) = u(idx,:) + p(idx,:);
      idx = idx + 1;
    end
end

scatter(u(:,1),u(:,2));

%numerical differencing
%reference coordinate
for ii = 1:N
    for jj = 1:N
        X = [v(ii) v(jj)];
        Fana = embedDefGrad(X,qf,qc,Xf, Xc);
        h = 0.001;
        Fnum = eye(2);
        for dim = 1:2
            Xp=X;
            Xp(dim) = Xp(dim) + h;
            up = embedDisp(Xp,qf,qc,Xf, Xc);
            Xm=X;
            Xm(dim) = Xm(dim) - h;
            um = embedDisp(Xm,qf,qc,Xf, Xc);
            Fnum(:,dim) = Fnum(:,dim)+(0.5/h)*(up-um)';
        end
        Fnum
        Fana
    end
end