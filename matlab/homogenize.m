function [U, C] = homogenize(nelx, nely, nelz, xPhys)
%% MATERIAL PROPERTIES
E0 = 1;
nu = 0.45;
maxIter = 100;
%% PREPARE FINITE ELEMENT ANALYSIS
dim = 3;
nele = nelx*nely*nelz;
ndof = dim*(nelx+1)*(nely+1)*(nelz+1);
KE = lk_H8(nu);
nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1)
nodenrs = zeros(nelx+1, nely+1, nelz+1);
for ii = 1:(nelx + 1)*(nely + 1)*(nelz + 1)
    tmp = (ii-1) / ((nelx + 1)*(nely + 1));
    zi = floor(tmp);
    jj = ii - 1 - zi * (nelx + 1)*(nely + 1);
    xi = floor(jj/(nely+1));
    yi = jj - xi * (nely+1);
    nodenrs(yi+1, xi+1 , zi+1) = ii;
end
nodeids = reshape(nodegrd(1:end-1,1:end-1),nely*nelx,1);
nodeidz = 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1);
nodeids = repmat(nodeids,size(nodeidz))+repmat(nodeidz,size(nodeids));
edofVec = 3*nodeids(:)+1;
%   7 8
%   5 6
%4 3
%1 2
edofMat = repmat(edofVec,1,24)+ ...
    repmat([0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1 ...
    3*(nely+1)*(nelx+1)+[0 1 2 3*nely + [3 4 5 0 1 2] -3 -2 -1]],nele,1);
iK = reshape(kron(edofMat,ones(24,1))',24*24*nele,1);
jK = reshape(kron(edofMat,ones(1,24))',24*24*nele,1);

%% PERIODIC BOUNDARY CONDITIONS
e0 = eye(6);
ufixed = zeros(24,6);
U = zeros(3*(nely+1)*(nelx+1)*(nelz+1),6);
alldofs = (1:3*(nely+1)*(nelx+1)*(nelz+1));
%corner vertices
%4 3
%1 2
n1 = [nodenrs(end, [1,end], 1), nodenrs(1, [end,1], 1)...
    nodenrs(end, [1,end], end), nodenrs(1, [end,1], end)];
d1 = reshape([ (3*n1-2) ; (3*n1-1); 3*n1 ] , 1, 24);
%lower left faces
%                  y2
%                 /|\
%                  |
%  n3_1             --->y1
%        n3_2
%
n3 = [reshape(nodenrs(2:end-1, 1, 1:end),   1, (nelx-1)*(nelz+1)) ...
reshape(nodenrs([1 end], 1, 2:end-1), 1, 2*(nelz-1)) ...
reshape(nodenrs(2:end-1, 1, 1:end),   1, (nelx-1)*(nelz+1)) ...
reshape(nodenrs([1 end], 1, 2:end-1), 1, 2*(nelz-1)) ...
reshape(nodenrs(2:end-1, 1, 1:end),   1, (nelx-1)*(nelz+1)) ...
reshape(nodenrs([1 end], 1, 2:end-1), 1, 2*(nelz-1))];

d3 = reshape([(2*n3-1);2*n3],1,2*(nelx+nely-2));
n4 = [nodenrs(2:end-1,end)',nodenrs(1,2:end-1)];
d4 = reshape([(2*n4-1);2*n4],1,2*(nelx+nely-2));
d2 = setdiff(alldofs,[d1,d3,d4]);
%ordering of ufixed must be consistent with n1 and d1
for j = 1:6
  ufixed(3:4,j) =[e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[nelx;0];
  ufixed(7:8,j) = [e0(1,j),e0(3,j)/2;e0(3,j)/2,e0(2,j)]*[0;nely];
  ufixed(5:6,j) = ufixed(3:4,j)+ufixed(7:8,j);
end
wfixed = [repmat(ufixed(3:4,:),nely-1,1); repmat(ufixed(7:8,:),nelx-1,1)];
change = 1;
loop = 0;

%% FE-ANALYSIS
sK = reshape(KE(:) * xPhys(:)', 64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
max(diag(K))

Kr = [K(d2,d2), K(d2,d3)+K(d2,d4); K(d3,d2)+K(d4,d2), K(d3,d3)+K(d4,d3)+K(d3,d4)+K(d4,d4)];
U(d1,:) = ufixed;
U([d2,d3],:) = Kr\(-[K(d2,d1); K(d3,d1)+K(d4,d1)]*ufixed-[K(d2,d4); K(d3,d4)+K(d4,d4)]*wfixed);
U(d4,:) = U(d3,:)+wfixed;
%% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
for i = 1:3
for j = 1:3
  U1 = U(:,i); U2 = U(:,j);
  qe{i,j} = reshape(sum((U1(edofMat)*KE).*U2(edofMat),2),nely,nelx)/(nelx*nely);
  Q(i,j) = sum(sum((xPhys).*qe{i,j}));
end
end
C=Q;
nu21 = Q(1,2)/Q(1,1)
nu12 = Q(2,1)/Q(2,2)
E1 = Q(1,1) * (1-nu12*nu21)
E2 = Q(2,2) * (1-nu12*nu21)
end

function [KE] = lk_H8(nu)
A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8;
    -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
k = 1/144*A'*[1; nu];

K1 = [k(1) k(2) k(2) k(3) k(5) k(5);
    k(2) k(1) k(2) k(4) k(6) k(7);
    k(2) k(2) k(1) k(4) k(7) k(6);
    k(3) k(4) k(4) k(1) k(8) k(8);
    k(5) k(6) k(7) k(8) k(1) k(2);
    k(5) k(7) k(6) k(8) k(2) k(1)];
K2 = [k(9)  k(8)  k(12) k(6)  k(4)  k(7);
    k(8)  k(9)  k(12) k(5)  k(3)  k(5);
    k(10) k(10) k(13) k(7)  k(4)  k(6);
    k(6)  k(5)  k(11) k(9)  k(2)  k(10);
    k(4)  k(3)  k(5)  k(2)  k(9)  k(12)
    k(11) k(4)  k(6)  k(12) k(10) k(13)];
K3 = [k(6)  k(7)  k(4)  k(9)  k(12) k(8);
    k(7)  k(6)  k(4)  k(10) k(13) k(10);
    k(5)  k(5)  k(3)  k(8)  k(12) k(9);
    k(9)  k(10) k(2)  k(6)  k(11) k(5);
    k(12) k(13) k(10) k(11) k(6)  k(4);
    k(2)  k(12) k(9)  k(4)  k(5)  k(3)];
K4 = [k(14) k(11) k(11) k(13) k(10) k(10);
    k(11) k(14) k(11) k(12) k(9)  k(8);
    k(11) k(11) k(14) k(12) k(8)  k(9);
    k(13) k(12) k(12) k(14) k(7)  k(7);
    k(10) k(9)  k(8)  k(7)  k(14) k(11);
    k(10) k(8)  k(9)  k(7)  k(11) k(14)];
K5 = [k(1) k(2)  k(8)  k(3) k(5)  k(4);
    k(2) k(1)  k(8)  k(4) k(6)  k(11);
    k(8) k(8)  k(1)  k(5) k(11) k(6);
    k(3) k(4)  k(5)  k(1) k(8)  k(2);
    k(5) k(6)  k(11) k(8) k(1)  k(8);
    k(4) k(11) k(6)  k(2) k(8)  k(1)];
K6 = [k(14) k(11) k(7)  k(13) k(10) k(12);
    k(11) k(14) k(7)  k(12) k(9)  k(2);
    k(7)  k(7)  k(14) k(10) k(2)  k(9);
    k(13) k(12) k(10) k(14) k(7)  k(11);
    k(10) k(9)  k(2)  k(7)  k(14) k(7);
    k(12) k(2)  k(9)  k(11) k(7)  k(14)];
KE = 1/((nu+1)*(1-2*nu))*...
    [ K1  K2  K3  K4;
    K2'  K5  K6  K3';
    K3' K6  K5' K2';
    K4  K3  K2  K1'];
end
