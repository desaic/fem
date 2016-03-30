%loads lhs matrices for homogenization
%d list of DOF for corner, interior, edge, face, 
%P Full stiffness matrix
%A modified stiffness matrix
%U U matrix for rhs
%W  rhs
function [d, K, A, U, W] = loadBC(dfile, Kfile, Afile, Ufile, Wfile)
d = loadArrs(dfile);
K = loadCSR(Kfile);
A = loadCSR(Afile);
U = loadCSR(Ufile);
W = loadCSR(Wfile);
end