function K = loadSparse(filename)
IN = fopen(filename);
N = fscanf(IN, '%d', 1);
nnz = fscanf(IN, '%d', 1);
M = fscanf(IN, '%f', [3 nnz]);
K=sparse(M(1,:), M(2,:), M(3,:), N,N);
fclose(IN);
end