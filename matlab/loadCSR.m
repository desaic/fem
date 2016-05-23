function K = loadCSR(filename)
IN = fopen(filename);
K=0;
N = fscanf(IN, '%d',[2 1]);
triplets = fscanf(IN, '%f', [3 inf]);
triplets = triplets';
N
size(triplets)
K=sparse(triplets(:,2)+1, triplets(:,1)+1, triplets(:,3),N(1), N(2));
fclose(IN);
end