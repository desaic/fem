function saveArrInt2d(A, filename)
FILE = fopen(filename,'w');
fprintf(FILE, '%d ', size(A));
fprintf(FILE,'\n');
for ii = 1:size(A,1)
  fprintf(FILE, '%d ', A(ii,:));
  fprintf(FILE,'\n');
end
fclose(FILE);
end