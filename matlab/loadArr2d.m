function a = loadArr2d(file)
    IN = fopen (file);
    N = fscanf(IN, '%d', 2);
    N = N';
    a = zeros(N(1),N(2));
    for ii = 1:N(1)
        for jj = 1:N(2)
           a(ii,jj) = fscanf(IN, '%f', 1);
        end
    end
    fclose(IN);
end