function a = loadArr3d(file)
    IN = fopen (file);
    N = fscanf(IN, '%d', 3);
    N = N';
    a = zeros(N(1),N(2),N(3));
    for ii = 1:N(1)
        for jj = 1:N(2)
            for kk = 1:N(3)
                a(ii,jj,kk) = fscanf(IN, '%f', 1);    
            end
        end
    end
    fclose(IN);
end