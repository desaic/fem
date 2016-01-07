function a = loadArr1d(file)
    IN = fopen (file);
    N = fscanf(IN, '%d', 1);
    N = N';
    a = zeros(N(1),1);
    for ii = 1:N(1)
      a(ii) = fscanf(IN, '%f', 1);
    end
    fclose(IN);
end