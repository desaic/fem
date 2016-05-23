function d = loadArrs(filename)
IN = fopen(filename, 'r');
nArr = fscanf(IN, '%d', 1);
d = cell(nArr,1);
for i = 1:nArr
    nEle = fscanf(IN,'%d', 1);
    d{i} = fscanf(IN, '%f', [nEle,1]) + 1;
end
fclose(IN);
end