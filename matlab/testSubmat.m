function testSubmat()
N = 36;
M = zeros(N,N);
cnt = 1;
for ii = 1:N
  j0 = mod(ii + 1 , 2)+1;
  for jj = j0:2: N
    M(jj, ii) = cnt;
    cnt = cnt + 1;
  end
end
M
block_size = 3;
subsize = 4;
n1 = zeros(subsize, 1);
n2 = zeros(subsize, 1);
for ii = 1: subsize
    n1(ii) = min((ii-1) * 4+1, floor(N / block_size) );
    n2(ii) = min((ii-1) * 4 + 2, floor(N / block_size));
end
n1e = expand(n1, block_size);
n2e = expand(n2, block_size);
sub = M(n1e, n2e)
end

function b = expand(a, block_size)
l = length(a);
b = zeros(length(a) * block_size,1);
for ii = 1:l
    for jj = 1:block_size
        b(3*(ii-1) + jj) = 3*(a(ii)-1) + jj;
    end
end
b
end