function [NC L] = mergeCells(A)
[ny nx] = size(A);
NC = 0;
L = A;
for ii = 1:ny
  startIdx = 0;
  endIdx = 1;
  prev = 0;
  for jj = 1:nx
    if(A(ii,jj) > 0)
      if(prev==0)
        NC = NC + 1;
        startIdx = jj;
        L(ii,jj) = NC;
      end
      endIdx = jj;
    else
      
    end
    prev = A(ii,jj);  
  end
end
end