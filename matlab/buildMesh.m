function [e, x] = buildMesh(sx,sy,sz)
hexV=[0,0,0;
0,0,1;
0,1,0;
0,1,1;
1,0,0;
1,0,1;
1,1,0;
1,1,1];
x = zeros((sx+1)*(sy+1)*(sz+1),3);
e = zeros(sx*sy*sz,8);
maxn = max(sx, sy);
maxn = max(maxn, sz);
dx = 1.0/maxn;
cnt = 1;
for ii=1:sx+1
    for jj=1:sy+1
      for kk=1:sz+1
        x(cnt,:) = dx * [ii, jj, kk];
        cnt = cnt+1;
      end
    end
end
cnt = 1;
for ii = 1:sx
    for jj = 1:sy
        for kk = 1:sz
            eidx = (ii-1)*sy*sz+(jj-1)*sz+kk;
            cnt = cnt+1;
            for ll = 1:8
                e(eidx,ll) = VX(ii+hexV(ll,1), jj+hexV(ll,2), kk+hexV(ll,3), sy,sz);
            end
        end
    end
end
end

function idx = VX(ii,jj,kk,sy,sz)
  idx = (ii-1)*(sz+1) * (sy+1) + (jj-1)*(sz+1) + (kk);
end