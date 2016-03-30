%draw displacement of a regular grid
%sx sy sz are sizes of the grid.
%i.e. sx+1 is number of points in 
%x direction.
function draw_disp(u, sx, sy, sz)
nEdges = 12 * sx * sy * sz;
coord = zeros(nEdges,2,3);
plot3(coord(:,:,1), coord(:,:,2), coord(:,:,3));
cnt = 1;
nvert = (sx+1) * (sy+1) * (sz+1);
x = zeros(nvert,3);
for ii = 1:sx+1
    for jj = 1:sy+1
        for kk = 1:sz+1
            x(linearIdx(ii,jj,kk,sx+1, sy+1, sz+1),:) = [ii,jj,kk];
        end
    end
end
end

function y = linearIdx(ii,jj,kk,si,sj,sk)
  y = (kk-1) * si*sj + (ii-1) * sj + jj;
end