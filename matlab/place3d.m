%@param input build volume (1280x800)
%@param p 2d print
%@param x0 integer coordinate for origin of print
function bv = place3d(bv0, p, x0)
bv = bv0;
ix = 800;
iy = 1280;
s = size(bv, 1);
sp = size(p);
zscale = 4;
h = zscale*sp(3);
if(s < h)
    bv{sp(3)} = zeros(ix,iy);
end
for ii = 1:h
  if(ii>s)
      bv{ii} = zeros(ix, iy);
  end
   bv{ii}(x0(1):x0(1)+sp(1)-1, x0(2):x0(2) + sp(2)-1) = p(:,:, ceil(ii/zscale));
end

end