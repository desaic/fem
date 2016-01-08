%@param input build volume (1280x800)
%@param p 2d print
%@param x0 integer coordinate for origin of print
%@param t thickness (25 microns per layer)
function bv = place2d(bv0, p, x0, t)
bv = bv0;
ix = 800;
iy = 1280;

s = size(bv, 1);
if(s<t)
    bv{t} = zeros(ix,iy);
end

for ii = 1:t
  if(ii>s)
      bv{ii} = zeros(ix, iy);
  end
  sp = size(p);
  bv{ii}(x0(1):x0(1)+sp(1)-1, x0(2):x0(2) + sp(2)-1) = p(:,:);
end

end