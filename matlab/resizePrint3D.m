function po = resizePrint3D(a, res)
r = 100;
if(nargin>=2)
 r = res;
end
thresh = 0.5;

nx = floor((size(a,2))^(1/3) + 0.5);
ny = nx;
nz = nx;
p = reshape(a, [nx ny nz]);
rep = 4;
po = zeros(r*rep, r*rep, r*rep);
for ir = 1:rep
for ii = 1:r
  tmp1 = nx*(ii-1);
  tmp = floor(tmp1/r);
  z = 1 + tmp;
  z=max(1,z);
  z = min(nz, z);
  if(z==2)
      z
  end
  layer = imresize(p(:,:, z) , [r r]);
  layer (layer >=thresh) = 1;
  layer (layer <thresh) = 0;
  layer = repmat(layer,4,4);
  po(:,:, (ir-1)*r + ii) = layer(:,:);
end
end
end