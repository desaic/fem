%@param a 2D array
%@param b 3D structure
%@param rep repeat how many layers
function v = replaceVoxels(a, b, rep)
  sa = size(a);
  sb = size(b);
  v = zeros([sa(1) * sb(1)  sa(2) * sb(2) rep * sb(3)]);
  for ix = 1:sa(1)
      for iy = 1:sa(2)
          if(a(ix,iy)<=0)
              continue
          end
          for iz = 1:rep
            v((ix-1) * sb(1) + 1 : ix* sb(1) , (iy-1) * sb(2) + 1 : iy * sb(2), (iz-1) * sb(3) + 1 : iz*sb(3)) = ...
                  b(:,:,:);
          end
      end
  end
end