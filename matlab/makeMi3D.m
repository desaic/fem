%@param t beam thickness
function v = makeMi3D(len, t)
  v = zeros(len, len, len);
  mid = ceil(len/2);
  for ii = 1:len
      for jj = 1:t
          j1 = jj - ceil(t/2);
          for kk = 1:t
              k1 = kk - ceil(t/2);
              v(ii, mid+j1, mid+k1) = 1;
              v(mid+j1, ii, mid+k1) = 1;
              v(mid+j1, mid+k1, ii) = 1;
          end
      end
  end
end