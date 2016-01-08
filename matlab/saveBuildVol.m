function saveBuildVol(bv, prefix)
N = size(bv, 1);
  for ii = 1:N
      filename = strcat(prefix, 'slice_', int2str(ii), '.png');
      imwrite(bv{ii}, filename);
  end
end