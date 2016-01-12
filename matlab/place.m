function platform = place(p, a, x, y)
  platform = p;
  platform (y:(y+size(a,1))-1, x:(x+size(a,2))-1) = a(:,:);
end