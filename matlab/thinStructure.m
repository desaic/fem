function b = thinStructure(a, scale, thresh)
sa = size(a);
b = imresize(a, scale*sa);
b=imgaussfilt(b, 2);
b(b>=thresh) = 1;
b(b<=thresh) = 0;
end