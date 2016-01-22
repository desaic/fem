function o = makeBin(im, size, thresh)
N = imresize(im, size);
o = N;
o(N>thresh) = 1;
o(N<thresh) = 0;
end
