function p = resizePrints(a, res, rep)
r = 100;
rp = 4;
if(nargin>=2)
 r = res;
end
if(nargin>=3)
    rp = rep;
end

nx = int32(sqrt(size(a,2)));
ny = nx;
p = reshape(a, [nx ny]);
p = imresize(p, [r r]);
thresh = 0.5;
p(p>=thresh) = 1;
p(p<thresh) = 0;
p = repmat(p,rp,rp);
end