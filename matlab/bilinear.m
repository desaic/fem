%@param x 2D point in natural cooridnates where the square is [-1,1]^2
%@return weights
% y     2      4
% ^ 
% |
%  -> x 1      3
function w = bilinear(x)
    w = zeros(4,1);
    w(1) = 0.25* (1-x(1)) * (1-x(2));
    w(2) = 0.25* (1-x(1)) * (1+x(2));
    w(3) = 0.25* (1+x(1)) * (1-x(2));
    w(4) = 0.25* (1+x(1)) * (1+x(2));
end