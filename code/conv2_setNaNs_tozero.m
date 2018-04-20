function Y = conv2_setNaNs_tozero(X,H)
% function Y = conv2_setNaNs_tozero(X,H)
% 
% 2-D convolution of matrix X with matrix H, after setting
% NaNs in X to zero.

if any(isnan(H(:)))
    error('No NANs allowed in second input matrix.');
end

X(isnan(X)) = 0;
Y = conv2(X,H,'same');