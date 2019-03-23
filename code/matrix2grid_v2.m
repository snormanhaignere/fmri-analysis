function G = matrix2grid_v2(X, G)

% Assumes spatial dimension is the last dimension of the matrix. Useful because
% this is the format returned by grid2matrix.

X = permute(X, [ndims(X), 1:ndims(X)-1]);
G = matrix2grid(X,G);